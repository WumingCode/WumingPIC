module app
  use iso_fortran_env, only: int64
  use mpi
  use wuming3d
  use wuming_utils
  use boundary_shock, bc__init        => boundary_shock__init,        &
                      bc__dfield      => boundary_shock__dfield,      &
                      bc__particle_yz => boundary_shock__particle_yz, &
                      bc__injection   => boundary_shock__injection,   &
                      bc__curre       => boundary_shock__curre,       &
                      bc__phi         => boundary_shock__phi,         &
                      bc__mom         => boundary_shock__mom
  implicit none
  private

  ! main simulation loop
  public:: app__main

  ! configuration file and initial parameter files
  character(len=*), parameter   :: config_default = 'config.json'
  character(len=:), allocatable :: config_string
  character(len=:), allocatable :: config
  character(len=*), parameter   :: param = 'init_param'

  ! read from "config" section
  logical                       :: restart
  character(len=128)            :: restart_file
  character(len=:), allocatable :: datadir
  real(8)                       :: max_elapsed
  integer                       :: max_it
  integer                       :: intvl_ptcl
  integer                       :: intvl_mom
  integer                       :: intvl_orb
  integer                       :: intvl_expand
  integer                       :: verbose

  ! read from "parameter" section
  integer :: num_process
  integer :: num_process_j
  integer :: n_ppc
  integer :: n_x
  integer :: n_x_ini
  integer :: n_y
  integer :: n_z
  real(8) :: u_inject
  real(8) :: mass_ratio
  real(8) :: sigma_e
  real(8) :: omega_pe
  real(8) :: v_the
  real(8) :: v_thi
  real(8) :: theta_bn
  real(8) :: phi_bn
  real(8) :: l_damp_ini

  integer :: nproc, nproc_j, nproc_k
  integer :: it0
  integer :: np
  integer :: n0
  integer :: nx, nxgs, nxge, nxs, nxe
  integer :: ny, nygs, nyge
  integer :: nz, nzgs, nzge
  integer :: mpierr

  integer, parameter :: ndim   = 7
  integer, parameter :: nsp    = 2

  ! OTHER CONSTANTS
  real(8), parameter :: c      = 1d0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 5.01d-1   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 1d0     !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: delx   = 1d0     !CELL WIDTH
  real(8), parameter :: pi     = 4d0*atan(1d0)

  ! TRACKING PARTICLES INITIALLY XRS <= X <= XRE ONLY WHEN NDIM=7
  real(8), parameter :: xrs   = 4.5d2
  real(8), parameter :: xre   = 5d2

  !
  ! main variables
  !
  integer, allocatable, public :: np2(:,:,:), cumcnt(:,:,:,:)
  real(8), allocatable, public :: uf(:,:,:,:)
  real(8), allocatable, public :: up(:,:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:,:)
  real(8), allocatable, public :: mom(:,:,:,:,:)
  real(8)                      :: r(nsp)
  real(8)                      :: q(nsp)
  real(8)                      :: delt
  real(8)                      :: b0
  real(8)                      :: u0, v0, gam0

contains
  !
  ! main simulation loop
  !
  subroutine app__main()
    implicit none
    integer :: it
    real(8) :: etime, etime0

    ! initialization
    call load_config()
    call init()

    ! current clock
    etime0 = get_etime()

    ! main loop
    do it = it0+1, max_it

      ! update
      call particle__solv(gp, up, uf, cumcnt, nxs, nxe)
      call bc__injection(gp, np2, nxs, nxe, u0)
      call field__fdtd_i(uf, up, gp, cumcnt, nxs, nxe, &
                         bc__dfield, bc__curre, bc__phi)
      call bc__particle_yz(gp, np2)
      call sort__bucket(up, gp, cumcnt, np2, nxs, nxe)

      ! injection
      call inject()

      ! expand box
      if( mod(it, intvl_expand) == 0 ) then
        call relocate()
      endif

      ! output entire particles
      if ( mod(it, intvl_ptcl) == 0 ) then
        call io__ptcl(up, uf, np2, it)
      endif

      ! ouput tracer particles
      if ( mod(it, intvl_orb) == 0 ) then
        call io__orb(up, uf, np2, it)
      endif

      ! output moments and electromagnetic fields
      if ( mod(it, intvl_mom) == 0 ) then
        call mom_calc__accl(gp, up, uf, cumcnt, nxs, nxe)
        call mom_calc__nvt(mom, gp, np2)
        call bc__mom(mom)
        call io__mom(mom, uf, it)
      endif

      ! check elapsed time
      etime = get_etime() - etime0
      if ( etime >= max_elapsed ) then
        ! save snapshoft for restart
        write(restart_file, '(i7.7, "_restart")') it
        call save_restart(up, uf, np2, nxs, nxe, it, restart_file)

        if ( nrank == 0 ) then
          write(0,'("*** Elapsed time limit exceeded ")')
          write(0,'("*** A snapshot ", a, " has been saved")') trim(restart_file)
        endif

        call finalize()
        stop
      endif

      if( verbose >= 1 .and. nrank == 0 ) then
        write(*,'("*** Time step: ", i7, " completed in ", e10.2, " sec.")') it, etime
      endif
    enddo

    ! save final state
    it = max_it + 1
    write(restart_file, '(i7.7, "_restart")') it
    call save_restart(up, uf, np2, nxs, nxe, it, restart_file)
    call finalize()

  end subroutine app__main

  !
  ! parse command line argument and load configuration file
  !
  subroutine load_config()
    implicit none

    logical :: status, found
    integer :: arg_count
    character(len=:), allocatable :: filename

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    ! check ndim
    if( ndim /= 7 ) then
      write(0,*) 'Error: ndim must be 7'
      stop
    endif

    !
    ! process command line to find a configuration file
    !
    arg_count = command_argument_count()

    if( arg_count == 0 ) then
      config = config_default
    else
      ! only the first argument is relevant
      call get_command_argument(1, config)
    endif

    ! check init file
    inquire(file=trim(config), exist=status)

    if( .not. status ) then
      write(0, '("Error: ", a, " does not exists")') trim(config)
      stop
    endif

    ! try loading init file
    call json%initialize()
    call file%initialize()
    call file%load(trim(config))

    if( file%failed() ) then
      write(0, '("Error: failed to load ", a)') trim(config)
      stop
    endif

    ! read "config" section
    call file%get(root)
    call json%get(root, 'config', p)

    call json%get(p, 'verbose', verbose)
    call json%get(p, 'datadir', datadir)
    call json%get(p, 'max_elapsed', max_elapsed)
    call json%get(p, 'max_it', max_it)
    call json%get(p, 'intvl_ptcl', intvl_ptcl)
    call json%get(p, 'intvl_mom', intvl_mom)
    call json%get(p, 'intvl_orb', intvl_orb)
    call json%get(p, 'intvl_expand', intvl_expand)

    ! make sure this is a directory
    datadir = trim(datadir) // '/'

    ! restart file
    call json%get(p, 'restart_file', filename, found)

    if ( found .and. filename /= '' ) then
      restart = .true.
      restart_file = filename
    endif

    ! read "parameter" section and initialize
    call json%get(root, 'parameter', p)
    call json%get(p, 'num_process', num_process)
    call json%get(p, 'num_process_j', num_process_j)
    call json%get(p, 'n_ppc', n_ppc)
    call json%get(p, 'n_x', n_x)
    call json%get(p, 'n_y', n_y)
    call json%get(p, 'n_z', n_z)
    call json%get(p, 'n_x_ini', n_x_ini)
    call json%get(p, 'u_inject', u_inject)
    call json%get(p, 'mass_ratio', mass_ratio)
    call json%get(p, 'sigma_e', sigma_e)
    call json%get(p, 'omega_pe', omega_pe)
    call json%get(p, 'v_the', v_the)
    call json%get(p, 'v_thi', v_thi)
    call json%get(p, 'theta_bn', theta_bn)
    call json%get(p, 'phi_bn', phi_bn)
    call json%get(p, 'l_damp_ini', l_damp_ini)

    nproc   = num_process
    nproc_j = num_process_j
    nproc_k = nproc/nproc_j
    nx      = n_x
    ny      = n_y
    nz      = n_z
    np      = n_ppc*nx*5
    n0      = n_ppc
    nxgs    = 2
    nxge    = nxgs+nx-1
    nygs    = 2
    nyge    = nygs+ny-1
    nzgs    = 2
    nzge    = nzgs+nz-1
    nxs     = nxgs
    nxe     = nxgs+n_x_ini

    ! degree to radian
    theta_bn = theta_bn*pi/1.8d2
    phi_bn   = phi_bn  *pi/1.8d2

    call json%serialize(root, config_string)
    call json%destroy()
    call file%destroy()

  end subroutine load_config

  !
  ! initialize simulation
  !
  subroutine init()
    implicit none
    integer :: isp, i, j, k
    real(8) :: wpe, wpi, wge, wgi, vte, vti

    ! MPI
    call mpi_set__init(nygs, nyge, nzgs, nzge, nproc, nproc_j, nproc_k)

    ! random number
    call init_random_seed()

    ! allocate memory and initialize everything by zero
    allocate(np2(nys:nye,nzs:nze,nsp))
    allocate(cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp))
    allocate(uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
    allocate(up(ndim,np,nys:nye,nzs:nze,nsp))
    allocate(gp(ndim,np,nys:nye,nzs:nze,nsp))
    allocate(mom(1:7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp))
!$OMP PARALLEL WORKSHARE
    np2(nys:nye,nzs:nze,1:nsp) = 0
    cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,1:nsp) = 0
    uf(1:6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2) = 0d0
    up(1:ndim,1:np,nys:nye,nzs:nze,1:nsp) = 0d0
    gp(1:ndim,1:np,nys:nye,nzs:nze,1:nsp) = 0d0
    mom(1:7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = 0d0
!$OMP END PARALLEL WORKSHARE

    ! set physical parameters
    delt = cfl*delx/c
    u0   =-abs(u_inject)
    gam0 = sqrt(1d0+(u0*u0)/(c*c))
    v0   = u0/gam0
    wpe  = omega_pe
    wge  = omega_pe*sqrt(sigma_e)
    wpi  = wpe/sqrt(mass_ratio)
    wgi  = wge/mass_ratio
    vte  = v_the
    vti  = v_thi
    r(1) = mass_ratio
    r(2) = 1.0d0
    q(1) =+sqrt(gam0*r(1)/(4*pi*n0))*wpi
    q(2) =-sqrt(gam0*r(2)/(4*pi*n0))*wpe
    b0   = r(1)*c/q(1)*wgi*gam0

    ! number of particles
!$OMP PARALLEL WORKSHARE    
    np2(nys:nye,nzs:nze,1:nsp) = n0*(nxe-nxs-1)
!$OMP END PARALLEL WORKSHARE
    if ( nrank == 0 ) then
      if ( n0*(nxge-nxgs-1) > np ) then
        write(0,*) 'Error: Too large number of particles'
        stop
      endif
    endif

    ! preparation of sort
    do isp = 1,nsp
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = nzs,nze
      do j = nys,nye
        cumcnt(nxs:nxs+1,j,k,isp) = 0
        do i = nxs+2, nxe
          cumcnt(i,j,k,isp) = cumcnt(i-1,j,k,isp) + n0
        enddo
        if ( cumcnt(nxe,j,k,isp) /= np2(j,k,isp) ) then
          write(0,*) 'Error: invalid values encounterd for cumcnt'
          stop
        endif
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

    ! initialize modules
    call bc__init(ndim,np,nsp,                                    &
                  nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze,  &
                  jup,jdown,kup,kdown,mnpi,mnpr,ncomw,nerr,nstat, &
                  delx,delt,c)
    call particle__init(ndim,np,nsp,                                   &
                        nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                        delx,delt,c,q,r)
    call field__init(ndim,np,nsp,                                   &
                     nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                     mnpr,ncomw,opsum,nerr,                         &
                     delx,delt,c,q,r,gfac)
    call sort__init(ndim,np,nsp, &
                    nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze)
    call io__init(ndim,np,nsp,                                   &
                  nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                  nproc,nproc_j,nproc_k,nrank,delx,delt,c,q,r,datadir)
    call mom_calc__init(ndim,np,nsp,                                   &
                        nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                        delx,delt,c,q,r)

    if ( restart ) then
      ! restart
      call io__input(gp,uf,np2,nxs,nxe,it0,restart_file)
      call sort__bucket(up,gp,cumcnt,np2,nxs,nxe)
    else
      ! output parameters and set initial condition
      call save_param(n0,wpe,wpi,wge,wgi,vti,vte,param)
      call set_initial_condition()
      it0 = 0
    endif

    ! copy
    gp = up

  end subroutine init

  !
  ! finalize simulation
  !
  subroutine finalize()
    implicit none

    call io__finalize()
    call MPI_Finalize(mpierr)

  end subroutine finalize

  !
  ! set initial condition for field and particles
  !
  subroutine set_initial_condition()
    implicit none
    integer :: i, j, k, ii, isp
    real(8) :: v1, gam1, gamp, sd(nsp)

    !
    ! electromagnetic field
    !
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
      do i=nxgs-2,nxge+2
        uf(1,i,j,k) = b0*cos(theta_bn)
        uf(2,i,j,k) = b0*sin(theta_bn)*cos(phi_bn)
        uf(3,i,j,k) = b0*sin(theta_bn)*sin(phi_bn)
        uf(4,i,j,k) = 0d0
        uf(5,i,j,k) =+vprofile(i*delx)*uf(3,i,j,k)/c
        uf(6,i,j,k) =-vprofile(i*delx)*uf(2,i,j,k)/c
      enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !
    ! particle position
    !
    isp = 1
!$OMP PARALLEL DO PRIVATE(ii,j,k)
    do k=nzs,nze
    do j=nys,nye
      do ii=1,np2(j,k,isp)
        up(1,ii,j,k,1) = (nxs+(nxe-nxs)*(ii-5d-1)/np2(j,k,isp))*delx
        up(2,ii,j,k,1) = (j+uniform_rand())*delx
        up(3,ii,j,k,1) = (k+uniform_rand())*delx
        up(1,ii,j,k,2) = up(1,ii,j,k,1)
        up(2,ii,j,k,2) = up(2,ii,j,k,1)
        up(3,ii,j,k,2) = up(3,ii,j,k,1)
      enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !
    ! particle velocity
    !
    sd(1) = v_thi
    sd(2) = v_the
    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(ii,j,k,v1,gam1,gamp)
      do k=nzs,nze
      do j=nys,nye
        do ii=1,np2(j,k,isp)
          ! Maxwellian in fluid rest frame
          up(4,ii,j,k,isp) = sd(isp)*normal_rand()
          up(5,ii,j,k,isp) = sd(isp)*normal_rand()
          up(6,ii,j,k,isp) = sd(isp)*normal_rand()

          ! Lorentz transform to lab frame
          v1   = vprofile(up(1,ii,j,k,isp))
          gam1 = 1d0/sqrt(1d0-(v1/c)**2)
          gamp = sqrt(1d0+(up(4,ii,j,k,isp)**2 &
                          +up(5,ii,j,k,isp)**2 &
                          +up(6,ii,j,k,isp)**2 )/c**2)
          up(4,ii,j,k,isp) = gam1*(up(4,ii,j,k,isp)+v1*gamp)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

    ! set particle IDs
    call set_particle_ids()

  end subroutine set_initial_condition

  !
  ! set initial particle IDs
  !
  subroutine set_particle_ids()
    implicit none
    integer :: isp, i, j, k

    integer(8) :: gcumsum(nproc+1,nsp), lcumsum(nys:nye,nzs:nze,nsp), pid

    if ( ndim /= 7) then
      return
    end if

    ! calculate the first particle IDs
    call get_global_cumsum(np2,gcumsum)

    do isp = 1, nsp
      lcumsum(nys,nzs,isp) = gcumsum(nrank+1,isp)
      do k=nzs,nze-1
        do j=nys,nye-1
          lcumsum(j+1,k,isp) = lcumsum(j,k,isp)+np2(j,k,isp)
        enddo
        lcumsum(nys,k+1,isp) = lcumsum(nye,k,isp)+np2(nye,k,isp)
      enddo

      do j=nys,nye-1
        lcumsum(j+1,nze,isp) = lcumsum(j,nze,isp)+np2(j,nze,isp)
      enddo
    enddo

    ! unique ID as 64bit integer (negative by default)
    do isp = 1, nsp
!$OMP PARALLEL DO PRIVATE(i,j,k,pid)
      do k=nzs,nze
      do j=nys,nye
        do i=1,np2(j,k,isp)
          pid = lcumsum(j,k,isp)+i
          up(7,i,j,k,isp) = transfer(-pid,1d0)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

    ! make particle ID positive for output
    do isp = 1, nsp
!$OMP PARALLEL DO PRIVATE(i,j,k,pid)
      do k=nzs,nze
      do j=nys,nye
        do i=1,np2(j,k,isp)
          if ( up(1,i,j,k,isp) >= xrs .and. up(1,i,j,k,isp) <= xre ) then
            pid = transfer(up(7,i,j,k,isp), 1_8)
            up(7,i,j,k,isp) = transfer(sign(pid, +1_8), 1.0_8)
          endif
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

  end subroutine set_particle_ids

  !
  ! save everything for restart
  !
  subroutine save_restart(up, uf, np2, nxs, nxe, it, restart_file)
    implicit none
    integer, intent(in)          :: np2(nys:nye,nzs:nze,nsp), nxs, nxe
    integer, intent(in)          :: it
    real(8), intent(in)          :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)          :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    character(len=*), intent(in) :: restart_file

    logical :: found
    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    call json%initialize()
    call file%initialize()
    call file%deserialize(config_string)
    call file%get(root)
    call json%get(root, 'config', p)
    call json%update(p, 'restart_file', trim(restart_file), found)

    ! write data to the disk
    call io__output(up, uf, np2, nxs, nxe, it, trim(restart_file))

    if ( nrank == 0 ) then
       call json%print(root, config)
    end if

    call json%destroy()
    call file%destroy()

  end subroutine save_restart

  !
  ! save parameters
  !
  subroutine save_param(n0, wpe, wpi, wge, wgi, vti, vte, filename)
    implicit none
    integer, intent(in)          :: n0
    real(8), intent(in)          :: wpe, wpi, wge, wgi, vti, vte
    character(len=*), intent(in) :: filename

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp
    integer :: fh

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    ! save deafult parameters
    call io__param(n0, wpe, wpi, wge, wgi, vti, vte, filename)

    ! save additional parameters
    datafile = trim(datadir) // trim(filename) // '.raw'
    jsonfile = trim(datadir) // trim(filename) // '.json'

    ! open json file
    call file%initialize()
    call json%initialize()
    call file%load(jsonfile)
    call file%get(root)

    ! open data file
    call mpiio_open_file(datafile, fh, disp, 'a')

    ! put attributes
    call json%get(root, 'attribute', p)

    call jsonio_put_attribute(json, p, u0, 'u0', disp, '')
    call mpiio_write_atomic(fh, disp, u0)

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine save_param


  subroutine relocate()
    implicit none
    integer :: isp, j, k, ii, ii1 ,ii2
    real(8) :: v1, gam1, gamp, sd(nsp)

    integer(8) :: gcumsum(nproc+1,nsp), nptotal(nsp), pid

    if(nxe==nxge) return

    ! expand box
    nxe = nxe+1

    ! get particle number for ID
    call get_global_cumsum(np2, gcumsum)
    nptotal(1:nsp) = gcumsum(nproc+1,1:nsp)

    !
    ! position
    !
!$OMP PARALLEL DO PRIVATE(ii,ii1,ii2,j,k)
    do k=nzs,nze
    do j=nys,nye
      do ii=1,n0
        ii1 = np2(j,k,1) + ii
        ii2 = np2(j,k,2) + ii

        up(1,ii1,j,k,1) = (nxe-1)*delx+(ii-5d-1)/n0*delx
        up(2,ii1,j,k,1) = (j+uniform_rand())*delx
        up(3,ii1,j,k,1) = (k+uniform_rand())*delx
        up(1,ii2,j,k,2) = up(1,ii1,j,k,1)
        up(2,ii2,j,k,2) = up(2,ii1,j,k,1)
        up(3,ii2,j,k,2) = up(3,ii1,j,k,1)
      enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !
    ! velocity
    !
    sd(1) = v_thi
    sd(2) = v_the
    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(ii,j,k,v1,gam1,gamp,pid)
      do k=nzs,nze
      do j=nys,nye
        do ii=np2(j,k,isp)+1,np2(j,k,isp)+n0
          ! Maxwellian in fluid rest frame
          up(4,ii,j,k,isp) = sd(isp)*normal_rand()
          up(5,ii,j,k,isp) = sd(isp)*normal_rand()
          up(6,ii,j,k,isp) = sd(isp)*normal_rand()

          ! Lorentz transform to lab frame
          v1   = vprofile(up(1,ii,j,k,isp))
          gam1 = 1d0/sqrt(1d0-(v1/c)**2)
          gamp = sqrt(1d0+(up(4,ii,j,k,isp)**2 &
                          +up(5,ii,j,k,isp)**2 &
                          +up(6,ii,j,k,isp)**2)/c**2)
          up(4,ii,j,k,isp) = gam1*(up(4,ii,j,k,isp)+v1*gamp)

          ! particle ID
          pid = ii-np2(j,k,isp)+((nyge-nygs+1)*(k-nzgs)+(j-nygs))*n0+nptotal(isp)
          up(7,ii,j,k,isp) = transfer(-pid, 1.0_8)
        enddo
        np2(j,k,isp)        = np2(j,k,isp)         +n0
        cumcnt(nxe,j,k,isp) = cumcnt(nxe-1,j,k,isp)+n0
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
      uf(2,nxe-1,j,k) = b0*sin(theta_bn)*cos(phi_bn)
      uf(3,nxe-1,j,k) = b0*sin(theta_bn)*sin(phi_bn)
      uf(5,nxe-1,j,k) =+v0*uf(3,nxe-1,j,k)/c
      uf(6,nxe-1,j,k) =-v0*uf(2,nxe-1,j,k)/c
      uf(2,nxe,j,k)   = b0*sin(theta_bn)*cos(phi_bn)
      uf(3,nxe,j,k)   = b0*sin(theta_bn)*sin(phi_bn)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine relocate

  !
  ! injection from the right-hand side boundary
  !
  subroutine inject()
    implicit none
    integer :: isp, ii, ii1, ii2, i, j, k
    real(8) :: v1, gam1, gamp, sd(nsp)

    real(8) :: pflux, x0, xinj
    integer :: nginj, ngmod, nginj_proc(nproc), index_proc(nproc)
    integer :: nlinj, nlmod, nlinj_grid(nys:nye,nzs:nze), index_grid((nye-nys+1)*(nze-nzs+1)), ig2(2)
    integer :: ncinj_proc(nproc+1), ncinj_grid(nys:nye,nzs:nze)
    integer(8) :: gcumsum(nproc+1,nsp), nptotal(nsp), pid

    !
    ! Determine number of particle injected into each cell with
    ! the following steps.
    !
    ! (1) Determine the total number of particles injected.
    !     (A fractional portion is taken into acccount by random numbers.)
    ! (2) Equally divide it among PEs. Reminders are added randomly.
    ! (3) In each PE, equally divide the number of particles for each
    !     cell. Reminders are added randomly.
    !

    ! * number of particles injected into the entire system
    pflux = n0*abs(v0)*delt*delx*(nyge-nygs+1)*(nzge-nzgs+1)
    nginj = int(pflux)
    if( uniform_rand() < pflux-int(pflux) ) then
      nginj = nginj+1
    endif

    ! * number of particles injected into the local domain
    ngmod = mod(nginj,nproc)
    do i =1,nproc
      nginj_proc(i) = nginj/nproc
      index_proc(i) = i
    enddo
    call shuffle(index_proc)
    do i=1,ngmod
      nginj_proc(index_proc(i)) = nginj_proc(index_proc(i))+1
    enddo

    ! send from root to other processes
    call MPI_Bcast(nginj_proc, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

    ! * number of particles injected into each cell
    nlinj = nginj_proc(nrank+1)
    nlmod = mod(nlinj, (nye-nys+1)*(nze-nzs+1))
    do k=nzs,nze
    do j=nys,nye
      nlinj_grid(j,k) = nlinj/((nye-nys+1)*(nze-nzs+1))
      index_grid(j-(nys-1)+(nye-nys+1)*(k-nzs)) = j-(nys-1)+(nye-nys+1)*(k-nzs)
    enddo
    enddo

    call shuffle(index_grid)
    do i=1,nlmod
      ig2(1) = nys+mod(index_grid(i)-1,nye-nys+1)
      ig2(2) = nzs+(index_grid(i)-1)/(nye-nys+1)
      nlinj_grid(ig2(1),ig2(2)) = nlinj_grid(ig2(1),ig2(2))+1
    enddo

    !
    ! The following steps are needed for assigning particle IDs.
    !

    ! get current number of particles in the system
    call get_global_cumsum(np2, gcumsum)
    nptotal(1:nsp) = gcumsum(nproc+1,1:nsp)

    ! cumulative sum of number of injection particles
    ncinj_proc(1) = 0
    do i = 1, nproc
      ncinj_proc(i+1) = ncinj_proc(i)+nginj_proc(i)
    enddo

    ncinj_grid(nys,nzs) = ncinj_proc(nrank+1)
    do k = nzs,nze-1
      do j = nys,nye-1
        ncinj_grid(j+1,k) = ncinj_grid(j,k)+nlinj_grid(j,k)
      enddo
      ncinj_grid(nys,k+1) = ncinj_grid(nye,k)+nlinj_grid(nye,k)
    enddo

    do j=nys,nye-1
      ncinj_grid(j+1,nze) = ncinj_grid(j,nze)+nlinj_grid(j,nze)
    enddo

    !
    ! position
    !
    x0 = abs(v0)*delt
!$OMP PARALLEL DO PRIVATE(ii,ii1,ii2,j,k)
    do k=nzs,nze
    do j=nys,nye
      do ii = 1, nlinj_grid(j,k)
        ii1 = np2(j,k,1)+ii
        ii2 = np2(j,k,2)+ii

        up(1,ii1,j,k,1) = nxe*delx+(ii-5d-1)/nlinj_grid(j,k)*x0
        up(2,ii1,j,k,1) = (j+uniform_rand())*delx
        up(3,ii1,j,k,1) = (k+uniform_rand())*delx
        up(1,ii2,j,k,2) = up(1,ii1,j,k,1)
        up(2,ii2,j,k,2) = up(2,ii1,j,k,1)
        up(3,ii2,j,k,2) = up(3,ii1,j,k,1)
      enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !
    ! velocity
    !
    sd(1) = v_thi
    sd(2) = v_the
    do isp = 1, nsp
!$OMP PARALLEL DO PRIVATE(ii,j,k,xinj,v1,gam1,gamp,pid)
      do k=nzs,nze
      do j=nys,nye
        do ii = np2(j,k,isp)+1, np2(j,k,isp)+nlinj_grid(j,k)
          ! Maxwellian in fluid rest frame
          up(4,ii,j,k,isp) = sd(isp)*normal_rand()
          up(5,ii,j,k,isp) = sd(isp)*normal_rand()
          up(6,ii,j,k,isp) = sd(isp)*normal_rand()

          ! injection (non-relativistic approximation)
          xinj = up(1,ii,j,k,isp)+(v0+up(4,ii,j,k,isp))*delt
          up(1,ii,j,k,isp) = xinj
!         if( xinj <= nxe*delx ) then
!         ! leave ux as is
!           up(1,ii,j,k,isp) = xinj
!         else
!         ! folding (x, ux)
!           up(4,ii,j,k,isp) =-up(4,ii,j,k,isp)
!           up(1,ii,j,k,isp) =-up(1,ii,j,k,isp) + &
!                 & 2*x0 + (v0 + up(4,ii,j,k,isp)) * delt
!         end if

          ! Lorentz transform to lab frame
          v1   = vprofile(up(1,ii,j,k,isp))
          gam1 = 1d0/sqrt(1d0-(v1/c)**2)
          gamp = sqrt(1d0+(up(4,ii,j,k,isp)**2 &
                          +up(5,ii,j,k,isp)**2 &
                          +up(6,ii,j,k,isp)**2)/c**2)
          up(4,ii,j,k,isp) = gam1*(up(4,ii,j,k,isp)+v1*gamp)

          ! particle ID
          pid = ii-np2(j,k,isp)+ncinj_grid(j,k)+nptotal(isp)
          up(7,ii,j,k,isp) = transfer(-pid, 1.0_8)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

    do isp = 1, nsp
!$OMP PARALLEL WORKSHARE
      np2(nys:nye,nzs:nze,isp)        = np2(nys:nye,nzs:nze,isp)       +nlinj_grid(nys:nye,nzs:nze)
      cumcnt(nxe,nys:nye,nzs:nze,isp) = cumcnt(nxe,nys:nye,nzs:nze,isp)+nlinj_grid(nys:nye,nzs:nze)
!$OMP END PARALLEL WORKSHARE
    enddo

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
      uf(2,nxe-1,j,k) = b0*sin(theta_bn)*cos(phi_bn)
      uf(3,nxe-1,j,k) = b0*sin(theta_bn)*sin(phi_bn)
      uf(5,nxe-1,j,k) =+v0*uf(3,nxe-1,j,k)/c
      uf(6,nxe-1,j,k) =-v0*uf(2,nxe-1,j,k)/c
      uf(2,nxe,j,k)   = b0*sin(theta_bn)*cos(phi_bn)
      uf(3,nxe,j,k)   = b0*sin(theta_bn)*sin(phi_bn)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine inject

  !
  ! get global cumulative sum of particle numbers
  !
  subroutine get_global_cumsum(np2,cumsum)
    implicit none
    integer, intent(in)       :: np2(nys:nye,nzs:nze,nsp)
    integer(8), intent(inout) :: cumsum(nproc+1,nsp)

    integer :: i, isp, mpierr
    integer(8) :: lcount(nsp), gcount(nsp,nproc)

    ! get number of particles for each proces
    do isp=1,nsp
      lcount(isp) = sum(np2(nys:nye,nzs:nze,isp))
    enddo
    call MPI_Allgather(lcount,nsp,MPI_INTEGER8,gcount,nsp,MPI_INTEGER8, &
                       MPI_COMM_WORLD,mpierr)

    ! calculate cumulative sum
    do isp=1,nsp
      cumsum(1,isp) = 0
      do i=1,nproc
        cumsum(i+1,isp) = cumsum(i,isp)+gcount(isp,i)
      enddo
    enddo

  end subroutine get_global_cumsum

  !
  ! initial velocity profile
  !
  function vprofile(x) result(y)
    implicit none
    real(8), intent(in) :: x
    real(8) :: y
    real(8) :: x0, xs

    x0 = l_damp_ini+nxgs*delx
    xs = l_damp_ini*1d-1
    y  = 5d-1*v0*(1d0+tanh((x-x0)/xs))

  end function vprofile

end module app
