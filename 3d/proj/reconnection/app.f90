module app
  use iso_fortran_env, only: int64
  use mpi
  use wuming3d
  use wuming_utils
  use boundary_reconnection, &
       & bc__init        => boundary_reconnection__init,        &
       & bc__dfield      => boundary_reconnection__dfield,      &
       & bc__particle_x  => boundary_reconnection__particle_x,  &
       & bc__particle_yz => boundary_reconnection__particle_yz, &
       & bc__curre       => boundary_reconnection__curre,       &
       & bc__phi         => boundary_reconnection__phi,         &
       & bc__mom         => boundary_reconnection__mom
  implicit none
  private

  ! main simulation loop
  public:: app__main

  ! configuration file and initial parameter files
  character(len=*), parameter   :: config_default = 'config.json'
  character(len=:), allocatable :: config_string
  character(len=:), allocatable :: config
  character(len=*), parameter   :: param = 'init_param'
  character(len=*), parameter   :: ehist = 'energy.dat'

  ! read from "config" section
  logical                       :: restart
  character(len=128)            :: restart_file
  character(len=:), allocatable :: datadir
  real(8)                       :: max_elapsed
  integer                       :: max_it
  integer                       :: intvl_ptcl
  integer                       :: intvl_mom
  integer                       :: intvl_orb
  integer                       :: verbose

  ! read from "parameter" section
  integer :: num_process
  integer :: num_process_j
  integer :: n_ppc
  integer :: n_x
  integer :: n_y
  integer :: n_z
  real(8) :: mass_ratio
  real(8) :: alpha
  real(8) :: rtemp
  real(8) :: lcs
  integer :: nbg
  integer :: ncs

  integer :: nproc, nproc_j, nproc_k
  integer :: it0
  integer :: np
  integer :: n0
  integer :: nx, nxgs, nxge, nxs, nxe
  integer :: ny, nygs, nyge !, nys, nye are defined in mpi_set.f90
  integer :: nz, nzgs, nzge
  integer :: mpierr

  integer, parameter :: ndim   = 7
  integer, parameter :: nsp    = 2
  integer, parameter :: nroot  = 0

  ! OTHER CONSTANTS
  real(8), parameter :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 0.5D0     !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

  !
  ! main variables
  !
  integer, allocatable, public :: np2(:,:,:), cumcnt(:,:,:,:)!, np_harris(:,:,:,:)
  real(8), allocatable, public :: uf(:,:,:,:)
  real(8), allocatable, public :: up(:,:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:,:)
  real(8), allocatable, public :: mom(:,:,:,:,:)
  real(8)                      :: r(nsp)
  real(8)                      :: q(nsp)
  real(8)                      :: delt
  real(8)                      :: b0, vte, vti
  real(8)                      :: x0, y0, z0

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
       call bc__particle_x(gp, np2, nxs, nxe)
       call field__fdtd_i(uf, up, gp, cumcnt, nxs, nxe, &
            & bc__dfield, bc__curre, bc__phi)
       call bc__particle_yz(gp, np2)
       call sort__bucket(up, gp, cumcnt, np2, nxs, nxe)
      
       ! output entire particles
       if ( mod(it, intvl_ptcl) == 0 ) then
          call io__ptcl(up, uf, np2, it)
       end if

       ! output tracer particles
       if ( mod(it, intvl_orb) == 0 ) then
          call io__orb(up, uf, np2, it)
       end if

       ! output moments and electromagnetic fields
       if ( mod(it, intvl_mom) == 0 ) then
          call mom_calc__accl(gp, up, uf, cumcnt, nxs, nxe)
          call mom_calc__nvt(mom, gp, np2)
          call bc__mom(mom)
          call io__mom(mom, uf, it)
          call energy_history(up, uf, np2, it)
       endif

       ! check elapsed time
       etime = get_etime() - etime0
       if ( etime >= max_elapsed ) then
          ! save snapshoft for restart
          write(restart_file, '(i7.7, "_restart")') it
          call save_restart(up, uf, np2, nxs, nxe, it, restart_file)

          if ( nrank == nroot ) then
             write(0,'("*** Elapsed time limit exceeded ")')
             write(0,'("*** A snapshot ", a, " has been saved")') &
                  & trim(restart_file)
          end if

          call finalize()
          stop
       endif

       if( verbose >= 1 .and. nrank == nroot ) then
          write(*,'("*** Time step: ", i7, " completed in ", e10.2, " sec.")') &
               & it, etime
       end if
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
    end if

    !
    ! process command line to find a configuration file
    !
    arg_count = command_argument_count()

    if( arg_count == 0 ) then
       config = config_default
    else
       ! only the first argument is relevant
       call get_command_argument(1, config)
    end if

    ! check init file
    inquire(file=trim(config), exist=status)

    if( .not. status ) then
       write(0, '("Error: ", a, " does not exists")') trim(config)
       stop
    end if

    ! try loading init file
    call json%initialize()
    call file%initialize()
    call file%load(trim(config))

    if( file%failed() ) then
       write(0, '("Error: failed to load ", a)') trim(config)
       stop
    end if

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
    call json%get(p, 'n_x', n_x)
    call json%get(p, 'n_y', n_y)
    call json%get(p, 'n_z', n_z)
    call json%get(p, 'mass_ratio', mass_ratio)
    call json%get(p, 'alpha', alpha)
    call json%get(p, 'rtemp', rtemp)
    call json%get(p, 'lcs', lcs)
    call json%get(p, 'nbg', nbg)
    call json%get(p, 'ncs', ncs)

    nproc = num_process
    nproc_j = num_process_j
    nproc_k = nproc/nproc_j
    nx    = n_x
    ny    = n_y
    nz    = n_z
    n0    = nbg + ncs
    np    = n0 * nx
    nxgs  = 2
    nxge  = nxgs + nx - 1
    nygs  = 2
    nyge  = nygs + ny - 1
    nzgs  = 2
    nzge  = nzgs + nz - 1
    nxs   = nxgs
    nxe   = nxge

    call json%serialize(root, config_string)
    call json%destroy()
    call file%destroy()

  end subroutine load_config

  !
  ! initialize simulation
  !
  subroutine init()
    implicit none
    integer :: isp, i, j, k, ii
    real(8) :: wpe, wpi, wge, wgi, ldb

    ! MPI
    call mpi_set__init(nygs,nyge,nzgs,nzge,nproc,nproc_j,nproc_k)
    !call mpi_set__init(nygs, nyge, nproc)

    ! random number
    call init_random_seed()
    
    ! allocate memory and initialize everything by zero
    allocate(np2(nys:nye,nzs:nze,nsp))
    allocate(cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp))
    allocate(uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
    allocate(up(ndim,np,nys:nye,nzs:nze,nsp))
    allocate(gp(ndim,np,nys:nye,nzs:nze,nsp))
    allocate(mom(1:7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp))
    np2    = 0
    cumcnt = 0
    uf     = 0
    up     = 0
    gp     = 0
    mom    = 0

    ! set physical parameters
    r(1) = mass_ratio
    r(2) = 1.0d0
    delt = cfl*delx/c
    ldb  = delx  ! Debye length = dx
    vte  = sqrt(rtemp)*c/(sqrt(1+rtemp)*alpha) !from the pressure balance
    vti  = vte*sqrt(r(2)/r(1))/sqrt(rtemp)
    wpe  = vte/ldb/sqrt(2.d0)
    wpi  = wpe*sqrt(r(2)/r(1))
    wge  = wpe/alpha
    wgi  = wge/mass_ratio
    q(1) =+sqrt(r(1) / (4*pi*n0/delx**3)) * wpi
    q(2) =-sqrt(r(2) / (4*pi*n0/delx**3)) * wpe
    b0   = r(1)*c/q(1) * wgi

    ! POSITION OF THE X-POINT
    x0  = 0.5*(nxge+nxgs)*delx
    y0  = 0.5*(nyge-nygs)*delx
    z0  = 0.5*(nzge-nzgs)*delx
    ! CURRENT SHEET THICKNESS
    lcs = lcs * c/wpi
   
    ! NUMBER OF PARTICLES IN EACH CELL IN Y
    np2(nys:nye,nzs:nze,1:nsp) = nbg*(nxge-nxgs) + ncs*2*lcs
   
    ! initialize modules
    call bc__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze, &
         & jup, jdown, kup, kdown, mnpi, mnpr, ncomw, nerr, nstat, delx, delt, c)
    call particle__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze, &
         & delx, delt, c, q, r)
    call field__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze, &
         & mnpr, ncomw, opsum, nerr, delx, delt, c, q, r, gfac)
    call sort__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze)
    call io__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze, &
         & nproc, nproc_j, nproc_k, nrank, delx, delt, c, q, r, datadir)
    call mom_calc__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze, &
         & delx, delt, c, q, r)

    if ( restart ) then
       ! restart
       call io__input(gp, uf, np2, nxs, nxe, it0, restart_file)
       call sort__bucket(up, gp, cumcnt, np2, nxs, nxe)
    else
       ! output parameters and set initial condition
       call save_param(n0, wpe, wpi, wge, wgi, vti, vte, param)

       call set_initial_condition()
       ! sort loaded particles
       call sort__bucket(gp, up, cumcnt, np2, nxs, nxe)
       ! copy
       up = gp

       it0 = 0
       call energy_history(up, uf, np2, it0)
    endif

  end subroutine init

  !
  ! finalize simulation
  !
  subroutine finalize()
    implicit none
    write (*,*) "Checkpoint!!!" 
    call io__finalize()
    call MPI_Finalize(mpierr)

  end subroutine finalize

  !
  ! set initial condition for field and particles
  !
  subroutine set_initial_condition()
    implicit none
    integer            :: i, j, k, ii, ibg, isp
    real(8), parameter :: e1 = 0.12d0  ! B1/B0: SZ recommend ~< O(0.5*nbg/ncs).
    real(8)            :: r1, r2
    real(8)            :: jz, density
    real(8)            :: b_harris, bx_pert, by_pert
    real(8)            :: x, y
    real(8)            :: f1, f2
    real(8)            :: sdi, sde  ! for Box-Mullter method

    ! ---------------- Utility functions ------------------
    ! magnetic field strength
    b_harris(x) = b0 * tanh((x-x0)/lcs)
    ! localized perturbation
    bx_pert(x,y) = +e1*b0 *((y-y0)/lcs) * exp(-((x-x0)**2+(y-y0)**2)/(2*lcs)**2)
    by_pert(x,y) = -e1*b0 *((x-x0)/lcs) * exp(-((x-x0)**2+(y-y0)**2)/(2*lcs)**2)
    ! density
    density(x) = ncs * cosh((x-x0)/lcs)**(-2) + nbg
    ! current jz_0 > 0, while jz_1 < 0
    jz(x,y)    = &
         +     b0/(4*pi*lcs) * cosh((x-x0)/lcs)**(-2) &
         -2*e1*b0/(4*pi*lcs) * ( 1.d0-((x-x0)**2+(y-y0)**2)/(2*lcs)**2 ) * exp(-((x-x0)**2+(y-y0)**2)/(2*lcs)**2)

    !
    ! electromagnetic field
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
      do i=nxs-2,nxe+2
         uf(1,i,j,k) = 0.0D0
         uf(1,i,j,k) = uf(1,i,j,k) + bx_pert(i*delx,j*delx)
         uf(2,i,j,k) = b_harris(i*delx)
         uf(2,i,j,k) = uf(2,i,j,k) + by_pert(i*delx,j*delx)
         uf(3,i,j,k) = 0.d0
         uf(4,i,j,k) = 0.d0
         uf(5,i,j,k) = 0.d0
         uf(6,i,j,k) = 0.d0
      enddo
    enddo
    enddo
    !$OMP END PARALLEL DO

    ! particle position
    !

    f1 =  1.d0 / ( (1.d0+rtemp) * q(1) )
    f2 = rtemp / ( (1.d0+rtemp) * q(2) )
    ! a factor of 1/sqrt(2) for Box-Muller method
    sdi = vti/sqrt(2.)
    sde = vte/sqrt(2.)
    isp = 1
    ibg = nbg*(nxge-nxgs-2)
    !$OMP PARALLEL DO PRIVATE(ii,j,k,r1)
    do k=nzs,nze
    do j=nys,nye
      ! Background density: Uniform distribution
      do ii=1,ibg
         up(1,ii,j,k,1) = (nxgs+1)*delx + uniform_rand()*delx*(nxge-nxgs-2)
         up(1,ii,j,k,2) = up(1,ii,j,k,1)
         up(2,ii,j,k,1) = dble(j)*delx + uniform_rand()*delx
         up(2,ii,j,k,2) = up(2,ii,j,k,1)
         up(3,ii,j,k,1) = dble(k)*delx + uniform_rand()*delx
         up(3,ii,j,k,2) = up(3,ii,j,k,1)
      enddo
      ! Current sheet: 2nd order logistic distribution in x
      do ii=ibg+1,np2(j,k,1)
         r1 = uniform_rand()
         r1 = (2.0*r1-1.0)*tanh(0.5*(nxge-nxgs-2)*delx/lcs)
         up(1,ii,j,k,1) = lcs * 0.5d0*(log(1.d0+r1)-log(1.d0-r1)) + x0
         up(1,ii,j,k,2) = up(1,ii,j,k,1)
         !r1 = uniform_rand()
         !up(2,ii,j,1) = dble(j)*delx + r1*delx
         !up(2,ii,j,2) = up(2,ii,j,1)
         up(2,ii,j,k,1) = dble(j)*delx + uniform_rand()*delx
         up(2,ii,j,k,2) = up(2,ii,j,k,1)
         up(3,ii,j,k,1) = dble(k)*delx + uniform_rand()*delx
         up(3,ii,j,k,2) = up(3,ii,j,k,1)
      enddo

      ! nonrelativistic Maxwellian 
      do ii=1,np2(j,k,1)
         up(4,ii,j,k,1) = sdi*normal_rand()
         up(5,ii,j,k,1) = sdi*normal_rand()
         up(6,ii,j,k,1) = sdi*normal_rand() + f1*jz(up(1,ii,j,k,1),up(2,ii,j,k,1))/density(up(1,ii,j,k,1))

         up(4,ii,j,k,2) = sde*normal_rand()
         up(5,ii,j,k,2) = sde*normal_rand()
         up(6,ii,j,k,2) = sde*normal_rand() + f2*jz(up(1,ii,j,k,2),up(2,ii,j,k,2))/density(up(1,ii,j,k,2))
      enddo

    enddo
    enddo
              
    !$OMP END PARALLEL DO

    ! set particle IDs
    call set_particle_ids()

  end subroutine set_initial_condition

  !
  ! set initial particle IDs
  !
  subroutine set_particle_ids() ! transported from weibel/app.f90
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

  end subroutine set_particle_ids

  !
  ! output energy history
  !
  subroutine energy_history(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)

    integer :: i, j, k, ii, isp, unit
    real(8) :: vene(nsp)
    real(8) :: efield, bfield, gam, u2
    real(8) :: energy_l(nsp+3), energy_g(nsp+3)

    ! open file
    if( nrank == 0 ) then
      if( it == 0 ) then
        open(newunit=unit,file=trim(datadir)//trim(ehist),status='replace')
      else
        open(newunit=unit,file=trim(datadir)//trim(ehist),status='old',position='append')
      endif
    endif

    ! initialize
    vene(1:nsp) = 0d0
    efield = 0d0
    bfield = 0d0

    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(ii,j,k,u2,gam) REDUCTION(+:vene)
      do k=nzs,nze
      do j=nys,nye
        do ii=1,np2(j,k,isp)
          u2 =  up(4,ii,j,k,isp)*up(4,ii,j,k,isp) &
               +up(5,ii,j,k,isp)*up(5,ii,j,k,isp) &
               +up(6,ii,j,k,isp)*up(6,ii,j,k,isp)
          gam = sqrt(1d0+u2/(c*c))
          vene(isp) = vene(isp)+r(isp)*(gam-1d0)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:bfield,efield)
    do k=nzs,nze
    do j=nys,nye
      do i=nxgs,nxge
        bfield = bfield+uf(1,i,j,k)*uf(1,i,j,k)+uf(2,i,j,k)*uf(2,i,j,k)+uf(3,i,j,k)*uf(3,i,j,k)
        efield = efield+uf(4,i,j,k)*uf(4,i,j,k)+uf(5,i,j,k)*uf(5,i,j,k)+uf(6,i,j,k)*uf(6,i,j,k)
      enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    do isp = 1, nsp
      energy_l(isp) = vene(isp)
    enddo
    energy_l(nsp+1) = efield / (8d0*pi)
    energy_l(nsp+2) = bfield / (8d0*pi)
    call MPI_Reduce(energy_l, energy_g, nsp+2, mnpr, opsum, 0, ncomw, nerr)

    if( nrank == 0 ) then
      ! time, particle1, particle2, efield, bfield, total
      energy_g(5) = sum(energy_g(1:4))
      write(unit, fmt='(f10.2, 5(1x, e12.5))') it*delt, &
            energy_g(1), energy_g(2), energy_g(3), energy_g(4), energy_g(5)
      close(unit)
    endif

  end subroutine energy_history

  !
  ! save everything for restart
  !
  subroutine save_restart(up, uf, np2, nxs, nxe, it, restart_file) ! transported from weibel/app.f90
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
    call io__output(up,uf,np2,nxs,nxe,it,trim(restart_file))

    if ( nrank == 0 ) then
       call json%print(root,config)
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

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine save_param


  !
  ! get global cumulative sum of particle numbers
  !
  subroutine get_global_cumsum(np2,cumsum) ! transported from weibel/app.f90
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


end module app
