module app
  use iso_fortran_env, only: int64
  use mpi
  use wuming2d
  use wuming_utils
  use boundary_periodic, &
       & bc__init        => boundary_periodic__init,       &
       & bc__dfield      => boundary_periodic__dfield,     &
       & bc__particle_x  => boundary_periodic__particle_x, &
       & bc__particle_y  => boundary_periodic__particle_y, &
       & bc__curre       => boundary_periodic__curre,      &
       & bc__phi         => boundary_periodic__phi,        &
       & bc__mom         => boundary_periodic__mom
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
  integer :: n_ppc
  integer :: n_x
  integer :: n_y
  real(8) :: mass_ratio
  real(8) :: sigma_e
  real(8) :: omega_pe
  real(8) :: v_the
  real(8) :: v_thi
  real(8) :: t_ani

  integer :: nproc
  integer :: it0
  integer :: np
  integer :: n0
  integer :: nx, nxgs, nxge, nxs, nxe
  integer :: ny, nygs, nyge !, nys, nye
  integer :: mpierr

  integer, parameter :: ndim   = 6
  integer, parameter :: nsp    = 2
  integer, parameter :: nroot  = 0

  ! OTHER CONSTANTS
  real(8), parameter :: c      = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac   = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl    = 1.0D0     !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: delx   = 1.0D0     !CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

  !
  ! main variables
  !
  integer, allocatable, public :: np2(:,:), cumcnt(:,:,:)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  real(8), allocatable, public :: mom(:,:,:,:)
  real(8)                      :: r(nsp)
  real(8)                      :: q(nsp)
  real(8)                      :: delt
  real(8)                      :: b0

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
       call field__fdtd_i(uf, up, gp, cumcnt, nxs, nxe, &
            & bc__dfield, bc__curre, bc__phi)
       call bc__particle_x(gp, np2)
       call bc__particle_y(gp, np2)
       call sort__bucket(up, gp, cumcnt, np2, nxs, nxe)

       ! output entire particles
       if ( mod(it, intvl_ptcl) == 0 ) then
          call io__ptcl(up, uf, np2, it)
       end if

       ! ouput tracer particles
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
    if( ndim /= 6 ) then
       write(0,*) 'Error: ndim must be 6'
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
    call json%get(p, 'n_ppc', n_ppc)
    call json%get(p, 'n_x', n_x)
    call json%get(p, 'n_y', n_y)
    call json%get(p, 'mass_ratio', mass_ratio)
    call json%get(p, 'sigma_e', sigma_e)
    call json%get(p, 'omega_pe', omega_pe)
    call json%get(p, 'v_the', v_the)
    call json%get(p, 'v_thi', v_thi)
    call json%get(p, 't_ani', t_ani)

    nproc = num_process
    nx    = n_x
    ny    = n_y
    np    = n_ppc * nx * 5
    n0    = n_ppc
    nxgs  = 2
    nxge  = nxgs + nx - 1
    nygs  = 2
    nyge  = nygs + ny - 1
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
    integer :: isp, i, j
    real(8) :: wpe, wpi, wge, wgi, vte, vti

    ! MPI
    call mpi_set__init(nygs, nyge, nproc)

    ! random number
    call init_random_seed()

    ! allocate memory and initialize everything by zero
    allocate(np2(nys:nye,nsp))
    allocate(cumcnt(nxgs:nxge+1,nys:nye,nsp))
    allocate(uf(6,nxgs-2:nxge+2,nys-2:nye+2))
    allocate(up(ndim,np,nys:nye,nsp))
    allocate(gp(ndim,np,nys:nye,nsp))
    allocate(mom(1:7,nxgs-1:nxge+1,nys-1:nye+1,1:nsp))
    np2    = 0
    cumcnt = 0
    uf     = 0
    up     = 0
    gp     = 0
    mom    = 0

    ! set physical parameters
    delt = cfl*delx/c
    wpe  = omega_pe
    wge  = omega_pe * sqrt(sigma_e)
    wpi  = wpe / sqrt(mass_ratio)
    wgi  = wge / mass_ratio
    vte  = v_the
    vti  = v_thi
    r(1) = mass_ratio
    r(2) = 1.0d0
    q(1) =+sqrt(r(1) / (4*pi*n0)) * wpi
    q(2) =-sqrt(r(2) / (4*pi*n0)) * wpe
    b0   = r(1)*c / q(1) * wgi

    ! number of particles
    np2(nys:nye,1:nsp) = n0*(nxge-nxgs+1)
    if ( nrank == nroot ) then
       if ( n0*(nxge-nxgs+1) > np ) then
          write(0,*) 'Error: Too large number of particles'
          stop
       endif
    endif

    ! preparation of sort
    do isp = 1, nsp
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j = nys, nye
          cumcnt(nxgs,j,isp) = 0
          do i = nxgs+1, nxge+1
             cumcnt(i,j,isp) = cumcnt(i-1,j,isp) + n0
          enddo
          if ( cumcnt(nxge+1,j,isp) /= np2(j,isp) ) then
             write(0,*) 'Error: invalid values encounterd for cumcnt'
             stop
          endif
       enddo
       !$OMP END PARALLEL DO
    enddo

    ! initialize modules
    call bc__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & nup, ndown, mnpi, mnpr, ncomw, nerr, nstat, delx, delt, c)
    call particle__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & delx, delt, c, q, r)
    call field__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & mnpr, ncomw, opsum, nerr, delx, delt, c, q, r, gfac)
    call sort__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye)
    call io__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & nproc, nrank, delx, delt, c, q, r, datadir)
    call mom_calc__init( &
         & ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye, &
         & delx, delt, c, q, r)

    if ( restart ) then
       ! restart
       call io__input(gp, uf, np2, nxs, nxe, it0, restart_file)
       call sort__bucket(up, gp, cumcnt, np2, nxs, nxe)
    else
       ! output parameters and set initial condition
       call save_param(n0, wpe, wpi, wge, wgi, vti, vte, param)
       call set_initial_condition()
       it0 = 0
       call energy_history(up, uf, np2, it0)
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
    integer :: i, j, ii, isp
    real(8) :: v1, gam1, gamp, sd(nsp)

    !
    ! electromagnetic field
    !
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=nys-2,nye+2
       do i=nxgs-2,nxge+2
          uf(1,i,j) = 0.0D0
          uf(2,i,j) = 0.0D0
          uf(3,i,j) = b0
          uf(4,i,j) = 0.0D0
          uf(5,i,j) = 0.0D0
          uf(6,i,j) = 0.0D0
       enddo
    enddo
    !$OMP END PARALLEL DO

    !
    ! particle position
    !
    isp = 1
    !$OMP PARALLEL DO PRIVATE(ii,j)
    do j=nys,nye
       do ii=1,np2(j,isp)
          up(1,ii,j,1) = (nxgs + (nxge-nxgs+1)*(ii - 0.5d0)/np2(j,isp)) * delx
          up(2,ii,j,1) = (j + uniform_rand()) * delx
          up(1,ii,j,2) = up(1,ii,j,1)
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !
    ! particle velocity
    !
    sd(1) = v_thi
    sd(2) = v_the
    do isp=1,nsp
       !$OMP PARALLEL DO PRIVATE(ii,j,v1,gam1,gamp)
       do j=nys,nye
          do ii=1,np2(j,isp)
             ! Maxwellian in fluid rest frame
             up(3,ii,j,isp) = sd(isp) * normal_rand()
             up(4,ii,j,isp) = sd(isp) * normal_rand()
             up(5,ii,j,isp) = t_ani * sd(isp) * normal_rand()
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
    integer :: isp, i, j

    integer(8) :: gcumsum(nproc+1,nsp), lcumsum(nys:nye+1,nsp), pid

    if( ndim /= 6 ) then
       return
    end if

    ! calculate the first particle IDs
    call get_global_cumsum(np2, gcumsum)

    do isp = 1, nsp
       lcumsum(nys,isp) = gcumsum(nrank+1,isp)
       do j = nys, nye
          lcumsum(j+1,isp) = lcumsum(j,isp) + np2(j,isp)
       end do
    end do

    ! unique ID as 64bit integer (negative by default)
    do isp = 1, nsp
       !$OMP PARALLEL DO PRIVATE(i,j,pid)
       do j = nys, nye
          do i = 1, np2(j,isp)
             pid = lcumsum(j,isp) + i
             up(6,i,j,isp) = transfer(-pid, 1.0_8)
          end do
       end do
       !$OMP END PARALLEL DO
    end do

  end subroutine set_particle_ids

  !
  ! output energy history
  !
  subroutine energy_history(up, uf, np2, it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    integer :: i, j, ii, isp, unit
    real(8) :: vene(nsp)
    real(8) :: efield, bfield, gam, u2
    real(8) :: energy_l(nsp+3), energy_g(nsp+3)

    ! open file
    if( nrank == 0 ) then
       if( it == 0 ) then
          open(newunit=unit, file=trim(datadir) // trim(ehist), &
               & status='replace')
       else
          open(newunit=unit, file=trim(datadir) // trim(ehist), &
               & status='old', position='append')
       end if
    endif

    ! initialize
    vene(1:nsp) = 0
    efield = 0
    bfield = 0

    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(ii,j,u2,gam) REDUCTION(+:vene)
       do j=nys,nye
       do ii=1,np2(j,isp)
          u2 =  up(3,ii,j,isp)*up(3,ii,j,isp) &
               +up(4,ii,j,isp)*up(4,ii,j,isp) &
               +up(5,ii,j,isp)*up(5,ii,j,isp)
          gam = sqrt(1+u2/(c*c))
          vene(isp) = vene(isp)+r(isp)*(gam-1)
       enddo
       enddo
!$OMP END PARALLEL DO
    enddo

!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:bfield,efield)
    do j=nys,nye
    do i=nxgs,nxge
       bfield = bfield+uf(1,i,j)*uf(1,i,j)+uf(2,i,j)*uf(2,i,j)+uf(3,i,j)*uf(3,i,j)
       efield = efield+uf(4,i,j)*uf(4,i,j)+uf(5,i,j)*uf(5,i,j)+uf(6,i,j)*uf(6,i,j)
    enddo
    enddo
!$OMP END PARALLEL DO

    do isp = 1, nsp
       energy_l(isp) = vene(isp)
    end do
    energy_l(nsp+1) = efield / (8*pi)
    energy_l(nsp+2) = bfield / (8*pi)
    call MPI_Reduce(energy_l, energy_g, nsp+2, mnpr, opsum, 0, ncomw, nerr)

    if( nrank == 0 ) then
       ! time, particle1, particle2, efield, bfield, total
       energy_g(nsp+3) = sum(energy_g(1:nsp+2))
       write(unit, fmt='(f8.2, 4(1x, e12.5))') it*delt, &
            & sum(energy_g(1:nsp)), energy_g(nsp+1), energy_g(nsp+2), energy_g(nsp+3)
       close(unit)
    endif

  end subroutine energy_history

  !
  ! save everything for restart
  !
  subroutine save_restart(up, uf, np2, nxs, nxe, it, restart_file)
    implicit none
    integer, intent(in)          :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in)          :: it
    real(8), intent(in)          :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in)          :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
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

    if ( nrank == nroot ) then
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
  subroutine get_global_cumsum(np2, cumsum)
    implicit none
    integer, intent(in)       :: np2(nys:nye,nsp)
    integer(8), intent(inout) :: cumsum(nproc+1,nsp)

    integer :: i, isp, mpierr
    integer(8) :: lcount(nsp), gcount(nsp, nproc)

    ! get number of particles for each proces
    lcount(1:nsp) = sum(np2(nys:nye,1:nsp), dim=1)
    call MPI_Allgather(lcount, nsp, MPI_INTEGER8, gcount, nsp, MPI_INTEGER8, &
         & MPI_COMM_WORLD, mpierr)

    ! calculate cumulative sum
    do isp = 1, nsp
       cumsum(1,isp) = 0
       do i = 1, nproc
          cumsum(i+1,isp) = cumsum(i,isp) + gcount(isp,i)
       end do
    end do

  end subroutine get_global_cumsum


end module app
