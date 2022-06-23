module paraio
  use iso_fortran_env, only: int64
  use mpi
  use jsonio
  use mpiio
  implicit none
  private

  public :: paraio__init
  public :: paraio__finalize
  public :: paraio__output
  public :: paraio__input
  public :: paraio__param
  public :: paraio__mom
  public :: paraio__ptcl
  public :: paraio__orb


  integer, parameter :: MOK = MPI_OFFSET_KIND ! assume to be the same as int64

  character(len=256), save :: dir
  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nproc, nrank
  real(8), save :: delx, delt, u0, c
  real(8), allocatable :: q(:), r(:)

  integer :: mpierr
  real(8), allocatable :: mpibuf1(:)
  integer, allocatable :: mpibuf2(:)

contains

  !
  ! initialize module
  !
  subroutine paraio__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
       & nproc_in,nrank_in,delx_in,delt_in,c_in,q_in,r_in,dir_in)
    implicit none
    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nproc_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

    integer :: psize, fsize, isize

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in
    nproc = nproc_in
    nrank = nrank_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in
    dir   = dir_in

    ! allocate MPI buffer
    psize = ndim*np*(nye-nys+5)*nsp
    fsize = 6*(nxge-nxgs+5)*(nye-nys+5)
    isize = (nye-nys+1)*nsp
    allocate(mpibuf1(max(psize, fsize)))
    allocate(mpibuf2(isize))

    is_init = .true.

  end subroutine paraio__init

  !
  ! finalize module
  !
  subroutine paraio__finalize()
    implicit none

    deallocate(mpibuf1)
    deallocate(mpibuf2)

  end subroutine paraio__finalize

  !
  ! output data for re-calculation
  !
  subroutine paraio__output(up,uf,np2,nxs,nxe,it,filename)
    implicit none
    integer, intent(in) :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in) :: it
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    character(len=*), intent(in) :: filename

    character(len=256) :: jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize, psize, goffset
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, gshape(5)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    ! filename
    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'w')
#ifdef _MPIIO_OPEN_CLOSE
    call mpiio_close_file(fh)
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'a')
#endif
    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    !
    ! attribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, MOK, 'dummy_attribute', disp, '')
    call mpiio_write_atomic(fh, disp, MOK)

    call jsonio_put_attribute(json, p, it, 'it', disp, '')
    call mpiio_write_atomic(fh, disp, it)

    call jsonio_put_attribute(json, p, nxs, 'nxs', disp, '')
    call mpiio_write_atomic(fh, disp, nxs)

    call jsonio_put_attribute(json, p, nxe, 'nxe', disp, '')
    call mpiio_write_atomic(fh, disp, nxe)

    call put_metadata(json, p, fh, disp)

    !
    ! dataset
    !

    call json%create_object(p, 'dataset')
    call json%add(root, p)

    ! particle
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    nd      = 5
    gshape  = (/ndim, np, nyl, nsp, nproc/)
    lsize   = size(up, kind=8)
    psize   = ndim
    gsize   = lsize * nproc
    goffset = lsize * nrank
    dsize   = gsize * 8
    desc    = 'particles'
    mpibuf1(1:lsize) = reshape(up, (/lsize/))
    call jsonio_put_metadata(json, p, 'up', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
         & 8, mpibuf1)

    nd      = 3
    gshape  = (/nyl, nsp, nproc, 0, 0/)
    lsize   = size(np2, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    dsize   = gsize * 4
    desc    = 'number of active particles'
    mpibuf2(1:lsize) = reshape(np2, (/lsize/))
    call jsonio_put_metadata(json, p, 'np2', 'i4', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
         & 4, mpibuf2)

    ! field: including ghost cells
    nxg = size(uf, 2) ! nxge - nxgs + 5
    nyl = size(uf, 3) ! nye  - nys  + 5
    nyg = nyl * nproc

    nd      = 4
    gshape  = (/6, nxg, nyl, nproc, 0/)
    lsize   = size(uf, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    dsize   = gsize * 8
    desc    = 'electromagnetic fields including ghost cells'
    mpibuf1(1:lsize) = reshape(uf, (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
         & 8, mpibuf1)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, trim(dir) // jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__output

  !
  ! input data for re-calculation
  !
  subroutine paraio__input(up,uf,np2,nxs,nxe,it,filename)
    implicit none
    real(8), intent(out)         :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out)         :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer, intent(out)         :: np2(nys:nye,nsp)
    integer, intent(out)         :: nxs
    integer, intent(out)         :: nxe
    integer, intent(out)         :: it
    character(len=*), intent(in) :: filename

    integer :: inp, indim, inxgs, inxge, inygs, inyge, insp, inproc

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp, dsize, lsize, psize, gsize, goffset
    integer :: fh, nd, gshape(5)

    type(json_core) :: json
    type(json_file) :: file
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'r')

    ! json file
    call json%initialize()
    call file%initialize()
    call file%load(trim(dir) // jsonfile)
    call file%get(root)

    !
    ! atttribute
    !
    call json%get(root, 'attribute', p)

    call jsonio_get_metadata(json, p, 'it', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, it)

    call jsonio_get_metadata(json, p, 'nxs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxs)

    call jsonio_get_metadata(json, p, 'nxe', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, nxe)

    call jsonio_get_metadata(json, p, 'ndim', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, indim)

    call jsonio_get_metadata(json, p, 'np', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inp)

    call jsonio_get_metadata(json, p, 'nxgs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inxgs)

    call jsonio_get_metadata(json, p, 'nxge', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inxge)

    call jsonio_get_metadata(json, p, 'nygs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inygs)

    call jsonio_get_metadata(json, p, 'nyge', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inyge)

    call jsonio_get_metadata(json, p, 'nsp', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, insp)

    call jsonio_get_metadata(json, p, 'nproc', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inproc)

    call jsonio_get_metadata(json, p, 'delx', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, delx)

    call jsonio_get_metadata(json, p, 'delt', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, delt)

    call jsonio_get_metadata(json, p, 'c', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, c)

    call jsonio_get_metadata(json, p, 'r', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, r)

    call jsonio_get_metadata(json, p, 'q', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, q)

    if( .not. ( &
         & inxgs == nxgs .and. &
         & inxge == nxge .and. &
         & inygs == nygs .and. &
         & inyge == nyge .and. &
         & inp == np .and. &
         & insp == nsp .and. &
         & inproc == nproc) ) then
       write(0,*) '*** Error: invalid input parameters ***'
       stop
    endif

    !
    ! dataset
    !
    call json%get(root, 'dataset', p)

    ! * particle

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'up', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 5                  .and. &
         &   gshape(1) == size(up, 1) .and. &
         &   gshape(2) == size(up, 2) .and. &
         &   gshape(3) == size(up, 3) .and. &
         &   gshape(4) == size(up, 4) .and. &
         &   gshape(5) == nproc ) ) then
       write(0, *) 'Fatal error in reading particle data'
       call MPI_Finalize(mpierr)
       stop
    end if

    ! read data
    lsize   = size(up, kind=8)
    psize   = ndim
    gsize   = lsize * nproc
    goffset = lsize * nrank
    call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
         & 8, mpibuf1)
    up = reshape(mpibuf1(1:lsize), shape(up))

    ! * particle number

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'np2', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 3                  .and. &
         &   gshape(1) == size(np2, 1) .and. &
         &   gshape(2) == size(np2, 2) .and. &
         &   gshape(3) == nproc ) ) then
       write(0, *) 'Fatal error in reading particle data'
       call MPI_Finalize(mpierr)
       stop
    end if

    ! read data
    lsize   = size(np2, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
         & 4, mpibuf2)
    np2 = reshape(mpibuf2(1:lsize), shape(np2))

    ! * field

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'uf', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 4                  .and. &
         &   gshape(1) == size(uf, 1) .and. &
         &   gshape(2) == size(uf, 2) .and. &
         &   gshape(3) == size(uf, 3) .and. &
         &   gshape(4) == nproc ) ) then
       write(0, *) 'Fatal error in reading particle data'
       call MPI_Finalize(mpierr)
       stop
    end if

    ! read data
    lsize   = size(uf, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
         & 8, mpibuf1)
    uf = reshape(mpibuf1(1:lsize), shape(uf))

    !
    ! finalize
    !

    ! close json
    call json%destroy()
    call file%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__input

  !
  ! output parameters
  !
  subroutine paraio__param(n0, wpe, wpi, wge, wgi, vti, vte, filename)
    implicit none
    integer, intent(in)          :: n0
    real(8), intent(in)          :: wpe, wpi, wge, wgi, vti, vte
    character(len=*), intent(in) :: filename

    character(len=256) :: jsonfile, datafile
    integer(int64) :: disp
    integer :: endian, fh, nx, ny
    real(8) :: vai, vae

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    vai = c * wgi/wpi
    vae = c * wge/wpe

    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    ! open data file
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'w')
#ifdef _MPIIO_OPEN_CLOSE
    call mpiio_close_file(fh)
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'a')
#endif
    !
    ! attribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, MOK, 'dummy_attribute', disp, '')
    call mpiio_write_atomic(fh, disp, MOK)

    nx = nxge - nxgs + 1
    call jsonio_put_attribute(json, p, nx, 'nx', disp, '')
    call mpiio_write_atomic(fh, disp, nx)

    ny = nyge - nygs + 1
    call jsonio_put_attribute(json, p, ny, 'ny', disp, '')
    call mpiio_write_atomic(fh, disp, ny)

    call jsonio_put_attribute(json, p, np, 'np', disp, '')
    call mpiio_write_atomic(fh, disp, np)

    call jsonio_put_attribute(json, p, delx, 'delx', disp, '')
    call mpiio_write_atomic(fh, disp, delx)

    call jsonio_put_attribute(json, p, delt, 'delt', disp, '')
    call mpiio_write_atomic(fh, disp, delt)

    call jsonio_put_attribute(json, p, c, 'c', disp, '')
    call mpiio_write_atomic(fh, disp, c)

    call jsonio_put_attribute(json, p, r, 'r', disp, '')
    call mpiio_write_atomic(fh, disp, r)

    call jsonio_put_attribute(json, p, q, 'q', disp, '')
    call mpiio_write_atomic(fh, disp, q)

    call jsonio_put_attribute(json, p, wpe, 'wpe', disp, '')
    call mpiio_write_atomic(fh, disp, wpe)

    call jsonio_put_attribute(json, p, wge, 'wge', disp, '')
    call mpiio_write_atomic(fh, disp, wge)

    call jsonio_put_attribute(json, p, wpi, 'wpi', disp, '')
    call mpiio_write_atomic(fh, disp, wpi)

    call jsonio_put_attribute(json, p, wgi, 'wgi', disp, '')
    call mpiio_write_atomic(fh, disp, wgi)

    call jsonio_put_attribute(json, p, vti, 'vti', disp, '')
    call mpiio_write_atomic(fh, disp, vti)

    call jsonio_put_attribute(json, p, vte, 'vte', disp, '')
    call mpiio_write_atomic(fh, disp, vte)

    call jsonio_put_attribute(json, p, vai, 'vai', disp, '')
    call mpiio_write_atomic(fh, disp, vai)

    call jsonio_put_attribute(json, p, vae, 'vae', disp, '')
    call mpiio_write_atomic(fh, disp, vae)

    call jsonio_put_attribute(json, p, n0, 'n0', disp, '')
    call mpiio_write_atomic(fh, disp, n0)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, trim(dir) // jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__param

  !
  ! output field and moment quantities
  !
  subroutine paraio__mom(mom,uf,it)
    implicit none
    integer, intent(in)    :: it
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: mom(7,nxgs-1:nxge+1,nys-1:nye+1,nsp)

    character(len=256) :: filename, jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize
    integer :: i, j
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, lshape(4), gshape(4), offset(4)
    real(8) :: tmp(1:6,nxgs:nxge,nys:nye)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    write(filename,'(i7.7, "_mom")') it
    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'w')
#ifdef _MPIIO_OPEN_CLOSE
    call mpiio_close_file(fh)
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'a')
#endif

    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    !
    ! atttribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, MOK, 'dummy_attribute', disp, '')
    call mpiio_write_atomic(fh, disp, MOK)

    call jsonio_put_attribute(json, p, it, 'it', disp, '')
    call mpiio_write_atomic(fh, disp, it)

    call put_metadata(json, p, fh, disp)

    !
    ! dataset
    !
    call json%create_object(p, 'dataset')
    call json%add(root, p)

    ! density
    nxg  = nxge - nxgs + 1
    nyg  = nyge - nygs + 1
    nyl  = nye  - nys  + 1

    ! density
    nd    = 3
    lshape = (/nxg, nyl, nsp, 0/)
    gshape = (/nxg, nyg, nsp, 0/)
    offset = (/0, nyl*nrank, 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'density'
    mpibuf1(1:lsize) = reshape(mom(1:1,nxgs:nxge,nys:nye,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'den', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! velocity
    mom(2,:,:,:) = mom(2,:,:,:) / mom(1,:,:,:)
    mom(3,:,:,:) = mom(3,:,:,:) / mom(1,:,:,:)
    mom(4,:,:,:) = mom(4,:,:,:) / mom(1,:,:,:)

    nd     = 4
    lshape = (/3, nxg, nyl, nsp/)
    gshape = (/3, nxg, nyg, nsp/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'velocity'
    mpibuf1(1:lsize) = reshape(mom(2:4,nxgs:nxge,nys:nye,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'vel', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! temperature
    mom(5,:,:,:) = mom(5,:,:,:) / mom(1,:,:,:)
    mom(6,:,:,:) = mom(6,:,:,:) / mom(1,:,:,:)
    mom(7,:,:,:) = mom(7,:,:,:) / mom(1,:,:,:)

    nd    = 4
    lshape = (/3, nxg, nyl, nsp/)
    gshape = (/3, nxg, nyg, nsp/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'temperature'
    mpibuf1(1:lsize) = reshape(mom(5:7,nxgs:nxge,nys:nye,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'temp', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! electromagnetic field on cell centers
    do j=nys,nye
       do i=nxgs,nxge
          tmp(1,i,j) = (uf(1,i,j)+uf(1,i,j+1)) / 2
          tmp(2,i,j) = (uf(2,i,j)+uf(2,i+1,j)) / 2
          tmp(3,i,j) = (uf(3,i,j)+uf(3,i+1,j)+uf(3,i+1,j)+uf(3,i+1,j+1)) / 4
          tmp(4,i,j) = (uf(4,i,j)+uf(4,i+1,j)) / 2
          tmp(5,i,j) = (uf(5,i,j)+uf(5,i,j+1)) / 2
          tmp(6,i,j) = uf(6,i,j)
       enddo
    enddo

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    mpibuf1(1:lsize) = reshape(tmp(1:6,nxgs:nxge,nys:nye), (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, trim(dir) // jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__mom

  !
  ! output all the active particles
  !
  subroutine paraio__ptcl(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    call write_particle(up, uf, np2, it, 0, '_ptcl')

  end subroutine paraio__ptcl

  !
  ! output tracer particles
  !
  subroutine paraio__orb(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    call write_particle(up, uf, np2, it, 1, '_orb')

  end subroutine paraio__orb

  !
  ! output particles
  !
  subroutine write_particle(up,uf,np2,it,mode,suffix)
    implicit none
    integer, intent(in)          :: it
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in)          :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer, intent(in)          :: mode
    character(len=*), intent(in) :: suffix

    character(len=256) :: filename, jsonfile, datafile, desc, name
    integer(int64) :: disp, dsize, lsize, psize, gsize, goffset
    integer(int64) :: cumsum(nproc+1,nsp), ip1, ip2
    integer(int64) :: npl, npg, npo, nd8, gshape8(2)
    integer :: isp, irank
    integer :: fh, endian, nxg, nyg, nyl
    integer :: nd, lshape(4), gshape(4), offset(4)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
       write(0,*) 'Initialize first by calling paraio__init()'
       stop
    endif

    ! filename
    write(filename,'(i7.7, a)') it, trim(suffix)
    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    ! open data file
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'w')
#ifdef _MPIIO_OPEN_CLOSE
    call mpiio_close_file(fh)
    call mpiio_open_file(trim(dir) // datafile, fh, disp, 'a')
#endif
    ! open json file
    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! metadata
    !
    endian = mpiio_get_endian_flag()
    call json%create_object(p, 'meta')
    call json%add(root, p)
    call json%add(p, 'endian', endian)
    call json%add(p, 'rawfile', trim(datafile))

    !
    ! atttribute
    !
    call json%create_object(p, 'attribute')
    call json%add(root, p)

    call jsonio_put_attribute(json, p, MOK, 'dummy_attribute', disp, '')
    call mpiio_write_atomic(fh, disp, MOK)

    call jsonio_put_attribute(json, p, it, 'it', disp, '')
    call mpiio_write_atomic(fh, disp, it)

    call put_metadata(json, p, fh, disp)

    !
    ! dataset
    !
    call json%create_object(p, 'dataset')
    call json%add(root, p)

    ! field
    nxg = nxge - nxgs + 1
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1

    nd     = 3
    lshape = (/6, nxg, nyl, 0/)
    gshape = (/6, nxg, nyg, 0/)
    offset = (/0, 0, nyl*nrank, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    mpibuf1(1:lsize) = reshape(uf(1:6,nxgs:nxge,nys:nye), (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
         & dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! particle
    call get_particle_count(up, np2, mpibuf1, cumsum, mode)

    irank = nrank + 1
    ip2   = 0
    do isp = 1, nsp
       write(desc, '("particle species #", i2.2)') isp
       write(name, '("up", i2.2)') isp

       npl     = cumsum(irank+1,isp) - cumsum(irank,isp)
       npg     = cumsum(nproc+1,isp)
       npo     = cumsum(irank,isp)
       ip1     = ip2 + 1
       ip2     = ip1 + ndim*npl - 1
       lsize   = ndim * npl
       psize   = ndim
       gsize   = ndim * npg
       goffset = ndim * npo
       dsize   = gsize * 8

       nd8 = 2
       gshape8(1) = ndim
       gshape8(2) = npg
       call jsonio_put_metadata(json, p, trim(name), 'f8', disp, &
            & dsize, nd8, gshape8, desc)
       call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
            & 8, transfer(mpibuf1(ip1:ip2), (/1/)))
    end do

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
       call json%print(root, trim(dir) // jsonfile)
    end if
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine write_particle

  !
  ! get particles distribution and pack to buffer
  !
  subroutine get_particle_count(up, np2, buf, cumsum, mode)
    implicit none
    integer, intent(in)       :: np2(nys:nye,nsp)
    real(8), intent(in)       :: up(ndim,np,nys:nye,nsp)
    real(8), intent(inout)    :: buf(:)
    integer(8), intent(inout) :: cumsum(nproc+1,nsp)
    integer, intent(in)       :: mode

    integer :: i, j, ip, jp, isp
    integer(8) :: lcount(nsp), gcount(nsp, nproc)
    integer(8) :: pid

    ! count number of particles and pack into buffer
    if ( mode == 0 ) then
       ! * mode 0: all the active particles
       ip = 1
       do isp = 1, nsp
          lcount(isp) = 0
          do j = nys, nye
             do i = 1, np2(j,isp)
                lcount(isp) = lcount(isp) + 1
                ! packing
                do jp = 1, ndim
                   buf(ip) = up(jp,i,j,isp)
                   ip = ip + 1
                end do
             enddo
          end do
       end do

    else if ( mode == 1 ) then
       ! * mode 1: tracer particles with positive IDs
       ip = 1
       do isp = 1, nsp
          lcount(isp) = 0
          do j = nys, nye
             do i = 1, np2(j,isp)
                ! get particle ID as 64bit integer
                pid = transfer(up(ndim,i,j,isp), 1_8)

                ! count positive
                if( pid > 0 ) then
                   lcount(isp) = lcount(isp) + 1
                   ! packing
                   do jp = 1, ndim
                      buf(ip) = up(jp,i,j,isp)
                      ip = ip + 1
                   end do
                endif
             enddo
          end do
       end do

    else
       ! error
       write(0,*) 'Error: invalid mode specified for get_particle_count'
       call MPI_Finalize(mpierr)
       stop
    end if

    call MPI_Allgather(lcount, nsp, MPI_INTEGER8, gcount, nsp, MPI_INTEGER8, &
         & MPI_COMM_WORLD, mpierr)

    ! calculate cumulative sum
    do isp = 1, nsp
       cumsum(1,isp) = 0
       do i = 1, nproc
          cumsum(i+1,isp) = cumsum(i,isp) + gcount(isp,i)
       end do
    end do

  end subroutine get_particle_count

  !
  ! put common metadata
  !
  subroutine put_metadata(json, p, file, disp)
    implicit none
    type(json_core), intent(inout)        :: json
    type(json_value), pointer, intent(in) :: p
    integer, intent(in)                   :: file
    integer(int64), intent(inout)         :: disp

    call jsonio_put_attribute(json, p, ndim, 'ndim', disp, '')
    call mpiio_write_atomic(file, disp, ndim)

    call jsonio_put_attribute(json, p, np, 'np', disp, '')
    call mpiio_write_atomic(file, disp, np)

    call jsonio_put_attribute(json, p, nxgs, 'nxgs', disp, '')
    call mpiio_write_atomic(file, disp, nxgs)

    call jsonio_put_attribute(json, p, nxge, 'nxge', disp, '')
    call mpiio_write_atomic(file, disp, nxge)

    call jsonio_put_attribute(json, p, nygs, 'nygs', disp, '')
    call mpiio_write_atomic(file, disp, nygs)

    call jsonio_put_attribute(json, p, nyge, 'nyge', disp, '')
    call mpiio_write_atomic(file, disp, nyge)

    call jsonio_put_attribute(json, p, nsp, 'nsp', disp, '')
    call mpiio_write_atomic(file, disp, nsp)

    call jsonio_put_attribute(json, p, nproc, 'nproc', disp, '')
    call mpiio_write_atomic(file, disp, nproc)

    call jsonio_put_attribute(json, p, delx, 'delx', disp, '')
    call mpiio_write_atomic(file, disp, delx)

    call jsonio_put_attribute(json, p, delt, 'delt', disp, '')
    call mpiio_write_atomic(file, disp, delt)

    call jsonio_put_attribute(json, p, c, 'c', disp, '')
    call mpiio_write_atomic(file, disp, c)

    call jsonio_put_attribute(json, p, r, 'r', disp, '')
    call mpiio_write_atomic(file, disp, r)

    call jsonio_put_attribute(json, p, q, 'q', disp, '')
    call mpiio_write_atomic(file, disp, q)

  end subroutine put_metadata

end module paraio
