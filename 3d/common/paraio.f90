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
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
  integer, save :: nproc, nproc_j, nproc_k, nrank
  real(8), save :: delx, delt, c
  real(8), allocatable :: q(:), r(:)

  integer :: mpierr
  real(8), allocatable :: mpibuf1(:)
  integer, allocatable :: mpibuf2(:)
  integer(int64), allocatable :: mpibuf3(:)

contains

  !
  ! initialize module
  !
  subroutine paraio__init(ndim_in,np_in,nsp_in,                                                        &
                          nxgs_in,nxge_in,nygs_in,nyge_in,nzgs_in,nzge_in,nys_in,nye_in,nzs_in,nze_in, &
                          nproc_in,nproc_j_in,nproc_k_in,nrank_in,delx_in,delt_in,c_in,q_in,r_in,dir_in)
    implicit none
    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nzgs_in, nzge_in, nys_in, nye_in, nzs_in, nze_in
    integer, intent(in) :: nproc_in, nproc_j_in, nproc_k_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

    integer :: psize, fsize, isize

    ndim    = ndim_in
    np      = np_in
    nsp     = nsp_in
    nxgs    = nxgs_in
    nxge    = nxge_in
    nygs    = nygs_in
    nyge    = nyge_in
    nzgs    = nzgs_in
    nzge    = nzge_in
    nys     = nys_in
    nye     = nye_in
    nzs     = nzs_in
    nze     = nze_in
    nproc   = nproc_in
    nproc_j = nproc_j_in
    nproc_k = nproc_k_in
    nrank   = nrank_in
    delx    = delx_in
    delt    = delt_in
    c       = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q       = q_in
    r       = r_in
    dir     = dir_in

    ! allocate MPI buffer
    psize = ndim*np*(nye-nys+5)*(nze-nzs+5)*nsp
    fsize = 6*(nxge-nxgs+5)*(nye-nys+5)*(nze-nzs+5)
    isize = (nye-nys+1)*(nze-nzs+1)*nsp
    allocate(mpibuf1(max(psize, fsize)))
    allocate(mpibuf2(isize))
    allocate(mpibuf3(nsp))

    is_init = .true.

  end subroutine paraio__init

  !
  ! finalize module
  !
  subroutine paraio__finalize()
    implicit none

    deallocate(mpibuf1)
    deallocate(mpibuf2)
    deallocate(mpibuf3)

  end subroutine paraio__finalize

  !
  ! output data for re-calculation
  !
  subroutine paraio__output(up,uf,np2,nxs,nxe,it,filename)
    implicit none
    integer, intent(in) :: np2(nys:nye,nzs:nze,nsp), nxs, nxe
    integer, intent(in) :: it
    real(8), intent(in) :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    character(len=*), intent(in) :: filename

    character(len=256) :: jsonfile, datafile, desc, name !, sharedtmp
    integer(int64) :: disp, dsize, lsize, gsize, psize, goffset
    integer(int64) :: cumsum(nproc+1,nsp), ip1, ip2, poffset(nsp)
    integer(int64) :: npl, npg, npo, nd8, gshape8(2)
    integer :: isp, irank
    integer :: fh, endian, nxg, nyg, nyl, nzg, nzl
    integer :: nd, gshape(6)

    type(json_core) :: json
    type(json_value), pointer :: root, p

    if ( .not. is_init ) then
      write(0,*) 'Initialize first by calling paraio__init()'
      stop
    endif

    ! filename
    datafile = trim(filename) // '.raw'
    jsonfile = trim(filename) // '.json'

    !
    ! Using shared temporary area in Fugaku
    !

    !call GET_ENVIRONMENT_VARIABLE("PJM_SHAREDTMP",sharedtmp)
    !call mpiio_open_file(trim(sharedtmp) // "/" // datafile, fh, disp, 'w')

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

    ! particle number
    nyg = nyge - nygs + 1
    nyl = nye  - nys  + 1
    nzg = nzge - nzgs + 1
    nzl = nze  - nzs  + 1

    nd      = 4
    gshape  = (/nyl, nzl, nsp, nproc, 0, 0/)
    lsize   = size(np2, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    dsize   = gsize * 4
    desc    = 'number of active particles'
    mpibuf2(1:lsize) = reshape(np2, (/lsize/))
    call jsonio_put_metadata(json, p, 'np2', 'i4', disp, &
                             dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
                                4, mpibuf2)

    ! particle
    call get_particle_count(up, np2, mpibuf1, cumsum, 0)

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
                               dsize, nd8, gshape8, desc)
      call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
                                  8, transfer(mpibuf1(ip1:ip2), (/1/)))

      poffset(isp) = goffset
    enddo

    ! particle offset
    nd      = 2
    gshape  = (/nsp, nproc, 0, 0, 0, 0/)
    lsize   = size(poffset, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    dsize   = gsize * 8
    desc    = 'particle offset'
    mpibuf3 = reshape(poffset, (/lsize/))
    call jsonio_put_metadata(json, p, 'poffset', 'i8', disp, &
                             dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
                                8, mpibuf3)

    ! field: including ghost cells
    nxg = size(uf, 2) ! nxge - nxgs + 5
    nyl = size(uf, 3) ! nye  - nys  + 5
    nyg = nyl * nproc_j
    nzl = size(uf, 4)
    nzg = nzl * nproc_k

    nd      = 5
    gshape  = (/6, nxg, nyl, nzl, nproc, 0/)
    lsize   = size(uf, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    dsize   = gsize * 8
    desc    = 'electromagnetic fields including ghost cells'
    mpibuf1(1:lsize) = reshape(uf, (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
                             dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
                                8, mpibuf1)

    !
    ! finalize
    !

    ! write json and close
    if( nrank == 0 ) then
      call json%print(root, trim(dir) // jsonfile)
    endif
    call json%destroy()

    ! close data file
    call mpiio_close_file(fh)

  end subroutine paraio__output

  !
  ! input data for re-calculation
  !
  subroutine paraio__input(up,uf,np2,nxs,nxe,it,filename)
    implicit none
    real(8), intent(out)         :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(out)         :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    integer, intent(out)         :: np2(nys:nye,nzs:nze,nsp)
    integer, intent(out)         :: nxs
    integer, intent(out)         :: nxe
    integer, intent(out)         :: it
    character(len=*), intent(in) :: filename

    integer :: inp, indim, inxgs, inxge, inygs, inyge, inzgs, inzge, insp, inproc
    integer :: i, j, k, isp, ip, jp

    character(len=256) :: jsonfile, datafile, name
    integer(int64) :: disp, dsize, lsize, psize, gsize, goffset, poffset(nsp), nd8, gshape8(2)
    integer :: fh, nd, gshape(6)

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

    call jsonio_get_metadata(json, p, 'nzgs', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inzgs)

    call jsonio_get_metadata(json, p, 'nzge', disp, dsize, nd, gshape)
    call mpiio_read_atomic(fh, disp, inzge)

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
         & inzgs == nzgs .and. &
         & inzge == nzge .and. &
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

    ! * particle number

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'np2', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 4                   .and. &
         &   gshape(1) == size(np2, 1) .and. &
         &   gshape(2) == size(np2, 2) .and. &
         &   gshape(3) == size(np2, 3) .and. &
         &   gshape(4) == nproc ) ) then
      write(0, *) 'Fatal error in reading particle number'
      call MPI_Finalize(mpierr)
      stop
    endif

    ! read data
    lsize   = size(np2, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
                               4, mpibuf2)
    np2 = reshape(mpibuf2(1:lsize), shape(np2))

    ! * particle offset

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'poffset', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 2                   .and. &
         &   gshape(1) == size(poffset) .and. &
         &   gshape(2) == nproc ) ) then
      write(0, *) 'Fatal error in reading particle offset'
      call MPI_Finalize(mpierr)
      stop
    endif

    ! read data
    lsize   = size(poffset, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
                               8, mpibuf3)
    poffset = reshape(mpibuf3(1:lsize), shape(poffset))

    ! * particle

    do isp=1,nsp
      ! read metadata and check
      write(name, '("up", i2.2)') isp
      call jsonio_get_metadata(json, p, trim(name), disp, dsize, nd8, gshape8)

      if ( .not. &
           & ( nd8 == 2           .and. &
           &   gshape8(1) == ndim ) ) then
        write(0, *) 'Fatal errorin reading particle data'
        call MPI_Finalize(mpierr)
        stop
      endif

      ! read data
      lsize   = ndim * sum(int(np2(:,:,isp), kind=8))
      psize   = 1
      gsize   = lsize * nproc
      goffset = poffset(isp)
      call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
                                 8, mpibuf1)


      ip = 1
      do k = nzs, nze
      do j = nys, nye
        do i = 1, np2(j,k,isp)
          ! unpacking
          do jp = 1, ndim
            up(jp,i,j,k,isp) = mpibuf1(ip)
            ip = ip + 1
          enddo
        enddo
      enddo
      enddo
    enddo

    ! * field

    ! read metadata and check
    call jsonio_get_metadata(json, p, 'uf', disp, dsize, nd, gshape)

    if ( .not. &
         & ( nd == 5                  .and. &
         &   gshape(1) == size(uf, 1) .and. &
         &   gshape(2) == size(uf, 2) .and. &
         &   gshape(3) == size(uf, 3) .and. &
         &   gshape(4) == size(uf, 4) .and. &
         &   gshape(5) == nproc ) ) then
      write(0, *) 'Fatal error in reading field data'
      call MPI_Finalize(mpierr)
      stop
    endif

    ! read data
    lsize   = size(uf, kind=8)
    psize   = 1
    gsize   = lsize * nproc
    goffset = lsize * nrank
    call mpiio_read_collective(fh, disp, gsize, lsize, psize, goffset, &
                               8, mpibuf1)
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
    integer :: endian, fh, nx, ny, nz
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

    nz = nzge - nzgs + 1
    call jsonio_put_attribute(json, p, nz, 'nz', disp, '')
    call mpiio_write_atomic(fh, disp, nz)

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
    endif
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
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(inout) :: mom(7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,nsp)

    character(len=256) :: filename, jsonfile, datafile, desc
    integer(int64) :: disp, dsize, lsize, gsize
    integer :: i, j, k
    integer :: fh, endian, nxg, nyg, nyl, nzg, nzl
    integer :: nd, lshape(5), gshape(5), offset(5)
    real(8) :: tmp(1:6,nxgs:nxge,nys:nye,nzs:nze)

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
    nzg  = nzge - nzgs + 1
    nzl  = nze  - nzs  + 1

    ! density
    nd    = 4
    lshape = (/nxg, nyl, nzl, nsp, 0/)
    gshape = (/nxg, nyg, nzg, nsp, 0/)
    offset = (/0, nyl*mod(nrank,nproc_j), nzl*(nrank/nproc_j), 0, 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'density'
    mpibuf1(1:lsize) = reshape(mom(1:1,nxgs:nxge,nys:nye,nzs:nze,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'den', 'f8', disp, &
                             dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! velocity
!$OMP PARALLEL WORKSHARE
    mom(2,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(2,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                        /mom(1,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    mom(3,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(3,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                        /mom(1,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    mom(4,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(4,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                        /mom(1,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
!$OMP END PARALLEL WORKSHARE

    nd     = 5
    lshape = (/3, nxg, nyl, nzl ,nsp/)
    gshape = (/3, nxg, nyg, nzg, nsp/)
    offset = (/0, 0, nyl*mod(nrank,nproc_j), nzl*(nrank/nproc_j), 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'velocity'
    mpibuf1(1:lsize) = reshape(mom(2:4,nxgs:nxge,nys:nye,nzs:nze,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'vel', 'f8', disp, &
                             dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! temperature
!$OMP PARALLEL WORKSHARE
    mom(5,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(5,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                        /mom(1,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    mom(6,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(6,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                        /mom(1,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    mom(7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                        /mom(1,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
!$OMP END PARALLEL WORKSHARE


    nd     = 5
    lshape = (/3, nxg, nyl, nzl, nsp/)
    gshape = (/3, nxg, nyg, nzg, nsp/)
    offset = (/0, 0, nyl*mod(nrank,nproc_j), nzl*(nrank/nproc_j), 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'temperature'
    mpibuf1(1:lsize) = reshape(mom(5:7,nxgs:nxge,nys:nye,nzs:nze,1:nsp), (/lsize/))
    call jsonio_put_metadata(json, p, 'temp', 'f8', disp, &
                             dsize, nd, gshape, desc)
    call mpiio_write_collective(fh, disp, nd, gshape, lshape, offset, mpibuf1)

    ! electromagnetic field on cell centers
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxgs,nxge
      tmp(1,i,j,k) = 2.5d-1*(+uf(1,i,j,k  )+uf(1,i,j+1,k  ) &
                             +uf(1,i,j,k+1)+uf(1,i,j+1,k+1))
      tmp(2,i,j,k) = 2.5d-1*(+uf(2,i,j,k  )+uf(2,i+1,j,k  ) &
                             +uf(2,i,j,k+1)+uf(2,i+1,j,k+1))
      tmp(3,i,j,k) = 2.5d-1*(+uf(3,i,j  ,k)+uf(3,i+1,j  ,k) &
                             +uf(3,i,j+1,k)+uf(3,i+1,j+1,k))
      tmp(4,i,j,k) = 5d-1*(+uf(4,i,j,k)+uf(4,i+1,j,k))
      tmp(5,i,j,k) = 5d-1*(+uf(5,i,j,k)+uf(5,i,j+1,k))
      tmp(6,i,j,k) = 5d-1*(+uf(6,i,j,k)+uf(6,i,j,k+1))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    nd     = 4
    lshape = (/6, nxg, nyl, nzl, 0/)
    gshape = (/6, nxg, nyg, nzg, 0/)
    offset = (/0, 0, nyl*mod(nrank,nproc_j), nzl*(nrank/nproc_j), 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    mpibuf1(1:lsize) = reshape(tmp(1:6,nxgs:nxge,nys:nye,nzs:nze), (/lsize/))
    call jsonio_put_metadata(json, p, 'uf', 'f8', disp, &
                             dsize, nd, gshape, desc)
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
    integer, intent(in) :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)

    call write_particle(up, uf, np2, it, 0, '_ptcl')

  end subroutine paraio__ptcl

  !
  ! output tracer particles
  !
  subroutine paraio__orb(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nzs+2)

    call write_particle(up, uf, np2, it, 1, '_orb')

  end subroutine paraio__orb

  !
  ! output particles
  !
  subroutine write_particle(up,uf,np2,it,mode,suffix)
    implicit none
    integer, intent(in)          :: it
    integer, intent(in)          :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in)          :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)          :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    integer, intent(in)          :: mode
    character(len=*), intent(in) :: suffix

    character(len=256) :: filename, jsonfile, datafile, desc, name
    integer(int64) :: disp, dsize, lsize, psize, gsize, goffset
    integer(int64) :: cumsum(nproc+1,nsp), ip1, ip2
    integer(int64) :: npl, npg, npo, nd8, gshape8(2)
    integer :: isp, irank
    integer :: fh, endian, nxg, nyg, nyl, nzg, nzl
    integer :: nd, lshape(5), gshape(5), offset(5)

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
    nzg = nzge - nzgs + 1
    nzl = nze  - nzs  + 1

    nd     = 4
    lshape = (/6, nxg, nyl, nzl, 0/)
    gshape = (/6, nxg, nyg, nzg, 0/)
    offset = (/0, 0, nyl*mod(nrank,nproc_j), nzl*(nrank/nproc_j), 0/)
    lsize  = product(lshape(1:nd))
    gsize  = product(gshape(1:nd))
    dsize  = gsize * 8
    desc   = 'electromagnetic field'
    mpibuf1(1:lsize) = reshape(uf(1:6,nxgs:nxge,nys:nye,nzs:nze), (/lsize/))
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
                                dsize, nd8, gshape8, desc)
       call mpiio_write_collective(fh, disp, gsize, lsize, psize, goffset, &
                                   8, transfer(mpibuf1(ip1:ip2), (/1/)))
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
    integer, intent(in)       :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in)       :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout)    :: buf(:)
    integer(8), intent(inout) :: cumsum(nproc+1,nsp)
    integer, intent(in)       :: mode

    integer :: i, j, k, ip, jp, isp
    integer(8) :: lcount(nsp), gcount(nsp, nproc)
    integer(8) :: pid

    ! count number of particles and pack into buffer
    if ( mode == 0 ) then
      ! * mode 0: all the active particles
      ip = 1
      do isp = 1, nsp
        lcount(isp) = 0
          do k = nzs, nze
          do j = nys, nye
            do i = 1, np2(j,k,isp)
              lcount(isp) = lcount(isp) + 1
              ! packing
              do jp = 1, ndim
                buf(ip) = up(jp,i,j,k,isp)
                ip = ip + 1
              enddo
            enddo
          enddo
          enddo
       enddo

    else if ( mode == 1 ) then
      ! * mode 1: tracer particles with positive IDs
      ip = 1
      do isp = 1, nsp
        lcount(isp) = 0
        do k = nzs, nze
        do j = nys, nye
          do i = 1, np2(j,k,isp)
            ! get particle ID as 64bit integer
            pid = transfer(up(ndim,i,j,k,isp), 1_8)

            ! count positive
            if( pid > 0 ) then
              lcount(isp) = lcount(isp) + 1
              ! packing
              do jp = 1, ndim
                buf(ip) = up(jp,i,j,k,isp)
                ip = ip + 1
              enddo
            endif
          enddo
        enddo
        enddo
      enddo

    else
      ! error
      write(0,*) 'Error: invalid mode specified for get_particle_count'
      call MPI_Finalize(mpierr)
      stop
    end if

    call MPI_Allgather(lcount, nsp, MPI_INTEGER8, gcount, nsp, MPI_INTEGER8, &
                       MPI_COMM_WORLD, mpierr)

    ! calculate cumulative sum
    do isp = 1, nsp
      cumsum(1,isp) = 0
      do i = 1, nproc
        cumsum(i+1,isp) = cumsum(i,isp) + gcount(isp,i)
      enddo
    enddo

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

    call jsonio_put_attribute(json, p, nzgs, 'nzgs', disp, '')
    call mpiio_write_atomic(file, disp, nzgs)

    call jsonio_put_attribute(json, p, nzge, 'nzge', disp, '')
    call mpiio_write_atomic(file, disp, nzge)

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
