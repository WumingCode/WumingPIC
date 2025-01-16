#ifdef USE_HDF5
module h5io
  use iso_fortran_env, only: int64
  use mpi
  use hdf5
  use h5util
  implicit none

  private

  public :: h5io__init
  public :: h5io__finalize
  public :: h5io__output
  public :: h5io__input
  public :: h5io__param
  public :: h5io__mom
  public :: h5io__ptcl
  public :: h5io__orb

  integer, parameter :: MAXDIM = 32
  integer :: hdferr

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nproc, nrank
  real(8), save :: delx, delt, u0, c
  real(8), allocatable :: q(:), r(:)
  character(len=256), save :: dir

  integer :: mpierr
  real(8), allocatable :: mpibuf(:)

contains

  !
  ! initialize module
  !
  subroutine h5io__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
                       nproc_in,nrank_in,                                                   &
                       delx_in,delt_in,c_in,q_in,r_in,dir_in)
    implicit none
    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nproc_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

    integer :: psize, fsize

    call h5util_init()

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
    allocate(mpibuf(max(psize, fsize)))

    is_init = .true.

  end subroutine h5io__init

  !
  ! finalize module
  !
  subroutine h5io__finalize()
    implicit none

    deallocate(mpibuf)
    call h5util_finalize()

  end subroutine h5io__finalize

  !
  ! output data for re-calculation
  !
  subroutine h5io__output(up,uf,np2,nxs,nxe,it,lflag)
    implicit none
    logical, intent(in) :: lflag
    integer, intent(in) :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in) :: it
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    character(len=256) :: filename

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(4), gdim(4), cdim(4), loff(4), goff(4)
    integer :: nd, nxg, nyg, nyl

    if(.not.is_init)then
       write(0,*) 'Initialize first by calling h5io__init()'
       stop
    endif

    !filename
    if(lflag)then
       write(filename,'(a,i7.7,a)')trim(dir),9999999,'.h5'
    else
       write(filename,'(a,i7.7,a)')trim(dir),it,'.h5'
    endif
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    !time & parameters
    call h5util_put_attribute(file_id, "it", it)
    call h5util_put_attribute(file_id, "ndim", ndim)
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "nxs", nxs)
    call h5util_put_attribute(file_id, "nxe", nxe)
    call h5util_put_attribute(file_id, "nxgs", nxgs)
    call h5util_put_attribute(file_id, "nxge", nxge)
    call h5util_put_attribute(file_id, "nygs", nygs)
    call h5util_put_attribute(file_id, "nyge", nyge)
    call h5util_put_attribute(file_id, "nsp", nsp)
    call h5util_put_attribute(file_id, "nproc", nproc)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)

    !particle
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 4
    ldim = (/ndim, np, nyl, nsp/)
    gdim = (/ndim, np, nyg, nsp/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)
    call h5util_create_dataset_chunk(file_id, "up", H5T_NATIVE_DOUBLE, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "up", nd, ldim, ldim, &
         & loff, goff, reshape(up, (/size(up)/)))

    nd   = 2
    ldim = (/nyl, nsp, 0, 0/)
    gdim = (/nyg, nsp, 0, 0/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/nyl*nrank, 0, 0, 0/)
    call h5util_create_dataset_chunk(file_id, "np2", H5T_NATIVE_INTEGER, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "np2", nd, ldim, ldim, &
         & loff, goff, reshape(np2, (/size(np2)/)))

    !field
    nxg  = nxge - nxgs + 5
    nyl  = nye  - nys  + 5
    nyg  = nyl * nproc

    nd   = 3
    ldim = (/6, nxg, nyl, 0/)
    gdim = (/6, nxg, nyg, 0/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)
    call h5util_create_dataset_chunk(file_id, "uf", H5T_NATIVE_DOUBLE, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, &
         & reshape(uf, (/6*nxg*nyl/)))

    call h5util_close_file(file_id)

  end subroutine h5io__output

  !
  ! input data for re-calculation
  !
  subroutine h5io__input(up,uf,np2,indim,nxs,nxe,it,file)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: np2(nys:nye,nsp), nxs, nxe, it, indim
    real(8), intent(out) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: inp, inxgs, inxge, inygs, inyge, insp, inproc

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(4), loff(4), goff(4)
    integer :: i, ii, nd, nxg, nyl, lsize

    if(.not.is_init)then
       write(0,*) 'Initialize first by calling h5io__init()'
       stop
    endif

    !open file
    call h5util_open_file(trim(dir)//trim(file), file_id)

    !time & parameters
    call h5util_get_attribute(file_id, "it", it)
    call h5util_get_attribute(file_id, "ndim", indim)
    call h5util_get_attribute(file_id, "np", inp)
    call h5util_get_attribute(file_id, "nxs", nxs)
    call h5util_get_attribute(file_id, "nxe", nxe)
    call h5util_get_attribute(file_id, "nxgs", inxgs)
    call h5util_get_attribute(file_id, "nxge", inxge)
    call h5util_get_attribute(file_id, "nygs", inygs)
    call h5util_get_attribute(file_id, "nyge", inyge)
    call h5util_get_attribute(file_id, "nsp", insp)
    call h5util_get_attribute(file_id, "nproc", inproc)
    call h5util_get_attribute(file_id, "delx", delx)
    call h5util_get_attribute(file_id, "delt", delt)
    call h5util_get_attribute(file_id, "c", c)
    call h5util_get_attribute(file_id, "m", r)
    call h5util_get_attribute(file_id, "q", q)

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

    !particle data
    nyl  = nye-nys+1

    nd   = 4
    ldim = (/indim, np, nyl, nsp/)
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)
    lsize = indim*np*nyl*nsp
    call h5util_read_dataset(file_id, "up", nd, ldim, ldim, loff, goff, mpibuf)
    up = reshape(mpibuf(1:lsize), ldim(1:4))

    nd   = 2
    ldim = (/nyl, nsp, 0, 0/)
    loff = (/0, 0, 0, 0/)
    goff = (/nyl*nrank, 0, 0, 0/)
    lsize = nyl*nsp
    call h5util_read_dataset(file_id, "np2", nd, ldim, ldim, loff, goff, mpibuf)
    np2 = reshape(mpibuf(1:lsize), ldim(1:2))

    !field data
    nxg  = nxge - nxgs + 5
    nyl  = nye  - nys  + 5

    nd   = 3
    ldim = (/6, nxg, nyl, 0/)
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)
    lsize = 6*nxg*nyl
    call h5util_read_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, mpibuf)
    uf = reshape(mpibuf(1:lsize), ldim(1:3))

    call h5util_close_file(file_id)

  end subroutine h5io__input

  !
  ! output parameters
  !
  subroutine h5io__param(n0,np2,temp,rtemp,fpe,fge,ls,file,nroot)
    implicit none
    integer, intent(in)          :: n0, nroot
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: temp, rtemp, fpe, fge, ls
    character(len=*), intent(in) :: file
    integer :: isp
    real(8) :: pi, fpi, fgi, vti, vte, vai, vae

    integer(hid_t) :: file_id
    character(len=256) :: filename

    if(.not.is_init)then
       write(0,*) 'Initialize first by calling h5io__init()'
       stop
    endif

    pi   = 4*atan(1.0D0)
    vti  = sqrt(2*temp/r(1))
    vte  = sqrt(2*temp*rtemp/r(2))
    vai  = abs(fge*r(1)*c/q(1))/sqrt(4*pi*r(1)*n0) * r(2)/r(1)
    vae  = abs(fge*r(2)*c/q(2))/sqrt(4*pi*r(2)*n0)
    fpi  = fpe * sqrt(r(2)/r(1))
    fgi  = fge * r(2)/r(1)

    filename = trim(dir)//trim(file)//trim(".h5")
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    call h5util_put_attribute(file_id, "nx", nxge-nxgs+1)
    call h5util_put_attribute(file_id, "ny", nyge-nygs+1)
    call h5util_put_attribute(file_id, "ls", ls)
    call h5util_put_attribute(file_id, "np2", np2(nys,:))
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "r", r)
    call h5util_put_attribute(file_id, "q", q)
    call h5util_put_attribute(file_id, "fpe", fpe)
    call h5util_put_attribute(file_id, "fge", fge)
    call h5util_put_attribute(file_id, "fpi", fpi)
    call h5util_put_attribute(file_id, "fgi", fgi)
    call h5util_put_attribute(file_id, "vai", vai)
    call h5util_put_attribute(file_id, "vae", vae)
    call h5util_put_attribute(file_id, "vte", vte)
    call h5util_put_attribute(file_id, "vti", vti)
    call h5util_put_attribute(file_id, "beta", (vti/vai)**2)
    call h5util_put_attribute(file_id, "rtemp", rtemp)
    call h5util_put_attribute(file_id, "rgi", vti/(fge*r(2)/r(1)))
    call h5util_put_attribute(file_id, "n0", n0)

    call h5util_close_file(file_id)

  end subroutine h5io__param

  !
  ! output field and moment quantities
  !
  subroutine h5io__mom(den,vel,temp,uf,it)
    implicit none
    integer, intent(in)    :: it
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nsp),    &
                              vel(nxgs-1:nxge+1,nys-1:nye+1,3,nsp),  &
                              temp(nxgs-1:nxge+1,nys-1:nye+1,3,nsp)
    integer :: i, j, isp
    real(8) :: tmp(1:6,nxgs:nxge,nys:nye)
    character(len=256) :: filename

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(4), gdim(4), cdim(4), loff(4), goff(4)
    integer :: nd, nxg, nyg, nyl


    if(.not.is_init)then
       write(0,*) 'Initialize first by calling h5io__init()'
       stop
    endif

    write(filename,'(a,i7.7,a)')trim(dir),it,'_mom.h5'
    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    !time & parameters
    call h5util_put_attribute(file_id, "it", it)
    call h5util_put_attribute(file_id, "ndim", ndim)
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "nxgs", nxgs)
    call h5util_put_attribute(file_id, "nxge", nxge)
    call h5util_put_attribute(file_id, "nygs", nygs)
    call h5util_put_attribute(file_id, "nyge", nyge)
    call h5util_put_attribute(file_id, "nsp", nsp)
    call h5util_put_attribute(file_id, "nproc", nproc)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)

    !density
    nxg  = nxge-nxgs+1
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 3
    ldim = (/nxg, nyl, nsp, 0/)
    gdim = (/nxg, nyg, nsp, 0/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, nyl*nrank, 0, 0/)
    call h5util_create_dataset_chunk(file_id, "den", H5T_NATIVE_REAL, &
        & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "den", nd, ldim, ldim, loff, goff, &
        & reshape(sngl(den(nxgs:nxge,nys:nye,1:nsp)), (/nxg*nyl*nsp/)))

    !velocity
    !$OMP PARALLEL DO PRIVATE(i,j,isp)
    do isp = 1, nsp
      do j = nys,nye
        do i = nxgs,nxge
          vel(i,j,1,isp) = vel(i,j,1,isp)/den(i,j,isp)
          vel(i,j,2,isp) = vel(i,j,2,isp)/den(i,j,isp)
          vel(i,j,3,isp) = vel(i,j,3,isp)/den(i,j,isp)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nd   = 4
    ldim = (/nxg, nyl, 3, nsp/)
    gdim = (/nxg, nyg, 3, nsp/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, nyl*nrank, 0, 0/)
    call h5util_create_dataset_chunk(file_id, "vel", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "vel", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(vel(nxgs:nxge,nys:nye,1:3,1:nsp)), (/nxg*nyl*3*nsp/)))

    !temperature
    !$OMP PARALLEL DO PRIVATE(i,j,isp)
    do isp = 1, nsp
      do j = nys,nye
        do i = nxgs,nxge
          temp(i,j,1,isp) = temp(i,j,1,isp)/den(i,j,isp)
          temp(i,j,2,isp) = temp(i,j,2,isp)/den(i,j,isp)
          temp(i,j,3,isp) = temp(i,j,3,isp)/den(i,j,isp)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nd   = 4
    ldim = (/nxg, nyl, 3, nsp/)
    gdim = (/nxg, nyg, 3, nsp/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, nyl*nrank, 0, 0/)
    call h5util_create_dataset_chunk(file_id, "temp", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "temp", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(temp(nxgs:nxge,nys:nye,1:3,1:nsp)), (/nxg*nyl*3*nsp/)))

    ! electromagnetic field on cell centers
    do j=nys,nye
       do i=nxgs,nxge
          tmp(1,i,j) = (uf(1,i,j)+uf(1,i,j+1)) / 2
          tmp(2,i,j) = (uf(2,i,j)+uf(2,i+1,j)) / 2
          tmp(3,i,j) = (uf(3,i,j)+uf(3,i+1,j)+uf(3,i+1,j)+uf(3,i+1,j+1)) / 2
          tmp(4,i,j) = (uf(4,i,j)+uf(4,i+1,j)) / 2
          tmp(5,i,j) = (uf(5,i,j)+uf(5,i,j+1)) / 2
          tmp(6,i,j) = uf(6,i,j)
       enddo
    enddo

    nd   = 3
    ldim = (/6, nxg, nyl, 0/)
    gdim = (/6, nxg, nyg, 0/)
    cdim = ldim
    loff = (/0, 0, 0, 0/)
    goff = (/0, 0, nyl*nrank, 0/)
    call h5util_create_dataset_chunk(file_id, "uf", H5T_NATIVE_REAL, &
         & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, &
         & reshape(sngl(tmp(1:6,nxgs:nxge,nys:nye)), (/6*nxg*nyl/)))

    call h5util_close_file(file_id)

  end subroutine h5io__mom

  !
  ! output all the active particles
  !
  subroutine h5io__ptcl(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    call write_particle(up, uf, np2, it, 0, '_ptcl')

  end subroutine h5io__ptcl

  !
  ! output tracer particles
  !
  subroutine h5io__orb(up,uf,np2,it)
    implicit none
    integer, intent(in) :: it
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)

    call write_particle(up, uf, np2, it, 1, '_orb')

  end subroutine h5io__orb

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

    integer :: isp, ii, i, j
    character(len=256) :: filename, name

    integer(hid_t) :: file_id
    integer(hsize_t) :: ldim(3), gdim(3), cdim(3), loff(3), goff(3)
    integer :: nd, nxg, nyg, nyl, nrec, npl, npg, npo, irank
    integer(int64) :: cumsum(nproc+1,nsp), ip1, ip2

    if(.not.is_init)then
       write(0,*) 'Initialize first by calling h5io__init()'
       stop
    endif

    !filename
    write(filename,'(a, i7.7, a)') trim(dir), it, trim(suffix) // '.h5'

    call h5util_create_file(filename)
    call h5util_open_file(filename, file_id)

    !time & parameters
    call h5util_put_attribute(file_id, "it", it)
    call h5util_put_attribute(file_id, "ndim", ndim)
    call h5util_put_attribute(file_id, "np", np)
    call h5util_put_attribute(file_id, "nxgs", nxgs)
    call h5util_put_attribute(file_id, "nxge", nxge)
    call h5util_put_attribute(file_id, "nygs", nygs)
    call h5util_put_attribute(file_id, "nyge", nyge)
    call h5util_put_attribute(file_id, "nsp", nsp)
    call h5util_put_attribute(file_id, "nproc", nproc)
    call h5util_put_attribute(file_id, "delx", delx)
    call h5util_put_attribute(file_id, "delt", delt)
    call h5util_put_attribute(file_id, "c", c)
    call h5util_put_attribute(file_id, "m", r)
    call h5util_put_attribute(file_id, "q", q)

    !field
    nxg  = nxge-nxgs+5
    nyg  = nyge-nygs+1
    nyl  = nye-nys+1

    nd   = 3
    ldim = (/6, nxg, nyl/)
    gdim = (/6, nxg, nyg/)
    cdim = ldim
    loff = (/0, 0, 0/)
    goff = (/0, 0, nyl*nrank/)
    call h5util_create_dataset_chunk(file_id, "uf", H5T_NATIVE_DOUBLE, &
        & nd, gdim, cdim)
    call h5util_write_dataset(file_id, "uf", nd, ldim, ldim, loff, goff, &
        & reshape(uf(1:6,nxgs-2:nxge+2,nys:nye), (/6*nxg*nyl/)))

    ! particle
    call get_particle_count(up, np2, mpibuf, cumsum, mode)

    irank = nrank + 1
    ip2   = 0
    do isp = 1, nsp
       write(name, '("up", i2.2)') isp

       npl    = cumsum(irank+1,isp) - cumsum(irank,isp)
       npg    = cumsum(nproc+1,isp)
       npo    = cumsum(irank,isp)
       ip1    = ip2 + 1
       ip2    = ip1 + ndim*npl - 1

       nd   = 2
       ldim = (/ndim, npl, 0/)
       gdim = (/ndim, npg, 0/)
       loff = (/0, 0, 0/)
       goff = (/0, npo, 0/)
       call h5util_create_dataset(file_id, trim(name), H5T_NATIVE_DOUBLE, &
            & nd, gdim)
       call h5util_write_dataset(file_id, trim(name), nd, ldim, ldim, &
            & loff, goff, mpibuf(ip1:ip2))
    end do

    call h5util_close_file(file_id)

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

end module h5io
#else

module h5io
  ! empty module
end module h5io

#endif
