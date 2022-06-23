module mpiio_test
  use mpi
  implicit none
  public

  character(len=*), parameter :: filename = 'mpiio_test.raw'
  integer, parameter :: ndim = 2
  integer, parameter :: nx   = 1024*4
  integer, parameter :: ny   = 1024*4
  integer, parameter :: nb   = 2
  integer, parameter :: nxs  = 1  + nb
  integer, parameter :: nxe  = nx + nb
  integer, parameter :: nys  = 1  + nb
  integer, parameter :: nye  = ny + nb
  integer, parameter :: zero(4) = (/0, 100, 200, 300/)

  integer :: file, comm, ierr, nproc, nrank, endian
  integer(MPI_OFFSET_KIND) :: disp

  integer :: mx, my, dims(ndim), coords(ndim)
  integer :: gshape(ndim), lshape(ndim), offset(ndim)

  integer :: charlen
  character(len=128) :: char
  integer(4) :: var_i4 = 1
  integer(8) :: var_i8 =-1
  real(4)    :: var_r4 = 4*atan(1.0_4)
  real(8)    :: var_r8 = 4*atan(1.0_8)
  integer(4) :: dat_i4(nx+2*nb,ny+2*nb)
  integer(8) :: dat_i8(nx+2*nb,ny+2*nb)
  real(4)    :: dat_r4(nx+2*nb,ny+2*nb)
  real(8)    :: dat_r8(nx+2*nb,ny+2*nb)
  integer(4) :: buf_i4(nx*ny)
  integer(8) :: buf_i8(nx*ny)
  real(4)    :: buf_r4(nx*ny)
  real(8)    :: buf_r8(nx*ny)


  interface filldata
     module procedure &
          & filldata_i4, &
          & filldata_i8, &
          & filldata_r4, &
          & filldata_r8
  end interface filldata

  interface testdata
     module procedure &
          & testdata_i4, &
          & testdata_i8, &
          & testdata_r4, &
          & testdata_r8
  end interface testdata

contains

  function formatstr(char, length) result(ret)
    implicit none
    character(len=*), intent(in) :: char
    integer, intent(in)          :: length
    character(len=length)        :: ret

    ret = char

  end function formatstr

  subroutine process_decomp(dims)
    implicit none
    integer, intent(inout) :: dims(ndim)

    dims(1) = 2
    dims(2) = nproc / 2

  end subroutine process_decomp

  subroutine filldata_i4(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    integer(4), intent(inout) :: data(:,:)
    integer, intent(in)       :: offset
    integer, intent(in)       :: nxs, nxe, nys, nye, nb
    integer, intent(in)       :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy

    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          data(ix,iy) = jx + jy*nx*dims(1) + offset
       end do
    end do

  end subroutine filldata_i4

  subroutine filldata_i8(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    integer(8), intent(inout) :: data(:,:)
    integer, intent(in)       :: offset
    integer, intent(in)       :: nxs, nxe, nys, nye, nb
    integer, intent(in)       :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy

    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          data(ix,iy) = jx + jy*nx*dims(1) + offset
       end do
    end do

  end subroutine filldata_i8

  subroutine filldata_r4(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    real(4), intent(inout) :: data(:,:)
    integer, intent(in)    :: offset
    integer, intent(in)    :: nxs, nxe, nys, nye, nb
    integer, intent(in)    :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy

    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          data(ix,iy) = jx + jy*nx*dims(1) + offset
       end do
    end do

  end subroutine filldata_r4

  subroutine filldata_r8(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    real(8), intent(inout) :: data(:,:)
    integer, intent(in)    :: offset
    integer, intent(in)    :: nxs, nxe, nys, nye, nb
    integer, intent(in)    :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy

    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          data(ix,iy) = jx + jy*nx*dims(1) + offset
       end do
    end do

  end subroutine filldata_r8

  subroutine testdata_i4(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    integer(4), intent(inout) :: data(:,:)
    integer, intent(in)       :: offset
    integer, intent(in)       :: nxs, nxe, nys, nye, nb
    integer, intent(in)       :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy, errcnt, sndbuf
    integer(4) :: element

    errcnt = 0
    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          element = jx + jy*nx*dims(1) + offset
          if( .not. data(ix,iy) == element ) then
             errcnt = errcnt + 1
          end if
       end do
    end do

    sndbuf = errcnt
    call MPI_Allreduce(sndbuf, errcnt, 1, MPI_INTEGER4, MPI_SUM, &
         & MPI_COMM_WORLD, ierr)

    if ( nrank == 0 .and. errcnt /= 0 ) then
          write(0, '(a, i8, a)') &
               & '# integer(4) read data failed with ', errcnt, ' errors'
    end if

  end subroutine testdata_i4

  subroutine testdata_i8(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    integer(8), intent(inout) :: data(:,:)
    integer, intent(in)       :: offset
    integer, intent(in)       :: nxs, nxe, nys, nye, nb
    integer, intent(in)       :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy, errcnt, sndbuf
    integer(8) :: element

    errcnt = 0
    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          element = jx + jy*nx*dims(1) + offset
          if( .not. data(ix,iy) == element ) then
             errcnt = errcnt + 1
          end if
       end do
    end do

    sndbuf = errcnt
    call MPI_Allreduce(sndbuf, errcnt, 1, MPI_INTEGER4, MPI_SUM, &
         & MPI_COMM_WORLD, ierr)

    if ( nrank == 0 .and. errcnt /= 0 ) then
          write(0, '(a, i8, a)') &
               & '# integer(8) read data failed with ', errcnt, ' errors'
    end if

  end subroutine testdata_i8

  subroutine testdata_r4(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    real(4), intent(inout)    :: data(:,:)
    integer, intent(in)       :: offset
    integer, intent(in)       :: nxs, nxe, nys, nye, nb
    integer, intent(in)       :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy, errcnt, sndbuf
    real(4) :: element

    errcnt = 0
    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          element = jx + jy*nx*dims(1) + offset
          if( .not. data(ix,iy) == element ) then
             errcnt = errcnt + 1
          end if
       end do
    end do

    sndbuf = errcnt
    call MPI_Allreduce(sndbuf, errcnt, 1, MPI_INTEGER4, MPI_SUM, &
         & MPI_COMM_WORLD, ierr)

    if ( nrank == 0 .and. errcnt /= 0 ) then
          write(0, '(a, i8, a)') &
               & '# real(4) read data failed with ', errcnt, ' errors'
    end if

  end subroutine testdata_r4

  subroutine testdata_r8(data, offset, nxs, nxe, nys, nye, nb, dims, coords)
    implicit none
    real(8), intent(inout)    :: data(:,:)
    integer, intent(in)       :: offset
    integer, intent(in)       :: nxs, nxe, nys, nye, nb
    integer, intent(in)       :: dims(ndim), coords(ndim)

    integer :: ix, iy, jx, jy, errcnt, sndbuf
    real(8) :: element

    errcnt = 0
    do iy = nys, nye
       jy = coords(2)*ny + iy-nys
       do ix = nxs, nxe
          jx = coords(1)*nx + ix-nxs
          element = jx + jy*nx*dims(1) + offset
          if( .not. data(ix,iy) == element ) then
             errcnt = errcnt + 1
          end if
       end do
    end do

    sndbuf = errcnt
    call MPI_Allreduce(sndbuf, errcnt, 1, MPI_INTEGER4, MPI_SUM, &
         & MPI_COMM_WORLD, ierr)

    if ( nrank == 0 .and. errcnt /= 0 ) then
          write(0, '(a, i8, a)') &
               & '# real(8) read data failed with ', errcnt, ' errors'
    end if

  end subroutine testdata_r8

end module mpiio_test
