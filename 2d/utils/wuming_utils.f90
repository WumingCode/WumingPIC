module wuming_utils
  use mpi
  implicit none
  private

  public :: get_etime
  public :: init_random_seed
  public :: uniform_rand
  public :: normal_rand
  public :: shuffle

  real(8), parameter :: pi = 4*atan(1.0d0)
  integer :: mpierr

contains
  !
  ! get elapsed time
  !
  function get_etime() result(y)
    implicit none
    real(8) :: y

    integer :: nrank

    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, mpierr)

    if ( nrank == 0 ) then
       y = MPI_Wtime()
    end if

    call MPI_Bcast(y, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpierr)

  end function get_etime

  !
  ! initialize raondom seed
  !
  subroutine init_random_seed()
    implicit none
    integer :: n
    integer, allocatable :: seed(:)

    integer :: nrank

    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, mpierr)

    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)*(nrank+1)
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_random_seed

  !
  ! uniform random number
  !
  function uniform_rand() result(y)
    implicit none
    real(8) :: y

    call random_number(y)

  end function uniform_rand

  !
  ! normal random number via Box-Muller method
  !
  function normal_rand() result(y)
    implicit none
    real(8) :: y
    real(8) :: xx(2), rr
    real(8), save :: yy
    logical, save :: has_saved

    if ( has_saved ) then
       y = yy
       has_saved = .false.
    else
       call random_number(xx)
       rr = sqrt(-2*log(1 - xx(1)) + 1.0d-30)
       yy = rr * cos(2*pi*xx(2))
       y  = rr * sin(2*pi*xx(2))
       has_saved = .true.
    end if

  end function normal_rand

  !
  ! quick sort
  !
  recursive subroutine qswap(key, idx, pivot, left, right)
    implicit none
    real(8), intent(inout) :: key(:)
    integer, intent(inout) :: idx(:)
    integer, intent(in)    :: pivot, left, right

    integer :: p, l, r, tmpidx
    real(8) :: pivkey

    pivkey = key( idx(pivot) )

    l = left
    r = right
    do
       do while(key( idx(l) ) < pivkey)
          l = l + 1
       end do

       do while(key( idx(r) ) > pivkey)
          r = r - 1
       end do

       if( l >= r ) then
          exit
       end if

       ! swap index
       tmpidx = idx(l)
       idx(l) = idx(r)
       idx(r) = tmpidx

       l = l + 1
       r = r - 1
    end do

    if( left < l-1 ) then
       p = left
       call qswap(key, idx, p, left, l-1)
    end if

    if( r+1 < right ) then
       p = r+1
       call qswap(key, idx, p, r+1, right)
    end if

  end subroutine qswap

  !
  ! quick sort for idx according to key
  !
  subroutine qsort(key, idx)
    implicit none
    real(8), intent(inout) :: key(:)
    integer, intent(inout) :: idx(:)

    call qswap(key, idx, 1, 1, size(key))

  end subroutine qsort

  !
  ! shuffle
  !
  subroutine shuffle(x)
    implicit none
    integer, intent(inout) :: x(:)

    integer :: i, idx(size(x)), val(size(x))
    real(8) :: key(size(x))

    do i = 1, size(x)
       idx(i) = i
       val(i) = x(i)
    end do

    call random_number(key)
    call qsort(key, idx)

    do i = 1, size(x)
       x(i) = val(idx(i))
    end do

  end subroutine shuffle

end module wuming_utils
