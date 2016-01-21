module mpi_set

  implicit none
  private

  public :: mpi_set__init

  include 'mpif.h'

  integer, public, parameter :: mnpi  = MPI_INTEGER
  integer, public, parameter :: mnpr  = MPI_DOUBLE_PRECISION
  integer, public, parameter :: mnpc  = MPI_CHARACTER
  integer, public, parameter :: opsum = MPI_SUM
  integer, public            :: nerr, ncomw, nstat(MPI_STATUS_SIZE)
  integer, public            :: nrank, nrank_i, nrank_j, nrank_k
  integer, public            :: iup, idown, jup, jdown, kup, kdown
  integer, public            :: nys, nye, nzs, nze


contains


  subroutine mpi_set__init(nygs,nyge,nzgs,nzge, &
                           nproc,nproc_i,nproc_j,nproc_k)

    integer, intent(in) :: nygs, nyge, nzgs, nzge
    integer, intent(in) :: nproc, nproc_i, nproc_j, nproc_k
    integer             :: i, j, k, irank, nsize 
    integer             :: ptable(-1:nproc_i,-1:nproc_j,-1:nproc_k)

    !*********** Initialization for MPI  ***************!
    call MPI_INIT(nerr)
    ncomw = MPI_COMM_WORLD
    call MPI_COMM_SIZE(ncomw, nsize, nerr)
    call MPI_COMM_RANK(ncomw, nrank, nerr)

    if(nsize /= nproc) then
       stop 'error in proc no.'
       call MPI_ABORT(ncomw, 9, nerr)
       call MPI_FINALIZE(nerr)
    endif
    if(nsize /= (nproc_i*nproc_j*nproc_k)) then
       stop 'error in proc no.'
       call MPI_ABORT(ncomw, 9, nerr)
       call MPI_FINALIZE(nerr)
    endif

    !rank table
    do k=-1,nproc_k
    do j=-1,nproc_j
    do i=-1,nproc_i
       ptable(i,j,k) = MPI_PROC_NULL
    enddo
    enddo
    enddo
    irank = 0
    do i=0,nproc_i-1
    do j=0,nproc_j-1
    do k=0,nproc_k-1
       ptable(i,j,k) = irank
       if(nrank == irank)then
          nrank_i = i
          nrank_j = j
          nrank_k = k
       endif
       irank = irank+1
    enddo
    enddo
    enddo

    call para_range(nys,nye,nygs,nyge,nproc_j,nrank_j)
    call para_range(nzs,nze,nzgs,nzge,nproc_k,nrank_k)

    !For MPI_SENDRECV
    iup   = ptable(nrank_i+1,nrank_j,nrank_k)
    idown = ptable(nrank_i-1,nrank_j,nrank_k)
    jup   = ptable(nrank_i,nrank_j+1,nrank_k)
    jdown = ptable(nrank_i,nrank_j-1,nrank_k)
    kup   = ptable(nrank_i,nrank_j,nrank_k+1)
    kdown = ptable(nrank_i,nrank_j,nrank_k-1)

    ! periodic boundary conditions in y and z directions
    if(nrank_j == 0) jdown = ptable(nrank_i,nproc_j-1,nrank_k)   
    if(nrank_k == 0) kdown = ptable(nrank_i,nrank_j,nproc_k-1)
    if(nrank_j == nproc_j-1) jup = ptable(nrank_i,0,nrank_k)
    if(nrank_k == nproc_k-1) kup = ptable(nrank_i,nrank_j,0)

  end subroutine mpi_set__init


  subroutine para_range(ns,ne,n1,n2,isize,irank)

    integer, intent(in)  :: n1, n2, isize, irank
    integer, intent(out) :: ns, ne
    integer              :: iwork1, iwork2

    !start and end of loop counter
    iwork1 = (n2-n1+1)/isize
    iwork2 = mod(n2-n1+1,isize)
    ns = irank*iwork1+n1+min(irank,iwork2)
    ne = ns+iwork1-1
    if(iwork2 > irank) ne = ne+1

  end subroutine para_range


end module mpi_set

