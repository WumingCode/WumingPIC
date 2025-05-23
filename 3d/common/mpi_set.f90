module mpi_set

  use mpi
  implicit none
  private
  public :: mpi_set__init, MPI_WTIME

  integer, public, parameter :: mnpi  = MPI_INTEGER
  integer, public, parameter :: mnpr  = MPI_DOUBLE_PRECISION
  integer, public, parameter :: mnpc  = MPI_CHARACTER
  integer, public, parameter :: opsum = MPI_SUM
  integer, public            :: nerr, ncomw, nstat(MPI_STATUS_SIZE)
  integer, public            :: nrank, nrank_j, nrank_k
  integer, public            :: jup, jdown, kup, kdown
  integer, public            :: nys, nye, nzs, nze


contains


  subroutine mpi_set__init(nygs,nyge,nzgs,nzge,nproc,nproc_j,nproc_k)

    integer, intent(in) :: nygs, nyge, nzgs, nzge
    integer, intent(in) :: nproc, nproc_j, nproc_k
    integer             :: j, k, irank, nsize
    integer             :: ptable(-1:nproc_j,-1:nproc_k)

    !*********** Initialization for MPI  ***************!
    call MPI_INIT(nerr)
    ncomw = MPI_COMM_WORLD
    call MPI_COMM_SIZE(ncomw, nsize, nerr)
    call MPI_COMM_RANK(ncomw, nrank, nerr)

    if(nsize /= nproc) then
       write(*,*) 'error in proc no.'
       call MPI_ABORT(ncomw, 9, nerr)
       call MPI_FINALIZE(nerr)
    endif
    if(nsize /= (nproc_j*nproc_k)) then
       write(*,*) 'error in proc no.'
       call MPI_ABORT(ncomw, 9, nerr)
       call MPI_FINALIZE(nerr)
    endif

    !rank table
    do k=-1,nproc_k
    do j=-1,nproc_j
       ptable(j,k) = MPI_PROC_NULL
    enddo
    enddo
    irank = 0
    do j=0,nproc_j-1
    do k=0,nproc_k-1
       ptable(j,k) = irank
       if(nrank == irank)then
          nrank_j = j
          nrank_k = k
       endif
       irank = irank+1
    enddo
    enddo

    call para_range(nys,nye,nygs,nyge,nproc_j,nrank_j)
    call para_range(nzs,nze,nzgs,nzge,nproc_k,nrank_k)

    !For MPI_SENDRECV
    jup   = ptable(nrank_j+1,nrank_k)
    jdown = ptable(nrank_j-1,nrank_k)
    kup   = ptable(nrank_j,nrank_k+1)
    kdown = ptable(nrank_j,nrank_k-1)

    ! periodic boundary conditions in y and z directions
    if(nrank_j == 0) jdown = ptable(nproc_j-1,nrank_k)
    if(nrank_k == 0) kdown = ptable(nrank_j,nproc_k-1)
    if(nrank_j == nproc_j-1) jup = ptable(0,nrank_k)
    if(nrank_k == nproc_k-1) kup = ptable(nrank_j,0)

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

