module const

  implicit none
  integer, parameter :: nx    = 8801     ! number of grid points in x
  integer, parameter :: ny    = 768      ! number of grid points in y
  integer, parameter :: nz    = 768      ! number of grid points in y
  integer, parameter :: nxgs  = 2         ! start point in x
  integer, parameter :: nxge  = nxgs+nx-1 ! end point
  integer, parameter :: nygs  = 2         ! start point in y
  integer, parameter :: nyge  = nygs+ny-1 ! end point
  integer, parameter :: nzgs  = 2         ! start point in z
  integer, parameter :: nzge  = nzgs+nz-1 ! end point
  integer, parameter :: np    = 100*nx     ! number of particles in each (y,z)
  integer, parameter :: nsp   = 2         ! number of particle species

  !number of processes
  integer, parameter :: nproc   = 9216
  integer, parameter :: nproc_i = 1 ! number of processes in x
  integer, parameter :: nproc_j = 96 ! number of processes in y
  integer, parameter :: nproc_k = nproc/(nproc_i*nproc_j) ! number of processes in z

end module

