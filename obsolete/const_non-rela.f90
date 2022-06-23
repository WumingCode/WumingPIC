module const

  implicit none

!!************************ NUMERICAL CONSTANTS ************************!!
  integer, parameter :: nx      = 8801      ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny      = 32       ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nz      = 32       ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs    = 2         ! START POINT IN X
  integer, parameter :: nxge    = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs    = 2         ! START POINT IN Y
  integer, parameter :: nyge    = nygs+ny-1 ! END POINT
  integer, parameter :: nzgs    = 2         ! START POINT IN Z
  integer, parameter :: nzge    = nzgs+nz-1 ! END POINT
  integer, parameter :: np      = 100*nx    ! MAX. NUMBER OF PARTICLES IN COLUMN AT (Y,Z)
  integer, parameter :: nsp     = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc   = 16      ! NUMBER OF PROCESSES
  integer, parameter :: nproc_i = 1         ! NUMBER OF PROCESSES IN X
  integer, parameter :: nproc_j = 4        ! NUMBER OF PROCESSES IN Y
  integer, parameter :: nproc_k = nproc/(nproc_i*nproc_j) ! NUMBER OF PROCESSES IN Z
  integer, parameter :: nroot   = 0         ! ROOT PROC. NUMBER

!! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs     = nxgs        ! START POINT IN X
  integer            :: nxe     = nxgs+nx*0.5 ! INITIAL SIZE

!! SETUP FOR MAIN PROGRAM & SUBROUTINES
  integer, parameter :: itmax   = 175000    !NUMBER OF ITERATION
  integer            :: it0     = 0	!0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1  = 175000    !INTERVAL FOR PARTICLES & FIELDS STORAGE          
  integer, parameter :: intvl2  = 175000    !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3  = 20        !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  integer, parameter :: intvl4  = 175000    !INTERVAL FOR RECORDING MOMENT DATA
  character(len=128) :: dir     = './'      !DIRECTORY FOR OUTPUT
  character(len=128) :: file9   = 'init_param.dat' !FILENAME OF INIT CONDITIONS
  real(8), parameter :: etlim   = 1.0*60.*60. !MAX. ELAPSE TIME IN SEC.

!! OTHER CONSTANTS
  real(8)            :: c       = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac    = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl     = 0.5D0     !CFL CONDITION FOR LIGHT WAVE
  real(8)            :: delx    = 1.0D0     !CELL WIDTH
  real(8), parameter :: rdbl    = 1.0       !DEBYE LENGTH / CELL WIDTH

!!************************ PHYSICAL CONSTANTS ************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!      mr : ION-TO-ELECTRON MASS RATIO
!!   alpha : wpe/wge = c/vth_e * sqrt(beta_e)
!!    beta : ION PLASMA BETA
!!   rtemp : Te/Ti
!!      ma : ALFVEN MACH NUMBER ~ (SIGMA*GAMMA)^(-1/2) FOR V0~C
!!   theta : SHOCK ANGLE(deg.)
!!   phi   : INCLINATION ANGLE(deg.) FROM X-Y PLANE FOR BY & BZ
  integer, parameter :: n0     = 20
  real(8), parameter :: mr     = 64.0D0
  real(8), parameter :: alpha  = 10.0D0, beta = 0.5D0, rtemp=1.0D0
  real(8), parameter :: ma     = 16.0D0
  real(8), parameter :: theta  = 70.0D0
  real(8), parameter :: phi    = 90.0D0

end module

