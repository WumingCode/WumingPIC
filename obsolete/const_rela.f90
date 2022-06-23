module const

  implicit none

!!************************ NUMERICAL CONSTANTS ************************!!
  integer, parameter :: nx      = 20001      ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny      = 1024      ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nz      = 1024       ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs    = 2         ! START POINT IN X
  integer, parameter :: nxge    = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs    = 2         ! START POINT IN Y
  integer, parameter :: nyge    = nygs+ny-1 ! END POINT
  integer, parameter :: nzgs    = 2         ! START POINT IN Z
  integer, parameter :: nzge    = nzgs+nz-1 ! END POINT
  integer, parameter :: np      = 30*nx    ! MAX. NUMBER OF PARTICLES IN COLUMN AT (Y,Z)
  integer, parameter :: nsp     = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: nproc   = 16384      ! NUMBER OF PROCESSES
  integer, parameter :: nproc_i = 1         ! NUMBER OF PROCESSES IN X
  integer, parameter :: nproc_j = 128        ! NUMBER OF PROCESSES IN Y
  integer, parameter :: nproc_k = nproc/(nproc_i*nproc_j) ! NUMBER OF PROCESSES IN Z
  integer, parameter :: nroot   = 0         ! ROOT PROC. NUMBER

!! SETUP FOR MOVING INJECTOR; INITIAL DOMAIN SIZE IN X
  integer            :: nxs     = nxgs        ! START POINT IN X
  integer            :: nxe     = nxgs+2500 ! INITIAL SIZE

!! SETUP FOR MAIN PROGRAM & SUBROUTINES
  integer, parameter :: itmax   = 20000    !NUMBER OF ITERATION
  integer            :: it0     = 0	!0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1  = 4000    !INTERVAL FOR PARTICLES & FIELDS STORAGE          
  integer, parameter :: intvl2  = 1         !INTERVAL FOR INJECTING PARTICLES
  integer, parameter :: intvl3  = 1         !INTERVAL FOR EXPANDING PHYSICAL REGION IN X
  integer, parameter :: intvl4  = 1000      !INTERVAL FOR RECORDING MOMENT DATA
  character(len=128) :: dir     = './'      !DIRECTORY FOR OUTPUT
  character(len=128) :: file9   = 'init_param.dat' !FILENAME OF INIT CONDITIONS
  real(8), parameter :: etlim   = 23.5*60.*60. !MAX. ELAPSE TIME IN SEC.

!! OTHER CONSTANTS
  real(8)            :: c       = 1.0D0     !SPEED OF LIGHT
  real(8), parameter :: gfac    = 0.501D0   !IMPLICITNESS FACTOR 0.501-0.505
  real(8), parameter :: cfl     = 1.0D0     !CFL CONDITION FOR LIGHT WAVE
  real(8)            :: delx    = 1.0D0     !CELL WIDTH
  real(8), parameter :: rdbl    = 20.0      !ELECTRON SKIN DEPTH / CELL WIDTH
  
!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!      mr : ION-TO-ELECTRON MASS RATIO
!!    gam0 : UPSTREAM BULK LORENTZ FACTOR
!!    sig0 : UPSTREAM SIGMA PARAMETER USING ION MASS
!!   rtemp : Te/Ti
!!     vt? : (ION/ELECTRON) THERMAL SPEED
!!   theta : SHOCK ANGLE
!!     phi : INCLINATION ANGLE FROM X-Y PLANE FOR BY & BZ
!!    ldmp : FOR INITIAL VELOCITY PROFILE POSITION NEAR X=0 IN UNIT OF C/WPI
  integer, parameter :: n0     = 10
  real(8), parameter :: mr     = 100.0D0
  real(8), parameter :: gam0   = 40.0D0
  real(8), parameter :: sig0   = 5.D-2
  real(8), parameter :: rtemp  = 1.0D0
  real(8), parameter :: vte    = 0.0D0
  real(8), parameter :: vti    = sqrt(1.0D0/(rtemp*mr))*vte
  real(8), parameter :: theta  = 90.0D0
  real(8), parameter :: phi    = 90.0D0
  real(8)            :: ldmp   = 10.0
!!TRACKING PARTICLES INITIALLY XRS <= X <= XRE ONLY WHEN NDIM=6
  real(8), parameter :: xrs   = 530.0D0
  real(8), parameter :: xre   = 540.0D0

end module

