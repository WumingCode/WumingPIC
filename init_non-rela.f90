module init

  use const
  use mpi_set

  implicit none

  private

  public :: init__set_param, init__inject, init__relocate

  integer, allocatable, public :: np2(:,:,:), cumcnt(:,:,:,:)
  real(8), public              :: delt, u0
  real(8), public              :: q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:,:)
  real(8), allocatable, public :: up(:,:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:,:)
  real(8), allocatable, public :: den(:,:,:,:),vel(:,:,:,:,:),temp(:,:,:,:,:)
  real(8), parameter           :: pi = 4.0D0*datan(1.0D0), deg2rad = pi/180.0D0
  real(8), save                :: v0, b0, vti, vte, gam0


contains

  
  subroutine init__set_param

    use fio,  only : fio__input, fio__param
    use sort, only : sort__bucket

    integer              :: n, i, j, k, isp
    integer, allocatable :: seed(:)
    real(8)              :: fgi, fpi, fpe, fge, va, rgi, rge, ldb
    character(len=128)   :: file11

!************** MPI settings  *******************!
    call mpi_set__init(nygs,nyge,nzgs,nzge, &
                       nproc,nproc_i,nproc_j,nproc_k)
!*********** End of MPI settings  ***************!

!*********** Memory Allocations  ****************!
    allocate(np2(nys:nye,nzs:nze,1:nsp))
    allocate(cumcnt(nxgs:nxge,nys:nye,nzs:nze,1:nsp))
    allocate(uf(1:6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
    allocate(up(1:6,np,nys:nye,nzs:nze,1:nsp))
    allocate(gp(1:6,np,nys:nye,nzs:nze,1:nsp))
    allocate(den(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp))
    allocate(vel(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp))
    allocate(temp(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp))

!$OMP WORKSHARE
    up(1:6,1:np,nys:nye,nzs:nze,1:nsp) = 0.0d0
    den(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = 0.0d0
    vel(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = 0.0d0
    temp(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = 0.0d0
!$OMP END WORKSHARE

!***************** End of  **********************!

!*********** Random seed *************!
    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)*(nrank+1)
    call random_seed(put=seed)
    deallocate(seed)
!***********   End of    *************!

!**** SETTING OTHER NUMERICAL & PHYSICAL CONSTANTS ****!
    r(2) = 1.0D0      ! ELECTRON MASS
    r(1) = r(2)*mr    ! ION MASS
    delt = cfl*delx/c ! TIME STEP SIZE
    ldb  = delx*rdbl
    fpe  = sqrt(beta*rtemp)*c/(sqrt(2.D0)*alpha*ldb)
    fge  = fpe/alpha
    fgi  = fge*r(2)/r(1)
    fpi  = fpe*sqrt(r(2)/r(1))
    va   = fge/fpe*c*sqrt(r(2)/r(1))
    rge  = alpha*ldb*sqrt(2.D0)
    rgi  = rge*sqrt(r(1)/r(2))/sqrt(rtemp)
    vte  = rge*fge
    vti  = vte*sqrt(r(2)/r(1))/sqrt(rtemp)

    !CHARGE
    q(1) = fpi*sqrt(r(1)/(4.0D0*pi*n0))
    q(2) = -q(1)

    !MAGNETIC FIELD STRENGTH
    b0 = fgi*r(1)*c/q(1)
    b0 = b0/sin(84.D0*deg2rad)

    !INJECTOR IS ON THE RIGHT-HAND-SIDE; MINUS SIGN IS NECESSARY
    v0   = -ma*va
    u0   = v0/sqrt(1.-(v0/c)**2)
    gam0 = sqrt(1.+u0**2/c**2)

    !INITIAL NUMBER OF PARTICLES IN COLUMN AT (Y, Z)
    np2(nys:nye,nzs:nze,1:nsp) = n0*(nxe-nxs-2)
    if(nrank == nroot)then
       if(n0*(nxge-nxgs) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif
 
    !PREPARATION FOR SORT
    do isp=1,nsp
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=nzs,nze
       do j=nys,nye
          cumcnt(nxs:nxs+1,j,k,isp) = 0
          do i=nxs+2,nxe-1
             cumcnt(i,j,k,isp) = cumcnt(i-1,j,k,isp)+n0
          enddo
          cumcnt(nxe,j,k,isp) = cumcnt(nxe-1,j,k,isp)
          if(cumcnt(nxe,j,k,isp) /= np2(j,k,isp))then
             write(*,*)'error in cumcnt'
             stop
          endif
       enddo
       enddo
!$OMP END PARALLEL DO
    enddo

    if(it0 /= 0)then
       !RESTART FROM THE PAST CALCULATION
       write(file11,'(i7.7,a,i5.5,a)')it0,'_rank=',nrank,'.dat'
       call fio__input(gp,uf,np2,c,q,r,delt,delx,it0,nxs,nxe,                &
                       nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze,np,nsp, &
                       nproc,nproc_i,nproc_j,nproc_k,nrank,                  &
                       dir,file11)
       call sort__bucket(up,gp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,np2,delx)

       return
    endif

    call init__loading

    if(nrank == nroot) &
         call fio__param(nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                         np,nsp,np2,n0,                                 &
                         c,q,r,0.5*r(1)*vti**2,rtemp,fpe,fge,           &
                         ldb,delt,delx,dir,file9)

  end subroutine init__set_param


  subroutine init__loading

    integer :: i, j, k, ii, isp
    real(8) :: sd, aa, bb, cc, gamp
    real(8) :: vfunc, x0, eps=1d-40

!INITIAL VELOCITY PROFILE
    vfunc(x0) = 0.5*v0*(1.d0+tanh( (x0-ldmp)/(delx+eps) ))
    
!!--- SETTING OF FIELDS ---!
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
    do i=nxgs-2,nxge+2
       uf(1,i,j,k) = b0*cos(theta*deg2rad)
       uf(2,i,j,k) = b0*sin(theta*deg2rad)*cos(phi*deg2rad)
       uf(3,i,j,k) = b0*sin(theta*deg2rad)*sin(phi*deg2rad)
       uf(4,i,j,k) = 0.0
       uf(5,i,j,k) =  vfunc(i*delx)*uf(3,i,j,k)/c
       uf(6,i,j,k) = -vfunc(i*delx)*uf(2,i,j,k)/c
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !PARTICLE POSITION
    isp = 1
!$OMP PARALLEL DO PRIVATE(ii,j,k,aa)
    do k=nzs,nze
    do j=nys,nye
       do ii=1,np2(j,k,isp)
          up(1,ii,j,k,1) = (nxs+1)*delx+(nxe-nxs-2)*delx*ii/(np2(j,k,isp)+1)
          up(1,ii,j,k,2) = up(1,ii,j,k,1)

          call random_number(aa)
          up(2,ii,j,k,1) = dble(j)*delx+delx*aa
          up(2,ii,j,k,2) = up(2,ii,j,k,1)

          call random_number(aa)
          up(3,ii,j,k,1) = dble(k)*delx+delx*aa
          up(3,ii,j,k,2) = up(3,ii,j,k,1)
       enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/sqrt(2.)
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.)
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,k,aa,bb,cc,gamp)
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)
             
             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)

             up(4,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*(2.*bb-1)
             up(5,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*2.*sqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(6,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*2.*sqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = sqrt(1.D0+(up(4,ii,j,k,isp)**2+up(5,ii,j,k,isp)**2+up(6,ii,j,k,isp)**2)/c**2)

             call random_number(cc)

             if(up(4,ii,j,k,isp)*v0 >= 0.)then
                up(4,ii,j,k,isp) = (+up(4,ii,j,k,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*up(4,ii,j,k,isp)/gamp))then
                   up(4,ii,j,k,isp) = (-up(4,ii,j,k,isp)+v0*gamp)*gam0
                else
                   up(4,ii,j,k,isp) = (+up(4,ii,j,k,isp)+v0*gamp)*gam0
                endif
             endif

          enddo
       enddo
       enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine init__loading


  subroutine init__relocate(it)

    integer, intent(in) :: it
    integer :: dn, isp, j, k, ii, ii2 ,ii3
    real(8) :: aa, bb, cc, sd, gamp, dx

    if(nxe==nxge) return

    dx  = v0*delt*mod(it,intvl2)/delx
    dn  = int(n0*abs(dx)+0.5)
    nxe  = nxe+1

    !PARTICLE POSITION
!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,k,aa)
    do k=nzs,nze
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,k,1)+ii
          ii3 = np2(j,k,2)+ii

          up(1,ii2,j,k,1) = (nxe-1.D0+dx)*delx-dx*delx*ii/(dn+1.D0)
          up(1,ii3,j,k,2) = up(1,ii2,j,k,1)

          call random_number(aa)
          up(2,ii2,j,k,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,k,2) = up(2,ii2,j,k,1)

          call random_number(aa)
          up(3,ii2,j,k,1) = dble(k)*delx+delx*aa
          up(3,ii3,j,k,2) = up(3,ii2,j,k,1)
       enddo
       do ii=1,n0
          ii2 = np2(j,k,1)+dn+ii
          ii3 = np2(j,k,2)+dn+ii

          up(1,ii2,j,k,1) = (nxe-1.)*delx+delx*ii/(n0+1.D0)
          up(1,ii3,j,k,2) = up(1,ii2,j,k,1)

          call random_number(aa)
          up(2,ii2,j,k,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,k,2) = up(2,ii2,j,k,1)

          call random_number(aa)
          up(3,ii2,j,k,1) = dble(k)*delx+delx*aa
          up(3,ii3,j,k,2) = up(3,ii2,j,k,1)
       enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/sqrt(2.)
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.)
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,k,aa,bb,cc,gamp)
       do k=nzs,nze
       do j=nys,nye
          do ii=np2(j,k,isp)+1,np2(j,k,isp)+dn+n0

             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)
             
             up(4,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*(2.*bb-1)
             up(5,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*2.*sqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(6,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*2.*sqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = sqrt(1.D0+(up(4,ii,j,k,isp)**2+up(5,ii,j,k,isp)**2+up(6,ii,j,k,isp)**2)/c**2)

             call random_number(cc)

             if(up(4,ii,j,k,isp)*v0 >= 0.)then
                up(4,ii,j,k,isp) = (+up(4,ii,j,k,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*up(4,ii,j,k,isp)/gamp))then
                   up(4,ii,j,k,isp) = (-up(4,ii,j,k,isp)+v0*gamp)*gam0
                else
                   up(4,ii,j,k,isp) = (+up(4,ii,j,k,isp)+v0*gamp)*gam0
                endif
             endif
          enddo
          np2(j,k,isp) = np2(j,k,isp)+dn+n0
          cumcnt(nxe-1,j,k,isp) = cumcnt(nxe-1,j,k,isp)+dn
          cumcnt(nxe,j,k,isp) = cumcnt(nxe-1,j,k,isp)+n0
       enddo
       enddo
!$OMP END PARALLEL DO

    enddo

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
       uf(2,nxe-1,j,k) = b0*sin(theta*deg2rad)*cos(phi*deg2rad)
       uf(3,nxe-1,j,k) = b0*sin(theta*deg2rad)*sin(phi*deg2rad)
       uf(5,nxe-1,j,k) = +v0*uf(3,nxe-1,j,k)/c
       uf(6,nxe-1,j,k) = -v0*uf(2,nxe-1,j,k)/c

       uf(2,nxe,j,k) = b0*sin(theta*deg2rad)*cos(phi*deg2rad)
       uf(3,nxe,j,k) = b0*sin(theta*deg2rad)*sin(phi*deg2rad)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine init__relocate


  subroutine init__inject(it)

    integer, intent(in) :: it
    integer :: isp, ii, ii2, ii3, j, k, dn
    real(8) :: sd, aa, bb, cc, dx, gamp

    !INJECT PARTICLES IN x=nxe-v0*dt~nxe*dt
    dx  = v0*delt*(it-max(int(it/intvl3)*intvl3,it-intvl2))/delx
    if(dx == 0.0D0) dx = v0*delt*min(intvl2,intvl3)/delx
    if(nxe == nxge) dx = v0*delt*intvl2/delx
    dn  = int(abs(n0*dx)+0.5)

!$OMP PARALLEL DO PRIVATE(ii,ii2,ii3,j,k,aa)
    do k=nzs,nze
    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,k,1)+ii
          ii3 = np2(j,k,2)+ii

          up(1,ii2,j,k,1) = (nxe-1)*delx+dx*(dn-ii+1.D0)/(dn+1.D0)
          up(1,ii3,j,k,2) = up(1,ii2,j,k,1)

          call random_number(aa)
          up(2,ii2,j,k,1) = dble(j)*delx+delx*aa
          up(2,ii3,j,k,2) = up(2,ii2,j,k,1)

          call random_number(aa)
          up(3,ii2,j,k,1) = dble(k)*delx+delx*aa
          up(3,ii3,j,k,2) = up(3,ii2,j,k,1)
       enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/sqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.0D0)
       endif

!$OMP PARALLEL DO PRIVATE(ii,j,k,aa,bb,cc,gamp)
       do k=nzs,nze
       do j=nys,nye
          do ii=np2(j,k,isp)+1,np2(j,k,isp)+dn

             aa = 0.0D0
             do while(aa==0.D0)
                call random_number(aa)
             enddo
             call random_number(bb)
             call random_number(cc)

             up(4,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*(2.*bb-1)
             up(5,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*2.*sqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(6,ii,j,k,isp) = sd*sqrt(-2.*dlog(aa))*2.*sqrt(bb*(1.-bb))*sin(2.*pi*cc)
             gamp = sqrt(1.D0+(up(4,ii,j,k,isp)**2+up(5,ii,j,k,isp)**2+up(6,ii,j,k,isp)**2)/c**2)

             call random_number(cc)

             if(up(4,ii,j,k,isp)*v0 >= 0.)then
                up(4,ii,j,k,isp) = (+up(4,ii,j,k,isp)+v0*gamp)*gam0
             else
                if(cc < (-v0*up(4,ii,j,k,isp)/gamp))then
                   up(4,ii,j,k,isp) = (-up(4,ii,j,k,isp)+v0*gamp)*gam0
                else
                   up(4,ii,j,k,isp) = (+up(4,ii,j,k,isp)+v0*gamp)*gam0
                endif
             endif
          enddo
       enddo
       enddo
!$OMP END PARALLEL DO
    enddo

    do isp=1,nsp
!$OMP WORKSHARE
       np2(nys:nye,nzs:nze,isp) = np2(nys:nye,nzs:nze,isp)+dn
       cumcnt(nxe,nys:nye,nzs:nze,isp) = cumcnt(nxe,nys:nye,nzs:nze,isp)+dn
!$OMP END WORKSHARE
    enddo

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
       uf(2,nxe-1,j,k) = b0*sin(theta*deg2rad)*cos(phi*deg2rad)
       uf(3,nxe-1,j,k) = b0*sin(theta*deg2rad)*sin(phi*deg2rad)
       uf(5,nxe-1,j,k) = v0*uf(3,nxe-1,j,k)/c
       uf(6,nxe-1,j,k) = -v0*uf(2,nxe-1,j,k)/c

       uf(2,nxe,j,k) = b0*sin(theta*deg2rad)*cos(phi*deg2rad)
       uf(3,nxe,j,k) = b0*sin(theta*deg2rad)*sin(phi*deg2rad)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine init__inject


end module init
