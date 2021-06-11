module field

  implicit none

  private

  public :: field__fdtd_i


contains

  
  subroutine field__fdtd_i(uf,up,gp,                                        &
                           nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                           jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                           q,c,delx,delt,gfac)

    use boundary, only : boundary__dfield, boundary__curre
 
    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze
    integer, intent(in)    :: np, nsp
    integer, intent(in)    :: cumcnt(nxgs:nxge,nys:nye,nzs:nze,nsp)
    integer, intent(in)    :: jup, jdown, kup, kdown, opsum, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: q(nsp), c, delx, delt, gfac
    real(8), intent(in)    :: gp(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)    :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    logical, save              :: lflag=.true.
    integer                    :: i, j, k, ieq
    real(8), parameter         :: pi = 4.0D0*datan(1.0D0)
    real(8)                    :: f1, f2, f3
    real(8)                    :: uj(3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), save, allocatable :: df(:,:,:,:), gkl(:,:,:,:)

    if(lflag)then
       allocate(df(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
       allocate(gkl(3,nxgs:nxge,nys:nye,nzs:nze))
!$OMP PARALLEL WORKSHARE
       df(1:6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2) = 0.0D0
!$OMP END PARALLEL WORKSHARE
!$OMP PARALLEL WORKSHARE
       gkl(1:3,nxgs:nxge,nys:nye,nzs:nze) = 0.0D0
!$OMP END PARALLEL WORKSHARE
       lflag=.false.
    endif

    call ele_cur(uj,up,gp,                                        &
                 nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                 q,delx,delt)

    call boundary__curre(uj,nxs,nxe,nys,nye,nzs,nze, &
                         jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    !calculation
    f1 = c*delt/delx
    f2 = gfac*f1*f1
    f3 = 4.0*pi*delx/c
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
       gkl(1,i,j,k) = +f2*(+uf(1,i,j,k-1)+uf(1,i,j-1,k)                &
                           +uf(1,i-1,j,k)-6.*uf(1,i,j,k)+uf(1,i+1,j,k) &
                           +uf(1,i,j+1,k)+uf(1,i,j,k+1)                &
                           +f3*(-uj(3,i,j-1,k)+uj(3,i,j,k)             &
                                +uj(2,i,j,k-1)-uj(2,i,j,k))            &
                          )                                            &
                      -f1*(-uf(6,i,j-1,k)+uf(6,i,j,k)                  &
                           +uf(5,i,j,k-1)-uf(5,i,j,k))
       gkl(2,i,j,k) = +f2*(+uf(2,i,j,k-1)+uf(2,i,j-1,k)                &
                           +uf(2,i-1,j,k)-6.*uf(2,i,j,k)+uf(2,i+1,j,k) &
                           +uf(2,i,j+1,k)+uf(2,i,j,k+1)                &
                           +f3*(-uj(1,i,j,k-1)+uj(1,i,j,k)             &
                                +uj(3,i-1,j,k)-uj(3,i,j,k))            &
                          )                                            &
                      -f1*(-uf(4,i,j,k-1)+uf(4,i,j,k)                  &
                           +uf(6,i-1,j,k)-uf(6,i,j,k))
       gkl(3,i,j,k) = +f2*(+uf(3,i,j,k-1)+uf(3,i,j-1,k)                &
                           +uf(3,i-1,j,k)-6.*uf(3,i,j,k)+uf(3,i+1,j,k) &
                           +uf(3,i,j+1,k)+uf(3,i,j,k+1)                &
                           +f3*(-uj(2,i-1,j,k)+uj(2,i,j,k)             &
                                +uj(1,i,j-1,k)-uj(1,i,j,k))            &
                          )                                            &
                      -f1*(-uf(5,i-1,j,k)+uf(5,i,j,k)                  &
                           +uf(4,i,j-1,k)-uf(4,i,j,k))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !solve  < bx, by & bz >
    call cgm(df,gkl,                                          &
             nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,               &
             jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
             c,delx,delt,gfac)

    call boundary__dfield(df,                                &
                          nxgs,nxge,nxs,nxe,nys,nye,nzs,nze, &
                          jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    !solve  < ex, ey & ez >
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
       df(4,i,j,k) = +f1*(+gfac*(-df(3,i,j,k)+df(3,i,j+1,k)    &
                                 +df(2,i,j,k)-df(2,i,j,k+1))   &
                          +     (-uf(3,i,j,k)+uf(3,i,j+1,k)    &
                                 +uf(2,i,j,k)-uf(2,i,j,k+1)) ) &
                     -4.*pi*delt*uj(1,i,j,k)

       df(5,i,j,k) = +f1*(+gfac*(-df(1,i,j,k)+df(1,i,j,k+1)    &
                                 +df(3,i,j,k)-df(3,i+1,j,k))   &
                          +     (-uf(1,i,j,k)+uf(1,i,j,k+1)    &
                                 +uf(3,i,j,k)-uf(3,i+1,j,k)) ) &
                     -4.*pi*delt*uj(2,i,j,k)

       df(6,i,j,k) = +f1*(+gfac*(-df(2,i,j,k)+df(2,i+1,j,k)    &
                                 +df(1,i,j,k)-df(1,i,j+1,k))   &
                          +     (-uf(2,i,j,k)+uf(2,i+1,j,k)    &
                                 +uf(1,i,j,k)-uf(1,i,j+1,k)) ) &
                     -4.*pi*delt*uj(3,i,j,k)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    call boundary__dfield(df,                                &
                          nxgs,nxge,nxs,nxe,nys,nye,nzs,nze, &
                          jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    !===== Update fields and particles ======
!$OMP PARALLEL DO PRIVATE(i,j,k,ieq)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       do ieq=1,6
          uf(ieq,i,j,k) = uf(ieq,i,j,k)+df(ieq,i,j,k)
       enddo
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine field__fdtd_i


  subroutine ele_cur(uj,up,gp,                                        &
                     nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                     q,delx,delt)

    integer, intent(in)  :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze, np, nsp 
    integer, intent(in)  :: cumcnt(nxgs:nxge,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: q(nsp), delx, delt
    real(8), intent(in)  :: gp(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(out) :: uj(3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2)

    integer            :: ii, i, j, k, isp, i2, inc, ip, jp, kp
    real(8), parameter :: fac = 1.D0/3.D0
    real(8)            :: idelx, idelt, dh
    real(8)            :: s1_1, s1_2, s1_3, smo_1, smo_2, smo_3
    real(8)            :: pjtmpx, pjtmpy, pjtmpz, dstmpx, dstmpy, dstmpz
    real(8)            :: s0x(-2:2), s0y(-2:2), s0z(-2:2), dsx(-2:2), dsy(-2:2), dsz(-2:2)
    real(8)            :: pjx(-2:2,-2:2,-2:2), pjy(-2:2,-2:2,-2:2), pjz(-2:2,-2:2,-2:2)

!$OMP PARALLEL WORKSHARE
    uj(1:3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2) = 0.D0
!$OMP END PARALLEL WORKSHARE

    idelt = 1.D0/delt
    idelx = 1.D0/delx

    !--------------Charge Conservation Method -------------!
    !---- Density Decomposition (Esirkepov, CPC, 2001) ----!
!$OMP PARALLEL DO PRIVATE(ii,i,j,k,isp,i2,inc,ip,jp,kp,              &
!$OMP                     s1_1,s1_2,s1_3,smo_1,smo_2,smo_3,dh,       &
!$OMP                     s0x,s0y,s0z,dsx,dsy,dsz,pjx,pjy,pjz,       &
!$OMP                     pjtmpx,pjtmpy,pjtmpz,dstmpx,dstmpy,dstmpz) &
!$OMP REDUCTION(+:uj)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe-1

!OCL UNROLL('FULL')
       pjx(-2:2,-2:2,-2:2) = 0.D0
!OCL UNROLL('FULL')
       pjy(-2:2,-2:2,-2:2) = 0.D0
!OCL UNROLL('FULL')
       pjz(-2:2,-2:2,-2:2) = 0.D0
 
       isp=1

       do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)
	
!          second order shape function
          dh = up(1,ii,j,k,isp)*idelx-0.5-i
          s0x(-2) = 0.D0
          s0x(-1) = 0.5*(0.5-dh)*(0.5-dh)
          s0x( 0) = 0.75-dh*dh
          s0x(+1) = 0.5*(0.5+dh)*(0.5+dh)
          s0x(+2) = 0.D0

          dh = up(2,ii,j,k,isp)*idelx-0.5-j
          s0y(-2) = 0.D0
          s0y(-1) = 0.5*(0.5-dh)*(0.5-dh)
          s0y( 0) = 0.75-dh*dh
          s0y(+1) = 0.5*(0.5+dh)*(0.5+dh)
          s0y(+2) = 0.D0

          dh = up(3,ii,j,k,isp)*idelx-0.5-k
          s0z(-2) = 0.D0
          s0z(-1) = 0.5*(0.5-dh)*(0.5-dh)
          s0z( 0) = 0.75-dh*dh
          s0z(+1) = 0.5*(0.5+dh)*(0.5+dh)
          s0z(+2) = 0.D0

          i2 = int(gp(1,ii,j,k,isp)*idelx)
          dh = gp(1,ii,j,k,isp)*idelx-0.5-i2
          inc = i2-i
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          dsx(-2) = s1_1*smo_1
          dsx(-1) = s1_1*smo_2+s1_2*smo_1
          dsx( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsx(+1) = s1_3*smo_2+s1_2*smo_3
          dsx(+2) = s1_3*smo_3

          i2 = int(gp(2,ii,j,k,isp)*idelx)
          dh = gp(2,ii,j,k,isp)*idelx-0.5-i2
          inc = i2-j
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          dsy(-2) = s1_1*smo_1
          dsy(-1) = s1_1*smo_2+s1_2*smo_1
          dsy( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsy(+1) = s1_3*smo_2+s1_2*smo_3
          dsy(+2) = s1_3*smo_3

          i2 = int(gp(3,ii,j,k,isp)*idelx)
          dh = gp(3,ii,j,k,isp)*idelx-0.5-i2
          inc = i2-k
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          dsz(-2) = s1_1*smo_1
          dsz(-1) = s1_1*smo_2+s1_2*smo_1
          dsz( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsz(+1) = s1_3*smo_2+s1_2*smo_3
          dsz(+2) = s1_3*smo_3

          dsx(-2:2) = dsx(-2:2)-s0x(-2:2)
          dsy(-2:2) = dsy(-2:2)-s0y(-2:2)
          dsz(-2:2) = dsz(-2:2)-s0z(-2:2)

!OCL UNROLL('FULL')
          do kp=-2,2
          do jp=-2,2
             pjtmpx = 0.D0
             pjtmpy = 0.D0
             pjtmpz = 0.D0
             dstmpx =  (s0y(jp)+0.5*dsy(jp))*s0z(kp) &
                      +(0.5*s0y(jp)+fac*dsy(jp))*dsz(kp)
             dstmpy =  (s0x(jp)+0.5*dsx(jp))*s0z(kp) &
                      +(0.5*s0x(jp)+fac*dsx(jp))*dsz(kp)
             dstmpz =  (s0x(jp)+0.5*dsx(jp))*s0y(kp) &
                      +(0.5*s0x(jp)+fac*dsx(jp))*dsy(kp)
          do ip=-2,1
             pjtmpx = pjtmpx-q(isp)*delx*idelt*dsx(ip)*dstmpx
             pjtmpy = pjtmpy-q(isp)*delx*idelt*dsy(ip)*dstmpy
             pjtmpz = pjtmpz-q(isp)*delx*idelt*dsz(ip)*dstmpz
             pjx(ip+1,jp,kp) = pjx(ip+1,jp,kp)+pjtmpx
             pjy(ip+1,jp,kp) = pjy(ip+1,jp,kp)+pjtmpy
             pjz(ip+1,jp,kp) = pjz(ip+1,jp,kp)+pjtmpz
          enddo
          enddo
          enddo

       enddo

       isp=2

       do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)
	
!          second order shape function
          dh = up(1,ii,j,k,isp)*idelx-0.5-i
          s0x(-2) = 0.D0
          s0x(-1) = 0.5*(0.5-dh)*(0.5-dh)
          s0x( 0) = 0.75-dh*dh
          s0x(+1) = 0.5*(0.5+dh)*(0.5+dh)
          s0x(+2) = 0.D0

          dh = up(2,ii,j,k,isp)*idelx-0.5-j
          s0y(-2) = 0.D0
          s0y(-1) = 0.5*(0.5-dh)*(0.5-dh)
          s0y( 0) = 0.75-dh*dh
          s0y(+1) = 0.5*(0.5+dh)*(0.5+dh)
          s0y(+2) = 0.D0

          dh = up(3,ii,j,k,isp)*idelx-0.5-k
          s0z(-2) = 0.D0
          s0z(-1) = 0.5*(0.5-dh)*(0.5-dh)
          s0z( 0) = 0.75-dh*dh
          s0z(+1) = 0.5*(0.5+dh)*(0.5+dh)
          s0z(+2) = 0.D0

          i2 = int(gp(1,ii,j,k,isp)*idelx)
          dh = gp(1,ii,j,k,isp)*idelx-0.5-i2
          inc = i2-i
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          dsx(-2) = s1_1*smo_1
          dsx(-1) = s1_1*smo_2+s1_2*smo_1
          dsx( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsx(+1) = s1_3*smo_2+s1_2*smo_3
          dsx(+2) = s1_3*smo_3

          i2 = int(gp(2,ii,j,k,isp)*idelx)
          dh = gp(2,ii,j,k,isp)*idelx-0.5-i2
          inc = i2-j
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          dsy(-2) = s1_1*smo_1
          dsy(-1) = s1_1*smo_2+s1_2*smo_1
          dsy( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsy(+1) = s1_3*smo_2+s1_2*smo_3
          dsy(+2) = s1_3*smo_3

          i2 = int(gp(3,ii,j,k,isp)*idelx)
          dh = gp(3,ii,j,k,isp)*idelx-0.5-i2
          inc = i2-k
          s1_1 = 0.5*(0.5-dh)*(0.5-dh)
          s1_2 = 0.75-dh*dh
          s1_3 = 0.5*(0.5+dh)*(0.5+dh)
          smo_1 = -(inc-abs(inc))*0.5+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*0.5+0
          dsz(-2) = s1_1*smo_1
          dsz(-1) = s1_1*smo_2+s1_2*smo_1
          dsz( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsz(+1) = s1_3*smo_2+s1_2*smo_3
          dsz(+2) = s1_3*smo_3

          dsx(-2:2) = dsx(-2:2)-s0x(-2:2)
          dsy(-2:2) = dsy(-2:2)-s0y(-2:2)
          dsz(-2:2) = dsz(-2:2)-s0z(-2:2)

!OCL UNROLL('FULL')
          do kp=-2,2
          do jp=-2,2
             pjtmpx = 0.D0
             pjtmpy = 0.D0
             pjtmpz = 0.D0
             dstmpx =  (s0y(jp)+0.5*dsy(jp))*s0z(kp) &
                      +(0.5*s0y(jp)+fac*dsy(jp))*dsz(kp)
             dstmpy =  (s0x(jp)+0.5*dsx(jp))*s0z(kp) &
                      +(0.5*s0x(jp)+fac*dsx(jp))*dsz(kp)
             dstmpz =  (s0x(jp)+0.5*dsx(jp))*s0y(kp) &
                      +(0.5*s0x(jp)+fac*dsx(jp))*dsy(kp)
          do ip=-2,1
             pjtmpx = pjtmpx-q(isp)*delx*idelt*dsx(ip)*dstmpx
             pjtmpy = pjtmpy-q(isp)*delx*idelt*dsy(ip)*dstmpy
             pjtmpz = pjtmpz-q(isp)*delx*idelt*dsz(ip)*dstmpz
             pjx(ip+1,jp,kp) = pjx(ip+1,jp,kp)+pjtmpx
             pjy(ip+1,jp,kp) = pjy(ip+1,jp,kp)+pjtmpy
             pjz(ip+1,jp,kp) = pjz(ip+1,jp,kp)+pjtmpz
          enddo
          enddo
          enddo

       enddo

       do kp=-2,2
       do jp=-2,2
       do ip=-2,2
          uj(1,i+ip,j+jp,k+kp) = uj(1,i+ip,j+jp,k+kp)+pjx(ip,jp,kp)
          uj(2,i+ip,j+jp,k+kp) = uj(2,i+ip,j+jp,k+kp)+pjy(jp,ip,kp)
          uj(3,i+ip,j+jp,k+kp) = uj(3,i+ip,j+jp,k+kp)+pjz(kp,ip,jp)
       enddo
       enddo
       enddo

    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine ele_cur


  subroutine cgm(db,gkl,                                          &
                 nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,               &
                 jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                 c,delx,delt,gfac)

    use boundary, only : boundary__phi

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !  #  this routine will be stopped after iteration number = ite_max
    !-----------------------------------------------------------------------

    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze
    integer, intent(in)    :: jup, jdown, kup, kdown, mnpr, opsum, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: c, delx, delt, gfac
    real(8), intent(in)    :: gkl(3,nxgs:nxge,nys:nye,nzs:nze)
    real(8), intent(inout) :: db(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)

    integer, parameter :: ite_max = 100 ! maximum number of interation
    integer            :: i, j, k, l, ite, bc
    real(8), parameter :: err = 1d-6 
    real(8)            :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g
    real(8)            :: phi(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: p(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: b(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: r(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: ap(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: bff_snd(2),bff_rcv(2)

    do l=1,3

       select case(l)
         case(1)
           bc = -1
         case(2,3)
           bc = 0
       endselect

       ! initial guess
       ite = 0
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sum)
       do k=nzs,nze
       do j=nys,nye
       do i=nxs,nxe+bc
          phi(i,j,k) = db(l,i,j,k)
          b(i,j,k) = f2*gkl(l,i,j,k)
          sum = sum+b(i,j,k)*b(i,j,k)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !****** boundary conditions of PHI *****!
       call boundary__phi(phi,                       &
                          nxs,nxe,nys,nye,nzs,nze,l, &
                          jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)
       !****** end of boundary conditions *****!

       f1 = 6.0+(delx/(c*delt*gfac))**2
       sumr = 0.0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sumr)
       do k=nzs,nze
       do j=nys,nye
       do i=nxs,nxe+bc
          r(i,j,k) = b(i,j,k)+phi(i,j,k-1)+phi(i,j-1,k)             &
                             +phi(i-1,j,k)-f1*phi(i,j,k)+phi(i+1,j,k) &
                             +phi(i,j+1,k)+phi(i,j,k+1)
          p(i,j,k) = r(i,j,k)
          sumr = sumr+r(i,j,k)*r(i,j,k)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(dsqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !****** boundary conditions of P *****!
             call boundary__phi(p,                        & 
                                nxs,nxe,nys,nye,nzs,nze,l,&
                                jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)
             !****** end of boundary conditions *****!
       
             sumr = 0.0
             sum2 = 0.0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sumr,sum2)
             do k=nzs,nze
             do j=nys,nye
             do i=nxs,nxe+bc
                ap(i,j,k) = -p(i,j,k-1)-p(i,j-1,k)             &
                            -p(i-1,j,k)+f1*p(i,j,k)-p(i+1,j,k) &
                            -p(i,j+1,k)-p(i,j,k+1)
                sumr = sumr+r(i,j,k)*r(i,j,k)
                sum2 = sum2+p(i,j,k)*ap(i,j,k)
             enddo
             enddo
             enddo
!$OMP END PARALLEL DO

             bff_snd(1) = sumr
             bff_snd(2) = sum2
             call MPI_ALLREDUCE(bff_snd,bff_rcv,2,mnpr,opsum,ncomw,nerr)
             sumr_g = bff_rcv(1)
             sum2_g = bff_rcv(2)

             av = sumr_g/sum2_g

!$OMP PARALLEL DO PRIVATE(i,j,k)
             do k=nzs,nze
             do j=nys,nye
             do i=nxs,nxe+bc
                phi(i,j,k) = phi(i,j,k)+av*p(i,j,k)
                r(i,j,k) = r(i,j,k)-av*ap(i,j,k)
             enddo
             enddo
             enddo
!$OMP END PARALLEL DO
             
             sum_g = dsqrt(sumr_g)

             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sum1)
             do k=nzs,nze
             do j=nys,nye
             do i=nxs,nxe+bc
                sum1 = sum1+r(i,j,k)*r(i,j,k)
             enddo
             enddo
             enddo
!$OMP END PARALLEL DO

             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g
             
!$OMP PARALLEL DO PRIVATE(i,j,k)
             do k=nzs,nze
             do j=nys,nye
             do i=nxs,nxe+bc
                p(i,j,k) = r(i,j,k)+bv*p(i,j,k)
             enddo
             enddo
             enddo
!$OMP END PARALLEL DO
             
          enddo
       endif

!$OMP PARALLEL WORKSHARE
       db(l,nxs:nxe+bc,nys:nye,nzs:nze) = phi(nxs:nxe+bc,nys:nye,nzs:nze)
!$OMP END PARALLEL WORKSHARE

    end do
    
  end subroutine cgm


end module field
