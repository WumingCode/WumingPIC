module field

  implicit none

  private

  public :: field__fdtd_i


contains

  
  subroutine field__fdtd_i(uf,up,gp,                                        &
                           nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                           jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                           q,c,delx,delt,gfac)

    use boundary, only : boundary__field, boundary__curre
 
    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze
    integer, intent(in)    :: np, nsp
    integer, intent(in)    :: cumcnt(nxgs:nxge,nys:nye,nzs:nze,nsp)
    integer, intent(in)    :: jup, jdown, kup, kdown, opsum, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: q(nsp), c, delx, delt, gfac
    real(8), intent(in)    :: gp(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    logical, save              :: lflag=.true.
    integer                    :: i, j, k, ieq
    real(8)                    :: pi, f1, f2, f3
    real(8)                    :: uj(3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2)
    real(8)                    :: gkl(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), save, allocatable :: gf(:,:,:,:)

    pi = 4.0*atan(1.0)

    if(lflag)then
       allocate(gf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
!$OMP PARALLEL WORKSHARE
       gf(1:6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2) = 0.0D0
!$OMP END PARALLEL WORKSHARE
       lflag=.false.
    endif

    call ele_cur(uj,up,gp,                                        &
                 nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                 q,delx,delt)

    call boundary__curre(uj,nxs,nxe,nys,nye,nzs,nze, &
                         jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    !calculation
    !< gkl(1:3) =  (c*delt)*rot(e) >
    !< gkl(4:6) =  (c*delt)*rot(b) - (4*pi*delt)*j >
    f1 = c*delt/delx
    f2 = 4.0*pi*delt
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
       gkl(1,i,j,k) = -f1*(+(-uf(6,i,j-1,k)+uf(6,i,j,k))-(-uf(5,i,j,k-1)+uf(5,i,j,k)))
       gkl(2,i,j,k) = -f1*(+(-uf(4,i,j,k-1)+uf(4,i,j,k))-(-uf(6,i-1,j,k)+uf(6,i,j,k)))
       gkl(3,i,j,k) = -f1*(+(-uf(5,i-1,j,k)+uf(5,i,j,k))-(-uf(4,i,j-1,k)+uf(4,i,j,k)))
       gkl(4,i,j,k) = +f1*(+(-uf(3,i,j,k)+uf(3,i,j+1,k))-(-uf(2,i,j,k)+uf(2,i,j,k+1)))-f2*uj(1,i,j,k)
       gkl(5,i,j,k) = +f1*(+(-uf(1,i,j,k)+uf(1,i,j,k+1))-(-uf(3,i,j,k)+uf(3,i+1,j,k)))-f2*uj(2,i,j,k)
       gkl(6,i,j,k) = +f1*(+(-uf(2,i,j,k)+uf(2,i+1,j,k))-(-uf(1,i,j,k)+uf(1,i,j+1,k)))-f2*uj(3,i,j,k)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    call boundary__field(gkl,                               &
                         nxgs,nxge,nxs,nxe,nys,nye,nzs,nze, &
                         jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    f3 = c*delt*gfac/delx
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
       gkl(1,i,j,k) = gkl(1,i,j,k)-f3*(+(-gkl(6,i,j-1,k)+gkl(6,i,j,k))-(-gkl(5,i,j,k-1)+gkl(5,i,j,k)))
       gkl(2,i,j,k) = gkl(2,i,j,k)-f3*(+(-gkl(4,i,j,k-1)+gkl(4,i,j,k))-(-gkl(6,i-1,j,k)+gkl(6,i,j,k)))
       gkl(3,i,j,k) = gkl(3,i,j,k)-f3*(+(-gkl(5,i-1,j,k)+gkl(5,i,j,k))-(-gkl(4,i,j-1,k)+gkl(4,i,j,k)))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !solve  < bx, by & bz >
    call cgm(gf,gkl,                                          &
             nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,               &
             jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
             c,delx,delt,gfac)
    call boundary__field(gf,                                &
                         nxgs,nxge,nxs,nxe,nys,nye,nzs,nze, &
                         jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    !solve  < ex, ey & ez >
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
       gf(4,i,j,k) = gkl(4,i,j,k)+f3*(+(-gf(3,i,j,k)+gf(3,i,j+1,k))-(-gf(2,i,j,k)+gf(2,i,j,k+1)))
       gf(5,i,j,k) = gkl(5,i,j,k)+f3*(+(-gf(1,i,j,k)+gf(1,i,j,k+1))-(-gf(3,i,j,k)+gf(3,i+1,j,k)))
       gf(6,i,j,k) = gkl(6,i,j,k)+f3*(+(-gf(2,i,j,k)+gf(2,i+1,j,k))-(-gf(1,i,j,k)+gf(1,i,j+1,k)))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    call boundary__field(gf,                                &
                         nxgs,nxge,nxs,nxe,nys,nye,nzs,nze, &
                         jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    !===== Update fields and particles ======
    
!$OMP PARALLEL DO PRIVATE(i,j,k,ieq)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       do ieq=1,6
          uf(ieq,i,j,k) = uf(ieq,i,j,k)+gf(ieq,i,j,k)
       enddo
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine field__fdtd_i


  subroutine ele_cur(uj,up,gp,                                        &
                     nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                     q,delx,delt)

    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze, np, nsp 
    integer, intent(in)    :: cumcnt(nxgs:nxge,nys:nye,nzs:nze,nsp)
    real(8), intent(in)    :: q(nsp), delx, delt
    real(8), intent(in)    :: gp(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(out)   :: uj(3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2)

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

       pjx(-2:2,-2:2,-2:2) = 0.D0
       pjy(-2:2,-2:2,-2:2) = 0.D0
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

          up(1,ii,j,k,isp) = gp(1,ii,j,k,isp)
          up(2,ii,j,k,isp) = gp(2,ii,j,k,isp)
          up(3,ii,j,k,isp) = gp(3,ii,j,k,isp)
          up(4,ii,j,k,isp) = gp(4,ii,j,k,isp)
          up(5,ii,j,k,isp) = gp(5,ii,j,k,isp)
          up(6,ii,j,k,isp) = gp(6,ii,j,k,isp)

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

          up(1,ii,j,k,isp) = gp(1,ii,j,k,isp)
          up(2,ii,j,k,isp) = gp(2,ii,j,k,isp)
          up(3,ii,j,k,isp) = gp(3,ii,j,k,isp)
          up(4,ii,j,k,isp) = gp(4,ii,j,k,isp)
          up(5,ii,j,k,isp) = gp(5,ii,j,k,isp)
          up(6,ii,j,k,isp) = gp(6,ii,j,k,isp)

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


  subroutine cgm(gb,gkl,                                          &
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
    real(8), intent(in)    :: gkl(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(inout) :: gb(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)

    integer, parameter :: ite_max = 100 ! maximum number of interation
    integer            :: i, j, k, l, ite, bc
    real(8), parameter :: err = 1d-6 
    real(8)            :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g
    real(8)            :: x(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1), b(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: r(nxs:nxe,nys:nye,nzs:nze), p(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: ap(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: bff_snd(2),bff_rcv(2)

    do l=1,3

       if(l == 1)then
          bc=-1
       else
          bc=0
       endif

       ! initial guess
       ite = 0
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sum)
       do k=nzs,nze
       do j=nys,nye
       do i=nxs,nxe+bc
          x(i,j,k) = gb(l,i,j,k)
          b(i,j,k) = f2*gkl(l,i,j,k)
          sum = sum+b(i,j,k)*b(i,j,k)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !****** boundary conditions of X *****!
       call boundary__phi(x,                       &
                          nxs,nxe,nys,nye,nzs,nze, &
                          jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

       select case(l)
       case(1)

!$OMP PARALLEL DO PRIVATE(j,k)
          do k=nzs-1,nze+1
          do j=nys-1,nye+1
             x(nxs-1,j,k) = -x(nxs  ,j,k)
             x(nxe  ,j,k) = -x(nxe-1,j,k)
          enddo
          enddo
!$OMP END PARALLEL DO

       case(2,3)

!$OMP PARALLEL DO PRIVATE(j,k)
          do k=nzs-1,nze+1
          do j=nys-1,nye+1
             x(nxs-1,j,k) = x(nxs+1,j,k)
             x(nxe+1,j,k) = x(nxe-1,j,k)
          enddo
          enddo
!$OMP END PARALLEL DO

       end select
       !****** end of boundary conditions *****!

       f1 = 6.0+(delx/(c*delt*gfac))**2
       sumr = 0.0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sumr)
       do k=nzs,nze
       do j=nys,nye
       do i=nxs,nxe+bc
          r(i,j,k) = b(i,j,k)+x(i,j,k-1)+x(i,j-1,k)             &
                             +x(i-1,j,k)-f1*x(i,j,k)+x(i+1,j,k) &
                             +x(i,j+1,k)+x(i,j,k+1)
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
             call boundary__phi(p,                       &
                                nxs,nxe,nys,nye,nzs,nze, &
                                jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

             select case(l)
             case(1)

!$OMP PARALLEL DO PRIVATE(j,k)
                do k=nzs-1,nze+1
                do j=nys-1,nye+1
                   p(nxs-1,j,k) = -p(nxs  ,j,k)
                   p(nxe  ,j,k) = -p(nxe-1,j,k)
                enddo
                enddo
!$OMP END PARALLEL DO

             case(2,3)

!$OMP PARALLEL DO PRIVATE(j,k)
                do k=nzs-1,nze+1
                do j=nys-1,nye+1
                   p(nxs-1,j,k) = p(nxs+1,j,k)
                   p(nxe+1,j,k) = p(nxe-1,j,k)
                enddo
                enddo
!$OMP END PARALLEL DO
             end select
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
                x(i,j,k) = x(i,j,k)+av* p(i,j,k)
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
       gb(l,nxs:nxe+bc,nys:nye,nzs:nze) = x(nxs:nxe+bc,nys:nye,nzs:nze)
!$OMP END PARALLEL WORKSHARE

    end do
    
  end subroutine cgm


end module field
