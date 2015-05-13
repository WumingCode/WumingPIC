module mom_calc

  implicit none

  private

  public :: mom_calc__nvt, mom_calc__accl


contains

  
  subroutine mom_calc__accl(gp,                                              &
                            nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                            c,q,r,delt,delx,                                 &
                            up,uf)

    integer, intent(in)  :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze, np, nsp
    integer, intent(in)  :: cumcnt(nxgs:nxge,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: c, q(nsp), r(nsp), delt, delx
    real(8), intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(out) :: gp(6,np,nys:nye,nzs:nze,nsp)

    integer            :: i,j, k, ii, isp, i0, j0, k0
    real(8)            :: tmpf(6,-1:2,-1:2,-1:2)
    real(8)            :: idelx, dh, fac1, fac1r, fac2, fac2r, gam, igam, txxx
    real(8)            :: bpx, bpy, bpz, epx, epy, epz
    real(8)            :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8)            :: s0xm, s0x, s0xp, s0ym, s0y, s0yp, s0zm, s0z, s0zp
    real(8)            :: shxm, shx, shxp, shym, shy, shyp, shzm, shz, shzp

    idelx = 1./delx

!$OMP PARALLEL DO PRIVATE(ii,i,j,k,isp,i0,j0,k0,                       &
!$OMP                     dh,gam,igam,fac1,fac2,txxx,fac1r,fac2r,tmpf, &
!$OMP                     s0xm,s0x,s0xp,s0ym,s0y,s0yp,s0zm,s0z,s0zp,   &
!$OMP                     shxm,shx,shxp,shym,shy,shyp,shzm,shz,shzp,   &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,                     &
!$OMP                     uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe-1

       isp = 1

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

       tmpf(1:6,-1:2,-1:2,-1:2) = uf(1:6,i-1:i+2,j-1:j+2,k-1:k+2)

       do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)

          !second order shape function
          i0 = int(up(1,ii,j,k,isp)*idelx+0.5)
          dh = up(1,ii,j,k,isp)*idelx-i0
          s0xm = 0.5*(0.5-dh)*(0.5-dh)
          s0x  = 0.75-dh*dh
          s0xp = 0.5*(0.5+dh)*(0.5+dh)

          j0 = int(up(2,ii,j,k,isp)*idelx+0.5)
          dh = up(2,ii,j,k,isp)*idelx-j0
          s0ym = 0.5*(0.5-dh)*(0.5-dh)
          s0y  = 0.75-dh*dh
          s0yp = 0.5*(0.5+dh)*(0.5+dh)

          k0 = int(up(3,ii,j,k,isp)*idelx+0.5)
          dh = up(3,ii,j,k,isp)*idelx-k0
          s0zm = 0.5*(0.5-dh)*(0.5-dh)
          s0z  = 0.75-dh*dh
          s0zp = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(1,ii,j,k,isp)*idelx-0.5-i
          shxm = 0.5*(0.5-dh)*(0.5-dh)
          shx  = 0.75-dh*dh
          shxp = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(2,ii,j,k,isp)*idelx-0.5-j
          shym = 0.5*(0.5-dh)*(0.5-dh)
          shy  = 0.75-dh*dh
          shyp = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(3,ii,j,k,isp)*idelx-0.5-k
          shzm = 0.5*(0.5-dh)*(0.5-dh)
          shz  = 0.75-dh*dh
          shzp = 0.5*(0.5+dh)*(0.5+dh)

          bpx = (+(+tmpf(1,-1,j0-j-1,k0-k-1)*shxm+tmpf(1,0,j0-j-1,k0-k-1)*shx+tmpf(1,+1,j0-j-1,k0-k-1)*shxp)*s0ym       &
                 +(+tmpf(1,-1,j0-j  ,k0-k-1)*shxm+tmpf(1,0,j0-j  ,k0-k-1)*shx+tmpf(1,+1,j0-j  ,k0-k-1)*shxp)*s0y        &
                 +(+tmpf(1,-1,j0-j+1,k0-k-1)*shxm+tmpf(1,0,j0-j+1,k0-k-1)*shx+tmpf(1,+1,j0-j+1,k0-k-1)*shxp)*s0yp)*s0zm &
               +(+(+tmpf(1,-1,j0-j-1,k0-k  )*shxm+tmpf(1,0,j0-j-1,k0-k  )*shx+tmpf(1,+1,j0-j-1,k0-k  )*shxp)*s0ym       &
                 +(+tmpf(1,-1,j0-j  ,k0-k  )*shxm+tmpf(1,0,j0-j  ,k0-k  )*shx+tmpf(1,+1,j0-j  ,k0-k  )*shxp)*s0y        &
                 +(+tmpf(1,-1,j0-j+1,k0-k  )*shxm+tmpf(1,0,j0-j+1,k0-k  )*shx+tmpf(1,+1,j0-j+1,k0-k  )*shxp)*s0yp)*s0z  &
               +(+(+tmpf(1,-1,j0-j-1,k0-k+1)*shxm+tmpf(1,0,j0-j-1,k0-k+1)*shx+tmpf(1,+1,j0-j-1,k0-k+1)*shxp)*s0ym       &
                 +(+tmpf(1,-1,j0-j  ,k0-k+1)*shxm+tmpf(1,0,j0-j  ,k0-k+1)*shx+tmpf(1,+1,j0-j  ,k0-k+1)*shxp)*s0y        &
                 +(+tmpf(1,-1,j0-j+1,k0-k+1)*shxm+tmpf(1,0,j0-j+1,k0-k+1)*shx+tmpf(1,+1,j0-j+1,k0-k+1)*shxp)*s0yp)*s0zp

          bpy = (+(+tmpf(2,i0-i-1,-1,k0-k-1)*s0xm+tmpf(2,i0-i,-1,k0-k-1)*s0x+tmpf(2,i0-i+1,-1,k0-k-1)*s0xp)*shym       &
                 +(+tmpf(2,i0-i-1, 0,k0-k-1)*s0xm+tmpf(2,i0-i, 0,k0-k-1)*s0x+tmpf(2,i0-i+1, 0,k0-k-1)*s0xp)*shy        &
                 +(+tmpf(2,i0-i-1,+1,k0-k-1)*s0xm+tmpf(2,i0-i,+1,k0-k-1)*s0x+tmpf(2,i0-i+1,+1,k0-k-1)*s0xp)*shyp)*s0zm &
               +(+(+tmpf(2,i0-i-1,-1,k0-k  )*s0xm+tmpf(2,i0-i,-1,k0-k  )*s0x+tmpf(2,i0-i+1,-1,k0-k  )*s0xp)*shym       &
                 +(+tmpf(2,i0-i-1, 0,k0-k  )*s0xm+tmpf(2,i0-i, 0,k0-k  )*s0x+tmpf(2,i0-i+1, 0,k0-k  )*s0xp)*shy        &
                 +(+tmpf(2,i0-i-1,+1,k0-k  )*s0xm+tmpf(2,i0-i,+1,k0-k  )*s0x+tmpf(2,i0-i+1,+1,k0-k  )*s0xp)*shyp)*s0z  &
               +(+(+tmpf(2,i0-i-1,-1,k0-k+1)*s0xm+tmpf(2,i0-i,-1,k0-k+1)*s0x+tmpf(2,i0-i+1,-1,k0-k+1)*s0xp)*shym       &
                 +(+tmpf(2,i0-i-1, 0,k0-k+1)*s0xm+tmpf(2,i0-i, 0,k0-k+1)*s0x+tmpf(2,i0-i+1, 0,k0-k+1)*s0xp)*shy        &
                 +(+tmpf(2,i0-i-1,+1,k0-k+1)*s0xm+tmpf(2,i0-i,+1,k0-k+1)*s0x+tmpf(2,i0-i+1,+1,k0-k+1)*s0xp)*shyp)*s0zp

          bpz = (+(+tmpf(3,i0-i-1,j0-j-1,-1)*s0xm+tmpf(3,i0-i,j0-j-1,-1)*s0x+tmpf(3,i0-i+1,j0-j-1,-1)*s0xp)*s0ym       &
                 +(+tmpf(3,i0-i-1,j0-j  ,-1)*s0xm+tmpf(3,i0-i,j0-j  ,-1)*s0x+tmpf(3,i0-i+1,j0-j  ,-1)*s0xp)*s0y        &
                 +(+tmpf(3,i0-i-1,j0-j+1,-1)*s0xm+tmpf(3,i0-i,j0-j+1,-1)*s0x+tmpf(3,i0-i+1,j0-j+1,-1)*s0xp)*s0yp)*shzm &
               +(+(+tmpf(3,i0-i-1,j0-j-1, 0)*s0xm+tmpf(3,i0-i,j0-j-1, 0)*s0x+tmpf(3,i0-i+1,j0-j-1, 0)*s0xp)*s0ym       &
                 +(+tmpf(3,i0-i-1,j0-j  , 0)*s0xm+tmpf(3,i0-i,j0-j  , 0)*s0x+tmpf(3,i0-i+1,j0-j  , 0)*s0xp)*s0y        &
                 +(+tmpf(3,i0-i-1,j0-j+1, 0)*s0xm+tmpf(3,i0-i,j0-j+1, 0)*s0x+tmpf(3,i0-i+1,j0-j+1, 0)*s0xp)*s0yp)*shz  &
               +(+(+tmpf(3,i0-i-1,j0-j-1,+1)*s0xm+tmpf(3,i0-i,j0-j-1,+1)*s0x+tmpf(3,i0-i+1,j0-j-1,+1)*s0xp)*s0ym       &
                 +(+tmpf(3,i0-i-1,j0-j  ,+1)*s0xm+tmpf(3,i0-i,j0-j  ,+1)*s0x+tmpf(3,i0-i+1,j0-j  ,+1)*s0xp)*s0y        &
                 +(+tmpf(3,i0-i-1,j0-j+1,+1)*s0xm+tmpf(3,i0-i,j0-j+1,+1)*s0x+tmpf(3,i0-i+1,j0-j+1,+1)*s0xp)*s0yp)*shzp

          epx = (+(+tmpf(4,i0-i-1,-1,-1)*s0xm+tmpf(4,i0-i,-1,-1)*s0x+tmpf(4,i0-i+1,-1,-1)*s0xp)*shym       &
                 +(+tmpf(4,i0-i-1, 0,-1)*s0xm+tmpf(4,i0-i, 0,-1)*s0x+tmpf(4,i0-i+1, 0,-1)*s0xp)*shy        &
                 +(+tmpf(4,i0-i-1,+1,-1)*s0xm+tmpf(4,i0-i,+1,-1)*s0x+tmpf(4,i0-i+1,+1,-1)*s0xp)*shyp)*shzm &
               +(+(+tmpf(4,i0-i-1,-1, 0)*s0xm+tmpf(4,i0-i,-1, 0)*s0x+tmpf(4,i0-i+1,-1, 0)*s0xp)*shym       &
                 +(+tmpf(4,i0-i-1, 0, 0)*s0xm+tmpf(4,i0-i, 0, 0)*s0x+tmpf(4,i0-i+1, 0, 0)*s0xp)*shy        &
                 +(+tmpf(4,i0-i-1,+1, 0)*s0xm+tmpf(4,i0-i,+1, 0)*s0x+tmpf(4,i0-i+1,+1, 0)*s0xp)*shyp)*shz  &
               +(+(+tmpf(4,i0-i-1,-1,+1)*s0xm+tmpf(4,i0-i,-1,+1)*s0x+tmpf(4,i0-i+1,-1,+1)*s0xp)*shym       &
                 +(+tmpf(4,i0-i-1, 0,+1)*s0xm+tmpf(4,i0-i, 0,+1)*s0x+tmpf(4,i0-i+1, 0,+1)*s0xp)*shy        &
                 +(+tmpf(4,i0-i-1,+1,+1)*s0xm+tmpf(4,i0-i,+1,+1)*s0x+tmpf(4,i0-i+1,+1,+1)*s0xp)*shyp)*shzp

          epy = (+(+tmpf(5,-1,j0-j-1,-1)*shxm+tmpf(5,0,j0-j-1,-1)*shx+tmpf(5,+1,j0-j-1,-1)*shxp)*s0ym       &
                 +(+tmpf(5,-1,j0-j  ,-1)*shxm+tmpf(5,0,j0-j  ,-1)*shx+tmpf(5,+1,j0-j  ,-1)*shxp)*s0y        &
                 +(+tmpf(5,-1,j0-j+1,-1)*shxm+tmpf(5,0,j0-j+1,-1)*shx+tmpf(5,+1,j0-j+1,-1)*shxp)*s0yp)*shzm &
               +(+(+tmpf(5,-1,j0-j-1, 0)*shxm+tmpf(5,0,j0-j-1, 0)*shx+tmpf(5,+1,j0-j-1, 0)*shxp)*s0ym       &
                 +(+tmpf(5,-1,j0-j  , 0)*shxm+tmpf(5,0,j0-j  , 0)*shx+tmpf(5,+1,j0-j  , 0)*shxp)*s0y        &
                 +(+tmpf(5,-1,j0-j+1, 0)*shxm+tmpf(5,0,j0-j+1, 0)*shx+tmpf(5,+1,j0-j+1, 0)*shxp)*s0yp)*shz  &
               +(+(+tmpf(5,-1,j0-j-1,+1)*shxm+tmpf(5,0,j0-j-1,+1)*shx+tmpf(5,+1,j0-j-1,+1)*shxp)*s0ym       &
                 +(+tmpf(5,-1,j0-j  ,+1)*shxm+tmpf(5,0,j0-j  ,+1)*shx+tmpf(5,+1,j0-j  ,+1)*shxp)*s0y        &
                 +(+tmpf(5,-1,j0-j+1,+1)*shxm+tmpf(5,0,j0-j+1,+1)*shx+tmpf(5,+1,j0-j+1,+1)*shxp)*s0yp)*shzp

          epz = (+(+tmpf(6,-1,-1,k0-k-1)*shxm+tmpf(6,0,-1,k0-k-1)*shx+tmpf(6,+1,-1,k0-k-1)*shxp)*shym       &
                 +(+tmpf(6,-1, 0,k0-k-1)*shxm+tmpf(6,0, 0,k0-k-1)*shx+tmpf(6,+1, 0,k0-k-1)*shxp)*shy        &
                 +(+tmpf(6,-1,+1,k0-k-1)*shxm+tmpf(6,0,+1,k0-k-1)*shx+tmpf(6,+1,+1,k0-k-1)*shxp)*shyp)*s0zm &
               +(+(+tmpf(6,-1,-1,k0-k  )*shxm+tmpf(6,0,-1,k0-k  )*shx+tmpf(6,+1,-1,k0-k  )*shxp)*shym       &
                 +(+tmpf(6,-1, 0,k0-k  )*shxm+tmpf(6,0, 0,k0-k  )*shx+tmpf(6,+1, 0,k0-k  )*shxp)*shy        &
                 +(+tmpf(6,-1,+1,k0-k  )*shxm+tmpf(6,0,+1,k0-k  )*shx+tmpf(6,+1,+1,k0-k  )*shxp)*shyp)*s0z  &
               +(+(+tmpf(6,-1,-1,k0-k+1)*shxm+tmpf(6,0,-1,k0-k+1)*shx+tmpf(6,+1,-1,k0-k+1)*shxp)*shym       &
                 +(+tmpf(6,-1, 0,k0-k+1)*shxm+tmpf(6,0, 0,k0-k+1)*shx+tmpf(6,+1, 0,k0-k+1)*shxp)*shy        &
                 +(+tmpf(6,-1,+1,k0-k+1)*shxm+tmpf(6,0,+1,k0-k+1)*shx+tmpf(6,+1,+1,k0-k+1)*shxp)*shyp)*s0zp

          !accel.
          uvm1 = up(4,ii,j,k,isp)+fac1*epx
          uvm2 = up(5,ii,j,k,isp)+fac1*epy
          uvm3 = up(6,ii,j,k,isp)+fac1*epz

          !rotate
          gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
          igam = 1./gam
          fac1r = fac1*igam
          fac2r = fac2/(gam+txxx*(bpx*bpx+bpy*bpy+bpz*bpz)*igam)

          uvm4 = uvm1+fac1r*(+uvm2*bpz-uvm3*bpy)
          uvm5 = uvm2+fac1r*(+uvm3*bpx-uvm1*bpz)
          uvm6 = uvm3+fac1r*(+uvm1*bpy-uvm2*bpx)

          uvm1 = uvm1+fac2r*(+uvm5*bpz-uvm6*bpy)
          uvm2 = uvm2+fac2r*(+uvm6*bpx-uvm4*bpz)
          uvm3 = uvm3+fac2r*(+uvm4*bpy-uvm5*bpx)

          !accel.
          gp(1,ii,j,k,isp) = up(1,ii,j,k,isp)
          gp(2,ii,j,k,isp) = up(2,ii,j,k,isp)
          gp(3,ii,j,k,isp) = up(3,ii,j,k,isp)
          gp(4,ii,j,k,isp) = uvm1+fac1*epx
          gp(5,ii,j,k,isp) = uvm2+fac1*epy
          gp(6,ii,j,k,isp) = uvm3+fac1*epz

       enddo

       isp = 2

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

       do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)

          !second order shape function
          i0 = int(up(1,ii,j,k,isp)*idelx+0.5)
          dh = up(1,ii,j,k,isp)*idelx-i0
          s0xm = 0.5*(0.5-dh)*(0.5-dh)
          s0x  = 0.75-dh*dh
          s0xp = 0.5*(0.5+dh)*(0.5+dh)

          j0 = int(up(2,ii,j,k,isp)*idelx+0.5)
          dh = up(2,ii,j,k,isp)*idelx-j0
          s0ym = 0.5*(0.5-dh)*(0.5-dh)
          s0y  = 0.75-dh*dh
          s0yp = 0.5*(0.5+dh)*(0.5+dh)

          k0 = int(up(3,ii,j,k,isp)*idelx+0.5)
          dh = up(3,ii,j,k,isp)*idelx-k0
          s0zm = 0.5*(0.5-dh)*(0.5-dh)
          s0z  = 0.75-dh*dh
          s0zp = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(1,ii,j,k,isp)*idelx-0.5-i
          shxm = 0.5*(0.5-dh)*(0.5-dh)
          shx  = 0.75-dh*dh
          shxp = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(2,ii,j,k,isp)*idelx-0.5-j
          shym = 0.5*(0.5-dh)*(0.5-dh)
          shy  = 0.75-dh*dh
          shyp = 0.5*(0.5+dh)*(0.5+dh)

          dh = up(3,ii,j,k,isp)*idelx-0.5-k
          shzm = 0.5*(0.5-dh)*(0.5-dh)
          shz  = 0.75-dh*dh
          shzp = 0.5*(0.5+dh)*(0.5+dh)

          bpx = (+(+tmpf(1,-1,j0-j-1,k0-k-1)*shxm+tmpf(1,0,j0-j-1,k0-k-1)*shx+tmpf(1,+1,j0-j-1,k0-k-1)*shxp)*s0ym       &
                 +(+tmpf(1,-1,j0-j  ,k0-k-1)*shxm+tmpf(1,0,j0-j  ,k0-k-1)*shx+tmpf(1,+1,j0-j  ,k0-k-1)*shxp)*s0y        &
                 +(+tmpf(1,-1,j0-j+1,k0-k-1)*shxm+tmpf(1,0,j0-j+1,k0-k-1)*shx+tmpf(1,+1,j0-j+1,k0-k-1)*shxp)*s0yp)*s0zm &
               +(+(+tmpf(1,-1,j0-j-1,k0-k  )*shxm+tmpf(1,0,j0-j-1,k0-k  )*shx+tmpf(1,+1,j0-j-1,k0-k  )*shxp)*s0ym       &
                 +(+tmpf(1,-1,j0-j  ,k0-k  )*shxm+tmpf(1,0,j0-j  ,k0-k  )*shx+tmpf(1,+1,j0-j  ,k0-k  )*shxp)*s0y        &
                 +(+tmpf(1,-1,j0-j+1,k0-k  )*shxm+tmpf(1,0,j0-j+1,k0-k  )*shx+tmpf(1,+1,j0-j+1,k0-k  )*shxp)*s0yp)*s0z  &
               +(+(+tmpf(1,-1,j0-j-1,k0-k+1)*shxm+tmpf(1,0,j0-j-1,k0-k+1)*shx+tmpf(1,+1,j0-j-1,k0-k+1)*shxp)*s0ym       &
                 +(+tmpf(1,-1,j0-j  ,k0-k+1)*shxm+tmpf(1,0,j0-j  ,k0-k+1)*shx+tmpf(1,+1,j0-j  ,k0-k+1)*shxp)*s0y        &
                 +(+tmpf(1,-1,j0-j+1,k0-k+1)*shxm+tmpf(1,0,j0-j+1,k0-k+1)*shx+tmpf(1,+1,j0-j+1,k0-k+1)*shxp)*s0yp)*s0zp

          bpy = (+(+tmpf(2,i0-i-1,-1,k0-k-1)*s0xm+tmpf(2,i0-i,-1,k0-k-1)*s0x+tmpf(2,i0-i+1,-1,k0-k-1)*s0xp)*shym       &
                 +(+tmpf(2,i0-i-1, 0,k0-k-1)*s0xm+tmpf(2,i0-i, 0,k0-k-1)*s0x+tmpf(2,i0-i+1, 0,k0-k-1)*s0xp)*shy        &
                 +(+tmpf(2,i0-i-1,+1,k0-k-1)*s0xm+tmpf(2,i0-i,+1,k0-k-1)*s0x+tmpf(2,i0-i+1,+1,k0-k-1)*s0xp)*shyp)*s0zm &
               +(+(+tmpf(2,i0-i-1,-1,k0-k  )*s0xm+tmpf(2,i0-i,-1,k0-k  )*s0x+tmpf(2,i0-i+1,-1,k0-k  )*s0xp)*shym       &
                 +(+tmpf(2,i0-i-1, 0,k0-k  )*s0xm+tmpf(2,i0-i, 0,k0-k  )*s0x+tmpf(2,i0-i+1, 0,k0-k  )*s0xp)*shy        &
                 +(+tmpf(2,i0-i-1,+1,k0-k  )*s0xm+tmpf(2,i0-i,+1,k0-k  )*s0x+tmpf(2,i0-i+1,+1,k0-k  )*s0xp)*shyp)*s0z  &
               +(+(+tmpf(2,i0-i-1,-1,k0-k+1)*s0xm+tmpf(2,i0-i,-1,k0-k+1)*s0x+tmpf(2,i0-i+1,-1,k0-k+1)*s0xp)*shym       &
                 +(+tmpf(2,i0-i-1, 0,k0-k+1)*s0xm+tmpf(2,i0-i, 0,k0-k+1)*s0x+tmpf(2,i0-i+1, 0,k0-k+1)*s0xp)*shy        &
                 +(+tmpf(2,i0-i-1,+1,k0-k+1)*s0xm+tmpf(2,i0-i,+1,k0-k+1)*s0x+tmpf(2,i0-i+1,+1,k0-k+1)*s0xp)*shyp)*s0zp

          bpz = (+(+tmpf(3,i0-i-1,j0-j-1,-1)*s0xm+tmpf(3,i0-i,j0-j-1,-1)*s0x+tmpf(3,i0-i+1,j0-j-1,-1)*s0xp)*s0ym       &
                 +(+tmpf(3,i0-i-1,j0-j  ,-1)*s0xm+tmpf(3,i0-i,j0-j  ,-1)*s0x+tmpf(3,i0-i+1,j0-j  ,-1)*s0xp)*s0y        &
                 +(+tmpf(3,i0-i-1,j0-j+1,-1)*s0xm+tmpf(3,i0-i,j0-j+1,-1)*s0x+tmpf(3,i0-i+1,j0-j+1,-1)*s0xp)*s0yp)*shzm &
               +(+(+tmpf(3,i0-i-1,j0-j-1, 0)*s0xm+tmpf(3,i0-i,j0-j-1, 0)*s0x+tmpf(3,i0-i+1,j0-j-1, 0)*s0xp)*s0ym       &
                 +(+tmpf(3,i0-i-1,j0-j  , 0)*s0xm+tmpf(3,i0-i,j0-j  , 0)*s0x+tmpf(3,i0-i+1,j0-j  , 0)*s0xp)*s0y        &
                 +(+tmpf(3,i0-i-1,j0-j+1, 0)*s0xm+tmpf(3,i0-i,j0-j+1, 0)*s0x+tmpf(3,i0-i+1,j0-j+1, 0)*s0xp)*s0yp)*shz  &
               +(+(+tmpf(3,i0-i-1,j0-j-1,+1)*s0xm+tmpf(3,i0-i,j0-j-1,+1)*s0x+tmpf(3,i0-i+1,j0-j-1,+1)*s0xp)*s0ym       &
                 +(+tmpf(3,i0-i-1,j0-j  ,+1)*s0xm+tmpf(3,i0-i,j0-j  ,+1)*s0x+tmpf(3,i0-i+1,j0-j  ,+1)*s0xp)*s0y        &
                 +(+tmpf(3,i0-i-1,j0-j+1,+1)*s0xm+tmpf(3,i0-i,j0-j+1,+1)*s0x+tmpf(3,i0-i+1,j0-j+1,+1)*s0xp)*s0yp)*shzp

          epx = (+(+tmpf(4,i0-i-1,-1,-1)*s0xm+tmpf(4,i0-i,-1,-1)*s0x+tmpf(4,i0-i+1,-1,-1)*s0xp)*shym       &
                 +(+tmpf(4,i0-i-1, 0,-1)*s0xm+tmpf(4,i0-i, 0,-1)*s0x+tmpf(4,i0-i+1, 0,-1)*s0xp)*shy        &
                 +(+tmpf(4,i0-i-1,+1,-1)*s0xm+tmpf(4,i0-i,+1,-1)*s0x+tmpf(4,i0-i+1,+1,-1)*s0xp)*shyp)*shzm &
               +(+(+tmpf(4,i0-i-1,-1, 0)*s0xm+tmpf(4,i0-i,-1, 0)*s0x+tmpf(4,i0-i+1,-1, 0)*s0xp)*shym       &
                 +(+tmpf(4,i0-i-1, 0, 0)*s0xm+tmpf(4,i0-i, 0, 0)*s0x+tmpf(4,i0-i+1, 0, 0)*s0xp)*shy        &
                 +(+tmpf(4,i0-i-1,+1, 0)*s0xm+tmpf(4,i0-i,+1, 0)*s0x+tmpf(4,i0-i+1,+1, 0)*s0xp)*shyp)*shz  &
               +(+(+tmpf(4,i0-i-1,-1,+1)*s0xm+tmpf(4,i0-i,-1,+1)*s0x+tmpf(4,i0-i+1,-1,+1)*s0xp)*shym       &
                 +(+tmpf(4,i0-i-1, 0,+1)*s0xm+tmpf(4,i0-i, 0,+1)*s0x+tmpf(4,i0-i+1, 0,+1)*s0xp)*shy        &
                 +(+tmpf(4,i0-i-1,+1,+1)*s0xm+tmpf(4,i0-i,+1,+1)*s0x+tmpf(4,i0-i+1,+1,+1)*s0xp)*shyp)*shzp

          epy = (+(+tmpf(5,-1,j0-j-1,-1)*shxm+tmpf(5,0,j0-j-1,-1)*shx+tmpf(5,+1,j0-j-1,-1)*shxp)*s0ym       &
                 +(+tmpf(5,-1,j0-j  ,-1)*shxm+tmpf(5,0,j0-j  ,-1)*shx+tmpf(5,+1,j0-j  ,-1)*shxp)*s0y        &
                 +(+tmpf(5,-1,j0-j+1,-1)*shxm+tmpf(5,0,j0-j+1,-1)*shx+tmpf(5,+1,j0-j+1,-1)*shxp)*s0yp)*shzm &
               +(+(+tmpf(5,-1,j0-j-1, 0)*shxm+tmpf(5,0,j0-j-1, 0)*shx+tmpf(5,+1,j0-j-1, 0)*shxp)*s0ym       &
                 +(+tmpf(5,-1,j0-j  , 0)*shxm+tmpf(5,0,j0-j  , 0)*shx+tmpf(5,+1,j0-j  , 0)*shxp)*s0y        &
                 +(+tmpf(5,-1,j0-j+1, 0)*shxm+tmpf(5,0,j0-j+1, 0)*shx+tmpf(5,+1,j0-j+1, 0)*shxp)*s0yp)*shz  &
               +(+(+tmpf(5,-1,j0-j-1,+1)*shxm+tmpf(5,0,j0-j-1,+1)*shx+tmpf(5,+1,j0-j-1,+1)*shxp)*s0ym       &
                 +(+tmpf(5,-1,j0-j  ,+1)*shxm+tmpf(5,0,j0-j  ,+1)*shx+tmpf(5,+1,j0-j  ,+1)*shxp)*s0y        &
                 +(+tmpf(5,-1,j0-j+1,+1)*shxm+tmpf(5,0,j0-j+1,+1)*shx+tmpf(5,+1,j0-j+1,+1)*shxp)*s0yp)*shzp

          epz = (+(+tmpf(6,-1,-1,k0-k-1)*shxm+tmpf(6,0,-1,k0-k-1)*shx+tmpf(6,+1,-1,k0-k-1)*shxp)*shym       &
                 +(+tmpf(6,-1, 0,k0-k-1)*shxm+tmpf(6,0, 0,k0-k-1)*shx+tmpf(6,+1, 0,k0-k-1)*shxp)*shy        &
                 +(+tmpf(6,-1,+1,k0-k-1)*shxm+tmpf(6,0,+1,k0-k-1)*shx+tmpf(6,+1,+1,k0-k-1)*shxp)*shyp)*s0zm &
               +(+(+tmpf(6,-1,-1,k0-k  )*shxm+tmpf(6,0,-1,k0-k  )*shx+tmpf(6,+1,-1,k0-k  )*shxp)*shym       &
                 +(+tmpf(6,-1, 0,k0-k  )*shxm+tmpf(6,0, 0,k0-k  )*shx+tmpf(6,+1, 0,k0-k  )*shxp)*shy        &
                 +(+tmpf(6,-1,+1,k0-k  )*shxm+tmpf(6,0,+1,k0-k  )*shx+tmpf(6,+1,+1,k0-k  )*shxp)*shyp)*s0z  &
               +(+(+tmpf(6,-1,-1,k0-k+1)*shxm+tmpf(6,0,-1,k0-k+1)*shx+tmpf(6,+1,-1,k0-k+1)*shxp)*shym       &
                 +(+tmpf(6,-1, 0,k0-k+1)*shxm+tmpf(6,0, 0,k0-k+1)*shx+tmpf(6,+1, 0,k0-k+1)*shxp)*shy        &
                 +(+tmpf(6,-1,+1,k0-k+1)*shxm+tmpf(6,0,+1,k0-k+1)*shx+tmpf(6,+1,+1,k0-k+1)*shxp)*shyp)*s0zp

          !accel.
          uvm1 = up(4,ii,j,k,isp)+fac1*epx
          uvm2 = up(5,ii,j,k,isp)+fac1*epy
          uvm3 = up(6,ii,j,k,isp)+fac1*epz

          !rotate
          gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
          igam = 1./gam
          fac1r = fac1*igam
          fac2r = fac2/(gam+txxx*(bpx*bpx+bpy*bpy+bpz*bpz)*igam)

          uvm4 = uvm1+fac1r*(+uvm2*bpz-uvm3*bpy)
          uvm5 = uvm2+fac1r*(+uvm3*bpx-uvm1*bpz)
          uvm6 = uvm3+fac1r*(+uvm1*bpy-uvm2*bpx)

          uvm1 = uvm1+fac2r*(+uvm5*bpz-uvm6*bpy)
          uvm2 = uvm2+fac2r*(+uvm6*bpx-uvm4*bpz)
          uvm3 = uvm3+fac2r*(+uvm4*bpy-uvm5*bpx)

          !accel.
          gp(1,ii,j,k,isp) = up(1,ii,j,k,isp)
          gp(2,ii,j,k,isp) = up(2,ii,j,k,isp)
          gp(3,ii,j,k,isp) = up(3,ii,j,k,isp)
          gp(4,ii,j,k,isp) = uvm1+fac1*epx
          gp(5,ii,j,k,isp) = uvm2+fac1*epy
          gp(6,ii,j,k,isp) = uvm3+fac1*epz
       enddo

    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine mom_calc__accl


  subroutine mom_calc__nvt(den,vel,temp,up,nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2,c)

    integer, intent(in)    :: nxgs, nxge, nys, nye, nzs, nze
    integer, intent(in)    :: np, nsp, np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in)    :: c
    real(8), intent(in)    :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    real(8), intent(inout) :: vel(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp)
    real(8), intent(inout) :: temp(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp)

    integer :: ii, ih, j, jh, k, kh, isp
    real(8) :: dx, dxm, dy, dym, dz, dzm, gam

!$OMP PARALLEL WORKSHARE
    den(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = 0.0D0
    vel(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = 0.0D0
    temp(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = 0.0D0
!$OMP END PARALLEL WORKSHARE

    !caluculate number density at (i+1/2, j+1/2, k+1/2)
    do isp=1,nsp
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)

             gam = 1./dsqrt(1.0+(+up(4,ii,j,k,isp)*up(4,ii,j,k,isp) &
                                 +up(5,ii,j,k,isp)*up(5,ii,j,k,isp) &
                                 +up(6,ii,j,k,isp)*up(6,ii,j,k,isp) &
                                )/(c*c))

             ih = floor(up(1,ii,j,k,isp)-0.5)
             jh = floor(up(2,ii,j,k,isp)-0.5)
             kh = floor(up(3,ii,j,k,isp)-0.5)

             dx = up(1,ii,j,k,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,k,isp)-0.5-jh
             dym = 1.-dy
             dz = up(3,ii,j,k,isp)-0.5-kh
             dzm = 1.-dz

             !N
             den(ih  ,jh  ,kh,isp) = den(ih  ,jh  ,kh,isp)+dxm*dym*dzm
             den(ih+1,jh  ,kh,isp) = den(ih+1,jh  ,kh,isp)+dx *dym*dzm
             den(ih  ,jh+1,kh,isp) = den(ih  ,jh+1,kh,isp)+dxm*dy *dzm
             den(ih+1,jh+1,kh,isp) = den(ih+1,jh+1,kh,isp)+dx *dy *dzm
             den(ih  ,jh  ,kh+1,isp) = den(ih  ,jh  ,kh+1,isp)+dxm*dym*dz
             den(ih+1,jh  ,kh+1,isp) = den(ih+1,jh  ,kh+1,isp)+dx *dym*dz
             den(ih  ,jh+1,kh+1,isp) = den(ih  ,jh+1,kh+1,isp)+dxm*dy *dz
             den(ih+1,jh+1,kh+1,isp) = den(ih+1,jh+1,kh+1,isp)+dx *dy *dz

             !Vx
             vel(ih  ,jh  ,kh  ,1,isp) = vel(ih  ,jh  ,kh  ,1,isp)+up(4,ii,j,k,isp)*gam*dxm*dym*dzm
             vel(ih+1,jh  ,kh  ,1,isp) = vel(ih+1,jh  ,kh  ,1,isp)+up(4,ii,j,k,isp)*gam*dx *dym*dzm
             vel(ih  ,jh+1,kh  ,1,isp) = vel(ih  ,jh+1,kh  ,1,isp)+up(4,ii,j,k,isp)*gam*dxm*dy *dzm
             vel(ih+1,jh+1,kh  ,1,isp) = vel(ih+1,jh+1,kh  ,1,isp)+up(4,ii,j,k,isp)*gam*dx *dy *dzm
             vel(ih  ,jh  ,kh+1,1,isp) = vel(ih  ,jh  ,kh+1,1,isp)+up(4,ii,j,k,isp)*gam*dxm*dym*dz
             vel(ih+1,jh  ,kh+1,1,isp) = vel(ih+1,jh  ,kh+1,1,isp)+up(4,ii,j,k,isp)*gam*dx *dym*dz
             vel(ih  ,jh+1,kh+1,1,isp) = vel(ih  ,jh+1,kh+1,1,isp)+up(4,ii,j,k,isp)*gam*dxm*dy *dz
             vel(ih+1,jh+1,kh+1,1,isp) = vel(ih+1,jh+1,kh+1,1,isp)+up(4,ii,j,k,isp)*gam*dx *dy *dz

             !Vy
             vel(ih  ,jh  ,kh  ,2,isp) = vel(ih  ,jh  ,kh  ,2,isp)+up(5,ii,j,k,isp)*gam*dxm*dym*dzm
             vel(ih+1,jh  ,kh  ,2,isp) = vel(ih+1,jh  ,kh  ,2,isp)+up(5,ii,j,k,isp)*gam*dx *dym*dzm
             vel(ih  ,jh+1,kh  ,2,isp) = vel(ih  ,jh+1,kh  ,2,isp)+up(5,ii,j,k,isp)*gam*dxm*dy *dzm
             vel(ih+1,jh+1,kh  ,2,isp) = vel(ih+1,jh+1,kh  ,2,isp)+up(5,ii,j,k,isp)*gam*dx *dy *dzm
             vel(ih  ,jh  ,kh+1,2,isp) = vel(ih  ,jh  ,kh+1,2,isp)+up(5,ii,j,k,isp)*gam*dxm*dym*dz
             vel(ih+1,jh  ,kh+1,2,isp) = vel(ih+1,jh  ,kh+1,2,isp)+up(5,ii,j,k,isp)*gam*dx *dym*dz
             vel(ih  ,jh+1,kh+1,2,isp) = vel(ih  ,jh+1,kh+1,2,isp)+up(5,ii,j,k,isp)*gam*dxm*dy *dz
             vel(ih+1,jh+1,kh+1,2,isp) = vel(ih+1,jh+1,kh+1,2,isp)+up(5,ii,j,k,isp)*gam*dx *dy *dz

             !Vz
             vel(ih  ,jh  ,kh  ,3,isp) = vel(ih  ,jh  ,kh  ,3,isp)+up(6,ii,j,k,isp)*gam*dxm*dym*dzm
             vel(ih+1,jh  ,kh  ,3,isp) = vel(ih+1,jh  ,kh  ,3,isp)+up(6,ii,j,k,isp)*gam*dx *dym*dzm
             vel(ih  ,jh+1,kh  ,3,isp) = vel(ih  ,jh+1,kh  ,3,isp)+up(6,ii,j,k,isp)*gam*dxm*dy *dzm
             vel(ih+1,jh+1,kh  ,3,isp) = vel(ih+1,jh+1,kh  ,3,isp)+up(6,ii,j,k,isp)*gam*dx *dy *dzm
             vel(ih  ,jh  ,kh+1,3,isp) = vel(ih  ,jh  ,kh+1,3,isp)+up(6,ii,j,k,isp)*gam*dxm*dym*dz
             vel(ih+1,jh  ,kh+1,3,isp) = vel(ih+1,jh  ,kh+1,3,isp)+up(6,ii,j,k,isp)*gam*dx *dym*dz
             vel(ih  ,jh+1,kh+1,3,isp) = vel(ih  ,jh+1,kh+1,3,isp)+up(6,ii,j,k,isp)*gam*dxm*dy *dz
             vel(ih+1,jh+1,kh+1,3,isp) = vel(ih+1,jh+1,kh+1,3,isp)+up(6,ii,j,k,isp)*gam*dx *dy *dz

             !Txx
             temp(ih  ,jh  ,kh  ,1,isp) = temp(ih  ,jh  ,kh  ,1,isp)+up(4,ii,j,k,isp)**2*gam*dxm*dym*dzm
             temp(ih+1,jh  ,kh  ,1,isp) = temp(ih+1,jh  ,kh  ,1,isp)+up(4,ii,j,k,isp)**2*gam*dx *dym*dzm
             temp(ih  ,jh+1,kh  ,1,isp) = temp(ih  ,jh+1,kh  ,1,isp)+up(4,ii,j,k,isp)**2*gam*dxm*dy *dzm
             temp(ih+1,jh+1,kh  ,1,isp) = temp(ih+1,jh+1,kh  ,1,isp)+up(4,ii,j,k,isp)**2*gam*dx *dy *dzm
             temp(ih  ,jh  ,kh+1,1,isp) = temp(ih  ,jh  ,kh+1,1,isp)+up(4,ii,j,k,isp)**2*gam*dxm*dym*dz
             temp(ih+1,jh  ,kh+1,1,isp) = temp(ih+1,jh  ,kh+1,1,isp)+up(4,ii,j,k,isp)**2*gam*dx *dym*dz
             temp(ih  ,jh+1,kh+1,1,isp) = temp(ih  ,jh+1,kh+1,1,isp)+up(4,ii,j,k,isp)**2*gam*dxm*dy *dz
             temp(ih+1,jh+1,kh+1,1,isp) = temp(ih+1,jh+1,kh+1,1,isp)+up(4,ii,j,k,isp)**2*gam*dx *dy *dz

             !Tyy
             temp(ih  ,jh  ,kh  ,2,isp) = temp(ih  ,jh  ,kh  ,2,isp)+up(5,ii,j,k,isp)**2*gam*dxm*dym*dzm
             temp(ih+1,jh  ,kh  ,2,isp) = temp(ih+1,jh  ,kh  ,2,isp)+up(5,ii,j,k,isp)**2*gam*dx *dym*dzm
             temp(ih  ,jh+1,kh  ,2,isp) = temp(ih  ,jh+1,kh  ,2,isp)+up(5,ii,j,k,isp)**2*gam*dxm*dy *dzm
             temp(ih+1,jh+1,kh  ,2,isp) = temp(ih+1,jh+1,kh  ,2,isp)+up(5,ii,j,k,isp)**2*gam*dx *dy *dzm
             temp(ih  ,jh  ,kh+1,2,isp) = temp(ih  ,jh  ,kh+1,2,isp)+up(5,ii,j,k,isp)**2*gam*dxm*dym*dz
             temp(ih+1,jh  ,kh+1,2,isp) = temp(ih+1,jh  ,kh+1,2,isp)+up(5,ii,j,k,isp)**2*gam*dx *dym*dz
             temp(ih  ,jh+1,kh+1,2,isp) = temp(ih  ,jh+1,kh+1,2,isp)+up(5,ii,j,k,isp)**2*gam*dxm*dy *dz
             temp(ih+1,jh+1,kh+1,2,isp) = temp(ih+1,jh+1,kh+1,2,isp)+up(5,ii,j,k,isp)**2*gam*dx *dy *dz

             !Tzz
             temp(ih  ,jh  ,kh  ,3,isp) = temp(ih  ,jh  ,kh  ,3,isp)+up(6,ii,j,k,isp)**2*gam*dxm*dym*dzm
             temp(ih+1,jh  ,kh  ,3,isp) = temp(ih+1,jh  ,kh  ,3,isp)+up(6,ii,j,k,isp)**2*gam*dx *dym*dzm
             temp(ih  ,jh+1,kh  ,3,isp) = temp(ih  ,jh+1,kh  ,3,isp)+up(6,ii,j,k,isp)**2*gam*dxm*dy *dzm
             temp(ih+1,jh+1,kh  ,3,isp) = temp(ih+1,jh+1,kh  ,3,isp)+up(6,ii,j,k,isp)**2*gam*dx *dy *dzm
             temp(ih  ,jh  ,kh+1,3,isp) = temp(ih  ,jh  ,kh+1,3,isp)+up(6,ii,j,k,isp)**2*gam*dxm*dym*dz
             temp(ih+1,jh  ,kh+1,3,isp) = temp(ih+1,jh  ,kh+1,3,isp)+up(6,ii,j,k,isp)**2*gam*dx *dym*dz
             temp(ih  ,jh+1,kh+1,3,isp) = temp(ih  ,jh+1,kh+1,3,isp)+up(6,ii,j,k,isp)**2*gam*dxm*dy *dz
             temp(ih+1,jh+1,kh+1,3,isp) = temp(ih+1,jh+1,kh+1,3,isp)+up(6,ii,j,k,isp)**2*gam*dx *dy *dz
          enddo
       enddo
       enddo
    enddo

  end subroutine mom_calc__nvt


end module mom_calc
