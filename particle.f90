module particle

  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(gp,                                              &
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
          gp(4,ii,j,k,isp) = uvm1+fac1*epx
          gp(5,ii,j,k,isp) = uvm2+fac1*epy
          gp(6,ii,j,k,isp) = uvm3+fac1*epz

          !move
          gam = 1./dsqrt(1.0+(+gp(4,ii,j,k,isp)*gp(4,ii,j,k,isp) &
                               +gp(5,ii,j,k,isp)*gp(5,ii,j,k,isp) &
                               +gp(6,ii,j,k,isp)*gp(6,ii,j,k,isp))/(c*c))

          gp(1,ii,j,k,isp) = up(1,ii,j,k,isp)+gp(4,ii,j,k,isp)*delt*gam
          gp(2,ii,j,k,isp) = up(2,ii,j,k,isp)+gp(5,ii,j,k,isp)*delt*gam
          gp(3,ii,j,k,isp) = up(3,ii,j,k,isp)+gp(6,ii,j,k,isp)*delt*gam

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
          gp(4,ii,j,k,isp) = uvm1+fac1*epx
          gp(5,ii,j,k,isp) = uvm2+fac1*epy
          gp(6,ii,j,k,isp) = uvm3+fac1*epz

          !move
          gam = 1./dsqrt(1.0+(+gp(4,ii,j,k,isp)*gp(4,ii,j,k,isp) &
                               +gp(5,ii,j,k,isp)*gp(5,ii,j,k,isp) &
                               +gp(6,ii,j,k,isp)*gp(6,ii,j,k,isp))/(c*c))

          gp(1,ii,j,k,isp) = up(1,ii,j,k,isp)+gp(4,ii,j,k,isp)*delt*gam
          gp(2,ii,j,k,isp) = up(2,ii,j,k,isp)+gp(5,ii,j,k,isp)*delt*gam
          gp(3,ii,j,k,isp) = up(3,ii,j,k,isp)+gp(6,ii,j,k,isp)*delt*gam

       enddo

    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine particle__solv


end module particle
