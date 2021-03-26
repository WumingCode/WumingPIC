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

    integer            :: i, j, k, ii, isp
    real(8)            :: tmpf(6,nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: idelx, dh, fac1, fac1r, fac2, fac2r, gam, igam, txxx
    real(8)            :: bpx, bpy, bpz, epx, epy, epz
    real(8)            :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8)            :: shxm, shx, shxp, shym, shy, shyp, shzm, shz, shzp

    idelx = 1./delx

    !fields at (i+1/2, j+1/2, k+1/2)
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs-1,nze+1
      do j=nys-1,nye+1
        do i=nxs-1,nxe+1
           tmpf(1,i,j,k) = 0.25*(+uf(1,i,j,k  )+uf(1,i,j+1,k  ) &
                                 +uf(1,i,j,k+1)+uf(1,i,j+1,k+1))
           tmpf(2,i,j,k) = 0.25*(+uf(2,i,j,k  )+uf(2,i+1,j,k  ) &
                                 +uf(2,i,j,k+1)+uf(2,i+1,j,k+1))
           tmpf(3,i,j,k) = 0.25*(+uf(3,i,j  ,k)+uf(3,i+1,j  ,k) &
                                 +uf(3,i,j+1,k)+uf(3,i+1,j+1,k))
           tmpf(4,i,j,k) = 0.5*(+uf(4,i,j,k)+uf(4,i+1,j  ,k  ))
           tmpf(5,i,j,k) = 0.5*(+uf(5,i,j,k)+uf(5,i  ,j+1,k  ))
           tmpf(6,i,j,k) = 0.5*(+uf(6,i,j,k)+uf(6,i  ,j  ,k+1))
        enddo
      enddo
    enddo
!OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,i,j,k,isp,                                &
!$OMP                     dh,gam,igam,fac1,fac2,txxx,fac1r,fac2r,      &
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

       do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)

          !second order shape function
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

          bpx = (+(+tmpf(1,i-1,j-1,k-1)*shxm+tmpf(1,i,j-1,k-1)*shx+tmpf(1,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(1,i-1,j  ,k-1)*shxm+tmpf(1,i,j  ,k-1)*shx+tmpf(1,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(1,i-1,j+1,k-1)*shxm+tmpf(1,i,j+1,k-1)*shx+tmpf(1,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(1,i-1,j-1,k  )*shxm+tmpf(1,i,j-1,k  )*shx+tmpf(1,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(1,i-1,j  ,k  )*shxm+tmpf(1,i,j  ,k  )*shx+tmpf(1,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(1,i-1,j+1,k  )*shxm+tmpf(1,i,j+1,k  )*shx+tmpf(1,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(1,i-1,j-1,k+1)*shxm+tmpf(1,i,j-1,k+1)*shx+tmpf(1,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(1,i-1,j  ,k+1)*shxm+tmpf(1,i,j  ,k+1)*shx+tmpf(1,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(1,i-1,j+1,k+1)*shxm+tmpf(1,i,j+1,k+1)*shx+tmpf(1,i+1,j+1,k+1)*shxp)*shyp)*shzp

          bpy = (+(+tmpf(2,i-1,j-1,k-1)*shxm+tmpf(2,i,j-1,k-1)*shx+tmpf(2,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(2,i-1,j  ,k-1)*shxm+tmpf(2,i,j  ,k-1)*shx+tmpf(2,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(2,i-1,j+1,k-1)*shxm+tmpf(2,i,j+1,k-1)*shx+tmpf(2,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(2,i-1,j-1,k  )*shxm+tmpf(2,i,j-1,k  )*shx+tmpf(2,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(2,i-1,j  ,k  )*shxm+tmpf(2,i,j  ,k  )*shx+tmpf(2,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(2,i-1,j+1,k  )*shxm+tmpf(2,i,j+1,k  )*shx+tmpf(2,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(2,i-1,j-1,k+1)*shxm+tmpf(2,i,j-1,k+1)*shx+tmpf(2,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(2,i-1,j  ,k+1)*shxm+tmpf(2,i,j  ,k+1)*shx+tmpf(2,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(2,i-1,j+1,k+1)*shxm+tmpf(2,i,j+1,k+1)*shx+tmpf(2,i+1,j+1,k+1)*shxp)*shyp)*shzp

          bpz = (+(+tmpf(3,i-1,j-1,k-1)*shxm+tmpf(3,i,j-1,k-1)*shx+tmpf(3,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(3,i-1,j  ,k-1)*shxm+tmpf(3,i,j  ,k-1)*shx+tmpf(3,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(3,i-1,j+1,k-1)*shxm+tmpf(3,i,j+1,k-1)*shx+tmpf(3,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(3,i-1,j-1,k  )*shxm+tmpf(3,i,j-1,k  )*shx+tmpf(3,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(3,i-1,j  ,k  )*shxm+tmpf(3,i,j  ,k  )*shx+tmpf(3,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(3,i-1,j+1,k  )*shxm+tmpf(3,i,j+1,k  )*shx+tmpf(3,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(3,i-1,j-1,k+1)*shxm+tmpf(3,i,j-1,k+1)*shx+tmpf(3,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(3,i-1,j  ,k+1)*shxm+tmpf(3,i,j  ,k+1)*shx+tmpf(3,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(3,i-1,j+1,k+1)*shxm+tmpf(3,i,j+1,k+1)*shx+tmpf(3,i+1,j+1,k+1)*shxp)*shyp)*shzp

          epx = (+(+tmpf(4,i-1,j-1,k-1)*shxm+tmpf(4,i,j-1,k-1)*shx+tmpf(4,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(4,i-1,j  ,k-1)*shxm+tmpf(4,i,j  ,k-1)*shx+tmpf(4,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(4,i-1,j+1,k-1)*shxm+tmpf(4,i,j+1,k-1)*shx+tmpf(4,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(4,i-1,j-1,k  )*shxm+tmpf(4,i,j-1,k  )*shx+tmpf(4,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(4,i-1,j  ,k  )*shxm+tmpf(4,i,j  ,k  )*shx+tmpf(4,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(4,i-1,j+1,k  )*shxm+tmpf(4,i,j+1,k  )*shx+tmpf(4,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(4,i-1,j-1,k+1)*shxm+tmpf(4,i,j-1,k+1)*shx+tmpf(4,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(4,i-1,j  ,k+1)*shxm+tmpf(4,i,j  ,k+1)*shx+tmpf(4,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(4,i-1,j+1,k+1)*shxm+tmpf(4,i,j+1,k+1)*shx+tmpf(4,i+1,j+1,k+1)*shxp)*shyp)*shzp

          epy = (+(+tmpf(5,i-1,j-1,k-1)*shxm+tmpf(5,i,j-1,k-1)*shx+tmpf(5,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(5,i-1,j  ,k-1)*shxm+tmpf(5,i,j  ,k-1)*shx+tmpf(5,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(5,i-1,j+1,k-1)*shxm+tmpf(5,i,j+1,k-1)*shx+tmpf(5,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(5,i-1,j-1,k  )*shxm+tmpf(5,i,j-1,k  )*shx+tmpf(5,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(5,i-1,j  ,k  )*shxm+tmpf(5,i,j  ,k  )*shx+tmpf(5,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(5,i-1,j+1,k  )*shxm+tmpf(5,i,j+1,k  )*shx+tmpf(5,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(5,i-1,j-1,k+1)*shxm+tmpf(5,i,j-1,k+1)*shx+tmpf(5,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(5,i-1,j  ,k+1)*shxm+tmpf(5,i,j  ,k+1)*shx+tmpf(5,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(5,i-1,j+1,k+1)*shxm+tmpf(5,i,j+1,k+1)*shx+tmpf(5,i+1,j+1,k+1)*shxp)*shyp)*shzp

          epz = (+(+tmpf(6,i-1,j-1,k-1)*shxm+tmpf(6,i,j-1,k-1)*shx+tmpf(6,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(6,i-1,j  ,k-1)*shxm+tmpf(6,i,j  ,k-1)*shx+tmpf(6,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(6,i-1,j+1,k-1)*shxm+tmpf(6,i,j+1,k-1)*shx+tmpf(6,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(6,i-1,j-1,k  )*shxm+tmpf(6,i,j-1,k  )*shx+tmpf(6,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(6,i-1,j  ,k  )*shxm+tmpf(6,i,j  ,k  )*shx+tmpf(6,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(6,i-1,j+1,k  )*shxm+tmpf(6,i,j+1,k  )*shx+tmpf(6,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(6,i-1,j-1,k+1)*shxm+tmpf(6,i,j-1,k+1)*shx+tmpf(6,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(6,i-1,j  ,k+1)*shxm+tmpf(6,i,j  ,k+1)*shx+tmpf(6,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(6,i-1,j+1,k+1)*shxm+tmpf(6,i,j+1,k+1)*shx+tmpf(6,i+1,j+1,k+1)*shxp)*shyp)*shzp

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

          bpx = (+(+tmpf(1,i-1,j-1,k-1)*shxm+tmpf(1,i,j-1,k-1)*shx+tmpf(1,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(1,i-1,j  ,k-1)*shxm+tmpf(1,i,j  ,k-1)*shx+tmpf(1,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(1,i-1,j+1,k-1)*shxm+tmpf(1,i,j+1,k-1)*shx+tmpf(1,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(1,i-1,j-1,k  )*shxm+tmpf(1,i,j-1,k  )*shx+tmpf(1,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(1,i-1,j  ,k  )*shxm+tmpf(1,i,j  ,k  )*shx+tmpf(1,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(1,i-1,j+1,k  )*shxm+tmpf(1,i,j+1,k  )*shx+tmpf(1,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(1,i-1,j-1,k+1)*shxm+tmpf(1,i,j-1,k+1)*shx+tmpf(1,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(1,i-1,j  ,k+1)*shxm+tmpf(1,i,j  ,k+1)*shx+tmpf(1,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(1,i-1,j+1,k+1)*shxm+tmpf(1,i,j+1,k+1)*shx+tmpf(1,i+1,j+1,k+1)*shxp)*shyp)*shzp

          bpy = (+(+tmpf(2,i-1,j-1,k-1)*shxm+tmpf(2,i,j-1,k-1)*shx+tmpf(2,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(2,i-1,j  ,k-1)*shxm+tmpf(2,i,j  ,k-1)*shx+tmpf(2,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(2,i-1,j+1,k-1)*shxm+tmpf(2,i,j+1,k-1)*shx+tmpf(2,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(2,i-1,j-1,k  )*shxm+tmpf(2,i,j-1,k  )*shx+tmpf(2,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(2,i-1,j  ,k  )*shxm+tmpf(2,i,j  ,k  )*shx+tmpf(2,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(2,i-1,j+1,k  )*shxm+tmpf(2,i,j+1,k  )*shx+tmpf(2,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(2,i-1,j-1,k+1)*shxm+tmpf(2,i,j-1,k+1)*shx+tmpf(2,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(2,i-1,j  ,k+1)*shxm+tmpf(2,i,j  ,k+1)*shx+tmpf(2,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(2,i-1,j+1,k+1)*shxm+tmpf(2,i,j+1,k+1)*shx+tmpf(2,i+1,j+1,k+1)*shxp)*shyp)*shzp

          bpz = (+(+tmpf(3,i-1,j-1,k-1)*shxm+tmpf(3,i,j-1,k-1)*shx+tmpf(3,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(3,i-1,j  ,k-1)*shxm+tmpf(3,i,j  ,k-1)*shx+tmpf(3,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(3,i-1,j+1,k-1)*shxm+tmpf(3,i,j+1,k-1)*shx+tmpf(3,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(3,i-1,j-1,k  )*shxm+tmpf(3,i,j-1,k  )*shx+tmpf(3,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(3,i-1,j  ,k  )*shxm+tmpf(3,i,j  ,k  )*shx+tmpf(3,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(3,i-1,j+1,k  )*shxm+tmpf(3,i,j+1,k  )*shx+tmpf(3,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(3,i-1,j-1,k+1)*shxm+tmpf(3,i,j-1,k+1)*shx+tmpf(3,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(3,i-1,j  ,k+1)*shxm+tmpf(3,i,j  ,k+1)*shx+tmpf(3,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(3,i-1,j+1,k+1)*shxm+tmpf(3,i,j+1,k+1)*shx+tmpf(3,i+1,j+1,k+1)*shxp)*shyp)*shzp

          epx = (+(+tmpf(4,i-1,j-1,k-1)*shxm+tmpf(4,i,j-1,k-1)*shx+tmpf(4,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(4,i-1,j  ,k-1)*shxm+tmpf(4,i,j  ,k-1)*shx+tmpf(4,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(4,i-1,j+1,k-1)*shxm+tmpf(4,i,j+1,k-1)*shx+tmpf(4,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(4,i-1,j-1,k  )*shxm+tmpf(4,i,j-1,k  )*shx+tmpf(4,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(4,i-1,j  ,k  )*shxm+tmpf(4,i,j  ,k  )*shx+tmpf(4,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(4,i-1,j+1,k  )*shxm+tmpf(4,i,j+1,k  )*shx+tmpf(4,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(4,i-1,j-1,k+1)*shxm+tmpf(4,i,j-1,k+1)*shx+tmpf(4,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(4,i-1,j  ,k+1)*shxm+tmpf(4,i,j  ,k+1)*shx+tmpf(4,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(4,i-1,j+1,k+1)*shxm+tmpf(4,i,j+1,k+1)*shx+tmpf(4,i+1,j+1,k+1)*shxp)*shyp)*shzp

          epy = (+(+tmpf(5,i-1,j-1,k-1)*shxm+tmpf(5,i,j-1,k-1)*shx+tmpf(5,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(5,i-1,j  ,k-1)*shxm+tmpf(5,i,j  ,k-1)*shx+tmpf(5,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(5,i-1,j+1,k-1)*shxm+tmpf(5,i,j+1,k-1)*shx+tmpf(5,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(5,i-1,j-1,k  )*shxm+tmpf(5,i,j-1,k  )*shx+tmpf(5,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(5,i-1,j  ,k  )*shxm+tmpf(5,i,j  ,k  )*shx+tmpf(5,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(5,i-1,j+1,k  )*shxm+tmpf(5,i,j+1,k  )*shx+tmpf(5,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(5,i-1,j-1,k+1)*shxm+tmpf(5,i,j-1,k+1)*shx+tmpf(5,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(5,i-1,j  ,k+1)*shxm+tmpf(5,i,j  ,k+1)*shx+tmpf(5,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(5,i-1,j+1,k+1)*shxm+tmpf(5,i,j+1,k+1)*shx+tmpf(5,i+1,j+1,k+1)*shxp)*shyp)*shzp

          epz = (+(+tmpf(6,i-1,j-1,k-1)*shxm+tmpf(6,i,j-1,k-1)*shx+tmpf(6,i+1,j-1,k-1)*shxp)*shym       &
                 +(+tmpf(6,i-1,j  ,k-1)*shxm+tmpf(6,i,j  ,k-1)*shx+tmpf(6,i+1,j  ,k-1)*shxp)*shy        &
                 +(+tmpf(6,i-1,j+1,k-1)*shxm+tmpf(6,i,j+1,k-1)*shx+tmpf(6,i+1,j+1,k-1)*shxp)*shyp)*shzm &
               +(+(+tmpf(6,i-1,j-1,k  )*shxm+tmpf(6,i,j-1,k  )*shx+tmpf(6,i+1,j-1,k  )*shxp)*shym       &
                 +(+tmpf(6,i-1,j  ,k  )*shxm+tmpf(6,i,j  ,k  )*shx+tmpf(6,i+1,j  ,k  )*shxp)*shy        &
                 +(+tmpf(6,i-1,j+1,k  )*shxm+tmpf(6,i,j+1,k  )*shx+tmpf(6,i+1,j+1,k  )*shxp)*shyp)*shz  &
               +(+(+tmpf(6,i-1,j-1,k+1)*shxm+tmpf(6,i,j-1,k+1)*shx+tmpf(6,i+1,j-1,k+1)*shxp)*shym       &
                 +(+tmpf(6,i-1,j  ,k+1)*shxm+tmpf(6,i,j  ,k+1)*shx+tmpf(6,i+1,j  ,k+1)*shxp)*shy        &
                 +(+tmpf(6,i-1,j+1,k+1)*shxm+tmpf(6,i,j+1,k+1)*shx+tmpf(6,i+1,j+1,k+1)*shxp)*shyp)*shzp

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
