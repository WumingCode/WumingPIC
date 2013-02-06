module particle

  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(gp,                                   &
                            nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2, &
                            c,q,r,delt,delx,                      &
                            up,uf)

    integer, intent(in)  :: nxgs, nxge, nys, nye, nzs, nze, np, nsp
    integer, intent(in)  :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: c, q(nsp), r(nsp), delt, delx
    real(8), intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(out) :: gp(6,np,nys:nye,nzs:nze,nsp)
    integer :: j, k, ii, isp, i0, j0, k0, ih
    real(8) :: idelx, dh, fac1, fac1r, fac2, fac2r, gam, igam, txxx
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8) :: s0xm, s0x, s0xp, s0ym, s0y, s0yp, s0zm, s0z, s0zp
    real(8) :: shxm, shx, shxp, shym, shy, shyp, shzm, shz, shzp

    idelx = 1./delx

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

!$OMP PARALLEL DO PRIVATE(ii,j,k,i0,j0,k0,ih,                        &
!$OMP                     dh,gam,igam,fac1r,fac2r,                   &
!$OMP                     s0xm,s0x,s0xp,s0ym,s0y,s0yp,s0zm,s0z,s0zp, &
!$OMP                     shxm,shx,shxp,shym,shy,shyp,shzm,shz,shzp, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,                   &
!$OMP                     uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)

             !second order shape function
             i0 = int(up(1,ii,j,k,isp)+0.5)
             dh = up(1,ii,j,k,isp)*idelx-i0
             s0xm = 0.5*(0.5-dh)*(0.5-dh)
             s0x  = 0.75-dh*dh
             s0xp = 0.5*(0.5+dh)*(0.5+dh)

             ih = int(up(1,ii,j,k,isp))
             dh = up(1,ii,j,k,isp)*idelx-0.5-ih
             shxm = 0.5*(0.5-dh)*(0.5-dh)
             shx  = 0.75-dh*dh
             shxp = 0.5*(0.5+dh)*(0.5+dh)

             j0 = int(up(2,ii,j,k,isp)+0.5)
             dh = up(2,ii,j,k,isp)*idelx-j0
             s0ym = 0.5*(0.5-dh)*(0.5-dh)
             s0y  = 0.75-dh*dh
             s0yp = 0.5*(0.5+dh)*(0.5+dh)

             dh = up(2,ii,j,k,isp)*idelx-0.5-j
             shym = 0.5*(0.5-dh)*(0.5-dh)
             shy  = 0.75-dh*dh
             shyp = 0.5*(0.5+dh)*(0.5+dh)

             k0 = int(up(3,ii,j,k,isp)+0.5)
             dh = up(3,ii,j,k,isp)*idelx-k0
             s0zm = 0.5*(0.5-dh)*(0.5-dh)
             s0z  = 0.75-dh*dh
             s0zp = 0.5*(0.5+dh)*(0.5+dh)

             dh = up(3,ii,j,k,isp)*idelx-0.5-k
             shzm = 0.5*(0.5-dh)*(0.5-dh)
             shz  = 0.75-dh*dh
             shzp = 0.5*(0.5+dh)*(0.5+dh)

             bpx = (+(+uf(1,ih-1,j0-1,k0-1)*shxm+uf(1,ih,j0-1,k0-1)*shx+uf(1,ih+1,j0-1,k0-1)*shxp)*s0ym       &
                    +(+uf(1,ih-1,j0  ,k0-1)*shxm+uf(1,ih,j0  ,k0-1)*shx+uf(1,ih+1,j0  ,k0-1)*shxp)*s0y        &
                    +(+uf(1,ih-1,j0+1,k0-1)*shxm+uf(1,ih,j0+1,k0-1)*shx+uf(1,ih+1,j0+1,k0-1)*shxp)*s0yp)*s0zm &
                  +(+(+uf(1,ih-1,j0-1,k0  )*shxm+uf(1,ih,j0-1,k0  )*shx+uf(1,ih+1,j0-1,k0  )*shxp)*s0ym       &
                    +(+uf(1,ih-1,j0  ,k0  )*shxm+uf(1,ih,j0  ,k0  )*shx+uf(1,ih+1,j0  ,k0  )*shxp)*s0y        &
                    +(+uf(1,ih-1,j0+1,k0  )*shxm+uf(1,ih,j0+1,k0  )*shx+uf(1,ih+1,j0+1,k0  )*shxp)*s0yp)*s0z  &
                  +(+(+uf(1,ih-1,j0-1,k0+1)*shxm+uf(1,ih,j0-1,k0+1)*shx+uf(1,ih+1,j0-1,k0+1)*shxp)*s0ym       &
                    +(+uf(1,ih-1,j0  ,k0+1)*shxm+uf(1,ih,j0  ,k0+1)*shx+uf(1,ih+1,j0  ,k0+1)*shxp)*s0y        &
                    +(+uf(1,ih-1,j0+1,k0+1)*shxm+uf(1,ih,j0+1,k0+1)*shx+uf(1,ih+1,j0+1,k0+1)*shxp)*s0yp)*s0zp

             bpy = (+(+uf(2,i0-1,j-1 ,k0-1)*s0xm+uf(2,i0,j-1 ,k0-1)*s0x+uf(2,i0+1,j-1 ,k0-1)*s0xp)*shym       &
                    +(+uf(2,i0-1,j   ,k0-1)*s0xm+uf(2,i0,j   ,k0-1)*s0x+uf(2,i0+1,j   ,k0-1)*s0xp)*shy        &
                    +(+uf(2,i0-1,j+1 ,k0-1)*s0xm+uf(2,i0,j+1 ,k0-1)*s0x+uf(2,i0+1,j+1 ,k0-1)*s0xp)*shyp)*s0zm &
                  +(+(+uf(2,i0-1,j-1 ,k0  )*s0xm+uf(2,i0,j-1 ,k0  )*s0x+uf(2,i0+1,j-1 ,k0  )*s0xp)*shym       &
                    +(+uf(2,i0-1,j   ,k0  )*s0xm+uf(2,i0,j   ,k0  )*s0x+uf(2,i0+1,j   ,k0  )*s0xp)*shy        &
                    +(+uf(2,i0-1,j+1 ,k0  )*s0xm+uf(2,i0,j+1 ,k0  )*s0x+uf(2,i0+1,j+1 ,k0  )*s0xp)*shyp)*s0z  &
                  +(+(+uf(2,i0-1,j-1 ,k0+1)*s0xm+uf(2,i0,j-1 ,k0+1)*s0x+uf(2,i0+1,j-1 ,k0+1)*s0xp)*shym       &
                    +(+uf(2,i0-1,j   ,k0+1)*s0xm+uf(2,i0,j   ,k0+1)*s0x+uf(2,i0+1,j   ,k0+1)*s0xp)*shy        &
                    +(+uf(2,i0-1,j+1 ,k0+1)*s0xm+uf(2,i0,j+1 ,k0+1)*s0x+uf(2,i0+1,j+1 ,k0+1)*s0xp)*shyp)*s0zp

             bpz = (+(+uf(3,i0-1,j0-1,k-1 )*s0xm+uf(3,i0,j0-1,k-1 )*s0x+uf(3,i0+1,j0-1,k-1 )*s0xp)*s0ym       &
                    +(+uf(3,i0-1,j0  ,k-1 )*s0xm+uf(3,i0,j0  ,k-1 )*s0x+uf(3,i0+1,j0  ,k-1 )*s0xp)*s0y        &
                    +(+uf(3,i0-1,j0+1,k-1 )*s0xm+uf(3,i0,j0+1,k-1 )*s0x+uf(3,i0+1,j0+1,k-1 )*s0xp)*s0yp)*shzm &
                  +(+(+uf(3,i0-1,j0-1,k   )*s0xm+uf(3,i0,j0-1,k   )*s0x+uf(3,i0+1,j0-1,k   )*s0xp)*s0ym       &
                    +(+uf(3,i0-1,j0  ,k   )*s0xm+uf(3,i0,j0  ,k   )*s0x+uf(3,i0+1,j0  ,k   )*s0xp)*s0y        &
                    +(+uf(3,i0-1,j0+1,k   )*s0xm+uf(3,i0,j0+1,k   )*s0x+uf(3,i0+1,j0+1,k   )*s0xp)*s0yp)*shz  &
                  +(+(+uf(3,i0-1,j0-1,k+1 )*s0xm+uf(3,i0,j0-1,k+1 )*s0x+uf(3,i0+1,j0-1,k+1 )*s0xp)*s0ym       &
                    +(+uf(3,i0-1,j0  ,k+1 )*s0xm+uf(3,i0,j0  ,k+1 )*s0x+uf(3,i0+1,j0  ,k+1 )*s0xp)*s0y        &
                    +(+uf(3,i0-1,j0+1,k+1 )*s0xm+uf(3,i0,j0+1,k+1 )*s0x+uf(3,i0+1,j0+1,k+1 )*s0xp)*s0yp)*shzp

             epx = (+(+uf(4,i0-1,j-1 ,k-1 )*s0xm+uf(4,i0,j-1 ,k-1 )*s0x+uf(4,i0+1,j-1 ,k-1 )*s0xp)*shym       &
                    +(+uf(4,i0-1,j   ,k-1 )*s0xm+uf(4,i0,j   ,k-1 )*s0x+uf(4,i0+1,j   ,k-1 )*s0xp)*shy        &
                    +(+uf(4,i0-1,j+1 ,k-1 )*s0xm+uf(4,i0,j+1 ,k-1 )*s0x+uf(4,i0+1,j+1 ,k-1 )*s0xp)*shyp)*shzm &
                  +(+(+uf(4,i0-1,j-1 ,k   )*s0xm+uf(4,i0,j-1 ,k   )*s0x+uf(4,i0+1,j-1 ,k   )*s0xp)*shym       &
                    +(+uf(4,i0-1,j   ,k   )*s0xm+uf(4,i0,j   ,k   )*s0x+uf(4,i0+1,j   ,k   )*s0xp)*shy        &
                    +(+uf(4,i0-1,j+1 ,k   )*s0xm+uf(4,i0,j+1 ,k   )*s0x+uf(4,i0+1,j+1 ,k   )*s0xp)*shyp)*shz  &
                  +(+(+uf(4,i0-1,j-1 ,k+1 )*s0xm+uf(4,i0,j-1 ,k+1 )*s0x+uf(4,i0+1,j-1 ,k+1 )*s0xp)*shym       &
                    +(+uf(4,i0-1,j   ,k+1 )*s0xm+uf(4,i0,j   ,k+1 )*s0x+uf(4,i0+1,j   ,k+1 )*s0xp)*shy        &
                    +(+uf(4,i0-1,j+1 ,k+1 )*s0xm+uf(4,i0,j+1 ,k+1 )*s0x+uf(4,i0+1,j+1 ,k+1 )*s0xp)*shyp)*shzp

             epy = (+(+uf(5,ih-1,j0-1,k-1 )*shxm+uf(5,ih,j0-1,k-1 )*shx+uf(5,ih+1,j0-1,k-1 )*shxp)*s0ym       &
                    +(+uf(5,ih-1,j0  ,k-1 )*shxm+uf(5,ih,j0  ,k-1 )*shx+uf(5,ih+1,j0  ,k-1 )*shxp)*s0y        &
                    +(+uf(5,ih-1,j0+1,k-1 )*shxm+uf(5,ih,j0+1,k-1 )*shx+uf(5,ih+1,j0+1,k-1 )*shxp)*s0yp)*shzm &
                  +(+(+uf(5,ih-1,j0-1,k   )*shxm+uf(5,ih,j0-1,k   )*shx+uf(5,ih+1,j0-1,k   )*shxp)*s0ym       &
                    +(+uf(5,ih-1,j0  ,k   )*shxm+uf(5,ih,j0  ,k   )*shx+uf(5,ih+1,j0  ,k   )*shxp)*s0y        &
                    +(+uf(5,ih-1,j0+1,k   )*shxm+uf(5,ih,j0+1,k   )*shx+uf(5,ih+1,j0+1,k   )*shxp)*s0yp)*shz  &
                  +(+(+uf(5,ih-1,j0-1,k+1 )*shxm+uf(5,ih,j0-1,k+1 )*shx+uf(5,ih+1,j0-1,k+1 )*shxp)*s0ym       &
                    +(+uf(5,ih-1,j0  ,k+1 )*shxm+uf(5,ih,j0  ,k+1 )*shx+uf(5,ih+1,j0  ,k+1 )*shxp)*s0y        &
                    +(+uf(5,ih-1,j0+1,k+1 )*shxm+uf(5,ih,j0+1,k+1 )*shx+uf(5,ih+1,j0+1,k+1 )*shxp)*s0yp)*shzp

             epz = (+(+uf(6,ih-1,j-1 ,k0-1)*shxm+uf(6,ih,j-1 ,k0-1)*shx+uf(6,ih+1,j-1 ,k0-1)*shxp)*shym       &
                    +(+uf(6,ih-1,j   ,k0-1)*shxm+uf(6,ih,j   ,k0-1)*shx+uf(6,ih+1,j   ,k0-1)*shxp)*shy        &
                    +(+uf(6,ih-1,j+1 ,k0-1)*shxm+uf(6,ih,j+1 ,k0-1)*shx+uf(6,ih+1,j+1 ,k0-1)*shxp)*shyp)*s0zm &
                  +(+(+uf(6,ih-1,j-1 ,k0  )*shxm+uf(6,ih,j-1 ,k0  )*shx+uf(6,ih+1,j-1 ,k0  )*shxp)*shym       &
                    +(+uf(6,ih-1,j   ,k0  )*shxm+uf(6,ih,j   ,k0  )*shx+uf(6,ih+1,j   ,k0  )*shxp)*shy        &
                    +(+uf(6,ih-1,j+1 ,k0  )*shxm+uf(6,ih,j+1 ,k0  )*shx+uf(6,ih+1,j+1 ,k0  )*shxp)*shyp)*s0z  &
                  +(+(+uf(6,ih-1,j-1 ,k0+1)*shxm+uf(6,ih,j-1 ,k0+1)*shx+uf(6,ih+1,j-1 ,k0+1)*shxp)*shym       &
                    +(+uf(6,ih-1,j   ,k0+1)*shxm+uf(6,ih,j   ,k0+1)*shx+uf(6,ih+1,j   ,k0+1)*shxp)*shy        &
                    +(+uf(6,ih-1,j+1 ,k0+1)*shxm+uf(6,ih,j+1 ,k0+1)*shx+uf(6,ih+1,j+1 ,k0+1)*shxp)*shyp)*s0zp

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
!$OMP END PARALLEL DO

    enddo

  end subroutine particle__solv


end module particle
