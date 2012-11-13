module particle

  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(gp,up,uf,                     &
                            np,nsp,np2,nxgs,nxge,nys,nye, &
                            c,q,r,delt,delx)

    integer, intent(in)  :: np, nxgs, nxge, nys, nye, nsp
    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(in)  :: c, q(nsp), r(nsp), delt, delx
    real(8), intent(out) :: gp(5,np,nys:nye,nsp)
    integer :: j, ii, isp, i0, j0, ih, ip, jp
    real(8) :: idelx, fac1, fac1r, fac2, fac2r, gam, igam, txxx, bt2
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8) :: s0(-1:1,2), sh(-1:1,2), dh

    idelx = 1./delx

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

!$OMP PARALLEL DO PRIVATE(ii,j,i0,j0,ih,ip,jp,s0,sh,dh,bt2,gam,igam,fac1r,fac2r, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
       do j=nys,nye
          do ii=1,np2(j,isp)

             i0 = int(up(1,ii,j,isp)+0.5)
             j0 = int(up(2,ii,j,isp)+0.5)
             ih = int(up(1,ii,j,isp))

             !second order shape function
             dh = up(1,ii,j,isp)*idelx-i0
             s0(-1,1) = 0.5*(0.5-dh)**2
             s0( 0,1) = 0.75-dh*dh
             s0(+1,1) = 0.5*(0.5+dh)**2

             dh = up(1,ii,j,isp)*idelx-0.5-ih
             sh(-1,1) = 0.5*(0.5-dh)**2
             sh( 0,1) = 0.75-dh*dh
             sh(+1,1) = 0.5*(0.5+dh)**2

             dh = up(2,ii,j,isp)*idelx-j0
             s0(-1,2) = 0.5*(0.5-dh)**2
             s0( 0,2) = 0.75-dh*dh
             s0(+1,2) = 0.5*(0.5+dh)**2

             dh = up(2,ii,j,isp)*idelx-0.5-j
             sh(-1,2) = 0.5*(0.5-dh)**2
             sh( 0,2) = 0.75-dh*dh
             sh(+1,2) = 0.5*(0.5+dh)**2

             bpx = +(+uf(1,ih-1,j0-1)*sh(-1,1)+uf(1,ih,j0-1)*sh(0,1)+uf(1,ih+1,j0-1)*sh(+1,1))*s0(-1,2) &
                   +(+uf(1,ih-1,j0  )*sh(-1,1)+uf(1,ih,j0  )*sh(0,1)+uf(1,ih+1,j0  )*sh(+1,1))*s0( 0,2) &
                   +(+uf(1,ih-1,j0+1)*sh(-1,1)+uf(1,ih,j0+1)*sh(0,1)+uf(1,ih+1,j0+1)*sh(+1,1))*s0(+1,2)

             bpy = +(+uf(2,i0-1,j-1 )*s0(-1,1)+uf(2,i0,j-1 )*s0(0,1)+uf(2,i0+1,j-1 )*s0(+1,1))*sh(-1,2) &
                   +(+uf(2,i0-1,j   )*s0(-1,1)+uf(2,i0,j   )*s0(0,1)+uf(2,i0+1,j   )*s0(+1,1))*sh( 0,2) &
                   +(+uf(2,i0-1,j+1 )*s0(-1,1)+uf(2,i0,j+1 )*s0(0,1)+uf(2,i0+1,j+1 )*s0(+1,1))*sh(+1,2)

             bpz = +(+uf(3,i0-1,j0-1)*s0(-1,1)+uf(3,i0,j0-1)*s0(0,1)+uf(3,i0+1,j0-1)*s0(+1,1))*s0(-1,2) &
                   +(+uf(3,i0-1,j0  )*s0(-1,1)+uf(3,i0,j0  )*s0(0,1)+uf(3,i0+1,j0  )*s0(+1,1))*s0( 0,2) &
                   +(+uf(3,i0-1,j0+1)*s0(-1,1)+uf(3,i0,j0+1)*s0(0,1)+uf(3,i0+1,j0+1)*s0(+1,1))*s0(+1,2)

             epx = +(+uf(4,i0-1,j-1 )*s0(-1,1)+uf(4,i0,j-1 )*s0(0,1)+uf(4,i0+1,j-1 )*s0(+1,1))*sh(-1,2) &
                   +(+uf(4,i0-1,j   )*s0(-1,1)+uf(4,i0,j   )*s0(0,1)+uf(4,i0+1,j   )*s0(+1,1))*sh( 0,2) &
                   +(+uf(4,i0-1,j+1 )*s0(-1,1)+uf(4,i0,j+1 )*s0(0,1)+uf(4,i0+1,j+1 )*s0(+1,1))*sh(+1,2)

             epy = +(+uf(5,ih-1,j0-1)*sh(-1,1)+uf(5,ih,j0-1)*sh(0,1)+uf(5,ih+1,j0-1)*sh(+1,1))*s0(-1,2) &
                   +(+uf(5,ih-1,j0  )*sh(-1,1)+uf(5,ih,j0  )*sh(0,1)+uf(5,ih+1,j0  )*sh(+1,1))*s0( 0,2) &
                   +(+uf(5,ih-1,j0+1)*sh(-1,1)+uf(5,ih,j0+1)*sh(0,1)+uf(5,ih+1,j0+1)*sh(+1,1))*s0(+1,2)

             epz = +(+uf(6,ih-1,j-1 )*sh(-1,1)+uf(6,ih,j-1 )*sh(0,1)+uf(6,ih+1,j-1 )*sh(+1,1))*sh(-1,2) &
                   +(+uf(6,ih-1,j   )*sh(-1,1)+uf(6,ih,j   )*sh(0,1)+uf(6,ih+1,j   )*sh(+1,1))*sh( 0,2) &
                   +(+uf(6,ih-1,j+1 )*sh(-1,1)+uf(6,ih,j+1 )*sh(0,1)+uf(6,ih+1,j+1 )*sh(+1,1))*sh(+1,2)

             bt2 = bpx*bpx+bpy*bpy+bpz*bpz

             !accel.
             uvm1 = up(3,ii,j,isp)+fac1*epx
             uvm2 = up(4,ii,j,isp)+fac1*epy
             uvm3 = up(5,ii,j,isp)+fac1*epz

             !rotate
             gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
             igam = 1./gam
             fac1r = fac1*igam
             fac2r = fac2/(gam+txxx*bt2*igam)

             uvm4 = uvm1+fac1r*(+uvm2*bpz-uvm3*bpy)
             uvm5 = uvm2+fac1r*(+uvm3*bpx-uvm1*bpz)
             uvm6 = uvm3+fac1r*(+uvm1*bpy-uvm2*bpx)

             uvm1 = uvm1+fac2r*(+uvm5*bpz-uvm6*bpy)
             uvm2 = uvm2+fac2r*(+uvm6*bpx-uvm4*bpz)
             uvm3 = uvm3+fac2r*(+uvm4*bpy-uvm5*bpx)

             !accel.
             gp(3,ii,j,isp) = uvm1+fac1*epx
             gp(4,ii,j,isp) = uvm2+fac1*epy
             gp(5,ii,j,isp) = uvm3+fac1*epz

             !move
             gam = 1./dsqrt(1.0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                 +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                 +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))
             gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
             gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam
          enddo
       enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine particle__solv


end module particle
