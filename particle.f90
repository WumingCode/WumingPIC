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
    integer :: j, k, ii, isp, i0, j0, k0, ih, ip, jp, kp
    real(8) :: idelx, fac1, fac1r, fac2, fac2r, gam, igam, txxx, bt2
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8) :: s0(-1:1,3), sh(-1:1,3), dh

    idelx = 1./delx

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

!$OMP PARALLEL DO PRIVATE(ii,j,k,i0,j0,k0,ih,ip,jp,kp,       &
!$OMP                     s0,sh,dh,bt2,gam,igam,fac1r,fac2r, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,           &
!$OMP                     uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)

             i0 = int(up(1,ii,j,k,isp)+0.5)
             j0 = int(up(2,ii,j,k,isp)+0.5)
             k0 = int(up(3,ii,j,k,isp)+0.5)
             ih = int(up(1,ii,j,k,isp))

             !second order shape function
             dh = up(1,ii,j,k,isp)*idelx-i0
             s0(-1,1) = 0.5*(0.5-dh)**2
             s0( 0,1) = 0.75-dh*dh
             s0(+1,1) = 0.5*(0.5+dh)**2

             dh = up(1,ii,j,k,isp)*idelx-0.5-ih
             sh(-1,1) = 0.5*(0.5-dh)**2
             sh( 0,1) = 0.75-dh*dh
             sh(+1,1) = 0.5*(0.5+dh)**2

             dh = up(2,ii,j,k,isp)*idelx-j0
             s0(-1,2) = 0.5*(0.5-dh)**2
             s0( 0,2) = 0.75-dh*dh
             s0(+1,2) = 0.5*(0.5+dh)**2

             dh = up(2,ii,j,k,isp)*idelx-0.5-j
             sh(-1,2) = 0.5*(0.5-dh)**2
             sh( 0,2) = 0.75-dh*dh
             sh(+1,2) = 0.5*(0.5+dh)**2

             dh = up(3,ii,j,k,isp)*idelx-k0
             s0(-1,3) = 0.5*(0.5-dh)**2
             s0( 0,3) = 0.75-dh*dh
             s0(+1,3) = 0.5*(0.5+dh)**2

             dh = up(3,ii,j,k,isp)*idelx-0.5-k
             sh(-1,3) = 0.5*(0.5-dh)**2
             sh( 0,3) = 0.75-dh*dh
             sh(+1,3) = 0.5*(0.5+dh)**2

             bpx = 0.D0
             bpy = 0.D0
             bpz = 0.D0
             epx = 0.D0
             epy = 0.D0
             epz = 0.D0

             do kp=-1,1
             do jp=-1,1
             do ip=-1,1
                bpx = bpx+uf(1,ih+ip,j0+jp,k0+kp)*sh(ip,1)*s0(jp,2)*s0(kp,3)
                bpy = bpy+uf(2,i0+ip,j +jp,k0+kp)*s0(ip,1)*sh(jp,2)*s0(kp,3)
                bpz = bpz+uf(3,i0+ip,j0+jp,k +kp)*s0(ip,1)*s0(jp,2)*sh(kp,3)
                epx = epx+uf(4,i0+ip,j +jp,k +kp)*s0(ip,1)*sh(jp,2)*sh(kp,3)
                epy = epy+uf(5,ih+ip,j0+jp,k +kp)*sh(ip,1)*s0(jp,2)*sh(kp,3)
                epz = epz+uf(6,ih+ip,j +jp,k0+kp)*sh(ip,1)*sh(jp,2)*s0(kp,3)
             enddo
             enddo
             enddo

             bt2 = bpx*bpx+bpy*bpy+bpz*bpz

             !accel.
             uvm1 = up(4,ii,j,k,isp)+fac1*epx
             uvm2 = up(5,ii,j,k,isp)+fac1*epy
             uvm3 = up(6,ii,j,k,isp)+fac1*epz

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
