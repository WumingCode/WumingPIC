module particle

  implicit none

  private

  public :: particle__init, particle__solv, particle__solv_vay

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
  real(8), save :: delx, delt, c, d_delx
  real(8), allocatable :: q(:), r(:)


contains


  subroutine particle__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nzgs_in,nzge_in,nys_in,nye_in,nzs_in,nze_in, &
                            delx_in,delt_in,c_in,q_in,r_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nzgs_in, nzge_in, nys_in, nye_in, nzs_in, nze_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nzgs  = nzgs_in
    nzge  = nzge_in
    nys   = nys_in
    nye   = nye_in
    nzs   = nzs_in
    nze   = nze_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in
    d_delx = 1d0/delx

    is_init = .true.

  end subroutine particle__init


  subroutine particle__solv(gp,up,uf,cumcnt,nxs,nxe)

    integer, intent(in)  :: nxs, nxe
    integer, intent(in)  :: cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(out) :: gp(ndim,np,nys:nye,nzs:nze,nsp)

    integer            :: i, j, k, ii, isp
    real(8)            :: tmpf(6,nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: dh, fac1, fac1r, fac2, fac2r, gam, igam, txxx
    real(8)            :: bpx, bpy, bpz, epx, epy, epz
    real(8)            :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8)            :: shxm, shx, shxp, shym, shy, shyp, shzm, shz, shzp

    if(.not.is_init)then
      write(6,*)'Initialize first by calling particle__init()'
      stop
    endif

    !FIELDS AT (I+1/2, J+1/2, K+1/2)
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs-1,nze+1
      do j=nys-1,nye+1
        do i=nxs-1,nxe+1
          tmpf(1,i,j,k) = 2.5d-1*(+uf(1,i,j,k  )+uf(1,i,j+1,k  ) &
                                  +uf(1,i,j,k+1)+uf(1,i,j+1,k+1))
          tmpf(2,i,j,k) = 2.5d-1*(+uf(2,i,j,k  )+uf(2,i+1,j,k  ) &
                                  +uf(2,i,j,k+1)+uf(2,i+1,j,k+1))
          tmpf(3,i,j,k) = 2.5d-1*(+uf(3,i,j  ,k)+uf(3,i+1,j  ,k) &
                                  +uf(3,i,j+1,k)+uf(3,i+1,j+1,k))
          tmpf(4,i,j,k) = 5d-1*(+uf(4,i,j,k)+uf(4,i+1,j  ,k  ))
          tmpf(5,i,j,k) = 5d-1*(+uf(5,i,j,k)+uf(5,i  ,j+1,k  ))
          tmpf(6,i,j,k) = 5d-1*(+uf(6,i,j,k)+uf(6,i  ,j  ,k+1))
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,i,j,k,isp,                              &
!$OMP                     dh,gam,igam,fac1,fac2,txxx,fac1r,fac2r,    &
!$OMP                     shxm,shx,shxp,shym,shy,shyp,shzm,shz,shzp, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,                   &
!$OMP                     uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe

      do isp=1,nsp

        fac1 = q(isp)/r(isp)*5d-1*delt
        txxx = fac1*fac1
        fac2 = q(isp)*delt/r(isp)

        do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)

            !second order shape function
            dh = up(1,ii,j,k,isp)*d_delx-5d-1-i
            shxm = 5d-1*(5d-1-dh)*(5d-1-dh)
            shx  = 7.5d-1-dh*dh
            shxp = 5d-1*(5d-1+dh)*(5d-1+dh)

            dh = up(2,ii,j,k,isp)*d_delx-5d-1-j
            shym = 5d-1*(5d-1-dh)*(5d-1-dh)
            shy  = 7.5d-1-dh*dh
            shyp = 5d-1*(5d-1+dh)*(5d-1+dh)

            dh = up(3,ii,j,k,isp)*d_delx-5d-1-k
            shzm = 5d-1*(5d-1-dh)*(5d-1-dh)
            shz  = 7.5d-1-dh*dh
            shzp = 5d-1*(5d-1+dh)*(5d-1+dh)

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
            gam = sqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
            igam = 1d0/gam
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
            gam = 1d0/sqrt(1d0+(+gp(4,ii,j,k,isp)*gp(4,ii,j,k,isp) &
                                +gp(5,ii,j,k,isp)*gp(5,ii,j,k,isp) &
                                +gp(6,ii,j,k,isp)*gp(6,ii,j,k,isp))/(c*c))

            gp(1,ii,j,k,isp) = up(1,ii,j,k,isp)+gp(4,ii,j,k,isp)*delt*gam
            gp(2,ii,j,k,isp) = up(2,ii,j,k,isp)+gp(5,ii,j,k,isp)*delt*gam
            gp(3,ii,j,k,isp) = up(3,ii,j,k,isp)+gp(6,ii,j,k,isp)*delt*gam
          enddo

      enddo

    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    if (ndim == 7) then
!$OMP PARALLEL WORKSHARE
    gp(7,:,:,:,:) = up(7,:,:,:,:)
!$OMP END PARALLEL WORKSHARE
    endif

  end subroutine particle__solv

  
  subroutine particle__solv_vay(gp,up,uf,cumcnt,nxs,nxe)

    ! 
    ! Vay solver [PoP 15, 056701 (2008)] originally written by Seiji Zenitani
    ! 
    integer, intent(in)  :: nxs, nxe
    integer, intent(in)  :: cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(out) :: gp(ndim,np,nys:nye,nzs:nze,nsp)

    integer            :: i, j, k, ii, isp
    real(8)            :: tmpf(6,nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: dh, fac1, fac1r, fac2, gam, gam2, txxx
    real(8)            :: taux, tauy, tauz, tau2, ua, sigma
    real(8)            :: bpx, bpy, bpz, epx, epy, epz
    real(8)            :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8)            :: shxm, shx, shxp, shym, shy, shyp, shzm, shz, shzp

    if(.not.is_init)then
      write(6,*)'Initialize first by calling particle__init()'
      stop
    endif

    !FIELDS AT (I+1/2, J+1/2, K+1/2)
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs-1,nze+1
      do j=nys-1,nye+1
        do i=nxs-1,nxe+1
          tmpf(1,i,j,k) = 2.5d-1*(+uf(1,i,j,k  )+uf(1,i,j+1,k  ) &
                                  +uf(1,i,j,k+1)+uf(1,i,j+1,k+1))
          tmpf(2,i,j,k) = 2.5d-1*(+uf(2,i,j,k  )+uf(2,i+1,j,k  ) &
                                  +uf(2,i,j,k+1)+uf(2,i+1,j,k+1))
          tmpf(3,i,j,k) = 2.5d-1*(+uf(3,i,j  ,k)+uf(3,i+1,j  ,k) &
                                  +uf(3,i,j+1,k)+uf(3,i+1,j+1,k))
          tmpf(4,i,j,k) = 5d-1*(+uf(4,i,j,k)+uf(4,i+1,j  ,k  ))
          tmpf(5,i,j,k) = 5d-1*(+uf(5,i,j,k)+uf(5,i  ,j+1,k  ))
          tmpf(6,i,j,k) = 5d-1*(+uf(6,i,j,k)+uf(6,i  ,j  ,k+1))
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,i,j,k,isp,                              &
!$OMP                     dh,gam,gam2,fac1,fac2,txxx,fac1r,          &
!$OMP                     taux,tauy,tauz,tau2,ua,sigma,              &
!$OMP                     shxm,shx,shxp,shym,shy,shyp,shzm,shz,shzp, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,                   &
!$OMP                     uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe

      do isp=1,nsp

        fac1 = q(isp)/r(isp)*5d-1*delt
        fac2 = q(isp)*delt/r(isp)

        do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)

            !second order shape function
            dh = up(1,ii,j,k,isp)*d_delx-5d-1-i
            shxm = 5d-1*(5d-1-dh)*(5d-1-dh)
            shx  = 7.5d-1-dh*dh
            shxp = 5d-1*(5d-1+dh)*(5d-1+dh)

            dh = up(2,ii,j,k,isp)*d_delx-5d-1-j
            shym = 5d-1*(5d-1-dh)*(5d-1-dh)
            shy  = 7.5d-1-dh*dh
            shyp = 5d-1*(5d-1+dh)*(5d-1+dh)

            dh = up(3,ii,j,k,isp)*d_delx-5d-1-k
            shzm = 5d-1*(5d-1-dh)*(5d-1-dh)
            shz  = 7.5d-1-dh*dh
            shzp = 5d-1*(5d-1+dh)*(5d-1+dh)

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
            uvm1 = up(4,ii,j,k,isp)
            uvm2 = up(5,ii,j,k,isp)
            uvm3 = up(6,ii,j,k,isp)

            !rotate
            gam = sqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
            fac1r = fac1/gam
            !uprime
            uvm4 = uvm1 + fac2*epx + fac1r*(+uvm2*bpz-uvm3*bpy)
            uvm5 = uvm2 + fac2*epy + fac1r*(+uvm3*bpx-uvm1*bpz)
            uvm6 = uvm3 + fac2*epz + fac1r*(+uvm1*bpy-uvm2*bpx)

            taux  = fac1*bpx/c
            tauy  = fac1*bpy/c
            tauz  = fac1*bpz/c
            tau2  = taux*taux+tauy*tauy+tauz*tauz
            ua    = (uvm4*taux+uvm5*tauy+uvm6*tauz)/c  ! u* = u'.tau/c
            sigma = 1.d0 + (uvm4*uvm4+uvm5*uvm5+uvm6*uvm6)/(c*c)-tau2
            gam2  = 0.5d0*(sigma+sqrt(sigma*sigma+4.0*(tau2+ua*ua)))
            gam   = sqrt(gam2)

            txxx  = 1.d0/(tau2+gam2)  ! s
            gp(4,ii,j,k,isp) = txxx*(gam2*uvm4+c*ua*taux+gam*(uvm5*tauz-uvm6*tauy))
            gp(5,ii,j,k,isp) = txxx*(gam2*uvm5+c*ua*tauy+gam*(uvm6*taux-uvm4*tauz))
            gp(6,ii,j,k,isp) = txxx*(gam2*uvm6+c*ua*tauz+gam*(uvm4*tauy-uvm5*taux))

            !move
            gam = 1d0/gam
            gp(1,ii,j,k,isp) = up(1,ii,j,k,isp)+gp(4,ii,j,k,isp)*delt*gam
            gp(2,ii,j,k,isp) = up(2,ii,j,k,isp)+gp(5,ii,j,k,isp)*delt*gam
            gp(3,ii,j,k,isp) = up(3,ii,j,k,isp)+gp(6,ii,j,k,isp)*delt*gam
          enddo

      enddo

    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    if (ndim == 7) then
!$OMP PARALLEL WORKSHARE
    gp(7,:,:,:,:) = up(7,:,:,:,:)
!$OMP END PARALLEL WORKSHARE
    endif

  end subroutine particle__solv_vay

end module particle
