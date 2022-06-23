module field
  use mpi
  implicit none

  private

  public :: field__init, field__fdtd_i

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
  integer, save :: mnpr, ncomw, opsum
  integer       :: nerr
  real(8), parameter :: pi = 4d0*atan(1d0)
  real(8), save :: delx, delt, c, gfac_in, d_delx, d_delt, gfac
  real(8), save :: f1, f2, f3, f4, f5
  real(8), allocatable :: q(:), r(:)


contains


  subroutine field__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nzgs_in,nzge_in,nys_in,nye_in,nzs_in,nze_in, &
                         mnpr_in,ncomw_in,opsum_in,nerr_in,                                                                &
                         delx_in,delt_in,c_in,q_in,r_in,gfac_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nzgs_in, nzge_in, nys_in, nye_in, nzs_in, nze_in
    integer, intent(in) :: mnpr_in, ncomw_in, opsum_in, nerr_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in), gfac_in

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
    mnpr  = mnpr_in
    ncomw = ncomw_in
    opsum = opsum_in
    nerr  = nerr_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in
    gfac  = gfac_in

    f1    = c*delt/delx
    f2    = gfac*f1*f1
    f3    = 4d0*pi*delx/c
    f4    = 6d0+(delx/(c*delt*gfac))**2
    f5    = (delx/(c*delt*gfac))**2
    d_delx = 1d0/delx
    d_delt = 1d0/delt

    is_init = .true.

  end subroutine field__init


  subroutine field__fdtd_i(uf,up,gp,cumcnt,nxs,nxe, &
                           set_boundary_dfield, &
                           set_boundary_curre, &
                           set_boundary_phi)

    interface
      ! set boundary for field
      subroutine set_boundary_dfield(df,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)
        real(8), intent(inout) :: df(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
        integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, nxgs, nxge
      end subroutine set_boundary_dfield

      ! set boundary for current
      subroutine set_boundary_curre(uj,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)
        real(8), intent(inout) :: uj(3,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
        integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, nxgs, nxge
      end subroutine set_boundary_curre

      ! set boundary for potential
      subroutine set_boundary_phi(phi,nxs,nxe,nys,nye,nzs,nze,l)
        real(8), intent(inout) :: phi(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
        integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, l
      end subroutine set_boundary_phi
    end interface

    integer, intent(in)    :: nxs, nxe
    integer, intent(in)    :: cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp)
    real(8), intent(in)    :: gp(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)    :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    logical, save          :: lflag=.true.
    integer                :: i, j, k, ieq
    real(8), save, allocatable :: df(:,:,:,:), gkl(:,:,:,:), uj(:,:,:,:)

    if(.not.is_init)then
      write(6,*)'Initialize first by calling field__init()'
      stop
    endif

    if(lflag)then
      allocate(df(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
      allocate(gkl(3,nxgs:nxge,nys:nye,nzs:nze))
      allocate(uj(3,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2))
!$OMP PARALLEL WORKSHARE
      df(1:6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2) = 0d0
!$OMP END PARALLEL WORKSHARE
!$OMP PARALLEL WORKSHARE
      gkl(1:3,nxgs:nxge,nys:nye,nzs:nze) = 0d0
!$OMP END PARALLEL WORKSHARE
!$OMP PARALLEL WORKSHARE
      uj(1:3,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2) = 0d0
!$OMP END PARALLEL WORKSHARE
      lflag=.false.
    endif

    call ele_cur(uj,up,gp,cumcnt,nxs,nxe)
    call set_boundary_curre(uj,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)

!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
      gkl(1,i,j,k) = +f2*(+uf(1,i,j,k-1)+uf(1,i,j-1,k)                 &
                          +uf(1,i-1,j,k)-6d0*uf(1,i,j,k)+uf(1,i+1,j,k) &
                          +uf(1,i,j+1,k)+uf(1,i,j,k+1)                 &
                          +f3*(-uj(3,i,j-1,k)+uj(3,i,j,k)              &
                               +uj(2,i,j,k-1)-uj(2,i,j,k))             &
                         )                                             &
                     -f1*(-uf(6,i,j-1,k)+uf(6,i,j,k)                   &
                          +uf(5,i,j,k-1)-uf(5,i,j,k))
      gkl(2,i,j,k) = +f2*(+uf(2,i,j,k-1)+uf(2,i,j-1,k)                 &
                          +uf(2,i-1,j,k)-6d0*uf(2,i,j,k)+uf(2,i+1,j,k) &
                          +uf(2,i,j+1,k)+uf(2,i,j,k+1)                 &
                          +f3*(-uj(1,i,j,k-1)+uj(1,i,j,k)              &
                               +uj(3,i-1,j,k)-uj(3,i,j,k))             &
                         )                                             &
                     -f1*(-uf(4,i,j,k-1)+uf(4,i,j,k)                   &
                          +uf(6,i-1,j,k)-uf(6,i,j,k))
      gkl(3,i,j,k) = +f2*(+uf(3,i,j,k-1)+uf(3,i,j-1,k)                 &
                          +uf(3,i-1,j,k)-6d0*uf(3,i,j,k)+uf(3,i+1,j,k) &
                          +uf(3,i,j+1,k)+uf(3,i,j,k+1)                 &
                          +f3*(-uj(2,i-1,j,k)+uj(2,i,j,k)              &
                               +uj(1,i,j-1,k)-uj(1,i,j,k))             &
                         )                                             &
                     -f1*(-uf(5,i-1,j,k)+uf(5,i,j,k)                   &
                          +uf(4,i,j-1,k)-uf(4,i,j,k))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !solve  < bx, by & bz >
    call cgm(df,gkl,nxs,nxe,set_boundary_phi)

    call set_boundary_dfield(df,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)

    !solve  < ex, ey & ez >
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe
      df(4,i,j,k) = +f1*(+gfac*(-df(3,i,j,k)+df(3,i,j+1,k)    &
                                +df(2,i,j,k)-df(2,i,j,k+1))   &
                         +     (-uf(3,i,j,k)+uf(3,i,j+1,k)    &
                                +uf(2,i,j,k)-uf(2,i,j,k+1)) ) &
                    -4d0*pi*delt*uj(1,i,j,k)

      df(5,i,j,k) = +f1*(+gfac*(-df(1,i,j,k)+df(1,i,j,k+1)    &
                                +df(3,i,j,k)-df(3,i+1,j,k))   &
                         +     (-uf(1,i,j,k)+uf(1,i,j,k+1)    &
                                +uf(3,i,j,k)-uf(3,i+1,j,k)) ) &
                    -4d0*pi*delt*uj(2,i,j,k)

      df(6,i,j,k) = +f1*(+gfac*(-df(2,i,j,k)+df(2,i+1,j,k)    &
                                +df(1,i,j,k)-df(1,i,j+1,k))   &
                         +     (-uf(2,i,j,k)+uf(2,i+1,j,k)    &
                                +uf(1,i,j,k)-uf(1,i,j+1,k)) ) &
                    -4d0*pi*delt*uj(3,i,j,k)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    call set_boundary_dfield(df,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)

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


  subroutine ele_cur(uj,up,gp,cumcnt,nxs,nxe)

    integer, intent(in)  :: nxs, nxe
    integer, intent(in)  :: cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: gp(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(out) :: uj(3,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)

    integer            :: ii, i, j, k, isp, i2, inc, ip, jp, kp
    real(8), parameter :: fac = 1d0/3d0
    real(8)            :: dh
    real(8)            :: s1_1, s1_2, s1_3, smo_1, smo_2, smo_3
    real(8)            :: pjtmp, dstmp, qdxdt
    real(8)            :: s0x(-2:2), s0y(-2:2), s0z(-2:2), dsx(-2:2), dsy(-2:2), dsz(-2:2)
    real(8)            :: pjx(-2:2,-2:2,-2:2), pjy(-2:2,-2:2,-2:2), pjz(-2:2,-2:2,-2:2)

!$OMP PARALLEL WORKSHARE
    uj(1:3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2) = 0d0
!$OMP END PARALLEL WORKSHARE

    !--------------Charge Conservation Method -------------!
    !---- Density Decomposition (Esirkepov, CPC, 2001) ----!
!$OMP PARALLEL DO PRIVATE(ii,i,j,k,isp,i2,inc,ip,jp,kp,              &
!$OMP                     s1_1,s1_2,s1_3,smo_1,smo_2,smo_3,dh,       &
!$OMP                     s0x,s0y,s0z,dsx,dsy,dsz,pjx,pjy,pjz,qdxdt, &
!$OMP                     pjtmp,dstmp) &
!$OMP REDUCTION(+:uj)
    do k=nzs,nze
    do j=nys,nye
    do i=nxs,nxe

      pjx(-2:2,-2:2,-2:2) = 0d0
      pjy(-2:2,-2:2,-2:2) = 0d0
      pjz(-2:2,-2:2,-2:2) = 0d0

      do isp=1,nsp

        qdxdt = q(isp)*delx*d_delt

!OCL LOOP_FISSION_TARGET(LS)
!OCL SWP_FREG_RATE(20)
        do ii=cumcnt(i,j,k,isp)+1,cumcnt(i+1,j,k,isp)

          !second order shape function
          dh = up(1,ii,j,k,isp)*d_delx-5d-1-i
          s0x(-2) = 0d0
          s0x(-1) = 5d-1*(5d-1-dh)*(5d-1-dh)
          s0x( 0) = 7.5d-1-dh*dh
          s0x(+1) = 5d-1*(5d-1+dh)*(5d-1+dh)
          s0x(+2) = 0d0

          dh = up(2,ii,j,k,isp)*d_delx-5d-1-j
          s0y(-2) = 0d0
          s0y(-1) = 5d-1*(5d-1-dh)*(5d-1-dh)
          s0y( 0) = 7.5d-1-dh*dh
          s0y(+1) = 5d-1*(5d-1+dh)*(5d-1+dh)
          s0y(+2) = 0d0

          dh = up(3,ii,j,k,isp)*d_delx-5d-1-k
          s0z(-2) = 0d0
          s0z(-1) = 5d-1*(5d-1-dh)*(5d-1-dh)
          s0z( 0) = 7.5d-1-dh*dh
          s0z(+1) = 5d-1*(5d-1+dh)*(5d-1+dh)
          s0z(+2) = 0d0

!OCL FISSION_POINT(1)

          i2 = int(gp(1,ii,j,k,isp)*d_delx)
          dh = gp(1,ii,j,k,isp)*d_delx-5d-1-i2
          inc = i2-i
          s1_1 = 5d-1*(5d-1-dh)*(5d-1-dh)
          s1_2 = 7.5d-1-dh*dh
          s1_3 = 5d-1*(5d-1+dh)*(5d-1+dh)
          smo_1 = -(inc-abs(inc))*5d-1+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*5d-1+0
          dsx(-2) = s1_1*smo_1
          dsx(-1) = s1_1*smo_2+s1_2*smo_1
          dsx( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsx(+1) = s1_3*smo_2+s1_2*smo_3
          dsx(+2) = s1_3*smo_3

          i2 = int(gp(2,ii,j,k,isp)*d_delx)
          dh = gp(2,ii,j,k,isp)*d_delx-5d-1-i2
          inc = i2-j
          s1_1 = 5d-1*(5d-1-dh)*(5d-1-dh)
          s1_2 = 7.5d-1-dh*dh
          s1_3 = 5d-1*(5d-1+dh)*(5d-1+dh)
          smo_1 = -(inc-abs(inc))*5d-1+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*5d-1+0
          dsy(-2) = s1_1*smo_1
          dsy(-1) = s1_1*smo_2+s1_2*smo_1
          dsy( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsy(+1) = s1_3*smo_2+s1_2*smo_3
          dsy(+2) = s1_3*smo_3

          i2 = int(gp(3,ii,j,k,isp)*d_delx)
          dh = gp(3,ii,j,k,isp)*d_delx-5d-1-i2
          inc = i2-k
          s1_1 = 5d-1*(5d-1-dh)*(5d-1-dh)
          s1_2 = 7.5d-1-dh*dh
          s1_3 = 5d-1*(5d-1+dh)*(5d-1+dh)
          smo_1 = -(inc-abs(inc))*5d-1+0
          smo_2 = -abs(inc)+1
          smo_3 = (inc+abs(inc))*5d-1+0
          dsz(-2) = s1_1*smo_1
          dsz(-1) = s1_1*smo_2+s1_2*smo_1
          dsz( 0) = s1_2*smo_2+s1_3*smo_1+s1_1*smo_3
          dsz(+1) = s1_3*smo_2+s1_2*smo_3
          dsz(+2) = s1_3*smo_3

          dsx(-2:2) = dsx(-2:2)-s0x(-2:2)
          dsy(-2:2) = dsy(-2:2)-s0y(-2:2)
          dsz(-2:2) = dsz(-2:2)-s0z(-2:2)

!OCL FISSION_POINT(1)

!OCL UNROLL('FULL')
          do kp=-2,2
!OCL FULLUNROLL_PRE_SIMD
          do jp=-2,2

            !pjx
            dstmp = ( (s0y(jp)+5d-1*dsy(jp))*s0z(kp) &
                     +(5d-1*s0y(jp)+fac*dsy(jp))*dsz(kp))*qdxdt

            pjtmp =      -dsx(-2)*dstmp
            pjx(-1,jp,kp) = pjx(-1,jp,kp)+pjtmp

            pjtmp = pjtmp-dsx(-1)*dstmp
            pjx( 0,jp,kp) = pjx( 0,jp,kp)+pjtmp

            pjtmp = pjtmp-dsx( 0)*dstmp
            pjx(+1,jp,kp) = pjx(+1,jp,kp)+pjtmp

            pjtmp = pjtmp-dsx(+1)*dstmp
            pjx(+2,jp,kp) = pjx(+2,jp,kp)+pjtmp

            !pjy
            dstmp = ( (s0x(jp)+5d-1*dsx(jp))*s0z(kp) &
                     +(5d-1*s0x(jp)+fac*dsx(jp))*dsz(kp))*qdxdt

            pjtmp =      -dsy(-2)*dstmp
            pjy(-1,jp,kp) = pjy(-1,jp,kp)+pjtmp

            pjtmp = pjtmp-dsy(-1)*dstmp
            pjy( 0,jp,kp) = pjy( 0,jp,kp)+pjtmp

            pjtmp = pjtmp-dsy( 0)*dstmp
            pjy(+1,jp,kp) = pjy(+1,jp,kp)+pjtmp

            pjtmp = pjtmp-dsy(+1)*dstmp
            pjy(+2,jp,kp) = pjy(+2,jp,kp)+pjtmp

            !pjz
            dstmp = ( (s0x(jp)+5d-1*dsx(jp))*s0y(kp) &
                     +(5d-1*s0x(jp)+fac*dsx(jp))*dsy(kp))*qdxdt

            pjtmp =      -dsz(-2)*dstmp
            pjz(-1,jp,kp) = pjz(-1,jp,kp)+pjtmp

            pjtmp = pjtmp-dsz(-1)*dstmp
            pjz( 0,jp,kp) = pjz( 0,jp,kp)+pjtmp

            pjtmp = pjtmp-dsz( 0)*dstmp
            pjz(+1,jp,kp) = pjz(+1,jp,kp)+pjtmp

            pjtmp = pjtmp-dsz(+1)*dstmp
            pjz(+2,jp,kp) = pjz(+2,jp,kp)+pjtmp
          enddo
          enddo

        enddo

      enddo

!OCL UNROLL('FULL')
      do kp=-2,2
!OCL UNROLL('FULL')
      do jp=-2,2
!OCL UNROLL('FULL')
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


  subroutine cgm(db,gkl,nxs,nxe,set_boundary_phi)

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method
    !  #  this routine will be stopped after iteration number = ite_max
    !-----------------------------------------------------------------------

    interface
      ! set boundary for potential
      subroutine set_boundary_phi(phi,nxs,nxe,nys,nye,nzs,nze,l)
        real(8), intent(inout) :: phi(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
        integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, l
      end subroutine set_boundary_phi
    end interface

    integer, intent(in)    :: nxs, nxe
    real(8), intent(in)    :: gkl(3,nxgs:nxge,nys:nye,nzs:nze)
    real(8), intent(inout) :: db(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    integer, parameter :: ite_max = 100 ! maximum number of interation
    integer            :: i, j, k, l, ite
    real(8), parameter :: err = 1d-6
    real(8)            :: eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g
    real(8)            :: phi(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1), p(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    real(8)            :: r(nxs:nxe,nys:nye,nzs:nze), b(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: ap(nxs:nxe,nys:nye,nzs:nze)
    real(8)            :: bff_snd(2),bff_rcv(2)

    do l=1,3

      ! initial guess
      ite = 0
      sum = 0d0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sum)
      do k=nzs,nze
      do j=nys,nye
      do i=nxs,nxe
        phi(i,j,k) = db(l,i,j,k)
        b(i,j,k) = f5*gkl(l,i,j,k)
        sum = sum+b(i,j,k)*b(i,j,k)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

      eps = sqrt(sum_g)*err

      call set_boundary_phi(phi,nxs,nxe,nys,nye,nzs,nze,l)

      sumr = 0d0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sumr)
      do k=nzs,nze
      do j=nys,nye
      do i=nxs,nxe
        r(i,j,k) = b(i,j,k)+phi(i,j,k-1)+phi(i,j-1,k)               &
                           +phi(i-1,j,k)-f4*phi(i,j,k)+phi(i+1,j,k) &
                           +phi(i,j+1,k)+phi(i,j,k+1)
        p(i,j,k) = r(i,j,k)
        sumr = sumr+r(i,j,k)*r(i,j,k)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

      if(sqrt(sumr_g) > eps)then

        do while(sum_g > eps)

          ite = ite+1

          call set_boundary_phi(p,nxs,nxe,nys,nye,nzs,nze,l)

          sumr = 0d0
          sum2 = 0d0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sumr,sum2)
          do k=nzs,nze
          do j=nys,nye
          do i=nxs,nxe
            ap(i,j,k) = -p(i,j,k-1)-p(i,j-1,k)             &
                        -p(i-1,j,k)+f4*p(i,j,k)-p(i+1,j,k) &
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
          do i=nxs,nxe
            phi(i,j,k) = phi(i,j,k)+av*p(i,j,k)
            r(i,j,k) = r(i,j,k)-av*ap(i,j,k)
          enddo
          enddo
          enddo
!$OMP END PARALLEL DO

          sum_g = sqrt(sumr_g)

          if(ite >= ite_max) then
            write(6,*)'********** stop at cgm after ite_max **********'
            stop
          endif

          sum1 = 0d0
!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:sum1)
          do k=nzs,nze
          do j=nys,nye
          do i=nxs,nxe
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
          do i=nxs,nxe
            p(i,j,k) = r(i,j,k)+bv*p(i,j,k)
          enddo
          enddo
          enddo
!$OMP END PARALLEL DO

        enddo
      endif

!$OMP PARALLEL WORKSHARE
      db(l,nxs:nxe,nys:nye,nzs:nze) = phi(nxs:nxe,nys:nye,nzs:nze)
!$OMP END PARALLEL WORKSHARE

    enddo

  end subroutine cgm


end module field
