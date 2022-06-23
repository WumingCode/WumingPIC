module sort

  implicit none

  private

  public :: sort__init, sort__bucket

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye


contains


  subroutine sort__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in

    is_init = .true.

  end subroutine sort__init


  subroutine sort__bucket(gp,up,cumcnt,np2,nxs,nxe)

    !BUCKET SORT FOR PARTICLES IN X
    integer, intent(in)    :: nxs, nxe
    integer, intent(in)    :: np2(nys:nye,nsp)
    integer, intent(out)   :: cumcnt(nxgs:nxge+1,nys:nye,nsp)
    real(8), intent(in)    :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out)   :: gp(ndim,np,nys:nye,nsp)
    integer                :: i, j, ii, isp
    integer                :: cnt(nxs:nxe), sum_cnt(nxs:nxe+1)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling sort__init()'
       stop
    endif

    do isp=1,nsp

!$OMP PARALLEL DO PRIVATE(ii,i,j,cnt,sum_cnt)
       do j=nys,nye

          cnt(nxs:nxe) = 0

          do ii=1,np2(j,isp)
             i = int(up(1,ii,j,isp))
             cnt(i) = cnt(i)+1
          enddo

          sum_cnt(nxs) = 0
          cumcnt(nxs,j,isp) = 0
          do i=nxs+1,nxe+1
             sum_cnt(i) = sum_cnt(i-1)+cnt(i-1)
             cumcnt(i,j,isp) = sum_cnt(i)
          enddo

          do ii=1,np2(j,isp)
             i = int(up(1,ii,j,isp))
             gp(1:ndim,sum_cnt(i)+1,j,isp) = up(1:ndim,ii,j,isp)
             sum_cnt(i) = sum_cnt(i)+1
          enddo

       enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine sort__bucket


end module sort
