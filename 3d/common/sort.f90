module sort

  implicit none

  private

  public :: sort__init, sort__bucket

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze


contains


  subroutine sort__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nzgs_in,nzge_in,nys_in,nye_in,nzs_in,nze_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nzgs_in, nzge_in, nys_in, nye_in, nzs_in, nze_in

    ndim = ndim_in
    np   = np_in
    nsp  = nsp_in
    nxgs = nxgs_in
    nxge = nxge_in
    nygs = nygs_in
    nyge = nyge_in
    nzgs = nzgs_in
    nzge = nzge_in
    nys  = nys_in
    nye  = nye_in
    nzs  = nzs_in
    nze  = nze_in

    is_init = .true.

  end subroutine sort__init


  subroutine sort__bucket(gp,up,cumcnt,np2,nxs,nxe)

    integer, intent(in)  :: nxs, nxe
    integer, intent(in)  :: np2(nys:nye,nzs:nze,nsp)
    integer, intent(out) :: cumcnt(nxgs:nxge+1,nys:nye,nzs:nze,nsp)
    real(8), intent(in)  :: up(ndim,np,nys:nye,nzs:nze,nsp)
    real(8), intent(out) :: gp(ndim,np,nys:nye,nzs:nze,nsp)
    integer              :: i, j, k, ii, isp
    integer              :: cnt(nxs:nxe), sum_cnt(nxs:nxe+1)

    if(.not.is_init)then
      write(6,*)'Initialize first by calling sort__init()'
      stop
    endif

    do isp=1,nsp

      !BUCKET SORT FOR PARTICLES IN X
!$OMP PARALLEL DO PRIVATE(ii,i,j,k,cnt,sum_cnt)
      do k=nzs,nze
      do j=nys,nye

        cnt(nxs:nxe) = 0

        do ii=1,np2(j,k,isp)
          i = int(up(1,ii,j,k,isp))
          cnt(i) = cnt(i)+1
        enddo

        sum_cnt(nxs) = 0
        cumcnt(nxs,j,k,isp) = 0
        do i=nxs+1,nxe+1
          sum_cnt(i) = sum_cnt(i-1)+cnt(i-1)
          cumcnt(i,j,k,isp) = sum_cnt(i)
        enddo

        do ii=1,np2(j,k,isp)
          i = int(up(1,ii,j,k,isp))
          gp(1:ndim,sum_cnt(i)+1,j,k,isp) = up(1:ndim,ii,j,k,isp)
          sum_cnt(i) = sum_cnt(i)+1
        enddo

      enddo
      enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine sort__bucket


end module sort
