module sort

  implicit none

  private

  public :: sort__bucket


contains


  subroutine sort__bucket(up,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)

    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze, np, nsp
    integer, intent(in)    :: np2(nys:nye,nzs:nze,nsp)
    integer, intent(out)   :: cumcnt(nxgs:nxge,nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: up(6,np,nys:nye,nzs:nze,nsp)
    logical, save              :: lflag=.true.
    integer                    :: i, j, k, ii, isp
    integer                    :: cnt(nxs:nxe-1), sum_cnt(nxs:nxe-1) 
    real(8), save, allocatable :: tmp(:,:,:,:)

    if(lflag)then
       allocate(tmp(1:6,np,nys:nye,nzs:nze))
       lflag=.false.
    endif

    do isp=1,nsp

      !BUCKET SORT FOR PARTICLES IN X
!$OMP PARALLEL DO PRIVATE(ii,i,j,k,cnt,sum_cnt)
       do k=nzs,nze
       do j=nys,nye
          
          cnt(nxs:nxe-1) = 0
          tmp(1:6,1:np2(j,k,isp),j,k) = up(1:6,1:np2(j,k,isp),j,k,isp)

          do ii=1,np2(j,k,isp)
             i = int(tmp(1,ii,j,k))
             cnt(i) = cnt(i)+1
          enddo

          sum_cnt(nxs) = 0
          cumcnt(nxs,j,k,isp) = 0
          do i=nxs+1,nxe-1
             sum_cnt(i) = sum_cnt(i-1)+cnt(i-1)
             cumcnt(i,j,k,isp) = sum_cnt(i)
          enddo
          cumcnt(nxe,j,k,isp) = np2(j,k,isp)

          do ii=1,np2(j,k,isp)
             i = int(tmp(1,ii,j,k))
             up(1:6,sum_cnt(i)+1,j,k,isp) = tmp(1:6,ii,j,k)
             sum_cnt(i) = sum_cnt(i)+1
          enddo

       enddo
       enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine sort__bucket


end module sort
