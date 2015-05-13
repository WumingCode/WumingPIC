module boundary

  implicit none

  private

  public :: boundary__field
  public :: boundary__particle_x
  public :: boundary__particle_yz
  public :: boundary__curre
  public :: boundary__phi
  public :: boundary__mom


contains


  subroutine boundary__particle_x(gp,                                 &
                                  nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)

    integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, np, nsp
    integer, intent(inout) :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: gp(6,np,nys:nye,nzs:nze,nsp)
    integer                :: j, k, ii, isp, ipos

    do isp=1,nsp

!$OMP PARALLEL DO PRIVATE(ii,j,k,ipos)
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)

             ipos = int(gp(1,ii,j,k,isp))

             if(ipos <= nxs-1)then
                gp(1,ii,j,k,isp) = 2.0*nxs-gp(1,ii,j,k,isp)
                gp(4,ii,j,k,isp) = -gp(4,ii,j,k,isp)
             else if(ipos >= nxe)then
                gp(1,ii,j,k,isp) = 2.0*nxe-gp(1,ii,j,k,isp)
                gp(4,ii,j,k,isp) = -gp(4,ii,j,k,isp)
             endif

          enddo
       enddo
       enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine boundary__particle_x


  subroutine boundary__particle_yz(up,                                  &
                                   nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                                   np,nsp,np2,                          &
                                   jup,jdown,kup,kdown,nstat,mnpi,mnpr,ncomw,nerr)

!$  use omp_lib
    integer, intent(in)    :: nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
    integer, intent(in)    :: np, nsp
    integer, intent(in)    :: jup, jdown, kup, kdown, mnpi, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    integer, intent(inout) :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: up(6,np,nys:nye,nzs:nze,nsp)

    logical, save              :: lflag=.true.
!$  integer(omp_lock_kind)     :: lck(nys-1:nye+1,nzs-1:nze+1)
    integer                    :: j, k, ii, iii, isp, jpos, kpos
    integer                    :: cnt(nys-1:nye+1,nzs-1:nze+1), cnt2(nys:nye,nzs:nze)
    integer                    :: cnt_tmp_j(nzs-1:nze+1), cnt_tmp_k(nys:nye), cnt_tmp
    integer                    :: ssize, rsize
    integer                    :: bff_cnt_j(nze-nzs+3)
    integer, save, allocatable :: flag(:,:,:)
    real(8), save, allocatable :: bff_ptcl(:,:,:)
    real(8), allocatable       :: ptcl_snd(:), ptcl_rcv(:)

    if(lflag)then
       allocate(flag(np,nys:nye,nzs:nze))
       allocate(bff_ptcl(np,nys-1:nye+1,nzs-1:nze+1))
       lflag=.false.
    endif

!$OMP PARALLEL DO PRIVATE(j,k)
!$    do k=nzs-1,nze+1
!$    do j=nys-1,nye+1
!$       call omp_init_lock(lck(j,k))
!$    enddo
!$    enddo
!$OMP END PARALLEL DO

    do isp=1,nsp

!$OMP PARALLEL

!$OMP WORKSHARE
       cnt(nys-1:nye+1,nzs-1:nze+1) = 0
       cnt2(nys:nye,nzs:nze) = 0
!$OMP END WORKSHARE

!$OMP DO PRIVATE(ii,j,k,jpos,kpos)
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)

             jpos = int(up(2,ii,j,k,isp))
             kpos = int(up(3,ii,j,k,isp))

             if(.not.(jpos == j .and. kpos == k))then

                if(jpos <= nygs-1)then
                   up(2,ii,j,k,isp) = up(2,ii,j,k,isp)+(nyge-nygs+1)
                else if(jpos >= nyge+1)then
                   up(2,ii,j,k,isp) = up(2,ii,j,k,isp)-(nyge-nygs+1)
                endif

                if(kpos <= nzgs-1)then
                   up(3,ii,j,k,isp) = up(3,ii,j,k,isp)+(nzge-nzgs+1)
                else if(kpos >= nzge+1)then
                   up(3,ii,j,k,isp) = up(3,ii,j,k,isp)-(nzge-nzgs+1)
                endif

!$              call omp_set_lock(lck(jpos,kpos))
                bff_ptcl(1+6*cnt(jpos,kpos),jpos,kpos) = up(1,ii,j,k,isp)
                bff_ptcl(2+6*cnt(jpos,kpos),jpos,kpos) = up(2,ii,j,k,isp)
                bff_ptcl(3+6*cnt(jpos,kpos),jpos,kpos) = up(3,ii,j,k,isp)
                bff_ptcl(4+6*cnt(jpos,kpos),jpos,kpos) = up(4,ii,j,k,isp)
                bff_ptcl(5+6*cnt(jpos,kpos),jpos,kpos) = up(5,ii,j,k,isp)
                bff_ptcl(6+6*cnt(jpos,kpos),jpos,kpos) = up(6,ii,j,k,isp)
                cnt(jpos,kpos) = cnt(jpos,kpos)+1
!$              call omp_unset_lock(lck(jpos,kpos))

                cnt2(j,k) = cnt2(j,k)+1
                flag(cnt2(j,k),j,k) = ii
             endif
          enddo
       enddo
       enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL 

       !!TRANSFER TO JDOWN and RECEIVE FROM JUP
       bff_cnt_j(1:nze-nzs+3) = cnt(nys-1,nzs-1:nze+1)
       call MPI_SENDRECV(bff_cnt_j(1)    ,nze-nzs+3,mnpi,jdown,100, &
                         cnt_tmp_j(nzs-1),nze-nzs+3,mnpi,jup  ,100, &
                         ncomw,nstat,nerr)
       ssize = 6*sum(bff_cnt_j(1:nze-nzs+3))
       rsize = 6*sum(cnt_tmp_j(nzs-1:nze+1))
       allocate(ptcl_snd(ssize))
       allocate(ptcl_rcv(rsize))

       j=nys-1
!$OMP PARALLEL DO PRIVATE(iii,ii,k)
       do k=nzs-1,nze+1
          do ii=1,cnt(j,k)
             iii = 6*(ii+sum(cnt(j,nzs-1:k-1))-1)

             ptcl_snd(1+iii) = bff_ptcl(1+6*(ii-1),j,k)
             ptcl_snd(2+iii) = bff_ptcl(2+6*(ii-1),j,k)
             ptcl_snd(3+iii) = bff_ptcl(3+6*(ii-1),j,k)
             ptcl_snd(4+iii) = bff_ptcl(4+6*(ii-1),j,k)
             ptcl_snd(5+iii) = bff_ptcl(5+6*(ii-1),j,k)
             ptcl_snd(6+iii) = bff_ptcl(6+6*(ii-1),j,k)
          enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,jdown,101, &
                         ptcl_rcv(1),rsize,mnpr,jup  ,101, &
                         ncomw,nstat,nerr)

       j=nye
!$OMP PARALLEL DO PRIVATE(iii,ii,k)
       do k=nzs-1,nze+1
          do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_j(k)
             iii = 6*(ii-1-cnt(j,k)+sum(cnt_tmp_j(nzs-1:k-1)))

             bff_ptcl(1+6*(ii-1),j,k) = ptcl_rcv(1+iii)
             bff_ptcl(2+6*(ii-1),j,k) = ptcl_rcv(2+iii)
             bff_ptcl(3+6*(ii-1),j,k) = ptcl_rcv(3+iii)
             bff_ptcl(4+6*(ii-1),j,k) = ptcl_rcv(4+iii)
             bff_ptcl(5+6*(ii-1),j,k) = ptcl_rcv(5+iii)
             bff_ptcl(6+6*(ii-1),j,k) = ptcl_rcv(6+iii)
          enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL WORKSHARE
       cnt(nye,nzs-1:nze+1) = cnt(nye,nzs-1:nze+1)+cnt_tmp_j(nzs-1:nze+1)
!$OMP END PARALLEL WORKSHARE

       deallocate(ptcl_snd)
       deallocate(ptcl_rcv)


       !!TRANSFER TO JUP and RECEIVE FROM JDOWN
       bff_cnt_j(1:nze-nzs+3) = cnt(nye+1,nzs-1:nze+1)
       call MPI_SENDRECV(bff_cnt_j(1)    ,nze-nzs+3,mnpi,jup  ,200, &
                         cnt_tmp_j(nzs-1),nze-nzs+3,mnpi,jdown,200, &
                         ncomw,nstat,nerr)
       ssize = 6*sum(bff_cnt_j(1:nze-nzs+3))
       rsize = 6*sum(cnt_tmp_j(nzs-1:nze+1))
       allocate(ptcl_snd(ssize))
       allocate(ptcl_rcv(rsize))

       j=nye+1
!$OMP PARALLEL DO PRIVATE(iii,ii,k)
       do k=nzs-1,nze+1
          do ii=1,cnt(j,k)
             iii = 6*(ii+sum(cnt(j,nzs-1:k-1))-1)

             ptcl_snd(1+iii) = bff_ptcl(1+6*(ii-1),j,k)
             ptcl_snd(2+iii) = bff_ptcl(2+6*(ii-1),j,k)
             ptcl_snd(3+iii) = bff_ptcl(3+6*(ii-1),j,k)
             ptcl_snd(4+iii) = bff_ptcl(4+6*(ii-1),j,k)
             ptcl_snd(5+iii) = bff_ptcl(5+6*(ii-1),j,k)
             ptcl_snd(6+iii) = bff_ptcl(6+6*(ii-1),j,k)
          enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,jup  ,201, &
                         ptcl_rcv(1),rsize,mnpr,jdown,201, &
                         ncomw,nstat,nerr)

       j=nys
!$OMP PARALLEL DO PRIVATE(iii,ii,k)
       do k=nzs-1,nze+1
          do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_j(k)
             iii = 6*(ii-1-cnt(j,k)+sum(cnt_tmp_j(nzs-1:k-1)))

             bff_ptcl(1+6*(ii-1),j,k) = ptcl_rcv(1+iii)
             bff_ptcl(2+6*(ii-1),j,k) = ptcl_rcv(2+iii)
             bff_ptcl(3+6*(ii-1),j,k) = ptcl_rcv(3+iii)
             bff_ptcl(4+6*(ii-1),j,k) = ptcl_rcv(4+iii)
             bff_ptcl(5+6*(ii-1),j,k) = ptcl_rcv(5+iii)
             bff_ptcl(6+6*(ii-1),j,k) = ptcl_rcv(6+iii)
          enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL WORKSHARE
       cnt(nys,nzs-1:nze+1) = cnt(nys,nzs-1:nze+1)+cnt_tmp_j(nzs-1:nze+1)
!$OMP END PARALLEL WORKSHARE

       deallocate(ptcl_snd)
       deallocate(ptcl_rcv)


       !!TRANSFER TO KDOWN and RECEIVE FROM KUP
       call MPI_SENDRECV(cnt(nys,nzs-1),nye-nys+1,mnpi,kdown,300, &
                         cnt_tmp_k(nys),nye-nys+1,mnpi,kup  ,300, &
                         ncomw,nstat,nerr)
       ssize = 6*sum(cnt(nys:nye,nzs-1))
       rsize = 6*sum(cnt_tmp_k(nys:nye))
       allocate(ptcl_snd(ssize))
       allocate(ptcl_rcv(rsize))

       k=nzs-1
!$OMP PARALLEL DO PRIVATE(iii,ii,j)
       do j=nys,nye
          do ii=1,cnt(j,k)
             iii = 6*(ii+sum(cnt(nys:j-1,k))-1)

             ptcl_snd(1+iii) = bff_ptcl(1+6*(ii-1),j,k)
             ptcl_snd(2+iii) = bff_ptcl(2+6*(ii-1),j,k)
             ptcl_snd(3+iii) = bff_ptcl(3+6*(ii-1),j,k)
             ptcl_snd(4+iii) = bff_ptcl(4+6*(ii-1),j,k)
             ptcl_snd(5+iii) = bff_ptcl(5+6*(ii-1),j,k)
             ptcl_snd(6+iii) = bff_ptcl(6+6*(ii-1),j,k)
          enddo
       enddo
!$OMP END PARALLEL DO

       call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,kdown,301, &
                         ptcl_rcv(1),rsize,mnpr,kup  ,301, &
                         ncomw,nstat,nerr)

       k=nze
!$OMP PARALLEL DO PRIVATE(iii,ii,j)
       do j=nys,nye
          do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_k(j)
             iii = 6*(ii-1-cnt(j,k)+sum(cnt_tmp_k(nys:j-1)))

             bff_ptcl(1+6*(ii-1),j,k) = ptcl_rcv(1+iii)
             bff_ptcl(2+6*(ii-1),j,k) = ptcl_rcv(2+iii)
             bff_ptcl(3+6*(ii-1),j,k) = ptcl_rcv(3+iii)
             bff_ptcl(4+6*(ii-1),j,k) = ptcl_rcv(4+iii)
             bff_ptcl(5+6*(ii-1),j,k) = ptcl_rcv(5+iii)
             bff_ptcl(6+6*(ii-1),j,k) = ptcl_rcv(6+iii)
          enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL WORKSHARE
       cnt(nys:nye,nze) = cnt(nys:nye,nze)+cnt_tmp_k(nys:nye)
!$OMP END PARALLEL WORKSHARE

       deallocate(ptcl_snd)
       deallocate(ptcl_rcv)


       !!TRANSFER TO KUP and RECEIVE FROM KDOWN
       call MPI_SENDRECV(cnt(nys,nze+1),nye-nys+1,mnpi,kup  ,400, &
                         cnt_tmp_k(nys),nye-nys+1,mnpi,kdown,400, &
                         ncomw,nstat,nerr)

       ssize = 6*sum(cnt(nys:nye,nze+1))
       rsize = 6*sum(cnt_tmp_k(nys:nye))
       allocate(ptcl_snd(ssize))
       allocate(ptcl_rcv(rsize))

       k=nze+1
!$OMP PARALLEL DO PRIVATE(iii,ii,j)
       do j=nys,nye
          do ii=1,cnt(j,k)
             iii = 6*(ii+sum(cnt(nys:j-1,k))-1)

             ptcl_snd(1+iii) = bff_ptcl(1+6*(ii-1),j,k)
             ptcl_snd(2+iii) = bff_ptcl(2+6*(ii-1),j,k)
             ptcl_snd(3+iii) = bff_ptcl(3+6*(ii-1),j,k)
             ptcl_snd(4+iii) = bff_ptcl(4+6*(ii-1),j,k)
             ptcl_snd(5+iii) = bff_ptcl(5+6*(ii-1),j,k)
             ptcl_snd(6+iii) = bff_ptcl(6+6*(ii-1),j,k)
          enddo
       enddo
!$OMP END PARALLEL DO
       call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,kup  ,401, &
                         ptcl_rcv(1),rsize,mnpr,kdown,401, &
                         ncomw,nstat,nerr)

       k=nzs
!$OMP PARALLEL DO PRIVATE(iii,ii,j)
       do j=nys,nye
          do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_k(j)
             iii = 6*(ii-1-cnt(j,k)+sum(cnt_tmp_k(nys:j-1)))

             bff_ptcl(1+6*(ii-1),j,k) = ptcl_rcv(1+iii)
             bff_ptcl(2+6*(ii-1),j,k) = ptcl_rcv(2+iii)
             bff_ptcl(3+6*(ii-1),j,k) = ptcl_rcv(3+iii)
             bff_ptcl(4+6*(ii-1),j,k) = ptcl_rcv(4+iii)
             bff_ptcl(5+6*(ii-1),j,k) = ptcl_rcv(5+iii)
             bff_ptcl(6+6*(ii-1),j,k) = ptcl_rcv(6+iii)
          enddo
       enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL WORKSHARE
       cnt(nys:nye,nzs) = cnt(nys:nye,nzs)+cnt_tmp_k(nys:nye)
!$OMP END PARALLEL WORKSHARE

       deallocate(ptcl_snd)
       deallocate(ptcl_rcv)


!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,j,k,cnt_tmp)
       do k=nzs,nze
       do j=nys,nye
          iii=0
          cnt_tmp = cnt2(j,k)
          loop1 :do ii=1,cnt2(j,k)
             if(cnt(j,k) == 0)then
                if(np2(j,k,isp) < flag(ii,j,k)) exit loop1
                do while(np2(j,k,isp) == flag(cnt_tmp,j,k))
                   np2(j,k,isp) = np2(j,k,isp)-1
                   if(np2(j,k,isp) < flag(ii,j,k)) exit loop1
                   cnt_tmp = cnt_tmp-1
                enddo
                up(1:6,flag(ii,j,k),j,k,isp) = up(1:6,np2(j,k,isp),j,k,isp)
                np2(j,k,isp) = np2(j,k,isp)-1
             else
                up(1,flag(ii,j,k),j,k,isp) = bff_ptcl(1+6*iii,j,k)
                up(2,flag(ii,j,k),j,k,isp) = bff_ptcl(2+6*iii,j,k)
                up(3,flag(ii,j,k),j,k,isp) = bff_ptcl(3+6*iii,j,k)
                up(4,flag(ii,j,k),j,k,isp) = bff_ptcl(4+6*iii,j,k)
                up(5,flag(ii,j,k),j,k,isp) = bff_ptcl(5+6*iii,j,k)
                up(6,flag(ii,j,k),j,k,isp) = bff_ptcl(6+6*iii,j,k)
                iii = iii+1
                cnt(j,k) = cnt(j,k)-1
             endif
          enddo loop1
          
          if(cnt(j,k) > 0)then
             do ii=1,cnt(j,k)
                up(1,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(6*iii+1+6*(ii-1),j,k)
                up(2,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(6*iii+2+6*(ii-1),j,k)
                up(3,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(6*iii+3+6*(ii-1),j,k)
                up(4,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(6*iii+4+6*(ii-1),j,k)
                up(5,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(6*iii+5+6*(ii-1),j,k)
                up(6,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(6*iii+6+6*(ii-1),j,k)
             enddo
          endif
       enddo
       enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(j,k)
       do k=nzs,nze
       do j=nys,nye
          np2(j,k,isp) = np2(j,k,isp)+cnt(j,k)
          if(np2(j,k,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(j,k,isp),j,k,isp
             stop
          endif
       enddo
       enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    enddo

!$OMP PARALLEL DO PRIVATE(j,k)
!$  do k=nzs-1,nze+1
!$  do j=nys-1,nye+1
!$     call omp_destroy_lock(lck(j,k))
!$  enddo
!$  enddo
!$OMP END PARALLEL DO

  end subroutine boundary__particle_yz


  subroutine boundary__field(uf,                                &
                             nxgs,nxge,nxs,nxe,nys,nye,nzs,nze, &
                             jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye, nzs, nze
    integer, intent(in)    :: jup, jdown, kup, kdown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    integer                :: i, j, k, ii
    real(8)                :: bff_snd_j(12*(nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_rcv_j(12*(nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_snd_k(12*(nxe-nxs+1)*(nye-nys+5))
    real(8)                :: bff_rcv_k(12*(nxe-nxs+1)*(nye-nys+5))

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(k-nzs))

       bff_snd_j(ii+1)  = uf(1,i,nys,k)
       bff_snd_j(ii+2)  = uf(2,i,nys,k)
       bff_snd_j(ii+3)  = uf(3,i,nys,k)
       bff_snd_j(ii+4)  = uf(4,i,nys,k)
       bff_snd_j(ii+5)  = uf(5,i,nys,k)
       bff_snd_j(ii+6)  = uf(6,i,nys,k)
       bff_snd_j(ii+7)  = uf(1,i,nys+1,k)
       bff_snd_j(ii+8)  = uf(2,i,nys+1,k)
       bff_snd_j(ii+9)  = uf(3,i,nys+1,k)
       bff_snd_j(ii+10) = uf(4,i,nys+1,k)
       bff_snd_j(ii+11) = uf(5,i,nys+1,k)
       bff_snd_j(ii+12) = uf(6,i,nys+1,k)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_j(1),12*(nxe-nxs+1)*(nze-nzs+1),mnpr,jdown,110, &
                      bff_rcv_j(1),12*(nxe-nxs+1)*(nze-nzs+1),mnpr,jup  ,110, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(k-nzs))
       uf(1,i,nye+1,k) = bff_rcv_j(ii+1)   
       uf(2,i,nye+1,k) = bff_rcv_j(ii+2)
       uf(3,i,nye+1,k) = bff_rcv_j(ii+3)
       uf(4,i,nye+1,k) = bff_rcv_j(ii+4)   
       uf(5,i,nye+1,k) = bff_rcv_j(ii+5)
       uf(6,i,nye+1,k) = bff_rcv_j(ii+6)
       uf(1,i,nye+2,k) = bff_rcv_j(ii+7)   
       uf(2,i,nye+2,k) = bff_rcv_j(ii+8)
       uf(3,i,nye+2,k) = bff_rcv_j(ii+9)
       uf(4,i,nye+2,k) = bff_rcv_j(ii+10)   
       uf(5,i,nye+2,k) = bff_rcv_j(ii+11)
       uf(6,i,nye+2,k) = bff_rcv_j(ii+12)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(k-nzs))

       bff_snd_j(ii+1)  = uf(1,i,nye-1,k)
       bff_snd_j(ii+2)  = uf(2,i,nye-1,k)
       bff_snd_j(ii+3)  = uf(3,i,nye-1,k)
       bff_snd_j(ii+4)  = uf(4,i,nye-1,k)
       bff_snd_j(ii+5)  = uf(5,i,nye-1,k)
       bff_snd_j(ii+6)  = uf(6,i,nye-1,k)
       bff_snd_j(ii+7)  = uf(1,i,nye,k)
       bff_snd_j(ii+8)  = uf(2,i,nye,k)
       bff_snd_j(ii+9)  = uf(3,i,nye,k)
       bff_snd_j(ii+10) = uf(4,i,nye,k)
       bff_snd_j(ii+11) = uf(5,i,nye,k)
       bff_snd_j(ii+12) = uf(6,i,nye,k)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_j(1),12*(nxe-nxs+1)*(nze-nzs+1),mnpr,jup  ,100, &
                      bff_rcv_j(1),12*(nxe-nxs+1)*(nze-nzs+1),mnpr,jdown,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(k-nzs))

       uf(1,i,nys-2,k) = bff_rcv_j(ii+1)
       uf(2,i,nys-2,k) = bff_rcv_j(ii+2)
       uf(3,i,nys-2,k) = bff_rcv_j(ii+3)
       uf(4,i,nys-2,k) = bff_rcv_j(ii+4)   
       uf(5,i,nys-2,k) = bff_rcv_j(ii+5)
       uf(6,i,nys-2,k) = bff_rcv_j(ii+6)
       uf(1,i,nys-1,k) = bff_rcv_j(ii+7)   
       uf(2,i,nys-1,k) = bff_rcv_j(ii+8)
       uf(3,i,nys-1,k) = bff_rcv_j(ii+9)
       uf(4,i,nys-1,k) = bff_rcv_j(ii+10)   
       uf(5,i,nys-1,k) = bff_rcv_j(ii+11)
       uf(6,i,nys-1,k) = bff_rcv_j(ii+12)
    enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-2,nye+2
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(j-(nys-2)))

       bff_snd_k(ii+1)  = uf(1,i,j,nzs)
       bff_snd_k(ii+2)  = uf(2,i,j,nzs)
       bff_snd_k(ii+3)  = uf(3,i,j,nzs)
       bff_snd_k(ii+4)  = uf(4,i,j,nzs)
       bff_snd_k(ii+5)  = uf(5,i,j,nzs)
       bff_snd_k(ii+6)  = uf(6,i,j,nzs)
       bff_snd_k(ii+7)  = uf(1,i,j,nzs+1)
       bff_snd_k(ii+8)  = uf(2,i,j,nzs+1)
       bff_snd_k(ii+9)  = uf(3,i,j,nzs+1)
       bff_snd_k(ii+10) = uf(4,i,j,nzs+1)
       bff_snd_k(ii+11) = uf(5,i,j,nzs+1)
       bff_snd_k(ii+12) = uf(6,i,j,nzs+1)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_k(1),12*(nxe-nxs+1)*(nye-nys+5),mnpr,kdown,210, &
                      bff_rcv_k(1),12*(nxe-nxs+1)*(nye-nys+5),mnpr,kup  ,210, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-2,nye+2
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(j-(nys-2)))

       uf(1,i,j,nze+1) = bff_rcv_k(ii+1)   
       uf(2,i,j,nze+1) = bff_rcv_k(ii+2)
       uf(3,i,j,nze+1) = bff_rcv_k(ii+3)
       uf(4,i,j,nze+1) = bff_rcv_k(ii+4)   
       uf(5,i,j,nze+1) = bff_rcv_k(ii+5)
       uf(6,i,j,nze+1) = bff_rcv_k(ii+6)
       uf(1,i,j,nze+2) = bff_rcv_k(ii+7)   
       uf(2,i,j,nze+2) = bff_rcv_k(ii+8)
       uf(3,i,j,nze+2) = bff_rcv_k(ii+9)
       uf(4,i,j,nze+2) = bff_rcv_k(ii+10)   
       uf(5,i,j,nze+2) = bff_rcv_k(ii+11)
       uf(6,i,j,nze+2) = bff_rcv_k(ii+12)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-2,nye+2
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(j-(nys-2)))

       bff_snd_k(ii+1)  = uf(1,i,j,nze-1)
       bff_snd_k(ii+2)  = uf(2,i,j,nze-1)
       bff_snd_k(ii+3)  = uf(3,i,j,nze-1)
       bff_snd_k(ii+4)  = uf(4,i,j,nze-1)
       bff_snd_k(ii+5)  = uf(5,i,j,nze-1)
       bff_snd_k(ii+6)  = uf(6,i,j,nze-1)
       bff_snd_k(ii+7)  = uf(1,i,j,nze)
       bff_snd_k(ii+8)  = uf(2,i,j,nze)
       bff_snd_k(ii+9)  = uf(3,i,j,nze)
       bff_snd_k(ii+10) = uf(4,i,j,nze)
       bff_snd_k(ii+11) = uf(5,i,j,nze)
       bff_snd_k(ii+12) = uf(6,i,j,nze)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_k(1),12*(nxe-nxs+1)*(nye-nys+5),mnpr,kup  ,200, &
                      bff_rcv_k(1),12*(nxe-nxs+1)*(nye-nys+5),mnpr,kdown,200, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-2,nye+2
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(j-(nys-2)))

       uf(1,i,j,nzs-2) = bff_rcv_k(ii+1)
       uf(2,i,j,nzs-2) = bff_rcv_k(ii+2)
       uf(3,i,j,nzs-2) = bff_rcv_k(ii+3)
       uf(4,i,j,nzs-2) = bff_rcv_k(ii+4)   
       uf(5,i,j,nzs-2) = bff_rcv_k(ii+5)
       uf(6,i,j,nzs-2) = bff_rcv_k(ii+6)
       uf(1,i,j,nzs-1) = bff_rcv_k(ii+7)   
       uf(2,i,j,nzs-1) = bff_rcv_k(ii+8)
       uf(3,i,j,nzs-1) = bff_rcv_k(ii+9)
       uf(4,i,j,nzs-1) = bff_rcv_k(ii+10)   
       uf(5,i,j,nzs-1) = bff_rcv_k(ii+11)
       uf(6,i,j,nzs-1) = bff_rcv_k(ii+12)
    enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
       uf(1,nxs-2,j,k) = +uf(1,nxs+1,j,k)
       uf(2,nxs-2,j,k) = +2.*uf(2,nxs,j,k)-uf(2,nxs+2,j,k)
       uf(3,nxs-2,j,k) = +2.*uf(3,nxs,j,k)-uf(3,nxs+2,j,k)
!       uf(4,nxs-2,j,k) = +uf(4,nxs+2,j,k)
       uf(5,nxs-2,j,k) = +uf(5,nxs+1,j,k)
       uf(6,nxs-2,j,k) = +uf(6,nxs+1,j,k)

       uf(1,nxs-1,j,k) = +uf(1,nxs  ,j,k)
       uf(2,nxs-1,j,k) = +2.*uf(2,nxs,j,k)-uf(2,nxs+1,j,k)
       uf(3,nxs-1,j,k) = +2.*uf(3,nxs,j,k)-uf(3,nxs+1,j,k)
!       uf(4,nxs-1,j,k) = +uf(4,nxs+1,j,k)
       uf(5,nxs-1,j,k) = +uf(5,nxs  ,j,k)
       uf(6,nxs-1,j,k) = +uf(6,nxs  ,j,k)

       uf(1,nxe  ,j,k) = +uf(1,nxe-1,j,k)
       uf(2,nxe+1,j,k) = +2.*uf(2,nxe,j,k)-uf(2,nxe-1,j,k)
       uf(3,nxe+1,j,k) = +2.*uf(3,nxe,j,k)-uf(3,nxe-1,j,k)
!       uf(4,nxe+1,j,k) = +uf(4,nxe-1,j,k)
       uf(5,nxe  ,j,k) = +uf(5,nxe-1,j,k)
       uf(6,nxe  ,j,k) = +uf(6,nxe-1,j,k)

       uf(1,nxe+1,j,k) = +uf(1,nxe-2,j,k)
       uf(2,nxe+2,j,k) = +2.*uf(2,nxe,j,k)-uf(2,nxe-2,j,k)
       uf(3,nxe+2,j,k) = +2.*uf(3,nxe,j,k)-uf(3,nxe-2,j,k)
!       uf(4,nxe+2,j,k) = +uf(4,nxe-2,j,k)
       uf(5,nxe+1,j,k) = +uf(5,nxe-2,j,k)
       uf(6,nxe+1,j,k) = +uf(6,nxe-2,j,k)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary__field


  subroutine boundary__curre(uj,nxs,nxe,nys,nye,nzs,nze, &
                             jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze
    integer, intent(in)    :: jup, jdown, kup, kdown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uj(3,nxs-2:nxe+2,nys-2:nye+2,nzs-2:nze+2)
    integer                :: i, j, k, ii
    real(8)                :: bff_snd_j(6*(nxe-nxs+5)*(nze-nzs+5))
    real(8)                :: bff_rcv_j(6*(nxe-nxs+5)*(nze-nzs+5))
    real(8)                :: bff_snd_k(6*(nxe-nxs+5)*(nye-nys+1))
    real(8)                :: bff_rcv_k(6*(nxe-nxs+5)*(nye-nys+1))

    
!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs-2,nze+2
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(k-(nzs-2)))

       bff_snd_j(ii+1) = uj(1,i,nys-2,k)
       bff_snd_j(ii+2) = uj(2,i,nys-2,k)
       bff_snd_j(ii+3) = uj(3,i,nys-2,k)
       bff_snd_j(ii+4) = uj(1,i,nys-1,k)
       bff_snd_j(ii+5) = uj(2,i,nys-1,k)
       bff_snd_j(ii+6) = uj(3,i,nys-1,k)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_j(1),6*(nxe-nxs+5)*(nze-nzs+5),mnpr,jdown,110, &
                      bff_rcv_j(1),6*(nxe-nxs+5)*(nze-nzs+5),mnpr,jup  ,110, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs-2,nze+2
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(k-(nzs-2)))

       uj(1,i,nye-1,k) = uj(1,i,nye-1,k)+bff_rcv_j(ii+1)
       uj(2,i,nye-1,k) = uj(2,i,nye-1,k)+bff_rcv_j(ii+2)
       uj(3,i,nye-1,k) = uj(3,i,nye-1,k)+bff_rcv_j(ii+3)
       uj(1,i,nye  ,k) = uj(1,i,nye  ,k)+bff_rcv_j(ii+4)
       uj(2,i,nye  ,k) = uj(2,i,nye  ,k)+bff_rcv_j(ii+5)
       uj(3,i,nye  ,k) = uj(3,i,nye  ,k)+bff_rcv_j(ii+6)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs-2,nze+2
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(k-(nzs-2)))

       bff_snd_j(ii+1) = uj(1,i,nye+1,k)
       bff_snd_j(ii+2) = uj(2,i,nye+1,k)
       bff_snd_j(ii+3) = uj(3,i,nye+1,k)
       bff_snd_j(ii+4) = uj(1,i,nye+2,k)
       bff_snd_j(ii+5) = uj(2,i,nye+2,k)
       bff_snd_j(ii+6) = uj(3,i,nye+2,k)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_j(1),6*(nxe-nxs+5)*(nze-nzs+5),mnpr,jup  ,100, &
                      bff_rcv_j(1),6*(nxe-nxs+5)*(nze-nzs+5),mnpr,jdown,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs-2,nze+2
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(k-(nzs-2)))

       uj(1,i,nys  ,k) = uj(1,i,nys  ,k)+bff_rcv_j(ii+1)
       uj(2,i,nys  ,k) = uj(2,i,nys  ,k)+bff_rcv_j(ii+2)
       uj(3,i,nys  ,k) = uj(3,i,nys  ,k)+bff_rcv_j(ii+3)
       uj(1,i,nys+1,k) = uj(1,i,nys+1,k)+bff_rcv_j(ii+4)
       uj(2,i,nys+1,k) = uj(2,i,nys+1,k)+bff_rcv_j(ii+5)
       uj(3,i,nys+1,k) = uj(3,i,nys+1,k)+bff_rcv_j(ii+6)
    enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys,nye
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(j-nys))

       bff_snd_k(ii+1) = uj(1,i,j,nzs-2)
       bff_snd_k(ii+2) = uj(2,i,j,nzs-2)
       bff_snd_k(ii+3) = uj(3,i,j,nzs-2)
       bff_snd_k(ii+4) = uj(1,i,j,nzs-1)
       bff_snd_k(ii+5) = uj(2,i,j,nzs-1)
       bff_snd_k(ii+6) = uj(3,i,j,nzs-1)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_k(1),6*(nxe-nxs+5)*(nye-nys+1),mnpr,kdown,210, &
                      bff_rcv_k(1),6*(nxe-nxs+5)*(nye-nys+1),mnpr,kup  ,210, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,j,ii)
    do j=nys,nye
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(j-nys))

       uj(1,i,j,nze-1) = uj(1,i,j,nze-1)+bff_rcv_k(ii+1)
       uj(2,i,j,nze-1) = uj(2,i,j,nze-1)+bff_rcv_k(ii+2)
       uj(3,i,j,nze-1) = uj(3,i,j,nze-1)+bff_rcv_k(ii+3)
       uj(1,i,j,nze  ) = uj(1,i,j,nze  )+bff_rcv_k(ii+4)
       uj(2,i,j,nze  ) = uj(2,i,j,nze  )+bff_rcv_k(ii+5)
       uj(3,i,j,nze  ) = uj(3,i,j,nze  )+bff_rcv_k(ii+6)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j,ii)
    do j=nys,nye
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(j-nys))

       bff_snd_k(ii+1) = uj(1,i,j,nze+1)
       bff_snd_k(ii+2) = uj(2,i,j,nze+1)
       bff_snd_k(ii+3) = uj(3,i,j,nze+1)
       bff_snd_k(ii+4) = uj(1,i,j,nze+2)
       bff_snd_k(ii+5) = uj(2,i,j,nze+2)
       bff_snd_k(ii+6) = uj(3,i,j,nze+2)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_k(1),6*(nxe-nxs+5)*(nye-nys+1),mnpr,kup  ,200, &
                      bff_rcv_k(1),6*(nxe-nxs+5)*(nye-nys+1),mnpr,kdown,200, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys,nye
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2)+(nxe-nxs+5)*(j-nys))

       uj(1,i,j,nzs  ) = uj(1,i,j,nzs  )+bff_rcv_k(ii+1)
       uj(2,i,j,nzs  ) = uj(2,i,j,nzs  )+bff_rcv_k(ii+2)
       uj(3,i,j,nzs  ) = uj(3,i,j,nzs  )+bff_rcv_k(ii+3)
       uj(1,i,j,nzs+1) = uj(1,i,j,nzs+1)+bff_rcv_k(ii+4)
       uj(2,i,j,nzs+1) = uj(2,i,j,nzs+1)+bff_rcv_k(ii+5)
       uj(3,i,j,nzs+1) = uj(3,i,j,nzs+1)+bff_rcv_k(ii+6)
    enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs,nze
    do j=nys,nye
!       uj(1,nxs+1,j,k) = uj(1,nxs+1,j,k)-uj(1,nxs-1,j,k)
!       uj(1,nxs+2,j,k) = uj(1,nxs+2,j,k)-uj(1,nxs-2,j,k)
       uj(1,nxs  ,j,k) = 2.*uj(1,nxs  ,j,k)
       uj(2,nxs  ,j,k) = uj(2,nxs  ,j,k)+uj(2,nxs-1,j,k)
       uj(3,nxs  ,j,k) = uj(3,nxs  ,j,k)+uj(3,nxs-1,j,k)
       uj(2,nxs+1,j,k) = uj(2,nxs+1,j,k)+uj(2,nxs-2,j,k)
       uj(3,nxs+1,j,k) = uj(3,nxs+1,j,k)+uj(3,nxs-2,j,k)

!       uj(1,nxe-2,j,k) = uj(1,nxe-2,j,k)-uj(1,nxe+2,j,k)
!       uj(1,nxe-1,j,k) = uj(1,nxe-1,j,k)-uj(1,nxe+1,j,k)
       uj(2,nxe-2,j,k) = uj(2,nxe-2,j,k)+uj(2,nxe+1,j,k)
       uj(3,nxe-2,j,k) = uj(3,nxe-2,j,k)+uj(3,nxe+1,j,k)
       uj(2,nxe-1,j,k) = uj(2,nxe-1,j,k)+uj(2,nxe  ,j,k)
       uj(3,nxe-1,j,k) = uj(3,nxe-1,j,k)+uj(3,nxe  ,j,k)
       uj(1,nxe  ,j,k) = 2.*uj(1,nxe  ,j,k)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary__curre


  subroutine boundary__phi(phi,                     &
                           nxs,nxe,nys,nye,nzs,nze, &
                           jup,jdown,kup,kdown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze
    integer, intent(in)    :: jup, jdown, kup, kdown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: phi(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    integer                :: i, j, k, ii
    real(8)                :: bff_snd_j((nxe-nxs+1)*(nze-nzs+1)), bff_rcv_j((nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_snd_k((nxe-nxs+1)*(nye-nys+3)), bff_rcv_k((nxe-nxs+1)*(nye-nys+3))

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(k-nzs)
       bff_snd_j(ii) = phi(i,nys,k)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_j(1),(nxe-nxs+1)*(nze-nzs+1),mnpr,jdown,101, &
                      bff_rcv_j(1),(nxe-nxs+1)*(nze-nzs+1),mnpr,jup  ,101, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(k-nzs)
       phi(i,nye+1,k) = bff_rcv_j(ii)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(k-nzs)
       bff_snd_j(ii) = phi(i,nye,k)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_j(1),(nxe-nxs+1)*(nze-nzs+1),mnpr,jup  ,100, &
                      bff_rcv_j(1),(nxe-nxs+1)*(nze-nzs+1),mnpr,jdown,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(k-nzs)
       phi(i,nys-1,k) = bff_rcv_j(ii)
    enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(j-(nys-1))
       bff_snd_k(ii) = phi(i,j,nzs)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_k(1),(nxe-nxs+1)*(nye-nys+3),mnpr,kdown,201, &
                      bff_rcv_k(1),(nxe-nxs+1)*(nye-nys+3),mnpr,kup  ,201, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(j-(nys-1))
       phi(i,j,nze+1) = bff_rcv_k(ii)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(j-(nys-1))
       bff_snd_k(ii) = phi(i,j,nze)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_k(1),(nxe-nxs+1)*(nye-nys+3),mnpr,kup  ,200, &
                      bff_rcv_k(1),(nxe-nxs+1)*(nye-nys+3),mnpr,kdown,200, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs,nxe
       ii = i-nxs+1+(nxe-nxs+1)*(j-(nys-1))
       phi(i,j,nzs-1) = bff_rcv_k(ii)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary__phi


  subroutine boundary__mom(den,vel,temp,nxgs,nxge,nys,nye,nzs,nze,nsp)

    integer, intent(in)    :: nxgs, nxge, nys, nye, nzs, nze, nsp
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,nsp)
    real(8), intent(inout) :: vel(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,3,nsp)
    real(8), intent(inout) :: temp(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,3,nsp)

    !reflective condition

!$OMP PARALLEL WORKSHARE
    den(nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:nsp) = den(nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                               +den(nxgs-1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    den(nxge-1,nys-1:nye+1,nzs-1:nze+1,1:nsp) = den(nxge-1,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                               +den(nxge  ,nys-1:nye+1,nzs-1:nze+1,1:nsp)

    vel(nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = vel(nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) &
                                                   +vel(nxgs-1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp)
    vel(nxge-1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = vel(nxge-1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) &
                                                   +vel(nxge  ,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp)

    temp(nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = temp(nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) &
                                                    +temp(nxgs-1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp)
    temp(nxge-1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) = temp(nxge-1,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp) &
                                                    +temp(nxge  ,nys-1:nye+1,nzs-1:nze+1,1:3,1:nsp)
!$OMP END PARALLEL WORKSHARE

  end subroutine boundary__mom


end module boundary
