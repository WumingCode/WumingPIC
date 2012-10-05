module boundary

  implicit none

  private

  public  :: boundary__field
  public  :: boundary__particle
  public  :: boundary__curre


contains


  subroutine boundary__particle(up,                                   &
                                np,nsp,np2,nygs,nyge,nxs,nxe,nys,nye, &
                                nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

!$  use omp_lib

    integer, intent(in)        :: np, nsp, nygs, nyge, nxs, nxe, nys, nye
    integer, intent(in)        :: nup, ndown, mnpi, mnpr, ncomw
    integer, intent(inout)     :: nerr, nstat(:)
    integer, intent(inout)     :: np2(nys:nye,nsp)
    real(8), intent(inout)     :: up(5,np,nys:nye,nsp)
    logical, save              :: lflag=.true.
!$  integer(omp_lock_kind)     :: lck(nys-1:nye+1)
    integer                    :: j, ii, iii, isp, ipos, jpos, ieq
    integer                    :: cnt(nys-1:nye+1), cnt2(nys:nye), cnt_tmp
    integer, save, allocatable :: flag(:,:)
    real(8), save, allocatable :: bff_ptcl(:,:)

    if(lflag)then
       allocate(flag(np,nys:nye))
       allocate(bff_ptcl(np,nys-1:nye+1))
       lflag=.false.
    endif

!$OMP PARALLEL DO PRIVATE(j)
!$    do j=nys-1,nye+1
!$       call omp_init_lock(lck(j))
!$    enddo
!$OMP END PARALLEL DO

    do isp=1,nsp

!$OMP PARALLEL

!$OMP WORKSHARE
       cnt(nys-1:nye+1) = 0
       cnt2(nys:nye) = 0
!$OMP END WORKSHARE

!$OMP DO PRIVATE(ii,j,ipos,jpos)
       do j=nys,nye
          do ii=1,np2(j,isp)

             ipos = int(up(1,ii,j,isp))
             jpos = int(up(2,ii,j,isp))

             if(ipos <= nxs-1)then
                up(1,ii,j,isp) = 2.0*nxs-up(1,ii,j,isp)
                up(3,ii,j,isp) = -up(3,ii,j,isp)
             else if(ipos >= nxe)then
                up(1,ii,j,isp) = 2.0*nxe-up(1,ii,j,isp)
                up(3,ii,j,isp) = -up(3,ii,j,isp)
             endif

             if(jpos /= j)then
                if(jpos <= nygs-1)then
                   up(2,ii,j,isp) = up(2,ii,j,isp)+(nyge-nygs+1)
                else if(jpos >= nyge+1)then
                   up(2,ii,j,isp) = up(2,ii,j,isp)-(nyge-nygs+1)
                endif

!$              call omp_set_lock(lck(jpos))
                bff_ptcl(1+5*cnt(jpos),jpos) = up(1,ii,j,isp)
                bff_ptcl(2+5*cnt(jpos),jpos) = up(2,ii,j,isp)
                bff_ptcl(3+5*cnt(jpos),jpos) = up(3,ii,j,isp)
                bff_ptcl(4+5*cnt(jpos),jpos) = up(4,ii,j,isp)
                bff_ptcl(5+5*cnt(jpos),jpos) = up(5,ii,j,isp)
                cnt(jpos) = cnt(jpos)+1
!$              call omp_unset_lock(lck(jpos))

                cnt2(j) = cnt2(j)+1
                flag(cnt2(j),j) = ii
             endif
          enddo
       enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL 

       !transfer to rank-1
       call MPI_SENDRECV(cnt(nys-1),1,mnpi,ndown,100, &
                         cnt_tmp   ,1,mnpi,nup  ,100, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1           ,nys-1),5*cnt(nys-1),mnpr,ndown,101, &
                         bff_ptcl(5*cnt(nye)+1,nye  ),5*cnt_tmp   ,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)
       cnt(nye) = cnt(nye)+cnt_tmp

       !transfer to rank+1
       call MPI_SENDRECV(cnt(nye+1),1,mnpi,nup  ,200, &
                         cnt_tmp   ,1,mnpi,ndown,200, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1           ,nye+1),5*cnt(nye+1),mnpr,nup  ,201, &
                         bff_ptcl(5*cnt(nys)+1,nys  ),5*cnt_tmp   ,mnpr,ndown,201, &
                         ncomw,nstat,nerr)
       cnt(nys) = cnt(nys)+cnt_tmp

!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,j,ieq,cnt_tmp)
       do j=nys,nye
          iii=0
          cnt_tmp = cnt2(j)
          loop1 :do ii=1,cnt2(j)
             if(cnt(j) == 0)then
                if(np2(j,isp) < flag(ii,j)) exit loop1
                do while(np2(j,isp) == flag(cnt_tmp,j))
                   np2(j,isp) = np2(j,isp)-1
                   if(np2(j,isp) < flag(ii,j)) exit loop1
                   cnt_tmp = cnt_tmp-1
                enddo
                do ieq=1,5
                   up(ieq,flag(ii,j),j,isp) = up(ieq,np2(j,isp),j,isp)
                enddo
                np2(j,isp) = np2(j,isp)-1
             else
                up(1,flag(ii,j),j,isp) = bff_ptcl(1+5*iii,j)
                up(2,flag(ii,j),j,isp) = bff_ptcl(2+5*iii,j)
                up(3,flag(ii,j),j,isp) = bff_ptcl(3+5*iii,j)
                up(4,flag(ii,j),j,isp) = bff_ptcl(4+5*iii,j)
                up(5,flag(ii,j),j,isp) = bff_ptcl(5+5*iii,j)
                iii = iii+1
                cnt(j) = cnt(j)-1
             endif
          enddo loop1
          
          if(cnt(j) > 0)then
             do ii=1,cnt(j)
                up(1,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+1+5*(ii-1),j)
                up(2,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+2+5*(ii-1),j)
                up(3,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+3+5*(ii-1),j)
                up(4,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+4+5*(ii-1),j)
                up(5,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+5+5*(ii-1),j)
             enddo
          endif
       enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(j)
       do j=nys,nye
          np2(j,isp) = np2(j,isp)+cnt(j)
          if(np2(j,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(j,isp),j,isp
             stop
          endif
       enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    enddo

!$OMP PARALLEL DO PRIVATE(j)
!$  do j=nys-1,nye+1
!$     call omp_destroy_lock(lck(j))
!$  enddo
!$OMP END PARALLEL DO

  end subroutine boundary__particle


  subroutine boundary__field(uf,                        &
                             nxgs,nxge,nxs,nxe,nys,nye, &
                             nup,ndown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxgs, nxge, nxs, nxe, nys, nye
    integer, intent(in)    :: nup, ndown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer                :: i, j, ii
    real(8)                :: bff_snd(12*(nxe-nxs+1)), bff_rcv(12*(nxe-nxs+1))

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i=nxs,nxe
       ii = 12*(i-nxs)
       bff_snd(ii+1)  = uf(1,i,nys)
       bff_snd(ii+2)  = uf(2,i,nys)
       bff_snd(ii+3)  = uf(3,i,nys)
       bff_snd(ii+4)  = uf(4,i,nys)
       bff_snd(ii+5)  = uf(5,i,nys)
       bff_snd(ii+6)  = uf(6,i,nys)
       bff_snd(ii+7)  = uf(1,i,nys+1)
       bff_snd(ii+8)  = uf(2,i,nys+1)
       bff_snd(ii+9)  = uf(3,i,nys+1)
       bff_snd(ii+10) = uf(4,i,nys+1)
       bff_snd(ii+11) = uf(5,i,nys+1)
       bff_snd(ii+12) = uf(6,i,nys+1)
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1),12*(nxe-nxs+1),mnpr,ndown,110, &
                      bff_rcv(1),12*(nxe-nxs+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i=nxs,nxe
       ii = 12*(i-nxs)
       uf(1,i,nye+1) = bff_rcv(ii+1)   
       uf(2,i,nye+1) = bff_rcv(ii+2)
       uf(3,i,nye+1) = bff_rcv(ii+3)
       uf(4,i,nye+1) = bff_rcv(ii+4)   
       uf(5,i,nye+1) = bff_rcv(ii+5)
       uf(6,i,nye+1) = bff_rcv(ii+6)
       uf(1,i,nye+2) = bff_rcv(ii+7)   
       uf(2,i,nye+2) = bff_rcv(ii+8)
       uf(3,i,nye+2) = bff_rcv(ii+9)
       uf(4,i,nye+2) = bff_rcv(ii+10)   
       uf(5,i,nye+2) = bff_rcv(ii+11)
       uf(6,i,nye+2) = bff_rcv(ii+12)
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii)
    do i=nxs,nxe
       ii = 12*(i-nxs)
       bff_snd(ii+1)  = uf(1,i,nye-1)
       bff_snd(ii+2)  = uf(2,i,nye-1)
       bff_snd(ii+3)  = uf(3,i,nye-1)
       bff_snd(ii+4)  = uf(4,i,nye-1)
       bff_snd(ii+5)  = uf(5,i,nye-1)
       bff_snd(ii+6)  = uf(6,i,nye-1)
       bff_snd(ii+7)  = uf(1,i,nye)
       bff_snd(ii+8)  = uf(2,i,nye)
       bff_snd(ii+9)  = uf(3,i,nye)
       bff_snd(ii+10) = uf(4,i,nye)
       bff_snd(ii+11) = uf(5,i,nye)
       bff_snd(ii+12) = uf(6,i,nye)
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1),12*(nxe-nxs+1),mnpr,nup  ,100, &
                      bff_rcv(1),12*(nxe-nxs+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i=nxs,nxe
       ii = 12*(i-nxs)
       uf(1,i,nys-2) = bff_rcv(ii+1)   
       uf(2,i,nys-2) = bff_rcv(ii+2)
       uf(3,i,nys-2) = bff_rcv(ii+3)
       uf(4,i,nys-2) = bff_rcv(ii+4)   
       uf(5,i,nys-2) = bff_rcv(ii+5)
       uf(6,i,nys-2) = bff_rcv(ii+6)
       uf(1,i,nys-1) = bff_rcv(ii+7)   
       uf(2,i,nys-1) = bff_rcv(ii+8)
       uf(3,i,nys-1) = bff_rcv(ii+9)
       uf(4,i,nys-1) = bff_rcv(ii+10)   
       uf(5,i,nys-1) = bff_rcv(ii+11)
       uf(6,i,nys-1) = bff_rcv(ii+12)
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(j)
    do j=nys-2,nye+2
       uf(1,nxs-2,j) = -uf(1,nxs+1,j)
       uf(2,nxs-2,j) = +uf(2,nxs+2,j)
       uf(3,nxs-2,j) = +uf(3,nxs+2,j)
       uf(4,nxs-2,j) = +uf(4,nxs+2,j)
       uf(5,nxs-2,j) = -uf(5,nxs+1,j)
       uf(6,nxs-2,j) = -uf(6,nxs+1,j)

       uf(1,nxs-1,j) = -uf(1,nxs  ,j)
       uf(2,nxs-1,j) = +uf(2,nxs+1,j)
       uf(3,nxs-1,j) = +uf(3,nxs+1,j)
       uf(4,nxs-1,j) = +uf(4,nxs+1,j)
       uf(5,nxs-1,j) = -uf(5,nxs  ,j)
       uf(6,nxs-1,j) = -uf(6,nxs  ,j)

       uf(1,nxe  ,j) = -uf(1,nxe-1,j)
       uf(2,nxe+1,j) = +uf(2,nxe-1,j)
       uf(3,nxe+1,j) = +uf(3,nxe-1,j)
       uf(4,nxe+1,j) = +uf(4,nxe-1,j)
       uf(5,nxe  ,j) = -uf(5,nxe-1,j)
       uf(6,nxe  ,j) = -uf(6,nxe-1,j)

       uf(1,nxe+1,j) = -uf(1,nxe-2,j)
       uf(2,nxe+2,j) = +uf(2,nxe-2,j)
       uf(3,nxe+2,j) = +uf(3,nxe-2,j)
       uf(4,nxe+2,j) = +uf(4,nxe-2,j)
       uf(5,nxe+1,j) = -uf(5,nxe-2,j)
       uf(6,nxe+1,j) = -uf(6,nxe-2,j)
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary__field


  subroutine boundary__curre(uj,nxs,nxe,nys,nye, &
                             nup,ndown,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nys, nye
    integer, intent(in)    :: nup, ndown, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    integer                :: i, j, ii
    real(8)                :: bff_rcv(6*(nxe-nxs+4+1)), bff_snd(6*(nxe-nxs+4+1))

    !send to rank-1
!$OMP PARALLEL DO PRIVATE(i,ii)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nys-2)
       bff_snd(ii+2) = uj(2,i,nys-2)
       bff_snd(ii+3) = uj(3,i,nys-2)
       bff_snd(ii+4) = uj(1,i,nys-1)
       bff_snd(ii+5) = uj(2,i,nys-1)
       bff_snd(ii+6) = uj(3,i,nys-1)
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,ndown,110, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nye-1) = uj(1,i,nye-1)+bff_rcv(ii+1)
       uj(2,i,nye-1) = uj(2,i,nye-1)+bff_rcv(ii+2)
       uj(3,i,nye-1) = uj(3,i,nye-1)+bff_rcv(ii+3)
       uj(1,i,nye  ) = uj(1,i,nye  )+bff_rcv(ii+4)
       uj(2,i,nye  ) = uj(2,i,nye  )+bff_rcv(ii+5)
       uj(3,i,nye  ) = uj(3,i,nye  )+bff_rcv(ii+6)
    enddo
!$OMP END DO NOWAIT

    !send to rank+1
!$OMP DO PRIVATE(i,ii)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nye+1)
       bff_snd(ii+2) = uj(2,i,nye+1)
       bff_snd(ii+3) = uj(3,i,nye+1)
       bff_snd(ii+4) = uj(1,i,nye+2)
       bff_snd(ii+5) = uj(2,i,nye+2)
       bff_snd(ii+6) = uj(3,i,nye+2)
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,nup  ,100, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nys  ) = uj(1,i,nys  )+bff_rcv(ii+1)
       uj(2,i,nys  ) = uj(2,i,nys  )+bff_rcv(ii+2)
       uj(3,i,nys  ) = uj(3,i,nys  )+bff_rcv(ii+3)
       uj(1,i,nys+1) = uj(1,i,nys+1)+bff_rcv(ii+4)
       uj(2,i,nys+1) = uj(2,i,nys+1)+bff_rcv(ii+5)
       uj(3,i,nys+1) = uj(3,i,nys+1)+bff_rcv(ii+6)
    enddo
!$OMP END PARALLEL DO

    !boundary condition in x
!$OMP PARALLEL DO PRIVATE(j)
    do j=nys,nye
       uj(1,nxs+1,j) = uj(1,nxs+1,j)-uj(1,nxs-1,j)
       uj(1,nxs+2,j) = uj(1,nxs+2,j)-uj(1,nxs-2,j)
       uj(2,nxs  ,j) = uj(2,nxs  ,j)+uj(2,nxs-1,j)
       uj(3,nxs  ,j) = uj(3,nxs  ,j)+uj(3,nxs-1,j)
       uj(2,nxs+1,j) = uj(2,nxs+1,j)+uj(2,nxs-2,j)
       uj(3,nxs+1,j) = uj(3,nxs+1,j)+uj(3,nxs-2,j)

       uj(1,nxe-2,j) = uj(1,nxe-2,j)-uj(1,nxe+2,j)
       uj(1,nxe-1,j) = uj(1,nxe-1,j)-uj(1,nxe+1,j)
       uj(2,nxe-2,j) = uj(2,nxe-2,j)+uj(2,nxe+1,j)
       uj(3,nxe-2,j) = uj(3,nxe-2,j)+uj(3,nxe+1,j)
       uj(2,nxe-1,j) = uj(2,nxe-1,j)+uj(2,nxe  ,j)
       uj(3,nxe-1,j) = uj(3,nxe-1,j)+uj(3,nxe  ,j)
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary__curre


end module boundary
