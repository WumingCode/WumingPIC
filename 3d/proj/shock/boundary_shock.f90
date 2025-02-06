module boundary_shock
  use mpi
  implicit none

  private

  public :: boundary_shock__init
  public :: boundary_shock__dfield
  public :: boundary_shock__particle_yz
  public :: boundary_shock__injection
  public :: boundary_shock__curre
  public :: boundary_shock__phi
  public :: boundary_shock__mom

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
  integer, save :: jup, jdown, kup, kdown, mnpi, mnpr, ncomw
  integer       :: nerr
  integer, allocatable :: nstat(:)
  real(8), save :: delx, delt, c, d_delx


contains

  subroutine boundary_shock__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nzgs_in,nzge_in,&
                                  nys_in,nye_in,nzs_in,nze_in, &
                                  jup_in,jdown_in,kup_in,kdown_in,mnpi_in,mnpr_in,ncomw_in,nerr_in,nstat_in, &
                                  delx_in,delt_in,c_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nzgs_in, nzge_in, nys_in, nye_in, nzs_in, nze_in
    integer, intent(in) :: jup_in, jdown_in, kup_in, kdown_in, mnpi_in, mnpr_in, ncomw_in, nerr_in, nstat_in(:)
    real(8), intent(in) :: delx_in, delt_in, c_in

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
    jup   = jup_in
    jdown = jdown_in
    kup   = kup_in
    kdown = kdown_in
    mnpi  = mnpi_in
    mnpr  = mnpr_in
    ncomw = ncomw_in
    nerr  = nerr_in
    allocate(nstat(size(nstat_in)))
    nstat = nstat_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in

    d_delx = 1d0/delx

    is_init = .true.


  end subroutine boundary_shock__init


  subroutine boundary_shock__particle_yz(up,np2)

!$  use omp_lib

    integer, intent(inout) :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(inout) :: up(ndim,np,nys:nye,nzs:nze,nsp)
    logical, save              :: lflag=.true.
!$  integer(omp_lock_kind)     :: lck(nys-1:nye+1,nzs-1:nze+1)
    integer                    :: j, k, ii, iii, isp, jpos, kpos, idim
    integer                    :: cnt(nys-1:nye+1,nzs-1:nze+1), cnt2(nys:nye,nzs:nze)
    integer                    :: cnt_tmp_j(nzs-1:nze+1), cnt_tmp_k(nys:nye), cnt_tmp
    integer                    :: ssize, rsize
    integer                    :: bff_cnt_j(nze-nzs+3)
    integer, save, allocatable :: flag(:,:,:)
    real(8), save, allocatable :: bff_ptcl(:,:,:)
    real(8), allocatable       :: ptcl_snd(:), ptcl_rcv(:)

    if(.not.is_init)then
      write(6,*)'Initialize first by calling boundary__init()'
      stop
    endif

    if(lflag)then
      allocate(flag(np,nys:nye,nzs:nze))
      allocate(bff_ptcl(ndim*np,nys-1:nye+1,nzs-1:nze+1))
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
!$OMP END WORKSHARE
!$OMP WORKSHARE
      cnt2(nys:nye,nzs:nze) = 0
!$OMP END WORKSHARE

!$OMP DO PRIVATE(ii,j,k,jpos,kpos,idim)
      do k=nzs,nze
      do j=nys,nye
        do ii=1,np2(j,k,isp)

          jpos = int(up(2,ii,j,k,isp)*d_delx)
          kpos = int(up(3,ii,j,k,isp)*d_delx)

          if(.not.(jpos == j .and. kpos == k))then

            if(jpos <= nygs-1)then
              up(2,ii,j,k,isp) = up(2,ii,j,k,isp)+(nyge-nygs+1)*delx
            else if(jpos >= nyge+1)then
              up(2,ii,j,k,isp) = up(2,ii,j,k,isp)-(nyge-nygs+1)*delx
            endif

            if(kpos <= nzgs-1)then
              up(3,ii,j,k,isp) = up(3,ii,j,k,isp)+(nzge-nzgs+1)*delx
            else if(kpos >= nzge+1)then
              up(3,ii,j,k,isp) = up(3,ii,j,k,isp)-(nzge-nzgs+1)*delx
            endif

!$          call omp_set_lock(lck(jpos,kpos))
            do idim=1,ndim
              bff_ptcl(idim+ndim*cnt(jpos,kpos),jpos,kpos) = up(idim,ii,j,k,isp)
            enddo
            cnt(jpos,kpos) = cnt(jpos,kpos)+1
!$          call omp_unset_lock(lck(jpos,kpos))

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
      ssize = ndim*sum(bff_cnt_j(1:nze-nzs+3))
      rsize = ndim*sum(cnt_tmp_j(nzs-1:nze+1))
      allocate(ptcl_snd(max(ssize,1)))
      allocate(ptcl_rcv(max(rsize,1)))

      j=nys-1
!$OMP PARALLEL DO PRIVATE(iii,ii,k,idim)
      do k=nzs-1,nze+1
        do ii=1,cnt(j,k)
          iii = ndim*(ii+sum(cnt(j,nzs-1:k-1))-1)
          do idim=1,ndim
            ptcl_snd(idim+iii) = bff_ptcl(idim+ndim*(ii-1),j,k)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,jdown,101, &
                        ptcl_rcv(1),rsize,mnpr,jup  ,101, &
                        ncomw,nstat,nerr)

      j=nye
!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,k,idim)
      do k=nzs-1,nze+1
        do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_j(k)
          iii = ndim*(ii-1-cnt(j,k)+sum(cnt_tmp_j(nzs-1:k-1)))
          do idim=1,ndim
            bff_ptcl(idim+ndim*(ii-1),j,k) = ptcl_rcv(idim+iii)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(k)
      do k=nzs-1,nze+1
        cnt(nye,k) = cnt(nye,k)+cnt_tmp_j(k)
      enddo
!$OMP END DO

!$OMP END PARALLEL

      deallocate(ptcl_snd)
      deallocate(ptcl_rcv)


      !!TRANSFER TO JUP and RECEIVE FROM JDOWN
      bff_cnt_j(1:nze-nzs+3) = cnt(nye+1,nzs-1:nze+1)
      call MPI_SENDRECV(bff_cnt_j(1)    ,nze-nzs+3,mnpi,jup  ,200, &
                        cnt_tmp_j(nzs-1),nze-nzs+3,mnpi,jdown,200, &
                        ncomw,nstat,nerr)
      ssize = ndim*sum(bff_cnt_j(1:nze-nzs+3))
      rsize = ndim*sum(cnt_tmp_j(nzs-1:nze+1))
      allocate(ptcl_snd(max(ssize,1)))
      allocate(ptcl_rcv(max(rsize,1)))

      j=nye+1
!$OMP PARALLEL DO PRIVATE(iii,ii,k,idim)
      do k=nzs-1,nze+1
        do ii=1,cnt(j,k)
          iii = ndim*(ii+sum(cnt(j,nzs-1:k-1))-1)
          do idim=1,ndim
            ptcl_snd(idim+iii) = bff_ptcl(idim+ndim*(ii-1),j,k)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,jup  ,201, &
                        ptcl_rcv(1),rsize,mnpr,jdown,201, &
                        ncomw,nstat,nerr)

      j=nys
!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,k,idim)
      do k=nzs-1,nze+1
        do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_j(k)
          iii = ndim*(ii-1-cnt(j,k)+sum(cnt_tmp_j(nzs-1:k-1)))
          do idim=1,ndim
            bff_ptcl(idim+ndim*(ii-1),j,k) = ptcl_rcv(idim+iii)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(k)
      do k=nzs-1,nze+1
        cnt(nys,k) = cnt(nys,k)+cnt_tmp_j(k)
      enddo
!$OMP END DO

!$OMP END PARALLEL

      deallocate(ptcl_snd)
      deallocate(ptcl_rcv)


      !!TRANSFER TO KDOWN and RECEIVE FROM KUP
      call MPI_SENDRECV(cnt(nys,nzs-1),nye-nys+1,mnpi,kdown,300, &
                        cnt_tmp_k(nys),nye-nys+1,mnpi,kup  ,300, &
                        ncomw,nstat,nerr)
      ssize = ndim*sum(cnt(nys:nye,nzs-1))
      rsize = ndim*sum(cnt_tmp_k(nys:nye))
      allocate(ptcl_snd(max(ssize,1)))
      allocate(ptcl_rcv(max(rsize,1)))

      k=nzs-1
!$OMP PARALLEL DO PRIVATE(iii,ii,j,idim)
      do j=nys,nye
        do ii=1,cnt(j,k)
          iii = 7*(ii+sum(cnt(nys:j-1,k))-1)
          do idim=1,ndim
            ptcl_snd(idim+iii) = bff_ptcl(idim+ndim*(ii-1),j,k)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,kdown,301, &
                        ptcl_rcv(1),rsize,mnpr,kup  ,301, &
                        ncomw,nstat,nerr)

      k=nze
!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,j,idim)
      do j=nys,nye
        do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_k(j)
          iii = ndim*(ii-1-cnt(j,k)+sum(cnt_tmp_k(nys:j-1)))
          do idim=1,ndim
            bff_ptcl(idim+ndim*(ii-1),j,k) = ptcl_rcv(idim+iii)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(j)
      do j=nys,nye
        cnt(j,nze) = cnt(j,nze)+cnt_tmp_k(j)
      enddo
!$OMP END DO

!$OMP END PARALLEL

      deallocate(ptcl_snd)
      deallocate(ptcl_rcv)


      !!TRANSFER TO KUP and RECEIVE FROM KDOWN
      call MPI_SENDRECV(cnt(nys,nze+1),nye-nys+1,mnpi,kup  ,400, &
                        cnt_tmp_k(nys),nye-nys+1,mnpi,kdown,400, &
                        ncomw,nstat,nerr)

      ssize = ndim*sum(cnt(nys:nye,nze+1))
      rsize = ndim*sum(cnt_tmp_k(nys:nye))
      allocate(ptcl_snd(max(ssize,1)))
      allocate(ptcl_rcv(max(rsize,1)))

      k=nze+1
!$OMP PARALLEL DO PRIVATE(iii,ii,j,idim)
      do j=nys,nye
        do ii=1,cnt(j,k)
          iii = 7*(ii+sum(cnt(nys:j-1,k))-1)
          do idim=1,ndim
            ptcl_snd(idim+iii) = bff_ptcl(idim+ndim*(ii-1),j,k)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      call MPI_SENDRECV(ptcl_snd(1),ssize,mnpr,kup  ,401, &
                        ptcl_rcv(1),rsize,mnpr,kdown,401, &
                        ncomw,nstat,nerr)

      k=nzs
!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,j,idim)
      do j=nys,nye
        do ii=cnt(j,k)+1,cnt(j,k)+cnt_tmp_k(j)
          iii = ndim*(ii-1-cnt(j,k)+sum(cnt_tmp_k(nys:j-1)))
          do idim=1,ndim
            bff_ptcl(idim+ndim*(ii-1),j,k) = ptcl_rcv(idim+iii)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(j)
      do j=nys,nye
        cnt(j,nzs) = cnt(j,nzs)+cnt_tmp_k(j)
      enddo
!$OMP END DO

!$OMP END PARALLEL

      deallocate(ptcl_snd)
      deallocate(ptcl_rcv)


!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,j,k,idim,cnt_tmp)
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
            do idim=1,ndim
              up(idim,flag(ii,j,k),j,k,isp) = up(idim,np2(j,k,isp),j,k,isp)
            enddo
            np2(j,k,isp) = np2(j,k,isp)-1
          else
            do idim=1,ndim
              up(idim,flag(ii,j,k),j,k,isp) = bff_ptcl(idim+ndim*iii,j,k)
            enddo
            iii = iii+1
            cnt(j,k) = cnt(j,k)-1
          endif
        enddo loop1

        if(cnt(j,k) > 0)then
          do ii=1,cnt(j,k)
            do idim=1,ndim
              up(idim,np2(j,k,isp)+ii,j,k,isp) = bff_ptcl(ndim*iii+idim+ndim*(ii-1),j,k)
            enddo
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

  end subroutine boundary_shock__particle_yz


  subroutine boundary_shock__injection(up,np2,nxs,nxe,u0)

    real(8), intent(inout) :: up(ndim,np,nys:nye,nzs:nze,nsp)
    integer, intent(in)    :: np2(nys:nye,nzs:nze,nsp)
    integer, intent(in)    :: nxs, nxe
    real(8), intent(in)    :: u0

    integer                :: j, k, ii, isp, ipos
    real(8)                :: xend

    if(.not.is_init)then
      write(6,*)'Initialize first by calling boundary_shock__init()'
      stop
    endif

    xend = nxe*delx+u0/sqrt(1d0+(u0*u0)/(c*c))*delt

    do isp=1,nsp

!$OMP PARALLEL DO PRIVATE(ii,j,k,ipos)
       do k=nzs,nze
       do j=nys,nye
          do ii=1,np2(j,k,isp)

             ipos = int(up(1,ii,j,k,isp)*d_delx)

             if(ipos < nxs+1)then
                up(1,ii,j,k,isp) = 2d0*(nxs+1)*delx-up(1,ii,j,k,isp)
                up(4,ii,j,k,isp) = -up(4,ii,j,k,isp)
                up(5,ii,j,k,isp) = -up(5,ii,j,k,isp)
                up(6,ii,j,k,isp) = -up(6,ii,j,k,isp)
             else if(up(1,ii,j,k,isp) > xend)then
                up(1,ii,j,k,isp) = 2d0*xend-up(1,ii,j,k,isp)
                up(4,ii,j,k,isp) = 2d0*u0-up(4,ii,j,k,isp)
                up(5,ii,j,k,isp) = -up(5,ii,j,k,isp)
                up(6,ii,j,k,isp) = -up(6,ii,j,k,isp)
             endif

          enddo
       enddo
       enddo
!$OMP END PARALLEL DO

    enddo

  end subroutine boundary_shock__injection


  subroutine boundary_shock__dfield(df,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)

    integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, nxgs, nxge
    real(8), intent(inout) :: df(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    integer                :: i, j, k, ii
    real(8)                :: bff_snd_j(12*(nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_rcv_j(12*(nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_snd_k(12*(nxe-nxs+1)*(nye-nys+5))
    real(8)                :: bff_rcv_k(12*(nxe-nxs+1)*(nye-nys+5))

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(k-nzs))

       bff_snd_j(ii+1)  = df(1,i,nys,k)
       bff_snd_j(ii+2)  = df(2,i,nys,k)
       bff_snd_j(ii+3)  = df(3,i,nys,k)
       bff_snd_j(ii+4)  = df(4,i,nys,k)
       bff_snd_j(ii+5)  = df(5,i,nys,k)
       bff_snd_j(ii+6)  = df(6,i,nys,k)
       bff_snd_j(ii+7)  = df(1,i,nys+1,k)
       bff_snd_j(ii+8)  = df(2,i,nys+1,k)
       bff_snd_j(ii+9)  = df(3,i,nys+1,k)
       bff_snd_j(ii+10) = df(4,i,nys+1,k)
       bff_snd_j(ii+11) = df(5,i,nys+1,k)
       bff_snd_j(ii+12) = df(6,i,nys+1,k)
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
       df(1,i,nye+1,k) = bff_rcv_j(ii+1)
       df(2,i,nye+1,k) = bff_rcv_j(ii+2)
       df(3,i,nye+1,k) = bff_rcv_j(ii+3)
       df(4,i,nye+1,k) = bff_rcv_j(ii+4)
       df(5,i,nye+1,k) = bff_rcv_j(ii+5)
       df(6,i,nye+1,k) = bff_rcv_j(ii+6)
       df(1,i,nye+2,k) = bff_rcv_j(ii+7)
       df(2,i,nye+2,k) = bff_rcv_j(ii+8)
       df(3,i,nye+2,k) = bff_rcv_j(ii+9)
       df(4,i,nye+2,k) = bff_rcv_j(ii+10)
       df(5,i,nye+2,k) = bff_rcv_j(ii+11)
       df(6,i,nye+2,k) = bff_rcv_j(ii+12)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(k-nzs))

       bff_snd_j(ii+1)  = df(1,i,nye-1,k)
       bff_snd_j(ii+2)  = df(2,i,nye-1,k)
       bff_snd_j(ii+3)  = df(3,i,nye-1,k)
       bff_snd_j(ii+4)  = df(4,i,nye-1,k)
       bff_snd_j(ii+5)  = df(5,i,nye-1,k)
       bff_snd_j(ii+6)  = df(6,i,nye-1,k)
       bff_snd_j(ii+7)  = df(1,i,nye,k)
       bff_snd_j(ii+8)  = df(2,i,nye,k)
       bff_snd_j(ii+9)  = df(3,i,nye,k)
       bff_snd_j(ii+10) = df(4,i,nye,k)
       bff_snd_j(ii+11) = df(5,i,nye,k)
       bff_snd_j(ii+12) = df(6,i,nye,k)
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

       df(1,i,nys-2,k) = bff_rcv_j(ii+1)
       df(2,i,nys-2,k) = bff_rcv_j(ii+2)
       df(3,i,nys-2,k) = bff_rcv_j(ii+3)
       df(4,i,nys-2,k) = bff_rcv_j(ii+4)
       df(5,i,nys-2,k) = bff_rcv_j(ii+5)
       df(6,i,nys-2,k) = bff_rcv_j(ii+6)
       df(1,i,nys-1,k) = bff_rcv_j(ii+7)
       df(2,i,nys-1,k) = bff_rcv_j(ii+8)
       df(3,i,nys-1,k) = bff_rcv_j(ii+9)
       df(4,i,nys-1,k) = bff_rcv_j(ii+10)
       df(5,i,nys-1,k) = bff_rcv_j(ii+11)
       df(6,i,nys-1,k) = bff_rcv_j(ii+12)
    enddo
    enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-2,nye+2
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(j-(nys-2)))

       bff_snd_k(ii+1)  = df(1,i,j,nzs)
       bff_snd_k(ii+2)  = df(2,i,j,nzs)
       bff_snd_k(ii+3)  = df(3,i,j,nzs)
       bff_snd_k(ii+4)  = df(4,i,j,nzs)
       bff_snd_k(ii+5)  = df(5,i,j,nzs)
       bff_snd_k(ii+6)  = df(6,i,j,nzs)
       bff_snd_k(ii+7)  = df(1,i,j,nzs+1)
       bff_snd_k(ii+8)  = df(2,i,j,nzs+1)
       bff_snd_k(ii+9)  = df(3,i,j,nzs+1)
       bff_snd_k(ii+10) = df(4,i,j,nzs+1)
       bff_snd_k(ii+11) = df(5,i,j,nzs+1)
       bff_snd_k(ii+12) = df(6,i,j,nzs+1)
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

       df(1,i,j,nze+1) = bff_rcv_k(ii+1)
       df(2,i,j,nze+1) = bff_rcv_k(ii+2)
       df(3,i,j,nze+1) = bff_rcv_k(ii+3)
       df(4,i,j,nze+1) = bff_rcv_k(ii+4)
       df(5,i,j,nze+1) = bff_rcv_k(ii+5)
       df(6,i,j,nze+1) = bff_rcv_k(ii+6)
       df(1,i,j,nze+2) = bff_rcv_k(ii+7)
       df(2,i,j,nze+2) = bff_rcv_k(ii+8)
       df(3,i,j,nze+2) = bff_rcv_k(ii+9)
       df(4,i,j,nze+2) = bff_rcv_k(ii+10)
       df(5,i,j,nze+2) = bff_rcv_k(ii+11)
       df(6,i,j,nze+2) = bff_rcv_k(ii+12)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-2,nye+2
    do i=nxs,nxe
       ii = 12*((i-nxs)+(nxe-nxs+1)*(j-(nys-2)))

       bff_snd_k(ii+1)  = df(1,i,j,nze-1)
       bff_snd_k(ii+2)  = df(2,i,j,nze-1)
       bff_snd_k(ii+3)  = df(3,i,j,nze-1)
       bff_snd_k(ii+4)  = df(4,i,j,nze-1)
       bff_snd_k(ii+5)  = df(5,i,j,nze-1)
       bff_snd_k(ii+6)  = df(6,i,j,nze-1)
       bff_snd_k(ii+7)  = df(1,i,j,nze)
       bff_snd_k(ii+8)  = df(2,i,j,nze)
       bff_snd_k(ii+9)  = df(3,i,j,nze)
       bff_snd_k(ii+10) = df(4,i,j,nze)
       bff_snd_k(ii+11) = df(5,i,j,nze)
       bff_snd_k(ii+12) = df(6,i,j,nze)
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

       df(1,i,j,nzs-2) = bff_rcv_k(ii+1)
       df(2,i,j,nzs-2) = bff_rcv_k(ii+2)
       df(3,i,j,nzs-2) = bff_rcv_k(ii+3)
       df(4,i,j,nzs-2) = bff_rcv_k(ii+4)
       df(5,i,j,nzs-2) = bff_rcv_k(ii+5)
       df(6,i,j,nzs-2) = bff_rcv_k(ii+6)
       df(1,i,j,nzs-1) = bff_rcv_k(ii+7)
       df(2,i,j,nzs-1) = bff_rcv_k(ii+8)
       df(3,i,j,nzs-1) = bff_rcv_k(ii+9)
       df(4,i,j,nzs-1) = bff_rcv_k(ii+10)
       df(5,i,j,nzs-1) = bff_rcv_k(ii+11)
       df(6,i,j,nzs-1) = bff_rcv_k(ii+12)
    enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-2,nze+2
    do j=nys-2,nye+2
       df(1  ,nxs-1,j,k) = -df(1,  nxs,  j,k)
       df(2:4,nxs-1,j,k) =  df(2:4,nxs+1,j,k)
       df(5:6,nxs-1,j,k) = -df(5:6,nxs  ,j,k)
       df(1  ,nxe+1,j,k) = 0d0
       df(2:4,nxe+1,j,k) = 0d0
       df(5:6,nxe+1,j,k) = 0d0
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary_shock__dfield


  subroutine boundary_shock__curre(uj,nxs,nxe,nys,nye,nzs,nze,nxgs,nxge)

    integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, nxgs, nxge
    real(8), intent(inout) :: uj(3,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    integer                :: i, j, k, ii
    real(8)                :: bff_snd_j(6*(nxe-nxs+5)*(nze-nzs+5))
    real(8)                :: bff_rcv_j(6*(nxe-nxs+5)*(nze-nzs+5))
    real(8)                :: bff_snd_k(6*(nxe-nxs+5)*(nye-nys+5))
    real(8)                :: bff_rcv_k(6*(nxe-nxs+5)*(nye-nys+5))

!--- in y-direction

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

!--- in z-direction

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

!--- update of nori-shiro --- !

!--- in y-direction

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(k-nzs))

       bff_snd_j(ii+1) = uj(1,i,nys,k)
       bff_snd_j(ii+2) = uj(2,i,nys,k)
       bff_snd_j(ii+3) = uj(3,i,nys,k)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_j(1),3*(nxe-nxs+5)*(nze-nzs+1),mnpr,jdown,110, &
                      bff_rcv_j(1),3*(nxe-nxs+5)*(nze-nzs+1),mnpr,jup  ,110, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(k-nzs))

       uj(1,i,nye+1,k) = bff_rcv_j(ii+1)
       uj(2,i,nye+1,k) = bff_rcv_j(ii+2)
       uj(3,i,nye+1,k) = bff_rcv_j(ii+3)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(k-nzs))

       bff_snd_j(ii+1) = uj(1,i,nye,k)
       bff_snd_j(ii+2) = uj(2,i,nye,k)
       bff_snd_j(ii+3) = uj(3,i,nye,k)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_j(1),3*(nxe-nxs+5)*(nze-nzs+1),mnpr,jup  ,100, &
                      bff_rcv_j(1),3*(nxe-nxs+5)*(nze-nzs+1),mnpr,jdown,100, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,k,ii)
    do k=nzs,nze
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(k-nzs))

       uj(1,i,nys-1,k) = bff_rcv_j(ii+1)
       uj(2,i,nys-1,k) = bff_rcv_j(ii+2)
       uj(3,i,nys-1,k) = bff_rcv_j(ii+3)
    enddo
    enddo
!$OMP END PARALLEL DO

!--- in z-direction

!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(j-(nys-1)))

       bff_snd_k(ii+1) = uj(1,i,j,nzs)
       bff_snd_k(ii+2) = uj(2,i,j,nzs)
       bff_snd_k(ii+3) = uj(3,i,j,nzs)
    enddo
    enddo
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd_k(1),3*(nxe-nxs+5)*(nye-nys+3),mnpr,kdown,210, &
                      bff_rcv_k(1),3*(nxe-nxs+5)*(nye-nys+3),mnpr,kup  ,210, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(j-(nys-1)))

       uj(1,i,j,nze+1) = bff_rcv_k(ii+1)
       uj(2,i,j,nze+1) = bff_rcv_k(ii+2)
       uj(3,i,j,nze+1) = bff_rcv_k(ii+3)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(j-(nys-1)))

       bff_snd_k(ii+1) = uj(1,i,j,nze)
       bff_snd_k(ii+2) = uj(2,i,j,nze)
       bff_snd_k(ii+3) = uj(3,i,j,nze)
    enddo
    enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd_k(1),3*(nxe-nxs+5)*(nye-nys+3),mnpr,kup  ,200, &
                      bff_rcv_k(1),3*(nxe-nxs+5)*(nye-nys+3),mnpr,kdown,200, &
                      ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,j,ii)
    do j=nys-1,nye+1
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2)+(nxe-nxs+5)*(j-(nys-1)))

       uj(1,i,j,nzs-1) = bff_rcv_k(ii+1)
       uj(2,i,j,nzs-1) = bff_rcv_k(ii+2)
       uj(3,i,j,nzs-1) = bff_rcv_k(ii+3)
    enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine boundary_shock__curre


  subroutine boundary_shock__phi(phi,nxs,nxe,nys,nye,nzs,nze,l)

    integer, intent(in)    :: nxs, nxe, nys, nye, nzs, nze, l
    real(8), intent(inout) :: phi(nxs-1:nxe+1,nys-1:nye+1,nzs-1:nze+1)
    integer                :: i, j, k, ii
    real(8)                :: bff_snd_j((nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_rcv_j((nxe-nxs+1)*(nze-nzs+1))
    real(8)                :: bff_snd_k((nxe-nxs+1)*(nye-nys+3))
    real(8)                :: bff_rcv_k((nxe-nxs+1)*(nye-nys+3))

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

    select case(l)
    case(1)

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-1,nze+1
    do j=nys-1,nye+1
       phi(nxs-1,j,k) = -phi(nxs,j,k)
       phi(nxe+1,j,k) = 0d0
    enddo
    enddo
!$OMP END PARALLEL DO

    case(2,3)

!$OMP PARALLEL DO PRIVATE(j,k)
    do k=nzs-1,nze+1
    do j=nys-1,nye+1
       phi(nxs-1,j,k) = phi(nxs+1,j,k)
       phi(nxe+1,j,k) = 0d0
    enddo
    enddo
!$OMP END PARALLEL DO

    end select

  end subroutine boundary_shock__phi


  subroutine boundary_shock__mom(mom)

    real(8), intent(inout) :: mom(7,nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,nsp)
    integer, parameter     :: nl = 7
    integer                :: i, j, k, l, ii, isp
    real(8)                :: bff_snd_j(nl*(nxge-nxgs+3)*(nze-nzs+3))
    real(8)                :: bff_rcv_j(nl*(nxge-nxgs+3)*(nze-nzs+3))
    real(8)                :: bff_snd_k(nl*(nxge-nxgs+3)*(nye-nys+1))
    real(8)                :: bff_rcv_k(nl*(nxge-nxgs+3)*(nye-nys+1))


!$OMP PARALLEL WORKSHARE
    mom(1:nl,nxgs,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(1:nl,nxgs  ,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                  +mom(1:nl,nxgs-1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
    mom(1:nl,nxge,nys-1:nye+1,nzs-1:nze+1,1:nsp) = mom(1:nl,nxge  ,nys-1:nye+1,nzs-1:nze+1,1:nsp) &
                                                  +mom(1:nl,nxge+1,nys-1:nye+1,nzs-1:nze+1,1:nsp)
!$OMP END PARALLEL WORKSHARE

    do isp = 1,nsp
!$OMP PARALLEL DO PRIVATE(i,k,l,ii)
      do k=nzs-1,nze+1
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(k-(nzs-1)))
        do l=1,nl
          bff_snd_j(ii+l) = mom(l,i,nys-1,k,isp)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_SENDRECV(bff_snd_j(1),nl*(nxge-nxgs+3)*(nze-nzs+3),mnpr,jdown,110, &
                        bff_rcv_j(1),nl*(nxge-nxgs+3)*(nze-nzs+3),mnpr,jup  ,110, &
                        ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,k,l,ii)
      do k=nzs-1,nze+1
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(k-(nzs-1)))
        do l=1,nl
          mom(l,i,nye,k,isp) = mom(l,i,nye,k,isp)+bff_rcv_j(ii+l)
        enddo
      enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,l,ii)
      do k=nzs-1,nze+1
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(k-(nzs-1)))
        do l=1,nl
          bff_snd_j(ii+l) = mom(l,i,nye+1,k,isp)
        enddo
      enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      call MPI_SENDRECV(bff_snd_j(1),nl*(nxge-nxgs+3)*(nze-nzs+3),mnpr,jup  ,100, &
                        bff_rcv_j(1),nl*(nxge-nxgs+3)*(nze-nzs+3),mnpr,jdown,100, &
                        ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,k,l,ii)
      do k=nzs-1,nze+1
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(k-(nzs-1)))
        do l=1,nl
          mom(l,i,nys,k,isp) = mom(l,i,nys,k,isp)+bff_rcv_j(ii+l)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j,l,ii)
      do j=nys,nye
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(j-nys))
        do l=1,nl
          bff_snd_k(ii+l) = mom(l,i,j,nzs-1,isp)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_SENDRECV(bff_snd_k(1),nl*(nxge-nxgs+3)*(nye-nys+1),mnpr,kdown,210, &
                        bff_rcv_k(1),nl*(nxge-nxgs+3)*(nye-nys+1),mnpr,kup  ,210, &
                        ncomw,nstat,nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,j,l,ii)
      do j=nys,nye
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(j-nys))
        do l=1,nl
          mom(l,i,j,nze,isp) = mom(l,i,j,nze,isp)+bff_rcv_k(ii+l)
        enddo
      enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,j,l,ii)
      do j=nys,nye
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(j-nys))
        do l=1,nl
          bff_snd_k(ii+l) = mom(l,i,j,nze+1,isp)
        enddo
      enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      call MPI_SENDRECV(bff_snd_k(1),nl*(nxge-nxgs+3)*(nye-nys+1),mnpr,kup  ,200, &
                        bff_rcv_k(1),nl*(nxge-nxgs+3)*(nye-nys+1),mnpr,kdown,200, &
                        ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,j,l,ii)
      do j=nys,nye
      do i=nxgs-1,nxge+1
        ii = nl*((i-(nxgs-1))+(nxge-nxgs+3)*(j-nys))
        do l=1,nl
          mom(l,i,j,nzs,isp) = mom(l,i,j,nzs,isp)+bff_rcv_k(ii+l)
        enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo

  end subroutine boundary_shock__mom



end module boundary_shock
