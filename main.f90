program main

  use const
  use mpi_set
  use init
!!$  use boundary
!!$  use fio
!!$  use particle
!!$  use field
  implicit none

  integer :: it=0
  real(8) :: etime, etlim, etime0, omp_get_wtime

!**********************************************************************c
!
!    two-dimensional electromagnetic plasma simulation code
!
!    written by M Hoshino,  ISAS, 1984/09/12
!    revised  1985/03/08  1985/04/05  1997/05/06
!    revised for CANS    (by Y. Matsumoto, STEL)  2004/06/22
!    re-written in F90   (by Y. Matsumoto, STEL)  2008/10/21
!    MPI parallelization (by Y. Matsumoto, STEL)  2009/4/1
!    2-D code            (by Y. Matsumoto, STEL)  2009/6/5
!    3-D code            (by Y. Matsumoto, Chiba-U) 2012/10/5
!
!**********************************************************************c

  !**** Maximum elapse time ****!o
  etlim = 48.*60.*60.-5.*60.
!!$  etlim = 20000.-20.*60.
  !Test runs
!!$  etlim = 15.*60.
!!$  !*****************************!
!!$  etime0 = omp_get_wtime()

  call init__set_param
!!$  call MPI_BCAST(etime0,1,mnpr,nroot,ncomw,nerr)
!!$
!!$  loop: do it=1,itmax-it0
!!$
!!$     if(nrank == nroot) etime = omp_get_wtime()
!!$
!!$     call MPI_BCAST(etime,1,mnpr,nroot,ncomw,nerr)
!!$
!!$     if(etime-etime0 >= etlim) then
!!$        call fio__output(up,uf,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,np2,nproc,nrank, &
!!$                         c,q,r,delt,delx,it-1,it0,dir,.true.)
!!$        if(nrank == nroot) write(*,*) '*** elapse time over ***',it,etime-etime0
!!$        exit loop
!!$     endif
!!$
!!$     call particle__solv(gp,up,uf,                     &
!!$                         np,nsp,np2,nxgs,nxge,nys,nye, &
!!$                         c,q,r,delt,delx)
!!$     call field__fdtd_i(uf,up,gp,                             &
!!$                        np,nsp,np2,nxgs,nxge,nxs,nxe,nys,nye, &
!!$                        q,c,delx,delt,gfac,                   &
!!$                        nup,ndown,mnpr,opsum,nstat,ncomw,nerr)
!!$     call boundary__particle(up,                                   &
!!$                             np,nsp,np2,nygs,nyge,nxs,nxe,nys,nye, &
!!$                             nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)
!!$
!!$     if(mod(it+it0,intvl3) == 0) call init__inject
!!$
!!$     if(mod(it+it0,intvl1) == 0)                                                             &
!!$          call fio__output(up,uf,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,np2,nproc,nrank, &
!!$                           c,q,r,delt,delx,it,it0,dir,.false.)
!!$     if(mod(it+it0,intvl4) == 0) call init__relocate
!!$
!!$  enddo loop
!!$
!!$  call MPI_FINALIZE(nerr)

end program main

