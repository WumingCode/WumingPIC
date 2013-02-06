program main

  use const
  use mpi_set
  use init
  use boundary
  use fio
  use particle
  use field

  implicit none

  integer :: it=0
  real(8) :: etime, etlim, etime0, omp_get_wtime

!**********************************************************************c
!
!    three-dimensional electromagnetic plasma simulation code
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
!  etlim = 48.*60.*60.-5.*60.
  !Test runs
  etlim = 5.*60.
  !*****************************!

!$  etime0 = omp_get_wtime()

  call init__set_param
  call MPI_BCAST(etime0,1,mnpr,nroot,ncomw,nerr)

  loop: do it=1,itmax-it0

!$     if(nrank == nroot) etime = omp_get_wtime()
     call MPI_BCAST(etime,1,mnpr,nroot,ncomw,nerr)

     if(etime-etime0 >= etlim) then
!        call fio__output(.true.,                                                &
!                         nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
!                         np,nsp,np2,it,it0,                                     &
!                         nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
!                         c,q,r,delt,delx,dir,                                   &
!                         up,uf)
        if(nrank == nroot) write(*,*) '*** elapse time over ***',it,etime-etime0
        exit loop
     endif

     call fapp_start("ptcl",1,1)
     call particle__solv(gp,                                   &
                         nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2, &
                         c,q,r,delt,delx,                      &
                         up,uf)
     call fapp_stop("ptcl",1,1)

     call fapp_start("bnd1",1,1)
     call boundary__particle_x(gp,                                 &
                                nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)
     call fapp_stop("bnd1",1,1)

     call fapp_start("fld",1,1)
     call field__fdtd_i(uf,up,gp,                                        &
                        nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,np2,    &
                        jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                        q,c,delx,delt,gfac)
     call fapp_stop("fld",1,1)

     call fapp_start("bnd2",1,1)
     call boundary__particle_y(up,                                  &
                               nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                               np,nsp,np2,                          &
                               jup,jdown,kup,kdown,nstat,mnpi,mnpr,ncomw,nerr)
     call fapp_stop("bnd2",1,1)

     if(mod(it+it0,intvl2) == 0) call init__inject

     if(mod(it+it0,intvl1) == 0)                                                  &
          call fio__output(.false.,                                               &
                           nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
                           np,nsp,np2,it,it0,                                     &
                           nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
                           c,q,r,delt,delx,dir,                                   &
                           up,uf)

     if(mod(it+it0,intvl3) == 0) call init__relocate

  enddo loop

  call MPI_FINALIZE(nerr)

end program main

