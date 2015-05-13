program main

  use const
  use mpi_set
  use init
  use boundary
  use fio
  use particle
  use field
  use sort, only : sort__bucket
  use mom_calc

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
!  etlim = 72.*60.*60.-10.*60.
  etlim = 24.*60.*60.-30.*60.
!  etlim = 7.*60.*60.+30.*60.
  !Test runs
!  etlim = 1.*60.*60.
!  etlim = 15.*60.
  !*****************************!

  etime0 = omp_get_wtime()

  call init__set_param
  call MPI_BCAST(etime0,1,mnpr,nroot,ncomw,nerr)

  loop: do it=1,itmax-it0,2

     if(nrank == nroot) etime = omp_get_wtime()
     call MPI_BCAST(etime,1,mnpr,nroot,ncomw,nerr)

     if(etime-etime0 >= etlim) then
        call fio__output(.true.,                                                &
                         nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
                         np,nsp,np2,it-1,it0,                                   &
                         nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
                         c,q,r,delt,delx,dir,                                   &
                         up,uf)
        call mom_calc__accl(gp,                                              &
                            nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                            c,q,r,0.5*delt,delx,                             &
                            up,uf)
        call mom_calc__nvt(den,vel,temp,gp,nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2,c)
        call boundary__mom(den,vel,temp,nxgs,nxge,nys,nye,nzs,nze,nsp)
        call fio__mom(den,vel,temp,uf,nxgs,nxge,nys,nye,nzs,nze,nsp,it-1+it0,nrank,dir)

        if(nrank == nroot) write(*,*) '*** elapse time over ***',it,etime-etime0
        exit loop
     endif

     call particle__solv(gp,                                              &
                         nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                         c,q,r,delt,delx,                                 &
                         up,uf)
     call boundary__particle_x(gp,                                 &
                               nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)
     call field__fdtd_i(uf,up,gp,                                       &
                        nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                        jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                        q,c,delx,delt,gfac)
     call boundary__particle_yz(up,                                  &
                                nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                                np,nsp,np2,                          &
                                jup,jdown,kup,kdown,nstat,mnpi,mnpr,ncomw,nerr)

     if(mod(it+it0,intvl2) == 0) call init__inject(up)
     if(mod(it+it0,intvl3) == 0) call init__relocate(up)

     call sort__bucket(gp,up,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)

     if(mod(it+it0,intvl1) == 0)then 
        call fio__output(.false.,                                               &
                         nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
                         np,nsp,np2,it,it0,                                     &
                         nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
                         c,q,r,delt,delx,dir,                                   &
                         gp,uf)
     endif

     if(mod(it+it0,intvl4) == 0)then 
        call mom_calc__accl(up,                                              &
                            nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                            c,q,r,0.5*delt,delx,                             &
                            gp,uf)
        call mom_calc__nvt(den,vel,temp,up,nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2,c)
        call boundary__mom(den,vel,temp,nxgs,nxge,nys,nye,nzs,nze,nsp)
        call fio__mom(den,vel,temp,uf,nxgs,nxge,nys,nye,nzs,nze,nsp,it+it0,nrank,dir)
     endif

     !it=it+1

     call particle__solv(up,                                              &
                         nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                         c,q,r,delt,delx,                                 &
                         gp,uf)
     call boundary__particle_x(up,                                 &
                               nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)
     call field__fdtd_i(uf,gp,up,                                        &
                        nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                        jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                        q,c,delx,delt,gfac)
     call boundary__particle_yz(gp,                                  &
                                nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                                np,nsp,np2,                          &
                                jup,jdown,kup,kdown,nstat,mnpi,mnpr,ncomw,nerr)

     if(mod(it+1+it0,intvl2) == 0) call init__inject(gp)
     if(mod(it+1+it0,intvl3) == 0) call init__relocate(gp)

     call sort__bucket(up,gp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,np2)

     if(mod(it+1+it0,intvl1) == 0)then
        call fio__output(.false.,                                               &
                         nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
                         np,nsp,np2,it+1,it0,                                   &
                         nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
                         c,q,r,delt,delx,dir,                                   &
                         up,uf)
     endif

     if(mod(it+1+it0,intvl4) == 0)then
        call mom_calc__accl(gp,                                              &
                            nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                            c,q,r,0.5*delt,delx,                             &
                            up,uf)
        call mom_calc__nvt(den,vel,temp,gp,nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2,c)
        call boundary__mom(den,vel,temp,nxgs,nxge,nys,nye,nzs,nze,nsp)
        call fio__mom(den,vel,temp,uf,nxgs,nxge,nys,nye,nzs,nze,nsp,it+1+it0,nrank,dir)
     endif

  enddo loop

  call MPI_FINALIZE(nerr)

end program main

