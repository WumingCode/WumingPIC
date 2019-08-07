program main

  use const
  use mpi_set
  use init
  use boundary, only : boundary__particle_yz, &
                       boundary__mom, boundary__particle_injection
  use fio, only : fio__mom, fio__output
  use particle
  use field
  use sort, only : sort__bucket
  use mom_calc
  use omp_lib

  implicit none

  integer :: it=0
  real(8) :: etime, etime0

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
!    3-D code w SIMD     (by Y. Matsumoto, Chiba-U) 2013/4/1
!
!**********************************************************************c

  etime0 = omp_get_wtime()

  call init__set_param
  call MPI_BCAST(etime0,1,mnpr,nroot,ncomw,nerr)

  loop: do it=1,itmax-it0

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
        if(nrank == nroot) write(*,*) '*** elapse time over ***',it-1+it0,etime-etime0
        exit loop
     endif

     call particle__solv(gp,                                              &
                         nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                         c,q,r,delt,delx,                                 &
                         up,uf)
     call boundary__particle_injection(gp,                                        &
                                       nxs,nxe,nys,nye,nzs,nze,np,nsp,np2,        &
                                       mod(it+it0-1,intvl2),mod(it+it0-1,intvl3), &
                                       delx,delt,u0,c)
     call field__fdtd_i(uf,up,gp,                                        &
                        nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                        jup,jdown,kup,kdown,mnpr,opsum,nstat,ncomw,nerr, &
                        q,c,delx,delt,gfac)
     call boundary__particle_yz(gp,                                  &
                                nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                                np,nsp,np2,delx,                     &
                                jup,jdown,kup,kdown,nstat,mnpi,mnpr,ncomw,nerr)

     call sort__bucket(up,gp,cumcnt,nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,np2,delx)

     if(mod(it+it0,intvl2) == 0) call init__inject(it+it0)
     if(mod(it+it0,intvl3) == 0) call init__relocate(it+it0)

     if(mod(it+it0,intvl1) == 0)then 
        call fio__output(.false.,                                               &
                         nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
                         np,nsp,np2,it,it0,                                     &
                         nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
                         c,q,r,delt,delx,dir,                                   &
                         up,uf)
     endif

     if(mod(it+it0,intvl4) == 0)then 
        call mom_calc__accl(gp,                                              &
                            nxgs,nxge,nxs,nxe,nys,nye,nzs,nze,np,nsp,cumcnt, &
                            c,q,r,0.5*delt,delx,                             &
                            up,uf)
        call mom_calc__nvt(den,vel,temp,gp,nxgs,nxge,nys,nye,nzs,nze,np,nsp,np2,c)
        call boundary__mom(den,vel,temp,nxgs,nxge,nys,nye,nzs,nze,nsp)
        call fio__mom(den,vel,temp,uf,nxgs,nxge,nys,nye,nzs,nze,nsp,it+it0,nrank,dir)
     endif

  enddo loop

  call MPI_FINALIZE(nerr)

end program main

