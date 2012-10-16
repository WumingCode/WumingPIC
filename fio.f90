module fio

  implicit none

  private

  public :: fio__output
  public :: fio__input
  public :: fio__param


contains


  subroutine fio__output(lflag,                                                 &
                         nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze, &
                         np,nsp,np2,it,it0,                                     &
                         nproc,nproc_i,nproc_j,nproc_k,nrank,                   &
                         c,q,r,delt,delx,dir,                                   &
                         up,uf)

    logical, intent(in) :: lflag
    integer, intent(in) :: nxgs, nxge, nygs, nyge, nzgs, nzge, nxs, nxe, nys, nye, nzs, nze
    integer, intent(in) :: np, nsp, np2(nys:nye,nzs:nze,nsp)
    integer, intent(in) :: it, it0
    integer, intent(in) :: nproc, nproc_i, nproc_j, nproc_k, nrank
    real(8), intent(in) :: c, q(nsp), r(nsp), delt, delx
    real(8), intent(in) :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    character(len=*), intent(in) :: dir
    integer :: it2
    character(len=256) :: filename

    it2=it+it0

    !filename
    if(lflag)then
       write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),9999999,'_rank=',nrank,'.dat'
    else
       write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it2,'_rank=',nrank,'.dat'
    endif
    open(100+nrank,file=filename,form='unformatted')

    !time & parameters
    write(100+nrank)it2,nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze
    write(100+nrank)np,nsp,nproc,nproc_i,nproc_j,nproc_k
    write(100+nrank)delt,delx,c
    write(100+nrank)np2
    write(100+nrank)q
    write(100+nrank)r

    !field data
    write(100+nrank)uf

    !particle data
    write(100+nrank)up

    close(100+nrank)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,c,q,r,delt,delx,it0,nxs,nxe,                &
                        nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze,np,nsp, &
                        nproc,nproc_i,nproc_j,nproc_k,nrank,                  &
                        dir,file)
    integer, intent(in)  :: nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
    integer, intent(in)  :: np, nsp
    integer, intent(in)  :: nproc, nproc_i, nproc_j, nproc_k, nrank
    character(len=*), intent(in) :: dir, file
    integer, intent(out) :: np2(nys:nye,nzs:nze,nsp), nxs, nxe, it0
    real(8), intent(out) :: up(6,np,nys:nye,nzs:nze,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(out) :: c, q(nsp), r(nsp), delt, delx
    integer :: inp, inxgs, inxge, inygs, inyge, inzgs, inzge, inys, inye, inzs, inze
    integer :: insp, inproc, inproc_i, inproc_j, inproc_k

    !filename
    open(101+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(101+nrank)it0,inxgs,inxge,inygs,inyge,inzgs,inzge,nxs,nxe,inys,inye,inzs,inze
    read(101+nrank)inp,insp,inproc,inproc_i,inproc_j,inproc_k
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or.(inygs /= nygs) .or. (inyge /= nyge)   &
        .or. (inzgs /= nzgs) .or. (inzge /= nzge) .or. (inys /= nys) .or. (inye /= nye) &
        .or. (inzs /= nzs) .or. (inze /= nze) .or. (inp /= np) .or. (insp /= nsp)       &
        .or. (inproc /= nproc) .or. (inproc_i /= nproc_i) .or. (inproc_j /= nproc_j)    & 
        .or. (inproc_k /= nproc_k) )then
       write(6,*) '** parameter mismatch **'
       stop
    endif

    read(101+nrank)delt,delx,c
    read(101+nrank)np2
    read(101+nrank)q
    read(101+nrank)r

    !field data
    read(101+nrank)uf

    !particle data
    read(101+nrank)up

    close(101+nrank)

  end subroutine fio__input


  subroutine fio__param(nxgs,nxge,nygs,nyge,nzgs,nzge,nys,nye,nzs,nze, &
                        np,nsp,np2,                                    &
                        c,q,r,n0,temp,rtemp,fpe,fge,                   &
                        ldb,delt,delx,dir,file)

    integer, intent(in)          :: nxgs, nxge, nygs, nyge, nzgs, nzge, nys, nye, nzs, nze
    integer, intent(in)          :: np, nsp 
    integer, intent(in)          :: np2(nys:nye,nzs:nze,nsp)
    real(8), intent(in)          :: c, q(nsp), r(nsp), n0, temp, rtemp, fpe, fge, ldb, delt, delx
    character(len=*), intent(in) :: dir, file
    integer :: isp
    real(8) :: pi, vti, vte, va

    pi = 4.0*atan(1.0)

    vti = sqrt(2.*temp/r(1))
    vte = sqrt(2.*temp*rtemp/r(2))
    va  = fge*r(2)*c/q(1)/sqrt(4.*pi*r(1)*n0)

    !filename
    open(9,file=trim(dir)//trim(file),status='unknown')

    write(9,610) nxge-nxgs+1,' x ',nyge-nygs+1, ' x ',nzge-nzgs+1, ldb
    write(9,620) (np2(nys,nzs,isp),isp=1,nsp),np
    write(9,630) delx,delt,c
    write(9,640) (r(isp),isp=1,nsp)
    write(9,650) (q(isp),isp=1,nsp)
    write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
    write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
    write(9,*)
610 format(' grid size, debye lngth ============> ',i6,a,i6,a,i6,f8.4)
620 format(' particle number in cell============> ',i8,i8,'/',i8)
630 format(' dx, dt, c =========================> ',f8.4,3x,f8.4,3x,f8.4)
640 format(' Mi, Me  ===========================> ',2(1p,e10.2,1x))
650 format(' Qi, Qe  ===========================> ',2(1p,e10.2,1x))
660 format(' Fpe, Fge, Fpi Fgi =================> ',4(1p,e10.2,1x))
670 format(' Va, Vi, Ve, beta, Te/Ti, rgi     ==> ',6(1p,e10.2,1x))
    close(9)

  end subroutine fio__param


end module fio
