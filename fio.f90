module fio

  implicit none

  private

  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__mom


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
       write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),9999999,'_rank=',nrank,'.dat'
    else
       write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it2,'_rank=',nrank,'.dat'
    endif
    open(200+nrank,file=filename,form='unformatted')

    !time & parameters
    write(200+nrank)it2,nxgs,nxge,nygs,nyge,nzgs,nzge,nxs,nxe,nys,nye,nzs,nze
    write(200+nrank)np,nsp,nproc,nproc_i,nproc_j,nproc_k
    write(200+nrank)delt,delx,c
    write(200+nrank)np2
    write(200+nrank)q
    write(200+nrank)r

    !field data
    write(200+nrank)uf

    !particle data
    write(200+nrank)up

    close(200+nrank)

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
    open(201+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(201+nrank)it0,inxgs,inxge,inygs,inyge,inzgs,inzge,nxs,nxe,inys,inye,inzs,inze
    read(201+nrank)inp,insp,inproc,inproc_i,inproc_j,inproc_k
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or.(inygs /= nygs) .or. (inyge /= nyge)   &
        .or. (inzgs /= nzgs) .or. (inzge /= nzge) .or. (inys /= nys) .or. (inye /= nye) &
        .or. (inzs /= nzs) .or. (inze /= nze) .or. (inp /= np) .or. (insp /= nsp)       &
        .or. (inproc /= nproc) .or. (inproc_i /= nproc_i) .or. (inproc_j /= nproc_j)    & 
        .or. (inproc_k /= nproc_k) )then
       write(6,*) '** parameter mismatch **'
       stop
    endif

    read(201+nrank)delt,delx,c
    read(201+nrank)np2
    read(201+nrank)q
    read(201+nrank)r

    !field data
    read(201+nrank)uf

    !particle data
    read(201+nrank)up

    close(201+nrank)

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


  subroutine fio__mom(den,vel,temp,uf,nxgs,nxge,nys,nye,nzs,nze,nsp,it0,irank,dir)

    integer, intent(in)    :: nxgs, nxge, nys, nye, nzs, nze, nsp, it0, irank
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2,nzs-2:nze+2)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,nsp),    &
                              vel(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,3,nsp),  &
                              temp(nxgs-1:nxge+1,nys-1:nye+1,nzs-1:nze+1,3,nsp)
    character(len=*), intent(in) :: dir
    integer            :: i, j, k
    real(8)            :: tmp(nxgs:nxge,nys-1:nye+1,nzs-1:nze+1,6)
    character(len=256) :: filename

    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_den_i_rank=',irank,'.dat'
    open(10,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_den_e_rank=',irank,'.dat'
    open(11,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_txx_i_rank=',irank,'.dat'
    open(12,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_tyy_i_rank=',irank,'.dat'
    open(13,file=filename,status='unknown',form='unformatted')	
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_tzz_i_rank=',irank,'.dat'
    open(14,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_txx_e_rank=',irank,'.dat'
    open(15,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_tyy_e_rank=',irank,'.dat'
    open(16,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_tzz_e_rank=',irank,'.dat'
    open(17,file=filename,status='unknown',form='unformatted')

    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_vx_i_rank=',irank,'.dat'
    open(18,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_vy_i_rank=',irank,'.dat'
    open(19,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_vz_i_rank=',irank,'.dat'
    open(20,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_vx_e_rank=',irank,'.dat'
    open(21,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_vy_e_rank=',irank,'.dat'
    open(22,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_vz_e_rank=',irank,'.dat'
    open(23,file=filename,status='unknown',form='unformatted')

    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_bx_rank=',irank,'.dat'
    open(24,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_by_rank=',irank,'.dat'
    open(25,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_bz_rank=',irank,'.dat'
    open(26,file=filename,status='unknown',form='unformatted')

    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_ex_rank=',irank,'.dat'
    open(27,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_ey_rank=',irank,'.dat'
    open(28,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i5.5,a)')trim(dir),it0,'_ez_rank=',irank,'.dat'
    open(29,file=filename,status='unknown',form='unformatted')

    !fields at (i+1/2, j+1/2, k+1/2)
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=nzs,nze
    do j=nys,nye
    do i=nxgs,nxge-1
       tmp(i,j,k,1) = 0.25*(+uf(1,i,j,k)  +uf(1,i,j+1,k) &
                            +uf(1,i,j,k+1)+uf(1,i,j+1,k+1))
       tmp(i,j,k,2) = 0.25*(+uf(2,i,j,k)  +uf(2,i+1,j,k) &
                            +uf(2,i,j,k+1)+uf(2,i+1,j,k+1))
       tmp(i,j,k,3) = 0.25*(+uf(3,i,j,k)  +uf(3,i+1,j,k) &
                            +uf(3,i,j+1,k)+uf(3,i+1,j+1,k))
       tmp(i,j,k,4) = 0.5*(+uf(4,i,j,k)+uf(4,i+1,j,k))
       tmp(i,j,k,5) = 0.5*(+uf(5,i,j,k)+uf(5,i,j+1,k))
       tmp(i,j,k,6) = 0.5*(+uf(6,i,j,k)+uf(6,i,j,k+1))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    write(10)sngl(den(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,1))
    write(11)sngl(den(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,2))
    write(12)sngl(temp(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,1,1))
    write(13)sngl(temp(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,2,1))
    write(14)sngl(temp(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,3,1))
    write(15)sngl(temp(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,1,2))
    write(16)sngl(temp(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,2,2))
    write(17)sngl(temp(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,3,2))
    write(18)sngl(vel(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,1,1))
    write(19)sngl(vel(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,2,1))
    write(20)sngl(vel(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,3,1))
    write(21)sngl(vel(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,1,2))
    write(22)sngl(vel(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,2,2))
    write(23)sngl(vel(nxgs:nxge-1,nys-1:nye+1,nzs-1:nze+1,3,2))
    write(24)sngl(tmp(nxgs:nxge-1,nys:nye,nzs:nze,1))
    write(25)sngl(tmp(nxgs:nxge-1,nys:nye,nzs:nze,2))
    write(26)sngl(tmp(nxgs:nxge-1,nys:nye,nzs:nze,3))
    write(27)sngl(tmp(nxgs:nxge-1,nys:nye,nzs:nze,4))
    write(28)sngl(tmp(nxgs:nxge-1,nys:nye,nzs:nze,5))
    write(29)sngl(tmp(nxgs:nxge-1,nys:nye,nzs:nze,6))
    
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)
    close(27)
    close(28)
    close(29)

  end subroutine fio__mom


end module fio
