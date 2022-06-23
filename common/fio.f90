module fio

  implicit none

  private

  public :: fio__init
  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__mom
  public :: fio__orb

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nproc, nrank
  real(8), save :: delx, delt, u0, c
  real(8), allocatable :: q(:), r(:)
  character(len=256), save :: dir


contains


  subroutine fio__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
                       nproc_in,nrank_in,                                                  &
                       delx_in,delt_in,c_in,q_in,r_in,dir_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nproc_in, nrank_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)
    character(len=*), intent(in) :: dir_in

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in
    nproc = nproc_in
    nrank = nrank_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in
    dir   = dir_in

    is_init = .true.

  end subroutine fio__init


  subroutine fio__output(up,uf,np2,nxs,nxe,it0,lflag)

    logical, intent(in) :: lflag
    integer, intent(in) :: np2(nys:nye,nsp), nxs, nxe
    integer, intent(in) :: it0
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    character(len=256) :: filename

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    if(lflag)then
       write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),9999999,'_rank=',nrank,'.dat'
    else
       write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_rank=',nrank,'.dat'
    endif
    open(200+nrank,file=filename,form='unformatted')

    !time & parameters
    write(200+nrank)it0,ndim,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,nproc,delt,delx,c
    write(200+nrank)np2
    write(200+nrank)q
    write(200+nrank)r

    !field data
    write(200+nrank)uf

    !particle data
    write(200+nrank)up

    close(200+nrank)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,indim,nxs,nxe,it0,file)

    character(len=*), intent(in) :: file
    integer, intent(out) :: np2(nys:nye,nsp), nxs, nxe, it0, indim
    real(8), intent(out) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: inp, inxgs, inxge, inygs, inyge, inys, inye, insp, inproc

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    open(201+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(201+nrank)it0,indim,inp,inxgs,inxge,inygs,inyge,nxs,nxe,inys,inye,insp,inproc,delt,delx,c
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or.(inygs /= nygs) .or. (inyge /= nyge) &
        .or. (inys /= nys) .or. (inye /= nye) .or. (inp /= np) .or. (insp /= nsp) &
        .or. (inproc /= nproc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif

    read(201+nrank)np2
    read(201+nrank)q
    read(201+nrank)r

    !field data
    read(201+nrank)uf

    !particle data
    read(201+nrank)up(1:indim,:,:,:)

    close(201+nrank)

  end subroutine fio__input


  subroutine fio__param(n0,np2,temp,rtemp,fpe,fge,ls,file,nroot)

    integer, intent(in)          :: n0, nroot
    integer, intent(in)          :: np2(nys:nye,nsp)
    real(8), intent(in)          :: temp, rtemp, fpe, fge, ls
    character(len=*), intent(in) :: file
    integer :: isp
    real(8) :: pi, vti, vte, va
    character(len=256) :: filename

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    pi = 4.0D0*datan(1.0D0)

    vti  = sqrt(2.*temp/r(1))
    vte  = sqrt(2.*temp*rtemp/r(2))
    va   = fge*r(2)*c/q(1)/sqrt(4.*pi*r(1)*n0)

    if(nrank == nroot)then

       !filename
       filename = trim(dir)//trim(file)//".dat"
       open(9,file=filename,status='unknown')

       write(9,610) nxge-nxgs+1,' x ',nyge-nygs+1, ls
       write(9,620) (np2(nys,isp),isp=1,nsp),np
       write(9,630) delx,delt,c
       write(9,640) (r(isp),isp=1,nsp)
       write(9,650) (q(isp),isp=1,nsp)
       write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
       write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
       write(9,*)
610    format(' grid size, electron skin depth ====>',i6,a,i6,f8.4)
620    format(' particle number in cell============> ',i8,i8,'/',i8)
630    format(' dx, dt, c =========================> ',f8.4,3x,f8.4,3x,f8.4)
640    format(' Mi, Me  ===========================> ',2(1p,e10.2,1x))
650    format(' Qi, Qe  ===========================> ',2(1p,e10.2,1x))
660    format(' Fpe, Fge, Fpi Fgi =================> ',4(1p,e10.2,1x))
670    format(' Va, Vi, Ve, beta, Te/Ti, rgi     ==> ',6(1p,e10.2,1x))
       close(9)

    endif

  end subroutine fio__param


  subroutine fio__mom(den,vel,temp,uf,it0)

    integer, intent(in)    :: it0
    real(8), intent(in)    :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nys-1:nye+1,nsp),    &
                              vel(nxgs-1:nxge+1,nys-1:nye+1,3,nsp),  &
                              temp(nxgs-1:nxge+1,nys-1:nye+1,3,nsp)
    integer :: i, j
    real(8) :: tmp(nxgs:nxge,nys:nye,1:6)
    character(len=256) :: filename


    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_den_i_rank=',nrank,'.dat'
    open(10,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_den_e_rank=',nrank,'.dat'
    open(11,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_txx_i_rank=',nrank,'.dat'
    open(12,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_tyy_i_rank=',nrank,'.dat'
    open(13,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_tzz_i_rank=',nrank,'.dat'
    open(14,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_txx_e_rank=',nrank,'.dat'
    open(15,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_tyy_e_rank=',nrank,'.dat'
    open(16,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_tzz_e_rank=',nrank,'.dat'
    open(17,file=filename,status='unknown',form='unformatted')

    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_vx_i_rank=',nrank,'.dat'
    open(18,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_vy_i_rank=',nrank,'.dat'
    open(19,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_vz_i_rank=',nrank,'.dat'
    open(20,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_vx_e_rank=',nrank,'.dat'
    open(21,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_vy_e_rank=',nrank,'.dat'
    open(22,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_vz_e_rank=',nrank,'.dat'
    open(23,file=filename,status='unknown',form='unformatted')

    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_bx_rank=',nrank,'.dat'
    open(24,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_by_rank=',nrank,'.dat'
    open(25,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_bz_rank=',nrank,'.dat'
    open(26,file=filename,status='unknown',form='unformatted')

    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_ex_rank=',nrank,'.dat'
    open(27,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_ey_rank=',nrank,'.dat'
    open(28,file=filename,status='unknown',form='unformatted')
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_ez_rank=',nrank,'.dat'
    open(29,file=filename,status='unknown',form='unformatted')


    do j=nys,nye
    do i=nxgs,nxge
       tmp(i,j,1) = 0.5*(uf(1,i,j)+uf(1,i,j+1))
       tmp(i,j,2) = 0.5*(uf(2,i,j)+uf(2,i+1,j))
       tmp(i,j,3) = 0.25*( uf(3,i,j)+uf(3,i+1,j)      &
                          +uf(3,i+1,j)+uf(3,i+1,j+1))
       tmp(i,j,4) = 0.5*(uf(4,i,j)+uf(4,i+1,j))
       tmp(i,j,5) = 0.5*(uf(5,i,j)+uf(5,i,j+1))
       tmp(i,j,6) = uf(6,i,j)
    enddo
    enddo

    write(10)sngl(den(nxgs:nxge,nys-1:nye+1,1))
    write(11)sngl(den(nxgs:nxge,nys-1:nye+1,2))
    write(12)sngl(temp(nxgs:nxge,nys-1:nye+1,1,1))
    write(13)sngl(temp(nxgs:nxge,nys-1:nye+1,2,1))
    write(14)sngl(temp(nxgs:nxge,nys-1:nye+1,3,1))
    write(15)sngl(temp(nxgs:nxge,nys-1:nye+1,1,2))
    write(16)sngl(temp(nxgs:nxge,nys-1:nye+1,2,2))
    write(17)sngl(temp(nxgs:nxge,nys-1:nye+1,3,2))
    write(18)sngl(vel(nxgs:nxge,nys-1:nye+1,1,1))
    write(19)sngl(vel(nxgs:nxge,nys-1:nye+1,2,1))
    write(20)sngl(vel(nxgs:nxge,nys-1:nye+1,3,1))
    write(21)sngl(vel(nxgs:nxge,nys-1:nye+1,1,2))
    write(22)sngl(vel(nxgs:nxge,nys-1:nye+1,2,2))
    write(23)sngl(vel(nxgs:nxge,nys-1:nye+1,3,2))
    write(24)sngl(tmp(nxgs:nxge,nys:nye,1))
    write(25)sngl(tmp(nxgs:nxge,nys:nye,2))
    write(26)sngl(tmp(nxgs:nxge,nys:nye,3))
    write(27)sngl(tmp(nxgs:nxge,nys:nye,4))
    write(28)sngl(tmp(nxgs:nxge,nys:nye,5))
    write(29)sngl(tmp(nxgs:nxge,nys:nye,6))

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


  subroutine fio__orb(up,uf,np2,it0)

    integer, intent(in) :: it0
    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: up(ndim,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer :: isp, ii, j
    character(len=256) :: filename

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_orb_i_rank=',nrank,'.dat'
    open(1000+nrank,file=filename)
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_orb_e_rank=',nrank,'.dat'
    open(2000+nrank,file=filename)
    write(filename,'(a,i7.7,a,i3.3,a)')trim(dir),it0,'_eb_rank=',nrank,'.dat'
    open(3000+nrank,file=filename,form='unformatted')

    do isp=1,nsp

!$OMP PARALLEL DO PRIVATE(ii,j)
       do j=nys,nye
          do ii=1,np2(j,isp)
             if(up(6,ii,j,isp) > 0.D0)then

                !write particle data
                if(isp==1)then
                   write(1000+nrank,'(5e17.7,f22.12)')up(1:5,ii,j,isp),up(6,ii,j,isp)
                else
                   write(2000+nrank,'(5e17.7,f22.12)')up(1:5,ii,j,isp),up(6,ii,j,isp)
                endif
             endif
          enddo
       enddo
!$OMP END PARALLEL DO

    enddo

    write(3000+nrank)sngl(uf(1:6,nxgs:nxge,nys:nye))

    close(1000+nrank)
    close(2000+nrank)
    close(3000+nrank)

  end subroutine fio__orb


end module fio
