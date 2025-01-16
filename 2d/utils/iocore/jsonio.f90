!
! JSON IO routines
!
! written by Takanobu Amano <amano@eps.s.u-tokyo.ac.jp>
!
module jsonio
  use json_module
  use iso_fortran_env, only: int64
  implicit none
  private

  public :: json_IK
  public :: json_RK
  public :: json_core
  public :: json_value
  public :: json_file
  public :: jsonio_check_error
  public :: jsonio_put_metadata
  public :: jsonio_put_attribute
  public :: jsonio_get_metadata
  public :: jsonio_get_attribute

  integer, parameter :: IK = json_IK
  integer, parameter :: RK = json_RK

  !
  ! check error
  !
  interface jsonio_check_error
     module procedure &
          & check_error_core, &
          & check_error_file
  end interface jsonio_check_error

  !
  ! put metadata
  !
  interface jsonio_put_metadata
     module procedure &
          & put_metadata_4, &
          & put_metadata_8
  end interface jsonio_put_metadata

  !
  ! put attribute; get metadata and data
  !
  interface jsonio_put_attribute
     module procedure &
          & put_scalar_i4, &
          & put_array_i4,  &
          & put_scalar_i8, &
          & put_array_i8,  &
          & put_scalar_r4, &
          & put_array_r4,  &
          & put_scalar_r8, &
          & put_array_r8
  end interface jsonio_put_attribute

  !
  ! get metadata
  !
  interface jsonio_get_metadata
     module procedure &
          & get_metadata_4, &
          & get_metadata_8
  end interface jsonio_get_metadata

  !
  ! get attribute; get metadata and data
  !
  interface jsonio_get_attribute
     module procedure &
          & get_scalar_i4, &
          & get_array_i4,  &
          & get_scalar_r4, &
          & get_array_r4,  &
          & get_scalar_r8, &
          & get_array_r8
  end interface jsonio_get_attribute

contains

  !
  ! put metadata for integer(4) ndim and dhspae
  !
  subroutine put_metadata_4(json, dst, name, datatype, disp, &
       & dsize, ndim, dshape, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    character(len=*), intent(in)   :: name
    character(len=*), intent(in)   :: datatype
    integer(int64), intent(in)     :: disp
    integer(int64), intent(in)     :: dsize
    integer(4), intent(in)         :: ndim
    integer(4), intent(in)         :: dshape(ndim)
    character(len=*), intent(in)   :: desc

    call put_metadata_ik(json, dst, name, datatype, disp, dsize, &
         & int(ndim, ik), int(dshape, ik), desc)

  end subroutine put_metadata_4

  !
  ! put metadata for integer(8) ndim and dshape
  !
  subroutine put_metadata_8(json, dst, name, datatype, disp, &
       & dsize, ndim, dshape, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    character(len=*), intent(in)   :: name
    character(len=*), intent(in)   :: datatype
    integer(int64), intent(in)     :: disp
    integer(int64), intent(in)     :: dsize
    integer(8), intent(in)         :: ndim
    integer(8), intent(in)         :: dshape(ndim)
    character(len=*), intent(in)   :: desc

    call put_metadata_ik(json, dst, name, datatype, disp, dsize, &
         & int(ndim, ik), int(dshape, ik), desc)

  end subroutine put_metadata_8

  !
  ! put metadata
  !
  subroutine put_metadata_ik(json, dst, name, datatype, disp, &
       & dsize, ndim, dshape, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    character(len=*), intent(in)   :: name
    character(len=*), intent(in)   :: datatype
    integer(int64), intent(in)     :: disp
    integer(int64), intent(in)     :: dsize
    integer(IK), intent(in)        :: ndim
    integer(IK), intent(in)        :: dshape(ndim)
    character(len=*), intent(in)   :: desc

    type(json_value), pointer :: obj

    call json%create_object(obj, name)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'datatype', datatype)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'offset', disp)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'size', dsize)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'ndim', ndim)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'shape', dshape)
    call jsonio_check_error(json, 'put_metadata')

    call json%add(obj, 'description', trim(desc))
    call jsonio_check_error(json, 'put_metadata')

    ! add to dst
    call json%add(dst, obj)
    call jsonio_check_error(json, 'put_metadata')

    nullify(obj)

  end subroutine put_metadata_ik

  ! put scalar of integer(4)
  subroutine put_scalar_i4(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(4), intent(in)         :: data
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    integer(IK) :: x

    x = data
    s = 4
    call jsonio_put_metadata(json, dst, name, 'i4', disp, s, &
         & 1, (/1/), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_i4

  ! put array of integer(4)
  subroutine put_array_i4(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(4), intent(in)         :: data(:)
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    integer(IK) :: x(size(data))

    x = data
    s = size(data) * 4
    call jsonio_put_metadata(json, dst, name, 'i4', disp, s, &
         & 1, shape(data), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_i4

  ! put scalar of integer(8)
  subroutine put_scalar_i8(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(8), intent(in)         :: data
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    integer(IK) :: x

    x = data
    s = 8
    call jsonio_put_metadata(json, dst, name, 'i8', disp, s, &
         & 1, (/1/), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_i8

  ! put array of integer(8)
  subroutine put_array_i8(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    integer(8), intent(in)         :: data(:)
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    integer(IK) :: x(size(data))

    x = data
    s = size(data) * 8
    call jsonio_put_metadata(json, dst, name, 'i8', disp, s, &
         & 1, shape(data), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_i8

  ! put scalar of real(4)
  subroutine put_scalar_r4(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(4), intent(in)            :: data
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    real(RK) :: x

    x = data
    s = 4
    call jsonio_put_metadata(json, dst, name, 'f4', disp, s, &
         & 1, (/1/), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_r4

  ! put array of real(4)
  subroutine put_array_r4(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(4), intent(in)            :: data(:)
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    real(RK) :: x(size(data))

    x = data
    s = size(data) * 4
    call jsonio_put_metadata(json, dst, name, 'f4', disp, s, &
         & 1, shape(data), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_r4

  ! put scalar of real(8)
  subroutine put_scalar_r8(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(8), intent(in)            :: data
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    real(RK) :: x

    x = data
    s = 8
    call jsonio_put_metadata(json, dst, name, 'f8', disp, s, &
         & 1, (/1/), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_scalar_r8

  ! put array of real(8)
  subroutine put_array_r8(json, dst, data, name, disp, desc)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: dst
    real(8), intent(in)            :: data(:)
    character(len=*), intent(in)   :: name
    integer(int64), intent(in)     :: disp
    character(len=*), intent(in)   :: desc

    integer(int64) :: s
    real(RK) :: x(size(data))

    x = data
    s = size(data) * 8
    call jsonio_put_metadata(json, dst, name, 'f8', disp, s, &
         & 1, shape(data), desc)
    call json%add_by_path(dst, name // '.data', x)
    call jsonio_check_error(json, 'put_attribute')

  end subroutine put_array_r8

  !
  ! get metadata for integer(4) ndim and dshape
  !
  subroutine get_metadata_4(json, src, name, disp, dsize, ndim, dshape)
    implicit none
    type(json_core), intent(inout)    :: json
    type(json_value), pointer         :: src
    character(len=*), intent(in)      :: name
    integer(int64), intent(out)       :: disp
    integer(int64), intent(out)       :: dsize
    integer(4), intent(out)           :: ndim
    integer(4), intent(out)           :: dshape(:)

    integer(IK) :: ndim_ik, dshape_ik(size(dshape))

    call get_metadata_ik(json, src, name, disp, dsize, ndim_ik, dshape_ik)
    ndim   = ndim_ik
    dshape = dshape_ik

  end subroutine get_metadata_4

  !
  ! get metadata for integer(8) ndim and dshape
  !
  subroutine get_metadata_8(json, src, name, disp, dsize, ndim, dshape)
    implicit none
    type(json_core), intent(inout)    :: json
    type(json_value), pointer         :: src
    character(len=*), intent(in)      :: name
    integer(int64), intent(out)       :: disp
    integer(int64), intent(out)       :: dsize
    integer(8), intent(out)           :: ndim
    integer(8), intent(out)           :: dshape(:)

    integer(IK) :: ndim_ik, dshape_ik(size(dshape))

    call get_metadata_ik(json, src, name, disp, dsize, ndim_ik, dshape_ik)
    ndim   = ndim_ik
    dshape = dshape_ik

  end subroutine get_metadata_8

  !
  ! get metadata
  !
  subroutine get_metadata_ik(json, src, name, disp, dsize, ndim, dshape)
    implicit none
    type(json_core), intent(inout)    :: json
    type(json_value), pointer         :: src
    character(len=*), intent(in)      :: name
    integer(int64), intent(out)       :: disp
    integer(int64), intent(out)       :: dsize
    integer(IK), intent(out)          :: ndim
    integer(IK), intent(out)          :: dshape(:)

    character(len=128) :: dataname
    integer :: i

    dshape = 0

    call json%get(src, name // '.offset', disp)
    call jsonio_check_error(json, 'get_metadata')

    call json%get(src, name // '.size', dsize)
    call jsonio_check_error(json, 'get_metadata')

    call json%get(src, name // '.ndim', ndim)
    call jsonio_check_error(json, 'get_metadata')

    do i = 1, ndim
       write(dataname,'(A,".shape(",i2,")")') name, i
       call json%get(src, trim(dataname), dshape(i))
       call jsonio_check_error(json, 'get_metadata')
    end do

  end subroutine get_metadata_ik

  ! get scalar of integer(4)
  subroutine get_scalar_i4(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    integer(4), intent(out)        :: data

    integer(IK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_i4

  ! get array of integer(4)
  subroutine get_array_i4(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    integer(4), intent(out)        :: data(:)

    integer(int64) :: dsize
    character(len=128) :: dataname
    integer(IK) :: i, x, ndim, dshape(1)

    call jsonio_get_metadata(json, src, name, disp, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_i4

  ! get scalar of integer(8)
  subroutine get_scalar_i8(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    integer(8), intent(out)        :: data

    integer(IK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_i8

  ! get array of integer(8)
  subroutine get_array_i8(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    integer(8), intent(out)        :: data(:)

    integer(int64) :: dsize
    character(len=128) :: dataname
    integer(IK) :: i, x, ndim, dshape(1)

    call jsonio_get_metadata(json, src, name, disp, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_i8

  ! get scalar of real(4)
  subroutine get_scalar_r4(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    real(4), intent(out)           :: data

    real(RK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_r4

  ! get array of real(4)
  subroutine get_array_r4(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    real(4), intent(out)           :: data(:)

    integer(int64) :: dsize
    character(len=128) :: dataname
    integer(IK) :: i, ndim, dshape(1)
    real(RK) :: x

    call jsonio_get_metadata(json, src, name, disp, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_r4

  ! get scalar of real(8)
  subroutine get_scalar_r8(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    real(8), intent(out)           :: data

    real(RK) :: x

    call json%get(src, name // '.data', x)
    call jsonio_check_error(json, 'get_attribute')
    data = x

  end subroutine get_scalar_r8

  ! get array of real(8)
  subroutine get_array_r8(json, src, name, disp, data)
    implicit none
    type(json_core), intent(inout) :: json
    type(json_value), pointer      :: src
    character(len=*), intent(in)   :: name
    integer(int64), intent(inout)  :: disp
    real(8), intent(out)           :: data(:)

    integer(int64) :: dsize
    character(len=128) :: dataname
    integer(IK) :: i, ndim, dshape(1)
    real(RK) :: x

    call jsonio_get_metadata(json, src, name, disp, dsize, ndim, dshape)

    do i = 1, size(data)
       write(dataname,'(A,".data(",i2,")")') name, i
       call json%get(src, trim(dataname), x)
       call jsonio_check_error(json, 'get_attribute')
       data(i) = x
    end do

  end subroutine get_array_r8

  ! error check code
  subroutine check_error_core(json, message)
    implicit none
    type(json_core), intent(inout) :: json
    character(len=*), intent(in)   :: message

    if( json%failed() ) then
       write(0, *) 'Some error detected during json operation: ', message
       call json%print_error_message(0)
    end if

  end subroutine check_error_core

  ! error check code
  subroutine check_error_file(json, message)
    implicit none
    type(json_file), intent(inout) :: json
    character(len=*), intent(in)   :: message

    if( json%failed() ) then
       write(0, *) 'Some error detected during json operation: ', message
       call json%print_error_message(0)
    end if

  end subroutine check_error_file

end module jsonio
