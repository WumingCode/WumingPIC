program jsonio_test
  use jsonio, IK => json_IK, RK => json_RK
  use iso_fortran_env, only: int64
  implicit none

  logical, parameter :: debug = .false.

  character(len=128), parameter :: jsonfile = 'jsonio_test.json'
  integer(int64), parameter :: large_int64 = (2_8)**32
  integer, parameter :: endian_flag = 1
  integer, parameter :: ndim = 2
  integer, parameter :: nx = 16
  integer, parameter :: ny = 32
  integer, parameter :: ns = 2
  integer, parameter :: shape_emf(ndim+1) = (/6_4, ny, nx/)
  integer, parameter :: shape_mom(ndim+1) = (/10_4, ny, nx/)
  real(8), parameter :: mass(ns) = (/1.0_8, 100.0_8/)
  real(8), parameter :: charge(ns) = (/-1.0_8, +1.0_8/)

  integer, parameter :: num_data_max = 32
  character(len=16)  :: json_name(num_data_max)
  character(len=16)  :: json_datatype(num_data_max)
  integer(int64)     :: json_offset(num_data_max)
  integer(int64)     :: json_size(num_data_max)
  integer(4)         :: json_ndim(num_data_max)
  integer(4)         :: json_shape(ndim+1,num_data_max)

  call write_json(jsonfile)
  call read_json(jsonfile)

contains

  ! store metadata
  subroutine store_metadata(n, name, datatype, offset, dsize, ndim, dshape)
    integer, intent(in)          :: n
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: datatype
    integer(int64), intent(in)   :: offset
    integer(int64), intent(in)   :: dsize
    integer(4), intent(in)       :: ndim
    integer(4), intent(in)       :: dshape(ndim)

    json_name(n)         = trim(name)
    json_datatype(n)     = trim(datatype)
    json_offset(n)       = offset
    json_size(n)         = dsize
    json_ndim(n)         = ndim
    json_shape(1:ndim,n) = dshape

  end subroutine store_metadata

  ! check metadata
  subroutine check_metadata(n, name, datatype, offset, dsize, ndim, dshape)
    integer, intent(in)          :: n
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: datatype
    integer(int64), intent(in)   :: offset
    integer(int64), intent(in)   :: dsize
    integer(4), intent(in)       :: ndim
    integer(4), intent(in)       :: dshape(ndim)

    integer :: i, errcnt = 0

    if( .not. trim(json_datatype(n)) == trim(datatype) ) then
       errcnt = errcnt + 1
    end if

    if( .not. json_offset(n) == offset ) then
       errcnt = errcnt + 1
    end if

    if( .not. json_size(n) == dsize ) then
       errcnt = errcnt + 1
    end if

    if( .not. json_ndim(n) == ndim ) then
       errcnt = errcnt + 1
    end if

    if( .not. all(json_shape(1:ndim,n) == dshape) ) then
       errcnt = errcnt + 1
    end if

    !
    ! report error
    !
    if( errcnt > 0 ) then
       write(0,'("*** Some errors detected for ", A)') trim(name)
       write(0,'(A10, " : ", A10, " <===> ", A10)') &
            & 'datatype', trim(json_datatype(n)), trim(datatype)
       write(0,'(A10, " : ", i10, " <===> ", i10)') &
            & 'offset', json_offset(n), offset
       write(0,'(A10, " : ", i10, " <===> ", i10)') &
            & 'size', json_size(n), dsize
       write(0,'(A10, " : ", i10, " <===> ", i10)') &
            & 'ndim', json_ndim(n), ndim
       ! shape
       write(0,'(A10, " : ")', advance='no') 'shape'
       write(0,'(A)', advance='no') '['
       do i = 1, ndim
          write(0,'(i6)', advance='no') json_shape(i,n)
       end do
       write(0,'(A)', advance='no') ' ] <===> ['
       do i = 1, ndim
          write(0,'(i6)', advance='no') dshape(i)
       end do
       write(0,'(A)', advance='no') ' ]'
       write(0,*)
    else if ( debug ) then
       write(0,'("No error detected for ", A)') trim(name)
    end if

  end subroutine check_metadata

  ! write a json file
  subroutine write_json(filename)
    implicit none
    character(len=*), intent(in) :: filename

    type(json_core) :: json
    type(json_value), pointer :: root, p

    integer(int64)     :: offset, dsize, large_int
    integer            :: numdata, endian, ndim, dshape(3)
    character(len=128) :: desc

    endian  = endian_flag
    numdata = 1

    call json%initialize()
    call json%create_object(root, 'root')

    !
    ! meta
    !
    call json%create_object(p, 'meta')
    call jsonio_check_error(json, 'write_json')

    call json%add(root, p)
    call jsonio_check_error(json, 'write_json')

    call json%add(p, 'endian', endian)
    call jsonio_check_error(json, 'write_json')

    call json%add(p, 'rawfile', 'test.raw')
    call jsonio_check_error(json, 'write_json')

    call json%add(p, 'info', 'some information')
    call jsonio_check_error(json, 'write_json')

    ! check for 64bit integer
    large_int = large_int64
    call json%add(p, 'large_int', large_int)

    nullify(p)

    !
    ! attributes
    !
    call json%create_object(p, 'attribute')
    call jsonio_check_error(json, 'write_json')

    call json%add(root, p)
    call jsonio_check_error(json, 'write_json')

    ! nx
    offset = 0
    dsize  = 4
    ndim   = 1
    dshape = (/1, 0, 0/)
    desc   = '# grid in x'
    call jsonio_put_attribute(json, p, nx, 'nx', offset, desc)
    call store_metadata(numdata, 'nx', 'i4', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    ! ny
    offset = offset + dsize
    dsize  = 4
    ndim   = 1
    dshape = (/1, 0, 0/)
    desc   = '# grid in y'
    call jsonio_put_attribute(json, p, ny, 'ny', offset, desc)
    call store_metadata(numdata, 'ny', 'i4', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    ! ns
    offset = offset + dsize
    dsize  = 4
    ndim   = 1
    dshape = (/1, 0, 0/)
    desc   = '# species'
    call jsonio_put_attribute(json, p, ns, 'ns', offset, desc)
    call store_metadata(numdata, 'ns', 'i4', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    ! mass
    offset = offset + dsize
    dsize  = 8 * 2
    ndim   = 1
    dshape = (/2, 0, 0/)
    desc   = 'mass'
    call jsonio_put_attribute(json, p, mass, 'mass', offset, desc)
    call store_metadata(numdata, 'mass', 'f8', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    ! charge
    offset = offset + dsize
    dsize  = 8 * 2
    ndim   = 1
    dshape = (/2, 0, 0/)
    desc   = 'charge'
    call jsonio_put_attribute(json, p, charge, 'charge', offset, desc)
    call store_metadata(numdata, 'charge', 'f8', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    nullify(p)

    !
    ! datasets
    !
    call json%create_object(p, 'dataset')
    call jsonio_check_error(json, 'write_json')

    call json%add(root, p)
    call jsonio_check_error(json, 'write_json')

    ! emf
    offset = offset + dsize
    dsize  = product(shape_emf) * 8
    ndim   = 3
    dshape = shape_emf
    desc   = 'electromagnetic fields'
    call jsonio_put_metadata(json, p, 'emf', 'f8', offset, &
         & dsize, ndim, dshape, desc)
    call store_metadata(numdata, 'emf', 'f8', offset, dsize, &
         & ndim, dshape)
    numdata = numdata + 1

    ! mom
    offset = offset + dsize
    dsize  = product(shape_mom) * 8
    ndim   = 3
    dshape = shape_mom
    desc   = 'moments'
    call jsonio_put_metadata(json, p, 'mom', 'f8', offset, &
         & dsize, ndim, dshape, desc)
    call store_metadata(numdata, 'mom', 'f8', offset, dsize, &
         & ndim, dshape)
    numdata = numdata + 1

    nullify(p)

    ! output
    call json%print(root, filename)

    nullify(root)
    call json%destroy()

  end subroutine write_json

  ! read and check a json file
  subroutine read_json(filename)
   implicit none
    character(len=*), intent(in) :: filename

    type(json_file) :: file
    type(json_core) :: json
    type(json_value), pointer :: root, p

    integer(int64)     :: offset, dsize, large_int
    integer            :: numdata, endian, ndim, dshape(3)
    integer            :: nx, ny, ns
    real(8)            :: mass(2), charge(2)

    numdata = 1

    call json%initialize()
    call jsonio_check_error(json, 'read_json')

    call file%initialize()
    call jsonio_check_error(file, 'read_json')

    call file%load(filename)
    call jsonio_check_error(file, 'read_json')

    ! get root
    call file%get(root)
    call jsonio_check_error(file, 'read_json')

    !
    ! meta
    !
    call json%get(root, 'meta', p)
    call json%get(p, 'endian', endian)
    call json%get(p, 'large_int', large_int)

    if( large_int /= large_int64) then
       write(0,*) 'Error in reading 64bit integer: ', large_int
    end if

    !
    ! attributes
    !
    call json%get(root, 'attribute', p)
    call jsonio_check_error(json, 'read_json')

    call jsonio_get_attribute(json, p, 'nx', offset, nx)
    call jsonio_get_metadata(json, p, 'nx', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'nx', 'i4', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    call jsonio_get_attribute(json, p, 'ny', offset, ny)
    call jsonio_get_metadata(json, p, 'ny', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'ny', 'i4', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    call jsonio_get_attribute(json, p, 'ns', offset, ns)
    call jsonio_get_metadata(json, p, 'ns', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'ns', 'i4', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    call jsonio_get_attribute(json, p, 'mass', offset, mass)
    call jsonio_get_metadata(json, p, 'mass', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'mass', 'f8', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    call jsonio_get_attribute(json, p, 'charge', offset, charge)
    call jsonio_get_metadata(json, p, 'charge', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'charge', 'f8', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    nullify(p)

    !
    ! datasets
    !
    call json%get(root, 'dataset', p)
    call jsonio_check_error(json, 'read_json')

    call jsonio_get_metadata(json, p, 'emf', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'emf', 'f8', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    call jsonio_get_metadata(json, p, 'mom', offset, dsize, ndim, dshape)
    call check_metadata(numdata, 'mom', 'f8', offset, dsize, ndim, dshape)
    numdata = numdata + 1

    nullify(root)
    nullify(p)

    call json%destroy()
    call file%destroy()

  end subroutine read_json

end program jsonio_test
