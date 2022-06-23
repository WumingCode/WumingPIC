#ifdef USE_HDF5
!
! Parallel HDF5 ultility routines
!
! written by Takanobu Amano <amano@eps.s.u-tokyo.ac.jp>
!
module h5util
  use hdf5
  use mpi
  implicit none
  private

  integer, parameter :: MAXDIM = 32
  integer :: hdferr

  ! initialization of hdf5 library
  interface h5util_init
     module procedure init
  end interface h5util_init

  ! finalization of hdf5 library
  interface h5util_finalize
     module procedure finalize
  end interface h5util_finalize

  ! create hdf5 file
  interface h5util_create_file
     module procedure create_file
  end interface h5util_create_file

  ! open existing file
  interface h5util_open_file
     module procedure open_file
  end interface h5util_open_file

  ! close file
  interface h5util_close_file
     module procedure close_file
  end interface h5util_close_file

  ! put attribute
  interface h5util_put_attribute
     module procedure &
          & put_attribute_char,&
          & put_attribute_i4, &
          & put_attribute_i8, &
          & put_attribute_r4, &
          & put_attribute_r8, &
          & put_attribute_arr_i4, &
          & put_attribute_arr_i8, &
          & put_attribute_arr_r4, &
          & put_attribute_arr_r8
  end interface h5util_put_attribute

  ! get attribute
  interface h5util_get_attribute
     module procedure &
          & get_attribute_char,   &
          & get_attribute_int,    &
          & get_attribute_r4,     &
          & get_attribute_r8,     &
          & get_attribute_arr_int,&
          & get_attribute_arr_r4, &
          & get_attribute_arr_r8
  end interface h5util_get_attribute

  ! write global dataset
  interface h5util_write_dataset
     module procedure &
          & write_dataset_i4, &
          & write_dataset_i8, &
          & write_dataset_r4, &
          & write_dataset_r8
  end interface h5util_write_dataset

  ! read global dataset
  interface h5util_read_dataset
     module procedure &
          & read_dataset_int,&
          & read_dataset_r4, &
          & read_dataset_r8
  end interface h5util_read_dataset

  ! create global dataset
  interface h5util_create_dataset
     module procedure create_dataset
  end interface h5util_create_dataset

  ! create global dataset with explicitly specified chunk
  interface h5util_create_dataset_chunk
     module procedure create_dataset_chunk
  end interface h5util_create_dataset_chunk

  ! extend last (i.e., time) dimension
  interface h5util_extend_dimension
     module procedure extend_dimension
  end interface h5util_extend_dimension

  ! get array dimension
  interface h5util_get_dimension
     module procedure get_dimension
  end interface h5util_get_dimension

  ! public routiens
  public :: h5util_init
  public :: h5util_finalize
  public :: h5util_create_file
  public :: h5util_open_file
  public :: h5util_close_file
  public :: h5util_put_attribute
  public :: h5util_get_attribute
  public :: h5util_write_dataset
  public :: h5util_read_dataset
  public :: h5util_create_dataset
  public :: h5util_create_dataset_chunk
  public :: h5util_extend_dimension
  public :: h5util_get_dimension

contains
  !
  ! initialize library
  !
  subroutine init()
    implicit none

    call h5open_f(hdferr)

  end subroutine init

  !
  ! finalize library
  !
  subroutine finalize()
    implicit none

    call h5close_f(hdferr)

  end subroutine finalize

  !
  ! create file (replace if already exists)
  !
  subroutine create_file(filename)
    implicit none
    character(len=*), intent(in) :: filename

    integer(hid_t) :: prop, file

    call h5pcreate_f(H5P_FILE_ACCESS_F, prop, hdferr)
    call h5pset_fapl_mpio_f(prop, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr, &
         & H5P_DEFAULT_F, prop)
    call h5pclose_f(prop, hdferr)
    call h5fclose_f(file, hdferr)

  end subroutine create_file

  !
  ! open file
  !
  subroutine open_file(filename, file)
    implicit none
    character(len=*), intent(in) :: filename
    integer(hid_t), intent(out)  :: file

    integer(hid_t) :: prop

    call h5pcreate_f(H5P_FILE_ACCESS_F, prop, hdferr)
    call h5pset_fapl_mpio_f(prop, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file, hdferr, prop)
    call h5pclose_f(prop, hdferr)

  end subroutine open_file

  !
  ! close file
  !
  subroutine close_file(file)
    implicit none
    integer(hid_t), intent(in)  :: file

    call h5fclose_f(file, hdferr)

  end subroutine close_file

  !
  ! put scalar character attribute
  !
  subroutine put_attribute_char(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)
    integer(hsize_t) :: length

    length = len(data)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, type, hdferr)
    call h5tset_size_f(type, length, hdferr)

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5tclose_f(type, hdferr)
    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_char

  !
  ! put scalar integer(4) attribute
  !
  subroutine put_attribute_i4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(4), intent(in)       :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_INTEGER

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_i4

  !
  ! put scalar integer(8) attribute
  !
  subroutine put_attribute_i8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(8), intent(in)       :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_INTEGER

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_i8

  !
  ! put scalar real(4) attribute
  !
  subroutine put_attribute_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(in)          :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_REAL

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_r4

  !
  ! put scalar real(8) attribute
  !
  subroutine put_attribute_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(in)          :: data

    integer(hid_t) :: scalar, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_DOUBLE

    call h5screate_f(H5S_SCALAR_F, scalar, hdferr)
    call h5acreate_f(dest, name, type, scalar, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(scalar, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_r8

  !
  ! put array integer(4) attribute
  !
  subroutine put_attribute_arr_i4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(4), intent(in)       :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_INTEGER
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_i4

  !
  ! put array integer(8) attribute
  !
  subroutine put_attribute_arr_i8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(8), intent(in)       :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_INTEGER
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_i8

  !
  ! put array integer attribute
  !
  subroutine put_attribute_arr_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(in)          :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_REAL
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_r4

  !
  ! put array integer attribute
  !
  subroutine put_attribute_arr_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(in)          :: data(:)

    integer(hid_t) :: array, attr, type
    integer(hsize_t) :: dims(1)

    type = H5T_NATIVE_DOUBLE
    dims(1) = size(data)

    call h5screate_simple_f(1, dims, array, hdferr)
    call h5acreate_f(dest, name, type, array, attr, hdferr, &
         & H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5awrite_f(attr, type, data, dims, hdferr)

    call h5sclose_f(array, hdferr)
    call h5aclose_f(attr, hdferr)

  end subroutine put_attribute_arr_r8

  !
  ! create dataset
  !
  subroutine create_dataset(dest, name, type, rank, dims, maxd)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: type
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in), optional :: maxd(rank)

    integer(hsize_t) :: max_dims(rank)

    if( present(maxd) ) then
       max_dims = maxd
    else
       max_dims = dims
    end if

    ! by default, assume chunk == dims
    call create_dataset_chunk(dest, name, type, rank, dims, dims, max_dims)

  end subroutine create_dataset

  !
  ! create dataset with explicitly specified chunk
  !
  subroutine create_dataset_chunk(dest, name, type, rank, dims, cdim, maxd)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: type
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: cdim(rank)
    integer(hsize_t), intent(in), optional :: maxd(rank)

    integer(hid_t) :: space, plist, dtype, dataset

    if( present(maxd) ) then
       call h5screate_simple_f(rank, dims, space, hdferr, maxd)
    else
       call h5screate_simple_f(rank, dims, space, hdferr)
    end if

    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, hdferr)
    call h5tcopy_f(type, dtype, hdferr)

    call h5pset_chunk_f(plist, rank, cdim, hdferr)
    call h5tset_order_f(dtype, H5T_ORDER_LE_F, hdferr)

    call h5dcreate_f(dest, name, dtype, space, dataset, hdferr, &
         & plist, H5P_DEFAULT_F, H5P_DEFAULT_F)

    call h5sclose_f(space, hdferr)
    call h5pclose_f(plist, hdferr)
    call h5tclose_f(dtype, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine create_dataset_chunk

  !
  ! write integer global array
  !
  subroutine write_dataset_i4(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    integer(4), intent(in)       :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_i4

  !
  ! write integer global array
  !
  subroutine write_dataset_i8(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    integer(8), intent(in)       :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_i8

  !
  ! write real(4) global array
  !
  subroutine write_dataset_r4(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(4), intent(in)          :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_REAL, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_r4

  !
  ! write real(8) global array
  !
  subroutine write_dataset_r8(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(8), intent(in)          :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! write
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, data, dims, hdferr, &
         & mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine write_dataset_r8

  !
  ! read scalar character attribute
  !
  subroutine get_attribute_char(dest, name, data)
    implicit none
    integer(hid_t), intent(in)    :: dest
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: data

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_CHARACTER, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_char

  !
  ! read scalar integer attribute
  !
  subroutine get_attribute_int(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(out)         :: data

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_INTEGER, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_int

  !
  ! read scalar real(4) attribute
  !
  subroutine get_attribute_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(out)         :: data

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_REAL, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_r4

  !
  ! read scalar real(8) attribute
  !
  subroutine get_attribute_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(out)         :: data

    integer(hid_t) :: attr, type
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_DOUBLE, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_r8

  !
  ! read array integer attribute
  !
  subroutine get_attribute_arr_int(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(out)         :: data(:)

    integer(hid_t) :: attr, type
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_INTEGER, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_arr_int

  !
  ! read array real(4) attribute
  !
  subroutine get_attribute_arr_r4(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(4), intent(out)          :: data(:)

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_REAL, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_arr_r4

  !
  ! read array real(8) attribute
  !
  subroutine get_attribute_arr_r8(dest, name, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    real(8), intent(out)         :: data(:)

    integer(hid_t) :: attr
    integer(hsize_t) :: dims(1)

    call h5aopen_f(dest, name, attr, hdferr, H5P_DEFAULT_F)

    call h5aread_f(attr, H5T_NATIVE_DOUBLE, data, dims, hdferr)

    call h5aclose_f(attr, hdferr)

  end subroutine get_attribute_arr_r8

  !
  ! read integer global array
  !
  subroutine read_dataset_int(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    integer, intent(out)         :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! read
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dread_f(dataset, H5T_NATIVE_INTEGER, data, dims, hdferr, mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine read_dataset_int

  !
  ! read real(4) global array
  !
  subroutine read_dataset_r4(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(4), intent(out)         :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! read
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dread_f(dataset, H5T_NATIVE_REAL, data, dims, hdferr, mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine read_dataset_r4

  !
  ! read real(8) global array
  !
  subroutine read_dataset_r8(dest, name, rank, dims, count, &
       & loffset, goffset, data)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name
    integer, intent(in)          :: rank
    integer(hsize_t), intent(in) :: dims(rank)
    integer(hsize_t), intent(in) :: goffset(rank)
    integer(hsize_t), intent(in) :: loffset(rank)
    integer(hsize_t), intent(in) :: count(rank)
    real(8), intent(out)         :: data(:)

    integer(hid_t) :: dataset, dspace, mspace, mpio

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, dspace, hdferr)
    call h5screate_simple_f(rank, dims, mspace, hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F, mpio, hdferr)

    ! destination hyperslab
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, goffset, count, hdferr)

    ! source hyperslab
    call h5sselect_hyperslab_f(mspace, H5S_SELECT_SET_F, loffset, count, hdferr)

    ! read
    call h5pset_dxpl_mpio_f(mpio, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call h5dread_f(dataset, H5T_NATIVE_DOUBLE, data, dims, hdferr, mspace, dspace, mpio)

    call h5pclose_f(mpio, hdferr)
    call h5sclose_f(dspace, hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine read_dataset_r8

  !
  ! extend last dimension by one
  !
  subroutine extend_dimension(dest, name)
    implicit none
    integer(hid_t), intent(in)   :: dest
    character(len=*), intent(in) :: name

    integer :: rank
    integer(hid_t) :: dataset, space
    integer(hsize_t) :: dims(MAXDIM), dimm(MAXDIM)

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, space, hdferr)
    call h5sget_simple_extent_ndims_f(space, rank, hdferr)
    call h5sget_simple_extent_dims_f(space, dims, dimm, hdferr)

    dims(rank) = dims(rank) + 1
    call h5dset_extent_f(dataset, dims, hdferr)

    call h5sclose_f(space, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine extend_dimension

  !
  ! get array dimension
  !
  subroutine get_dimension(dest, name, ndim, dims)
    implicit none
    integer(hid_t), intent(in)    :: dest
    character(len=*), intent(in)  :: name
    integer, intent(out)          :: ndim
    integer(hsize_t), intent(out) :: dims(:)

    integer(hid_t) :: dataset, space
    integer(hsize_t) :: dimm(MAXDIM)

    call h5dopen_f(dest, name, dataset, hdferr, H5P_DEFAULT_F)
    call h5dget_space_f(dataset, space, hdferr)
    call h5sget_simple_extent_ndims_f(space, ndim, hdferr)

    if( ndim <= size(dims) ) then
       ! get dimensions
       call h5sget_simple_extent_dims_f(space, dims, dimm, hdferr)
    else
       ! show error message but try to continue
       write(0) 'Error in subroutine get_dimension'
    end if

    call h5sclose_f(space, hdferr)
    call h5dclose_f(dataset, hdferr)

  end subroutine get_dimension

end module h5util

#else

module h5util
  ! empty module
end module h5util

#endif
