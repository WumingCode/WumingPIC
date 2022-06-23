program mpiio_test2
  use mpi
  use mpiio
  use mpiio_test
  implicit none

  integer :: r_nproc, r_nrank, r_dims(ndim), r_nx, r_ny

  call MPI_Init(ierr)
  call MPI_Comm_size(mpi_comm_world, nproc, ierr)
  call MPI_Comm_rank(mpi_comm_world, nrank, ierr)

  call process_decomp(dims)

  if ( nproc /= product(dims) ) then
     ! exit with error
     if ( nrank == 0 ) then
        print *, 'error !'
     end if

     call MPI_Finalize(ierr)
     stop
  end if

  call MPI_Cart_create(mpi_comm_world, ndim, dims, &
       & (/.true., .true./), .false., comm, ierr)
  call MPI_Cart_coords(comm, nrank, ndim, coords, ierr)

  ! open
  call mpiio_open_file(filename, file, disp, 'r')

  !
  ! atomic read
  !
  call mpiio_read_atomic(file, disp, endian)
  call mpiio_read_atomic(file, disp, r_nproc)
  call mpiio_read_atomic(file, disp, r_dims)
  call mpiio_read_atomic(file, disp, r_nx)
  call mpiio_read_atomic(file, disp, r_ny)

  call mpiio_read_atomic(file, disp, var_i4)
  call mpiio_read_atomic(file, disp, var_i8)
  call mpiio_read_atomic(file, disp, var_r4)
  call mpiio_read_atomic(file, disp, var_r8)
  call mpiio_read_atomic(file, disp, charlen)
  call mpiio_read_atomic(file, disp, char, charlen)

  ! check endian
  if ( .not. endian == mpiio_get_endian_flag() ) then
     write(0, *) 'endian does not match !'
     call MPI_Finalize(ierr)
     stop
  end if

  ! check MPI process etc are ok
  if(    nproc /= r_nproc .or. &
       & dims(1) /= r_dims(1) .or. &
       & dims(2) /= r_dims(2) .or. &
       & nx /= r_nx .or. &
       & ny /= r_ny ) then
     write(0, *) 'MPI processes or data size are wrong !'
     call MPI_Finalize(ierr)
     stop
  end if

  if( nrank == 0 ) then
     mx = r_nx * r_dims(1)
     my = r_ny * r_dims(2)
     write(*, '(a, " : ", i8)') formatstr('# endian flag', 30), endian
     write(*, '(a, " : ", i8)') formatstr('# MPI process', 30), r_nproc
     write(*, '(a, " : ", i8)') formatstr('# MPI process in x', 30), r_dims(1)
     write(*, '(a, " : ", i8)') formatstr('# MPI process in y', 30), r_dims(2)
     write(*, '(a, " : ", i8)') formatstr('# local grid in x', 30), r_nx
     write(*, '(a, " : ", i8)') formatstr('# local grid in y', 30), r_ny
     write(*, '(a, " : ", i8)') formatstr('# global grid in x', 30), mx
     write(*, '(a, " : ", i8)') formatstr('# global grid in y', 30), my
     write(*, '(a, " : ", i8)') formatstr('# integer(4)', 30), var_i4
     write(*, '(a, " : ", i8)') formatstr('# integer(8)', 30), var_i8
     write(*, '(a, " : ", f10.4)') formatstr('# real(4)', 30), var_r4
     write(*, '(a, " : ", f10.4)') formatstr('# real(8)', 30), var_r8
     write(*, '(a, " : ", i8)') formatstr('# character(*) len', 30), charlen
     write(*, '(a, " : ", a)') formatstr('# character(*)', 30), trim(char)
  end if

  !
  ! collective read
  !
  gshape(1) = nx * dims(1)
  gshape(2) = ny * dims(2)
  lshape(1) = nx
  lshape(2) = ny
  offset(1) = nx * coords(1)
  offset(2) = ny * coords(2)

  ! integer(4)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for integer(4)', 30), disp
  end if
  call mpiio_read_collective(file, disp, ndim, gshape, lshape, offset, buf_i4)
  dat_i4(nxs:nxe,nys:nye) = reshape(buf_i4, (/nx, ny/))

  ! integer(8)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for integer(8)', 30), disp
  end if
  call mpiio_read_collective(file, disp, ndim, gshape, lshape, offset, buf_i8)
  dat_i8(nxs:nxe,nys:nye) = reshape(buf_i8, (/nx, ny/))

  ! real(4)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for real(4)', 30), disp
  end if
  call mpiio_read_collective(file, disp, ndim, gshape, lshape, offset, buf_r4)
  dat_r4(nxs:nxe,nys:nye) = reshape(buf_r4, (/nx, ny/))

  ! real(8)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for real(8)', 30), disp
  end if
  call mpiio_read_collective(file, disp, ndim, gshape, lshape, offset, buf_r8)
  dat_r8(nxs:nxe,nys:nye) = reshape(buf_r8, (/nx, ny/))

  ! check data validity
  call testdata(dat_i4, zero(1), nxs, nxe, nys, nye, nb, dims, coords)
  call testdata(dat_i8, zero(2), nxs, nxe, nys, nye, nb, dims, coords)
  call testdata(dat_r4, zero(3), nxs, nxe, nys, nye, nb, dims, coords)
  call testdata(dat_r8, zero(4), nxs, nxe, nys, nye, nb, dims, coords)

  call mpiio_close_file(file)

  call MPI_Finalize(ierr)

end program mpiio_test2
