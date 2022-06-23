program mpiio_test1
  use mpi
  use mpiio
  use mpiio_test
  implicit none

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
  call mpiio_open_file(filename, file, disp, 'w')

  !
  ! atomic write
  !
  endian = mpiio_get_endian_flag()
  call mpiio_write_atomic(file, disp, endian)
  call mpiio_write_atomic(file, disp, nproc)
  call mpiio_write_atomic(file, disp, dims)
  call mpiio_write_atomic(file, disp, nx)
  call mpiio_write_atomic(file, disp, ny)

  call mpiio_write_atomic(file, disp, var_i4)
  call mpiio_write_atomic(file, disp, var_i8)
  call mpiio_write_atomic(file, disp, var_r4)
  call mpiio_write_atomic(file, disp, var_r8)

  char = '*** Character Description ***'
  charlen = len_trim(char)
  call mpiio_write_atomic(file, disp, charlen)
  call mpiio_write_atomic(file, disp, char)

  ! close once
  call mpiio_close_file(file)


  if( nrank == 0 ) then
     mx = nx * dims(1)
     my = ny * dims(2)
     write(*, '(a, " : ", i8)') formatstr('# endian flag', 30), endian
     write(*, '(a, " : ", i8)') formatstr('# MPI process', 30), nproc
     write(*, '(a, " : ", i8)') formatstr('# MPI process in x', 30), dims(1)
     write(*, '(a, " : ", i8)') formatstr('# MPI process in y', 30), dims(2)
     write(*, '(a, " : ", i8)') formatstr('# local grid in x', 30), nx
     write(*, '(a, " : ", i8)') formatstr('# local grid in y', 30), ny
     write(*, '(a, " : ", i8)') formatstr('# global grid in x', 30), mx
     write(*, '(a, " : ", i8)') formatstr('# global grid in y', 30), my
     write(*, '(a, " : ", i8)') formatstr('# integer(4)', 30), var_i4
     write(*, '(a, " : ", i8)') formatstr('# integer(8)', 30), var_i8
     write(*, '(a, " : ", f10.4)') formatstr('# real(4)', 30), var_r4
     write(*, '(a, " : ", f10.4)') formatstr('# real(8)', 30), var_r8
     write(*, '(a, " : ", i8)') formatstr('# character(*) len', 30), charlen
     write(*, '(a, " : ", a)') formatstr('# character(*)', 30), trim(char)
  end if

  ! now open with append mode
  call mpiio_open_file(filename, file, disp, 'a')

  !
  ! collective write
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
  call filldata(dat_i4, zero(1), nxs, nxe, nys, nye, nb, dims, coords)
  buf_i4 = reshape(dat_i4(nxs:nxe,nys:nye), (/nx*ny/))
  call mpiio_write_collective(file, disp, ndim, gshape, lshape, offset, buf_i4)

  ! integer(8)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for integer(8)', 30), disp
  end if
  call filldata(dat_i8, zero(2), nxs, nxe, nys, nye, nb, dims, coords)
  buf_i8 = reshape(dat_i8(nxs:nxe,nys:nye), (/nx*ny/))
  call mpiio_write_collective(file, disp, ndim, gshape, lshape, offset, buf_i8)

  ! real(4)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for real(4)', 30), disp
  end if
  call filldata(dat_r4, zero(3), nxs, nxe, nys, nye, nb, dims, coords)
  buf_r4 = reshape(dat_r4(nxs:nxe,nys:nye), (/nx*ny/))
  call mpiio_write_collective(file, disp, ndim, gshape, lshape, offset, buf_r4)

  ! real(8)
  if( nrank == 0 ) then
     write(*, '(a, " : ", i16)') formatstr('# offset for real(8)', 30), disp
  end if
  call filldata(dat_r8, zero(4), nxs, nxe, nys, nye, nb, dims, coords)
  buf_r8 = reshape(dat_r8(nxs:nxe,nys:nye), (/nx*ny/))
  call mpiio_write_collective(file, disp, ndim, gshape, lshape, offset, buf_r8)

  call mpiio_close_file(file)

  call MPI_Finalize(ierr)

end program mpiio_test1
