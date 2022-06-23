module wuming2d
  use json_module
  use jsonio
  use mpiio
  use field
  use particle
  use mom_calc
  use sort
  use mpi_set
  use fio
  !
  ! select default I/O routines
  !
#ifdef USE_HDF5
  use h5io, &
       & io__init     => h5io__init,     &
       & io__finalize => h5io__finalize, &
       & io__param    => h5io__param,    &
       & io__input    => h5io__input,    &
       & io__output   => h5io__output,   &
       & io__mom      => h5io__mom,      &
       & io__ptcl     => h5io__ptcl,     &
       & io__orb      => h5io__orb
  use paraio
#else
  use h5io
  use paraio, &
       & io__init     => paraio__init,     &
       & io__finalize => paraio__finalize, &
       & io__param    => paraio__param,    &
       & io__input    => paraio__input,    &
       & io__output   => paraio__output,   &
       & io__mom      => paraio__mom,      &
       & io__ptcl     => paraio__ptcl,     &
       & io__orb      => paraio__orb
#endif
  implicit none

end module wuming2d
