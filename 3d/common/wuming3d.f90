module wuming3d
  use json_module
  use jsonio
  use mpiio
  use field
  use particle
  use mom_calc
  use sort
  use mpi_set
  use paraio, &
       & io__init     => paraio__init,     &
       & io__finalize => paraio__finalize, &
       & io__param    => paraio__param,    &
       & io__input    => paraio__input,    &
       & io__output   => paraio__output,   &
       & io__mom      => paraio__mom,      &
       & io__ptcl     => paraio__ptcl,     &
       & io__orb      => paraio__orb

  implicit none

end module wuming3d
