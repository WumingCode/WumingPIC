# Collisionless Shock

## Configuration Parameters
In addition to common configuration parameters described in `README.md`  
in the top directory, following parameters can be specified.

- `config`
  - `intvl_expand`  
     Interval of time step for expanding the box size in x direction.
- `parameter` :
  - `num_process`  
     Number of MPI processes. This must be the same as the number specified  
     in execution with `mpiexec` or `mpirun`.
  - `n_ppc`  
     Number of particles per cell in the upstream region.
  - `n_x`  
     Number of grids in x direction.
  - `n_y`  
     Number of grids in y direction.
  - `n_x_ini`  
     Initial number of grids in x direction. The simulation box will expand  
     every `intvl_expand` step.
  - `u_injection`  
     Injection four-velocity from the boundary.
  - `mass_ratio`  
     Ion to electron mass ratio: m_i/m_e.
  - `sigma_e`  
     Ratio of electron cyclotron-to-plasma frequency squared.
  - `omega_pe`  
     Proper plasma frequency in the upstream.
  - `v_thi`  
     Thermal velocity of ions in the upstream.
  - `v_the`  
     Thermal velocity of electrons in the upstream.
  - `theta_bn`  
     Magnetic obliquity in the upstream or the shock angle in degrees.
  - `phi_bn`  
     Magnetic field orientation in the upstream with respect to the  
     simulation plane in degrees. The in-plane and out-of-plane magnetic  
     field configurations correspond to 0 and 90, respectively.
  - `l_damp_ini`  
     Length scale for the initial velocity profile in unit of grid size.  
     This is to smooth out initial disturbances and nothing physical.

## Checking Results
For a quicklook, run the script `plot_sample.py` as

```bash
$ python plot_sample.py
```
