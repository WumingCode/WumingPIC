# Magnetic Reconnection
To minimize the workload imbalance, we set the configuration of the background magnetic field in the y-direction contrary to usual setups.

## Configuration Parameters
In addition to common configuration parameters described in `README.md`
in the top directory, following parameters can be specified.

- `parameter` :
  - `num_process`  
     Number of MPI processes. This must be the same as the number specified
     in execusion with `mpiexec` or `mpirun`.
  - `num_process_j`  
     Number of MPI processes in y direction. Number of MPI processes in z direction is
     determined by `num_process`/`num_process_j`.
  - `n_x`  
     Number of grid in x direction.
  - `n_y`  
     Number of grid in y direction.
  - `n_z`  
     Number of grid in z direction.
  - `mass_ratio`  
     Ion to electron mass ratio: m_i/m_e.
  - `alpha`  
     Ratio of the electron plasma-to-cyclotron frequency.
  - `rtemp`  
     Ion-to-Electron temperature ratio.
  - `lcs`  
     Current sheet thickness in the unit of the ion inertia length
  - `nbg`  
     Number of particles per cell per species for the background population
  - `ncs`  
     Number of particles per cell per species for the current sheet population

## Checking Results
For a quicklook, run the script `plot_sample.py` as

```bash
$ python plot_sample.py
```
