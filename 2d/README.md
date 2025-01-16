# Wuming PIC2D
[![DOI](https://zenodo.org/badge/377835665.svg)](https://zenodo.org/badge/latestdoi/377835665)

Two-dimentional, special relativistic, electromagnetic particle-in-cell simulation code for general puposes in space and astrophysical plasmas.

## Features
* Solves the Vlasov-Maxwell equations by the particle-in-cell method
* Buneman-Boris method for the equation of motion of particles.
* Implicit FDTD scheme for the Maxwell equation (Hoshino, ApJ, 2013)
* Esirkepov's charge conservation scheme for the current deposit with the 2nd-order shape function (Esirkepov, CPC, 2001)
* Written in Fortran 90/95
* Hybrid parallelization by MPI and OpenMP
   - 1D domain decomposition in the x direction
* SIMD optimization and efficient cache usage
* MPI-IO raw data output with JSON-based metadata
* Python scripts for HDF5 format convertor and quicklook

## Requirements
* Fortran compiler
  - We prepare setting files for Makefile for different compilers, including Intel Fortran, GCC-Fortran, Fujitsu compiler
  - Because of the requirement of JSON-Fortran library used in this code (https://github.com/jacobwilliams/json-fortran), **GCC-Fortran's version must be greater than 4.9**.

* MPI library
  - We have tested with MPICH and Open MPI
  - The code works with the vender's MPI library on Fujitsu's supercomputer systems
   
* Python [OPTIONAL]
  - The code generates raw binary data and corresponding JSON-based metadata files. A Python script in each physics directory converts them to HDF5 files as a post process
  - A sample python script is prepared for quick look of the results

## Installation
```bash
$ git clone git@github.com:WumingCode/WumingPIC2D.git
```

## Code structure
``` 
WumingPIC2D
├── Makefile
│
├── README.md
│
├── common
│   └── common files of PIC algorithms
│
├── common.mk
│
├── compiler-fujitsu.mk
│
├── compiler-gcc.mk
│
├── compiler-intel.mk
│
├── include
│   └── directory for module files
│
├── lib
│   └── directory for common library
│
├── proj
│   ├── shock
│   │   └── collsion-less shock simulation setup files and scripts for post process
│   ├── weibel
│   │   └── Weibel instability simulation setup files and scripts for post process
│   └── reconnection
│       └── magnetic reconnection simulation setup files and scripts for post process
│
├── python
│   └── json2hdf5.py - A python script to convert JSON files to HDF5 metadata
│
└── utils
    └── utility files for MPI-IO and JSON output
```

## Preparation
1. Move to the installed directory.  

   ```bash
   $ cd ./WumingPIC2D
   ```

2. Copy one of `comiler-*.mk` files depending on your compiler environment to `compiler.mk`.  
   For instance, copy `compiler-gcc.mk` if you are using gfortran.

   ```bash
   $ cp compiler-gcc.mk compiler.mk
   ```

3. Make a common library.  
   Compile via

   ```bash
   $ make
   ```

   and make sure `libwuming*.a` are genererated in the library directory `lib/`.  
   You are now ready for executing a specific physics problem.

## Physics Problems

Following physics problem setups are available at present.
* [Weibel instability](proj/weibel/README.md)
* [Collision-less shock](proj/shock/README.md)
* [Magnetic reconnection](proj/reconnection/README.md)

### How to run
Go to one of the physics problem directories `proj/*` and make an executable `main.out`.  
For instance,

```bash
$ cd proj/weibel
$ make
```

will create an executable for Weibel instability. This will read parameters from a configuration file in JSON format. You may copy a sample configuration file `config_sample.json` to `config.json`:

```bash
$ cp config_sample.json config.json
```

and edit it as you like. By default, running the code via

```bash
$ mpiexec -np 4 ./main.out
```

will try to read `config.json`.  

If you want, you may explicitly specify the filename with a command line argument:

```bash
$ mpiexec -np 4 ./main.out somethingelse.json
```

in which case `sometingelse.json` will be read.

### Configuration Parameters
Configuration parameters that are common for all physics problems are as follows.

- `config`
  - `verbose`  
     Print verbose messages if >= 1.
  - `datadir`  
     Data directory to which all the output will be saved.
  - `max_elapsed`  
     Maximum elapsed time. The code will stop automatically when the elapsed  
     time goes beyond this limit. A snapshot will be saved for restart.
  - `max_it`  
     Maximum number of iteration step.
  - `intvl_ptcl`  
     Interval of time step for entire particle data output.
  - `intvl_mom`  
     Interval of time step for for moment data output.
  - `intvl_orb`  
     Interval of time step for tracer particle data output.
  - `restart_file`  
     Snapshot file to be read during initialzation. If this is specified,  
     the code will start from this state. Otherwise, it will start from the  
     initial condition. Note that this parameter will be overwritten by the  
     code automatically when it finishes.

For specific problem-dependent parameters, please refer `README.md` in each of physics problem directories.


### How it works
The code will produce a lot of files in a specified directory. The output should
appear in pairs of `*.json` and `*.raw` files. A `.json` file in JSON format
describes meta data and how the actual simulation data are stored in a `.raw`
file, which contains raw data in binary format.

The JSON files can be processed to generate HDF5 format files for data analysis
via a script `json2hdf5.py` which is located in `python/` directory.
For instance in the working directory,

```bash
$ python ../../python/json2hdf5.py *.json
```

will process all JSON files in the current directory and generate HDF5 format
files for each JSON format file. You can use the generated HDF5 files for data
analysis.

When the code running time goes beyond a given elapsed time limit, it will save
a snapshot data and stop running. By specifying a snapshot in the configuration
file, you may restart the run.

By default, the configuration file will be overwritten by the code so that you
can restart the run next time with the exactly same command.

For instance, if you run the code via

```bash
$ mpiexec -np 4 ./main.out
```

and the elapsed time limit is reached, it will overwrite the configuration file
`config.json` to properly set `restart_file` option.
In the next time, you may run again via

```bash
$ mpiexec -np 4 ./main.out
```

then the previous snapshot data will be read automatically.

## Support
Join Slack workspace via https://join.slack.com/t/wumingpic/shared_invite/zt-xlm8cixg-NOV33dyorO1Whc4~FcVJ0g .

## Credits
WumingPIC2D code uses 
* [JSON-Fortran](https://github.com/jacobwilliams/json-fortran) API for reading/writing JSON files from Fortran.
* [Amano's MPI-IO, JSON, HDF5 utitlity files](https://github.com/amanotk)

## License
WumingPIC2D code is distributed under [the MIT license](LICENSE.txt).

## Cite as
[![DOI](https://zenodo.org/badge/377835665.svg)](https://zenodo.org/badge/latestdoi/377835665)

Cite https://zenodo.org/doi/10.5281/zenodo.10990575, which represents all versions and will always resolve to the latest release.
