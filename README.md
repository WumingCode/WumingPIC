# Wuming PIC
<!--[![DOI](https://zenodo.org/badge/377835665.svg)](https://zenodo.org/badge/latestdoi/377835665)-->

Two- and Three-dimentional, special relativistic, electromagnetic particle-in-cell simulation code for general puposes in space and astrophysical plasmas.

![shock](/sample.png)

## Features
* Solves the Vlasov-Maxwell equations by the particle-in-cell method
* Buneman-Boris method for the equation of motion of particles.
* Implicit FDTD scheme for the Maxwell equation (Hoshino, ApJ, 2013)
* Esirkepov's charge conservation scheme for the current deposit with the 2nd-order shape function (Esirkepov, CPC, 2001)
* Written in Fortran 90/95
* Hybrid parallelization by MPI and OpenMP
   - 1D/2D domain decomposition in the y-/y-z- directions.
* SIMD optimization and efficient cache usage
* MPI-IO raw data output with JSON-based metadata
* Python scripts for HDF5 format convertor and quicklook

## Requirements
* Fortran compiler
  - We prepare setting files for Makefile for different compilers, including Cray compiler, Fujitsu compiler, GCC-Fortran, Intel Fortran, NVIDIA HPC Fortran
  - Because of the requirement of JSON-Fortran library used in this code (https://github.com/jacobwilliams/json-fortran), **GCC-Fortran's version must be greater than 4.9**.

* MPI library
  - We have tested with MPICH and Open MPI
  - The code works with the vender's MPI library on Fujitsu's supercomputer systems
   
* Python [OPTIONAL]
  - The code generates raw binary data and corresponding JSON-based metadata files. A Python script in each physics directory converts them to HDF5 files as a post process
  - A sample python script is prepared for quick look of the results
  - [Here](pyvista_demo.rst) is a demonstration of 3D visualazation using PyVista library

## Installation
```bash
$ git clone git@github.com:WumingCode/WumingPIC.git
```

## Code structure
``` 
WumingPIC/
├── 2d
│   ├── common
│   ├── include
│   ├── lib
│   └── proj
│       ├── reconnection
│       ├── shock
│       └── weibel
├── 3d
│   ├── common
│   ├── include
│   ├── lib
│   └── proj
│       ├── reconnection
│       ├── shock
│       └── weibel
│
├── CITATION.cff
│
├── LICENSE.txt
│
├── Makefile
│
├── README.md
│
├── common.mk
│
├── compiler-cray.mk
│
├── compiler-fujitsu.mk
│
├── compiler-gcc.mk
│
├── compiler-intel.mk
│
├── compiler-nvidia.mk
│
├── include
│   └── directory for module files
├── lib
│   └── directory for common library
├── python
│   ├── json2hdf5.py - A python script to convert JSON files to HDF5 metadata
│   └── jsoncheck.py - For checking reading json files
└── utils - some utility files for I/O functions
    ├── iocore
    └── json
```
In each 2d or 3d code, files are organized as follows

```
{2d,3d}
├── Makefile
│
├── common
│   └── common files of PIC algorithms
│
├── common{2d,3d}.mk
│
├── include
│   └── directory for module files
│
├── lib
│   └── directory for common library
│
└── proj
    ├── shock
    │   └── collsion-less shock simulation setup files and scripts for post process
    ├── weibel
    │   └── Weibel instability simulation setup files and scripts for post process
    └── reconnection
        └── magnetic reconnection simulation setup files and scripts for post process
```

## Preparation
1. Move to the installed directory.

   ```bash
   $ cd ./WumingPIC
   ```

2. Copy one of `comiler-*.mk` files depending on your compiler environment to `compiler.mk`.
   For instance, copy `compiler-gcc.mk` if you are using gfortran.

   ```bash
   $ cp compiler-gcc.mk compiler.mk
   ```

3. Make common libraries.
   Compile via

   ```bash
   $ make
   ```

   and make sure `libwuming*.a` are genererated in the library directories of `lib/` and `{2d,3d}/lib/`.
   You are now ready for executing a specific physics problem.

## Physics Problems

Following physics problem setups are available at present.
* Weibel instability ([2d](2d/proj/weibel/README.md), [3d](3d/proj/weibel/README.md))
* Collision-less shock ([2d](2d/proj/shock/README.md), [3d](3d/proj/shock/README.md))
* Magnetic reconnection ([2d](2d/proj/reconnection/README.md), [3d](3d/proj/reconnection/README.md))

### How to run
Go to one of the physics problem directories `{2d,3d}/proj/*` and make an executable `main.out`.
For instance,

```bash
$ cd 2d/proj/weibel
$ make
```

will create an executable for 2D Weibel instability. This will read parameters from a configuration file in JSON format. You may copy a sample configuration file `config_sample.json` to `config.json`:

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
$ python ../../../python/json2hdf5.py *.json
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
WumingPIC code uses 
* [JSON-Fortran](https://github.com/jacobwilliams/json-fortran) API for reading/writing JSON files from Fortran.
* [Amano's MPI-IO, JSON, HDF5 utitlity files](https://github.com/amanotk)

## License
WumingPIC code is distributed under [the MIT license](LICENSE.txt).

## Cite as
<!--[![DOI](https://zenodo.org/badge/377835665.svg)](https://zenodo.org/badge/latestdoi/377835665)-->

Cite https://zenodo.org/doi/10.5281/zenodo.10990575, which represents all versions and will always resolve to the latest release.
