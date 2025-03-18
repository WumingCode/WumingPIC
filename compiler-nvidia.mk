# -*- Makefile -*-

# HDF5
#HDF5DIR = $(shell spack location -i hdf5 %gcc)

# compilers and arguments
AR      = ar
CC      = mpicc -fast -mp -mcmodel=medium
FC      = mpif90 -fast -mp -mcmodel=medium -Mbackslash
FCFLAGS = -cpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)

