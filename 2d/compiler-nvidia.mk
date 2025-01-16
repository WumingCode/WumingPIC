# -*- Makefile -*-

# HDF5
#HDF5DIR = $(shell spack location -i hdf5 %gcc)

# compilers and arguments
AR      = ar
CC      = mpicc -Mcuda -O3 -mp -mcmodel=medium
FC      = mpif90 -Mcuda -O3 -mp -mcmodel=medium -Mbackslash
FCFLAGS = -cpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib

