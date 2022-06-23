# -*- Makefile -*-

# HDF5
#HDF5DIR = $(shell spack location -i hdf5 %intel)

# compilers and arguments
AR      = ar
CC      = mpiicc -qopenmp
FC      = mpiifort -qopenmp
FCFLAGS = -fpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib
