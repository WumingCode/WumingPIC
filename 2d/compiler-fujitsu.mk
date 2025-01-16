# -*- Makefile -*-

# HDF5
#HDF5DIR = $(shell spack location -i hdf5 %fj /7cxah5)

# compilers and arguments
AR      = ar
CC      = mpifccpx -Kfast,openmp
FC      = mpifrtpx -Kfast,openmp
FCFLAGS = -Nalloc_assign -Cpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib
