# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpicc -fopenmp -mcmodel=medium
FC      = mpif90 -fopenmp -mcmodel=medium
FCFLAGS = -cpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib

