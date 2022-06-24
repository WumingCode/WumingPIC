# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpicc -fopenmp -mcmodel=medium
FC      = mpif90 -fopenmp -mcmodel=medium
FCFLAGS = -cpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
