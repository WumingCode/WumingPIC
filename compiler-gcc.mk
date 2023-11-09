# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpicc -O3 -fopenmp -mcmodel=medium
FC      = mpif90 -O3 -fopenmp -mcmodel=medium
FCFLAGS = -cpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
