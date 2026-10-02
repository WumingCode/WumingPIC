# -*- Makefile -*-

# compilers and arguments
AR      = ar
FC      = mpif90 -O3 -fopenmp -mcmodel=medium
FCFLAGS = -cpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
