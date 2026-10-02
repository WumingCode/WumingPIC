# -*- Makefile -*-

# compilers and arguments
AR      = ar
FC      = mpiifort -qopenmp
FCFLAGS = -fpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
