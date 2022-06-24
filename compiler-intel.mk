# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpiicc -qopenmp
FC      = mpiifort -qopenmp
FCFLAGS = -fpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
