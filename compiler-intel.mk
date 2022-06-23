# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpiicc -qopenmp
FC      = mpiifort -qopenmp
FCFLAGS = -fpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib
