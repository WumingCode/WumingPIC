# -*- Makefile -*-

# base directory
BASEDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# include compilers
include $(BASEDIR)/compiler.mk

# default
.PHONY : all
.PHONY : clean

.SUFFIXES :
.SUFFIXES : .o .f90

%.o: %.f90
	$(FC) -c $(FCFLAGS) $< -o $@

%.o: %.F90
	$(FC) -c $(FCFLAGS) $< -o $@


# Wuming
WM_INCLUDE    = $(BASEDIR)/include
WM_LIB        = $(BASEDIR)/lib
