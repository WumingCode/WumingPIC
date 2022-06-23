# -*- Makefile -*-

# base directory
BASEDIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# include compilers
include $(BASEDIR)/compiler.mk

# add HDF5 options if HDF5DIR is specified
ifdef HDF5DIR
	CC      += -DUSE_HDF5
	FC      += -DUSE_HDF5
	FCFLAGS += -I$(HDF5DIR)/include
	LDFLAGS += -L$(HDF5DIR)/lib
	HDF5LIB  = -lhdf5 -lhdf5_fortran
endif

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
WM_INCLUDE2D  = $(BASEDIR)/2d/include
WM_INCLUDE3D  = $(BASEDIR)/3d/include
WM_LIB        = $(BASEDIR)/lib
WM_LIB2D      = $(BASEDIR)/2d/lib
WM_LIB3D      = $(BASEDIR)/3d/lib
