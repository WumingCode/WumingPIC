# -*- Makefile -*-
include ../../common.mk
include ../common3d.mk

SRCS   = \
	wuming3d.f90 boundary_periodic.f90 field.f90 particle.f90 mom_calc.f90 \
    sort.f90 mpi_set.f90 paraio.f90
OBJS   = $(SRCS:%.f90=%.o)
MODS   = $(SRCS:%.f90=%.mod)

default: $(OBJS)
	$(AR) r $(WM_LIB3D)/libwuming3d_common.a $(OBJS)
	cp $(MODS) $(WM_INCLUDE3D)

clean :
	rm -f $(TARGET) $(OBJS) *.i *.mod *.out

# Dependencies
h5io.o: mpi_set.o
wuming3d.o: field.o particle.o mom_calc.o sort.o \
	mpi_set.o paraio.o
