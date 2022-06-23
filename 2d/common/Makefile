# -*- Makefile -*-
include ../../common.mk
include ../common2d.mk

SRCS   = \
	wuming2d.f90 boundary_periodic.f90 field.f90 particle.f90 mom_calc.f90 \
    sort.f90 mpi_set.f90 fio.f90 h5io.f90 paraio.f90
OBJS   = $(SRCS:%.f90=%.o)
MODS   = $(SRCS:%.f90=%.mod)

default: $(OBJS)
	$(AR) r $(WM_LIB2D)/libwuming2d_common.a $(OBJS)
	cp $(MODS) $(WM_INCLUDE2D)

clean :
	rm -f $(TARGET) $(OBJS) *.i *.mod *.out

# Dependencies
h5io.o: mpi_set.o
wuming2d.o: field.o particle.o mom_calc.o sort.o \
	mpi_set.o fio.o h5io.o paraio.o
