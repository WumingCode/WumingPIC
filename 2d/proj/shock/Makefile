# -*- Makefile -*-
include ../../../common.mk
include ../../common2d.mk

SRCS   = main.f90 app.f90 boundary_shock.f90
OBJS   = $(SRCS:%.f90=%.o)
TARGET = main.out

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS) -lwuming2d_common -lwuming_utils $(HDF5LIB)

clean :
	rm -f $(TARGET) $(OBJS) *.i *.mod *.out

# Dependencies
main.o: app.o
app.o: boundary_shock.o $(WM_LIB2D)/libwuming2d_common.a $(WM_LIB)/libwuming_utils.a
