# -*- Makefile -*-
include ../../../common.mk
include ../../common3d.mk

SRCS   = main.f90 app.f90
OBJS   = $(SRCS:%.f90=%.o)
TARGET = main.out

$(TARGET): $(OBJS)
	         $(FC) -o $@ $^ $(LDFLAGS) -lwuming3d_common -lwuming_utils

clean:
	    rm -f $(TARGET) $(OBJS) *.i *.mod *.out

main.o: app.o
app.o: $(WM_LIB3D)/libwuming3d_common.a $(WM_LIB)/libwuming_utils.a
