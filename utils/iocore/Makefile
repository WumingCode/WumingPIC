# -*- Makefile -*-
include ../../common.mk

SRCS    = mpiio.f90 jsonio.f90 h5util.f90
OBJS    = $(SRCS:%.f90=%.o)
MODS    = $(SRCS:%.f90=%.mod)
EXES    = jsonio_test mpiio_test1 mpiio_test2

default: $(OBJS)
	$(AR) r $(WM_LIB)/libwuming_utils.a $(OBJS)
	cp $(MODS) $(WM_INCLUDE)

test: $(EXES)
	python unittest.py

mpiio_test1: mpiio_test1.o mpiio_test.o mpiio.o
	$(FC) -o $@ $^

mpiio_test2: mpiio_test2.o mpiio_test.o mpiio.o
	$(FC) -o $@ $^

jsonio_test: jsonio.o jsonio_test.o
	$(FC) -o $@ $^ $(LDFLAGS) -lwuming_utils

clean :
	rm -f $(EXES) $(OBJS) *.o *.i *.mod *.out mpiio_test.raw mpiio_test.h5

# Dependencies
jsonio_test.o: jsonio.o
mpiio_test1.o: mpiio.o mpiio_test.o
mpiio_test2.o: mpiio.o mpiio_test.o
