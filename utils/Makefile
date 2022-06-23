# -*- Makefile -*-
include ../common.mk

SUBDIRS = json iocore
SRCS    = wuming_utils.f90
OBJS    = $(SRCS:%.f90=%.o)
MODS    = $(SRCS:%.f90=%.mod)

default: $(OBJS) subdirs
	$(AR) r $(WM_LIB)/libwuming_utils.a $(OBJS)
	cp $(MODS) $(WM_INCLUDE)

subdirs:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

clean :
	rm -f $(OBJS) $(TARGET) *.i *.mod *.out
	# clean subdirectories
	for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done

json: iocore
