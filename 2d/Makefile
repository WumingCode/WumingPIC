# -*- Makefile -*-
include common.mk

SUBDIRS = utils common

default:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

clean :
	rm -f $(OBJS) $(WM_INCLUDE)/*.mod $(WM_LIB)/*.a *.i *.mod *.out
	# clean subdirectories
	for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done
