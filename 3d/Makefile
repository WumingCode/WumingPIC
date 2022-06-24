# -*- Makefile -*-
include ../common.mk

SUBDIRS = common

default:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

clean:
	rm -f $(WM_LIB3D)/*.a $(WM_INCLUDE3D)/*.mod
	# clean subdirectories
	for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done
