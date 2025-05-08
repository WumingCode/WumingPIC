# -*- Makefile -*-
include common.mk

SUBDIRS = utils 2d 3d
LIBS    = $(WM_LIB)/libwuming_utils.a

default: $(LIBS) 2d 3d

$(LIBS):
	$(MAKE) -C utils

2d: $(LIBS)
	$(MAKE) -C 2d

3d: $(LIBS)
	$(MAKE) -C 3d

clean:
	rm -f $(WM_LIB)/*.a $(WM_INCLUDE)/*.mod
	# clean subdirectories
	for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done
