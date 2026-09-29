# -*- Makefile -*-

# compilers and arguments
AR      = ar
FC      = ftn -h omp
FCFLAGS = -dynamic -em -ef -eZ -eF -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
