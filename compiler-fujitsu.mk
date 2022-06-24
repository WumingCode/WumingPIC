# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpifccpx -Kfast,openmp
FC      = mpifrtpx -Kfast,openmp
FCFLAGS = -Nalloc_assign -Cpp -I$(WM_INCLUDE)
LDFLAGS = -L$(WM_LIB)
