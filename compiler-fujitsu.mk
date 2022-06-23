# -*- Makefile -*-

# compilers and arguments
AR      = ar
CC      = mpifccpx -Kfast,openmp
FC      = mpifrtpx -Kfast,openmp
FCFLAGS = -Nalloc_assign -Cpp -I$(BASEDIR)/include
LDFLAGS = -L$(BASEDIR)/lib
