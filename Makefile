EXEDIR = ./
FC = mpif90
FFLAGS = -O3 -openmp
OBJS = fio.o particle.o field.o boundary.o mpi_set.o const.o init.o main.o sort.o mom_calc.o

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	$(FC) -c $< $(FFLAGS)

all: 
	@echo 'make test or make 1 or make 2 or make 3'

test: $(OBJS)
	$(FC) -o $(EXEDIR)em3dsk_$@.out $(FFLAGS) $(OBJS)

1: $(OBJS)
	$(FC) -o $(EXEDIR)em3dsk$@.out $(FFLAGS) $(OBJS)

2: $(OBJS)
	$(FC) -o $(EXEDIR)em3dsk$@.out $(FFLAGS) $(OBJS)

3: $(OBJS)
	$(FC) -o $(EXEDIR)em3dsk$@.out $(FFLAGS) $(OBJS)

# Dependencies
field.o : boundary.o
main.o : init.o const.o mpi_set.o boundary.o fio.o particle.o field.o sort.o mom_calc.o
init.o : const.o mpi_set.o boundary.o fio.o sort.o

clean :
	rm -f $(OBJS) $(TARGET) *.mod *.out
