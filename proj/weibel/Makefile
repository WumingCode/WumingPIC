include ../../common.mk

SRCS   = main.f90 app.f90
OBJS   = $(SRCS:%.f90=%.o)
TARGET = main.out

$(TARGET): $(OBJS)
	         $(FC) -o $@ $^ $(LDFLAGS) -lwuming2d_common -lwuming_utils

clean:
	    rm -f $(TARGET) $(OBJS) *.i *.mod *.out

main.o: app.o
