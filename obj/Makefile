# -*- Makefile -*-
SHELL=/bin/sh
CC=mpicc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -fexpensive-optimizations -msse3 -I ../inc -fpic -DOPENMP -fopenmp -D_MPI
LDLIBS=-lm

# list of source files
SRC=prototypes.h datastructures.h velverlet.c force.c utilities.c output.c ljmd.c
VPATH=../src ../inc
OBJ=$(SRC:%.c=%.o)

default: ../ljmd.x ../obj/libljmd.so

clean:
	rm -f *.mod *.o ../ljmd.x ../obj/libljmd.so

# linker rule
../ljmd.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

../obj/libljmd.so: utilities.o force.o velverlet.o output.o
	$(CC) -shared -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
utilities.o: prototypes.h datastructures.h
force.o: prototypes.h datastructures.h
velverlet.o: prototypes.h datastructures.h
output.o: prototypes.h datastructures.h
ljmd.o: force.o utilities.o output.o
