# -*- Makefile -*-

SHELL=/bin/sh

CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

INCLUDE= -I ../src -I integration -I Parallel-OPenmp

# list of source files
SRC=prototypes.h datastructures.h velverlet.c force.c utilities.c output.c test_integration.c test_kinetic.c

VPATH=../src integration
OBJ=$(SRC:%.c=%.o)

default: integration/test_integration.x  Parallel-OPenmp/test_kinetic.x

clean:
	rm -f *.mod *.o integration/test_integration.x

# linker rule
integration/test_integration.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS) $(INCLUDE)

Parallel-OPenmp/test_kinetic.x: $(OBJ)
        $(CC) -o $@ $(CFLAGS) $^ $(LDLIBS) $(INCLUDE)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $< $(INCLUDE)

# dependencies
utilities.o: prototypes.h datastructures.h
force.o: prototypes.h datastructures.h
velverlet.o: prototypes.h datastructures.h
output.o: prototypes.h datastructures.h
test_integration.o: prototypes.h datastructures.h velverlet.c force.c utilities.c output.c test_integration.c
