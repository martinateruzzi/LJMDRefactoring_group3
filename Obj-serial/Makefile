# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -I ../inc -fpic
LDLIBS=-lm

# list of source files
SRC=prototypes.h datastructures.h velverlet.c force.c utilities.c output.c ljmd.c
VPATH=../src ../inc
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.x ../libljmd.so

clean:
	rm -f *.mod *.o ../ljmd-serial.x

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

../libljmd.so: utilities.o force.o velverlet.o output.o
	$(CC) -shared -o $@ $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
utilities.o: prototypes.h datastructures.h
force.o: prototypes.h datastructures.h
velverlet.o: prototypes.h datastructures.h
output.o: prototypes.h datastructures.h
ljmd.o: force.o utilities.o output.o
