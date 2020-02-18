#-*- Makefile -*-
SHELL=/bin/sh
CC=mpicc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm
OFLAGS=-D_MPI
# list of source files
SRC=ljmd.c \
    force.c \
    utilities.c \
    output.c \
    velverlet.c \

IPATH=../inc/ 
#VPATH=../src/new 

OBJ=$(SRC:%.c=%.o)
EXE=main_mpi.x

default: ../$(EXE)

clean:
	rm -f *.mod *.o ../main_mpi.x

# linker rule
../$(EXE): $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c 
	$(CC) -c $(CFLAGS) $(OFLAGS) $< -I $(IPATH) 

# dependencies
%.o: $(IPATH)%.h
