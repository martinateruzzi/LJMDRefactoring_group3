SHELL=/bin/bash

TST := $(wildcard *.c)
OBJ := $(TST:.c=.o)
EXE := $(TST:.c=.x)

SDIR := ../../src
IDIR := ../../inc
ODIR := ../../obj

CC := mpicc
CFLAGS := -Wall -g -O3 -ffast-math -fomit-frame-pointer -fopenmp -I $(IDIR) -D_MPI
LDLIBS := -lm -lljmd -fopenmp
LDFLAGS :=  -L $(ODIR)

all: $(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.x: %.o 
	$(CC) $< -o $@ $(LDFLAGS) $(LDLIBS) 

.PHONY: clean run

run: $(EXE)
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(ODIR)
	./$(EXE) < argon_003.inp

clean:
	rm *.x *.frc 
