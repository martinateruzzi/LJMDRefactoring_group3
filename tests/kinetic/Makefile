SHELL=/bin/bash

TST := $(wildcard *.c)
OBJ := $(TST:.c=.o)
EXE := $(TST:.c=.x)

SDIR := ../../src
IDIR := ../../inc
ODIR := ../../obj

LD_PATH := ${LD_LIBRARY_PATH}

CFLAGS := -Wall -g -O3 -ffast-math -fomit-frame-pointer -DOPENMP -fopenmp  -I $(IDIR)
CC := mpicc
LDLIBS := -lm -lljmd -fopenmp
LDFLAGS :=  -L $(ODIR)

all: $(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.x: %.o 
	$(CC) $< -o $@ $(LDFLAGS) $(LDLIBS) 

.PHONY: clean
clean:
	rm *.x *.xyz

run: test_force.x
	@$(eval export LD_LIBRARY_PATH=$(LD_PATH):$(ODIR))
	./test_kinetic.x < argon_003.inp
	cmp argon_003.dat argon_003_serial.reference || exit 1
	rm -f argon_003.dat     

test: 
	$(eval export LD_LIBRARY_PATH=$(LD_PATH):$(ODIR))
	@echo ${LD_LIBRARY_PATH}
	@echo "LD_PATH" $(LD_PATH)



