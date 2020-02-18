SHELL=/bin/bash

TST := $(wildcard *.c)
OBJ := $(TST:.c=.o)
EXE := $(TST:.c=.x)

SDIR := ../../src
IDIR := ../../inc
ODIR := ../../obj

LD_PATH := ${LD_LIBRARY_PATH}


CC := gcc
CFLAGS := -Wall -g -O3 -ffast-math -fomit-frame-pointer -DOPENMP -fopenmp  -I $(IDIR) 
LDLIBS := -lm -lljmd
LDFLAGS :=  -L $(ODIR)

all: $(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.x: %.o 
	$(CC) $< -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY: clean

clean:
	rm -f *.x *.xyz

run: test_integration.x
	@$(eval export LD_LIBRARY_PATH=$(LD_PATH):$(ODIR))
	./test_integration.x < argon_003.inp
	cmp argon_003.dat argon_003_serial.reference || exit 1
	rm -f argon_003.dat	

test: 
	$(eval export LD_LIBRARY_PATH=$(LD_PATH):$(ODIR))
	@echo ${LD_LIBRARY_PATH}
	@echo "LD_PATH" $(LD_PATH)





