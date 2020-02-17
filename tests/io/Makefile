SHELL=/bin/bash

TST := $(wildcard *.c)
OBJ := $(TST:.c=.o)
EXE := $(TST:.c=.x)

SDIR := ../../src
IDIR := ../../inc
ODIR := ../../obj

CC := gcc
CFLAGS := -Wall -g -O3 -ffast-math -fomit-frame-pointer  -I $(IDIR)
LDLIBS := -lm -lljmd
LDFLAGS :=  -L $(ODIR)

all: $(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.x: %.o 
	$(CC) $(LDFLAGS) $(LDLIBS) $< -o $@

.PHONY: clean

clean:
	rm *.x
