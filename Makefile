# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=obj/%.o)
############################################

default: serial

serial: 
	$(MAKE) $(MFLAGS) -C obj

clean:
	$(MAKE) $(MFLAGS) -C obj clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C tests clean

check: serial
	$(MAKE) $(MFLAGS) -C examples check

test : default
	$(MAKE) -C tests

.PHONY: test default
