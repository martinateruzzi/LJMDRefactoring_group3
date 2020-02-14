
all: force integration io kinetic
	(cd force; make all)
	(cd integration; make all)
	(cd io; make all)
	(cd Parallel-OPenmp; make all)

force:
	$(MAKE) -C force

integration:
	$(MAKE) -C integration

kinetic:
	$(MAKE) -C kinetic

io:
	$(MAKE) -C io

.PHONY: force integration kinetic io
