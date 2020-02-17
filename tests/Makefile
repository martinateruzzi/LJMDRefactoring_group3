
all: force integration io kinetic

force:
	$(MAKE) -C force

integration:
	$(MAKE) -C integration

kinetic:
	$(MAKE) -C kinetic

io:
	$(MAKE) -C io

.PHONY: force integration kinetic io clean

clean:
	$(MAKE) $(MFLAGS) -C force
	$(MAKE) $(MFLAGS) -C integration
	$(MAKE) $(MFLAGS) -C kinetic
	$(MAKE) $(MFLAGS) -C io
