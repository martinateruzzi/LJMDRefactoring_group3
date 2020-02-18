
all: force integration io kinetic

force:
	$(MAKE) -C force

integration:
	#$(MAKE) -C integration
	$(MAKE) -C integration run

kinetic:
	$(MAKE) -C kinetic

io:
	$(MAKE) -C io

.PHONY: force integration kinetic io clean

clean:
	$(MAKE) $(MFLAGS) -C force clean
	$(MAKE) $(MFLAGS) -C integration clean
	$(MAKE) $(MFLAGS) -C kinetic clean
	$(MAKE) $(MFLAGS) -C io clean
