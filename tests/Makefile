
all: force integration io kinetic

force:
	#$(MAKE) -C force
	$(MAKE) -C force run
integration:
	#$(MAKE) -C integration
	$(MAKE) -C integration run

kinetic:
	#$(MAKE) -C kinetic
	$(MAKE) -C kinetic run
io:
	#$(MAKE) -C io
	$(MAKE) -C io run
.PHONY: force integration kinetic io clean

clean:
	$(MAKE) $(MFLAGS) -C force clean
	$(MAKE) $(MFLAGS) -C integration clean
	$(MAKE) $(MFLAGS) -C kinetic clean
	$(MAKE) $(MFLAGS) -C io clean
