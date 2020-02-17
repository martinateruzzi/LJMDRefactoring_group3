
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
	$(MAKE) clean -C force
	$(MAKE) clean -C integration
	$(MAKE) clean -C kinetic
	$(MAKE) clean -C io
