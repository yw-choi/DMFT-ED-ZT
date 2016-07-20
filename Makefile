.PHONY: all clean

default: all

all:
	(cd src; $(MAKE); cd ..; \
	 if [ ! -d bin ]; then mkdir bin; fi;\
	 cp src/dmft-ed.x bin/;\
	 cp utils/spectral_ftn.py bin/)	
clean:
	(cd src; $(MAKE) clean)
