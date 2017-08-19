.PHONY: docs clean distclean examples

default: libraries

all: libraries

clean:
	make -Csrc clean
	make -Cexamples clean
	-rm *~
	-rm */*~

distclean: clean
	-rm bin/*
	-rm lib/*
	-rm -Rf doc/*

libraries:
	make -Csrc all


dist: distclean
	tar -cf ../smctc.tar *
	bzip2 ../smctc.tar
	tar -cf ../smctc.tar *
	zip -r ../smctc.zip *
