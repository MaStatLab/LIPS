include Makefile.in

C = main.cc rng.cc history.cc simfunctions.cc
H = smctc.hh sampler.hh particle.hh moveset.hh history.hh rng.hh

PFC = modelsel_pf.cc ## modelsel_pf_funcs.cc
PFO = LIPS.so ## modelsel_pf_funcs.o
PFH = modelsel_pf_funcs.hh

CCFLAGS := $(CCFLAGS) -Ismctc/include -Lsmctc/lib
LFLAGS := -Ismctc/include -Lsmctc/lib $(LFLAGS)

all: smctc pf_R

.PHONY: pf_R smctc clean

pf_R: $(PFC) $(PFH)
	R CMD SHLIB -o $(PFO) $(PFC) $(CCFLAGS) -lsmctc -lgsl -lgslcblas

smctc:
	+$(MAKE) -Csmctc

clean:
	-rm *.o *.so
	-rm *~

	make -Csmctc clean
