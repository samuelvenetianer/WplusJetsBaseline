# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2020 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all: MyPythia8Simul


# PYTHIA library.
$(PYTHIA):
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

MyPythia8Simul: $(PYTHIA) MyPythia8Simul.o ANA_utils.o TruthPart.o TruthJets.o 
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(ROOT_USE),truetruetrue)
	$(CXX) MyPythia8Simul.o ANA_utils.o TruthPart.o TruthJets.o -o $@ -w $(FASTJET3_INCLUDE) $(HEPMC2_INCLUDE) $(CXX_COMMON)\
        $(FASTJET3_LIB) -lfastjettools $(HEPMC2_LIB) $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
#MyPythia8Simul: $(PYTHIA) MyPythia8Simul.o ANA_utils.o TruthPart.o TruthJets.o 
#ifeq ($(FASTJET3_USE)$(ROOT_USE),truetrue)
#	$(CXX) MyPythia8Simul.o ANA_utils.o TruthPart.o TruthJets.o -o $@ -w $(FASTJET3_INCLUDE) $(CXX_COMMON)\
#        $(FASTJET3_LIB) -lfastjettools $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
else
        @echo "Error: $@ requires ROOT and FASTJET3"
endif



MyPythia8Simul.o: $(PYTHIA) MyPythia8Simul.cc MyPythia8Simul.h ANA_utils.h StdArg.hpp TruthPart.h TruthJets.h 
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(ROOT_USE),truetruetrue)
	${CXX} -c ${CXXFLAGS} MyPythia8Simul.cc -w $(HEPMC2_INCLUDE) $(FASTJET3_INCLUDE) $(CXX_COMMON)\
        $(FASTJET3_LIB) -lfastjettools $(HEPMC2_LIB) $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
#MyPythia8Simul.o: $(PYTHIA) MyPythia8Simul.cc MyPythia8Simul.h ANA_utils.h StdArg.hpp TruthPart.h TruthJets.h 
#ifeq ($(FASTJET3_USE)$(ROOT_USE),truetrue)
#	${CXX} -c ${CXXFLAGS} MyPythia8Simul.cc -w $(FASTJET3_INCLUDE) $(CXX_COMMON)\
#        $(FASTJET3_LIB) -lfastjettools $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
else
        @echo "Error: $@ requires ROOT and FASTJET3"
endif



ANA_utils.o: ANA_utils.cc ANA_utils.h
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(ROOT_USE),truetruetrue)
	${CXX} -c ${CXXFLAGS} ANA_utils.cc -w $(HEPMC2_INCLUDE) $(FASTJET3_INCLUDE) $(CXX_COMMON)\
        $(FASTJET3_LIB) -lfastjettools $(HEPMC2_LIB) $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
#ANA_utils.o: ANA_utils.cc ANA_utils.h
#ifeq ($(FASTJET3_USE)$(ROOT_USE),truetrue)
#	${CXX} -c ${CXXFLAGS} ANA_utils.cc -w $(FASTJET3_INCLUDE) $(CXX_COMMON)\
#        $(FASTJET3_LIB) -lfastjettools $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
else
        @echo "Error: $@ requires ROOT, FASTJET3 and HEPMC2"
endif


#   Compile the physics objets classes
#   ----------------------------------

TruthPart.o:  TruthPart.cc TruthPart.h 
	${CXX} -c ${CXXFLAGS} TruthPart.cc $(CXX_COMMON)

TruthJets.o:  TruthJets.cc TruthJets.h 
	${CXX} -c ${CXXFLAGS} TruthJets.cc $(CXX_COMMON)




# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat;\
	rm -f weakbosons.lhe; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;\
	rm -f *.log; rm -f *plot.py; rm -f *.pcm; \
	rm -f *.o;
