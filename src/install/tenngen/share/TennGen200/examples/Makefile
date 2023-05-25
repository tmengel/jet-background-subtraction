
################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libtenngen200.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -ltenngen200 -ldl
TENNGEN=$(PREFIX_LIB)/libtenngen200$(LIB_SUFFIX)
CXX_PYTHIA= $(PYTHIA8_INCLUDE) $(PYTHIA8_LIB) -lpythia8 -ldl

################################################################################
# RULES: Definition of the rules used to build the TENNGEN examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	$(info Usage: make tenngen-)

# TENNGEN library.
$(TENNGEN):
	$(error Error: TennGen must be built, please run "make"\
                in the top TennGen directory)
$(ROOT):
	$(error Error: TennGen requires root to be built)

# Examples without external dependencies.
tenngen-200: $(TENNGEN) $(AUAU200) tenngen200.cpp 
	$(CXX) tenngen200.cpp -o $@ $(CXX_COMMON) $(ROOT_LIB)\
	 `$(ROOT_CONFIG) --cflags --glibs`

tenngen-pythia: $(TENNGEN) tenngenpythia.cpp
ifeq ($(ROOT_USE)$(PYTHIA8_USE),truetrue)
	$(CXX) tenngenpythia.cpp -o $@ -w $(CXX_COMMON) $(CXX_PYTHIA) $(ROOT_LIB)\
	 `$(ROOT_CONFIG) --cflags --glibs`
else
	$(error Error: $@ requires ROOT and PYTHIA)
endif

tenngen-jets: $(TENNGEN) tenngenJets.cpp
ifeq ($(FASTJET3_USE)$(ROOT_USE)$(PYTHIA_USE),truetruetrue)
	$(CXX) $@.cpp -o $@ -w $(CXX_COMMON) $(FASTJET3_INCLUDE) $(FASTJET3_LIB)\
	 $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
else
	$(error Error: $@ requires FASTJET3, ROOT and PYTHIA)
endif


clean:
	@rm -f tenngen-*;\
	rm -f *.root;
	rm -rf output-*
