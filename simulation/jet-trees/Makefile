
-include Makefile.inc


BINS = $(patsubst %.C, %, $(wildcard *.C))
		


CXX = g++
CCFLAGS = $(shell root-config --cflags) #-ggdb 

CXX_PYTHIA= $(PYTHIA8_INCLUDE) $(PYTHIA8_LIB) #-lpythia8 -ldl
CXX_COMMON:=$(CXX_COMMON) $(GZIP_LIB)
# CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -ltenngen200 -ldl


LD = g++
LDFLAGS = 

LIBS = $(shell root-config --libs) 
CXX_ROOT = $(ROOT_LIB) $(shell root-config --cflags --glibs --libs) 
default : $(BINS)
# $(FASTJET3_INCLUDE) $(FASTJET3_LIB)
$(BINS): % : %.C 
	@echo -n "Building $@ ... "
	$(CXX)  $(CCFLAGS) $<  -o $@ -w  $(CXX_COMMON) $(shell root-config --cflags --glibs --libs) $(FASTJET3_INCLUDE) $(FASTJET3_LIB) $(shell $(FASTJET3_CONFIG)/fastjet-config --cxxflags --libs)
	@echo "Done"


clean: 
	rm -f *.o	$(BINS)

jobclean:
	rm *.o* *.e*
