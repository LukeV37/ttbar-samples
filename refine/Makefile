CPP=g++
LD=g++
# the flag -D__USE_XOPEN2K8 is needed for slc6 32 bit version, see
# https://twiki.cern.ch/twiki/bin/view/Atlas/RPMCompatSLC6
#SPECIALFLAGS=-O2 -D__USE_XOPEN2K8
SPECIALFLAGS=-O2
#SPECIALFLAGS=-g -O0
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

CFLAGS = $(SPECIALFLAGS) -I. -I./RootCore/include
LFLAGS = -L. -L./RootCore/lib

RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC=d_ana.cpp d.cpp

%.o: %.cpp d.h 
	@$(CPP) $(RCXX) -c $<

all: $(SRC:.cpp=.o)
	@$(LD) $(SRC:.cpp=.o) $(RLXX) -o d_ana

clean:
	@rm -f *~ *.o
	@rm -f d_ana
	@rm -f a1.root
