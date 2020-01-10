#
# Makefile for Neutrino Generator Test Programs
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include
INCLUDES += -I$(GENIE)/include/GENIE/

INCLUDES += -I$(MYGENIE)

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
MYGENIE_LIBS := -L$(MYGENIE)/lib -lGenieDriversIceCube -lHelperClasses

#LIBRARIES  := $(LHAPDF_LIBRARIES) $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)  $(MYGENIE_LIBS)


TGT =	hepevt_writer \
	gicevgen

all: $(TGT)

hepevt_writer: FORCE
	$(CXX) $(CXXFLAGS) -c hepevt_writer.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) hepevt_writer.o $(LIBRARIES) -o hepevt_writer

gicevgen: FORCE
	$(CXX) $(CXXFLAGS) -c gICEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gICEvGen.o $(LIBRARIES) -o gicevgen


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) gicevgen
	$(RM) hepevt_writer



FORCE:

# DO NOT DELETE
