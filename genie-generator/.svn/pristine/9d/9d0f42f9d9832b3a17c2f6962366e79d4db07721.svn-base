#
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include
INCLUDES += -I$(GENIE)/include/GENIE/



PACKAGE           = HelperClasses
DICTIONARY        = _ROOT_DICT_HelperClasses
LIBNAME           = libHelperClasses
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
#install : install-inc install-lib

# Include standard package makefile targets
#
#include $(GENIE)/src/make/Make.std-package-targets

# this is a copy of what is in the previous file to deal with the fact
# we are working outside GENIE src code
#
PACKAGE_LIB_PATH              = $(MYGENIE)/lib
PACKAGE_LIB                   = $(LIBNAME)$(GLIBVERSION).$(DllSuf)
PACKAGE_LIB_WITH_PATH         = $(PACKAGE_LIB_PATH)/$(LIBNAME)$(GLIBVERSION).$(DllSuf)
PACKAGE_LIBNOVRS_WITH_PATH    = $(PACKAGE_LIB_PATH)/$(LIBNAME).$(DllSuf)
PACKAGE_LIBMACLINK_WITH_PATH  = $(PACKAGE_LIB_PATH)/$(LIBNAME).$(DllLinkSuf)
# PACKAGE_LIBNOVRS_WITH_IPATH   = $(PACKAGE_LIB_PATH)/$(LIBNAME).$(DllSuf)
# PACKAGE_LIBMACLINK_WITH_IPATH = $(PACKAGE_LIB_PATH)/$(LIBNAME).$(DllLinkSuf)

PACKAGE_ALL_SOURCES          := $(wildcard *.cxx)
DICTGEN_HEADERS              := $(addsuffix .h, $(basename $(PACKAGE_ALL_SOURCES)))
DICTGEN_OBJECTS              := $(addsuffix .$(ObjSuf), $(basename $(DICTGEN_HEADERS)))
PACKAGE_ALL_SOURCES          := $(PACKAGE_ALL_SOURCES) $(DICTIONARY).cc $(wildcard *.cpp)
PACKAGE_ALL_OBJECTS           = $(addsuffix .$(ObjSuf), $(basename $(PACKAGE_ALL_SOURCES)))
PACKAGE_SOURCES              := $(wildcard *.cxx)
PACKAGE_OBJECTS               = $(addsuffix .$(ObjSuf), $(basename $(PACKAGE_SOURCES)))

LIB_DEPEND=
ifeq ($(strip $(DICTIONARY)),undefined)
  LIB_DEPEND = $(PACKAGE_OBJECTS) 
else
  LIB_DEPEND = $(PACKAGE_ALL_OBJECTS)
endif

rootcint: $(DICTIONARY).o

$(DICTIONARY).cc: LinkDef.h $(DICTGEN_OBJECTS)
	$(RM) $(DICTIONARY).*
	@echo "[package: $(PACKAGE)] Generating ROOT dictionary ..."
	rootcint -f $(DICTIONARY).cc -c -I$(ROOTSYS)/include -I$(GENIE)/src -I$(GENIE)/src/NuValidator -I./ $(DICTGEN_HEADERS) LinkDef.h

lib: $(LIB_DEPEND)
	@echo "[package: $(PACKAGE)] Generating shared library ..."
ifeq ($(strip $(DICTIONARY)),undefined)
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PACKAGE_OBJECTS) $(EXTRA_EXT_LIBS) $(OutPutOpt) $(PACKAGE_LIB_WITH_PATH)
else
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PACKAGE_ALL_OBJECTS) $(EXTRA_EXT_LIBS) $(OutPutOpt) $(PACKAGE_LIB_WITH_PATH)
endif
#
#
lib-link: FORCE
	@echo "[package: $(PACKAGE)] Creating symbolic link to shared library ..."
ifeq ($(strip $(GOPT_ENABLE_DYLIBVERSION)),YES)
	$(SYMLINK) $(PACKAGE_LIB) $(PACKAGE_LIBNOVRS_WITH_PATH)
endif
ifeq ($(strip $(ARCH)),macosx)
	$(SYMLINK) $(PACKAGE_LIB) $(PACKAGE_LIBMACLINK_WITH_PATH)
endif


#
#
purge: FORCE
	$(RM) $(PACKAGE_OBJECTS) *~ core a.out Makefile.temp
 

#
#
clean: FORCE   
	$(RM) *.o *~ core a.out $(PACKAGE_LIB_PATH)/$(LIBNAME)* $(DICTIONARY).* Makefile.temp




FORCE:

