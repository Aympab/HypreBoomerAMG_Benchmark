# Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

########################################################################
# Compiler and external dependences
########################################################################
CC        = mpicc
F77       = mpif77
CXX       = mpicxx
F90       = mpifort
#HYPRE_DIR = ../hypre
#HYPRE_DIR = /home/christophe/Works/Projects/hypre/src/hypre
#HYPRE_DIR = /gext/christophe.bovet/hypre/src/hypre
HYPRE_DIR = bin_hypre/hypre_cuda/hypre
#HYPRE_DIR = bin_hypre/hypre/hypre
########################################################################
# Compiling and linking options
########################################################################
COPTS     = -g -Wall
CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)
FOPTS     = -g
FINCLUDES = $(CINCLUDES)
FFLAGS    = $(FOPTS) $(FINCLUDES)
CXXOPTS   = $(COPTS) -Wno-deprecated
CXXINCLUDES = $(CINCLUDES) -I..
CXXDEFS   = $(CDEFS)
IFLAGS_BXX = 
CXXFLAGS  = $(CXXOPTS) $(CXXINCLUDES) $(CXXDEFS) $(IFLAGS_BXX)
IF90FLAGS = 
F90FLAGS = $(FFLAGS) $(IF90FLAGS)

LINKOPTS  = $(COPTS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++ 
LFLAGS_B =\
 -L${HYPRE_DIR}/lib\
 -lbHYPREClient-C\
 -lbHYPREClient-CX\
 -lbHYPREClient-F\
 -lbHYPRE\
 -lsidl -ldl -lxml2
LFLAGS77 = $(LFLAGS)
LFLAGS90 =

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .f .cxx .f90

.c.o:
	$(CC) $(CFLAGS) -c  $< 
.f.o:
	$(F77) $(FFLAGS) -c $<
.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

########################################################################
# List of all programs to be compiled
########################################################################
ALLPROGS = tp original

all: $(ALLPROGS)

default: all


########################################################################
# Example 5
########################################################################
tp: tp.o
	$(CC) -o $@ $^ $(LFLAGS)

original: original.o
	$(CC) -o $@ $^ $(LFLAGS)


########################################################################
# Clean up
########################################################################
clean:
	rm -f $(ALLPROGS:=.o)
	cd vis; make clean
distclean: clean
	rm -f $(ALLPROGS) $(ALLPROGS:=*~)
	rm -fr README*
