# /******************************************************************************
#  *   Copyright (C) 2014  Juan Antonio Garcia Martin , Peter Clote, Ivan Dotu  *
#  *                                                                            *
#  *  This program is free software: you can redistribute it and/or modify      *
#  *  it under the terms of the GNU General Public License as published by      *
#  *  the Free Software Foundation, either version 3 of the License, or         *
#  *  (at your option) any later version.                                       *
#  *                                                                            *
#  *  This program is distributed in the hope that it will be useful,           *
#  *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#  *  GNU General Public License for more details.                              *
#  *                                                                            *
#  *  You should have received a copy of the GNU General Public License         *
#  *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
#  ******************************************************************************/

ORTOOLS_DIR=$(HOME)/OR-Tools-src/or-tools-7.3
ORTOOLS_DEPENDENCIES_INCLUDES=$(ORTOOLS_DIR)/dependencies/install/include
ORTOOLS_SRC=$(ORTOOLS_DIR)/ortools
ORTOOLS_LIB=$(ORTOOLS_DIR)/lib
ORTOOLS_DEPENDENCIES_LIB=$(ORTOOLS_DIR)/dependencies/install/lib
ORTOOLS_SRC_GEN=$(ORTOOLS_SRC)/gen

VIENNA_FILES = data_structures.h  energy_const.h  fold_vars.h  libRNA.a  viennaC.cc  viennaC.h

############## Operating-System and Processor Architecture Detection #########
# This section attempts to auto-detect the OS. This is used only for convenience
# and allows the same Makefile to be used on multiple operating systems
# without modification.
# You can circumvent auto-detection by setting these environment variables:
#   RNA_MAKE_OS   -- the operating system name (e.g. Linux, Windows, Mac)
#   RNA_MAKE_ARCH -- the target architecture (e.g. x86 or x86_64)
# (These can be set as environment variables or on the MAKE command-line. 
#    e.g. `make all RNA_MAKE_OS=Linux`)
    ifneq (${RNA_MAKE_OS},) 
      #if RNA_MAKE_OS is NOT blank, use it as the OS
      OPSYSTEM=$(RNA_MAKE_OS)
    else ifeq (${OPSYSTEM},)
      # IF both RNA_MAKE_OS and OPSYSTEM are blank, use the `uname` command 
      #   (if available) to determine the OS.
      #   Replace 'UNKNOWN' with default OS if desired. 
      OPSYSTEM_RAW:=$(shell uname -s 2>/dev/null || echo UNKNOWN) 
      # Perform some replacements to normalize the output of uname on various systems.
      # OS_REPLACEMENTS= CYGWIN%=Windows MSYS%=Windows Darwin=Mac GNU%=Linux
      OPSYSTEM := $(OPSYSTEM_RAW:CYGWIN%=Windows)
      OPSYSTEM := $(OPSYSTEM:MSYS%=Windows)
      OPSYSTEM := $(OPSYSTEM:Darwin=Mac)
      OPSYSTEM := $(OPSYSTEM:GNU%=Linux)
      $(if $(DEBUG), $(info Make: Operating System: $(OPSYSTEM)))
      export OPSYSTEM #make it available for recursive calls to make, so auto-detection is performed only once.
    endif


CXX     	= g++
###CXXFLAGS        = -fPIC -std=c++0x -O4 -DNDEBUG -I$(ORTOOLS_SRC) -I$(ORTOOLS_DEPENDENCIES_INCLUDES) -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DFOLD_DEBUG
CXXFLAGS        = -fPIC -std=c++0x -O4 -DNDEBUG -I. -I$(ORTOOLS_DIR) -I$(ORTOOLS_SRC) -I$(ORTOOLS_DEPENDENCIES_INCLUDES) -I$(ORTOOLS_SRC_GEN) -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP 
    ifeq (${OPSYSTEM},Linux)
      ############ LINUX ##########################################################
      LDFLAGS         = -Wl,-rpath,$(ORTOOLS_DEPENDENCIES_LIB) -Wl,-rpath,$(ORTOOLS_LIB) -L$(ORTOOLS_LIB) -L$(ORTOOLS_DEPENDENCIES_LIB)  -Wl,-Bdynamic -lortools -lCbc -lCbcSolver -lClp -lCoinUtils -lOsiClp -lOsiCommonTests -lOsi -lCgl -lprotobuf -lgflags -lglog  -lz -lrt -lpthread -lstdc++ -lgcc -lm -L./ViennaRNA -l RNA -L./RNAstructure -l RNA -l HybridRNA 
    else ifeq (${OPSYSTEM},Mac)
      ############ MAC ############################################################
      LDFLAGS         = -lortools -lCbc -lCbcSolver -lClp -lCoinUtils -lOsiClp -lOsiCommonTests -lOsi -lCgl -lprotobuf -lgflags -lglog -Wl,-rpath,$(ORTOOLS_DEPENDENCIES_LIB) -Wl,-rpath $(ORTOOLS_LIB) -L$(ORTOOLS_LIB) -L$(ORTOOLS_DEPENDENCIES_LIB) -lz -lpthread -lstdc++ -lm -L./ViennaRNA -l RNA -L./RNAstructure -l RNA -l HybridRNA 
    else ifeq (,${OPSYSTEM})
      ############ NO OS ########################################################
	  $(error No Operating system defined!!!)
    else ifeq (${OPSYSTEM},Windows)
      ############ WINDOWS ########################################################
      LDFLAGS         = -Wl,-Bstatic -lconstraint_solver -lglop -llinear_solver -lalgorithms -lshortestpaths -lsat -lgraph -lbase -lCbc -lCbcSolver -lClp -lCoinUtils -lOsiClp -lOsiCommonTests -lOsi -lCgl -lprotobuf -lgflags -Wl,-rpath $(ORTOOLS_LIB) -L$(ORTOOLS_LIB) -L$(ORTOOLS_DEPENDENCIES_LIB) -L./ViennaRNA -l RNA -Wl,-Bdynamic -lz -lrt -lpthread -lstdc++ -lgcc -lm  -L./RNAstructure -l RNA -l HybridRNA 
   else
	  ############ UNKNOWN OS ###################################################
      $(error Unknown Operating system defined: $(OPSYSTEM))
    endif

all: checkDirectories RNAiFold

checkDirectories:
	@if [ ! -d "${ORTOOLS_DIR}" ]; then \
		echo "CANNOT BUILD RNAiFold: The OR-Tools directory ${ORTOOLS_DIR} does not exist. Check the ORTOOLS_DIR variable in Makefile."; \
		exit 1; \
	fi; \
	if [ ! -d "./RNAstructure" ]; then \
		echo "CANNOT BUILD RNAiFold: Compiled RNAstructure libraries must be copied to RNAstructure directory."; \
		exit 1; \
	fi; \
	if [ ! -d "./ViennaRNA" ]; then \
		echo "CANNOT BUILD RNAiFold: Compiled ViennaRNA source files must be copied to ViennaRNA directory."; \
		exit 1; \
	fi; \
	if [ ! -f "./ViennaRNA/libRNA.a" ]; then \
		echo "CANNOT BUILD RNAiFold: Compiled ViennaRNA source files must be copied to ViennaRNA directory."; \
		exit 1; \
	fi;

RNAiFold: mfe.o misc.o ifold.o cphelix.o strtree.o strtreegroup.o value_heuristic.o  vienna_constraint.o vienna_plugin.o rnastructure_plugin.o rna_plugin.o restart_monitor.o diversity_measures.o minenergy_constraint.o minensdef_constraint.o helix_constraint.o local_constraint.o energy_constraint.o ensdef_constraint.o aaconstraint.o
	$(CXX) $(CXXFLAGS) mfe.o misc.o ifold.o cphelix.o strtree.o strtreegroup.o value_heuristic.o  vienna_constraint.o vienna_plugin.o rnastructure_plugin.o rna_plugin.o restart_monitor.o diversity_measures.o minenergy_constraint.o minensdef_constraint.o helix_constraint.o local_constraint.o energy_constraint.o ensdef_constraint.o aaconstraint.o $(LDFLAGS) -o RNAiFold

ifold.o: ifold.cc
	$(CXX) $(CXXFLAGS) -c ifold.cc -o ifold.o

mfe.o: mfe.cc
	$(CXX) $(CXXFLAGS) -c mfe.cc -o mfe.o

strtree.o: strtree.cc
	$(CXX) $(CXXFLAGS) -c strtree.cc -o strtree.o

cphelix.o: cphelix.cc
	$(CXX) $(CXXFLAGS) -c cphelix.cc -o cphelix.o

strtreegroup.o: strtreegroup.cc
	$(CXX) $(CXXFLAGS) -c strtreegroup.cc -o strtreegroup.o

value_heuristic.o: value_heuristic.cc
	$(CXX) $(CXXFLAGS) -c value_heuristic.cc -o value_heuristic.o

misc.o: misc.c
	$(CXX) $(CXXFLAGS) -c misc.c -o misc.o

vienna_constraint.o: vienna_constraint.cc
	$(CXX) $(CXXFLAGS) -c vienna_constraint.cc -o vienna_constraint.o

vienna_plugin.o: vienna_plugin.cc
	$(CXX) $(CXXFLAGS) -c vienna_plugin.cc -o vienna_plugin.o

rnastructure_plugin.o: rnastructure_plugin.cc
	$(CXX) $(CXXFLAGS) -c rnastructure_plugin.cc -o rnastructure_plugin.o

rna_plugin.o: rna_plugin.cc
	$(CXX) $(CXXFLAGS) -c rna_plugin.cc -o rna_plugin.o

restart_monitor.o: restart_monitor.cc
	$(CXX) $(CXXFLAGS) -c restart_monitor.cc -o restart_monitor.o

diversity_measures.o: diversity_measures.cc
	$(CXX) $(CXXFLAGS) -c diversity_measures.cc -o diversity_measures.o

minenergy_constraint.o: minenergy_constraint.cc
	$(CXX) $(CXXFLAGS) -c minenergy_constraint.cc -o minenergy_constraint.o

minensdef_constraint.o: minensdef_constraint.cc
	$(CXX) $(CXXFLAGS) -c minensdef_constraint.cc -o minensdef_constraint.o

helix_constraint.o: helix_constraint.cc
	$(CXX) $(CXXFLAGS) -c helix_constraint.cc -o helix_constraint.o

local_constraint.o: local_constraint.cc
	$(CXX) $(CXXFLAGS) -c local_constraint.cc -o local_constraint.o

energy_constraint.o: energy_constraint.cc
	$(CXX) $(CXXFLAGS) -c energy_constraint.cc -o energy_constraint.o

ensdef_constraint.o: ensdef_constraint.cc
	$(CXX) $(CXXFLAGS) -c ensdef_constraint.cc -o ensdef_constraint.o

aaconstraint.o: aaconstraint.cc
	$(CXX) $(CXXFLAGS) -c aaconstraint.cc -o aaconstraint.o

clean:
	rm -rf RNAiFold *.o 

