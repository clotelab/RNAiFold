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

ORTOOLS_DIR=$(HOME)/OR-Tools-src/or-tools-read-only
ORTOOLS_DEPENDENCIES_INCLUDES=$(ORTOOLS_DIR)/dependencies/install/include
ORTOOLS_SRC=$(ORTOOLS_DIR)/src
ORTOOLS_LIB=$(ORTOOLS_DIR)/lib
ORTOOLS_DEPENDENCIES_LIB=$(ORTOOLS_DIR)/dependencies/install/lib

VIENNA_FILES = data_structures.h  energy_const.h  fold_vars.h  libRNA.a  viennaC.cc  viennaC.h

CXX     	= g++
CXXFLAGS        = -fPIC -std=c++0x -O4 -DNDEBUG -I$(ORTOOLS_SRC) -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP
LDFLAGS         = -Wl,-Bstatic -lconstraint_solver -lglop -llinear_solver -lalgorithms -lshortestpaths -lsat -lgraph -lbase -lutil -lCbc -lCbcSolver -lClp -lCoinUtils -lOsiClp -lOsiCommonTests -lOsi -lCgl -lprotobuf -lgflags -Wl,-rpath $(ORTOOLS_LIB) -L$(ORTOOLS_LIB) -L$(ORTOOLS_DEPENDENCIES_LIB) -L./ViennaRNA -l RNA -Wl,-Bdynamic -lz -lrt -lpthread -lstdc++ -lgcc -lm  -L./RNAstructure -l RNA -l HybridRNA 

all: RNAiFold

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

