/******************************************************************************
 *   Copyright (C) 2014  Juan Antonio Garcia Martin , Peter Clote, Ivan Dotu  *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 ******************************************************************************/


#include <vector>
#include "strtree.h"

#ifndef _STR_TREE_GROUP
#define _STR_TREE_GROUP
#include <vector>
#include <map>
#include "constraint_solver/constraint_solver.h"
#include "util/string_array.h"

using namespace std;

namespace operations_research {
	
	// Prints the instance of the Frequency Assignment Problem
	class StrTreeGroup{
		public:			
			StrTreeGroup(vector<StrTree*> *str_trees);
			int getStrId(int index);
			int getHelixPos(int index);
			CPHelix* getHelix(int index);
			int getNumHelices();
			// void optimizeValueHeuristic(int strategy);
			
			std::vector<std::vector<int>> optimizeValueHeuristic(int n, std::vector<std::vector<int> > BPO,std::vector<std::vector<int> > BPC,std::vector<int*> str_int,std::vector<double> temps);
			std::vector<std::vector<int>> optimizeUPHeuristic(int n,std::vector<std::vector<int> > BPO,std::vector<std::vector<int> > BPC,std::vector<int*> str_int);
			vector<int64> optimizeOrder(int strategy);
			vector<size_t> getOrderIndex();
			~StrTreeGroup();
		private:
			int nHelices;
			vector<StrTree*>* str_trees_;
			vector<int>  helixIndex;
			vector<int>  posIndex;
			vector<int64> order_;
			vector<size_t> orderIndex;
	};

	
	
}
#endif
