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

#include "cphelix.h"
#include <vector>

#ifndef _STR_TREE
#define _STR_TREE

#define TH_RNAIFOLD 1
#define TH_RNAIFOLD_DANGLES 2

#define MAX_LEFT_MISMATCH 2
#define MAX_RIGHT_MISMATCH 2

namespace operations_research {
	
	// Prints the instance of the Frequency Assignment Problem
	class StrTree{
		public:
			~StrTree();
			StrTree(int tree_id, int* str_int,int n, std::vector<int> BPO, std::vector<int> BPC, int treeType, int cutPoint);
			StrTree(int* str_int,int n, std::vector<int> BPO, std::vector<int> BPC, int treeType, int cutPoint);

			void createRNAiFoldHelices(int* str_int,int n, std::vector<int> BPO, std::vector<int> BPC, bool includeDangles, int cutPoint);

			int getTreeId();
			void setTreeId(int tree_id);
			
			std::vector<CPHelix*> getHelices();
			std::vector<CPHelix* > getSortedHelices();
			void showTree();

		private:
			int consecutive(int aX, int aY, int bX, int bY);
			std::vector<CPHelix*> helices;
			int tree_id_;
	};

	inline bool compByLevel(CPHelix* a, CPHelix* b){
		return a->getLevel() > b->getLevel();
	}	
	
	
}
#endif
