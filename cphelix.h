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
#include "energy_constant.h"

#ifndef _CP_HELIX
#define _CP_HELIX
namespace operations_research {	
	class CPHelix{
		public:
			~CPHelix();
			CPHelix(int _id, int _i, int _j);

			int getId();
			void setId(int _id);
			int getI();
			void setI(int _i);
			int getJ();
			void setJ(int _j);
			int getParent();
			void setParent(int _parent);
			int getLevel();
			void setLevel(int _level);

			int getPosLeft();
			void setPosLeft(int npos);
			int getPosRight();
			void setPosRight(int npos);

			
			void addInside(int h);
			void removeInside(int h);
			std::vector<int> getInside();
			
			int contains(CPHelix* h);
			void setHelixPositions(std::vector<CPHelix*> helices);
			
			std::vector<int> getPositions();
			void addPositions(std::vector<int> newPositions);

			int getPositionsIntersect(std::vector<int> newPositions);
						
			void adaptToRemoveOf(int h);
			int getStackingEnergy(int* str_int, int* hairpin37, int bp_energy);

			void setPositionIndexAndType(int* str_int, std::vector<int> BPO, std::vector<int> BPC);
			
			std::vector<int> getBPs();
			void addBPs(std::vector<int> newBPs);

			
			std::vector<int> getClosingBPs();
			void addClosingBPs(std::vector<int> newCBPs);
					
			std::vector<int> getBPpositions();
			int getBPintersect(std::vector<int> newBPpositions);

			std::vector<int> getUPs();	
			void addUPs(std::vector<int> newUPs);
			std::vector<int> getOrderedUPs(int* str_int);
		
			bool isStack(int bp);
						
			void toString();


		private:
			int id;
			int i;
			int j;
			int parent;
			int level;
			int pos_left;
			int pos_right;			
			std::vector<int> inside;
			std::vector<int> positions;
			std::vector<int> bps;
			std::vector<int> cbps;
			std::vector<int> ups;
			std::vector<int> stackBPs;
			std::vector<int> bpPositions;

	};
}
#endif
