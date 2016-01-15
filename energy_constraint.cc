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

#include "energy_constraint.h"
#include "misc.h"
#include <cstring>

namespace operations_research {
	void EnergyConstraint::InitialPropagate() {
	    int allBound =1;
		vector<int64> freeVars;
		for (int i = 0; i < size(); ++i) {
			if (!vars_[i]->Bound()) {
				allBound=0;
				return;
			}
		}
		if(allBound){
			char seq[n+1];
			double energy=0;

			for(int i = 0; i < size(); ++i) {
				seq[i]=ToNucl(vars_[i]->Value());
			}
			seq[n]='\0';

			v->setSequence(seq);

			energy=v->energyOfStruct();
			
//			cout << seq << endl<<trgStr_<<endl <<"Energy "<< energy << endl;
			if(max_energy_!= NO_ENERGY_LIMIT && (float)max_energy_<(float)energy){
				solver()->Fail();
			}
			if(min_energy_!= NO_ENERGY_LIMIT && (float)min_energy_>(float)energy){
				solver()->Fail();
			}
			return;
		}
	}
}
