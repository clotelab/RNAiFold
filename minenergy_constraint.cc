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

#include "minenergy_constraint.h"
#include "misc.h"
#include <cstring>

namespace operations_research {
	void MinEnergyConstraint::InitialPropagate() {
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
			if(trgStr_==""){
				if(cutPoint_==-1){
					energy=v->fold();
				}
				else{
					energy=v->cofold();
				}
			}
			else{
				energy=v->energyOfStruct();
			}
			
			//cout << "***************************Setting maximum to "<<energy << " in max: " <<vMaxMFE->Max()<< " - min: " <<vMaxMFE->Min() << endl;

			vMaxMFE->SetValue((energy-0.005)*100);
			//cout << "***************************Set maximum to "<<energy << " in max: " <<vMaxMFE->Max()<< " - min: " <<vMaxMFE->Min() << endl;
			return;
		}
	}
}
