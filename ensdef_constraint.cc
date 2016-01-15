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

#include "ensdef_constraint.h"
#include "diversity_measures.h"
#include "misc.h"
#include <cstring>

namespace operations_research {
	void EnsDefConstraint::InitialPropagate() {
	    int allBound =1;
		vector<int64> freeVars;
		for (int i = 1; i < size(); ++i) {
			if (!vars_[i]->Bound()) {
				allBound=0;
				return;
			}
		}
		if(allBound){
			char seq[n+1];
			double ens_defect=0;

			for(int i = 0; i < size(); ++i) {
				seq[i]=ToNucl(vars_[i]->Value());
			}
			seq[n]='\0';

			double **outBpPr;
			v->setSequence(seq);

			if(cutPoint_==-1){
				outBpPr=v->basePairProbsMatrix();
			}
			else{
				outBpPr=v->basePairProbsMatrixCofold();					
			}
			
			ens_defect=ensembleDefect(outBpPr,n,bpList_);
			
//			cout << seq <<endl;
//			for(int i=0; i< n; i++){
//				cout << bpList_[i] << " ";
//			}
//			cout << endl;
//			cout << "Ensemble defect constraint: " << ens_defect << " (Limits: " << min_ensdef_ <<" - " <<max_ensdef_ << ")"<< endl;
			
			if(max_ensdef_!= NO_ED_LIMIT && max_ensdef_<ens_defect){
				solver()->Fail();
				free(outBpPr);
				return;
			}
			if(min_ensdef_!= NO_ED_LIMIT && min_ensdef_>ens_defect){
				solver()->Fail();
				free(outBpPr);
				return;
			}

			free(outBpPr);
			return;
		}
	}
}
