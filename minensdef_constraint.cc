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

#include "minensdef_constraint.h"
#include "diversity_measures.h"
#include "misc.h"
#include <cstring>

namespace operations_research {
	void MinEnsDefConstraint::InitialPropagate() {
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
			int freeBP =0;
			if(bpList_==NULL){
				freeBP=1;
				if(cutPoint_==-1){
					v->fold();
				}
				else{
					v->cofold();					
				}
				bpList_=v->getBasePairs();
			}

			if(cutPoint_==-1){
				outBpPr=v->basePairProbsMatrix();
			}
			else{
				outBpPr=v->basePairProbsMatrixCofold();					
			}
			
			ens_defect=ensembleDefect(outBpPr,n,bpList_);
			
			if(freeBP==1){
				free(bpList_);
				bpList_=NULL;
			}
			
//			cout << "Ensemble defect constraint: " << ens_defect << " (MAX: " << maxEnsDef << ")"<< endl;
			vMaxEnsDef->SetValue((ens_defect+(0.5/ED_PRECISION))*ED_PRECISION);


			free(outBpPr);
			return;
		}
	}
}
