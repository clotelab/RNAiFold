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

#include <string>
#include <vector>
#include <map>
#include <limits.h>
#include <ctime>


#include "base/logging.h"
#include "util/string_array.h"
#include "constraint_solver/constraint_solver.h"
#include "constraint_solver/constraint_solveri.h"

#include "rna_plugin.h"


namespace operations_research {
// This constraint maintains:
// 
// 
class ViennaConstraint : public Constraint {
	public:

	// This constructor does not take any ownership on its arguments.
	ViennaConstraint(Solver* const s, const std::vector<IntVar*>& vars, const std::vector<int64>& structure, int dangles, std::string rnaLib, std::string energyModel, double foldTemp, int cutPoint, int posLeft, int posRight, double upperMFE, double lowerMFE) : Constraint(s), vars_(vars), structure_(structure), cutPoint_(cutPoint), posLeft_(posLeft), posRight_(posRight), upperMFE_((float)upperMFE), lowerMFE_((float)lowerMFE) {
		int i;
		n=size();
		v = newRNAplugin(rnaLib,n,dangles);
		v->setEnergyModel(energyModel);
		v->setTemperature(foldTemp);
		if(cutPoint!= -1){
			v->setCutPoint(cutPoint_);
		}
		targetStr = (char *) malloc(sizeof(char)*(n+1));
		for(i=1;i<=n;i++){
			if(structure[i]<0){
					targetStr[i-1]='.';
			}
			else if(structure[i]<i){
				targetStr[i-1]=')';
			}
			else{
				targetStr[i-1]='(';
			}
		}  
		targetStr[n]='\0';
/*		// TIME FOR RESTART
		nRestarts=0;
		time(&timeLastFail);*/

	}

	virtual ~ViennaConstraint() {
		free(targetStr);
		delete v;
	}

	// Adds observers (named Demon) to variable events. These demons are
	// responsible for implementing the propagation algorithm of the
	// constraint.
	virtual void Post() {
		// Create a demon 'global_demon' that will bind events on
		// variables to the calling of the 'InitialPropagate()' method. As
		// this method is expensive, 'global_demon' has a low priority. As
		// such, InitialPropagate will be called after all normal demons
		// and constraints have reached a fixed point. Note
		// that ownership of the 'global_demon' belongs to the solver.

		//Demon* const global_demon = solver()->MakeDelayedConstraintInitialPropagateCallback(this);
		Demon* const global_demon = solver()->MakeConstraintInitialPropagateCallback(this);

		// Attach to all variables.
		for (int i = 0; i < size(); ++i) {
			vars_[i]->WhenBound(global_demon);
		}
	}


	protected:
		std::vector<IntVar*> vars_;
		const std::vector<int64> structure_;
		int cutPoint_;
		int posLeft_;
		int posRight_;
		int64 size() const { return vars_.size(); };
		RNAPlugin* v;
		int n;
		char* targetStr;
		double upperMFE_;
		double lowerMFE_;
/*		// TIME FOR RESTART
		int nRestarts;
		time_t timeLastFail;*/
		
};

class ViennaConstraintUndet : public ViennaConstraint {
	public:
		ViennaConstraintUndet(Solver* const s, const std::vector<IntVar*>& vars, const std::vector<int64>& structure, int dangles, std::string rnaLib, std::string energyModel, double foldTemp, int cutPoint, int posLeft, int posRight, double upperMFE, double lowerMFE) : ViennaConstraint(s,vars,structure,dangles,rnaLib,energyModel,foldTemp,cutPoint, posLeft, posRight,upperMFE, lowerMFE) {}
		// This is the main propagation method.
		//
		// It scans all variables are bound, storing those variables with more than one element in its domain
		//
		// If all variables are bound, calculates the MFE structure and compares with the target structure
		//
		// If the structure is not correct and there are no more possible values it sends the fail signal
		//
		// I not, it removes the assigned value from the domain of one the free variables

		virtual void InitialPropagate();

};

class ViennaConstraintDet : public ViennaConstraint {
	public:
		ViennaConstraintDet(Solver* const s, const std::vector<IntVar*>& vars, const std::vector<int64>& structure, int dangles, std::string rnaLib, std::string energyModel, double foldTemp, int cutPoint, int posLeft, int posRight, double upperMFE, double lowerMFE) : ViennaConstraint(s,vars,structure,dangles,rnaLib,energyModel,foldTemp,cutPoint, posLeft, posRight,upperMFE, lowerMFE) {}	
		// This is the main propagation method.
		//
		// It scans all variables are bound, storing those variables with more than one element in its domain
		//
		// If all variables are bound, calculates the MFE structure and compares with the target structure
		//
		// If the structure is not correct and there are no more possible values it sends the fail signal
		//
		// I not, it removes the assigned value from the domain of one the free variables
		virtual void InitialPropagate();
};


}
