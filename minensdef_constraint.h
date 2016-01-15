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

#include "energy_constant.h"
#include "rna_plugin.h"


namespace operations_research {
// This constraint maintains:
// 
// 
class MinEnsDefConstraint : public Constraint {
	public:

	// This constructor does not take any ownership on its arguments.
	MinEnsDefConstraint(Solver* const s, const std::vector<IntVar*>& vars, int dangles, std::string rnaLib, std::string energyModel, int* bpList, double foldTemp, int cutPoint) : Constraint(s), bpList_(bpList), cutPoint_(cutPoint), vars_(vars) {
		n=size();
		v = newRNAplugin(rnaLib,n,dangles);
		v->setEnergyModel(energyModel);
		v->setTemperature(foldTemp);
		if(cutPoint!= -1){
			v->setCutPoint(cutPoint_);
		}
		vMaxEnsDef = s->MakeIntVar(0,n*ED_PRECISION, StringPrintf("MaxEnsDef"));
		obj = s->MakeMinimize(vMaxEnsDef,1);
		
		
	}

	virtual ~MinEnsDefConstraint() {
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

		Demon* const global_demon = solver()->MakeDelayedConstraintInitialPropagateCallback(this);
		//Demon* const global_demon = solver()->MakeConstraintInitialPropagateCallback(this);

		// Attach to all variables.
		for (int i = 0; i < size(); ++i) {
			vars_[i]->WhenBound(global_demon);
		}
	}

	virtual void InitialPropagate();

	public:
		OptimizeVar* obj;
		IntVar* vMaxEnsDef;
		
	protected:
		int* bpList_;
		int cutPoint_;
		std::vector<IntVar*> vars_;
		int64 size() const { return vars_.size(); };
		RNAPlugin* v;
		int n;
};


}
