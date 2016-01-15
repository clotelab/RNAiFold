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

#include <iostream>
#include "value_heuristic.h"

namespace operations_research {
	
	// FirstUnboundSelector functions
	IntVar* FirstUnboundSelector::Select(Solver* const s, int64* id) {
		for (int i = first_; i < vars_.size(); ++i) {
			IntVar* const var = vars_[i];
			if (!var->Bound()) {
				s->SaveAndSetValue(&first_, i);
				*id = i;
				return var;
			}
		}
		s->SaveAndSetValue(&first_, static_cast<int>(vars_.size()));
		*id = vars_.size();
		return nullptr;
	}

	int64 RandomValueSelector::Select(const IntVar* const v, int64 id){
		const uint64 span = v->Max() - v->Min() + 1;
		Solver* const s = v->solver();		
		for (;;) {
			const int64 value = v->Min() + s->Rand64(span);
			if (v->Contains(value)) {
				return value;
			}
		}
		// Never reached, remove warnings
		return v->Max();
	}
	
	int64 UpValueSelector::Select(const IntVar* const v, int64 id) {
		switch(v->Size()){
			case 4:
			case 2:
				return v->Min();
				break;
			case 3:
			case 1:
			default:
				return v->Max();
				break;
		}
		// Never reached, remove warnings
		return v->Max();
	}

	int64 UpValueSelectorContains::Select(const IntVar* const v, int64 id) {
		for(int i=0;i<updomain.size();i++){
			if(v->Contains(updomain[i])){
				return updomain[i];
			}
		}
		// Never reached, remove warnings
		return v->Max();
	}

	int64 UpValueSelectorContainsRand::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		int randVal;
		for(int i=0;i<updomain.size();i++){
			if(v->Contains(updomain[i])){
				retVal=updomain[i];
				randVal=rand()%100;
				if (randVal < threshold_){
					return retVal;
				}
			}
		}
		return retVal;
		// Never reached, remove warnings
		return v->Max();
	}
	
	int64 UpValueSelectorClosing::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		int randVal;
		for(int i=0;i<updomainClosing.size();i++){
			if(v->Contains(updomainClosing[i])){
				retVal=updomainClosing[i];
				randVal=rand()%100;
				if (randVal < threshold_){
					return retVal;
				}
			}
		}
		return retVal;
		// Never reached, remove warnings
		return v->Max();
	}
	
	int64 BpValueSelector::Select(const IntVar* const v, int64 id) {
		for(int i=0;i<bpdomain.size();i++){
			if(v->Contains(bpdomain[i])){
				return bpdomain[i];
			}
		}
		// Never reached, remove warnings
		return v->Max();
	}

	int64 BpValueSelectorRand::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		int randVal;
		for(int i=0;i<bpdomain.size();i++){
			if(v->Contains(bpdomain[i])){
				retVal=bpdomain[i];
				randVal=((double) rand() / (RAND_MAX));
				if (randVal < threshold_){
					return retVal;
				}
			}
		}
		return retVal;
	}

	int64 BpValueSelectorMultiRand::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		int randVal;
		for(int i=0;i<bpdomainMulti.at(type_).size();i++){
			if(v->Contains(bpdomainMulti.at(type_)[i])){
				retVal=bpdomainMulti.at(type_)[i];
				randVal=rand()%100;
				if (randVal < threshold_){
					return retVal;
				}
			}
		}
		return retVal;
	}
	
	int64 BpValueSelectorBoltz::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		int randVal;
		for(int i=0;i<bpdomain.size();i++){
			if(v->Contains(bpdomain[i])){
				retVal=bpdomain[i];
				randVal=rand()%100;
				if(randVal<=stackProbabilities[i]){
					return retVal;
				}
			}
		}
		return retVal;
	}
	int64 StackValueSelectorContains::Select(const IntVar* const v, int64 id) {
		vector<int>* currentDomain;
		if(vars_[id-1]->Bound()){
			currentDomain = &stackDomain[vars_[id-1]->Value()];
		}
		else{
			currentDomain = &stackDomain[0];
		}

		for(int i=0;i<currentDomain->size();i++){
			if(v->Contains((*currentDomain)[i])){
				return (*currentDomain)[i];
			}
		}
		// Never reached, remove warnings
		return v->Max();
	}

	int64 StackValueSelectorContainsRand::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		double randVal;
				
		vector<int> currentDomain;
		if(vars_[id-1]->Bound()){
			currentDomain = stackDomain[vars_[id-1]->Value()];
		}
		else{
			currentDomain = stackDomain[0];
		}
		for(int i=0;i<currentDomain.size();i++){
			if(v->Contains(currentDomain[i])){
				retVal=currentDomain[i];
				randVal=rand()%100;
				if(randVal<threshold_){
					return retVal;
				}
			}
		}
		
		return retVal;
	}

	int64 StackValueSelectorContainsBoltz::Select(const IntVar* const v, int64 id) {
		int64 retVal;
		double randVal;
		
		vector<int> currentDomain;
		vector<double> currentProbs;
		if(vars_[id-1]->Bound()){
			currentDomain = stackDomain[vars_[id-1]->Value()];
			currentProbs = stackProbabilities[vars_[id-1]->Value()];
		}
		else{
			currentDomain = stackDomain[0];
			currentProbs = stackProbabilities[0];			
		}
		for(int i=0;i<currentDomain.size();i++){
			if(v->Contains(currentDomain[i])){
				retVal=currentDomain[i];
				randVal=((double) rand() / (RAND_MAX));
				if(randVal<=currentProbs[i]){
					return retVal;
				}
			}
		}
		
		return retVal;
	}



	// VariableAssignmentSelector functions	
	std::string VariableAssignmentSelector::DebugString() const {
		return var_selector_->DebugString() + "_" + var_selector_->VarDebugString() + "\n Value Selectors: " + std::to_string(value_selectors_.size());
		//return var_selector_->DebugString() + "_" + value_selectors_[0]->DebugString() + var_selector_->VarDebugString();
	}
	
	
	// AssignOneVariableValue functions	

	AssignOneVariableValue::AssignOneVariableValue(IntVar* const v, int64 val) : var_(v), value_(val) {}

	std::string AssignOneVariableValue::DebugString() const {
		return StringPrintf("[%s == %" GG_LL_FORMAT "d]", var_->DebugString().c_str(),	value_);
	}

	void AssignOneVariableValue::Apply(Solver* const s) { var_->SetValue(value_); }

	void AssignOneVariableValue::Refute(Solver* const s) {	var_->RemoveValue(value_);}


	// BaseAssignVariables functions	


	BaseAssignVariables::~BaseAssignVariables() {}

	Decision* BaseAssignVariables::Next(Solver* const s) {
		int64 id = 0;
		// Assign first variables fixed in the restart
		if(!fixedVars.empty()){
			IntVar* const variable = fixedVars.back();
			fixedVars.pop_back();
			int64 const val = fixedValues.back();
			fixedValues.pop_back();			
			//printf("Assigning position %d\n",variable->Size());
			
			return s->RevAlloc(new AssignOneVariableValue(variable, val));
		}
		IntVar* const var = selector_->SelectVariable(s, &id);
		//printf("Variable Selected %d\n",id);
		if (nullptr != var) {
			const int64 value = selector_->SelectValue(var, id);
			//printf("Value Selected %d\n\n",value);			
			return s->RevAlloc(new AssignOneVariableValue(var, value));
		}
		return nullptr;
	}

	std::string BaseAssignVariables::DebugString() const {
		return selector_->DebugString();
	}

	
	// MakePhase	

	BaseAssignVariables* MakePhase(Solver* s, const std::vector<IntVar*>& vars, vector<int> value_selector_type, int upthreshold, int bpthreshold) {
		//VariableSelector
		FirstUnboundSelector* var_selector = s->RevAlloc(new FirstUnboundSelector(vars));
		if(upthreshold !=100 || bpthreshold !=100){
			long seed = time(NULL);
			srand(seed);
		}
		vector<ValueSelector*> value_selectors;
		value_selectors.push_back(s->RevAlloc(new MinValueSelector));
		value_selectors.push_back(s->RevAlloc(new UpValueSelector));
		if(upthreshold ==100){
			value_selectors.push_back(s->RevAlloc(new UpValueSelectorContains));
		}
		else{
			value_selectors.push_back(s->RevAlloc(new UpValueSelectorContainsRand(upthreshold)));
		}
		if(bpthreshold ==100){
			value_selectors.push_back(s->RevAlloc(new BpValueSelector));
			value_selectors.push_back(s->RevAlloc(new StackValueSelectorContains(vars)));
		}
		else if(bpthreshold ==0){
			value_selectors.push_back(s->RevAlloc(new BpValueSelectorBoltz));
			value_selectors.push_back(s->RevAlloc(new StackValueSelectorContainsBoltz(vars)));
		}
		else {
			value_selectors.push_back(s->RevAlloc(new BpValueSelectorRand(bpthreshold)));
			value_selectors.push_back(s->RevAlloc(new StackValueSelectorContainsRand(vars, bpthreshold)));
		}
		value_selectors.push_back(s->RevAlloc(new RandomValueSelector));		

		// Value selector for multiple structures
		value_selectors.push_back(s->RevAlloc(new BpValueSelectorMultiRand(bpthreshold,0)));
		value_selectors.push_back(s->RevAlloc(new BpValueSelectorMultiRand(bpthreshold,1)));
		value_selectors.push_back(s->RevAlloc(new BpValueSelectorMultiRand(bpthreshold,2)));
		value_selectors.push_back(s->RevAlloc(new BpValueSelectorMultiRand(bpthreshold,3)));
		value_selectors.push_back(s->RevAlloc(new BpValueSelectorMultiRand(bpthreshold,4)));
		value_selectors.push_back(s->RevAlloc(new BpValueSelectorMultiRand(bpthreshold,5)));
		value_selectors.push_back(s->RevAlloc(new UpValueSelectorClosing(upthreshold)));



		// BaseVariableAssignmentSelector
		VariableAssignmentSelector* selector = s->RevAlloc(new VariableAssignmentSelector(var_selector, value_selectors, value_selector_type));
		
		return s->RevAlloc(new BaseAssignVariables(selector));
	}	

	
// OTHER CLASSES (REMOVE)

// AssignOneNtValue

	AssignOneNtValue::AssignOneNtValue(IntVar* const v, int64 val)
	    : var_(v), value_(val) {}

	std::string AssignOneNtValue::DebugString() const {
	  return StringPrintf("[%s == %" GG_LL_FORMAT "d]", var_->DebugString().c_str(),
		              value_);
	}

	void AssignOneNtValue::Apply(Solver* const s) { var_->SetValue(value_); }

	void AssignOneNtValue::Refute(Solver* const s) {
	  var_->RemoveValue(value_);
	}


}
