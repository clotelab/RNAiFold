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
#include <map>
#include <algorithm>

#ifndef _VALUE_HEURISTIC
#define _VALUE_HEURISTIC

#include "base/logging.h"
#include "util/string_array.h"
#include "constraint_solver/constraint_solver.h"

#define NO_RAND_ASSIGN -600

#define VH_MIN_VALUE 0
#define VH_UP_CASE 1
#define VH_UP_CONTAINS 2
#define VH_BP 3
#define VH_BP_STACK 4 // Order by stacking energy
#define VH_RANDOM 5
#define VH_MULTI_BP_0 6  // Unpaired in other structure 
#define VH_MULTI_BP_1 7  // Same base pair in other structure 
#define VH_MULTI_BP_2 8  //  Different base pair in other structure //GC-CG-GU-UG-AU-UA //-6, 6,-11,11,-9,9 // BPType == 2
#define VH_MULTI_BP_3 9  //  Opening unpaired in other structure //GC-GU-UG-UA-CG-AU //-6,-11,11,9,6,-9 // BPType == 3
#define VH_MULTI_BP_4 10 //  Closing unpaired in other structure //CG-UG-GU-AU-GC-UA //6,11,-11,-9,-6, 9 // BPType == 4
#define VH_MULTI_BP_5 11 //  Opening or closing unpaired in other structure but adjacent paired with the same position //UG-GU-AU-UA-CG-GC //11,-11,-9, 9, 6, -6 // BPType == 5
#define VH_MULTI_UP_0 VH_UP_CONTAINS  //  Unpaired position is also unpaired in other structures // A-U-C-G
#define VH_MULTI_UP_1 VH_RANDOM //  Unpaired position is a normal base pair another structures // RAndom
#define VH_MULTI_UP_2 12 //  Unpaired position is a closing base pair in other structure // G-C-U-A
#define VH_THREEP_UP 13 //  Unpaired position is a three prime dangle // Avoid pairing with the 5'
#define VH_BP_STACK_INV 14 // Inverse to stacking energy 

struct IdxCompare
{
    const std::vector<int>& target;

    IdxCompare(const std::vector<int>& target): target(target) {}

    bool operator()(int a, int b) const { return target[a] < target[b]; }
};

namespace operations_research {
	// Base value domains
	const std::vector<int> updomain {0,3,1,2};
	const std::vector<int> bpdomain {-6,6,-9,9,-11,11}; //GC=-6,CG=6, AU=-9, UA=9, GU=-11, UG=11
	
	// Variables for multiple structures assignment
	const std::vector<std::vector<int>> bpdomainMulti {{11,-9,-11,9,6,-6}, // Type 0 - Inverse of normal order
	                                                   {-6,6,-9,9,-11,11}, // Type 1 - Normal order
	                                                   {-6,6,-11,11,-9,9}, // Type 2
	                                                   {-6,-11,11,-9,6,9}, // Type 3
	                                                   {6,11,-11,-9,-6,9}, // Type 4
	                                                   {11,-11,-9,9,6,-6}};// Type 5
	const std::vector<int> updomainClosing {2,1,3,0};


	// Variables for energy based assignment
	const std::map<int, vector<int> > baseStack {{  0,{0,0,200,200,200,200}},
	                                             { -6,{-240,-330,-210,-210,-210,-140}},
	                                             {  6,{-330,-340,-220,-240,-250,-150}},
	                                             { -9,{-210,-220,-110,-90,-140,-60}},
	                                             {  9,{-210,-240,-90,-130,-130,-100}},
	                                             {-11,{-210,-250,-140,-130,-130,-50}},
	                                             { 11,{-140,-150,-60,-100,-50,30}}};

	const std::map<int, vector<int> > baseStackInv {{  0,{0,0,-200,-200,-200,-200}},
	                                                { -6,{240,330,210,210,210,140}},
	                                                {  6,{330,340,220,240,250,150}},
	                                                { -9,{210,220,110,90,140,60}},
	                                                {  9,{210,240,90,130,130,100}},
	                                                {-11,{210,250,140,130,130,50}},
	                                                { 11,{140,150,60,100,50,-30}}};
	                                                
	const std::vector<int> baseEnergyMulti {-290,-223,-190,-137,-103,-10};


// replaced inherit VariableSelector

	class FirstUnboundSelector : public BaseObject {
		public:
			explicit FirstUnboundSelector(const std::vector<IntVar*>& vars)	: vars_(vars), first_(0) {}
			virtual ~FirstUnboundSelector() {}
			virtual IntVar* Select(Solver* const s, int64* id);
			void Accept(ModelVisitor* const visitor) const {
				visitor->BeginVisitExtension(ModelVisitor::kVariableGroupExtension);
				visitor->VisitIntegerVariableArrayArgument(ModelVisitor::kVarsArgument, vars_);
				visitor->EndVisitExtension(ModelVisitor::kVariableGroupExtension);
			}

			virtual std::string DebugString() const { return "ChooseFirstUnbound"; }

			std::string VarDebugString() const {
				return StringPrintf("(%s)", JoinDebugStringPtr(vars_, ", ").c_str());
			}

		protected:
			const std::vector<IntVar*> vars_;
		private:
			int first_;

	};


	// ----- Value selectors -----

	class ValueSelector : public BaseObject {
		public:
			ValueSelector() {}
			virtual ~ValueSelector() {}
			virtual int64 Select(const IntVar* const v, int64 id) = 0;
	};

	class MinValueSelector : public ValueSelector {
		public:
			MinValueSelector() {}
			virtual ~MinValueSelector() {}
			virtual int64 Select(const IntVar* const v, int64 id) { return v->Min(); }
			std::string DebugString() const { return "AssignMin"; }
	};

	class RandomValueSelector : public ValueSelector {
		public:
			RandomValueSelector() {}
			virtual ~RandomValueSelector() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignRandom"; }
	};

	class UpValueSelector : public ValueSelector {
		public:
			UpValueSelector() {}
			virtual ~UpValueSelector() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignOrderedUp"; }
	};

	class UpValueSelectorContains : public ValueSelector {
		public:
			UpValueSelectorContains() {}
			virtual ~UpValueSelectorContains() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignOrderedUpContains"; }
	};

	class UpValueSelectorContainsRand : public ValueSelector {
		public:
			UpValueSelectorContainsRand(int threshold) : threshold_(threshold) {}
			virtual ~UpValueSelectorContainsRand() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignOrderedUpContainsRand"; }
		protected:
			int threshold_;
	};

	class UpValueSelectorClosing : public ValueSelector {
		public:
			UpValueSelectorClosing(int threshold) : threshold_(threshold) {}
			virtual ~UpValueSelectorClosing() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignOrderedUpClosing"; }
		protected:
			int threshold_;
	};
	
	class BpValueSelector : public ValueSelector {
		public:
			BpValueSelector() {}
			virtual ~BpValueSelector() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignBP"; }

	};

	class BpValueSelectorRand : public ValueSelector {
		public:
			BpValueSelectorRand(int threshold) : threshold_(threshold) {}
			virtual ~BpValueSelectorRand() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignBPRand"; }
		protected:
			int threshold_;

	};

	class BpValueSelectorMultiRand : public ValueSelector {
		public:
			BpValueSelectorMultiRand(int threshold, int type) : threshold_(threshold), type_(type) {}
			virtual ~BpValueSelectorMultiRand() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignBPMultiRand"; }
		protected:
			int threshold_;
			int type_;

	};

		
	class BpValueSelectorBoltz : public ValueSelector {
		public:
			BpValueSelectorBoltz() {
				stackProbabilities= vector<double> {0.4638510717,0.8651534066,0.25,0.3333333333,0.5,1};
			}
			virtual ~BpValueSelectorBoltz() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "AssignBPRand"; }
		protected:
			vector<double> stackProbabilities;						

	};

	class StackValueSelectorContains : public ValueSelector {
		public:
			StackValueSelectorContains(const std::vector<IntVar*>& vars) : vars_(vars) {
				stackDomain[0] = vector<int> {-6,6,-9,9,-11,11};
				stackDomain[-6] = vector<int> {-6,6,-9,9,-11,11};
				stackDomain[6] = vector<int> {-6,6,-11,9,-9,11};
				stackDomain[-9] = vector<int> {-6,6,-11,-9,9,11};
				stackDomain[9] = vector<int> {-6,6,9,-11,11,-9};
				stackDomain[-11] = vector<int> {-6,6,-9,9,11,-11};
				stackDomain[11] = vector<int> {-6,6,9,-9,-11,11};	
			}
			virtual ~StackValueSelectorContains() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "StackValueSelectorContains"; }

		protected:
			const std::vector<IntVar*> vars_;
			std::map<int, vector<int> > stackDomain;
	};

	class StackValueSelectorContainsBoltz : public ValueSelector {
		public:
			StackValueSelectorContainsBoltz(const std::vector<IntVar*>& vars) : vars_(vars) {
				stackDomain[0] = vector<int> {-6,6,-9,9,-11,11};
				stackDomain[-6] = vector<int> {-6,6,-9,9,-11,11};
				stackDomain[6] = vector<int> {-6,6,-11,9,-9,11};
				stackDomain[-9] = vector<int> {-6,6,-11,-9,9,11};
				stackDomain[9] = vector<int> {-6,6,9,-11,11,-9};
				stackDomain[-11] = vector<int> {-6,6,-9,9,11,-11};
				stackDomain[11] = vector<int> {-6,6,9,-9,-11,11};	


				stackProbabilities[0] = vector<double> {0.4638510717,0.8651534066,0.25,0.3333333333,0.5,1};
				stackProbabilities[-6] = vector<double> {0.4051334918,0.5790458752,0.230863874,0.4152255624,0.8351439873,1};
				stackProbabilities[6] = vector<double> {0.5861354829,0.3288129468,0.301098345,0.4308164717,0.7569025635,1};
				stackProbabilities[-9] = vector<double> {0.4020874815,0.5717645813,0.4288200348,0.4614292616,0.6193427597,1};
				stackProbabilities[9] = vector<double> {0.4670415686,0.5385986108,0.3187579885,0.4679071212,0.5404745369,1};
				stackProbabilities[-11] = vector<double> {0.5335731621,0.5977888617,0.4773469243,0.7765241403,0.948853341,1};
				stackProbabilities[11] = vector<double> {0.3599718785,0.4781937122,0.4788857635,0.4802157077,0.7855024678,1};	

			}
			virtual ~StackValueSelectorContainsBoltz() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "StackValueSelectorContainsBoltz"; }

		protected:
			const std::vector<IntVar*> vars_;
			std::map<int, vector<int> > stackDomain;
			std::map<int, vector<double> > stackProbabilities;			
	};

	class StackValueSelectorContainsRand : public ValueSelector {
		public:
			StackValueSelectorContainsRand(const std::vector<IntVar*>& vars,int threshold) : vars_(vars),threshold_(threshold) {
				stackDomain[0] = vector<int> {-6,6,-9,9,-11,11};
				stackDomain[-6] = vector<int> {-6,6,-9,9,-11,11};
				stackDomain[6] = vector<int> {-6,6,-11,9,-9,11};
				stackDomain[-9] = vector<int> {-6,6,-11,-9,9,11};
				stackDomain[9] = vector<int> {-6,6,9,-11,11,-9};
				stackDomain[-11] = vector<int> {-6,6,-9,9,11,-11};
				stackDomain[11] = vector<int> {-6,6,9,-9,-11,11};				
			}
			virtual ~StackValueSelectorContainsRand() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "StackValueSelectorContainsRand"; }

		protected:
			const std::vector<IntVar*> vars_;
			std::map<int, vector<int> > stackDomain;
			int threshold_;
	};


	class ThreePrimeValueSelector : public ValueSelector {
		public:
			ThreePrimeValueSelector(const std::vector<IntVar*>& vars,int threshold) : vars_(vars),threshold_(threshold) {
//				fivePrimeDom[0] = vector<int> {0,2,1,3};
//				fivePrimeDom[1] = vector<int> {1,0,3,2};
//				fivePrimeDom[2] = vector<int> {2,0,3,1};
//				fivePrimeDom[3] = vector<int> {3,1,2,0};
				fivePrimeDom[0] = vector<int> {0,2,3,1};
				fivePrimeDom[1] = vector<int> {1,0,3,2};
				fivePrimeDom[2] = vector<int> {2,0,3,1};
				fivePrimeDom[3] = vector<int> {3,2,0,1};

				fivePrimeDom[4] = updomain; 
			}
			virtual ~ThreePrimeValueSelector() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "ThreePrimeValueSelector"; }

		protected:
			const std::vector<IntVar*> vars_;
			std::map<int, vector<int> > fivePrimeDom;
			int threshold_;
	};


	class SimplifiedBPenergyHeuristic : public ValueSelector {
		public:
			SimplifiedBPenergyHeuristic(const std::vector<IntVar*>& vars,int threshold, const vector<int>& value_selector_type) : vars_(vars),threshold_(threshold),value_selector_type_(value_selector_type) {
				nvars_ = vars.size();
				if(random_.size() == 0){
					orderArray_.resize(nvars_);
					random_.resize(nvars_);
					baseOrder_.resize(nvars_);
					for(int i=0; i<nvars_;i++){
						switch(value_selector_type_.at(i)){
							case VH_BP:
								baseOrder_[i] = bpdomain;
								break;
							case VH_BP_STACK_INV:
								baseOrder_[i] = bpdomain;
								break;
							case VH_BP_STACK:
								baseOrder_[i] = bpdomain;
								break;
							case VH_MULTI_BP_0:
								baseOrder_[i] = bpdomainMulti.at(0);
								break;
							case VH_MULTI_BP_1:
								baseOrder_[i] = bpdomainMulti.at(1);
								break;
							case VH_MULTI_BP_2:
								baseOrder_[i] = bpdomainMulti.at(2);
								break;
							case VH_MULTI_BP_3:
								baseOrder_[i] = bpdomainMulti.at(3);
								break;
							case VH_MULTI_BP_4:
								baseOrder_[i] = bpdomainMulti.at(4);
								break;
							case VH_MULTI_BP_5:
								baseOrder_[i] = bpdomainMulti.at(5);
								break;						
						}
					} 
					for(int i=0; i<nvars_;i++){
						orderArray_[i]=baseOrder_[i];
						random_[i].resize(bpdomain.size(),NO_RAND_ASSIGN);
					}
				}
			}
			virtual ~SimplifiedBPenergyHeuristic() {}
			virtual int64 Select(const IntVar* const v, int64 id);
			std::string DebugString() const { return "SimplifiedBPenergyHeuristic"; }

		protected:
			static std::vector<std::vector<int> > orderArray_;
			static std::vector<std::vector<int> > random_;
			static std::vector<std::vector<int> > baseOrder_;	
		
			const std::vector<IntVar*> vars_;
			const std::vector<int> value_selector_type_;
			int threshold_;
			int nvars_;
	};

	
// replaced inherit public BaseVariableAssignmentSelector
	class VariableAssignmentSelector : public BaseObject  {
		public:
			VariableAssignmentSelector(FirstUnboundSelector* const var_selector, vector<ValueSelector*> const value_selectors, vector<int> const value_selector_type) : var_selector_(var_selector), value_selectors_(value_selectors), value_selector_type_(value_selector_type) {}
			virtual ~VariableAssignmentSelector() {}
			virtual int64 SelectValue(const IntVar* const var, int64 id) {
				//printf("Variable %d -- ValueSelector %s -- DomSize %d\n",id,value_selectors_[value_selector_type_[id]]->DebugString().c_str(), var->Size());
				return value_selectors_[value_selector_type_[id]]->Select(var, id);
			}
			virtual IntVar* SelectVariable(Solver* const s, int64* id) {
				return var_selector_->Select(s, id);
			}
			std::string DebugString() const;

			virtual void Accept(ModelVisitor* const visitor) const {
				var_selector_->Accept(visitor);
			}

		private:
			FirstUnboundSelector* const var_selector_;
			vector<ValueSelector*> const value_selectors_;
			vector<int> const value_selector_type_;
	};






	class AssignOneVariableValue : public Decision {
		public:
			AssignOneVariableValue(IntVar* const v, int64 val);
			virtual ~AssignOneVariableValue() {}
			virtual void Apply(Solver* const s);
			virtual void Refute(Solver* const s);
			virtual std::string DebugString() const;
			virtual void Accept(DecisionVisitor* const visitor) const {
				visitor->VisitSetVariableValue(var_, value_);
			}

		private:
			IntVar* const var_;
			int64 value_;
	};






	class BaseAssignVariables : public DecisionBuilder {
		public:
			explicit BaseAssignVariables(VariableAssignmentSelector* const selector) : selector_(selector) {}
			virtual ~BaseAssignVariables();
			virtual Decision* Next(Solver* const s);
			virtual std::string DebugString() const;
			virtual void Accept(ModelVisitor* const visitor) const {
				selector_->Accept(visitor);
			}
			void addFixedPosition(IntVar* position, int64 value)  {
				fixedVars.push_back(position);
				fixedValues.push_back(value);
				
			}			
			void clearFixedPositions()  {
				fixedVars.clear();
				fixedValues.clear();				
			}			

		protected:
			VariableAssignmentSelector* const selector_;
		private:
			vector<IntVar*> fixedVars;
			vector<int64> fixedValues;			
	};

	BaseAssignVariables* MakePhase(Solver* s, const std::vector<IntVar*>& vars, vector<int> value_selector_type, int upthreshold, int bpthreshold);






	// OTHER CLASSES (REMOVE)
	class AssignOneNtValue : public Decision {
		public:
			AssignOneNtValue(IntVar* const v, int64 val);
			virtual ~AssignOneNtValue() {}
			virtual void Apply(Solver* const s);
			virtual void Refute(Solver* const s);
			virtual std::string DebugString() const;
			virtual void Accept(DecisionVisitor* const visitor) const {
				visitor->VisitSetVariableValue(var_, value_);
			}

		private:
			IntVar* const var_;
			int64 value_;
	};

	class AssignNtFromLastAssignment : public DecisionBuilder {
		public:
			AssignNtFromLastAssignment(const Assignment* const assignment,
							DecisionBuilder* const db,
							const std::vector<IntVar*>& vars)
			: assignment_(assignment), db_(db), vars_(vars), iter_(0) {}

			~AssignNtFromLastAssignment() {}

			Decision* Next(Solver* const s) {
				if (iter_ < vars_.size()) {
					IntVar* const var = vars_[iter_++];
					return s->RevAlloc(new AssignOneNtValue(var, assignment_->Value(var)));
				} else {
					return db_->Next(s);
				}
			}

			virtual void Accept(ModelVisitor* const visitor) const {
				visitor->BeginVisitExtension(ModelVisitor::kVariableGroupExtension);
				visitor->VisitIntegerVariableArrayArgument(ModelVisitor::kVarsArgument,
										   vars_);
				visitor->EndVisitExtension(ModelVisitor::kVariableGroupExtension);
			}

		private:
			const Assignment* const assignment_;
			DecisionBuilder* const db_;
			const std::vector<IntVar*> vars_;
			int iter_;
	};


}
#endif
