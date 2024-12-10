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

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include "strtreegroup.h"
#include "misc.h"
#include "value_heuristic.h"


namespace operations_research {
	StrTreeGroup::~StrTreeGroup(){};
	StrTreeGroup::StrTreeGroup(vector<StrTree*> *str_trees) : str_trees_(str_trees) {
		nHelices=0;
		for(int i=0;i<str_trees_->size();i++){
			nHelices += str_trees_->at(i)->getHelices().size();
		}
		helixIndex.assign(nHelices, 0);
		posIndex.assign(nHelices, 0);
		order_.assign(nHelices, 0);
		for(int i=0;i<nHelices;i++){
			int str1=-1;
			int pos1=-1;
			int sum=0;
			for(int k=0;k<str_trees->size();k++){
				if(str1 == -1 && i-sum<str_trees_->at(k)->getHelices().size()){
					str1=k;
					pos1=i-sum;
					break;
				}
				sum+=str_trees_->at(k)->getHelices().size();
			}
			helixIndex[i]=str1;
			posIndex[i]=pos1;			
		}
	}

	int StrTreeGroup::getStrId(int index){
		return helixIndex[index];
	}
	
	int StrTreeGroup::getHelixPos(int index){
		return posIndex[index];
	}
	CPHelix* StrTreeGroup::getHelix(int index){
		return str_trees_->at(helixIndex[index])->getHelices().at(posIndex[index]);
	}
	int StrTreeGroup::getNumHelices(){
		return nHelices;
	}
	std::vector<std::vector<int>> StrTreeGroup::optimizeValueHeuristic(int n,std::vector<std::vector<int> > BPO,std::vector<std::vector<int> > BPC,std::vector<int*> str_int,std::vector<double> temps){ //, int []bp, int []BPO, int []BPC, int []UP, int []bpint strategy
		std::vector<std::vector<int>> BPtypes(str_int.size());
		if(str_int.size()==1){
			return BPtypes;
		}
		for(int h=0; h<str_int.size();h++){
			BPtypes.at(h) = std::vector<int>(BPO.at(h).size(),VH_MULTI_BP_0);
		}
		for(int h=0; h<str_int.size();h++){
//			cout << "Str "<< h << " - BPO:";
//			for(int i=0;i< BPO.at(h).size();i++){
//				cout << "\t" << BPO.at(h)[i];
//			}
//			cout << endl << "Str "<< h << " - BPC:";
//			for(int i=0;i< BPC.at(h).size();i++){
//				cout << "\t" << BPC.at(h)[i];
//			}
//			cout << endl << "Str "<< h << " -Type:";
			for(int i=0;i< BPO.at(h).size();i++){				
				int newBPtype = 0;
				if((BPO.at(h)[i] != 1 && str_int.at(h)[BPO.at(h)[i]-1] == -1 && str_int.at(h)[BPO.at(h)[i]+1] == -1) ||  
				   (BPC.at(h)[i] != n && str_int.at(h)[BPC.at(h)[i]-1] == -1 && str_int.at(h)[BPC.at(h)[i]+1] == -1)){ // If one on the paired bases has unpaired positions at both sides
					newBPtype = VH_MULTI_BP_1;
				}
				else{
					for(int h2=0; h2<str_int.size(); h2++){
						if(h!=h2){
							if(str_int.at(h2)[BPO.at(h)[i]] == -1 && str_int.at(h2)[BPC.at(h)[i]] == -1){	// Unpaired in structure 2
								if(temps.at(h)>temps.at(h2)){
									newBPtype = VH_MULTI_BP_1;
								}
								else{
									newBPtype = VH_MULTI_BP_0;
								}
							}
							else if(str_int.at(h2)[BPO.at(h)[i]] == BPC.at(h)[i]){ // Same base pair in structure 2
								newBPtype = VH_MULTI_BP_1;
							}
							else if(str_int.at(h2)[BPO.at(h)[i]] != -1 && str_int.at(h2)[BPC.at(h)[i]] != -1){ // Different base pair in structure 2
								newBPtype = VH_MULTI_BP_2;
							}
							else if(str_int.at(h2)[BPO.at(h)[i]] == -1){ // Opening unpaired in structure 2 
								if((BPO.at(h)[i] != 1 && str_int.at(h2)[BPO.at(h)[i]-1] == BPC.at(h)[i]) || (str_int.at(h2)[BPO.at(h)[i]+1] == BPC.at(h)[i])){
									newBPtype = VH_MULTI_BP_5;
								}
								else{
									newBPtype = VH_MULTI_BP_3;
								}
							}
							else if(str_int.at(h2)[BPC.at(h)[i]] == -1){ // Closing unpaired in structure 2 
								if((BPC.at(h)[i] != n && str_int.at(h2)[BPC.at(h)[i]+1] == BPO.at(h)[i]) || (str_int.at(h2)[BPC.at(h)[i]-1] == BPO.at(h)[i])){
									newBPtype = VH_MULTI_BP_5;
								}
								else{
									newBPtype = VH_MULTI_BP_4;
								}
							}
							else{
								cout << "WARNING: Not a type!" << endl;
							}
						}
					}
				}
				BPtypes.at(h)[i] = newBPtype;
//				cout << "\t" << newBPtype;
			}			

		}
		return BPtypes;

	}

	std::vector<std::vector<int>> StrTreeGroup::optimizeUPHeuristic(int n,std::vector<std::vector<int> > BPO,std::vector<std::vector<int> > BPC,std::vector<int*> str_int){ 
		std::vector<std::vector<int>> UPtypes(str_int.size(),std::vector<int>(n+1, VH_MULTI_UP_1));
		if(str_int.size()>1){
			for(int h=0; h<str_int.size();h++){
//				cout << "Str"<< h<<" - UP types ";

				for(int i=1;i<=n;i++){
					int UPtype = VH_MULTI_UP_1;
					if(str_int.at(h)[i] == -1){
						for(int h2=0; h2<str_int.size(); h2++){
							if(h!=h2){
								if(str_int.at(h2)[i] == -1){
									UPtype = VH_MULTI_UP_0;
								}
								else{
									vector<CPHelix*> newH = str_trees_->at(h2)->getHelices();
									for (int j=0; j< newH.size();j++){
										std::vector<int> cbps = newH.at(j)->getClosingBPs();
										for(int k=0;k<cbps.size();k++){
											if(i==BPO.at(h2)[cbps[k]] || i==BPC.at(h2)[cbps[k]]){
												UPtype = VH_MULTI_UP_2;
											}
										}
									}
								}
							}
						}
					}
					UPtypes.at(h)[i] = UPtype;
//					cout << UPtype<<" ";
				}
//					cout << endl;
			}
		}
		return UPtypes;


	}
	
	vector<int64> StrTreeGroup::optimizeOrder(int strategy){
		if(str_trees_->size()==1){
			vector<CPHelix* > sortedHelices = str_trees_->at(0)->getSortedHelices();
			for(int i=0;i<nHelices;i++){
				order_[sortedHelices[i]->getId()]=i;
			}
		}
		else if (strategy == HH_NO_OVERLAP){
			std::vector<CPHelix* > sortedHelices;

			/* // Heuristic -1 Not used
			int countHelices=0;
			for(int i=0;i<str_trees_->size();i++){
				std::vector<CPHelix* > tmpHelices = str_trees_->at(i)->getSortedHelices();
				sortedHelices.insert(sortedHelices.end(), tmpHelices.begin(), tmpHelices.end());
				for(int j=0;j<str_trees_->at(i)->getHelices().size();j++){
					order_[sortedHelices[j]->getId()+countHelices]=j+countHelices;
				}
				countHelices = str_trees_->at(i)->getHelices().size();
			}

			// End Heuristic -1 */
			//Heuristic 0			
			std::vector<int> y;
			for(int i=0;i<nHelices;i++){
				y.push_back(i);
				sortedHelices.push_back(getHelix(i));
			}
			//std::sort(y.begin(), y.end(), HelixCompare(sortedHelices));
			std::sort(y.begin(), y.end(), LeavesDistanceCompare(sortedHelices));
			

			for(int i=0;i<nHelices;i++){
				order_[y[i]]= i;
			}
			// End Heuristic 0
		}
		else {
			vector<vector<int> > includes(nHelices, std::vector<int>(nHelices, 0));
			vector<vector<int> > overlap(nHelices, std::vector<int>(nHelices, 0));
			vector<int>  degree(nHelices, 0);
			for(int i=0;i<nHelices;i++){
				for(int j=i+1;j<nHelices;j++){
					if(helixIndex[i]!=helixIndex[j]){
						CPHelix* helix1 = getHelix(i);
						CPHelix* helix2 = getHelix(j);
						if(helix2->getI()==helix1->getI() && helix2->getJ()==helix1->getJ()){
							if(str_trees_->at(helixIndex[i])->getTemp() > str_trees_->at(helixIndex[j])->getTemp()){
								includes[j][i]=1;
							}
							else{
								includes[i][j]=1;
							}
						}

						if((helix1->getI()==helix2->getJ() || helix1->getJ()==helix2->getI()) || // Helixes ovelap at last or first position
						   (helix1->getI()>helix2->getI() && helix1->getI()<helix2->getJ() && helix1->getJ()>helix2->getJ()) ||  // Start of helix A included in helix B
						   (helix1->getJ()>helix2->getI() && helix1->getJ()<helix2->getJ() && helix1->getI()<helix2->getI()) ||   // End of helix A included in helix B
						   (helix1->getI()>=helix2->getI() && helix1->getJ()<=helix2->getJ()) || 
						   (helix1->getI()<=helix2->getI() && helix1->getJ()>=helix2->getJ())){
							if((helix1->getI()<helix2->getI() && helix1->getJ()>=helix2->getJ()) || 
							   (helix1->getI()<=helix2->getI() && helix1->getJ()>helix2->getJ())){
								includes[i][j]=1;
							}
							else if((helix2->getI()<helix1->getI() && helix2->getJ()>=helix1->getJ()) || 
								   (helix2->getI()<=helix1->getI() && helix2->getJ()>helix1->getJ())){
								includes[j][i]=1;						
							}
							
							
							switch(strategy){
								// 1 If helices overlap at any position
								case HH_OVERLAP_SIMPLE:
									overlap[i][j]=1;
									overlap[j][i]=1;
									break;
								// Overlap value is the number of overlapping positions involved in a base pair
								case HH_OVERLAP_BP:
									overlap[i][j]=helix1->getBPintersect(helix2->getBPpositions());
									overlap[j][i]=overlap[i][j];
									break;
								// Overlap value is the number of overlapping positions, no matter if they are base pair or unpaired positions
								case HH_OVERLAP_POSITIONS:
									overlap[i][j]=helix1->getPositionsIntersect(helix2->getPositions());
									overlap[j][i]=overlap[i][j];
									break;
								// Overlap value is the percentage of overlapping positions involved in a base pair over the total number of positions involved in a base pair in each helix 
								case HH_OVERLAP_BP_PERCENT:
									int intersect=helix1->getBPintersect(helix2->getBPpositions())*100;
									overlap[i][j]=helix1->getBPpositions().size() == 0 ? 0 : (intersect/helix1->getBPpositions().size());
									overlap[j][i]=helix2->getBPpositions().size() == 0 ? 0 : (intersect/helix2->getBPpositions().size());
									break;
							}
							
							degree[i]+=overlap[i][j];
							degree[j]+=overlap[j][i];
						}
					}
				}
			}

			
//			cout << "Overlap" << endl;
//			for(int j=0;j<nHelices;j++){
//				cout <<"\t"<< j;
//			}
//			cout << endl;
//			for(int i=0;i<nHelices;i++){
//				cout << i;
//				for(int j=0;j<nHelices;j++){
//					cout << "\t" << overlap[i][j] ;
//				}
//				cout << endl;
//			}

//			cout << "Degree" << endl;
//			for(int j=0;j<nHelices;j++){
//				cout <<"\t"<< degree[j];
//			}
//			cout << endl;


			vector<vector<int> > degreedist(nHelices, std::vector<int>(nHelices, 0));
			for(int i=0;i<nHelices;i++){
				for(int j=i+1;j<nHelices;j++){
					if(degree[i]<degree[j]){
						degreedist[i][j]=degree[j]-degree[i];
					}
					else if(degree[j]<degree[i]){
						degreedist[j][i]=degree[i]-degree[j];
					}
				}
			}
			
//			cout << "Degreedist" << endl;
//			for(int j=0;j<nHelices;j++){
//				cout <<"\t"<< j;
//			}
//			cout << endl;
//			for(int i=0;i<nHelices;i++){
//				cout << i;
//				for(int j=0;j<nHelices;j++){
//					cout << "\t" << degreedist[i][j] ;
//				}
//				cout << endl;
//			}

//			cout << "Includes" << endl;
//			for(int j=0;j<nHelices;j++){
//				cout <<"\t"<< j;
//			}
//			cout << endl;
//			for(int i=0;i<nHelices;i++){
//				cout << i;
//				for(int j=0;j<nHelices;j++){
//					cout << "\t" << includes[i][j] ;
//				}
//				cout << endl;
//			}


			Solver orderSolver("Helix order solver");
			vector<IntVar *> order = vector<IntVar *> (nHelices);
			vector<IntVar *> diff = vector<IntVar *> (nHelices*nHelices);
			for(int i=0; i<nHelices;i++){
				order[i] = orderSolver.MakeIntVar(0,nHelices-1,absl::StrFormat("Order_%03d", i));
				for(int j=0; j<nHelices;j++){
					diff[(i*nHelices)+j]= orderSolver.MakeIntVar(0,degreedist[i][j],absl::StrFormat("Diff_%03d_%03d", i,j));
				}
			}
			
			orderSolver.AddConstraint(orderSolver.MakeAllDifferent(order));
			for(int i=0; i<nHelices;i++){
				if(getHelix(i)->getParent()!= -1){
					orderSolver.AddConstraint(orderSolver.MakeLess(order[i],order[i-posIndex[i]+getHelix(i)->getParent()]));
	//				cout << "Parent" << i << "<" << i-posIndex[i]+getHelix(i)->getParent() << endl;
				}
				for(int j=0; j<nHelices;j++){
					if(includes[i][j] ==1){
	//					cout << "Include" << j << "<" << i << endl;
						orderSolver.AddConstraint(orderSolver.MakeLess(order[j],order[i]));					
					}
					orderSolver.AddConstraint(orderSolver.MakeEquality(diff[(i*nHelices)+j],orderSolver.MakeProd(orderSolver.MakeIsLessOrEqualVar(order[i],order[j]),degreedist[i][j])));
				}
			}
			IntVar* minimizeVar = orderSolver.MakeSum(diff)->Var();
			OptimizeVar* minimizeDiff= orderSolver.MakeMinimize(minimizeVar,1);
			std::vector<SearchMonitor*> monitors;
			monitors.push_back(minimizeDiff);
			// UNCOMMENT TO PRINT SEARCH TREE
			// SearchMonitor* const cpviz = orderSolver.MakeTreeMonitor(order, "configurationOrder.xml","treeOrder.xml","visualizationOrder.xml");
			// monitors.push_back(cpviz);
			monitors.push_back(orderSolver.MakeTimeLimit(30000));			
			DecisionBuilder* db = orderSolver.MakePhase (order,Solver::CHOOSE_FIRST_UNBOUND,Solver::ASSIGN_MIN_VALUE);
			orderSolver.NewSearch(db,monitors);
			while (orderSolver.NextSolution()){
				for(int i=0;i<nHelices;i++){
					order_[i]=order[i]->Value();
				}
//				cout << "Order:";
//				for(int i=0;i<nHelices;i++){
//					cout <<"\t" << order_[i];
//				}
//				cout <<endl;
//				for(int i=0;i<nHelices;i++){
//					cout << i;
//					for(int j=0;j<nHelices;j++){
//						cout << "\t" << diff[(i*nHelices)+j]->Value();
//					}
//					cout << endl;
//				}

//				cout << "Diff" << endl;
//				for(int j=0;j<nHelices;j++){
//					cout <<"\t"<< j;
//				}
//				cout <<endl;

//				cout << endl << "SumDiff= " << minimizeVar->Value() << endl;
				
			}
			
			orderSolver.EndSearch();
		}

		orderIndex = sort_indexes(order_);

		return order_;
	}
	vector<size_t> StrTreeGroup::getOrderIndex(){
		return orderIndex;
	}
}


