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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "energy_constant.h"
#include "vienna_constraint.h"
#include "minenergy_constraint.h"
#include "minensdef_constraint.h"
#include "energy_constraint.h"
#include "ensdef_constraint.h"
#include "restart_monitor.h"
#include "ifold.h"
#include "diversity_measures.h"


using namespace std;

namespace operations_research {


	IFold::~IFold(){};
	IFold::IFold(Solver* const _solver,int _n) : 
	             solver(_solver),n(_n){
		
		
		vSeq = vector<IntVar *> (n+1);
		vSeqGolomb = vector<IntVar *> (n+1);		
		dangles_=1;
		
		energyModel_=getExecPath("");
		energyModel_.append("/");
		energyModel_.append(TURNER_04_FILE);

		rnaLib_ = VIENNA_LIB;
			
	}

	IFold::IFold(Solver* const _solver,int _n, int dangles, std::string rnaLib, std::string energyModel) : 
	             solver(_solver),n(_n), dangles_(dangles), rnaLib_(rnaLib), energyModel_(energyModel){
		vSeq = vector<IntVar *> (n+1);
		vSeqGolomb = vector<IntVar *> (n+1);		
	}

	void IFold::InitSequenceDomain(char* sequence){
		int i;

		IUPAC['N'] = vector<int> {2,0,1,3}; //All
		IUPAC['A'] = vector<int> {0}; //A
		IUPAC['C'] = vector<int> {1}; //C
		IUPAC['G'] = vector<int> {2}; //G
		IUPAC['U'] = vector<int> {3}; //T
		IUPAC['V'] = vector<int> {2,0,1}; //G,A,C
		IUPAC['B'] = vector<int> {2,3,1}; //G,U,C
		IUPAC['H'] = vector<int> {0,3,1}; //A,U,C
		IUPAC['D'] = vector<int> {2,0,3}; //G,A,U
		IUPAC['K'] = vector<int> {2,3}; //G,U
		IUPAC['S'] = vector<int> {2,1}; //G,C
		IUPAC['W'] = vector<int> {0,3}; //A,U
		IUPAC['M'] = vector<int> {0,1};//A,C
		IUPAC['Y'] = vector<int> {1,3}; //C,T
		IUPAC['R'] = vector<int> {0,2}; //A,G        
		
		if(strlen(sequence) == 0){
			for(i = 0; i < n; i++){
				sequence[i] = 'N';		
			}
		}
		sequence[n] = '\0';

		vSeq[0] = solver->MakeIntConst(-1);
		vSeqGolomb[0] = solver->MakeIntConst(-1);		

		for(i = 1; i <= n; i++){
			vSeq[i] = solver->MakeIntVar(IUPAC[sequence[i-1]], StringPrintf("S_%03d", i));		
			vSeqGolomb[i] =solver->MakeElement(golomb,vSeq[i])->Var();
		}	
		
	}
	
	void IFold::InitDomains(char* sequence){
		InitSequenceDomain(sequence);
	}
	
    void IFold::AddAminoAcidConstraint(AAConstraint* aaConstraint){

		int nCodons= aaConstraint->getLength();		
		if(nCodons>0){
			codonPositions.push_back(std::make_pair(aaConstraint->getStartPos(),nCodons));
			IntVar* const12 = solver->MakeIntConst(12);				
			
			int currentCodons=vCodons.size();
			
			for(int i=0;i<nCodons;i++){

				vCodons.push_back(solver->MakeIntVar(aaConstraint->getDomain(i), StringPrintf("AA_%03d", i+currentCodons)));	

				// Constraint to transform nucleotides sequence to codons ((vSeq[i*3]-vSeq[(i*3)+1] + ((vSeq[i*3]==vSeq[(i*3)+1])*vSeq[i*3]) + 12)+(vSeq[(i*3)+2] * 25))
				IntVar* equalSecondFirst = solver->MakeIsEqualVar(vSeqGolomb[(i*3)+aaConstraint->getStartPos()],vSeqGolomb[((i*3)+1)+aaConstraint->getStartPos()]);						
				
				solver->AddConstraint(solver->MakeEquality(vCodons[currentCodons+i],
				solver->MakeSum(
					solver->MakeSum(
						solver->MakeDifference(vSeqGolomb[(i*3)+aaConstraint->getStartPos()],vSeqGolomb[((i*3)+1)+aaConstraint->getStartPos()])->Var(),
						solver->MakeSum(solver->MakeProd(equalSecondFirst,vSeqGolomb[(i*3)+aaConstraint->getStartPos()])->Var(),const12)->Var())->Var(), 
					solver->MakeProd(vSeqGolomb[((i*3)+2)+aaConstraint->getStartPos()],25)->Var())->Var()
				));

				// Contraints for BLOSUM score maximization
				if(aaConstraint->maximizeScore() && aaConstraint->getTarget().length() > 0){
					maxBlosumValue+=aaConstraint->getMaxBlosumValue(i);
					vAaSeqSimilarity.push_back(solver->MakeElement(aaBlosum[aaConstraint->getBlosumIndex(i)],vCodons[currentCodons+i])->Var());
				}				
			}
		}
	}

	void IFold::AddNucleotideRangeConstraints(double minGCcont, double maxGCcont, int minAU, int maxAU, int minGC, int maxGC, int minGU, int maxGU,std::vector<std::tuple<int,int,int>> listMinA,std::vector<std::tuple<int,int,int>> listMaxA,std::vector<std::tuple<int,int,int>> listMinC,std::vector<std::tuple<int,int,int>> listMaxC,std::vector<std::tuple<int,int,int>> listMinG,std::vector<std::tuple<int,int,int>> listMaxG,std::vector<std::tuple<int,int,int>> listMinU,std::vector<std::tuple<int,int,int>> listMaxU,std::vector<std::tuple<int,int,int>> listConsA,std::vector<std::tuple<int,int,int>> listConsC,std::vector<std::tuple<int,int,int>> listConsG,std::vector<std::tuple<int,int,int>> listConsU){
		int i,j;
		// GC content contraints
		int intMinGCcont = (minGCcont == 0) ? 0 : (int) ceil((n*minGCcont)/100.0);
		int intMaxGCcont = (maxGCcont == 100) ? n : (int) floor((n*maxGCcont)/100.0);
		if(intMinGCcont>0 || intMaxGCcont < n){
			vector<IntVar*> gcCont = vector<IntVar *> (n);
			for(i = 1; i <= n; i++){
				gcCont[i-1] = solver->MakeElement(isGC,vSeq[i])->Var();		
			}
			solver->AddConstraint(solver->MakeBetweenCt(solver->MakeSum(gcCont)->Var(),intMinGCcont,intMaxGCcont));
		}
		
		// Contraints for number of base pairs of each type (applicable to all target structures)
		for(j=0;j<vBPs.size();j++){
			int intMinAU = (minAU == 0) ? 0 : minAU;
			int intMaxAU = (maxAU == -1) ? vBPs.at(j).size(): maxAU;
			if(intMinAU > 0 || intMaxAU < vBPs.at(j).size()){
				vector<IntVar*> numAU = vector<IntVar *> (vBPs.at(j).size());
				for(i = 0; i <vBPs.at(j).size() ; i++){
					numAU[i] = solver->MakeIsEqualCstVar(solver->MakeAbs(vBPs.at(j)[i]),9);		
				}
				solver->AddConstraint(solver->MakeBetweenCt(solver->MakeSum(numAU)->Var(),intMinAU,intMaxAU));
			}

			int intMinGC = (minGC == 0) ? 0 : minGC;
			int intMaxGC = (maxGC == -1) ? vBPs.at(j).size(): maxGC;
			if(intMinGC > 0 || intMaxGC < vBPs.at(j).size()){
				vector<IntVar*> numGC = vector<IntVar *> (vBPs.at(j).size());
				for(i = 0; i <vBPs.at(j).size() ; i++){
					numGC[i] = solver->MakeIsEqualCstVar(solver->MakeAbs(vBPs.at(j)[i]),6);		
				}
				solver->AddConstraint(solver->MakeBetweenCt(solver->MakeSum(numGC)->Var(),intMinGC,intMaxGC));
			}
			
			int intMinGU = (minGU == 0) ? 0 : minGU;
			int intMaxGU = (maxGU == -1) ? vBPs.at(j).size(): maxGU;
			if(intMinGU > 0 || intMaxGU < vBPs.at(j).size()){
				vector<IntVar*> numGU = vector<IntVar *> (vBPs.at(j).size());
				for(i = 0; i <vBPs.at(j).size() ; i++){
					numGU[i] = solver->MakeIsEqualCstVar(solver->MakeAbs(vBPs.at(j)[i]),11);		
				}
				solver->AddConstraint(solver->MakeBetweenCt(solver->MakeSum(numGU)->Var(),intMinGU,intMaxGU));
			}
		}		
		
		// Constraints for maximum and minimum nucleotides
		if(listMinA.size()>0 || listMaxA.size()>0){
			vector<IntVar*> numA = vector<IntVar *> (n);
			for(i = 1; i<=n ; i++){
				numA[i-1] = solver->MakeIsEqualCstVar(vSeq[i],0);		
			}
			for(j=0;j<listMinA.size();j++){
				solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(numA.begin () + get<1>(listMinA[j])-1, numA.begin()+ get<2>(listMinA[j])))->Var(), get<0>(listMinA[j])));
			}			
			for(j=0;j<listMaxA.size();j++){
				solver->AddConstraint(solver->MakeLessOrEqual(solver->MakeSum(std::vector<IntVar*>(numA.begin () + get<1>(listMaxA[j])-1, numA.begin()+ get<2>(listMaxA[j])))->Var(), get<0>(listMaxA[j])));
			}			
		}
		
		if(listMinC.size()>0 || listMaxC.size()>0){
			vector<IntVar*> numC = vector<IntVar *> (n);
			for(i = 1; i<=n ; i++){
				numC[i-1] = solver->MakeIsEqualCstVar(vSeq[i],1);		
			}
			for(j=0;j<listMinC.size();j++){
				solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(numC.begin () + get<1>(listMinC[j])-1, numC.begin()+ get<2>(listMinC[j])))->Var(), get<0>(listMinC[j])));
			}			
			for(j=0;j<listMaxC.size();j++){
				solver->AddConstraint(solver->MakeLessOrEqual(solver->MakeSum(std::vector<IntVar*>(numC.begin () + get<1>(listMaxC[j])-1, numC.begin()+ get<2>(listMaxC[j])))->Var(), get<0>(listMaxC[j])));
			}			
		}
		
		if(listMinG.size()>0 || listMaxG.size()>0){
			vector<IntVar*> numG = vector<IntVar *> (n);
			for(i = 1; i<=n ; i++){
				numG[i-1] = solver->MakeIsEqualCstVar(vSeq[i],2);		
			}
			for(j=0;j<listMinG.size();j++){
				solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(numG.begin () + get<1>(listMinG[j])-1, numG.begin()+ get<2>(listMinG[j])))->Var(), get<0>(listMinG[j])));
			}			
			for(j=0;j<listMaxG.size();j++){
				solver->AddConstraint(solver->MakeLessOrEqual(solver->MakeSum(std::vector<IntVar*>(numG.begin () + get<1>(listMaxG[j])-1, numG.begin()+ get<2>(listMaxG[j])))->Var(), get<0>(listMaxG[j])));
			}			
		}

		if(listMinG.size()>0 || listMaxG.size()>0){
			vector<IntVar*> numG = vector<IntVar *> (n);
			for(i = 1; i<=n ; i++){
				numG[i-1] = solver->MakeIsEqualCstVar(vSeq[i],3);		
			}
			for(j=0;j<listMinG.size();j++){
				solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(numG.begin () + get<1>(listMinG[j])-1, numG.begin()+ get<2>(listMinG[j])))->Var(), get<0>(listMinG[j])));
			}			
			for(j=0;j<listMaxG.size();j++){
				solver->AddConstraint(solver->MakeLessOrEqual(solver->MakeSum(std::vector<IntVar*>(numG.begin () + get<1>(listMaxU[j])-1, numG.begin()+ get<2>(listMaxU[j])))->Var(), get<0>(listMaxG[j])));
			}			
		}
		
		// Constraints for maximum consecutive nucleotides		
		if(listConsA.size() > 0){
			vector<IntVar*> vNoA = vector<IntVar *> (n);
			for(i = 1; i <= n; i++){
				vNoA[i-1] = solver->MakeIsDifferentCstVar(vSeq[i],0);
			}
			for(j=0;j<listConsA.size();j++){
				for(i=(get<1>(listConsA[j]))-1;i<get<2>(listConsA[j])-get<0>(listConsA[j])-1;i++){
					solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(vNoA.begin () + i, vNoA.begin()+ i+get<0>(listConsA[j])+1))->Var(),1));
				}
			}
		}

		if(listConsC.size() > 0){
			vector<IntVar*> vNoC = vector<IntVar *> (n);
			for(i = 1; i <= n; i++){
				vNoC[i-1] = solver->MakeIsDifferentCstVar(vSeq[i],1);
			}
			for(j=0;j<listConsC.size();j++){
				for(i=(get<1>(listConsC[j]))-1;i<get<2>(listConsC[j])-get<0>(listConsC[j])-1;i++){
					solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(vNoC.begin () + i, vNoC.begin()+ i+get<0>(listConsC[j])+1))->Var(),1));
				}
			}
		}

		if(listConsG.size() > 0){
			vector<IntVar*> vNoG = vector<IntVar *> (n);
			for(i = 1; i <= n; i++){
				vNoG[i-1] = solver->MakeIsDifferentCstVar(vSeq[i],2);
			}
			for(j=0;j<listConsG.size();j++){
				for(i=(get<1>(listConsG[j]))-1;i<get<2>(listConsG[j])-get<0>(listConsG[j])-1;i++){
					solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(vNoG.begin () + i, vNoG.begin()+ i+get<0>(listConsG[j])+1))->Var(),1));
				}
			}
		}

		if(listConsU.size() > 0){
			vector<IntVar*> vNoU = vector<IntVar *> (n);
			for(i = 1; i <= n; i++){
				vNoU[i-1] = solver->MakeIsDifferentCstVar(vSeq[i],3);
			}
			for(j=0;j<listConsU.size();j++){
				for(i=(get<1>(listConsU[j]))-1;i<get<2>(listConsU[j])-get<0>(listConsU[j])-1;i++){
					solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSum(std::vector<IntVar*>(vNoU.begin () + i, vNoU.begin()+ i+get<0>(listConsU[j])+1))->Var(),1));
				}
			}
		}
	}

	void IFold::makeBp(int str_index){
		int nBPs=0;
		vector<IntVar *> str_vBPs;
		vector<int> str_BPO;
		vector<int> str_BPC;
		
		for (int i=1; i<=n; i++){
			if(str_int.at(str_index)[i]>i){
				str_BPO.push_back(i);
				str_BPC.push_back(str_int.at(str_index)[i]);
				str_vBPs.push_back(solver->MakeIntVar(bpdoms, StringPrintf("BP_%03d_%03d", str_index,nBPs)));
				// cout << "BP[" << nBPs << "] = Seq["<< i << "] - Seq["<< str_int.at(str_index)[i] << "]" << endl;
				IntVar* bpShift = solver->MakeIntConst(11);
				
				// REPLACED BY GOLOMB
				//solver->AddConstraint(solver->MakeEquality(str_vBPs[nBPs], solver->MakeDifference(solver->MakeElement(golomb,vSeq[i]),solver->MakeElement(golomb,vSeq[str_int.at(str_index)[i]]))->Var()));
				solver->AddConstraint(solver->MakeEquality(str_vBPs[nBPs], solver->MakeDifference(vSeqGolomb[i],vSeqGolomb[str_int.at(str_index)[i]])->Var()));

				solver->AddConstraint(solver->MakeEquality(vSeq[i], solver->MakeElement(bpToUpO,solver->MakeSum(str_vBPs[nBPs],bpShift)->Var())));
				solver->AddConstraint(solver->MakeEquality(vSeq[str_int.at(str_index)[i]], solver->MakeElement(bpToUpC,solver->MakeSum(str_vBPs[nBPs],bpShift)->Var())));				

				solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSquare(str_vBPs[nBPs]),36));
				nBPs++;				
			}
		}
		vBPs.push_back(str_vBPs);
		BPO.push_back(str_BPO);
		BPC.push_back(str_BPC);
	}

	void IFold::AddCompatibilityConstraint(int* comp_str_int){
		if(comp_str_int!= NULL){
			int nBPComp=0;
			vector<IntVar *> vBPComp;
			for (int i=1; i<=n; i++){
				if(comp_str_int[i]>i){
					vBPComp.push_back(solver->MakeIntVar(bpdoms, StringPrintf("BPComp_%03d", nBPComp)));
					// cout << "Compatible BP[" << nBPComp << "] = Seq["<< i << "] - Seq["<< comp_str_int[i] << "]" << endl;

					// REPLACED BY GOLOMB
					// solver->AddConstraint(solver->MakeEquality(vBPComp[nBPComp], solver->MakeDifference(solver->MakeElement(golomb,vSeq[i]),solver->MakeElement(golomb,vSeq[comp_str_int[i]]))->Var()));
					solver->AddConstraint(solver->MakeEquality(vBPComp[nBPComp], solver->MakeDifference(vSeqGolomb[i],vSeqGolomb[comp_str_int[i]])->Var()));


					// REMOVE
					//IntVar* bpShift = solver->MakeIntConst(11);
					//solver->AddConstraint(solver->MakeEquality(vSeq[i], solver->MakeElement(bpToUpO,solver->MakeSum(vBPComp[nBPComp],bpShift)->Var()))); // ??
					//solver->AddConstraint(solver->MakeEquality(vSeq[comp_str_int[i]], solver->MakeElement(bpToUpC,solver->MakeSum(vBPComp[nBPComp],bpShift)->Var())));	//??
					//solver->AddConstraint(solver->MakeGreaterOrEqual(solver->MakeSquare(vBPComp[nBPComp]),36));

					nBPComp++;
					
				}
			}
		}

	}
		
	void IFold::AddIncompatibilityConstraint(vector<pair<int,int>> vIncompBP){
			for(int i=0;i< vIncompBP.size();i++){
					IntVar * vBPIncomp = solver->MakeIntVar(nobpdoms, StringPrintf("IncompBP_%03d", i));
					// REPLACED BY GOLOMB
					//solver->AddConstraint(solver->MakeEquality(vBPIncomp, solver->MakeDifference(solver->MakeElement(golomb,vSeq[vIncompBP[i].first]),solver->MakeElement(golomb,vSeq[vIncompBP[i].second]))->Var()));
					solver->AddConstraint(solver->MakeEquality(vBPIncomp, solver->MakeDifference(vSeqGolomb[vIncompBP[i].first],vSeqGolomb[vIncompBP[i].second])->Var()));

					// REMOVE
					//IntVar* bpShift = solver->MakeIntConst(11);	
					//solver->AddConstraint(solver->MakeEquality(vSeq[vIncompBP[i].first], solver->MakeElement(bpToUpO,solver->MakeSum(vBPIncomp,bpShift)->Var()))); // ??
					//solver->AddConstraint(solver->MakeEquality(vSeq[vIncompBP[i].second], solver->MakeElement(bpToUpC,solver->MakeSum(vBPIncomp,bpShift)->Var())));	//??
					//solver->AddConstraint(solver->MakeLess(solver->MakeSquare(vBPIncomp),36));
			}
	}

	void IFold::AddStructureConstraints(int* _str_int, int* _str_int_undet, int includeDangles, double foldTemp, int cutPoint, int MFEstructure, int showHelices, std::vector<HelixCstr> helixCstrs){
		int strIndex=str_int.size();
		str_int.push_back(_str_int);
		str_int_undet.push_back(_str_int_undet);
		trgFoldTemps.push_back(foldTemp);

		int i;
		cutPoint_=cutPoint;	
//		if(strIndex==0){
			// Base pair constraints
			makeBp(strIndex);
//		}
//		else{
//			AddCompatibilityConstraint(str_int.at(strIndex));
//		}

		// Create tree of helices
		if(includeDangles){
			str_tree.push_back(new StrTree(strIndex,str_int.at(strIndex),n,BPO[strIndex],BPC[strIndex],TH_RNAIFOLD_DANGLES, cutPoint, foldTemp));
		}
		else{
			str_tree.push_back(new StrTree(strIndex,str_int.at(strIndex),n,BPO[strIndex],BPC[strIndex],TH_RNAIFOLD, cutPoint, foldTemp));
		}
		if(showHelices){
			str_tree[str_tree.size()-1]->showTree();
		}
		
		// Add MFE structure, energy and ensemble defect contraints
		for(i = 0; i < str_tree.at(strIndex)->getHelices().size(); i++){
			double minEnergy = NO_ENERGY_LIMIT;
			double maxEnergy = NO_ENERGY_LIMIT;
			double minED = NO_ED_LIMIT;
			double maxED = NO_ED_LIMIT;
			for(int j = 0; j < helixCstrs.size(); j++){
				if(helixCstrs[j].getTreeId()==strIndex && helixCstrs[j].getHelixId()==i){
					switch(helixCstrs[j].getLimitType()){
						case HC_ENERGY_MIN: minEnergy = helixCstrs[j].getLimitValue();
										 break;
						case HC_ENERGY_MAX: maxEnergy = helixCstrs[j].getLimitValue();
										 break;
						case HC_ED_MIN:  minED = helixCstrs[j].getLimitValue();
										 break;
						case HC_ED_MAX:  maxED = helixCstrs[j].getLimitValue();
										 break;
						default: break;
					}
				}
			}
			
			vector<int64> structure(str_tree.at(strIndex)->getHelices()[i]->getJ()-str_tree.at(strIndex)->getHelices()[i]->getI()+2);
			structure[0]=str_tree.at(strIndex)->getHelices()[i]->getJ()-str_tree.at(strIndex)->getHelices()[i]->getI()+1;
			int strCutPoint = (cutPoint==-1 || cutPoint<str_tree.at(strIndex)->getHelices()[i]->getI() || cutPoint>str_tree.at(strIndex)->getHelices()[i]->getJ()-1)? -1 :(cutPoint-str_tree.at(strIndex)->getHelices()[i]->getI()+2);
			if(str_int_undet.at(strIndex) == NULL){
				for(int j = str_tree.at(strIndex)->getHelices()[i]->getI(); j <= str_tree.at(strIndex)->getHelices()[i]->getJ(); j++){
					if(str_int.at(strIndex)[j]==-1){
						structure[j-str_tree.at(strIndex)->getHelices()[i]->getI()+1]=-1;
					}
					else{
						structure[j-str_tree.at(strIndex)->getHelices()[i]->getI()+1]=str_int.at(strIndex)[j]-str_tree.at(strIndex)->getHelices()[i]->getI()+1;
					}
				}
				// MFE constraint
				if(MFEstructure){
					solver->AddConstraint(solver->RevAlloc(new ViennaConstraintDet(solver, std::vector<IntVar*>(vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getI(), 
																										  vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getJ()+1),
																					 structure,dangles_,rnaLib_,energyModel_,foldTemp, strCutPoint, str_tree.at(strIndex)->getHelices()[i]->getPosLeft()+1, str_tree.at(strIndex)->getHelices()[i]->getJ()-str_tree.at(strIndex)->getHelices()[i]->getI()-str_tree.at(strIndex)->getHelices()[i]->getPosRight()+1,
																					 str_tree.at(strIndex)->getHelices()[i]->getTryLeft() == 1 ? vSeq.at(str_tree.at(strIndex)->getHelices()[i]->getI()-1) : NULL, str_tree.at(strIndex)->getHelices()[i]->getTryRight() == 1 ? vSeq.at(str_tree.at(strIndex)->getHelices()[i]->getJ()+1) : NULL,
																					 maxEnergy, minEnergy)));
				}
				else{
					if(minEnergy != NO_ENERGY_LIMIT || maxEnergy != NO_ENERGY_LIMIT){
						solver->AddConstraint(solver->RevAlloc(new EnergyConstraint(solver, std::vector<IntVar*>(vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getI(), vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getJ()+1), dangles_,rnaLib_,energyModel_,structure, foldTemp, strCutPoint, minEnergy, maxEnergy)));
					}
				}
				if(minED != NO_ED_LIMIT || maxED != NO_ED_LIMIT){
					solver->AddConstraint(solver->RevAlloc(new EnsDefConstraint(solver, std::vector<IntVar*>(vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getI(), vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getJ()+1), dangles_,rnaLib_,energyModel_,structure, foldTemp, strCutPoint, minED, maxED, str_tree.at(strIndex)->getHelices()[i]->getI())));
				}
			}
			else{
				for(int j = str_tree.at(strIndex)->getHelices()[i]->getI(); j <= str_tree.at(strIndex)->getHelices()[i]->getJ(); j++){
					if(str_int_undet.at(strIndex)[j]<0){
						structure[j-str_tree.at(strIndex)->getHelices()[i]->getI()+1]=str_int_undet.at(strIndex)[j];
					}
					else{
						structure[j-str_tree.at(strIndex)->getHelices()[i]->getI()+1]=str_int_undet.at(strIndex)[j]-str_tree.at(strIndex)->getHelices()[i]->getI()+1;
					}
				}
				// MFE constraint
				if(MFEstructure){					
					solver->AddConstraint(solver->RevAlloc(new ViennaConstraintUndet(solver, std::vector<IntVar*>(vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getI(), 
																										  vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getJ()+1),
																					 structure,dangles_,rnaLib_,energyModel_,foldTemp, strCutPoint, str_tree.at(strIndex)->getHelices()[i]->getPosLeft()+1, str_tree.at(strIndex)->getHelices()[i]->getJ()-str_tree.at(strIndex)->getHelices()[i]->getI()-str_tree.at(strIndex)->getHelices()[i]->getPosRight()+1,
																					 str_tree.at(strIndex)->getHelices()[i]->getTryLeft() == 1 ? vSeq.at(str_tree.at(strIndex)->getHelices()[i]->getI()-1) : NULL, str_tree.at(strIndex)->getHelices()[i]->getTryRight() == 1 ? vSeq.at(str_tree.at(strIndex)->getHelices()[i]->getJ()+1) : NULL,
																					 maxEnergy, minEnergy)));
				}
				else{
					if(minEnergy != NO_ENERGY_LIMIT || maxEnergy != NO_ENERGY_LIMIT){
						solver->AddConstraint(solver->RevAlloc(new EnergyConstraint(solver, std::vector<IntVar*>(vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getI(), vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getJ()+1), dangles_,rnaLib_,energyModel_,structure, foldTemp, strCutPoint, minEnergy, maxEnergy)));
					}
				}
				if(minED != NO_ED_LIMIT || maxED != NO_ED_LIMIT){
					solver->AddConstraint(solver->RevAlloc(new EnsDefConstraint(solver, std::vector<IntVar*>(vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getI(), vSeq.begin () + str_tree.at(strIndex)->getHelices()[i]->getJ()+1), dangles_,rnaLib_,energyModel_,structure, foldTemp, strCutPoint, minED, maxED, str_tree.at(strIndex)->getHelices()[i]->getI())));
				}
			}
		}
		
	}
	void IFold::AddLocalStructureConstraints(double foldTemp, int cutPoint, std::vector<LocalCstr> localCstrs){
		std::vector<int> visited;
		for(int i = 0; i < localCstrs.size(); i++){
			if(find(visited.begin(), visited.end(), i) == visited.end()){
				double minEnergy = NO_ENERGY_LIMIT;
				double maxEnergy = NO_ENERGY_LIMIT;
				double minED = NO_ED_LIMIT;
				double maxED = NO_ED_LIMIT;
				int MFEactive = 0;
				visited.push_back(i);
				switch(localCstrs[i].getLimitType()){
					case LC_ENERGY_MIN: minEnergy = localCstrs[i].getLimitValue();
									 break;
					case LC_ENERGY_MAX: maxEnergy = localCstrs[i].getLimitValue();
									 break;
					case LC_ED_MIN:  minED = localCstrs[i].getLimitValue();
									 break;
					case LC_ED_MAX:  maxED = localCstrs[i].getLimitValue();
									 break;
					case LC_MFE:     MFEactive = 1;
									 break;
					default: break;
				}
				
				for(int j = i+1; j < localCstrs.size(); j++){
					if(localCstrs[i].isPairedConstraint(localCstrs[j]) && find(visited.begin(), visited.end(), j) == visited.end()){
						visited.push_back(j);
						switch(localCstrs[j].getLimitType()){
							case LC_ENERGY_MIN: minEnergy = localCstrs[j].getLimitValue();
											 break;
							case LC_ENERGY_MAX: maxEnergy = localCstrs[j].getLimitValue();
											 break;
							case LC_ED_MIN:  minED = localCstrs[j].getLimitValue();
											 break;
							case LC_ED_MAX:  maxED = localCstrs[j].getLimitValue();
											 break;
							case LC_MFE:     MFEactive = 1;
											 break;
							default: break;
						}
					}
				}
				
				// Add compatibility constraints for the local structure (solutions must be compatible with the target local structure)
				vector<IntVar *> vBPComp;
				int nBPComp=0;
				for (int j=1; j<=localCstrs[i].getStructureStr().length(); j++){
					if(localCstrs[i].getStructureArr()[j]>j){
						vBPComp.push_back(solver->MakeIntVar(bpdoms, StringPrintf("LocalBPComp_%03d", nBPComp)));
						solver->AddConstraint(solver->MakeEquality(vBPComp[nBPComp], solver->MakeDifference(vSeqGolomb[localCstrs[i].getStartPos()+j-1],vSeqGolomb[localCstrs[i].getStartPos()+localCstrs[i].getStructureArr()[j]-1])->Var()));
						//cout << "BP Constraint" << (localCstrs[i].getStartPos()+j-1) <<"-"<<(localCstrs[i].getStartPos()+localCstrs[i].getStructureArr()[j]-1) << "  " << localCstrs[i].getStructureArr()[j]<< "-"<< j <<endl;
						nBPComp++;
					}
				}
				
				// Adjust cutpoint to the local structure
				int strCutPoint = (cutPoint==-1 || cutPoint<localCstrs[i].getStartPos() || cutPoint>localCstrs[i].getStartPos() + localCstrs[i].getStructureStr().length())? -1 :(cutPoint-localCstrs[i].getStartPos()+2);

				// If MFE constraint is active energy constraints are included in the Vienna Constraint
				if(MFEactive){
					if(localCstrs[i].hasUndet()){
						solver->AddConstraint(solver->RevAlloc(new ViennaConstraintUndet(solver, std::vector<IntVar*>(vSeq.begin () + localCstrs[i].getStartPos(), 
																											  vSeq.begin () + localCstrs[i].getStartPos() + localCstrs[i].getStructureStr().length()),
																						 localCstrs[i].getStructureArr(),dangles_,rnaLib_,energyModel_,foldTemp, strCutPoint, 1, localCstrs[i].getStructureStr().length(),NULL,NULL,maxEnergy, minEnergy)));
					}
					else{
						solver->AddConstraint(solver->RevAlloc(new ViennaConstraintDet(solver, std::vector<IntVar*>(vSeq.begin () + localCstrs[i].getStartPos(), 
																											  vSeq.begin () + localCstrs[i].getStartPos() + localCstrs[i].getStructureStr().length()),
																						 localCstrs[i].getStructureArr(),dangles_,rnaLib_,energyModel_,foldTemp, strCutPoint, 1, localCstrs[i].getStructureStr().length(),NULL,NULL,maxEnergy, minEnergy)));
					}
				}
				// If MFE constraint is inactive energy constraint is independent
				else{
					if(localCstrs[i].hasUndet()){
						cout << "Undefined pairing positions are allowed only if MFE contraint is active, otherwise there is no target structure to compute the energy for."<< endl;
						exit(1);
					}
					if(minEnergy != NO_ENERGY_LIMIT || maxEnergy != NO_ENERGY_LIMIT){
						solver->AddConstraint(solver->RevAlloc(new EnergyConstraint(solver, std::vector<IntVar*>(vSeq.begin () + localCstrs[i].getStartPos(), vSeq.begin () + localCstrs[i].getStartPos() + localCstrs[i].getStructureStr().length()), dangles_,rnaLib_,energyModel_,localCstrs[i].getStructureArr(), foldTemp, strCutPoint, minEnergy, maxEnergy)));
					}
				}
				
				// Ad ensemble defect constraints
				if(minED != NO_ED_LIMIT || maxED != NO_ED_LIMIT){
					solver->AddConstraint(solver->RevAlloc(new EnsDefConstraint(solver, std::vector<IntVar*>(vSeq.begin () + localCstrs[i].getStartPos(), vSeq.begin () + localCstrs[i].getStartPos() + localCstrs[i].getStructureStr().length()), dangles_,rnaLib_,energyModel_,localCstrs[i].getStructureArr(), foldTemp, strCutPoint, minED, maxED, 1)));
				}

			}			
		}		
		return;
	}
	
	void IFold::SetSearchHeuristic(int helixHeuristic, int varHeuristic, int randomAssignment, int upthreshold, int bpthreshold){
		tree_group = new StrTreeGroup(&str_tree);

		vector<int64> order = tree_group->optimizeOrder(helixHeuristic);
//		cout << "Final order: " << endl;
//		for(int i=0;i<order.size();i++){
//			cout << "\t" <<order[i];
//		}
//		cout << endl;
//		string a;
//		cin >> a;

		vector<vector<int> > BPtypes = tree_group->optimizeValueHeuristic(n,BPO,BPC,str_int, trgFoldTemps);
		vector<vector<int> > UPtypes = tree_group->optimizeUPHeuristic(n,BPO,BPC,str_int);


//		if(BPtypes.size()>1){
//			varHeuristic = SH_MULTI_STR;
//		}

		upthreshold_=upthreshold;
		bpthreshold_=bpthreshold;
		switch(varHeuristic){
			// No search heuristic, worst case
			case SH_NONE:{
				for(int i=1;i<=n;i++){
					vars.push_back(vSeq[i]);
					var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);
					var_str.push_back(0);
				}
			}
			break;
			
			//	Search by tree levels			
			case SH_BOTTOM_TO_TOP:{
				vector<size_t> orderIndex = tree_group->getOrderIndex();
				
				// Get helices ordered by levels, from highest to lowest
//				vector<CPHelix* > sortedHelices = str_tree.at(tree_group->getStrId(i))->getSortedHelices();

				for(size_t i :orderIndex){
					vector<int> bps = tree_group->getHelix(i)->getBPs();
					//Add base pairing positions from inside to outside
					for(int j=bps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(tree_group->getStrId(i))[bps[j]]);
//						cout << "Adding vBP["<<bps[j]<<"]"<< endl;
						if(tree_group->getHelix(i)->isStack(bps[j])){
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP_STACK);
						}
						else{
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);	
						}
						var_str.push_back(tree_group->getStrId(i));
					}

					vector<int> ups = tree_group->getHelix(i)->getUPs();
					//Add unpaired positions
					for(int j=0; j<ups.size();j++){
						vars.push_back(vSeq[ups[j]]);
//						cout << "Adding vSeq["<<ups[j]<<"]"<< endl;						
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);
						var_str.push_back(tree_group->getStrId(i));
						
					}


					vector<int> cbps = tree_group->getHelix(i)->getClosingBPs();
					//Add closing base pairs
					for(int j=cbps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(tree_group->getStrId(i))[cbps[j]]);
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);	
						var_str.push_back(tree_group->getStrId(i));
					}



/*
					vector<int> positions=tree_group->getHelix(i)->getPositions();

					//Add base pairing positions from inside to outside
					for(int j=positions.size()-1; j>=0;j--){
						if(str_int[positions[j]]!=-1 && str_int[positions[j]]>positions[j]){
							vars.push_back(vSeq[positions[j]]);
							vars.push_back(vSeq[str_int[positions[j]]]);
						}
					}

					//Add unpaired positions
					for(int j=0; j<positions.size();j++){
						if(str_int[positions[j]]==-1){
							vars.push_back(vSeq[positions[j]]);
						}
					}
*/
				}
			}		
			break;



			//	Search for multiple structures
			case SH_MULTI_STR:{
				vector<size_t> orderIndex = tree_group->getOrderIndex();
				for(size_t i :orderIndex){
					vector<int> bps = tree_group->getHelix(i)->getBPs();
					//Add base pairing positions from inside to outside
					for(int j=bps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(tree_group->getStrId(i))[bps[j]]);
						if(BPtypes.at(tree_group->getStrId(i))[bps[j]] == VH_MULTI_BP_1){
							if(tree_group->getHelix(i)->isStack(bps[j])){
								var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP_STACK);
							}
							else{
								var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);	
							}
						}
						else  if(BPtypes.at(tree_group->getStrId(i))[bps[j]] == VH_MULTI_BP_0){
							if(tree_group->getHelix(i)->isStack(bps[j])){
								var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP_STACK_INV);
							}
							else{
								var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_MULTI_BP_0);	
							}
						}
						else{
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : BPtypes.at(tree_group->getStrId(i))[bps[j]]);
						}
						var_str.push_back(tree_group->getStrId(i));
					}



					vector<int> ups = tree_group->getHelix(i)->getUPs();
					int danglePos = tree_group->getHelix(i)->getPosLeft()+tree_group->getHelix(i)->getPosRight();
					//Add unpaired positions
					for(int j=0; j<ups.size()-danglePos;j++){
						vars.push_back(vSeq[ups[j]]);
//						cout << "Adding vSeq["<<ups[j]<<"]"<< endl;
						var_types.push_back(randomAssignment != 0 ?	VH_RANDOM : UPtypes.at(tree_group->getStrId(i))[ups[j]]);
						var_str.push_back(tree_group->getStrId(i));
						
					}

					//Add dangles
					for(int j=ups.size()-danglePos; j<ups.size();j++){
						vars.push_back(vSeq[ups[j]]);
//						cout << "Adding dangle vSeq["<<ups[j]<<"]"<< endl;
						// If assigning the three prime
						if(danglePos == 2 && j==ups.size()-1){
							var_types.push_back(randomAssignment != 0 ?	VH_RANDOM : VH_THREEP_UP);
						}
						else{
							var_types.push_back(randomAssignment != 0 ?	VH_RANDOM : UPtypes.at(tree_group->getStrId(i))[ups[j]]);
						}
						var_str.push_back(tree_group->getStrId(i));
						
					}
					vector<int> cbps = tree_group->getHelix(i)->getClosingBPs();
					//Add closing base pairs
					for(int j=cbps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(tree_group->getStrId(i))[cbps[j]]);
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_MULTI_BP_1);
						var_str.push_back(tree_group->getStrId(i));
					}

				}
			}		
			break;









			
			//	Search by tree levels			
			case SH_BOTTOM_TO_TOP_BAK:{

				// Get helices ordered by levels, from highest to lowest
				vector<CPHelix* > sortedHelices = str_tree.at(0)->getSortedHelices();

				for(int i=0; i <sortedHelices.size();i++){					
					vector<int> bps = sortedHelices[i]->getBPs();
					//Add base pairing positions from inside to outside
					for(int j=bps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(0)[bps[j]]);
//						cout << "Adding vBP["<<bps[j]<<"]"<< endl;
						if(sortedHelices[i]->isStack(bps[j])){
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP_STACK);
						}
						else{
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);	
						}
						var_str.push_back(tree_group->getStrId(0));
					}

					vector<int> cbps = sortedHelices[i]->getClosingBPs();
					//Add closing base pairs
					for(int j=cbps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(0)[cbps[j]]);
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);	
						var_str.push_back(tree_group->getStrId(0));
					}


					vector<int> ups = sortedHelices[i]->getUPs();
					//Add unpaired positions
					for(int j=0; j<ups.size();j++){
						vars.push_back(vSeq[ups[j]]);
//						cout << "Adding vSeq["<<ups[j]<<"]"<< endl;						
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);
						var_str.push_back(tree_group->getStrId(0));
						
					}

/*
					vector<int> positions=sortedHelices[i]->getPositions();

					//Add base pairing positions from inside to outside
					for(int j=positions.size()-1; j>=0;j--){
						if(str_int[positions[j]]!=-1 && str_int[positions[j]]>positions[j]){
							vars.push_back(vSeq[positions[j]]);
							vars.push_back(vSeq[str_int[positions[j]]]);
						}
					}

					//Add unpaired positions
					for(int j=0; j<positions.size();j++){
						if(str_int[positions[j]]==-1){
							vars.push_back(vSeq[positions[j]]);
						}
					}
*/
				}
			}		
			break;

			case SH_BOTTOM_TO_TOP_UP:{
				// Get helices ordered by levels, from highest to lowest
				vector<CPHelix* > sortedHelices = str_tree.at(0)->getSortedHelices();

				for(int i=0; i <sortedHelices.size();i++){
					vector<int> bps = sortedHelices[i]->getBPs();
					//Add base pairing positions from inside to outside
					for(int j=bps.size()-1; j>=0;j--){
						vars.push_back(vBPs.at(0)[bps[j]]);
						//cout << "Adding vBP["<<bps[j]<<"]"<< endl;
						if(sortedHelices[i]->isStack(bps[j])){
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP_STACK);
						}
						else{
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);	
						}
						var_str.push_back(tree_group->getStrId(0));
					}
					vector<int> ups = sortedHelices[i]->getUPs();

					//Add unpaired positions
					for(int j=0; j<ups.size();j++){
						vars.push_back(vSeq[ups[j]]);
						//cout << "Adding vSeq["<<ups[j]<<"]"<< endl;						
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);
						var_str.push_back(tree_group->getStrId(0));
					}
				}
			}
			break;


			// Default, from inner to outer stems
			case SH_IN_TO_OUT:
			default:				
				vector<bool> used(n+1,false);
				for (int d = SIGMA+1; d <= n; d++){
					for(int i=1; (i+d<=n); i++){
						int j=i+d;
						if(str_int.at(0)[i]==j){
							for(int l=0;l<BPO.at(0).size();l++){
								if(BPO.at(0)[l]==i){
									vars.push_back(vBPs.at(0)[l]);
									if(str_int.at(0)[i+1]==j-1){
										var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP_STACK);
									}
									else{
										var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_BP);
									}
									var_str.push_back(tree_group->getStrId(0));
									
									used[i] = true;
									used[j] = true;
									for(int k=i+1;k<j;k++){
										if(!used[k]&&str_int.at(0)[k]==-1){
											vars.push_back(vSeq[k]);
											var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);
											var_str.push_back(tree_group->getStrId(0));
											used[k] = true;
										}
									}
								}
							}
/*							vars.push_back(vSeq[i]);
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);

							vars.push_back(vSeq[j]);
							var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);							
							used[i] = true;
							used[j] = true;
							for(int k=i+1;k<j;k++){
								if(!used[k]&&str_int[k]==-1){
									vars.push_back(vSeq[k]);
									var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);									
									used[k] = true;
								}
							}*/
						}
					}
				}
				// Set remaining positions (external loops)
				for(int i=1;i<=n;i++){
					if(!used[i]){
						vars.push_back(vSeq[i]);
						var_types.push_back(randomAssignment != 0 ? VH_RANDOM : VH_UP_CONTAINS);
						var_str.push_back(tree_group->getStrId(0));
					}
				}
			          
			break;
		}
		
	}
	void IFold::Search(int maxSolutions, int64 timeLimit, int minimizeMFE, int minimizeEnsDef, int LNS, int LNSunchangedRestarts, int LNSrestartTime, int showMeasures){
		int i;

		BaseAssignVariables* const db = MakePhase(solver,vars,var_types, upthreshold_, bpthreshold_);
		std::vector<SearchMonitor*> monitors;
		SolutionCollector* collector;
		
		// UNCOMMENT TO ACTIVATE LOG MONITOR
		//SearchMonitor* const log = solver->MakeSearchLog(100000);
		//monitors.push_back(log);

		//SearchMonitor* const log = solver->MakeSearchTrace("EEEEEEEEEE");
		//monitors.push_back(log);

		// UNCOMMENT TO PRINT SEARCH TREE
		//SearchMonitor* const cpviz = solver->MakeTreeMonitor(vars, "tree.xml","visualization.xml");
		//monitors.push_back(cpviz);

		// Variables for diversity measures
		vector<RNAPlugin*> measuresVienna;
		vector<string> targetStructures;
		vector<int*> str_int_0index;
		
		for(i=0;i<str_int.size();i++){
			string targetStr="";
			measuresVienna.push_back(newRNAplugin(rnaLib_,n,dangles_));
			measuresVienna.at(i)->setEnergyModel(energyModel_);
			measuresVienna.at(i)->setCutPoint(cutPoint_);
			measuresVienna.at(i)->setTemperature(trgFoldTemps.at(i));
			str_int_0index.push_back((int*) malloc(sizeof(int)*(n)));
			targetStructures.push_back(targetStr);
		}

		//if(showMeasures>0){
			for(int j=0;j<str_int.size();j++){
				for(i=1;i<=n;i++){
					if(str_int.at(j)[i]==-1){
						str_int_0index.at(j)[i-1]=-1;
						targetStructures.at(j).append(".");
					}
					else{
						str_int_0index.at(j)[i-1]=str_int.at(j)[i]-1;
						if(str_int.at(j)[i]>i){
							targetStructures.at(j).append("(");
						}
						else{
							targetStructures.at(j).append(")");
						}
					}
				} 
				measuresVienna.at(j)->setTestStructure(targetStructures.at(j));
			}
		//}
		
		IntVar* blosumScore;

		// Maximize blosum similarity monitor
		if(maxBlosumValue > 0){
			blosumScore = solver->MakeSum(vAaSeqSimilarity)->Var();
			OptimizeVar* opBlosumScore;
			opBlosumScore = solver->MakeMaximize(blosumScore,1);
			monitors.push_back(opBlosumScore);
		}

		// LNS restart monitor 
		SearchMonitor* restartMonitor;
		if(LNS>0){
			restartMonitor=MakeLNSRestart(solver,db, vSeq, str_int, dangles_,rnaLib_,energyModel_,trgFoldTemps, LNSrestartTime, LNSunchangedRestarts,tree_group);
			monitors.push_back(restartMonitor);
		}


		// Time limit monitor
		SearchLimit* limit = NULL;
		//clock_t begin_search = clock();

		if(timeLimit>0){
			// Search limit.
			limit = solver->MakeTimeLimit(timeLimit);
			monitors.push_back(limit);
		}

		// Energy minimization constraint  (for all structures)
		if(minimizeMFE==1){
			for(i =0; i< str_int.size();i++){
				MinEnergyConstraint* mfeConst  = new MinEnergyConstraint(solver, std::vector<IntVar*>(vSeq.begin () +1,vSeq.end()),dangles_,rnaLib_,energyModel_,((str_int_undet.at(i)==NULL) ? targetStructures.at(i) : "") ,trgFoldTemps.at(i), cutPoint_);
				solver->AddConstraint(solver->RevAlloc(mfeConst));
				monitors.push_back(mfeConst->obj);
			}
		}

		// Ensemble defect minimization constraint  (for all structures)
		if(minimizeEnsDef==1){
			for(i =0; i< str_int.size();i++){
				MinEnsDefConstraint* minEdConst  = new MinEnsDefConstraint(solver, std::vector<IntVar*>(vSeq.begin () +1,vSeq.end()),dangles_,rnaLib_,energyModel_,((str_int_undet.at(i)==NULL) ? str_int_0index.at(i) : NULL), trgFoldTemps.at(i), cutPoint_);
				solver->AddConstraint(solver->RevAlloc(minEdConst));
				monitors.push_back(minEdConst->obj);
			}
		}

		int nSolutions =0;
		solver->NewSearch(db,monitors);

		while (solver->NextSolution() && (maxSolutions==0 || nSolutions<maxSolutions)) {
			nSolutions++;
			
			
			string outSequence ="";
			for(i=1;i<=n;i++)
				outSequence+=ToNucl(vSeq[i]->Value());
				
			cout << outSequence << endl;

			if(vCodons.size()>0){
				int codonIndex=0;
				for(i=0;i<codonPositions.size();i++){
					cout << "Amino acid sequence at "<< codonPositions[i].first<<": ";
					for(int j=0;j<codonPositions[i].second;j++){
						cout << codonAaDict.at(vCodons[codonIndex]->Value());
						codonIndex++;
					}
					cout << endl;						
				}
			}
			
			if(maxBlosumValue > 0){
				cout << "BLOSUM score:" << blosumScore->Value() << " of " << maxBlosumValue << endl;
			}

			// Print structural diversity measures
			if(showMeasures>0){
				for(int i=0;i<str_int.size();i++){
					cout << getNtContent(str_int.at(i), vSeqGolomb, n) << endl;
					double energy =0;
					measuresVienna.at(i)->setSequence(outSequence);
					int* mfeStr=NULL;
					if(str_int_undet.at(i)!=NULL){
						if(cutPoint_==1){
							energy=measuresVienna.at(i)->fold();
						}
						else{
							energy=measuresVienna.at(i)->cofold();
						}
						measuresVienna.at(i)->setTestStructure(measuresVienna.at(i)->getStructure());
						cout << measuresVienna.at(i)->getStructure() << endl;
						mfeStr= measuresVienna.at(i)->getBasePairs();
					}
					else{
						energy=measuresVienna.at(i)->energyOfStruct();
					}
				
					
					 

					double **outBpPr;
					if(cutPoint_==-1){
						outBpPr=measuresVienna.at(i)->basePairProbsMatrix();
					}
					else{
						outBpPr=measuresVienna.at(i)->basePairProbsMatrixCofold();					
					}
					cout << "Energy of the structure:" << energy << endl;
					cout << "Probability of MFE structure:" << measuresVienna.at(i)->probOfStruct() << endl;
					cout << "Expected pointwise entropy:" << expectedPointwiseEntropy(outBpPr,n) <<endl;
					cout << "Morgan-Higgs structural Diversity:" << morganHiggsStructuralDiversity(outBpPr,n) <<endl;
					cout << "Vienna structural Diversity:" << viennaStructuralDiversity(outBpPr,n) <<endl;
					cout << "Expected base pair distance:" << expectedBpDistance(outBpPr,n,mfeStr==NULL ? str_int_0index.at(i): mfeStr) <<endl;
					cout << "Ensemble defect:" << ensembleDefect(outBpPr,n,mfeStr==NULL ? str_int_0index.at(i): mfeStr) <<endl;

					for(int j=0;j<n;j++){
						free(outBpPr[j]);
					}
					free(outBpPr);
					if(mfeStr!=NULL){
						free(mfeStr);
					}
				}
			}
		} 
				
		solver->EndSearch();	  


		if(limit!=NULL && limit->crossed()){
			if (nSolutions == 0) {
				cout << "No solution found." << endl;
			}
		}
		else{
			if (nSolutions == 0) {
				cout << "No possible solution." << endl;				
			}
			else if (maxSolutions==0 || nSolutions < maxSolutions) {
				if(maxBlosumValue == 0){
					cout << "All possible solutions found!"<< endl;
				}
				else{
					cout << "Optimal solution found!"<< endl;
				}
			} 
		}
		
	}	


// PRIVATE FUNCTIONS


	vector<CPHelix* > IFold::getHelices(){
		return str_tree.at(0)->getHelices();
	}

	vector<IntVar* > IFold::getSearchVars(){
		return vars;
	}



	/**************************************** ENERGY CONSTRAINT FUNCTIONS FOR FIXED POSITIONS  ************************************/

	bool IFold::isFixed(int i, int j){
		for(int pos=i;pos<=j;pos++){
			if(!vSeq[pos]->Bound()){
				return false;
			}

		}
		return true;
	}

	int64 IFold::GetTableValue(const std::vector<int64>& table, const std::vector<int64>& indexes){
		int64 index=0;
		for(int i=0;i< indexes.size();i++){
			index+=indexes[i]*coeffs[i];
		}
		return table[index];
	}

}
