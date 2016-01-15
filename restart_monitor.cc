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

#include "restart_monitor.h"
#include "misc.h"

namespace operations_research {
	void LNSRestart::FixPositions(){
		char seq[size()];
		int isAssigned[size()];
		for(int i = 1; i < size(); i++) {
			if(vSeq_[i]->Bound()){
				seq[i-1]=ToNucl(vSeq_[i]->Value());
				isAssigned[i-1]=1;
			}
			else{
				seq[i-1]='N';
				isAssigned[i-1]=0;
			}
		}
		seq[size()-1]='\0';
//		cout << seq << endl;

		int* mfeStr[viennaPlugin.size()];
		for(int i=0;i<viennaPlugin.size();i++)	{
			viennaPlugin.at(i)->setSequence(seq);		
			viennaPlugin.at(i)->fold();
			mfeStr[i] = viennaPlugin.at(i)->getBasePairs1Index();
//			cout<< viennaPlugin.at(i)->getStructure() << endl;

		}

		for(int i = 1; i < size(); i++) {			
			if(isAssigned[i-1] == 1){
				int isCorrect = 1;
				for(int j=0;j<viennaPlugin.size();j++){
					//	if(mfeStr[j][i] == -1 || mfeStr[i] !=str_int_[i]){
					if(mfeStr[j][i] !=str_int_.at(j)[i]){
						isCorrect = 0;
						break;
					}					
				}
				if(isCorrect==1){
					fixedPositions[i]=vSeq_[i]->Value();
				}
				else{
					fixedPositions[i]=NOT_ASSIGNED;
				}
			}
			else{
				fixedPositions[i]=NOT_ASSIGNED;
			}
			
		}

		// SAVE ASSIGNMENT AND RESET RESTART COUNTER
		int nDiff=0;
		for(int i = 1; i < size(); i++) {
			if(lastAssignment[i]!=(fixedPositions[i]== NOT_ASSIGNED)){
				lastAssignment[i]=(fixedPositions[i]== NOT_ASSIGNED);
				nRestarts=1;
				nDiff++;
			}
		}

//		ASSIGNMENT BY VALUES
//		for(int i = 1; i < size(); i++) {
//			if(lastAssignment[i]!=fixedPositions[i]){
//				lastAssignment[i]=fixedPositions[i];
//				nRestarts=1;
//				nDiff++;
//			}
//		}

		if(nDiff>0){
			threshold_= initThreshold;
			//cout << "Diff " << nDiff << endl;
		}
		


		// PRINT FIXED  POSITIONS
//		for(int i = 1; i < size(); i++) {
//			if(fixedPositions[i]!=NOT_ASSIGNED){
//				seq[i-1]=ToNucl(fixedPositions[i]);
//			}
//			else{
//				seq[i-1]='N';
//			}
//		}
//		cout << seq << endl<< endl ;

	}
	
	LNSRestart* MakeLNSRestart(Solver *s,  BaseAssignVariables* myDB, const std::vector<IntVar*>& vSeq, vector<int*> str_int, int dangles, std::string rnaLib, std::string energyModel, vector<double> temperatures, double timeLimit, int maxRestarts, StrTreeGroup* treeGroup) {
		return s->RevAlloc(new LNSRestart(s, myDB, vSeq, str_int, dangles,rnaLib,energyModel, temperatures, timeLimit,maxRestarts,treeGroup));
	}
	
	
	void HelixLNSRestart::FixPositions(){
		char seq[size()];
		int isAssigned[size()];
		for(int i = 1; i < size(); i++) {
			if(vSeq_[i]->Bound()){
				seq[i-1]=ToNucl(vSeq_[i]->Value());
				isAssigned[i-1]=1;
			}
			else{
				seq[i-1]='N';
				isAssigned[i-1]=0;
			}
		}
		seq[size()-1]='\0';
//		cout << seq << endl;
		string strSeq(seq);
//		cout<< strSeq << endl;
		int* mfeStr[nHelices];
		bool solvedHelices[nHelices];
		bool boundHelices[nHelices];
		for(int i=0;i<nHelices;i++){
			solvedHelices[i]=false;
			boundHelices[i]=true;
			int iIndex=treeGroup_->getHelix(i)->getI();
			int treeIndex=treeGroup_->getStrId(i);
			string subSeq = strSeq.substr(iIndex-1,treeGroup_->getHelix(i)->getJ()-iIndex+1);
			if(subSeq.find("N")==std::string::npos){
				viennaPlugin.at(i)->setSequence(subSeq.c_str());
				viennaPlugin.at(i)->fold();
				mfeStr[i] = viennaPlugin.at(i)->getBasePairs1Index();

				bool helixSolved = true;
				for(int j = iIndex; j <= treeGroup_->getHelix(i)->getJ(); j++) {
					if((mfeStr[i][j-iIndex+1] == -1 && str_int_.at(treeIndex)[j]!=-1) || (mfeStr[i][j-iIndex+1] != -1 && mfeStr[i][j-iIndex+1] != str_int_.at(treeIndex)[j]-iIndex+1)){
//						cout << "No solved " <<j<< " " <<j-iIndex+1<< " " << mfeStr[i][j-iIndex+1]<<" " <<  str_int_.at(treeIndex)[j] << " " << str_int_.at(treeIndex)[j]-iIndex+1 <<  endl;
						helixSolved=false;
						break;
					}
				}
				solvedHelices[i]=helixSolved;
//				cout<< subSeq << endl << viennaPlugin.at(i)->getStructure() << endl;
			}
			else{
				boundHelices[i]=false;
			}
			
		}
		
		fixedPositions.assign(size(),NOT_ASSIGNED);
		for(int i=0;i<nHelices;i++){
			if(boundHelices[i] && solvedHelices[i]){
				for(int j = treeGroup_->getHelix(i)->getI(); j <= treeGroup_->getHelix(i)->getJ(); j++){
					fixedPositions[j]=vSeq_[j]->Value();
				}
			}
		}


		// SAVE ASSIGNMENT AND RESET RESTART COUNTER
		int nDiff=0;
		
		for(int i = 1; i < size(); i++) {
			if(lastAssignment[i]!=(fixedPositions[i]== NOT_ASSIGNED)){
				lastAssignment[i]=(fixedPositions[i]== NOT_ASSIGNED);
				nRestarts=1;
				nDiff++;
			}
		}
//			ASSIGNMENT BY VALUES
//		for(int i = 1; i < size(); i++) {
//			if(lastAssignment[i]!=fixedPositions[i]){
//				lastAssignment[i]=fixedPositions[i];
//				nRestarts=0;
//				nDiff++;
//			}
//		}

		if(nDiff>0){
			threshold_= initThreshold;
			//cout << "Diff " << nDiff << endl;
		}
		


		// PRINT FIXED  POSITIONS
//		for(int i = 1; i < size(); i++) {
//			if(fixedPositions[i]!=NOT_ASSIGNED){
//				seq[i-1]=ToNucl(fixedPositions[i]);
//			}
//			else{
//				seq[i-1]='N';
//			}
//		}
//		cout << seq << endl<< endl ;
			
	}
	
	HelixLNSRestart* MakeHelixLNSRestart(Solver *s,  BaseAssignVariables* myDB, const std::vector<IntVar*>& vSeq, vector<int*> str_int, int dangles, std::string rnaLib, std::string energyModel, vector<double> temperatures, double timeLimit, int maxRestarts, StrTreeGroup* treeGroup) {
		return s->RevAlloc(new HelixLNSRestart(s, myDB, vSeq, str_int, dangles,rnaLib,energyModel, temperatures, timeLimit,maxRestarts,treeGroup));
	}	
}
