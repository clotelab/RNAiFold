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

#include "energy_constant.h"
#include "strtree.h"

extern int hairpin37[31];
int max_bp_energy = -330;

namespace operations_research {
	StrTree::~StrTree(){};
	StrTree::StrTree(int tree_id, int* str_int,int n, std::vector<int> BPO, std::vector<int> BPC, int treeType, int cutPoint, double temp) : tree_id_(tree_id), temp_(temp) {
		if(treeType==TH_RNAIFOLD){
			createRNAiFoldHelices(str_int,n,BPO,BPC,false, cutPoint);
		}
		else{
			createRNAiFoldHelices(str_int,n,BPO,BPC,true, cutPoint);
		}
		// DEBUG 
		//for(int i = 0; i < helices.size(); i++){	
		//	helices[i]->toString();
		//}

	}

	//******************************** BEGIN RNAIFOLD HELICES *********************************************//

	void StrTree::createRNAiFoldHelices(int* str_int,int n, std::vector<int> BPO, std::vector<int> BPC, bool includeDangles, int cutPoint){
		int helixId = 0;
		int i,j;
		int pos,c;

		int lbp=-1;
		int level;
		int bps = 0;

		// Store full sequence as first helix
		i = 1;
		j = n;
	
		// Create root helix 
		level = 0;
		helices.push_back(new CPHelix(helixId,i,j));
		helices[helixId]->setLevel(level);

		// Initialize variables and start new helix
		helixId++;

		for(pos=1; pos<=n; pos++){
			if(pos < str_int[pos]){
				if(lbp == -1){	
					// Start first helix	
					i = pos;
					j = str_int[pos];
				}else if(!consecutive(lbp, str_int[lbp],pos,str_int[pos])){
						// create new helix
						helices.push_back(new CPHelix(helixId,i,j));

						// Find Parent and set level
						level = 1;
						for (c=helixId-1; c>=0;c--){
							if(helices[c]->contains(helices[helixId])){
								level = helices[helices[c]->getId()]->getLevel()+1; 
								helices[helices[c]->getId()]->addInside(helixId);
								helices[helixId]->setParent(helices[c]->getId());				
								break;  // Only in the nearest helix
							}
						}
						helices[helixId]->setLevel(level);

			
						// Initialize variables and increase helix counter
						helixId++;
						i = pos;
						j= str_int[pos];
						bps =0;
				}
				lbp = pos;
				bps++;
			}
		}
	
		if(lbp != -1){	
				helices.push_back(new CPHelix(helixId,i,j));

				// Find Parent and set level
				level = 1;
				for (c=helixId-1; c>=0;c--){
					if(helices[c]->contains(helices[helixId])){
						level = helices[helices[c]->getId()]->getLevel()+1; 
						helices[helices[c]->getId()]->addInside(helixId);
						helices[helixId]->setParent(helices[c]->getId());				
						break;  // Only in the nearest helix
					}
				}
				helices[helixId]->setLevel(level);

	
				// Initialize variables and increase helix counter
				helixId++;
				i = pos;
				j= str_int[pos];
				bps =0;
		}
		if(includeDangles){
			// Add dangling positions except in the following cases:
			//   1.- Opening base pair is at the first position of the structure/strand
			//   2.- Closing base pair is at the last position of the structure/strand
			//   3.- Dangling position is paired (MAX_LEFT_MISMATCH)
			for(i=helixId-1;i>0;i--){
				if(helices[i]->getI()>1 &&  (str_int[helices[i]->getI()-1]==-1) && (cutPoint!=helices[i]->getI()-1)){
					helices[i]->setPosLeft(1);
					helices[i]->setI(helices[i]->getI()-1);
				}
				if(helices[i]->getJ()<n &&  (str_int[helices[i]->getJ()+1]==-1) && (cutPoint!=helices[i]->getJ()+1)){
					helices[i]->setPosRight(1);
					helices[i]->setJ(helices[i]->getJ()+1);
				}			
			}
		}
		else{
			for(i=helixId-1;i>0;i--){
				if(helices[i]->getI()>1 &&  (str_int[helices[i]->getI()-1]==-1) && (cutPoint!=helices[i]->getI()-1)){
					helices[i]->setTryLeft(1);
				}
				if(helices[i]->getJ()<n &&  (str_int[helices[i]->getJ()+1]==-1) && (cutPoint!=helices[i]->getJ()+1)){
					helices[i]->setTryRight(1);
				}			
			}
		}


		for(i=0; i<helices.size(); i++){
			helices[i]->setHelixPositions(helices);
			helices[i]->setPositionIndexAndType(str_int, BPO,BPC);
			// helices[i]->toString();
		}
		initLeavesDistance();
		
		//*********************** POST-PROCESS ******************************//
		// Post-processing for fusion helixes, find helices with pattern (.XXXXXX) or (XXXXXX.) or 
		int post_process=1;
		if(post_process){
			for(i=helixId-1;i>0;i--){
				int firstOpening;
				int lastClosing;
				if(includeDangles){
					firstOpening=(str_int[helices[i]->getI()] == -1) ? helices[i]->getI()+1 : helices[i]->getI();
					lastClosing =(str_int[helices[i]->getJ()] == -1) ? helices[i]->getJ()-1 : helices[i]->getJ();
				}
				else{
					firstOpening= helices[i]->getI();
					lastClosing = helices[i]->getJ();
					
				}

				if(str_int[firstOpening+1]==-1 || str_int[lastClosing-1]==-1 || (lastClosing-firstOpening == 6) 
					 || (helices[i]->getInside().size()==0 && helices[i]->getStackingEnergy(str_int, hairpin37,max_bp_energy)>0)){ //|| helixes[i].getStructure(bp).equals("((...))")){){
					// Remove helix i
					CPHelix* parent = helices[helices[i]->getParent()];
					
					/// EQUAL TO RNAiFold 1 
					parent->addPositions(helices[i]->getPositions());
					parent->addBPs(helices[i]->getBPs());
					parent->addClosingBPs(helices[i]->getClosingBPs());
					parent->addUPs(helices[i]->getUPs());

					
					std::vector<int> inside = helices[i]->getInside();
					parent->removeInside(i);
					for(j=0;j<inside.size();j++){
						helices[inside[j]]->setParent(parent->getId());
						helices[inside[j]]->setLevel(parent->getLevel()+1);
						parent->addInside(inside[j]);
					}
					for(j=0;j<helixId;j++){
						helices[j]->adaptToRemoveOf(i);
						//  helices[j]->toString();
					}
					delete helices[i]; 					
					helices.erase(helices.begin()+i);
					helixId--;
					
				}
			}
		}
		
		return;

	}


	//********************************** END RNAIFOLD HELICES *********************************************//

	int StrTree::getTreeId(){
		return tree_id_;
	}

	void StrTree::setTreeId(int tree_id){
		tree_id_=tree_id;
		return;
	}
	
	double StrTree::getTemp(){
		return temp_;
	}
	
	std::vector<CPHelix* > StrTree::getHelices(){
		return helices;
	}

	std::vector<CPHelix* > StrTree::getSortedHelices(){
		std::vector<CPHelix* > sortedHelices(helices);
		
		sort(sortedHelices.begin(), sortedHelices.end(), compByLevel);
		return sortedHelices;
	}

	void StrTree::initLeavesDistance(){
		for(int i=helices.size()-1;i>=0;i--){
			if(helices[i]->getInside().size() == 0){
				helices[i]->setLeavesDistance(0);
				int tmpDistance=0;
				int tmpHelix = helices[i]->getParent();
				while (tmpHelix != -1 && (helices[tmpHelix]->getLeavesDistance() == -1 || helices[tmpHelix]->getLeavesDistance() > tmpDistance+1)){
					tmpDistance+=1;
					helices[tmpHelix]->setLeavesDistance(tmpDistance);					
					tmpHelix = helices[tmpHelix]->getParent();
				}
			}
		}
	}	

	int StrTree::consecutive(int aX, int aY, int bX, int bY){
		int i,j;
		for(i=1; i<=MAX_LEFT_MISMATCH+1; i++){
			for(j=1; j<=MAX_RIGHT_MISMATCH+1; j++){
				if((aX+i==bX)&&(aY-j==bY))
					return 1;   
			}
		}
		return 0;
	}
	void StrTree::showTree(){
		printf("Struct.\tHelix\tStart-End positions\n");
		for(int i=0; i<helices.size(); i++){
			printf("%d\t%d\t",tree_id_, getHelices()[i]->getId());
			for(int j=0;j<getHelices()[i]->getLevel();j++){
				printf("  ");
			}
			printf("%d-%d %d\n",getHelices()[i]->getI(),getHelices()[i]->getJ(), getHelices()[i]->getLeavesDistance());
			//getHelices()[i]->toString();
		}
	}
}


