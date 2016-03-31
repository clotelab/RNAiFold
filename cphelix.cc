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
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "cphelix.h"

namespace operations_research {
	CPHelix::~CPHelix(){};
	CPHelix::CPHelix(int _id, int _i, int _j) : 
		id(_id),
		i(_i),	    
		j(_j){
		parent=-1;
		level=0;
		leaves_distance=-1;

		pos_left=0;
		pos_right=0;
		try_left=0;
		try_right=0;
	}
	
	
	int CPHelix::getId(){
		return id;
	}
	void CPHelix::setId(int _id){
		id=_id;
		return;
	}
	int CPHelix::getI(){
		return i;
	}
	void CPHelix::setI(int _i){
		i=_i;
		return;
	}
	int CPHelix::getJ(){
		return j;
	}
	void CPHelix::setJ(int _j){
		j=_j;
		return;
	}
	int CPHelix::getParent(){
		return parent;
	}
	void CPHelix::setParent(int _parent){
		parent=_parent;
		return;
	}
	int CPHelix::getLevel(){
		return level;
	}
	void CPHelix::setLevel(int _level){
		level=_level;
		return;
	}

	int CPHelix::getLeavesDistance(){
		return leaves_distance;
	}
	void CPHelix::setLeavesDistance(int distance){
		leaves_distance=distance;
		return;
	}

	int CPHelix::getPosLeft(){
		return pos_left;
	}
	void CPHelix::setPosLeft(int npos){
		pos_left=npos;
		return;
	}
	int CPHelix::getPosRight(){
		return pos_right;
	}
	void CPHelix::setPosRight(int npos){
		pos_right=npos;
		return;
	}

	int CPHelix::getTryLeft(){
		return try_left;
	}
	void CPHelix::setTryLeft(int new_val){
		try_left=new_val;
		return;
	}
	int CPHelix::getTryRight(){
		return try_right;
	}
	
	void CPHelix::setTryRight(int new_val){
		try_right=new_val;
		return;
	}

	void CPHelix::addInside(int h){
		inside.push_back(h);
	}
	
	void CPHelix::removeInside(int h){
		for(int c=0;c<inside.size();c++){
			if(inside[c]==h){
				inside.erase(inside.begin()+c);
				break;
			}
		}
	}

	std::vector<int> CPHelix::getInside(){
		return inside;
	}
	int CPHelix::contains(CPHelix* h){
		if(i<=h->getI() && j>=h->getJ()){
			return 1;
		}
		else{
			return 0;
		}
	}

	void CPHelix::setHelixPositions(std::vector<CPHelix*> helices){
		int c,son;

		for(c=i; c<=j; c++){
			bool isInside=true;
			for(son=0; son<inside.size(); son++){
				if(helices[inside[son]]->getI()<=c && helices[inside[son]]->getJ()>=c){
					isInside=false;
				}
			}
			if(isInside==true){
				positions.push_back(c);
			}
		}
	}	
	
	std::vector<int> CPHelix::getPositions(){
		return positions;
	}

	void CPHelix::addPositions(std::vector<int> newPositions){
		positions.insert(positions.end(), newPositions.begin(), newPositions.end());
	}

	int CPHelix::getPositionsIntersect(std::vector<int> newPositions){
		std::vector<int> intersect;
		
		set_intersection(positions.begin(),positions.end(),newPositions.begin(),newPositions.end(),back_inserter(intersect)); 
		
		return intersect.size();
	}
	
	void CPHelix::adaptToRemoveOf(int h){
		if(id>h){
			id--;
		}
		if(parent>h){			
			parent--;
		}
		for(int c=0;c<inside.size();c++){
			if(inside[c]>h){			
				inside[c]--;
			}
		}
		
	}
	
	int CPHelix::getStackingEnergy(int* str_int, int* hairpin37, int bp_energy){
		int nBPs=0;
		int lastOpen=i;
		int firstClose=j;
		for(int c=i; c<=j;c++){
			if(str_int[c]>c){
				nBPs++;
				lastOpen=c;
				firstClose=str_int[c];
			}
		} 
		int loopSize = firstClose-lastOpen-1;
//		if(((nBPs-1)*(bp_energy) + ((loopSize==4) ? 0 : hairpin37[std::min(firstClose-lastOpen-1,30)])) >0){
//			printf("Energy not passed %d %d - %d %d - %d %d\n",lastOpen,firstClose, nBPs-1,firstClose-lastOpen,(nBPs-1)*(bp_energy), hairpin37[std::min(firstClose-lastOpen-1,30)]);
//		}
		return ((nBPs-1)*(bp_energy) + ((loopSize==4) ? 0 : hairpin37[std::min(firstClose-lastOpen-1,30)]));
	}
	
	void CPHelix::setPositionIndexAndType(int* str_int, std::vector<int> BPO, std::vector<int> BPC){
		int i,j;
		bpPositions.clear();
		for(i=0; i< positions.size();i++){
			if(str_int[positions[i]]==-1 || positions[i]<this->i+this->pos_left || positions[i]>this->j-this->pos_right ){				
				ups.push_back(positions[i]);
			}
			else if(str_int[positions[i]]>positions[i]){
				for(j=0; j< BPO.size();j++){
					if(BPO[j]==positions[i]){
						if(BPO[j]==this->getI()+this->getPosLeft()){
							cbps.push_back(j);
							bpPositions.push_back(BPO[j]);
							bpPositions.push_back(BPC[j]);
						}
						else{
							bps.push_back(j);
							bpPositions.push_back(BPO[j]);
							bpPositions.push_back(BPC[j]);

							if(str_int[positions[i]+1]==str_int[positions[i]]-1){
								stackBPs.push_back(j);
							}
						}
					}
				}				
			}			
		}
    
		//STORE UP POSITIONS IN THE OPTIMAL ORDER FOR SEARCH (same strategy as RNAiFold 1)
		std::vector<int> orderedUPs= getOrderedUPs(str_int);
		ups.swap(orderedUPs);
		
	    sort(bpPositions.begin(), bpPositions.end());
	    sort(positions.begin(), positions.end());

		
	}
	
	std::vector<int> CPHelix::getBPs(){
		return bps;
	}

	void CPHelix::addBPs(std::vector<int> newBPs){
		bps.insert(bps.end(), newBPs.begin(), newBPs.end());
	}
	

	std::vector<int> CPHelix::getClosingBPs(){
		return cbps;
	}

	void CPHelix::addClosingBPs(std::vector<int> newCBPs){
		cbps.insert(cbps.end(), newCBPs.begin(), newCBPs.end());
		return;
	}
	
	std::vector<int> CPHelix::getBPpositions(){
		return bpPositions;
	}

	int CPHelix::getBPintersect(std::vector<int> newBPpositions){
		std::vector<int> intersect;
		
		set_intersection(bpPositions.begin(),bpPositions.end(),newBPpositions.begin(),newBPpositions.end(),back_inserter(intersect)); 
		
		return intersect.size();
	}
	
	std::vector<int> CPHelix::getUPs(){
		return ups;
	}

	void CPHelix::addUPs(std::vector<int> newUPs){
		ups.insert(ups.end(), newUPs.begin(), newUPs.end());
	}

	std::vector<int> CPHelix::getOrderedUPs(int* str_int){
		int levelUP=0;
		int maxLevel=0;	
		int maxStretch=0;	

		int stretch=0;	
		int nUPs=0;
		std::vector<int> levelUPs(ups.size());
		std::vector<int> stretchUPs(ups.size());		
		std::vector<int> currentStretch;		
		std::vector<int> orderedUPs;
		std::vector<int> cleanVector;

		for(int k=i;k<=j;k++){
			if(k<i+pos_left || k>j-pos_right || str_int[k]==-1){
				stretch++;
				if(stretch>maxStretch){
					maxStretch=stretch;
				}

				if(std::find(positions.begin(), positions.end(), k) != positions.end()) {
					levelUPs[nUPs]=levelUP;
					currentStretch.push_back(nUPs);
					stretchUPs[nUPs]=1;

					nUPs++;	
				}

			}
			else if(str_int[k]>k){
				levelUP++;
				if(currentStretch.size()!=0){
					for(int m=0; m<currentStretch.size();m++){
						stretchUPs[currentStretch[m]]=stretch;
					}
					currentStretch.clear();
				}
				stretch=0;
				if(levelUP>maxLevel){
					maxLevel=levelUP;
				}
			}
			else {
				levelUP--;
				if(currentStretch.size()!=0){
					for(int m=0; m<currentStretch.size();m++){
						stretchUPs[currentStretch[m]]=stretch;
					}
					currentStretch.clear();
				}
				stretch=0;	
			}
			if(k==j){
					for(int m=0; m<currentStretch.size();m++){
						stretchUPs[currentStretch[m]]=stretch;
					}
			}
		}
		
		if(nUPs != ups.size()){
			printf ("Error assigning levels %d %lu\n",nUPs,ups.size());
			exit(0);
		}

		nUPs=0;
		for(int l=maxLevel;l>=0;l--){
			for(int m=maxStretch;m>=0;m--){	
				if(l==0 && 	stretchUPs.size()>0 && m == stretchUPs[ups.size()-1]){	
					for(int k=ups.size()-1;k>=0;k--){
						if(levelUPs[k]==l && stretchUPs[k]==m){
							orderedUPs.push_back(ups[k]);
							nUPs++;
						}
					}
				}
				else{
					for(int k=0;k<ups.size();k++){
						if(levelUPs[k]==l && stretchUPs[k]==m){
							orderedUPs.push_back(ups[k]);
							nUPs++;
						}
					}
				}
			}
		}		

		if(nUPs != ups.size()){
			printf ("Error assigning levels\n");
			exit(0);
		}

		 
		return orderedUPs;
	}


	bool CPHelix::isStack(int bp){
		if(std::find(stackBPs.begin(), stackBPs.end(), bp) != stackBPs.end()) {
			return true;
		}
		else{
			return false;
		}
	}

	void CPHelix::toString(){
		int c;
		printf("Helix:%d --> %d-%d  Level:%d Parent: %d Contains:",id,i,j,level, parent); 
		for(c=0;c<inside.size();c++){
			printf(" %d",inside[c]);
		}

		printf(" - Positions:"); 
		for(c=0;c<positions.size();c++){
			printf(" %d",positions[c]);
		}

		printf(" - Bps:"); 
		for(c=0;c<bps.size();c++){
			if(std::find(stackBPs.begin(), stackBPs.end(), bps[c]) != stackBPs.end()) {
				printf(" %d*",bps[c]);				
			}
			else{
				printf(" %d",bps[c]);
			}
		}

		printf(" - Ups:"); 
		for(c=0;c<ups.size();c++){
			printf(" %d",ups[c]);
		}

		printf("\n"); 
	}


}
