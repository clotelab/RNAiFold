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

#include "misc.h"
#include "aaconstraint.h"
#include <algorithm>
#include <iostream>


using namespace std;


AAConstraint::~AAConstraint(){};
AAConstraint::AAConstraint(int startPos, std::string  aaTarget, std::string  aaSeq, int maxBlosumScore, int aaSimilCstr, int strLen) : 
	startPos_(startPos),
	aaTarget_(aaTarget),
	aaSeq_(aaSeq),
	maxBlosumScore_(maxBlosumScore),
	aaSimilCstr_(aaSimilCstr),
	strLen_(strLen),
	errCode_(AA_ERR_OK){

	if(aaSimilCstr_<AA_MIN_BLOSUM || aaSimilCstr_ > AA_SIMIL_CLASS2){
		errCode_=AA_ERR_SIMILARITY;
		return;

	}

	if(maxBlosumScore<0 || maxBlosumScore > 1){
		errCode_=AA_ERR_MAXIMIZE;
		return;

	}

	if(startPos_<1 ||  startPos_>strLen_-2){
		errCode_=AA_ERR_STARTPOS;
		return;
	}

	
	// If there is no amino acid sequence constraints but there is an amino acid sequence target, the target sequence is used as the sequence constraints
	if(aaTarget_.length()>0 && aaSeq_.length()==0){
		if(aaSimilCstr==AA_SIMIL_CLASS1 || aaSimilCstr==AA_SIMIL_CLASS2){ // Current classifications
			aaSeq_= getAAConstraintsByClassification(aaTarget_,aaSimilCstr_);
		}
		else if(aaSimilCstr_==AA_SIMIL_NONE && maxBlosumScore == 1){
			aaSeq_.assign(aaTarget_.length(),'X');
		}
		else{
			aaSeq_.assign(aaTarget_);
		}
	}
	// The target sequence is only for blosum score maximization. 
	if(maxBlosumScore==0){
		aaTarget_.clear();
	}

	// Validate target sequence
	if(aaTarget_.length()>0){
		if((aaTarget_.length()*3)+(startPos_-1) > strLen_){
			errCode_=AA_ERR_LENGTH_OVERFLOW;
			return;
		}
		else{
			for(int j=0; j<aaTarget_.length();j++){
				if(std::find(validAaIUPAC.begin(), validAaIUPAC.end(), aaTarget_[j]) == validAaIUPAC.end()) {
					errCode_=AA_ERR_WRONG_TARGET;
					return;
				}
			}
		}
	}
	if(aaSeq_.length()>0){
		for(int j=0; j<aaSeq_.length();j++){
			if(std::find(validAaExtIUPAC.begin(), validAaExtIUPAC.end(), aaSeq_[j]) == validAaExtIUPAC.end()) {
				errCode_=AA_ERR_WRONG_SEQUENCE;
				return;
			}
		}

		// Parse sequence constraints and store as a list of amino acids 
		int bracketStack=0;
		std::string currentString;
		bool isGroup=false;
		for(int i=0;i<aaSeq_.length();i++){
			if(aaSeq_.at(i)=='['){
				bracketStack++;
				isGroup=true;
				currentString.clear();
			}
			else if(aaSeq_.at(i)==']'){
				bracketStack--;
				allowedAAs.push_back(currentString);
				isGroup=false;
			}
			else if(isGroup){
				currentString += aaSeq_.at(i);				
			}
			else {
				allowedAAs.push_back(std::string(1,aaSeq_.at(i)));				
			}
			if(bracketStack!=0 && bracketStack!=1){
				errCode_=AA_ERR_WRONG_BRACKETS;
				return;				
			}
		}
		if(bracketStack!=0){
			errCode_=AA_ERR_WRONG_BRACKETS;
			return;				
		}
		if(aaTarget_.length()>0 && allowedAAs.size() != aaTarget_.length()){
			errCode_=AA_ERR_DIFFERENT_LENGTH;
			return;
		}
		
		if((allowedAAs.size()*3)+(startPos_-1) > strLen_){
			errCode_=AA_ERR_LENGTH_OVERFLOW;
			return;
		}
	}
	size_=allowedAAs.size();

	if(aaSimilCstr_<AA_SIMIL_NONE || maxBlosumScore_==1){
		blosumHelper = new Blosum();
	}

	// Create codon domains
	if(aaSimilCstr_>=AA_SIMIL_NONE){		
		for(int i=0;i<allowedAAs.size();i++){
			std::vector<int> newCodonDomain;			
			for(int j=0;j<allowedAAs.at(i).length();j++){
				std::vector<int> toInsert = aaIUPAC.at(allowedAAs.at(i).at(j));
				for(int k =0; k<toInsert.size();k++){
					if(newCodonDomain.empty() || std::find(newCodonDomain.begin(), newCodonDomain.end(), toInsert[k]) == newCodonDomain.end()){
						newCodonDomain.push_back(toInsert[k]);
					}
				}				
			}
			if(newCodonDomain.empty()){
				//No valid values, define an unfeasible domain
				newCodonDomain.push_back(-1);
			}
			codonDomains.push_back(newCodonDomain);
		}
	}
	else{
		for(int i=0;i<allowedAAs.size();i++){
			std::vector<int> newCodonDomain;			
			for(int j=0;j<allowedAAs.at(i).length();j++){
				std::vector<char> compatibleAA = blosumHelper->getAtDistance(allowedAAs.at(i).at(j), aaSimilCstr_);
				for(int l =0; l<compatibleAA.size();l++){
					std::vector<int> toInsert = aaIUPAC.at(compatibleAA[l]);
					for(int k =0; k<toInsert.size();k++){
						if(newCodonDomain.empty() || std::find(newCodonDomain.begin(), newCodonDomain.end(), toInsert[k]) == newCodonDomain.end()){
							newCodonDomain.push_back(toInsert[k]);
						}
					}
				}
			}
			if(newCodonDomain.empty()){
				//No valid values, define an unfeasible domain
				newCodonDomain.push_back(-1);
			}
			codonDomains.push_back(newCodonDomain);

			
		}
	}
}

std::string AAConstraint::getTarget(){
	return aaTarget_;
}
char AAConstraint::getTrgConstAt(int pos){
	return aaTarget_.at(pos);
}

std::string AAConstraint::getAllowedAaAt(int pos){
	return allowedAAs.at(pos);
}

std::vector<int> AAConstraint::getDomain(int pos){
	return codonDomains.at(pos);
}
int AAConstraint::getSimilarity(){
	return aaSimilCstr_;
}
bool AAConstraint::maximizeScore(){
	return (maxBlosumScore_==1);
}

int AAConstraint::getStartPos(){
	return startPos_;
}

int AAConstraint::getLength(){
	return size_;
}

int AAConstraint::isValid(){
	return errCode_;
}
			
void AAConstraint::toString(){
	printf("Size %d, Target: %s\n",size_, aaTarget_.c_str());
	printf("Amino acid constraints: \n");
	for(int i=0;i<allowedAAs.size();i++){
		printf("%s: ",allowedAAs[i].c_str());
		for(int j=0;j<codonDomains[i].size();j++){
			printf(" %d",codonDomains.at(i).at(j));
		}
		printf("\n");
	}
	return;
}

int AAConstraint::getMaxBlosumValue(int pos){
	return blosumHelper->getBlosumValue(aaTarget_[pos],aaTarget_[pos]);
}

int AAConstraint::getBlosumIndex(int pos){
	return blosumHelper->getIndex(aaTarget_[pos]);
}


std::string AAConstraint::getErrorMessage(){
	switch(errCode_){
		case AA_ERR_LENGTH_OVERFLOW:
			return "Amino acid constraint is longer than 3*(target structure length)!\n";
		case AA_ERR_DIFFERENT_LENGTH:
			return "Both amino acid constraint and target amino acid sequences provided but lengths differ!\n";
		case AA_ERR_WRONG_TARGET:
			return "Invalid values in amino acid target constraint!\n";
		case AA_ERR_WRONG_SEQUENCE:
			return "Invalid values in amino acid sequence constraint!\n";
		case AA_ERR_WRONG_BRACKETS:
			return "Brackets are not balanced in amino acid sequence constraint!\n";
		case AA_ERR_STARTPOS:
			return "Starting position for amino acid constraint must be between 1 and target structure length-3!\n";
		case AA_ERR_SIMILARITY:
			return "Blosum similarity threshold must range between -4 and 7!\n";
		case AA_ERR_MAXIMIZE:
			return "Valid values for MaxBlosumScore are 0 (disabled) and 1 (enabled))!\n";

		default: return "";
	}
}

// Auxiliary functions 

std::string getAAConstraintsByClassification(std::string aaTarget, int aaSimilCstr){
		std::string returnString =""; 
		for(int i =0; i< aaTarget.length();i++){
				if ( aaClassification.at(aaSimilCstr).find(aaTarget.at(i)) == aaClassification.at(aaSimilCstr).end() ) {
					returnString+=aaTarget.at(i);
				}
				else{
					returnString+= aaClassification.at(aaSimilCstr).at(aaTarget.at(i));
				}
		}
		return returnString;
}
void parseAAconstraints(std::vector<AAConstraint*> *aaConstraints, std::string aaTarget, std::string aaConst,std::string aaStartPos, int maxBlosumScore, int aaSimilCstr, int strLen){
	int i;
	vector<string> aaConstList;
	vector<string> aaTargetList;
	
	vector<int> aaStartPosList;

	if(aaTarget.length()>0){
		split(aaTarget,FIELD_DELIMITER, aaTargetList);		
	}

	if(aaConst.length()>0){
		split(aaConst,FIELD_DELIMITER, aaConstList);
	}
	
	int maxAArules = max(aaTargetList.size(),aaConstList.size());	
	if(aaStartPos.length()>0){
		vector<string> aaStartPosListStr;
		split(aaStartPos,FIELD_DELIMITER, aaStartPosListStr);
		if(aaStartPosListStr.size()>maxAArules){
			cout << "There are more starting positions for amino acid constraints specified than amino acid contraints!"  << endl;
			exit(1);
		}
		for(i=0; i<aaStartPosListStr.size();i++){
			if(aaStartPosListStr[i] ==""){
				aaStartPosList.push_back(1);
			}
			else if(!is_number(aaStartPosListStr[i])){
				cout << "Amino acid constraint starting positions must be numeric!" << endl;
				exit(1);
			}
			else{
				aaStartPosList.push_back(atoi(aaStartPosListStr[i].c_str()));
			}
		}
	}
	
	// Set all amino acid constraint vectors to the same size

	if(aaTargetList.size()<maxAArules){
			for(i=aaTargetList.size();i<maxAArules;i++){
				aaTargetList.push_back("");
			}
	}	
	if(aaConstList.size()<maxAArules){
			for(i=aaConstList.size();i<maxAArules;i++){
				aaConstList.push_back("");
			}
	}	
	if(aaStartPosList.size()<maxAArules){
			for(i=aaStartPosList.size();i<maxAArules;i++){
				aaStartPosList.push_back(1);
			}
	}	

	for(int i=0; i<maxAArules;i++){
		aaConstraints->push_back(new AAConstraint(aaStartPosList[i],aaTargetList[i], aaConstList[i], maxBlosumScore,aaSimilCstr, strLen));
//		if(aaConstraints->at(aaConstraints->size()-1)->isValid()==AA_ERR_OK){
//			aaConstraints->at(i)->toString();
//		}
	}	
}
