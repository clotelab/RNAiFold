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
#include "local_constraint.h"
#include "misc.h"
#include <iostream>


using namespace std; 
LocalCstr::~LocalCstr(){};

LocalCstr::LocalCstr(std::string structure, int limitType, double limitValue) : LocalCstr(1, structure, limitType, limitValue){}

LocalCstr::LocalCstr(int start_pos, std::string structure, int limitType, double limitValue) : 
	start_pos_(start_pos),
	structure_(structure),
	limitType_(limitType),
	limitValue_(limitValue), 
	has_undet_(false){
		if(validSecondaryStructure(structure_, true)<0){
			limitType_=LC_WRONG_STRUCTURE;
		}
		else{
			int* tmp_int_str_=make_BasePair_Table(structure_, true);
			for(int i =0; i<=structure_.length();i++){
				if(tmp_int_str_[i]==-2) has_undet_= true;
				int_str_.push_back(tmp_int_str_[i]);
			}
			free(tmp_int_str_);
		}
	}

LocalCstr::LocalCstr(std::string structure, const std::string limitType, double limitValue) : LocalCstr(1, structure, limitType, limitValue){}

LocalCstr::LocalCstr(int start_pos, std::string structure, const std::string limitType, double limitValue) : 
	start_pos_(start_pos),
	structure_(structure),
	limitValue_(limitValue), 
	has_undet_(false){
	std::map<const std::string,int>::iterator it = limitTypes_.find(limitType);
	if(it != limitTypes_.end()){
		limitType_ = it->second;
	}
	else{
		cout << "Invalid helix constraint type ("<<limitType<< ")"<<endl;
		limitType_=LC_WRONG_TYPE;
	}
	if(validSecondaryStructure(structure_, true)<0){
		limitType_=LC_WRONG_STRUCTURE;
	}
	else{
		findCutPoint(&structure_);
		int* tmp_int_str_=make_BasePair_Table(structure_, true);
		for(int i =0; i<=structure_.length();i++){
			if(tmp_int_str_[i]==-2) has_undet_= true;
			int_str_.push_back(tmp_int_str_[i]);
		}
		free(tmp_int_str_);
	}
}	

int LocalCstr::getStartPos(){
	return start_pos_;
}
void LocalCstr::setStartPos(int id){
	start_pos_=id;
	return;
}

std::string LocalCstr::getStructureStr(){
	return structure_;
}
void LocalCstr::setStructureStr(std::string structure){
	structure=structure;
	return;
}

std::vector<int64> LocalCstr::getStructureArr(){
	return int_str_;
}
void LocalCstr::setStructureArr(std::vector<int64> int_str){
	int_str_=int_str;
	return;
}

int LocalCstr::getLimitType(){
	return limitType_;
}
void LocalCstr::setLimitType(int limitType){
	limitType_=limitType;
	return;
}

double LocalCstr::getLimitValue(){
	return limitValue_;
}
void LocalCstr::setLimitValue(double limitValue){
	limitType_=limitValue;
	return;
}

bool LocalCstr::hasUndet(){
		return has_undet_;
}

bool LocalCstr::isPairedConstraint(LocalCstr c){
	if (this->getStructureStr().compare(c.getStructureStr())==0 && this->getStartPos()== c.getStartPos()){
	   return true;
	}
	return false;
}

std::string LocalCstr::printTypes(){
	std::string outString="";
	for (std::map<const std::string,int>::iterator it = limitTypes_.begin(); it != limitTypes_.end(); ++it){		
		outString+=it->first;
		outString+=" ";		
	}
	return outString;
}
