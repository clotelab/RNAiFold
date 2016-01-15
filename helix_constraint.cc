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
#include "helix_constraint.h"
#include <iostream>

using namespace std; 
HelixCstr::~HelixCstr(){};
HelixCstr::HelixCstr(int helix_id, int limitType, double limitValue) : 
	tree_id_(0),
	helix_id_(helix_id),
	limitType_(limitType),
	limitValue_(limitValue){}

HelixCstr::HelixCstr(int tree_id, int helix_id, int limitType, double limitValue) : 
	tree_id_(tree_id),
	helix_id_(helix_id),
	limitType_(limitType),
	limitValue_(limitValue){}

HelixCstr::HelixCstr(int helix_id, const std::string limitType, double limitValue) : 
	tree_id_(0),
	helix_id_(helix_id),
	limitValue_(limitValue){
	std::map<const std::string,int>::iterator it = limitTypes_.find(limitType);
	if(it != limitTypes_.end()){
		limitType_ = it->second;
	}
	else{
		cout << "Invalid helix constraint type ("<<limitType<< ")"<<endl;
		limitType_=HC_WRONG_TYPE;
	}
}

HelixCstr::HelixCstr(int tree_id, int helix_id, const std::string limitType, double limitValue) : 
	tree_id_(tree_id),
	helix_id_(helix_id),
	limitValue_(limitValue){
	std::map<const std::string,int>::iterator it = limitTypes_.find(limitType);
	if(it != limitTypes_.end()){
		limitType_ = it->second;
	}
	else{
		cout << "Invalid helix constraint type ("<<limitType<< ")"<<endl;
		limitType_=HC_WRONG_TYPE;
	}
}	

int HelixCstr::getTreeId(){
	return tree_id_;
}
void HelixCstr::setTreeId(int id){
	tree_id_=id;
	return;
}

int HelixCstr::getHelixId(){
	return helix_id_;
}
void HelixCstr::setHelixId(int id){
	helix_id_=id;
	return;
}

int HelixCstr::getLimitType(){
	return limitType_;
}
void HelixCstr::setLimitType(int limitType){
	limitType_=limitType;
	return;
}

double HelixCstr::getLimitValue(){
	return limitValue_;
}
void HelixCstr::setLimitValue(double limitValue){
	limitType_=limitValue;
	return;
}
std::string HelixCstr::printTypes(){
	std::string outString="";
	for (std::map<const std::string,int>::iterator it = limitTypes_.begin(); it != limitTypes_.end(); ++it){		
		outString+=it->first;
		outString+=" ";		
	}
	return outString;
}
