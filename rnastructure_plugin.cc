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
#include <stdlib.h>
#include <cstring>
#include <math.h>

#include "rnastructure_plugin.h"


using namespace std;

RNAstructurePlugin::RNAstructurePlugin(int size){
	_winSize =
		( size > 1200 ) ? 20 :
		( size > 800 ) ? 15 :
		( size > 500 ) ? 11 :
		( size > 300 ) ? 7 :
		( size > 120 ) ? 5 :
		( size > 50 ) ? 3 :
		2;
	_temp = 310.15;
	_maxLoop=30;
	_cutpoint = -1;

}
RNAstructurePlugin::RNAstructurePlugin(int size,int dangles){
	_winSize =
		( size > 1200 ) ? 20 :
		( size > 800 ) ? 15 :
		( size > 500 ) ? 11 :
		( size > 300 ) ? 7 :
		( size > 120 ) ? 5 :
		( size > 50 ) ? 3 :
		2;
	_temp = 310.15;
	_maxLoop=30;
	_cutpoint = -1;
	_d = dangles;

}
RNAstructurePlugin::~RNAstructurePlugin(){
	
};

void RNAstructurePlugin::initMe(std::string name, int d) 
{
	_name = name;
	_winSize = 70;
	_maxLoop=30;
	_temp = 310.15;
	_cutpoint = -1;
	_d = d;
}

void RNAstructurePlugin::setSequence(std::string seq)
{
	_seq = seq;
}

//deprecated
void RNAstructurePlugin::setStructure(std::string str)
{
	_str = str;
}

void RNAstructurePlugin::setTestStructure(std::string str)
{
	_tstr = str;
	int i, j,hx;
	int n = str.length();
	vector<int> stack(n);
	_vtstr.assign(n,-1);

	hx=0;

	for (i=0; i<n; i++){
		switch (str[i]){
		case '.':
		     _vtstr[i]= -1;
		     break;
		case '(':
		     stack[hx++]=i;
		     break;
		case ')':
		     j = stack[--hx];
		     if (hx<0){
				//printf("%s \n",str);
				printf("Unbalanced brackets in make_BasePair_Table\n");
				exit(0);
		     }
		     _vtstr[i]=j;
		     _vtstr[j]=i;
		     break;
		}
	}
	if (hx!=0){
//		printf("%s \n",str);
		printf("Unbalanced brackets in make_BasePair_Table\n");
		exit(0);
	}
	return;	
}

void RNAstructurePlugin::setWindowSize(int ws){
	_winSize = ws;
}
void RNAstructurePlugin::setDangles(int d){
  _d = d;
}
void RNAstructurePlugin::setTemperature(double t){
	_temp = 273.15+t;
}

void RNAstructurePlugin::setCutPoint(int p){
	if(p!=-1){
		_cutpoint = p-1;
	}
}


double RNAstructurePlugin::fold(){
	RNA* strand = new RNA(_seq.c_str(), true);

	strand->SetTemperature(_temp);
	int errorCode = strand->FoldSingleStrand( 10, 1, _winSize, "", _maxLoop, true);
	_structure =  strand->GetStructure();
	if(errorCode==0){
		_str.clear();
		for (int j=1;j<=_structure->GetSequenceLength();j++) {
			if (_structure->GetPair(j,1)>j) _str.append("(");
			else if (_structure->GetPair(j,1)==0) _str.append(".");
			else _str.append(")");
		}
		//int e = _structure->GetEnergy(1);
		//return e;
		double e= strand->CalculateFreeEnergy(1);
		delete strand;
		return e;
	}
	delete strand;
	return -1;
}

std::string RNAstructurePlugin::getStructure(){
        return _str;
}

double RNAstructurePlugin::cofold(){
	if(_cutpoint==-1){
		return fold();
	}
	string _seq1=_seq.substr(0,_cutpoint);
	string _seq2=_seq.substr(_cutpoint);
	
	HybridRNA* strand = new HybridRNA(_seq1.c_str(),_seq2.c_str(), true);

	strand->SetTemperature(_temp);
	int errorCode = strand->FoldBimolecular( 10, 0, _winSize, "", false);
	
	if(errorCode==0){
		_structure=strand->GetStructure();
		_str.clear();
		for (int j=1;j<=_structure->GetSequenceLength();j++) {
			if(j<=_seq1.length() || j>_seq1.length()+3){
				if (_structure->GetPair(j,1)>j) _str.append("(");
				else if (_structure->GetPair(j,1)==0) _str.append(".");
				else _str.append(")");
			}
		}
		
		int e = _structure->GetEnergy(1);
		delete strand;
		return e/10.;
	}
	delete strand;
	return -1;
}

double RNAstructurePlugin::energyOfStruct(){
	double e=0;
	
	if(_cutpoint==-1){
		RNA* strand = new RNA(_seq.c_str(), true);
		strand->SetTemperature(_temp);

		strand->GetStructure()->AddStructure();
		for(int i=0; i<_vtstr.size();i++){
			if(_vtstr[i]>i){
				strand->SpecifyPair(i+1,_vtstr[i]+1,1);
			}
		}
		_structure =  strand->GetStructure();
		e=strand->CalculateFreeEnergy(1);
		delete strand;
	}
	else{
		string _seq1=_seq.substr(0,_cutpoint);
		string _seq2=_seq.substr(_cutpoint);
		HybridRNA* strand = new HybridRNA(_seq1.c_str(),_seq2.c_str(), true);
		strand->SetTemperature(_temp);
		strand->SetupBimolecular();
		strand->GetStructure()->AddStructure();
		for(int i=0; i<_vtstr.size();i++){
			if(_vtstr[i]>i){
				strand->SpecifyPair((i>=_seq1.length())? i+4: i+1,(_vtstr[i]>=_seq1.length())? _vtstr[i]+4: _vtstr[i]+1,1);
			}
		}
		e=strand->CalculateFreeEnergy(1);
		delete strand;
	}


	return e;
}

double RNAstructurePlugin::energyOfEnsemble(){
	if(_cutpoint==-1){
		RNA* strand = new RNA(_seq.c_str(), true);
		strand->SetTemperature(_temp);
		strand->PartitionFunction();
		double ensEnergy = strand->GetEnsembleEnergy();
		delete strand;
		return ensEnergy;
	}
	else{
		string _seq1=_seq.substr(0,_cutpoint);
		string _seq2=_seq.substr(_cutpoint);
		HybridRNA* strand = new HybridRNA(_seq1.c_str(),_seq2.c_str(), true);
		strand->SetTemperature(_temp);
		strand->SetupBimolecular();
		strand->PartitionFunction();
		double ensEnergy = strand->GetEnsembleEnergy();
		delete strand;
		return ensEnergy;
	}
}

double RNAstructurePlugin::probOfStruct(){
	double kT = (_temp+273.15)*1.98717/1000.; // in Kcal 

	if(_cutpoint!=-1){
		string _seq1=_seq.substr(0,_cutpoint);
		string _seq2=_seq.substr(_cutpoint);
		
		HybridRNA* strand = new HybridRNA(_seq1.c_str(),_seq2.c_str(), true);

		strand->SetTemperature(_temp);
		int errorCode = strand->FoldBimolecular( 10, 0, _winSize, "", false);
		
		if(errorCode==0){
			int e = strand->GetStructure()->GetEnergy(1);
			strand->PartitionFunction();
			double prob=exp(-(e/10.)/kT)/exp(-(strand->GetEnsembleEnergy())/kT);
			delete strand;
			return prob;
		}		
	}
	else{
		RNA* strand = new RNA(_seq.c_str(), true);
		int length = strand->GetSequenceLength();
		_winSize =
			( length > 1200 ) ? 20 :
			( length > 800 ) ? 15 :
			( length > 500 ) ? 11 :
			( length > 300 ) ? 7 :
			( length > 120 ) ? 5 :
			( length > 50 ) ? 3 :
			2;	
		strand->SetTemperature(_temp);
		int errorCode = strand->FoldSingleStrand( 10, 1, _winSize, "", _maxLoop, true);
		if(errorCode==0){
			double e= strand->CalculateFreeEnergy(1);
			strand->PartitionFunction();
			double prob=exp(-(e)/kT)/exp(-(strand->GetEnsembleEnergy())/kT);
			delete strand;
			return prob;
		}
	}
	return 0;
}

int* RNAstructurePlugin::getBasePairs(){
	int i;
	int n = _seq.size();
	int stck[n];
	int pointer = 0;
	int *bprs;
	bprs = (int *) malloc(sizeof(int)*(n));  
	for(i=0;i<n;i++){
		if(_str[i]=='('){
			stck[pointer]=i;
			pointer++;
		}
		if(_str[i]==')'){
			bprs[i]=stck[pointer-1];
			bprs[stck[pointer-1]]=i;
			pointer--;
		}
		if(_str[i]=='.')
			bprs[i]=-1;
	}
	return bprs;
}

int* RNAstructurePlugin::getBasePairs1Index(){
	int i;
	int n = _seq.size();
	int stck[n];
	int pointer = 0;
	int *bprs;
	bprs = (int *) malloc(sizeof(int)*(n+1));  
	bprs[0]=n;
	for(i=0;i<n;i++){
		if(_str[i]=='('){
			stck[pointer]=i;
			pointer++;
		}
		if(_str[i]==')'){
			bprs[i+1]=stck[pointer-1]+1;
			bprs[stck[pointer-1]+1]=i+1;
			pointer--;
		}
		if(_str[i]=='.')
			bprs[i+1]=-1;
	}
	return bprs;
}

void RNAstructurePlugin::fillBasePairs1Index(int *bprs){
	int i;
	int n = _seq.size();
	int stck[n];
	int pointer = 0;
	bprs[0]=n;
	for(i=0;i<n;i++){
		if(_str[i]=='('){
			stck[pointer]=i;
			pointer++;
		}
		if(_str[i]==')'){
			bprs[i+1]=stck[pointer-1]+1;
			bprs[stck[pointer-1]+1]=i+1;
			pointer--;
		}
		if(_str[i]=='.')
			bprs[i+1]=-1;
	}
	return;
}

double* RNAstructurePlugin::basePairProbsArray(){
	int i,j;
	int n = _seq.size();
	double* da = (double *) malloc(sizeof(double)*(n*n));  
	float sums = 0;
	//update_pf_params(n);
	RNA* strand = new RNA(_seq.c_str(), true);
	strand->SetTemperature(_temp);
	strand->PartitionFunction();
	for(i=0;i<n;i++){
		for(j=i+1;j<n;j++){
			da[i*n+j]=strand->GetPairProbability(i+1,j+1);
			da[j*n+i]=da[i*n+j];
			sums+=da[i*n+j];
		}
		sums = 0;
		for(j=0;j<n;j++)
			if(j!=i)
				sums+=da[i*n+j];
			da[i*n+i] =(1-sums);
	}
	delete strand;
	return da;
}

double**  RNAstructurePlugin::basePairProbsMatrix(){
	int i,j;
	int n = _seq.size();
	double** da = (double **) malloc(sizeof(double*)*(n*n));  
	for(i=0;i<n;i++){
		da[i] = (double *) malloc(sizeof(double)*(n+1));  
	}

	float sums = 0;
	RNA* strand = new RNA(_seq.c_str(), true);
	strand->SetTemperature(_temp);
	strand->PartitionFunction();
	for(i=0;i<n;i++){
		for(j=i+1;j<n;j++){
			da[i][j]=strand->GetPairProbability(i+1,j+1);
			da[j][i]=strand->GetPairProbability(i+1,j+1);
			sums+=da[i][j];
		}
		sums = 0;
		da[i][i] =0;
		for(j=0;j<n;j++)
			if(j!=i)
				sums+=da[i][j];
			da[i][n] =(1-sums);
	}
	delete strand;
	return da;
}


double* RNAstructurePlugin::basePairProbsArrayCofold(){
	int i,j;
	int n = _seq.size();
	double* da = (double *) malloc(sizeof(double)*(n*n));  
	float sums = 0;
	string _seq1=_seq.substr(0,_cutpoint);
	string _seq2=_seq.substr(_cutpoint);
	
	HybridRNA* strand = new HybridRNA(_seq1.c_str(),_seq2.c_str(), true);
	strand->SetTemperature(_temp);
	strand->PartitionFunctionBimolecular("");
	for(i=0;i<n;i++){
		int iindx = (i>=_seq1.length())? i+3: i;
		for(j=i+1;j<n;j++){
			int jindx = (j>=_seq1.length())? j+3: j;
			da[i*n+j]=strand->GetPairProbability(iindx+1,jindx+1);
			da[j*n+i]=da[i*n+j];
			sums+=da[i*n+j];
		}
		sums = 0;
		for(j=0;j<n;j++)
			if(j!=i)
				sums+=da[i*n+j];
			da[i*n+i] =(1-sums);
	}
	delete strand;	
	return da;
}

double**  RNAstructurePlugin::basePairProbsMatrixCofold(){
	int i,j;
	int n = _seq.size();
	double** da = (double **) malloc(sizeof(double*)*(n*n));  
	for(i=0;i<n;i++){
		da[i] = (double *) malloc(sizeof(double)*(n+1));  
	}

	float sums = 0;
	string _seq1=_seq.substr(0,_cutpoint);
	string _seq2=_seq.substr(_cutpoint);
	
	HybridRNA* strand = new HybridRNA(_seq1.c_str(),_seq2.c_str(), true);
	strand->SetTemperature(_temp);
	strand->PartitionFunctionBimolecular("");
	for(i=0;i<n;i++){
		int iindx = (i>=_seq1.length())? i+3: i;
		for(j=i+1;j<n;j++){
			int jindx = (j>=_seq1.length())? j+3: j;
			da[i][j]=strand->GetPairProbability(iindx+1,jindx+1);
			da[j][i]=strand->GetPairProbability(iindx+1,jindx+1);
			sums+=da[i][j];
		}
		sums = 0;
		da[i][i] =0;
		for(j=0;j<n;j++)
			if(j!=i)
				sums+=da[i][j];
			da[i][n] =(1-sums);
	}
	delete strand;

	return da;
}


void RNAstructurePlugin::setEnergyModel(std::string energy_model){
	setenv("DATAPATH", energy_model.c_str(), 1);
	return;
}

void RNAstructurePlugin::freeMe() {
	return;
}





