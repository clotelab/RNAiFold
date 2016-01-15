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
 
 #include "vienna_plugin.h"
 #include "ViennaRNA/fold_vars.h"
 
 extern "C" {
	float fold(char* sequence, char* structure);
	float cofold(char* sequence, char* structure);
	void initialize_fold(int length);
	void initialize_cofold(int length);
	float pf_fold(char* sequence, char* structure);
	void init_pf_fold(int length);
	void init_pf_foldLP(int length);
	float energy_of_structure(char* sequence, char* structure, int verbosity_level);
	void free_arrays();
	void free_co_arrays();
	void free_pf_arrays();
	void free_pf_arraysLP(void);
	void update_pf_params(int length);
	void init_co_pf_fold(int length);
	cofoldF co_pf_fold(char* sequence, char* structure);
	void free_co_pf_arrays();
	void read_parameter_file(const char fname[]);
 }
 
 using namespace std;
 
 ViennaPlugin::ViennaPlugin(int size){
 	_winSize = 70;
 	_temp = 37.0;
 	_cstr = new char[size+1];
 	_ctstr = new char[size+1];
 	_cseq = new char[size+1];
 	_cutpoint = -1;
 
 }
 ViennaPlugin::ViennaPlugin(int size,int dangles){
 	_winSize = 70;
 	_temp = 37.0;
 	_cstr = new char[size+1];
 	_ctstr = new char[size+1];
 	_cseq = new char[size+1];
 	_cutpoint = -1;
 	_d = dangles;

 }
 ViennaPlugin::~ViennaPlugin(){
 	delete[] _cstr;
 	delete[] _ctstr;
 	delete[] _cseq;
 	
 };
 
 void ViennaPlugin::initMe(std::string name, int d) 
 {
 	_name = name;
 	_winSize = 70;
 	_temp = 37.0;
 	_cstr = NULL;
 	_ctstr = NULL;
 	_cseq = NULL;
 	_cutpoint = -1;
 	_d = d;
   //      cout << "initializing vienna(" << this << "," << _name << ")" << endl;
 }
 
 void ViennaPlugin::setSequence(std::string seq)
 {
 	_seq = seq;
 /*	if(_cseq == NULL || (strlen(_cseq) != _seq.size()+1)){
 		if(_cseq != NULL){
 		        free(_cseq);
 		}
 		_cseq = new char[_seq.size()+1];
 	}
 
 	if(_cstr == NULL || (strlen(_cstr) != _seq.size()+1)){
 		if(_cstr != NULL){
 		        free(_cstr);
 		}
 	        _cstr = new char[_seq.size()+1];
 	}*/
 
 	strcpy(_cseq, _seq.c_str());
 }
 
 //deprecated
 void ViennaPlugin::setStructure(std::string str)
 {
 	_str = str;
 
 /*	if(_cstr == NULL || (strlen(_cstr) != _str.size()+1)){
 		if(_cstr != NULL){
 		        free(_cstr);
 		}
 		_cstr = new char[_str.size()+1];
 	}*/
 
 	strcpy(_cstr, _str.c_str());
 }
 
 void ViennaPlugin::setTestStructure(std::string str)
 {
        _tstr = str;
 /*	if(_ctstr == NULL || (strlen(_ctstr) != _tstr.size()+1)){
 		if(_ctstr != NULL){
 		        free(_ctstr);
 		}
 	        _ctstr = new char[_tstr.size()+1];
 	}*/
        strcpy(_ctstr, _tstr.c_str());
 }
 
 double ViennaPlugin::fold(){
        dangles = _d;
 	temperature =_temp;
	initialize_fold(_seq.size());
 	float e = ::fold(_cseq,_cstr);
 	_str = string(_cstr);
        free_arrays();
 	return e;
 }
 
 std::string ViennaPlugin::getStructure(){
        return _str;
 }
 
 double ViennaPlugin::energyOfStruct(){
 	dangles = _d;
 	temperature = _temp;
 	if(_cutpoint!=-1){
 		cut_point = _cutpoint+1;
 		//initialize_cofold(_seq.size());	
 	}
 	else{
 		//initialize_fold(_seq.size());
 	}
 	float e = energy_of_structure(_cseq,_ctstr,0);
 	return e;
 }
 
 double ViennaPlugin::energyOfEnsemble(){
 	dangles = _d;
 	temperature = _temp;
 	float e;
 	if(_cutpoint!=-1){
 		cut_point = _cutpoint;
        init_co_pf_fold(_seq.size());
 		cofoldF co = co_pf_fold(_cseq,NULL);
 		e=co.FAB;
         free_pf_arrays();
 	}
 	else{
        init_pf_fold(_seq.size());
 		e = pf_fold(_cseq,NULL);
         free_pf_arrays();
 	}
 	return e;
 }
 
 double ViennaPlugin::probOfStruct(){
 	dangles = _d;
 	temperature = _temp;
 	int n = _seq.size();
 	double kT = (_temp+273.15)*1.98717/1000.; /* in Kcal */
 
 	if(_cutpoint!=-1){
 		cut_point = _cutpoint;
		initialize_cofold(n);
 		float mfe = ::cofold(_cseq,_cstr);
 		free_co_arrays();
 
		init_co_pf_fold(n);
 		cofoldF e = co_pf_fold(_cseq,NULL);
 		free_co_pf_arrays();
 		float prob = exp(-mfe/kT)/exp(-e.FAB/kT);
 		return prob;
 
 	}
 	else{
 		float mfe = energy_of_structure(_cseq,_ctstr,0);
		init_pf_fold(n);
 		float e = pf_fold(_cseq,NULL);
 		free_pf_arrays();
 		float prob = exp(-mfe/kT)/exp(-e/kT);
 		return prob;
 	}
 }
 
 int* ViennaPlugin::getBasePairs(){
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
 
 int* ViennaPlugin::getBasePairs1Index(){
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
 
 double* ViennaPlugin::basePairProbsArray(){
 	int i,j;
 	int n = _seq.size();
 	double* da = (double *) malloc(sizeof(double)*(n*n));  
 	float sums = 0;
 	dangles = _d;
 	temperature = _temp;
 	//update_pf_params(n);
	init_pf_fold(n);
 	pf_fold(_cseq,_cstr);
 	for(i=0;i<n;i++){
 		for(j=i+1;j<n;j++){
			da[i*n+j]=pr[iindx[i+1] - (j+1)];
			da[j*n+i]=pr[iindx[i+1] - (j+1)];
			sums+=pr[iindx[i+1]-(j+1)];
 		}
 		sums = 0;
 		for(j=0;j<n;j++)
 			if(j!=i)
 				sums+=da[i*n+j];
 			da[i*n+i] =(1-sums);
 	}
 
 	free_pf_arrays();
 	return da;
 }
 
 double**  ViennaPlugin::basePairProbsMatrix(){
 	int i,j;
 	int n = _seq.size();
 	double** da = (double **) malloc(sizeof(double*)*(n*n));  
 	for(i=0;i<n;i++){
 		da[i] = (double *) malloc(sizeof(double)*(n+1));  
 	}
 	float sums = 0;
 	dangles = _d;
 	temperature = _temp;
 	//update_pf_params(n);
	init_pf_fold(n);
 	pf_fold(_cseq,_cstr);
 	for(i=0;i<n;i++){
 		for(j=i+1;j<n;j++){
			da[i][j]=pr[iindx[i+1] - (j+1)];
			da[j][i]=pr[iindx[i+1] - (j+1)];
			sums+=pr[iindx[i+1]-(j+1)];
 		}
 		sums = 0;
 		da[i][i] =0;
 		for(j=0;j<n;j++)
 			if(j!=i)
 				sums+=da[i][j];
 			da[i][n] =(1-sums);
 	}
 
 	free_pf_arrays();
 	return da;
 }
 
 
 double* ViennaPlugin::basePairProbsArrayCofold(){
 	int i,j;
 	int n = _seq.size();
 	double* da = (double *) malloc(sizeof(double)*(n*n));  
 	float sums = 0;
 	dangles = _d;
 	temperature = _temp;
 	//update_pf_params(n);
	init_co_pf_fold(n);
 	co_pf_fold(_cseq,_cstr);
 	for(i=0;i<n;i++){
 		for(j=i+1;j<n;j++){
			da[i*n+j]=pr[iindx[i+1] - (j+1)];
			da[j*n+i]=pr[iindx[i+1] - (j+1)];
			sums+=pr[iindx[i+1]-(j+1)];
 		}
 		sums = 0;
 		for(j=0;j<n;j++)
 			if(j!=i)
 				sums+=da[i*n+j];
 			da[i*n+i] =(1-sums);
 	}
 
 	free_co_pf_arrays();
 	return da;
 }
 
 double**  ViennaPlugin::basePairProbsMatrixCofold(){
 	int i,j;
 	int n = _seq.size();
 	double** da = (double **) malloc(sizeof(double*)*(n*n));  
 	for(i=0;i<n;i++){
 		da[i] = (double *) malloc(sizeof(double)*(n+1));  
 	}
 	float sums = 0;
 	dangles = _d;
 	temperature = _temp;
 	//update_pf_params(n);
	init_co_pf_fold(n);
 	co_pf_fold(_cseq,_cstr);
 	for(i=0;i<n;i++){
 		for(j=i+1;j<n;j++){
			da[i][j]=pr[iindx[i+1] - (j+1)];
			da[j][i]=pr[iindx[i+1] - (j+1)];
			sums+=pr[iindx[i+1]-(j+1)];
 		}
 		sums = 0;
 		da[i][i] =0;
 		for(j=0;j<n;j++)
 			if(j!=i)
 				sums+=da[i][j];
 			da[i][n] =(1-sums);
 	}
 
 	free_co_pf_arrays();
 	return da;
 }
 //CometFloatArray ViennaPlugin::localBasePairProbsArray(){
   /*      int i,j;
         int n = _seq.size();
         double prob;
         double sums;
         int pos;
         bool exit;
         struct plist *dp2 = NULL;
         struct plist **dpp2 = &dp2;
         CometFloatArray dal(n*n);
         double **pup = new double*[n+1];
         pup[0] = new double[1];
         pup[0][0] = 20;
         dangles = _d;
 	temperature = _temp;
         init_pf_foldLP(n);  
         struct plist* lbpprs = pfl_fold(_cseq,_winSize,_winSize,0,pup,dpp2,NULL,NULL);
         exit = false;
         pos = 0;
         while(!exit){
           i = lbpprs[pos].i - 1;
           j = lbpprs[pos].j - 1;
           prob = (double)lbpprs[pos].p;
           if((i==j)&&(i==-1))
             exit = true;
           else{
             dal.set(i*n+j,prob);
             dal.set(j*n+i,prob);
             pos++;
           }
         }
         //q_i
         for(i=0;i<n;i++)
           dal.set(i*n+i,pup[i+1][1]);
 
         delete[] pup;
         free(lbpprs);
         //scale
         for(i=0;i<n;i++){
           sums=0;
           for(j=0;j<n;j++){
             sums+=dal[i*n+j];
           }
           for(j=0;j<n;j++){
             dal.set(i*n+j,(dal[i*n+j]/sums));
           }
         }
         free_pf_arraysLP();
         return dal;*/
 //}
 
 //CometFloatArray ViennaPlugin::mixedBasePairProbsArray(){
 /*        int n = _seq.size();
         CometFloatArray mda(n*n);
         CometFloatArray da = basePairProbsArray();
         CometFloatArray dal = localBasePairProbsArray();
         double pstar[n];
         double sum = 0;
         for(int i=0;i<n;i++){
           sum = 0;
           for(int j=max(0,i-_winSize);j<min(n,i+_winSize);j++)
             sum += da[i*n+j];
           pstar[i] = sum;
         }         
         for(int i=0;i<n;i++)
           for(int j=i; j<n;j++){
             if(j-i>=_winSize){
               mda.set(i*n+j,da[i*n+j]);
               mda.set(j*n+i,mda[i*n+j]);
             }else{
                 mda.set(i*n+j,(pstar[i]*dal[i*n+j]));
                 mda.set(j*n+i,mda[i*n+j]);
             }
           }
 
         return mda;*/
 //}
 
 
 void ViennaPlugin::setWindowSize(int ws){
 	_winSize = ws;
 }
 
 void ViennaPlugin::setTemperature(double t){
 	_temp = t;
 }
 
 void ViennaPlugin::setCutPoint(int p){
  _cutpoint = p;
 }
 
 double ViennaPlugin::cofold(){
        cut_point = _cutpoint;
        dangles = _d;
 	temperature = _temp;
        initialize_cofold(_seq.size());
        float e = ::cofold(_cseq,_cstr);
        _str = string(_cstr);
 	free_co_arrays();
        return e;
 }
 
 void ViennaPlugin::freeMe() {
 	if(_cseq != NULL){
 		free(_cseq);
 	}
 	if(_cstr != NULL){
 	        free(_cstr);
 	}
 	if(_ctstr != NULL){
 		free(_ctstr);
 	}
 }
 
 void ViennaPlugin::setDangles(int d){
  _d = d;
 }
 
 void ViennaPlugin::setEnergyModel(string energy_model){
 	if(energy_model.find(TURNER_04_FILE)==string::npos){
 		read_parameter_file(energy_model.c_str());
 	}
 }
