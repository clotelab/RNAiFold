
#if !defined (RNA_LIBRARY)
#define RNA_LIBRARY


#include "structure.h"

#include <cmath>
#include <fstream>
#include "stackclass.h"
#include "stackstruct.h"



struct datatable //this structure contains all the info read from thermodynamic
							//data files
{
 short int poppen [5],maxpen,eparam[11],dangle[6][6][6][3],inter[31],bulge[31],
		hairpin[31],stack[6][6][6][6],tstkh[6][6][6][6],tstki[6][6][6][6],
		tloop[maxtloop+1][2],numoftloops,iloop22[6][6][6][6][6][6][6][6],
		iloop21[6][6][6][6][6][6][6],iloop11[6][6][6][6][6][6],
      coax[6][6][6][6],tstackcoax[6][6][6][6],coaxstack[6][6][6][6],
      tstack[6][6][6][6],tstkm[6][6][6][6],auend,gubonus,cint,cslope,c3,
      efn2a,efn2b,efn2c,triloop[maxtloop+1][2],numoftriloops,init,mlasym,strain,numofhexaloops,
	  singlecbulge,tstki23[6][6][6][6],tstki1n[6][6][6][6];
 int hexaloop[maxtloop+1][2];
 float RT;

float prelog;
datatable();


};






integersize ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data);
//this function calculates flush coaxial stacking
//it requires sequence in i,j,ip, and jp
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data); 
//this funtion calculates an intervening mismatch coaxial stack
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data); 
//this funtion calculates an intervening mismatch coaxial stack
integersize ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data);
int decon1(int x);//used by ergmulti to find a nucleotide from a base pair
int decon2(int x);//used by ergmulti to find a nucleotide from a base pair
integersize ergmulti(int st, int ip, structure *ct, datatable *data, bool simplemb);
//calculate the multi branch loop free energy for a loop starting at nuc ip
//	in structure number st of ct
integersize ergexterior(int st, structure *ct, datatable *data);
//calculate the exterior loop free energy in structure number ip



//write is used to write data to a save file
void write(std::ofstream *out,short *i);
void write(std::ofstream *out,bool *i);
void write(std::ofstream *out,int *i);
void write(std::ofstream *out,char *i);
void write(std::ofstream *out,string *i);
void write(std::ofstream *out,float *i);
void write(std::ofstream *out,double *i);
void writesinglechar(std::ofstream *out,char *i); 

//read is used to read data from a save file
void read(std::ifstream *out,short *i);
void read(std::ifstream *out,bool *i);
void read(std::ifstream *out,int *i);
void read(std::ifstream *out,char *i);
void read(std::ifstream *out,string *i);
void read(std::ifstream *out,float *i);
void read(std::ifstream *out,double *i);
void readsinglechar(std::ifstream *out,char *i);

int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, char *hexaloop,
	  char *tstacki23, char *tstacki1n, datatable* data);
		//gets thermodynamic data from data files

void push(stackstruct *stack,int a,int b,int c,int d);//push info onto the stack
void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz);
														//pull info from the stack

integersize erg1(int i,int j,int ip,int jp,structure *ct,datatable *data);
		//calculates energy of stacked base pairs
integersize erg2(int i,int j,int ip,int jp,structure *ct,datatable *data,char a,
	char b);
		//calculates energy of a bulge/internal loop
integersize erg2in(int i,int j,int ip,int jp,structure *ct, datatable *data,char a,
	char b);
		//calculates the energy of an interior part of an internal loop (includes asymmetry)
integersize erg2ex(int i,int j,int size,structure *ct, datatable *data);
		//calculates the energy of an exterior part of an internal loop (only has length and terminal stack components)
integersize erg3(int i,int j,structure *ct,datatable *data,char dbl);
		//calculates energy of a hairpin loop
integersize erg4(int i,int j,int ip,int jp,structure *ct,datatable *data,
	bool lfce);
		//calculates energy of a dangling base

inline integersize SHAPEend(int i, structure *ct);//calculate the SHAPE pseudo energy for a single
									//paired nucleotide

//this function calculates whether a terminal pair i,j requires the end penalty
inline integersize penalty(int i,int j,structure* ct, datatable *data) {
	integersize energy;
	if (ct->numseq[i]==4||ct->numseq[j]==4)
   	energy= data->auend;
	else energy= 0;//no end penalty
	
	return (energy/*+SHAPEend(i,ct)+SHAPEend(j,ct)*/);



}

integersize penalty2(int i, int j, datatable *data);
integersize ergcoax(int i, int j, int ip, int jp, int k, structure *ct, datatable *data);
	//returns the free energy of coaxial stacking of i-j onto ip-jp


void dG_T (float T, datatable &data, datatable &dhdata, datatable &dg);
	//Determine the free energy parameters at temperature T and store in dg 
		//from the free energy parameters at 37 degrees C (data) and the
		//temperature independent enathalpies (dhdata).

//Given the absolute temperature, T, return the free energy for a paramater as determined from the 
//free energy at 37 degrees C (dG) and the enthalpy change (assumed to be temperature independent).
inline short Tscale(float T,short dG, short dH)
{
	if (dG==INFINITE_ENERGY) return INFINITE_ENERGY; //keep infinity==infinity
	return dH - (short)floor((float)(dH-dG)*T/T37inK+0.5);
}


//When adding the included and excluded fragments for a structure with pair i-j, the SHAPE free energy needs correction
//so that it is not counted twice.


//Write the structure class-specific items in a save file
void writestructuresave(ofstream *out, structure *ct);

//Read the structure class-specific items in a save file
void openstructuresave(ifstream *out, structure *ct);



#endif


