#include <string>
#include <stdlib.h>
#include <vector>
#if !defined(STRUCTURE_H)
#define STRUCTURE_H


#include <vector>
#include <string>
using namespace std;

#include "defines.h"

#ifdef EXTENDED_DOUBLE
	#include "extended_double.h" //inlcude code for extended double if needed
#endif//defined EXTENDED_DOUBLE

//This is a depracted array size that must be removed:
#define maxforce 3000

//Single structure is a wrapper for the information associated with just a single structure, i.e. pairs, energy, and labels.
struct singlestructure {

	//This function sizes the vectors to the approprate size:
	singlestructure(int sequencelength);

	//keep track of the basepairs
	vector<int> basepr;

	//keep track of the energy of the structure, if available
	int energy;

	//keep a string, from a ct file or sequence file with the sequence desription
	string ctlabel;


};


//! structure Class.
/*!
	The structure class provides a class for handling sequences and structures.
*/



//////////////////////////////////////////////////////////////////////
class structure //this structure contains all the info for a structure
{
/*<<<<<<< structure.h
public:
int numofstructures;//number of structures 
int numofbases;//number of nucleotides in sequence
short int pair[maxforce][2],npair,nforbid,forbid[maxforce][2];//arrays to hold lists of forced pairs or pairs forbidden
short int *numseq,*hnumber;
int **basepr;
int *energy;//[maxstructures+1];
char **ctlabel;//[maxstructures+1][ctheaderlength];
short int ndbl, dbl[maxforce];
int inter[3],allocatedstructures;
short int nnopair,*nopair,nmod,mod[maxforce];
int nopairmax;//maximum number of nucleotides not allowed to pair
short int ngu,gu[maxgu];
char *nucs;
bool intermolecular,allocated,templated,stacking;
bool **tem;//tem stores template information as to whether a pair is allowed
bool SHAPEFileRead,distsread;
structure(int structures = maxstructures+1);
~structure();
void allocate(int size = maxbases);
void allocatestructure();
void deletestructure();
void checknopair();//make sure that the nopair array is large enough to add more items
void checknumberofstructures();//check to make sure there is room for one more structure
void allocatetem();//allocate space in **tem 
short int min_gu, min_g_or_u;//NMR-derived constraint variables
short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
//regional NMR constraints:
short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
//microarray type constraints:
short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
void sort();//sort the structures so that the lowest free energy structure is in position one
			//NOTE: This function sorts the the list of base pairs according to energy.  The energies (*energy) and basepairs (**basepr)
			//		arrays are correctly sorted, but ct labels (**ctlabel) is not sorted.
double CalculatePseudoEnergy(double data, std::string modifier, double, double);
double Gammadist(double data, double shape, double loc, double scale);
double Potential(double data, std::vector< std::vector<double> > params, double kT);
void ReadProbabilisticPotentialParams();//Read chemical modifier distributions from file
void ReadSHAPE(const char *filename, std::string modifier="SHAPE", bool calculate=true, bool nosum=false);//Read SHAPE reactivity data from a file
void ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file
void ReadOffset(const char *SSOffset, const char *DSOffset);//Read Free Energy Offset Files.
void ReadExperimentalPairBonus(const char *filename, double const experimentalOffset = 0.0, double const experimentalScaling = 1.0 );
bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
int maxdistance;//maximum distance between nucs in base pairs
double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
double **EX;// double array that contains experimental bonuses/penalties
bool shaped;//keeps track of whether SHAPE data was loaded
bool experimentalPairBonusExists;//keeps track of whether experimental bonus data was loaded
bool ssoffset;//keeps track of whether a single stranded offset was read from disk
double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability
//SINGLE STRANDED SHAPE ENERGY VARIABLES AND FUNCTIONS
double *SHAPEss; //short int array that contains SHAPE data for single-stranded segments
double SHAPEslope_ss, SHAPEintercept_ss; //values of the slope and intercept for SHAPE data modifying single stranded loop stability
short int **SHAPEss_region;  //2-d short int array containing energy values for hairpin loop combinations
//Parameters for distributions
std::vector< std::vector<double> > SHAPE_params;
std::vector< std::vector<double> > DMS_params;
std::vector< std::vector<double> > CMCT_params;
int SHAPEss_calc(int index_i, int index_j);  //Returns pseudoenergy term for a hairpin loop using single stranded SHAPE data
short int SHAPEss_give_value(int index);  //Returns the single stranded SHAPE pseudo energy for a given nucleotide
||||||| ../../RNAstructure-before-fork/src/structure.h
public:
int numofstructures;//number of structures 
int numofbases;//number of nucleotides in sequence
short int pair[maxforce][2],npair,nforbid,forbid[maxforce][2];//arrays to hold lists of forced pairs or pairs forbidden
short int *numseq,*hnumber;
int **basepr;
int *energy;//[maxstructures+1];
char **ctlabel;//[maxstructures+1][ctheaderlength];
short int ndbl, dbl[maxforce];
int inter[3],allocatedstructures;
short int nnopair,*nopair,nmod,mod[maxforce];
int nopairmax;//maximum number of nucleotides not allowed to pair
short int ngu,gu[maxgu];
char *nucs;
bool intermolecular,allocated,templated,stacking;
bool **tem;//tem stores template information as to whether a pair is allowed
structure(int structures = maxstructures+1);
~structure();
void allocate(int size = maxbases);
void allocatestructure();
void deletestructure();
void checknopair();//make sure that the nopair array is large enough to add more items
void checknumberofstructures();//check to make sure there is room for one more structure
void allocatetem();//allocate space in **tem 
short int min_gu, min_g_or_u;//NMR-derived constraint variables
short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
//regional NMR constraints:
short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
//microarray type constraints:
short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
void sort();//sort the structures so that the lowest free energy structure is in position one
			//NOTE: This function sorts the the list of base pairs according to energy.  The energies (*energy) and basepairs (**basepr)
			//		arrays are correctly sorted, but ct labels (**ctlabel) is not sorted.
void ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file
void ReadSHAPE(const char *filename, bool calculate=true);//Read SHAPE reactivity data from a file
void ReadOffset(const char *SSOffset, const char *DSOffset);//Read Free Energy Offset Files.
void ReadExperimentalPairBonus(const char *filename, double const experimentalOffset = 0.0, double const experimentalScaling = 1.0 );
bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
int maxdistance;//maximum distance between nucs in base pairs
double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
double **EX;// double array that contains experimental bonuses/penalties
bool shaped;//keeps track of whether SHAPE data was loaded
bool experimentalPairBonusExists;//keeps track of whether experimental bonus data was loaded
bool ssoffset;//keeps track of whether a single stranded offset was read from disk
double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability
//SINGLE STRANDED SHAPE ENERGY VARIABLES AND FUNCTIONS
double *SHAPEss; //short int array that contains SHAPE data for single-stranded segments
double SHAPEslope_ss, SHAPEintercept_ss; //values of the slope and intercept for SHAPE data modifying single stranded loop stability
short int **SHAPEss_region;  //2-d short int array containing energy values for hairpin loop combinations
int SHAPEss_calc(int index_i, int index_j);  //Returns pseudoenergy term for a hairpin loop using single stranded SHAPE data
short int SHAPEss_give_value(int index);  //Returns the single stranded SHAPE pseudo energy for a given nucleotide
=======*/
	public:



		//******************************
		//Constructor:
		//******************************

		//!Constructor.
		//!	\param sructures is an int that specifies how many structures should be anticipated.  This sets up an initial memory allocation, but this can expand as needed.
		structure(int structures = maxstructures+1);

		//!Destructor.
		~structure();

		//*********************************
		//Get and receive sequence and structure information:
		//*********************************

		//! Get the label for structure numer structurenumber.

		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \return a string that gives the label.
		string GetCtLabel(int structurenumber);

		//! Get the energy for structure numer structurenumber.

		//! This function requires that an energy calculation has been performed.  It does not do an eenrgy calculation.
		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \return an int that gives the energy in kcal/mol*conversionfactor, where conversionfactor is set in /src/defines.h.
		int GetEnergy(int structurenumber);
		
		//! Get the number of structures stored.

		//! \return An integer that is the number of structures encoded.
		int GetNumberofStructures();

		//! Get the pairing partner for i in structure structurenumber.

		//! \param i is the nucleotide index, which is one-indexed.
		//! \param structurenumber is the structure number, which is one-indexed.
		//! \return The pairing partner, as a nucleotide index.
		int GetPair(int i, int structurenumber=1);

		//! Get the label associated with the sequence.

		//! \return The string read from the sequence file.
		string GetSequenceLabel();

		//! Get the length of the sequence.

		//! \return An integer that is the sequence length.
		inline int GetSequenceLength() {
			//Return the value of numofbases:
			return numofbases;
		}
		
		//! Remove the pair at index i for structure number structurenumber.

		//! If i is paired to j, the pairing for j is also removed.
		//! \param i is the nucleotide for whic pairing is to be removed.  This is one-indexed.
		//! \param structurenumber is an int that provides from which structure the pair should be removed.  This is one-indexed.
		void RemovePair(int i, int structurenumber=1);

		//! Set the label for a structure, using a string.

		//! \param label is a string that will be strored.
		//! \param structurenumber is the index to which structure will hold the label.  This is one-indexed.
		void SetCtLabel(string label, int structurenumber);
		

		//! Set the label for a structure,using a pointer to cstring.

		//! \param label is a pointer to char that provides the label.
		//! \param structurenumber is the index to which structure will hold the label.
		void SetCtLabel(char *label, int structurenumber);

		//! Set the energy for structure numer structurenumber.

		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \param energy is an int that sets the energy in kcal/mol*conversionfactor, where conversionfactor is set in /src/defines.h.
		void SetEnergy(int structurenumber, int energy);

		
		//Set the pairing partner for i in structure structurenumber.

		//! This function sets nucleotide i paired to j in structure number structurenumber.
		//! \param i is the nucleotide index of the first pairing partner, which is one-indexed.
		//! \param j is the nucleotide index of the second pairing partner.
		//! \param structurenumber is the structure number, which is one-indexed.
		void SetPair(int i, int j, int structurenumber=1);

		//! Set the label from a sequence, using a string.

		//! \param label is a string that will be stored.
		void SetSequenceLabel(string label);


		//********************************
		// Get and Set constraint information
		//*********************************
		

		//! Allocate a bool array with information about whether any given pair is allowed.

		//! This function must be called after reading a sequence and before any template information is read or written.
		//! This mechanism is orthogonal to the functions AddForbiddenPair, GetForbiddenPair5, and GetForbiddenPair3.  It exists for the convenience of coding functiopns that need to forbid a large number of pairs.
		//! The memory use is cleaned up in the destructor.
		void allocatetem();


		//! Add a nucleotide to the list of those that must pair.

		//! \param i is an int that indicates the nucleotide position, one indexed. 
		void AddDouble(int i);

		//! Add a pair of nucleotides to the list of those not allowed to form.

		//! \param i is an int that indicates the 5' nucleotide position, one indexed.
		//! \param j is an int that indicates the 3' nucleotide position, one indexed.
		void AddForbiddenPair(int i, int j);

		//! Add a nucleotide to the list of Us in GU pairs.

		//! Note that there is no error checking.  This nucleotide must be a U.
		//! \param i is an int that indicates the nucleotide position, one indexed.
		void AddGUPair(int i);

		//! Add a nucleotide to the list of those accessible to traditional chemical modification.

		//! These nucleotides can only be unpaired, at the end of a helix, in a GU pair, or adjacent to a GU pair.
		//! \param i is an int that indicates the nucleotide position, one indexed.
		void AddModified(int i);

		//! Add a pair of nucleotides to the list of those that must form.

		//! Note that there is no error checking.  This should be an allowed base pair.
		//! \param i is an int that indicates the 5' nucleotide position, one indexed.
		//! \param j is an int that indicates the 3' nucleotide position, one indexed.
		void AddPair(int i, int j);

		//! Add a nucleotide to the list of those not able to pair.

		//! \param i is an int that indicates the nucleotide position, one indexed. 
		void AddSingle(int i);

		//!	Indicate if pairing distance is limited for structrure prediction methods.
//>>>>>>> ../../RNAstructure-object/src/structure.h

		//! \return A bool that is true if the pairing distance has a limit.
		inline bool DistanceLimited() {

			return limitdistance;

		}

		//! Get a nucleotide that must be base paired.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofDoubles()-1, inclusive.
		int GetDouble(int i);

		//! Get a nucleotide that must not be in a specific pair.

		//! \return An int that gives the nucleotide position for the 5' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofForbiddenPairs()-1, inclusive.
		int GetForbiddenPair5(int i);

		//! Get a nucleotide that must not be in a specific pair.

		//! \return An int that gives the nucleotide position for the 3' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofForbiddenPairs()-1, inclusive.
		int GetForbiddenPair3(int i);

		//! Get a nucleotide that must be a U in a GU pair.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofGU()-1, inclusive.
		int GetGUpair(int i);

		//! Get a nucleotide that is accessible to chemical modification.

		//! These nucleotides can only be unpaired, at the end of a helix, in a GU pair, or adjacent to a GU pair.
		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofMofidied()-1, inclusive.
		int GetModified(int i);

		//! Get the number of nucleotides forced to be double-stranded.

		//! \return An int that is the number of nucleotides constrained to be double stranded.
		int GetNumberofDoubles();

		//! Get the number of pairs that are forbidden.

		//! \return An int that is the number of forbidden pairs.
		int GetNumberofForbiddenPairs();
		
		//! Get the number of Us forced to be in GU pairs.

		//! \return An int that is the number of nucleotides constrained to be in GU pairs.
		int GetNumberofGU();
		
		//! Get the number of nucleotides that are accessible to chemical modification.

		//! \return An int that is the number of nucleotides constrained to be chemically modified.
		int GetNumberofModified();
		
		//! Get the number of nucleotides forced to be single-stranded.

		//! \return An int that is the number of nucleotides constrained to be single stranded.
		int GetNumberofSingles();
		
		//! Get the number of pairs that are constrained to occur.

		//! \return An int that is the number of forced pairs.
		int GetNumberofPairs();
		
		//! Get a nucleotide that must be in a specific pair.

		//! \return An int that gives the nucleotide position for the 5' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofPairs()-1, inclusive.
		int GetPair5(int i);

		//! Get a nucleotide that must be in a specific pair.

		//! \return An int that gives the nucleotide position for the 3' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofPairs()-1, inclusive.
		int GetPair3(int i);

		//!	Provide the maximum distance between nucleotides that can pair.

		//! \return An int that is the maximum distance in sequence for nucleotides that can pair.
		inline int GetPairingDistanceLimit() {

			return maxdistance;

		}

		//! Get a nucleotide that must be single stranded.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofSingles()-1, inclusive.
		int GetSingle(int i);
		

		//! Reset, i.e. remove, all constraints
		void RemoveConstraints();


		//! Set a maximum pairing distance between nucleotides that can pair.  

		//! For nucleotides i and j, if |j-i|>=ct->maxdistance, the nucleotides cannot pair.
		//! This also sets a bool so that DistanceLimited() will return true.
		//! \param maxdistance is an int that will be the maximum distance.
		void SetPairingDistance(int maxdistance);


		//*********************************
		//Functions for disk I/O
		//*********************************

		//! Write a ct file to disk.

		//!	if append is set to true, the ct information is appended to an existing file, if the file exists.
		//! \param ctoutfile is a const char pointer to a Null-terminated cstring that provides a filename.
		//! \param append is a bool that indicates if these structures should be appended to the end of the file.  The default, falase, is to overwrite any existing file.
		void ctout (const char *ctoutfile, bool append=false);		


		//! Open a CT File.

		//! This opens a ct file and stores all the information in this instance of structure.
		//! A non-zero return indicates and error in reading the file.  The value of the return is the line number that contains an error.
		//! \param ctfile is a pointer to a Null-terminated cstring that gives the filename, including any necessary path information.
		//! \return An int that gives the linenumber of an error or -1 if the file was not found.
		long openct(const char *ctfile);


		//! Open a sequence file.

		//! This function works on both .seq and FASTA files.

		//! \param seqfile is a const char pointer to a cstring that gives the filename, including any path information.
		//! \return An int that indicates an error state: 1 on no error and 0 on error.
		int openseq (const char *seqfile);

		//! Write a .bracket file (Vienna Format).

		//! The file will be decipherable only if there are no pseudoknots in the structure.  There is no error checking on this.
		//! \param filename is a const char pointer to a Null-terminated cstring that provides a filename.
		void writedotbracket(const char *filename);

		//*******************************
		//Functions that act on whole structures
		//*******************************

		//! Add another empty structure to the list of singlestructures.
		// DHM: Remember if this is structure 1 to set the label from some sequence label.!
		void AddStructure();

		//! Remove all pairs from a structure, i.e. make it a clean slate for a new set of pairs
		
		//! \param struturenumber is an index to the structure to be cleaned.
		void CleanStructure(int structurenumber);

		

		//! Remove the last structure.
		void RemoveLastStructure();

		//! Remove the structure at structurenumber.
		
		//! If the last structure is being removed, it is more efficient to use RemoveLastStructure();
		//! \param structurenumber is an int that is the index to which structure should be removed.  This is one indexed.
		void RemoveStructure(int structurenumber);

		//********************************
		//Additional functions
		//********************************

		//! Sort structures by energy.

		//! This function sorts structures in energy from lowest to highest.
		//! It is important that the structure energies be present by structure prediction or by efn2.
		void sort();


		//! Find problem in the set of structures.

		//! This function is for debugging.  It checks each pair in each structure to look for inconsistencies. 
		//! \return true when an inconstency with pairing is found and false otherwise.
		bool ProblemwithStructures();

		
		void ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file
		//void ReadSHAPE(const char *filename, bool calculate=true);//Read SHAPE reactivity data from a file
		void ReadSHAPE(const char *filename, std::string modifier="SHAPE", bool calculate=true, bool nosum=false);//Read SHAPE reactivity data from a file
		void ReadOffset(const char *SSOffset, const char *DSOffset);//Read Free Energy Offset Files.
		void ReadExperimentalPairBonus(const char *filename, double const experimentalOffset = 0.0, double const experimentalScaling = 1.0 );

		

		double **constant;//constant is used to hold an array of equilibrium constants.  In partition function calculations, 
					//the equilibrium constant is multiplied by constant[j][i] when the i-j pair is formed. 
					//NOTE: The use of constant is NOT orthogonal to using chemical modification data.  They cannot
					//both be used at once.
		void allocateconstant();//Function to allocate memory for constant array.
		bool SHAPEFileRead;


	//private:
		
		
		
		

		string sequencelabel;//a label that was read from disk along with a sequence

		
		short int *numseq,*hnumber;
		
		int inter[3],allocatedstructures;
		char *nucs;
		bool intermolecular,allocated,templated,stacking;
		bool **tem;//tem stores template information as to whether a pair is allowed

		void allocate(int size = maxbases);
		void allocatestructure(int structures);
		
		
		
		 
		short int min_gu, min_g_or_u;//NMR-derived constraint variables
		short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
		//regional NMR constraints:
		short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
		short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
		//microarray type constraints:
		short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
		bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
		
		
		double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
		double **EX;// double array that contains experimental bonuses/penalties
		bool shaped;//keeps track of whether SHAPE data was loaded
		bool experimentalPairBonusExists;//keeps track of whether experimental bonus data was loaded
		bool ssoffset;//keeps track of whether a single stranded offset was read from disk
		double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability
		//SINGLE STRANDED SHAPE ENERGY VARIABLES AND FUNCTIONS
		double *SHAPEss; //short int array that contains SHAPE data for single-stranded segments
		double SHAPEslope_ss, SHAPEintercept_ss; //values of the slope and intercept for SHAPE data modifying single stranded loop stability
		short int **SHAPEss_region;  //2-d short int array containing energy values for hairpin loop combinations
		int SHAPEss_calc(int index_i, int index_j);  //Returns pseudoenergy term for a hairpin loop using single stranded SHAPE data
		short int SHAPEss_give_value(int index);  //Returns the single stranded SHAPE pseudo energy for a given nucleotide
		double CalculatePseudoEnergy(double data, std::string modifier, double, double);
		double Gammadist(double data, double shape, double loc, double scale);
		double Potential(double data, std::vector< std::vector<double> > params, double kT);
		void ReadProbabilisticPotentialParams();//Read chemical modifier distributions from file
		//Parameters for distributions
		std::vector< std::vector<double> > SHAPE_params;
		std::vector< std::vector<double> > DMS_params;
		std::vector< std::vector<double> > CMCT_params;
		bool distsread;//keep track if the distruibution files have been read from disk.


	private:
		
		int numofbases;//number of nucleotides in sequence
		bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
		int maxdistance;//maximum distance between nucs in base pairs
		
		vector<singlestructure> arrayofstructures;//This holds an array of structures, i.e. base pairing information and comments
			
		//variables for holding folding constraints:
		vector<int> doublestranded; //nucleotides that must be double stranded
		vector<int> singlestranded; //nucleotides that must be single stranded
		vector<int> GUpair; //Us in GU pairs
		vector<int>	modified; //nucleotides accessible to tradictional chemical modification agents
		vector<int> pair5; //5' partner in forced pair
		vector<int> pair3; //3' partner in forced pair
		vector<int> forbid5; //5' partner in a forbidden pair
		vector<int> forbid3; //3' partner in a forbidden pair

		

};

short int tonumi(char *base); //converts base to a numeric


char *tobase (int i);//convert a numeric value for a base to the familiar
								//character


void tonum(char *base,structure *ct,int count); //converts base to a numeric



int ecompare(const void *i, const void *j);

void swap(int *a,int *b);//Swap two variables
void swap(short int *a,short int *b);//swap two variables
void swap(float *a,float *b);//swap two variables
void swap(double *a,double *b);//swap two variables



#endif
