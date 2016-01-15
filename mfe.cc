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
#include <fstream>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <sstream>

#include "energy_constant.h"

#include "mfe.h"
#include "ifold.h"


//minimum unpaired nts in a loop

DEFINE_string(InputFile, "", "File with input data");
DEFINE_string(RNAscdstr, "", "Input RNA secondary structure");
DEFINE_int64(TimeLimit, DEFAULT_TIME_LIMIT, "Search time limit in seconds");
DEFINE_int32(MAXsol, DEFAULT_MAX_SOLUTIONS, "Input MAX solutions");

DEFINE_string(RNAseqcon, "", "Input RNA sequence constraints");
DEFINE_string(AAseqcon, "", "List of amino acid sequence constraints");
DEFINE_string(AAtarget, "", "List of target amino acid sequence for blosum similarity score maximization");
DEFINE_string(AAstartPos, "", "List of starting positions for amino acid constraints");

DEFINE_int32(AAsimilCstr, 5, "Blosum similarity threshold (-4 to 4)(default 5 forces aminoacids to be equal to those specified) or allow similar amino acids based on a specific classification (6- Classification 1, 7-Classification 2 (see manual)");
DEFINE_int32(MaxBlosumScore, 0, "Maximize blosum score of target sequences (default disabled)");


DEFINE_string(RNAcompstr, "", "Compatible RNA secondary structure");
DEFINE_string(IncompBP, "", "List of incompatible base pairs");
DEFINE_string(temp, DEFAULT_TEMP, "Target folding temperature(s)");


DEFINE_int32(IncludeDangles, 1, "Include dangling positions when creating helices"); 
DEFINE_int32(HelixHeuristic, HH_OVERLAP_BP, "Helix ordering heuristic for the search 1-Simple overlap 2-Base pair overlap 3-Total overlap (default 2)");
DEFINE_int32(VarHeuristic, SH_BOTTOM_TO_TOP, "Variable heuristic 1-Helices bottom to top 2-In to out ");
DEFINE_int32(UPthreshold, DEFAULT_UP_THRESHOLD, "Probablility of selecting the next UP assignment in value heuristic");
DEFINE_int32(BPthreshold, DEFAULT_BP_THRESHOLD, "Probablility of selecting the next BP assignment in value heuristic");
DEFINE_int32(RandomAssignment, 0, "Activate random value heuristic");
DEFINE_int32(dangles, DEFAULT_DANGLES, "Dangling treatment");

DEFINE_string(RNAlibrary, DEFAULT_RNA_LIBRARY, "RNA library for folding and computing energy values: Allowed models are: Vienna (ViennaRNA package) RNAstructure (Mathews' Lab). Default RNA library is ViennaRNA (Vienna)");
DEFINE_string(EnergyModel, DEFAULT_ENERGY_MODEL, "Energy model for ViennaRNa library. Allowed models are: 2004 (Turner '04), 1999 (Turner '99), 2007 (Andronescu '07) . Default energy model is Turner '04 (2004)");

DEFINE_double(minGCcont, 0, "Minimum GC content");
DEFINE_double(maxGCcont, 100.0, "Maximum GC content");
DEFINE_int32(minAU, 0, "Minimum number of AU base pairs");
DEFINE_int32(maxAU, -1, "Maximum number of AU base pairs");
DEFINE_int32(minGC, 0, "Minimum  number of GC base pairs");
DEFINE_int32(maxGC, -1, "Maximum number of GC base pairs");
DEFINE_int32(minGU, 0, "Minimum  number of GU base pairs");
DEFINE_int32(maxGU, -1, "Maximum number of GU base pairs");

DEFINE_string(minA, "", "List of minimum number of As in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(maxA, "", "List of maximum number of As in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(minC, "", "List of minimum number of Cs in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(maxC, "", "List of maximum number of Cs in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(minG, "", "List of minimum number of Gs in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(maxG, "", "List of maximum number of Gs in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(minU, "", "List of minimum number of Us in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(maxU, "", "List of maximum number of Us in the full sequence (N) or in a specific range (N StartPos EndPos)");


DEFINE_string(consA, "", "List of maximum number of consecutive As in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(consC, "", "List of maximum number of consecutive Cs in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(consG, "", "List of maximum number of consecutive Gs in the full sequence (N) or in a specific range (N StartPos EndPos)");
DEFINE_string(consU, "", "List of maximum number of consecutive Us in the full sequence (N) or in a specific range (N StartPos EndPos)");

DEFINE_int32(MFEstructure, 1, "Sequence MFE structure(s) must be the target structure(s)");
DEFINE_double(MaxMFE, NO_ENERGY_LIMIT, "Maximum free energy allowed for the sequence folded into its minimum free energy structure");
DEFINE_double(MinMFE, NO_ENERGY_LIMIT, "Minimum free energy allowed for the sequence folded into its minimum free energy structure");
DEFINE_int32(MinimizeMFE, 0, "Minimize free energy of the MFE structure: print only sequences whose free energy when folded into their MFE structure is equal or lower than the previous one");
DEFINE_int32(MinimizeEnsDef, 0, "Minimize ensemble defect for the target structure: print only sequences whose ensemble defect for the target structure is equal or lower than the previous one");

DEFINE_int32(ShowHelices, 0, "Show helix identifiers for local constraints");
DEFINE_string(HelixCstrs, "", "Helix local constraints");
DEFINE_string(LocalCstrs, "", "Local constraints");


DEFINE_int32(LNS, 0, "Activate Large Neighborhood Search");
DEFINE_int32(LNSunchangedRestarts, DEFAULT_LNS_UNCHANGED_RESTARTS, "Maximum consecutive restarts in LNS without changes on fixed positions");
DEFINE_int32(LNSrestartTime, DEFAULT_LNS_RESTART_TIME, "Maximum consecutive restarts in LNS without changes on fixed positions");
DEFINE_int32(LNStimeMultiplier, DEFAULT_LNS_RESTART_TIME_MULTIPLIER, "Multiplier for restart time in milliseconds (Time=length*multiplier)");

DEFINE_int32(ShowMeasures, 1, "Show structural diversity measures");

using namespace std;

/* BEGIN OR-TOOLS FUNCTION */

namespace operations_research {

	void IfoldCp(std::vector<int*> strs_int, std::vector<int*> int_strs_undet, int n, int maxSolutions, int64 timeLimit, char* sequence, vector<AAConstraint*> aaConstraints, int helixHeuristic, int varHeuristic, int randomAssignment, int upthreshold, int bpthreshold, int includeDangles, int dangles, std::string rnaLib, std::string energyModel,std::vector<double> foldTemps, double minGCcont, double maxGCcont, int minAU, int maxAU, int minGC, int maxGC, int minGU, int maxGU,std::vector<std::tuple<int,int,int>> listMinA,std::vector<std::tuple<int,int,int>> listMaxA,std::vector<std::tuple<int,int,int>> listMinC,std::vector<std::tuple<int,int,int>> listMaxC,std::vector<std::tuple<int,int,int>> listMinG,std::vector<std::tuple<int,int,int>> listMaxG,std::vector<std::tuple<int,int,int>> listMinU,std::vector<std::tuple<int,int,int>> listMaxU,std::vector<std::tuple<int,int,int>> listConsA,std::vector<std::tuple<int,int,int>> listConsC,std::vector<std::tuple<int,int,int>> listConsG,std::vector<std::tuple<int,int,int>> listConsU, int MFEstructure, int minimizeMFE, int minimizeEnsDef, int* comp_str_int, vector<pair<int,int>> vIncompBP, int showHelices, std::vector<HelixCstr> helixCstrs, std::vector<LocalCstr> localCstrs, int LNS, int LNSunchangedRestarts, int LNSrestartTime, int showMeasures, int cutPoint){
		int i;

		// Constraint programming engine
		Solver solver("Inverse Folding CP");
		clock_t start_init = clock();
//		cout << "Creating iFold"<< endl;
		IFold ifold = IFold(&solver, n, dangles,rnaLib, energyModel);
//		cout << "Initializing domains"<< endl;
		ifold.InitDomains(sequence);

//		cout << "Adding structure constraints"<< endl;
		for(i=0;i<strs_int.size();i++){
			ifold.AddStructureConstraints(strs_int[i],int_strs_undet[i],includeDangles,foldTemps[i],cutPoint, MFEstructure, showHelices, helixCstrs);
		}
		// If showHelices option id active finish after printing the trees
		if(showHelices){
			return;
		}

//		cout << "Adding local structure constraints"<< endl;		
		ifold.AddLocalStructureConstraints(foldTemps[0],cutPoint,localCstrs);
				
//		cout << "Adding amino acid constraints"<< endl;
		for(i=0;i<aaConstraints.size();i++){
			ifold.AddAminoAcidConstraint(aaConstraints[i]);
		}
//		cout << "Adding nucleotide range constraints"<< endl;
		ifold.AddNucleotideRangeConstraints(minGCcont,maxGCcont,minAU,maxAU,minGC,maxGC,minGU,maxGU,listMinA,listMaxA,listMinC,listMaxC,listMinG,listMaxG,listMinU,listMaxU,listConsA,listConsC,listConsG,listConsU);

//		cout << "Adding compatibile structure constraint"<< endl;
		ifold.AddCompatibilityConstraint(comp_str_int);

//		cout << "Adding incompatibile base pairs constraint"<< endl;
		ifold.AddIncompatibilityConstraint(vIncompBP);

//		cout << "Selecting search heuristic"<< endl;
		ifold.SetSearchHeuristic(helixHeuristic,varHeuristic, randomAssignment, upthreshold, bpthreshold);
		clock_t end_init = clock();
		float seconds = (float)(end_init - start_init) / CLOCKS_PER_SEC;
		cout << "Init time: " << seconds << endl;

//		cout << "Starting search"<< endl;			
		ifold.Search(maxSolutions, timeLimit*1000, minimizeMFE, minimizeEnsDef,LNS,LNSunchangedRestarts, LNSrestartTime, showMeasures);
		clock_t end_search = clock();
		seconds = (float)(end_search - end_init) / CLOCKS_PER_SEC;
		cout << "Search time: " << seconds << endl;
		return;
	}

}   // namespace operations_research

/* END OR-TOOLS FUNCTION */


static string GetProgramUsage(const char *argv0) {
	std::ostringstream oss;
	oss << endl << RNAIFOLD_NAME<<" " << RNAIFOLD_VERSION << ": Constraint programming software for complete RNA inverse folding." <<endl;
	oss << "Copyright (C) 2014 Ivan Dotú, Juan Antonio García Martín, Peter Clote " << endl<< endl;
	oss << RNAIFOLD_NAME << " finds one or more sequences whose MFE structure is the target structure." <<endl<< endl;
	oss << "INPUT:" << endl << "  Only a target secondary structure is required, which can be provided either as a parameter by command line flags or inside an input file." << endl;
	oss << "   " << argv0 << " -RNAscdstr <SECONDARY_STRUCTURE>" << endl;
	oss << "   " << argv0 << " -InputFile <INPUT_FILE>" << endl << endl;
	oss << "   " << "Design constraints can be also provided as command line flags or inside the input file using the appropriate label preceeded by the \"pound\" symbol (\"#\") and writing the value in the next line. " << endl;
	oss << "   Valid labels/flags are: " << endl;
	vector<google::CommandLineFlagInfo> flags;
	google::GetAllFlags(&flags);           // flags are sorted by filename, then flagname

	for (vector<google::CommandLineFlagInfo>::const_iterator flag = flags.begin();flag != flags.end();++flag) {
		if (flag->filename == "mfe.cc") {
			oss << "     -"<< flag->name << ": (" << flag->type <<") " << flag->description << ". DEFAULT:"<< flag->default_value << endl; 
//			oss << DescribeOneFlag(*flag).c_str();
		}
	}
	oss << endl<< "   Input File format:" << endl;
	oss << "     Input file must contain a valid secondary structure, all the other fields are optional, "<< RNAIFOLD_NAME<<" input file format is:" << endl;
	oss << "       > Fasta comment" << endl;
	oss << "       Target structure" << endl;
	oss << "       Sequence constraints" << endl;
	oss << "       # Parameter" << endl;
	oss << "       Parameter value" << endl;
		
	oss << endl<< "OUTPUT:" << endl;
	oss <<  "  Three possible types of results can be returned:" << endl;
	oss <<  "    - Solution found: For each solution found the following information is displayed." << endl;
	oss <<  "        Sequence." << endl;
	oss <<  "        GC content and the number of base pairs of each type (strong, weak and wobble)." << endl;
	oss <<  "        Amino acid sequence and blosum score if blosum maximization is active." << endl;
	oss <<  "        Free energy of the structure in kcal/mol" << endl;
	oss <<  "        Additional measures" << endl;
	oss <<  "    - No solution found: If search time is reached and no solution has been found within this time limit." << endl;
	oss <<  "    - No possible solution: If the target structure (with specified constraints) has no solution and the time limit has not been reached." << endl;
	
	oss << endl<< "EXAMPLES:" << endl;
	oss <<  "  "<< argv0 << " -InputFile examples/tRNA.fas" << endl;	
	oss <<  "  "<< argv0 << " -RNAscdstr '(((((.(..((((.........)))).(((((.......))))).....(((((.......)))))).))))).' -RNAseqcon 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCUGGCUCG' -MAXsol 5 -dangles 2 -minGC 10 -maxGC 20 -LNS 1" << endl;
	
	
	return oss.str();
}

// MAIN function
int main(int argc, char *argv[]){
	google::SetUsageMessage(GetProgramUsage(argv[0]));
	google::SetVersionString(RNAIFOLD_VERSION);
	char sequence[MAXSIZE];
	vector<string> structures;
	vector<int*> int_strs;
	vector<int*> int_strs_undet;
	vector<double> foldTemps;

	int* comp_str_int = NULL;	

	int i,len;

	google::ParseCommandLineFlags(&argc, &argv, true);
	string inputFile = FLAGS_InputFile;
	string strs = FLAGS_RNAscdstr;

	int maxSolutions = FLAGS_MAXsol;
	int64 timeLimit = FLAGS_TimeLimit;

	
	string seqConst = FLAGS_RNAseqcon;
	std::transform(seqConst.begin(), seqConst.end(),seqConst.begin(), ::toupper);
	
	string aaConst = FLAGS_AAseqcon;	
	std::transform(aaConst.begin(), aaConst.end(),aaConst.begin(), ::toupper);	
		

	string aaTarget = FLAGS_AAtarget;	
	std::transform(aaTarget.begin(), aaTarget.end(),aaTarget.begin(), ::toupper);	
	
	string aaStartPos = FLAGS_AAstartPos;


	int aaSimilCstr = FLAGS_AAsimilCstr;
	int maxBlosumScore = FLAGS_MaxBlosumScore;
	
	string compStr = FLAGS_RNAcompstr;
	string incompBP = FLAGS_IncompBP;


	int includeDangles = FLAGS_IncludeDangles;

	int helixHeuristic = FLAGS_HelixHeuristic;
	int varHeuristic = FLAGS_VarHeuristic;
	
	int randomAssignment = FLAGS_RandomAssignment;
	int upthreshold = FLAGS_UPthreshold;
	int bpthreshold = FLAGS_BPthreshold;

	string rnaLib =FLAGS_RNAlibrary;
	string energyModel = FLAGS_EnergyModel;

	int dangles = FLAGS_dangles;
	string strTemps = FLAGS_temp;	
	
	double minGCcont = FLAGS_minGCcont;
	double maxGCcont = FLAGS_maxGCcont;
	int minAU =  FLAGS_minAU;
	int maxAU =  FLAGS_maxAU;
	int minGC =  FLAGS_minGC;
	int maxGC =  FLAGS_maxGC;
	int minGU =  FLAGS_minGU;
	int maxGU =  FLAGS_maxGU;

	string minA =  FLAGS_minA;
	string maxA =  FLAGS_maxA;
	string minC =  FLAGS_minC;
	string maxC =  FLAGS_maxC;
	string minG =  FLAGS_minG;
	string maxG =  FLAGS_maxG;
	string minU =  FLAGS_minU;
	string maxU =  FLAGS_maxU;

	string consA =  FLAGS_consA;
	string consC =  FLAGS_consC;
	string consG =  FLAGS_consG;
	string consU =  FLAGS_consU;

	int MFEstructure = FLAGS_MFEstructure;

	double upperMFE =  FLAGS_MaxMFE;
	double lowerMFE =  FLAGS_MinMFE;
	int minimizeMFE = FLAGS_MinimizeMFE;
	int minimizeEnsDef = FLAGS_MinimizeEnsDef;
	
	int showHelices = FLAGS_ShowHelices;
	string helixCstrsStr = FLAGS_HelixCstrs;
	string localCstrsStr = FLAGS_LocalCstrs;
	
	int LNS = FLAGS_LNS;
	int LNSunchangedRestarts = FLAGS_LNSunchangedRestarts;
	int LNSrestartTime = FLAGS_LNSrestartTime;
	int LNStimeMultiplier = FLAGS_LNStimeMultiplier;

	int showMeasures = FLAGS_ShowMeasures;

	// If input file exists, read parameters from file
	if(inputFile!=""){
		if(!readDatafromFile(inputFile,&strs,&maxSolutions,&timeLimit,&seqConst,&aaConst,&aaSimilCstr,&aaTarget, &aaStartPos, &maxBlosumScore, &compStr,&incompBP,&helixHeuristic,&varHeuristic,&randomAssignment,&upthreshold,&includeDangles,&bpthreshold,&dangles,&rnaLib,&energyModel,&strTemps,&minGCcont,&maxGCcont,&minAU,&maxAU,&minGC,&maxGC,&minGU,&maxGU,&minA,&maxA,&minC,&maxC,&minG,&maxG,&minU,&maxU,&consA,&consC,&consG,&consU,&MFEstructure,&upperMFE,&lowerMFE,&minimizeMFE,&minimizeEnsDef,&showHelices,&helixCstrsStr,&localCstrsStr,&LNS,&LNSunchangedRestarts, &LNSrestartTime, &LNStimeMultiplier)){
			cout << "Cannot open "<< inputFile << endl;
			exit(1);
		}
	}	

	// Check that RNA structure exists
	if(strs.empty()){
		cout << google::ProgramUsage() <<  endl;	
		exit(1);
	}
	else{
		split(strs,STRUCTURE_DELIMITER, structures);
	}

	if(!strTemps.empty()){
		vector<string> tempList;
		split(strTemps,FIELD_DELIMITER,tempList);
		try{
			for(i=0; i<tempList.size();i++){
				foldTemps.push_back(stod(tempList[i]));
			}
		}
		catch (exception e){
			cout << google::ProgramUsage() <<  endl;	
			exit(1);
		}
	}
	
	
	// Find cut point for co-fold. Check all '&' symbols are at the same positions. Store cut point position in a variable and remove '&' symbol from structures and sequences
	int cutPoint=findAndRemoveCutPoints(&structures,&seqConst,&compStr,&incompBP);	
	if(cutPoint==-2){
			cout << "Indicator '&' of multiple structures is not at the same position in target structure and constraints" << endl;
			exit(1);
	}
	
	// Parse amino acid constraints 
	vector<AAConstraint*> aaConstraints;
	parseAAconstraints(&aaConstraints, aaTarget, aaConst,aaStartPos, maxBlosumScore, aaSimilCstr,structures[0].length());
	
	
	// Parse helix constraints
	std::vector<HelixCstr> helixCstrs;
	
	if(!helixCstrsStr.empty()){
		vector<string> vHelixCstrsStr;
		split(helixCstrsStr,FIELD_DELIMITER, vHelixCstrsStr);
		for(i=0; i<vHelixCstrsStr.size();i++){
			std::vector<string> helixCstrsFields;
			split(vHelixCstrsStr[i],' ',helixCstrsFields);
			if(helixCstrsFields.size()<3){
				cout << "Invalid syntax for for helix constraint \""<< vHelixCstrsStr[i]<<"\"" << endl <<" Correct syntax is: ID TYPE VALUE" <<  endl;	
				exit(1);
			}
			else{
				if(helixCstrsFields.size()==3){
					try{
						HelixCstr tmpCstr(stoi(helixCstrsFields[0]), helixCstrsFields[1], stod(helixCstrsFields[2]));
						helixCstrs.push_back(tmpCstr);
					}
					catch (exception e){
						cout << "Invalid syntax for for helix constraint \""<< vHelixCstrsStr[i]<<"\"" << endl <<" Correct syntax is: HELIX_ID TYPE VALUE" <<  endl;	
						exit(1);
					}
				}
				else if(helixCstrsFields.size()==4){
					try{
						HelixCstr tmpCstr(stoi(helixCstrsFields[0]), stoi(helixCstrsFields[1]), helixCstrsFields[2], stod(helixCstrsFields[3]));
						helixCstrs.push_back(tmpCstr);
					}
					catch (exception e){
						cout << "Invalid syntax for for helix constraint \""<< vHelixCstrsStr[i]<<"\"" << endl <<" Correct syntax is: STRUCTURE_ID HELIX_ID TYPE VALUE" <<  endl;	
						exit(1);
					}
				}
				else{
					cout << "Invalid syntax for for helix constraint \""<< vHelixCstrsStr[i]<<"\"" << endl <<" Accepted syntaxes are: "<< endl <<"  STRUCTURE_ID HELIX_ID TYPE VALUE"<< endl <<"  HELIX_ID TYPE VALUE" <<  endl;	
					exit(1);
				}
			}
		}
	}

	// Parse local constraints
	std::vector<LocalCstr> localCstrs;
	
	if(!localCstrsStr.empty()){
		vector<string> vLocalCstrsStr;
		split(localCstrsStr,STRUCTURE_DELIMITER, vLocalCstrsStr);
		for(i=0; i<vLocalCstrsStr.size();i++){
			std::vector<string> localCstrsFields;
			split(vLocalCstrsStr[i],' ',localCstrsFields);
			if(localCstrsFields.size()<3){
				cout << "Invalid syntax for local constraint \""<< vLocalCstrsStr[i]<<"\"" << endl <<" Accepted syntaxes are: "<< endl <<"  START_POS STRUCTURE TYPE VALUE"<< endl <<"  STRUCTURE TYPE VALUE" <<  endl;	
				exit(1);
			}
			else{
				if(localCstrsFields.size()==3){
					try{
						LocalCstr tmpCstr(localCstrsFields[0], localCstrsFields[1], stod(localCstrsFields[2]));
						localCstrs.push_back(tmpCstr);
					}
					catch (exception e){
						cout << "Invalid syntax for for local constraint \""<< vLocalCstrsStr[i]<<"\"" << endl <<" Correct syntax is: STRUCTURE TYPE VALUE" <<  endl;	
						exit(1);
					}
				}
				else if(localCstrsFields.size()==4){
					try{
						LocalCstr tmpCstr(stoi(localCstrsFields[0]), localCstrsFields[1], localCstrsFields[2], stod(localCstrsFields[3]));
						localCstrs.push_back(tmpCstr);
					}
					catch (exception e){
						cout << "Invalid syntax for for local constraint \""<< vLocalCstrsStr[i]<<"\"" << endl <<" Correct syntax is: START_POS STRUCTURE TYPE VALUE" <<  endl;	
						exit(1);
					}
				}
				else{
					cout << "Invalid syntax for for local constraint \""<< vLocalCstrsStr[i]<<"\"" << endl <<" Accepted syntaxes are: "<< endl <<"  START_POS STRUCTURE TYPE VALUE"<< endl <<"  STRUCTURE TYPE VALUE" <<  endl;	
					exit(1);
				}
			}
		}
	}	
	// Input parameter validations
	if(!checkInputParameters(structures,maxSolutions,timeLimit,seqConst,aaConstraints,compStr,incompBP,helixHeuristic,varHeuristic,randomAssignment,upthreshold,bpthreshold,includeDangles,dangles,rnaLib,energyModel,foldTemps,minGCcont,maxGCcont,minAU,maxAU,minGC,maxGC,minGU,maxGU,MFEstructure, upperMFE,lowerMFE,minimizeMFE,minimizeEnsDef,showHelices, helixCstrs, localCstrs, LNS, LNSunchangedRestarts, LNSrestartTime, LNStimeMultiplier)){	
		exit(1);		
	}
	// If there are more structures than temperatures assume that the target temperature for the remaining structures is the first
	if(foldTemps.size() < structures.size()){
		for(int i=foldTemps.size(); i<structures.size() ;i++){
			foldTemps.push_back(foldTemps[0]);
		}	
	}
	
	// Set energy parameter file paths
	string base_path=getExecPath(string(argv[0]));
	if(rnaLib.compare(RNASTRUCTURE_LIB)==0){
		energyModel= base_path;
		energyModel.append("/");
		energyModel.append(RNASTRUCTURE_DIR);
		//cout<< "Energy Model: "<< energyModel<< endl;
	}
	else{
		if(energyModel.compare(TURNER_04_CODE)==0){
			energyModel= base_path;
			energyModel.append("/");
			energyModel.append(TURNER_04_FILE);
		}
		else if(energyModel.compare(TURNER_99_CODE)==0){
			energyModel= base_path;
			energyModel.append("/");
			energyModel.append(TURNER_99_FILE);
		}
		else if(energyModel.compare(ANDRONESCU_07_CODE)==0){
			energyModel= base_path;
			energyModel.append("/");
			energyModel.append(ANDRONESCU_07_FILE);
		}
	}

	// Convert maxMFE and minMFE to helix constraints (a constraint for all root helices)
	if(lowerMFE != NO_ENERGY_LIMIT){
		for(i=0; i< structures.size(); i++){
			HelixCstr tmpCstr(i, 0, HC_ENERGY_MIN_STR, lowerMFE);
			helixCstrs.push_back(tmpCstr);
		}
	}
	if(upperMFE != NO_ENERGY_LIMIT){
		for(i=0; i< structures.size(); i++){
			HelixCstr tmpCstr(i, 0, HC_ENERGY_MAX_STR, upperMFE);
			helixCstrs.push_back(tmpCstr);
		}
	}

	len = structures[0].length();


	//Parse and check nucleotide constraints
	vector<tuple<int,int,int>> listConsA;
	if(!parseTriple(consA, len, "consA", &listConsA)){
		exit(0);
	}
	vector<tuple<int,int,int>> listConsC;
	if(!parseTriple(consC, len, "consC", &listConsC)){
		exit(0);
	}
	vector<tuple<int,int,int>> listConsG;
	if(!parseTriple(consG, len, "consG", &listConsG)){
		exit(0);
	}
	vector<tuple<int,int,int>> listConsU;
	if(!parseTriple(consU, len, "consU", &listConsU)){
		exit(0);
	}	
	vector<tuple<int,int,int>> listMinA;
	if(!parseTriple(minA, len, "minA", &listMinA)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMaxA;
	if(!parseTriple(maxA, len, "maxA", &listMaxA)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMinC;
	if(!parseTriple(minC, len, "minC", &listMinC)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMaxC;
	if(!parseTriple(maxC, len, "maxC", &listMaxC)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMinG;
	if(!parseTriple(minG, len, "minG", &listMinG)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMaxG;
	if(!parseTriple(maxG, len, "maxG", &listMaxG)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMinU;
	if(!parseTriple(minU, len, "minU", &listMinU)){
		exit(0);
	}
	vector<tuple<int,int,int>> listMaxU;
	if(!parseTriple(maxU, len, "maxU", &listMaxU)){
		exit(0);
	}

	
	// Transform structure to base pair array
	for(i=0; i<structures.size(); i++){
		int_strs.push_back(make_BasePair_Table(structures[i], false));
		// Inthe structure contains undertermined positions create the base pair array
		if (structures[i].find(',') != std::string::npos){
			int_strs_undet.push_back(make_BasePair_Table(structures[i], true));
		}
		else{
			int_strs_undet.push_back(NULL);
		}
		// Print structure and base pair array
		cout << structures[i] << endl;
		//for(j=1; j<=len;j++){
		//	printf("%d\t",int_strs[i][j]);
		//}
		//printf("\n"); 

	}


	// COPY COMPATIBLE STRUCTURE
	if(compStr.compare("")!=0){
		if(compStr.length()!=len){
			cout << "Compatible structure and target structure length differ!" << endl;
			exit(1);
		}
		char compStructure[len];		
		strcpy(compStructure, compStr.c_str());
		comp_str_int = make_BasePair_Table(compStructure, false);
	}
	
	
	// Split incompatible base pairs into tokens 
	vector<pair<int,int>> vIncompBP;	
	if(incompBP.compare("")!=0){
		if(validSecondaryStructure(incompBP,false)!= -1){
			vector<int> openingPos;
			for(i=0; i<incompBP.length();i++){
				switch (incompBP[i]){
					case '(':
						openingPos.push_back(i);
						break;
					case ')':						
						vIncompBP.push_back(std::make_pair(openingPos.back()+1,i+1));
						openingPos.pop_back();
						break;
					default:
						break;
				}
			}

		}
		else{
			vector<string> listIncompBP;	

			char incompatibleBP[incompBP.length()];
			strcpy(incompatibleBP, incompBP.c_str());

			
			char* item;
			item = strtok (incompatibleBP,",");
			while (item != NULL)
			{
				listIncompBP.push_back(item);
				item = strtok (NULL, ",");
			}	


			for(i=0; i<listIncompBP.size();i++){
				char newItem[listIncompBP[i].length()];
				strcpy(newItem,listIncompBP[i].c_str());
				char* openBP;
				openBP = strtok (newItem," ");
				if(openBP==NULL){
					cout << "Wrong syntax in incompatible base pair ("<< newItem << ") - Syntax is (OpenBP CloseBP, ...)"<< endl; 
					exit(1);				

				}
				else{
					if(strcmp(openBP,"P")==0){
							openBP = strtok (NULL," ");	
							if(openBP!=NULL && is_number(openBP)){
								char* beginCloseBP = strtok (NULL," ");
								if(beginCloseBP!=NULL && is_number(beginCloseBP)){
									char* endCloseBP= strtok (NULL," ");
									if(endCloseBP!=NULL && is_number(endCloseBP)){
										int openBPList=atoi(openBP);
										int beginCloseBPList=atoi(beginCloseBP);
										int endCloseBPList=atoi(endCloseBP);										
										if(openBPList <=0 || beginCloseBPList+endCloseBPList-1>len || openBPList >= beginCloseBPList){
											cout << "Closing position (max: "<< len<<") must be higher than opening position (min: 1) in incompatible base pair ("<< newItem << ")" << endl; 
											exit(1);
										}
										else{
											for(int j=beginCloseBPList; j<beginCloseBPList+endCloseBPList;j++){
												pair<int,int> newIncomp = std::make_pair(openBPList,j);
												vIncompBP.push_back(newIncomp);
											}
										}
									}
								}
							}
							
					}
					else if(is_number(openBP)){
							char* closeBP = strtok (NULL," ");
							if(closeBP!=NULL && is_number(closeBP)){
								pair<int,int> newIncomp = std::make_pair(atoi(openBP),atoi(closeBP));
								if(newIncomp.first<=0 || newIncomp.second> len || newIncomp.first >= newIncomp.second){
									cout << "Closing position (max: "<< len<<") must be higher than opening position (min: 1) in incompatible base pair ("<< newItem << ")" << endl; 
									exit(1);
								}
								vIncompBP.push_back(newIncomp);	
							}							
					}
				}
			}
		}

// 		Show list of incompatible base pairs
//		for(i=0; i<vIncompBP.size();i++){
//			cout << vIncompBP[i].first << "-" <<vIncompBP[i].second << endl;
//		}

	}
	
	
	// Compute LNS restart time if the multiplier factor has been assigned
	if(LNStimeMultiplier>0){
		LNSrestartTime = (LNStimeMultiplier*len)/1000;
	}	

	// Read sequence constraints	
	sequence[0]='\0';
	if(seqConst!=""){
		if(seqConst.size() != len){
			cout << "Sequence constraint and target structure differ!" << endl;
		}
		else{
			strcpy(sequence, seqConst.c_str());
		}
	}
	
	// CALL IFOLD
	clock_t start = clock();
	operations_research::IfoldCp(int_strs,int_strs_undet,len, maxSolutions, timeLimit, sequence, aaConstraints, helixHeuristic,varHeuristic, randomAssignment, upthreshold, bpthreshold,includeDangles,dangles,rnaLib,energyModel,foldTemps, minGCcont, maxGCcont, minAU, maxAU, minGC, maxGC, minGU, maxGU,listMinA,listMaxA,listMinC,listMaxC,listMinG,listMaxG,listMinU,listMaxU,listConsA,listConsC,listConsG,listConsU, MFEstructure,minimizeMFE, minimizeEnsDef,comp_str_int, vIncompBP,showHelices,helixCstrs,localCstrs,LNS,LNSunchangedRestarts,LNSrestartTime,showMeasures,cutPoint);
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed " << seconds << endl;
	for(i=0;i<int_strs.size();i++){
		free(int_strs[i]);
		if(int_strs_undet[i]!=NULL){
			free(int_strs_undet[i]);
		}
	}
	
	return 0;
}

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}


bool checkInputParameters(std::vector<std::string> structures, int maxSolutions, int64 timeLimit, std::string seqConst,std::vector<AAConstraint*> aaConstraints, std::string compStr, std::string incompBP, int helixHeuristic, int varHeuristic, int randomAssignment, int upthreshold, int bpthreshold, int includeDangles, int dangles, std::string rnaLib, std::string energyModel,std::vector<double> foldTemps, double minGCcont, double maxGCcont, int minAU,int maxAU,int minGC,int maxGC,int minGU,int maxGU, int MFEstructure, double upperMFE, double lowerMFE, int minimizeMFE, int minimizeEnsDef,int showHelices, std::vector<HelixCstr> helixCstrs, std::vector<LocalCstr> localCstrs,int LNS, int LNSunchangedRestarts, int LNSrestartTime, int LNStimeMultiplier){
	std::string errorMessage ="";
	bool retVal = true;
	int numBasePairs = 0;
	int strLen=0;
	// SECONDARY STRUCTURE 
	for(int i=0; i< structures.size(); i++){
		if(structures[i].length()==0){
			errorMessage += "No target structure\n";
			retVal = false;
			cout << errorMessage;
			return retVal;
		}
		else{
			numBasePairs= validSecondaryStructure(structures[i],true);
			if(i==0){
				strLen = structures[i].length();
			}
			else{
				if(structures[i].length()!=strLen){
					errorMessage += "Secondary target structure and main target structure length differ!\n";
					retVal = false;
				}
			}

			if(numBasePairs==-1){
				errorMessage += "Unbalanced parenthesis in target structure\n";
				retVal = false;
			}
		}
	}
	
	// MAX SOLUTIONS
	if(maxSolutions < 0){
			errorMessage += "Number of solutions must be positive! (Complete search - 0)\n";		
			retVal = false;
	}
	
	// TIME LIMIT
	if(timeLimit < 0){	
			errorMessage += "Time limit must be positive! (No limit - 0)\n";		
			retVal = false;
	}

	// SEQUENCE CONSTRAINT 
	if(seqConst.length()>0){
		if(seqConst.length() != strLen){
				errorMessage += "Sequence constraint and target structure length differ!\n";
				retVal = false;
		}
		else{
			vector<char> validIUPAC = {'N','A','C','G','U','V','B','H','D','K','S','W','M','Y','R'}; 			
			for(int i=0; i<seqConst.length();i++){
				if(std::find(validIUPAC.begin(), validIUPAC.end(), seqConst[i]) == validIUPAC.end()) {
					errorMessage += "Invalid values in sequence constraint!\n";
					retVal = false;
					break;
				}
			}
		}
	}
	
	// AMINO ACID CONSTRAINTS AND BLOSUM MAXIMIZATION
	for(int i=0; i< aaConstraints.size();i++){
		if(aaConstraints[i]->isValid()!=AA_ERR_OK){
			errorMessage += aaConstraints[i]->getErrorMessage();
			retVal = false;
		}
	}


	// COMPATIBLE STRUCTURE 
	if(compStr.length()!=0){
		if(compStr.length()!=strLen){
			errorMessage += "Compatible structure and target structure length differ!\n";
			retVal = false;
		}
		else{
			if(validSecondaryStructure(compStr,false)==-1){
				errorMessage += "Unbalanced parenthesis in compatible structure\n";
				retVal = false;
			}
		}
	}
	
	// INCOMPATIBLE BASE PAIRS 
	if(incompBP.length()!=0){
		if(incompBP.length() !=strLen || validSecondaryStructure(incompBP,false)==-1){
			vector<string> listIncompBP;	
			
			char incompatibleBP[incompBP.length()];
			strcpy(incompatibleBP, incompBP.c_str());
			
			char* item;
			item = strtok (incompatibleBP,",");
			while (item != NULL)
			{
				listIncompBP.push_back(item);
				item = strtok (NULL, ",");
			}	
			int validList=1;
			for(int i=0; i<listIncompBP.size();i++){
				char newItem[listIncompBP[i].length()];
				strcpy(newItem,listIncompBP[i].c_str());
				char* openBP;
				openBP = strtok (newItem," ");
				if(openBP==NULL){
					validList=0;
					break;
				}
				else{
					if(strcmp(openBP,"P")==0){
							openBP = strtok (NULL," ");	
							if(openBP!=NULL && is_number(openBP)){
								char* beginCloseBP = strtok (NULL," ");
								if(beginCloseBP!=NULL && is_number(beginCloseBP)){
									char* endCloseBP= strtok (NULL," ");
									if(endCloseBP!=NULL && is_number(endCloseBP)){
										int openBPList=atoi(openBP);
										int beginCloseBPList=atoi(beginCloseBP);
										int endCloseBPList=atoi(endCloseBP);										
										if(openBPList <=0 || beginCloseBPList+endCloseBPList-1>strLen || openBPList >= beginCloseBPList){
												errorMessage += "Closing position  must be higher than opening position (min: 1) in incompatible base pair\n"; 
												validList=0;
												break;

										}										
									}
									else{
										validList=0;
										break;
									}									
								}
								else{
									validList=0;
									break;
								}								
							}
							else{
								validList=0;
								break;
							}		
							
					}
					else if(is_number(openBP)){
							char* closeBP = strtok (NULL," ");
							if(closeBP!=NULL && is_number(closeBP)){
								int openBPInt=atoi(openBP);
								int closeBPInt=atoi(closeBP);
								
								if(openBPInt<=0 || closeBPInt> strLen || openBPInt >= closeBPInt){
									errorMessage += "Closing position must be higher than opening position (min: 1) in incompatible base pair\n"; 
									validList=0;
									break;
								}
								
							}							
							else{
								validList=0;
								break;
							}		
					}
					else{
						validList=0;
						break;
					}
				}
			}
			if(validList==0){
				errorMessage += "Wrong syntax in incompatible base pair list - Syntax is (OpenBP CloseBP, ...)\n";
				retVal = false;		
			}
		}
	}
	
	
	// INCLUDE DANGLING POSITIONS
	if(includeDangles<0 || includeDangles > 1){
		errorMessage += "Invalid IncludeDangles flag value. [0|1]\n";
		retVal = false;		
	}

	//  HELIX HEURISTIC
	if(helixHeuristic<1 || helixHeuristic > 4){
		errorMessage += "Invalid helix heuristic. Allowed values(1 - SIMPLE_OVERLAP, 2 - BASE_PAIR_OVERLAP, 3 - TOTAL_OVERLAP, 2 - BASE_PAIR_PERCENT_OVERLAP\n";
		retVal = false;		
	}

	//  VARIABLE HEURISTIC
	if(varHeuristic<0 || varHeuristic > 3){
		errorMessage += "Invalid variable heuristic. Allowed values(0 - NONE, 1 - IN_TO_OUT, 2 - BOTTOM_TO_TOP, 3 - BOTTOM_TOP_UP\n";
		retVal = false;		
	}
	
	//  VALUE HEURISTICS
	if(randomAssignment<0 || randomAssignment > 1){
		errorMessage += "Invalid random assignment flag value. [0|1]\n";
		retVal = false;		
	}

	if(upthreshold<1 || upthreshold > 100){
		errorMessage += "UPthreshold must be an integer value between 1 and 100\n";
		retVal = false;		
	}
	if(bpthreshold<1 || bpthreshold > 100){
		errorMessage += "BPthreshold must be an integer value between 1 and 100\n";
		retVal = false;		
	}

	//  FOLDING PARAMETERS (DANGLES, ENERGY MODEL AND TEMPERATURE )
	if(dangles<0 || dangles > 3){
		errorMessage += "Invalid dangling treatment (0 to 3)\n";
		retVal = false;		
	}
	
	if(rnaLib.compare(VIENNA_LIB) != 0 && rnaLib.compare(RNASTRUCTURE_LIB) != 0){
		errorMessage += "Invalid RNA library. Allowed RNA libraries are \"Vienna\" and \"RNAstructure\"\n";
		retVal = false;		
	}
	else{	
		if(rnaLib.compare(VIENNA_LIB) == 0){
			if(energyModel.compare(TURNER_04_CODE) != 0 && energyModel.compare(TURNER_99_CODE) != 0 && energyModel.compare(ANDRONESCU_07_CODE) != 0){
				errorMessage += "Invalid energy model: Allowed models are: ";
				errorMessage += TURNER_04_CODE;
				errorMessage += " (Turner '04), ";
				errorMessage += TURNER_99_CODE;
				errorMessage += " (Turner '99), ";
				errorMessage += ANDRONESCU_07_CODE;
				errorMessage += " (Andronescu '07) \n";
				retVal = false;		
			}
		}
		else{
			energyModel=RNASTRUCTURE_DIR;
		}
	}
	
	for(int i=0; i< foldTemps.size();i++){
		if(foldTemps[i]<-273){
			errorMessage += "Invalid folding temperature (min: 273)\n";
			retVal = false;		
		}
	}
	if(foldTemps.size() > structures.size()){
		errorMessage += "More target folding temperatures than target structures provided\n";
		retVal = false;
	}
	
	
	//  GC CONTENT 
	if(minGCcont<0 || minGCcont > 100){
		errorMessage += "Minimum GC content must be between 0 and 100\n";
		retVal = false;		
	}

	if(maxGCcont<0 || maxGCcont > 100){
		errorMessage += "Maximum GC content must be between 0 and 100\n";
		retVal = false;		
	}
	if(minGCcont>maxGCcont){
		errorMessage += "Minimum GC content can not be higher than maximum GC content\n";
		retVal = false;		
	
	}
	
	// NUCLEOTIDE CONTENT 
	if(minAU<0 || minAU > numBasePairs){
		errorMessage += "Minimum number of AU must be between 0 and number of base pairs\n";
		retVal = false;		
	}
	if(maxAU<-1 || maxAU > numBasePairs){
		errorMessage += "Maximum number of AU must be between 0 and number of base pairs (-1 = no limit)\n";
		retVal = false;		
	}
	if(maxAU!= -1 && maxAU<minAU){
		errorMessage += "Minimum number of AU content can not be higher than maximum\n";
		retVal = false;		
	}

	if(minGC<0 || minGC > numBasePairs){
		errorMessage += "Minimum number of GC must be between 0 and number of base pairs\n";
		retVal = false;		
	}
	if(maxGC<-1 || maxGC > numBasePairs){
		errorMessage += "Maximum number of GC must be between 0 and number of base pairs (-1 = no limit)\n";
		retVal = false;		
	}
	if(maxGC!= -1 && maxGC<minGC){
		errorMessage += "Minimum number of GC content can not be higher than maximum\n";
		retVal = false;		
	}

	if(minGU<0 || minGU > numBasePairs){
		errorMessage += "Minimum number of GU must be between 0 and number of base pairs\n";
		retVal = false;		
	}
	if(maxGU<-1 || maxGU > numBasePairs){
		errorMessage += "Maximum number of GU must be between 0 and number of base pairs (-1 = no limit)\n";
		retVal = false;		
	}
	if(maxGU!= -1 && maxGU<minGU){
		errorMessage += "Minimum number of GU content can not be higher than maximum\n";
		retVal = false;		
	}

	//if(upperMFE> NO_ENERGY_LIMIT && upperMFE>0){
	if(upperMFE> NO_ENERGY_LIMIT){
		//errorMessage += "Maximum free energy allowed for the sequence folded into its minimum free energy structure must be 0 or less (1000 disabled)\n";
		errorMessage += "Maximum free energy allowed for the sequence folded into its minimum free energy structure must be less than 1000\n";
		retVal = false;		
	}

	//if(lowerMFE> NO_ENERGY_LIMIT && lowerMFE>0){
	if(lowerMFE> NO_ENERGY_LIMIT){
		//errorMessage += "Minimum free energy allowed for the sequence folded into its minimum free energy structure must be 0 or less (1000 disabled)\n";
		errorMessage += "Minimum free energy allowed for the sequence folded into its minimum free energy structure must be less than 1000\n";		
		retVal = false;		
	}
	if(lowerMFE != NO_ENERGY_LIMIT && upperMFE!=NO_ENERGY_LIMIT && upperMFE<lowerMFE){
		errorMessage += "Lower bound is higher than upper bound in free energy allowed for the sequence folded into its minimum free energy structure\n";
		retVal = false;			
	}
	
	// Minimum free energy constraint
	if(MFEstructure<0 || MFEstructure > 1){
		errorMessage += "Valid values for MFEstructure are 0 (disabled) and 1 (enabled))!\n";
		retVal = false;
	}	

	// Minimization constraints
	if(minimizeMFE<0 || minimizeMFE > 1){
		errorMessage += "Valid values for MinimizeMFE are 0 (disabled) and 1 (enabled))!\n";
		retVal = false;
	}

	if(minimizeEnsDef<0 || minimizeEnsDef > 1){
		errorMessage += "Valid values for MinimizeEnsDef are 0 (disabled) and 1 (enabled))!\n";
		retVal = false;
	}	

	// Helix constraints
	if(showHelices<0 || showHelices > 1){
		errorMessage += "Valid values for ShowHelices are 0 (disabled) and 1 (enabled))!\n";
		retVal = false;
	}
	
	for(int i=0; i<helixCstrs.size();i++){
		int helixType=helixCstrs[i].getLimitType();
		if(helixType == HC_WRONG_TYPE){
			errorMessage += "Invalid helix constraint type: Valid types are (";
			errorMessage += helixCstrs[i].printTypes();
			errorMessage += ")!\n";
			retVal = false;			
		}
		else if((helixType == HC_ED_MIN || helixType == HC_ED_MAX) && helixCstrs[i].getLimitValue() < 0){
			errorMessage += "Invalid value for helix ensemble defect constraint, ensemble defect must be positive!\n";
			retVal = false;			
		}
		else if(helixCstrs[i].getTreeId()>=structures.size()){
			//errorMessage += "Invalid structure in helix constraint ("<< helixCstrs[i].getTreeId() <<"). There are only "<< structures.size() << " structures!\n";
			errorMessage += "Invalid structure in helix constraint. Tree identifier is higher than the number of structures!\n";
			retVal = false;			
		}
	}
	
	// Local constraints
	for(int i=0; i<localCstrs.size();i++){
		int helixType=localCstrs[i].getLimitType();
		if(helixType == LC_WRONG_TYPE){
			errorMessage += "Invalid local constraint type: Valid types are (";
			errorMessage += localCstrs[i].printTypes();
			errorMessage += ")!\n";
			retVal = false;			
		}
		else if(helixType == LC_WRONG_STRUCTURE){
			errorMessage += "Invalid structure in local constraints!\n";
			retVal = false;			
		}
		else if((helixType == LC_ED_MIN || helixType == LC_ED_MAX) && localCstrs[i].getLimitValue() < 0){
			errorMessage += "Invalid value for local ensemble defect constraint, ensemble defect must be positive!\n";
			retVal = false;			
		}
		else if(localCstrs[i].getStartPos()+localCstrs[i].getStructureStr().length()-1>strLen){
			errorMessage += "Invalid helix constraint. Local constraint structure overflows target structure!\n";
			retVal = false;			
		}
	}
	
		
	// LNS FLAG
	if(LNS<0 || LNS > 1){
		errorMessage += "Invalid random LNS flag value. [0|1]\n";
		retVal = false;		
	}
	else if(LNS ==1){
		if(LNSunchangedRestarts<0){
			errorMessage += "Number of consecutive restarts in LNS without changes must be positive.\n";
			retVal = false;		
		}
		if(LNSrestartTime<1){
			errorMessage += "Restarts time for LNS must be positive.\n";
			retVal = false;		
		}
		if(LNStimeMultiplier<0){
			errorMessage += "Restarts time multiplier for LNS must be positive.\n";
			retVal = false;		
		}
	}
			
	if(!retVal){
		cout << errorMessage;
	}
	
	return retVal;
}


int findAndRemoveCutPoints(std::vector<std::string>* structures,std::string* seqConst,std::string* compStr,std::string* incompBP){
	int cutPoint=findCutPoint(&structures->at(0));
	
	if(cutPoint!=-1){
		if(seqConst->length()!=0){
			if(findCutPoint(seqConst)!=cutPoint){
				return -2;
			}
		}
		if(compStr->length()!=0){
			if(findCutPoint(compStr)!=cutPoint){
				return -2;
			}
		}
		if(incompBP->length()!=0){
			if(validSecondaryStructure(*incompBP, false)>=0){
				if(findCutPoint(incompBP)!=cutPoint){
					return -2;
				}
			}
		}
		for(int i=1; i<structures->size();i++){
			if(structures->at(i).length()!=0){
				if(findCutPoint(&structures->at(i))!=cutPoint){
					return -2;
				}
			}
		}

	}
	return cutPoint;

}

bool readDatafromFile(std::string inputFile, std::string* strs, int* maxSolutions, int64* timeLimit, std::string* seqConst,std::string* aaConst, int* aaSimilCstr,std::string* aaTarget, std::string* aaStartPos, int* maxBlosumScore, std::string* compStr, std::string* incompBP, int* helixHeuristic, int* varHeuristic, int* randomAssignment, int* upthreshold, int* bpthreshold, int* includeDangles, int* dangles, std::string* rnaLib, std::string* energyModel, string* strTemps, double* minGCcont, double* maxGCcont, int* minAU,int* maxAU,int* minGC,int* maxGC,int* minGU,int* maxGU,std::string* minA,std::string* maxA,std::string* minC,std::string* maxC,std::string* minG,std::string* maxG,std::string* minU,std::string* maxU,std::string* consA,std::string* consC,std::string* consG,std::string* consU, int* MFEstructure, double* upperMFE, double* lowerMFE, int* minimizeMFE, int* minimizeEnsDef, int* showHelices, std::string* helixCstrsStr, std::string* localCstrsStr, int* LNS, int* LNSunchangedRestarts, int* LNSrestartTime, int* LNStimeMultiplier){
	string line;
	ifstream infile;
	infile.open (inputFile);
	vector<char> validIUPAC = {'N','A','C','G','U','V','B','H','D','K','S','W','M','Y','R'}; 
	string currentParameter = ""; 
	if(infile.fail()){
		return false;
	}
	while(!infile.eof()) {
		getline(infile,line); 		
		if(line.length()>0){
			if(line[0]=='#' && line.length()>1){
				currentParameter.assign(line.substr(1));
				trim(currentParameter);				
			}
			else{
				if(currentParameter.length()==0){
					// Single target secondary structure
					if(validSecondaryStructure(line,true)>=0){
						strs->assign(line);
					}
					// Multiple target secondary structures
					else if(line.find(STRUCTURE_DELIMITER) != std::string::npos){
						bool isStr=true;
						std::vector<std::string> decomposedLine;
						split(line,STRUCTURE_DELIMITER, decomposedLine);
						for(int i=0; i<decomposedLine.size();i++){
							if(validSecondaryStructure(decomposedLine[i],true)<0){
								isStr=false;
							}
						}
						if(isStr){
							strs->assign(line);
						}						
					}
					// Sequence constraints
					else{
						bool isSeq = true;
						for(int i=0; i<line.length();i++){
							if(std::find(validIUPAC.begin(), validIUPAC.end(), line[i]) == validIUPAC.end()) {
								isSeq = false;
								break;
							}
						}
						if(isSeq){
							seqConst->assign(line);
						}
					}
				}
				else{
					trim(line);
					if(currentParameter=="RNAscdstr"){
						strs->assign(line);
					}
					else if(currentParameter=="RNAseqcon"){
						seqConst->assign(line);
					}
					else if(currentParameter=="MAXsol"){
						*maxSolutions=atoi(line.c_str());
					}
					else if(currentParameter=="TimeLimit"){
						*timeLimit=atoi(line.c_str());
					}
					else if(currentParameter=="temp"){
						strTemps->assign(line);
					}
					else if(currentParameter=="dangles"){
						*dangles=atoi(line.c_str());
					}
					else if(currentParameter=="RNAlibrary"){
						rnaLib->assign(line);
					}
					else if(currentParameter=="EnergyModel"){
						energyModel->assign(line);
					}
					else if(currentParameter=="AAseqcon"){
						aaConst->assign(line);
					}
					else if(currentParameter=="AAsimilCstr"){
						*aaSimilCstr=atoi(line.c_str());
					}
					else if(currentParameter=="AAtarget"){
						aaTarget->assign(line);
					}					
					else if(currentParameter=="AAstartPos"){
						aaStartPos->assign(line);
					}
					else if(currentParameter=="MaxBlosumScore"){
						*maxBlosumScore=atoi(line.c_str());
					}
					else if(currentParameter=="RNAcompstr"){
						compStr->assign(line);						
					}
					else if(currentParameter=="IncompBP"){
						incompBP->assign(line);
					}
					else if(currentParameter=="minGCcont"){
						*minGCcont=atof(line.c_str());
					}
					else if(currentParameter=="maxGCcont"){
						*maxGCcont=atof(line.c_str());
					}
					else if(currentParameter=="minAU"){
						*minAU=atoi(line.c_str());
					}
					else if(currentParameter=="maxAU"){
						*maxAU=atoi(line.c_str());
					}
					else if(currentParameter=="minGC"){
						*minGC=atoi(line.c_str());
					}
					else if(currentParameter=="maxGC"){
						*maxGC=atoi(line.c_str());
					}
					else if(currentParameter=="minGU"){
						*minGU=atoi(line.c_str());
					}
					else if(currentParameter=="maxGU"){
						*maxGU=atoi(line.c_str());
					}
					else if(currentParameter=="minA"){
						minA->assign(line);
					}
					else if(currentParameter=="maxA"){
						maxA->assign(line);
					}
					else if(currentParameter=="minC"){
						minC->assign(line);
					}
					else if(currentParameter=="maxC"){
						maxC->assign(line);
					}
					else if(currentParameter=="minG"){
						minG->assign(line);
					}
					else if(currentParameter=="maxG"){
						maxG->assign(line);
					}
					else if(currentParameter=="minU"){
						minU->assign(line);
					}
					else if(currentParameter=="maxU"){
						maxU->assign(line);
					}
					else if(currentParameter=="consA"){
						consA->assign(line);
					}
					else if(currentParameter=="consC"){
						consC->assign(line);
					}
					else if(currentParameter=="consG"){
						consG->assign(line);
					}
					else if(currentParameter=="consU"){
						consU->assign(line);
					}
					else if(currentParameter=="MFEstructure"){
						*MFEstructure=atoi(line.c_str());
					}
					else if(currentParameter=="MaxMFE"){
						*upperMFE=atof(line.c_str());
					}
					else if(currentParameter=="MinMFE"){
						*lowerMFE=atof(line.c_str());
					}
					else if(currentParameter=="MinimizeMFE"){
						*minimizeMFE=atoi(line.c_str());
					}
					else if(currentParameter=="MinimizeEnsDef"){
						*minimizeEnsDef=atoi(line.c_str());
					}
					else if(currentParameter=="ShowHelices"){
						*showHelices=atoi(line.c_str());
					}
					else if(currentParameter=="HelixCstrs"){
						helixCstrsStr->assign(line);
					}
					else if(currentParameter=="LocalCstrs"){
						localCstrsStr->assign(line);
					}
					else if(currentParameter=="VarHeuristic"){
						*varHeuristic=atoi(line.c_str());
					}
					else if(currentParameter=="HelixHeuristic"){
						*helixHeuristic=atoi(line.c_str());
					}
					else if(currentParameter=="LNS"){
						*LNS=atoi(line.c_str());
					}
					

					//int randomAssignment = FLAGS_RandomAssignment;
					//int upthreshold = FLAGS_UPthreshold;
					//int bpthreshold = FLAGS_BPthreshold;
					//int includeDangles = FLAGS_BPthreshold;
					//int LNS = FLAGS_LNS;
					//int LNSunchangedRestarts = FLAGS_LNSunchangedRestarts;
					//int LNSrestartTime = FLAGS_LNSrestartTime;
					//int LNStimeMultiplier = FLAGS_LNStimeMultiplier;

					currentParameter="";

				}
			}
		}
	}
	infile.close();
	return true;
}

