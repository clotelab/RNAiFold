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

#include <string.h>
#include <vector>

#include "helix_constraint.h"
#include "local_constraint.h"
#include "aaconstraint.h"

// OR-Tools libraries 
#include "base/logging.h"
#include "constraint_solver/constraint_solver.h"

#define RNAIFOLD_VERSION "3.0"
#define RNAIFOLD_NAME "RNAiFold"

#define DEFAULT_MAX_SOLUTIONS 8
#define DEFAULT_TIME_LIMIT 600

#define DEFAULT_BP_THRESHOLD 85
#define DEFAULT_UP_THRESHOLD 100

#define DEFAULT_TEMP "37"
#define DEFAULT_DANGLES 2
#define DEFAULT_ENERGY_MODEL TURNER_04_CODE
#define DEFAULT_RNA_LIBRARY VIENNA_LIB

#define DEFAULT_LNS_UNCHANGED_RESTARTS 10
#define DEFAULT_LNS_RESTART_TIME 15
#define DEFAULT_LNS_RESTART_TIME_MULTIPLIER 0 // 50

#define MAXSIZE 3000

namespace operations_research {
	void IfoldCp(std::vector<int*> strs_int, std::vector<int*> int_strs_undet, int n, int maxSolutions, int64 timeLimit, char* sequence, std::vector<AAConstraint*> aaConstraints, int helixHeuristic, int varHeuristic, int randomAssignment, int upthreshold, int bpthreshold, int includeDangles, int dangles, std::string rnaLib, std::string energyModel, std::vector<double> foldTemps, double minGCcont, double maxGCcont, int minAU, int maxAU, int minGC, int maxGC, int minGU, int maxGU,std::vector<std::tuple<int,int,int>> listMinA,std::vector<std::tuple<int,int,int>> listMaxA,std::vector<std::tuple<int,int,int>> listMinC,std::vector<std::tuple<int,int,int>> listMaxC,std::vector<std::tuple<int,int,int>> listMinG,std::vector<std::tuple<int,int,int>> listMaxG,std::vector<std::tuple<int,int,int>> listMinU,std::vector<std::tuple<int,int,int>> listMaxU,std::vector<std::tuple<int,int,int>> listConsA,std::vector<std::tuple<int,int,int>> listConsC,std::vector<std::tuple<int,int,int>> listConsG,std::vector<std::tuple<int,int,int>> listConsU, int MFEstructure, int minimizeMFE, int minimizeEnsDef, int* comp_str_int, vector<pair<int,int>> vIncompBP, int showHelices, std::vector<HelixCstr> helixCstrs, std::vector<LocalCstr> localCstrs, int LNS, int LNSunchangedRestarts, int LNSrestartTime, int showMeasures, int cutPoint);
}

bool checkInputParameters(std::vector<std::string> structures, int maxSolutions, int64 timeLimit, std::string seqConst,std::vector<AAConstraint*> aaConstraints,std::string compStr, std::string incompBP, int helixHeuristic, int varHeuristic, int randomAssignment, int upthreshold, int bpthreshold, int includeDangles , int dangles, std::string rnaLib, std::string energyModel, std::vector<double> foldTemps, double minGCcont, double maxGCcont, int minAU,int maxAU,int minGC,int maxGC,int minGU,int maxGU,int MFEstructure,double upperMFE, double lowerMFE, int minimizeMFE, int minimizeEnsDef,int showHelices, std::vector<HelixCstr> helixCstrs, std::vector<LocalCstr> localCstrs, int LNS, int LNSunchangedRestarts, int LNSrestartTime, int LNStimeMultiplier);
int findAndRemoveCutPoints(std::vector<std::string>* structures,std::string* seqConst,std::string* compStr,std::string* incompBP);


bool readDatafromFile(std::string inputFile, std::string *strs, int* maxSolutions, int64* timeLimit, std::string* seqConst,std::string* aaConst, int* aaSimilCstr, std::string* aaTarget, std::string* aaStartPos, int* maxBlosumScore, std::string* compStr, std::string* incompBP, int* helixHeuristic, int* varHeuristic, int* randomAssignment, int* upthreshold, int* bpthreshold, int *includeDangles , int* dangles, std::string* rnaLib, std::string* energyModel, std::string* strTemps, double* minGCcont, double* maxGCcont, int* minAU,int* maxAU,int* minGC,int* maxGC,int* minGU,int* maxGU,std::string* minA,std::string* maxA,std::string* minC,std::string* maxC,std::string* minG,std::string* maxG,std::string* minU,std::string* maxU,std::string* consA,std::string* consC,std::string* consG,std::string* consU,int* MFEstructure, double* upperMFE, double* lowerMFE, int* minimizeMFE, int* minimizeEnsDef, int* showHelices, std::string* helixCstrsStr, std::string* localCstrsStr, int* LNS, int* LNSunchangedRestarts, int* LNSrestartTime, int* LNStimeMultiplier);
