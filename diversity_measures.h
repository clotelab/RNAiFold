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
#include <string.h>
#include "constraint_solver/constraint_solver.h"

using namespace std; 

namespace operations_research {
	std::string getNtContent(int* str, vector<IntVar *> seq, int n);
}
double getBoltzmanProbability(double );
double expectedPointwiseEntropy(double ** , int);
double morganHiggsStructuralDiversity (double ** , int);
double viennaStructuralDiversity (double ** , int);
double expectedBpDistance(double ** , int , int* bpList);
double ensembleDefect(double ** , int n , int* bpList);


