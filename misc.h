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

#include <string>
#include <vector>
#include <tuple>
#ifndef _MISC_LIB
#define _MISC_LIB
/* declaration of functions for initializaing, transformation */

#define STRUCTURE_DELIMITER '|'
#define FIELD_DELIMITER ','

extern int d2[13];
extern int simple[13];
extern int d1[29];
extern int cbp[29];
extern int guau[29];
extern double R;
extern double T;

int* make_BasePair_Table(std::string structure,bool withUndet);
char ToNucl(int trans);
int Nt2Int(char nucl);
	
int * IntSeq(char * sequence);
int * TransformSeq(char* sequence); // converts sequence to list of int through dictionary G->1, A->3, C->7, U->12
int GetIndex(int i, int j); // returns value in  d1 through index from formula (i-j)+12+(i==j)[d2[i]]
void PrintStr(char *ptr); // prints string given ptr, end character must be \0
void PrintInt(int *ptr, int len);
int CountChar(char * str, char c);// counts the number of char c in str
int **Allocate2DMatrix(int a, int b);// Allocate 2D Matrix of axb
int **Allocate2Dmatrix(int a, int b);
double **Allocate2DMatrixDouble(int a, int b);// Allocate 2D double Matrix of axb
int **BPBetweenBool(int** Matrix, int * bp, int len);//fills in matrix such that matrix[i][j] returns 1 if BP between i,j, otherwise 0
int * GetBPList(char * str,int len);
char GetNucl(int trans);
int linearize2D(int row, int column);
bool is_number(const std::string& s);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
int parseTriple(std::string strTuple, int maxValue, std::string name, std::vector<std::tuple<int,int,int>>* tTuple);
int validSecondaryStructure(std::string secStr, bool acceptCommas);
int findCutPoint(std::string* secStr);
std::string getExecPath(std::string argv0);

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
	// initialize original index locations
	std::vector<size_t> idx(v.size());
	for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
	   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}
#endif
