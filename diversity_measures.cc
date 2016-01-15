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
#include <math.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include "diversity_measures.h"

using namespace std;
namespace operations_research {
	string getNtContent(int* str, vector<IntVar *> seq, int n){
		std::ostringstream oss;
		int AUnum = 0;
		int GCnum = 0;
		int GUnum = 0;
		double GplusC =0;
		double GCcont =0;
		for (int i=1;i<=n;i++){
			if(seq[i]->Value()==1 || seq[i]->Value()==7){
				GplusC++;
			}		
			if(str[i]>i){			
				switch(abs(seq[i]->Value()-seq[str[i]]->Value())){
					case 6: GCnum ++;
							break;
					case 9: AUnum ++;
							break;				
					case 11: GUnum ++;
							 break;				
					default: break;					
				}
			}
		}
		GCcont = GplusC/n;

		oss << "GC content: " <<std::setprecision(2)<< GCcont << " - AUs: " << AUnum << " - GCs: " << GCnum << " - GUs: " << GUnum;

		return oss.str();
	}
}
double expectedPointwiseEntropy(double ** bppr, int n)
{
	double h=0;
	double h_rna = 0;
	for (int i=0;i<n;i++)
		for (int j=0;j<n+1;j++)
		{
			if(bppr[i][j]>0)
				h += bppr[i][j]*log(bppr[i][j]);
		}
		h_rna = -h/n;
	return h_rna;
}

double morganHiggsStructuralDiversity (double **bppr, int n)
{
	double sum=0;
	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++){
			if(i!=j){
				sum += pow(bppr[i][j],2);
			}
		}
	for (int i=0;i<n;i++)
		sum += pow(bppr[i][n],2);
	return n-sum; 	
}

double viennaStructuralDiversity (double **bppr, int n)
{
	double sum = 0;
	for (int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			sum += bppr[i][j]*(1-bppr[i][j]);
	return sum;
}

double expectedBpDistance(double **bppr, int n, int* bpList)
{
	double sum = 0;
	for (int i=0;i<n;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			if(bpList[i]!=j )
			  sum += bppr[i][j];
			else 
			  sum += 1 - bppr[i][j];
		}
	}

	return sum;
}

double ensembleDefect(double **bppr, int n , int* bpList)
{
	double sum = 0;
	for(int i=0;i<n;i++)
	{
		if (bpList[i] == -1)
			sum += bppr[i][n];
		else
			sum += bppr[i][bpList[i]];
	}
	return n-sum;
}

int * getBasePairs(char * str, int n)
{
	vector <int> stack;
	int * bp = new int[n];
	for(int i=0; i<n;i++)
	{
		if (str[i] == '(' )
			stack.push_back(i);
		else if (str[i] == ')')
		{
			int m = stack.back();
			bp[i] = m;
			bp[m] = i;
			stack.pop_back();
		}
		else
			bp[i] = -1;
	}
		
	return bp;
}

