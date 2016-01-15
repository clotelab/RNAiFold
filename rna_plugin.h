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

using namespace std;
#ifndef _RNA_PLUGIN
#define _RNA_PLUGIN

#include "energy_constant.h"

class RNAPlugin{
	public:
		virtual void initMe(std::string name, int d) = 0;
		virtual void setSequence(std::string seq) = 0;
		virtual void setStructure(std::string str) = 0;
		virtual void setTestStructure(std::string str) = 0;
		virtual double fold() = 0;
		virtual string getStructure() = 0;
		virtual double energyOfStruct() = 0; 
		virtual double energyOfEnsemble() = 0; 
		virtual double probOfStruct() = 0;
		virtual int* getBasePairs() = 0;
		virtual int* getBasePairs1Index() = 0;
		virtual double* basePairProbsArray() = 0;  
		virtual double** basePairProbsMatrix() = 0;
		virtual double* basePairProbsArrayCofold() = 0;
		virtual double** basePairProbsMatrixCofold() = 0;
		virtual void setWindowSize(int ws) = 0;
		virtual void setTemperature(double t) = 0;
		virtual void setCutPoint(int p) = 0;
		virtual double cofold() = 0;
		virtual void freeMe() = 0;
		virtual void setDangles(int d) = 0;
		virtual void setEnergyModel(string energy_model) = 0;

		virtual ~RNAPlugin(){};
};

RNAPlugin* newRNAplugin(std::string rnaLib, int size, int dangles);
#endif
