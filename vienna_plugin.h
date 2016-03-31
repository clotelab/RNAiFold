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
 #ifndef _VIENNA_PLUGIN
 #define _VIENNA_PLUGIN
 #include "rna_plugin.h"
 
 class ViennaPlugin: public RNAPlugin{
 	public:
 		ViennaPlugin(){};
 		ViennaPlugin(int size);
 		ViennaPlugin(int size,int dangles);
 		~ViennaPlugin();
 	
 		void initMe(std::string name, int d);
 		void setSequence(std::string seq);
 		void setStructure(std::string str);
 		void setTestStructure(std::string str);
 		double fold();
 		string getStructure();
 		double energyOfStruct(); 
 		double energyOfEnsemble(); 
 		double probOfStruct();
 		int* getBasePairs();
 		int* getBasePairs1Index();
		void fillBasePairs1Index(int *bprs);
 		double* basePairProbsArray();  
 		double** basePairProbsMatrix();
 		double* basePairProbsArrayCofold();
 		double** basePairProbsMatrixCofold();
 		void setWindowSize(int ws);
 		void setTemperature(double t);
 		void setCutPoint(int p);
 		double cofold();
 		void freeMe();
 		void setDangles(int d);
 		void setEnergyModel(string energy_model);
 	private:
 		string _name;
 		string _seq;
 		string _str;
 		string _tstr;
 		char *_cseq;
 		char *_cstr;
 		char *_ctstr;
 		int _winSize;
 		double _temp;
 		int _cutpoint;
 		int _d;

 };
 
 #endif
