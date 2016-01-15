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

#include "vienna_constraint.h"
#include "misc.h"
#include <cstring>

namespace operations_research {
	void ViennaConstraintDet::InitialPropagate() {
	    int allBound =1;
//	    printf("Initial propagate\n");
		vector<int64> freeVars;
		for (int i = 0; i < size(); ++i) {
			if (!vars_[i]->Bound()) {
				allBound=0;
//				for(int i = 0; i < size(); ++i) {
//					if (!vars_[i]->Bound()) {
//						cout <<"N";
//					}
//					else{
//						cout << ToNucl(vars_[i]->Value());
//					}
//				}
//				cout << endl << targetStr <<endl;
				
				return;
			}
		}
		if(allBound){
			int* mfeStr;
			char seq[n+1];
			double energy=0;

			for(int i = 0; i < size(); ++i) {
				seq[i]=ToNucl(vars_[i]->Value());
			}
			seq[n]='\0';


			v->setSequence(seq);
			if(cutPoint_==-1){
				energy=v->fold();
			}
			else{
				energy=v->cofold();
			}
			
//			string a;
//			cout << "**" << endl;
//			cout << seq<<endl;
//			cout << "**" << endl;
//			cin >> a;


			mfeStr = v->getBasePairs1Index();

//			for(int i = 1; i <= size(); ++i) {
//				cout << structure_[i]<< "\t";
//			}
//			cout << endl;
//			for(int i = 1; i <= size(); ++i) {
//				cout << mfeStr[i]<< "\t";
//			}
//			cout << endl;
			
			for (int i = posLeft_; i <= posRight_; ++i) {
				if (mfeStr[i] != structure_[i]) {
//					cout << endl << "Fail " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl;
//					string a;
//					cin >> a;*/

//					printf("position %d --> %d != %d\n", i,mfeStr[i],structure_[i]);
					free(mfeStr);
					solver()->Fail();
					return;
				}
			}

// 			cout << endl << "OKOK " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl << "Energy:" << energy <<endl;
//			string a;
//			cin >> a;
			free(mfeStr);

			if(upperMFE_< NO_ENERGY_LIMIT && (float)energy>(float)upperMFE_){
				solver()->Fail();
			}
			else if(lowerMFE_< NO_ENERGY_LIMIT && (float)energy<(float)lowerMFE_){
				solver()->Fail();
			}

			return;
		}
	}

	void ViennaConstraintUndet::InitialPropagate() {
	    int allBound =1;
//	    printf("Initial propagate\n");
		vector<int64> freeVars;
		for (int i = 0; i < size(); ++i) {
			if (!vars_[i]->Bound()) {
				allBound=0;
//				for(int i = 0; i < size(); ++i) {
//					if (!vars_[i]->Bound()) {
//						cout <<"N";
//					}
//					else{
//						cout << ToNucl(vars_[i]->Value());
//					}
//				}
//				cout << endl << targetStr <<endl;
				
				return;
			}
		}
		if(allBound){
			int* mfeStr;
			char seq[n+1];
			double energy=0;

			for(int i = 0; i < size(); ++i) {
				seq[i]=ToNucl(vars_[i]->Value());
			}
			seq[n]='\0';


			v->setSequence(seq);
			if(cutPoint_==-1){
				energy=v->fold();
			}
			else{
				energy=v->cofold();
			}
			
//			string a;
//			cout << "**" << endl;
//			cout << seq<<endl;
//			cout << "**" << endl;
//			cin >> a;

//			for(int i = 1; i <= size(); ++i) {
//				cout << structure_[i]<< "\t";
//			}
//			cout << endl;
//			for(int i = 1; i <= size(); ++i) {
//				cout << mfeStr[i]<< "\t";
//			}
//			cout << endl;

			mfeStr = v->getBasePairs1Index();

			for (int i = posLeft_; i <= posRight_; ++i) {
				if (structure_[i] != -2 && mfeStr[i] != structure_[i]) {
//					cout << endl << "Fail " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl;
//					string a;
//					cin >> a;

//					printf("position %d --> %d != %d\n", i,mfeStr[i],structure_[i]);

					free(mfeStr);
					solver()->Fail();
					return;
				}
			}

// 			cout << endl << "OKOK " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl << "Energy:" << energy <<endl;
//			string a;
//			cin >> a;
			free(mfeStr);

			if(upperMFE_< NO_ENERGY_LIMIT && (float)energy>(float)upperMFE_){
					solver()->Fail();
			}
			else if(lowerMFE_<NO_ENERGY_LIMIT && (float)energy<(float)lowerMFE_){
					solver()->Fail();
			}

			return;
		}
	}


}
