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


			v->fillBasePairs1Index(mfeStr);
			
//			for(int i = 1; i <= size(); ++i) {
//				cout << structure_[i]<< "\t";
//			}
//			cout << endl;
//			for(int i = 1; i <= size(); ++i) {
//				cout << mfeStr[i]<< "\t";
//			}
//			cout << endl;
			bool continueSearch=false;
						
//			for (int i = posLeft_; i <= posRight_; ++i) {  // CHECK ONLY POSITIONS INTERNAL HELIX 
			for (int i = 0; i < size(); ++i) {
				if (mfeStr[i] != structure_[i]) {
					if(!checkSides){
//						cout << endl << "Fail " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl;
//						string a;
//						cin >> a;*/

//						printf("position %d --> %d != %d\n", i,mfeStr[i],structure_[i]);
						solver()->Fail();
						return;						
					}
					else{
						continueSearch=true;
						break;
					}
				}
			}

			if(!continueSearch){
// 				cout << endl << "OKOK " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl << "Energy:" << energy <<endl;
//				string a;
//				cin >> a;

				if(upperMFE_< NO_ENERGY_LIMIT && (float)energy>(float)upperMFE_){
					solver()->Fail();
				}
				else if(lowerMFE_< NO_ENERGY_LIMIT && (float)energy<(float)lowerMFE_){
					solver()->Fail();
				}

				return;
			}
			else{
				string seqSides(nSides,'N');
				for(int i = 0; i < size(); ++i) {
					seqSides[i+(varLeft_==NULL? 0: 1)]=ToNucl(vars_[i]->Value());
				}
				vector<int64> leftDom; 
				vector<int64> rightDom; 
				if(varLeft_!=NULL){
					for (leftIt_->Init(); leftIt_->Ok(); leftIt_->Next()) {
						leftDom.push_back(leftIt_->Value());						
					}
				}
				if(varRight_!=NULL){
					for (rightIt_->Init(); rightIt_->Ok(); rightIt_->Next()) {
						rightDom.push_back(rightIt_->Value());						
					}
				}
				for(int l = 0; l < (varLeft_==NULL ? 1: leftDom.size()); l++) {
					if(varLeft_!=NULL){
						seqSides[0]=ToNucl(leftDom[l]);
					}
					for(int r = 0; r < (varRight_==NULL ? 1: rightDom.size()); r++) {
						if(varRight_!=NULL){
							seqSides[nSides-1]=ToNucl(rightDom[r]);
						}
						vSides->setSequence(seqSides.c_str());
						if(cutPoint_==-1){
							energy=vSides->fold();
						}
						else{
							energy=vSides->cofold();
						}

						vSides->fillBasePairs1Index(mfeStrSides);
						bool isSolution=false;
						// Check sides
						if((varLeft_==NULL || mfeStrSides[1] == -1) && (varRight_==NULL || mfeStrSides[nSides] == -1)){
							isSolution=true;
							for (int i = posLeft_; i <= posRight_; ++i) {
								if (mfeStrSides[i+(varLeft_==NULL ? 0: 1)] != structure_[i]) {
//									cout << endl << "Fail " << l << " " << r<<endl << seqSides<< endl << vSides->getStructure() << endl << (varLeft_==NULL ? "" : " ") << targetStr << (varRight_==NULL ? "" : " ") << endl;
//									string a;
//									cin >> a;*/

//									printf("position %d --> %d != %d\n", i,mfeStrSides[i+(varLeft_==NULL ? 0: 1)],structure_[i]);
									isSolution=false;
									break;
								}
							}
						}
						
						if(isSolution){
							if((upperMFE_ == NO_ENERGY_LIMIT || (float)energy<=(float)upperMFE_) && (lowerMFE_== NO_ENERGY_LIMIT || (float)energy>=(float)lowerMFE_)){
//								cout << endl << "OK " << l << " " << r<<endl << seqSides<< endl << vSides->getStructure() << endl << (varLeft_==NULL ? "" : " ") << targetStr << (varRight_==NULL ? "" : " ") << endl;
								
								return;
							}
						}
					}
				}
				solver()->Fail();
				return;
			}
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

			v->fillBasePairs1Index(mfeStr);

			bool continueSearch=false;

			for (int i = posLeft_; i <= posRight_; ++i) {
				if (structure_[i] != -2 && mfeStr[i] != structure_[i]) {
//					cout << endl << "Fail " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl;
//					string a;
//					cin >> a;

//					printf("position %d --> %d != %d\n", i,mfeStr[i],structure_[i]);

					if(!checkSides){
						solver()->Fail();
						return;						
					}
					else{
						continueSearch=true;
						break;
					}
				}
			}

			if(!continueSearch){
// 				cout << endl << "OKOK " << endl << seq<< endl << v->getStructure() << endl << targetStr << endl << "Energy:" << energy <<endl;
//				string a;
//				cin >> a;
				if(upperMFE_< NO_ENERGY_LIMIT && (float)energy>(float)upperMFE_){
						solver()->Fail();
				}
				else if(lowerMFE_<NO_ENERGY_LIMIT && (float)energy<(float)lowerMFE_){
						solver()->Fail();
				}

				return;
			}
			else{
				string seqSides(nSides,'N');
				for(int i = 0; i < size(); ++i) {
					seqSides[i+(varLeft_==NULL? 0: 1)]=ToNucl(vars_[i]->Value());
				}
				vector<int64> leftDom; 
				vector<int64> rightDom; 
				if(varLeft_!=NULL){
					for (leftIt_->Init(); leftIt_->Ok(); leftIt_->Next()) {
						leftDom.push_back(leftIt_->Value());						
					}
				}
				if(varRight_!=NULL){
					for (rightIt_->Init(); rightIt_->Ok(); rightIt_->Next()) {
						rightDom.push_back(rightIt_->Value());						
					}
				}
				for(int l = 0; l < (varLeft_==NULL ? 1: leftDom.size()); l++) {
					if(varLeft_!=NULL){
						seqSides[0]=ToNucl(leftDom[l]);
					}
					for(int r = 0; r < (varRight_==NULL ? 1: rightDom.size()); r++) {
						if(varRight_!=NULL){
							seqSides[nSides-1]=ToNucl(rightDom[r]);
						}
						vSides->setSequence(seqSides.c_str());
						if(cutPoint_==-1){
							energy=vSides->fold();
						}
						else{
							energy=vSides->cofold();
						}

						vSides->fillBasePairs1Index(mfeStrSides);
						bool isSolution=false;
						// Check sides
						if((varLeft_==NULL || mfeStrSides[1] == -1) && (varRight_==NULL || mfeStrSides[nSides] == -1)){
							isSolution=true;
							for (int i = posLeft_; i <= posRight_; ++i) {
								if (structure_[i] != -2 && mfeStrSides[i+(varLeft_==NULL ? 0: 1)] != structure_[i]) {
				//					cout << endl << "Fail " << l << " " << r<<endl << seqSides<< endl << vSides->getStructure() << endl << (varLeft_==NULL ? "" : " ") << targetStr << (varRight_==NULL ? "" : " ") << endl;
				//					string a;
				//					cin >> a;*/

				//					printf("position %d --> %d != %d\n", i,mfeStrSides[i+(varLeft_==NULL ? 0: 1)],structure_[i]);
									isSolution=false;
									break;
								}
							}
						}
						
						if(isSolution){
							if((upperMFE_ == NO_ENERGY_LIMIT || (float)energy<=(float)upperMFE_) && (lowerMFE_== NO_ENERGY_LIMIT || (float)energy>=(float)lowerMFE_)){
								return;
							}
						}
					}
				}
				solver()->Fail();
				return;
			}
			
		}
	}


}
