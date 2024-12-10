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

#include "base/logging.h"
#include "util/string_array.h"
#include "constraint_solver/constraint_solver.h"

#include <ctime>
#include "strtreegroup.h"
#include "value_heuristic.h"
#include "vienna_plugin.h"
#define NOT_ASSIGNED -2

using namespace std;
namespace operations_research {
	class LNSRestart : public SearchMonitor {
	public:
		LNSRestart(Solver* s, BaseAssignVariables* myDB_, const std::vector<IntVar*>& vSeq, vector<int*> str_int, int dangles, std::string rnaLib, std::string energyModel, vector<double> temperatures, double timeLimit, int maxRestarts, StrTreeGroup* treeGroup)
		: SearchMonitor(s), myDB(myDB_), vSeq_(vSeq), str_int_(str_int), dangles_(dangles), timeLimit_(timeLimit), maxRestarts_(maxRestarts), treeGroup_(treeGroup), nRestarts(0) { //, current_fails_(0) {
			CHECK_GE(maxRestarts, 1);
			CHECK_GE(timeLimit, 1);			
			fixedPositions.assign(size(), NOT_ASSIGNED);
			lastAssignment.assign(size(), true);
//			ASSIGNMENT BY VALUES
//			lastAssignment.assign(size(), NOT_ASSIGNED);
			time(&timeLastFail);
			for(int i=0; i< str_int.size(); i++){
				RNAPlugin* newViennaPlugin = newRNAplugin(rnaLib,size()-1,dangles_);
				newViennaPlugin->setEnergyModel(energyModel);
				newViennaPlugin->setTemperature(temperatures[i]);
				viennaPlugin.push_back(newViennaPlugin);
			}
			threshold_= initThreshold;
		}

		virtual ~LNSRestart() {
			for(int i=0; i< viennaPlugin.size();i++){
				delete(viennaPlugin.at(i));
			}
		}

		virtual void BeginFail() {
			double diffSeconds;
			time_t currentTime;
			time(&currentTime);
			diffSeconds=difftime(currentTime,timeLastFail);
			if(diffSeconds>timeLimit_){
				//cout << "More than "<< timeLimit_ << " seconds since last restart: " << diffSeconds<< " - Restarts: "<< nRestarts << endl;
				time(&timeLastFail);
				if(nRestarts<maxRestarts_){
					myDB->clearFixedPositions();
					FixPositions();
//					cout << "Threshold" << threshold_ << endl;
					for (int i = 1; i < size(); ++i) {
						if(fixedPositions[i]!=NOT_ASSIGNED){
							double randVal = rand()%100;
							if(randVal<threshold_){							
								myDB->addFixedPosition(vSeq_[i],fixedPositions[i]);
							}
						}
					}
					threshold_-=thresholdDecrease;
				}
				else{
//					cout << " Full restart: threshold = " << threshold_ << endl;
					threshold_= initThreshold;
					myDB->clearFixedPositions();
					nRestarts=1;
				}

				nRestarts++;
				solver()->RestartCurrentSearch(); 
			}
		}

		virtual std::string DebugString() const {
			return absl::StrFormat("LNSRestart(%f %d)", timeLimit_, maxRestarts_);
		}
/*		virtual void EnterSearch (){
			printf("ENTERING SEARCH\n");
		} */
		virtual void RestartSearch (){
//			printf("RESTARTING SEARCH\n");
//			for (int i = 1; i < size(); ++i) {
//				printf("%d ",fixedPositions[i]);
//			}
//			printf("\n");

		}
		void FixPositions();

		void FullRestart(){
			myDB->clearFixedPositions();
			for (int i = 1; i < size(); ++i) {
				fixedPositions[i]=NOT_ASSIGNED;
			}
			nRestarts=0;
			threshold_=initThreshold;
			nRestarts++;
			solver()->RestartCurrentSearch(); 
		}
				
		int GetNumRestarts(){
			return nRestarts;
		}
	private:
		const double initThreshold=90;
		const double thresholdDecrease=10;
		BaseAssignVariables* myDB;
		const std::vector<IntVar*> vSeq_;
		vector<int*> str_int_;

		vector<RNAPlugin*> viennaPlugin;
		int dangles_;
		double threshold_;
		const double timeLimit_;
		const int maxRestarts_;
		StrTreeGroup* treeGroup_;
		
		int nRestarts;

		int64 size() const { return vSeq_.size(); };
		std::vector<int64> fixedPositions;
		std::vector<bool> lastAssignment;
//		ASSIGNMENT BY VALUES
//		std::vector<int64> lastAssignment;

		time_t timeLastFail;
		
	};
	
	
	LNSRestart* MakeLNSRestart(Solver *s,  BaseAssignVariables* myDB, const std::vector<IntVar*>& vSeq, vector<int*> str_int, int dangles, std::string rnaLib, std::string energyModel, vector<double> temperatures, double timeLimit, int maxRestarts, StrTreeGroup* treeGroup);
	
	
	
	class HelixLNSRestart : public SearchMonitor {
	public:
		HelixLNSRestart(Solver* s, BaseAssignVariables* myDB_, const std::vector<IntVar*>& vSeq, vector<int*> str_int, int dangles, std::string rnaLib, std::string energyModel, vector<double> temperatures, double timeLimit, int maxRestarts, StrTreeGroup* treeGroup)
		: SearchMonitor(s), myDB(myDB_), vSeq_(vSeq), str_int_(str_int), dangles_(dangles), timeLimit_(timeLimit), maxRestarts_(maxRestarts), treeGroup_(treeGroup), nRestarts(0) { //, current_fails_(0) {
			CHECK_GE(maxRestarts, 1);
			CHECK_GE(timeLimit, 1);			
			fixedPositions.assign(size(), NOT_ASSIGNED);
			lastAssignment.assign(size(), true);
			time(&timeLastFail);
			nHelices = treeGroup_->getNumHelices();
			threshold_= initThreshold;
			for(int i = 0; i < nHelices; i++){
				RNAPlugin* newViennaPlugin = newRNAplugin(rnaLib,treeGroup_->getHelix(i)->getJ() - treeGroup_->getHelix(i)->getI()+1,dangles_);
				newViennaPlugin->setEnergyModel(energyModel);
				newViennaPlugin->setTemperature(temperatures[treeGroup_->getStrId(i)]);
				viennaPlugin.push_back(newViennaPlugin);

				
			}
		}

		virtual ~HelixLNSRestart() {
			for(int i=0; i< viennaPlugin.size();i++){
				delete(viennaPlugin.at(i));
			}
		}

		virtual void BeginFail() {
			double diffSeconds;
			time_t currentTime;
			time(&currentTime);
			diffSeconds=difftime(currentTime,timeLastFail);
			if(diffSeconds>timeLimit_){
//				cout << "More than "<< timeLimit_ << " seconds since last restart: " << diffSeconds<< " - Restarts: "<< nRestarts << endl;
				time(&timeLastFail);
				if(nRestarts<maxRestarts_){
					myDB->clearFixedPositions();
					FixPositions();
//					cout << "Threshold" << threshold_ << endl;
					for (int i = 1; i < size(); ++i) {
						if(fixedPositions[i]!=NOT_ASSIGNED){
							double randVal = rand()%100;
							if(randVal<threshold_){							
								myDB->addFixedPosition(vSeq_[i],fixedPositions[i]);
							}
						}
					}
					threshold_-=thresholdDecrease;
				}
				else{
//					cout << " Full restart: threshold = " << threshold_ << endl;
					threshold_= initThreshold;
					myDB->clearFixedPositions();
					nRestarts=0;
				}

				nRestarts++;
				solver()->RestartCurrentSearch(); 
			}
		}

		virtual std::string DebugString() const {
			return absl::StrFormat("HelixLNSRestart(%f %d)", timeLimit_, maxRestarts_);
		}
/*		virtual void EnterSearch (){
			printf("ENTERING SEARCH\n");
		} */
		virtual void RestartSearch (){
//			printf("RESTARTING SEARCH\n");
//			for (int i = 1; i < size(); ++i) {
//				printf("%d ",fixedPositions[i]);
//			}
//			printf("\n");

		}
		void FixPositions();

		void FullRestart(){
			myDB->clearFixedPositions();
			for (int i = 1; i < size(); ++i) {
				fixedPositions[i]=NOT_ASSIGNED;
			}
			nRestarts=0;
			threshold_=initThreshold;
			nRestarts++;
			solver()->RestartCurrentSearch(); 
		}
				
		int GetNumRestarts(){
			return nRestarts;
		}
	private:
		const double initThreshold=90;
		const double thresholdDecrease=5;
		BaseAssignVariables* myDB;
		const std::vector<IntVar*> vSeq_;
		vector<int*> str_int_;

		vector<RNAPlugin*> viennaPlugin;
		int dangles_;
		double threshold_;
		const double timeLimit_;
		const int maxRestarts_;
		StrTreeGroup* treeGroup_;
		int nHelices;
		int nRestarts;

		int64 size() const { return vSeq_.size(); };
		std::vector<int64> fixedPositions;
		std::vector<bool> lastAssignment;

		time_t timeLastFail;
		
	};

	HelixLNSRestart* MakeHelixLNSRestart(Solver *s,  BaseAssignVariables* myDB, const std::vector<IntVar*>& vSeq, vector<int*> str_int, int dangles, std::string rnaLib, std::string energyModel, vector<double> temperatures, double timeLimit, int maxRestarts, StrTreeGroup* treeGroup);
	
}  // namespace

