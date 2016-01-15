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

#include <vector>
#include <string>
#include <map>
#include <string>
#include "constraint_solver/constraint_solver.h"

#ifndef _LOCAL_CSTR
#define _LOCAL_CSTR

#define LC_WRONG_STRUCTURE -2
#define LC_WRONG_TYPE -1
#define LC_ENERGY_MIN 1
#define LC_ENERGY_MAX 2
#define LC_ED_MIN 3
#define LC_ED_MAX 4
#define LC_MFE 5

#define LC_ENERGY_MIN_STR "MinEnergy"
#define LC_ENERGY_MAX_STR "MaxEnergy" 
#define LC_ED_MIN_STR "MinED"
#define LC_ED_MAX_STR "MaxED"
#define LC_MFE_STR "MFE"

class LocalCstr{
	public:
		~LocalCstr();
		LocalCstr(std::string structure, int limitType, double limitValue);
		LocalCstr(std::string structure, const std::string limitType, double limitValue);
		LocalCstr(int start_pos, std::string structure, int limitType, double limitValue);
		LocalCstr(int start_pos, std::string structure, const std::string limitType, double limitValue);
		int getStartPos();
		void setStartPos(int id);
		std::string getStructureStr();
		void setStructureStr(std::string structure);
		std::vector<int64> getStructureArr();
		void setStructureArr(std::vector<int64> int_str);
		int getLimitType();
		void setLimitType(int limitType);
		double getLimitValue();
		void setLimitValue(double limitValue);
		bool hasUndet();
		bool isPairedConstraint(LocalCstr c);
		std::string printTypes();
	private:
		int start_pos_;
		std::string structure_;
		std::vector<int64> int_str_;
		int limitType_;
		double limitValue_;
		std::map<const std::string, int> limitTypes_ = {{LC_ENERGY_MIN_STR,LC_ENERGY_MIN},{LC_ENERGY_MAX_STR,LC_ENERGY_MAX},{LC_ED_MIN_STR,LC_ED_MIN},{LC_ED_MAX_STR,LC_ED_MAX},{LC_MFE_STR,LC_MFE}};
		bool has_undet_;
		
};
#endif

