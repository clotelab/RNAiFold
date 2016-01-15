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

#ifndef _HELIX_CSTR
#define _HELIX_CSTR

#define HC_WRONG_TYPE -1
#define HC_ENERGY_MIN 1
#define HC_ENERGY_MAX 2
#define HC_ED_MIN 3
#define HC_ED_MAX 4

#define HC_ENERGY_MIN_STR "MFEmin"
#define HC_ENERGY_MAX_STR "MFEmax" 
#define HC_ED_MIN_STR "EDmin"
#define HC_ED_MAX_STR "EDmax"

class HelixCstr{
	public:
		~HelixCstr();
		HelixCstr(int helix_id, int limitType, double limitValue);
		HelixCstr(int tree_id, int helix_id, int limitType, double limitValue);
		HelixCstr(int helix_id, const std::string limitType, double limitValue);
		HelixCstr(int tree_id, int helix_id, const std::string limitType, double limitValue);
		int getTreeId();
		void setTreeId(int id);
		int getHelixId();
		void setHelixId(int id);
		int getLimitType();
		void setLimitType(int limitType);
		double getLimitValue();
		void setLimitValue(double limitValue);
		std::string printTypes();
	private:
		int tree_id_;
		int helix_id_;
		int limitType_;
		double limitValue_;
		std::map<const std::string, int> limitTypes_ = {{HC_ENERGY_MIN_STR,HC_ENERGY_MIN},{HC_ENERGY_MAX_STR,HC_ENERGY_MAX},{HC_ED_MIN_STR,HC_ED_MIN},{HC_ED_MAX_STR,HC_ED_MAX}};
		
};
#endif

