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
#ifndef _ENERGY_CONSTANT

#define SIGMA 3

#define NO_ENERGY_LIMIT 1000
#define ENERGY_PRECISION 100

#define NO_ED_LIMIT -1
#define ED_PRECISION 10000


#define VIENNA_LIB "Vienna"
#define RNASTRUCTURE_LIB "RNAstructure"

#define VIENNA_LIB "Vienna"
#define RNASTRUCTURE_LIB "RNAstructure"

#define TURNER_04_CODE "2004"
#define TURNER_99_CODE "1999"
#define ANDRONESCU_07_CODE "2007"


#define TURNER_04_FILE "ViennaRNA/energy_par_files/rna_turner2004.par"  // Energy parameters file for Turner '04 model 
#define TURNER_99_FILE "ViennaRNA/energy_par_files/rna_turner1999.par"  // Energy parameters file for Turner '99 model 
#define ANDRONESCU_07_FILE "ViennaRNA/energy_par_files/rna_andronescu2007.par"  // Energy parameters file for Andronescu '07 model 
#define RNASTRUCTURE_DIR "RNAstructure/data_tables"    // Energy parameter files directory for RNAstructure model 

#define HH_OVERLAP_SIMPLE 1
#define HH_OVERLAP_BP 2
#define HH_OVERLAP_POSITIONS 3
#define HH_OVERLAP_BP_PERCENT 4

#define SH_NONE 0
#define SH_IN_TO_OUT 1
#define SH_BOTTOM_TO_TOP 2
#define SH_BOTTOM_TO_TOP_UP 3
#define SH_BOTTOM_TO_TOP_BAK 4
#define SH_MULTI_STR 5


#endif
