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
#include <map>
#include <string>
#include "energy_constant.h"
#include "blosum.h"

#ifndef _AA_CONSTRAINT
#define _AA_CONSTRAINT


#define AA_MIN_BLOSUM -4
#define AA_SIMIL_NONE 5
#define AA_SIMIL_CLASS1 6
#define AA_SIMIL_CLASS2 7


#define AA_ERR_OK 0
#define AA_ERR_LENGTH_OVERFLOW 1
#define AA_ERR_DIFFERENT_LENGTH 2
#define AA_ERR_WRONG_TARGET 3
#define AA_ERR_WRONG_SEQUENCE 4
#define AA_ERR_WRONG_BRACKETS 5
#define AA_ERR_STARTPOS 6
#define AA_ERR_SIMILARITY 7
#define AA_ERR_MAXIMIZE 8

const std::vector<char> validAaIUPAC =    {'A','E','C','D','F','G','H','I','K','M','P','L','N','Q','R','S','T','V','W','Y','*','U','X'};
const std::vector<char> validAaExtIUPAC = {'A','B','E','C','D','F','G','H','I','K','M','P','L','N','Q','R','S','T','V','W','Y','Z','*','U','X','+','-','=','&','@','^','$','%','!','~','?','[',']'};
const std::map<int,std::map<char, char > > aaClassification={{AA_SIMIL_CLASS1,{{'S', '&'}, {'T', '&'}, {'Y', '&'}, {'N', '&'}, {'Q', '&'}, {'C', '&'},                                     // Polar (S,T,Y,N,Q,C)
                                                                               {'G', '='}, {'A', '='}, {'P', '='}, {'V', '='}, {'L', '='}, {'I', '='}, {'M', '='}, {'W', '='}, {'F', '='}, // Hydrophobic (G,A,P,V,L,I,M,W,F)
                                                                               {'D', '-'}, {'E', '-'},                                                                                     // Negatively charged (D,E)
                                                                               {'H', '+'}, {'K', '+'}, {'R', '+'}                                                                          // Positively charged (H,K,R)
      	                                                                      }
      	                                                    },
      													    {AA_SIMIL_CLASS2,{{'F', '@'}, {'W', '@'}, {'Y', '@'},                                                                          // Aromatic (F,W,Y)
                                                                              {'G', '^'}, {'A', '^'}, {'V', '^'}, {'L', '^'}, {'I', '^'},                                                  // Aliphatic (G,A,V,L,I) 
                                                                              {'S', '$'}, {'C', '$'}, {'T', '$'}, {'M', '$'}, {'U', '$'},                                                  // Hydroxyl or Sulfur/Selenium-containing (S,C,T,M,U)
                                                                              {'D', '%'}, {'E', '%'}, {'N', '%'}, {'Q', '%'},                                                              // Acidic and their Amide (D,E,N,Q)
                                                                              {'H', '+'}, {'K', '+'}, {'R', '+'},                                                                          // Basic (H,K,R)
                                                                              {'P', 'P'}                                                                                                   // Cyclic == P
      	                                                                     }
      	                                                   }};
const std::map<char, std::vector<int> >  aaIUPAC = {{'A', std::vector<int> {306,181,81,31}}, 
                                                    {'E', std::vector<int> {85,35}}, 
                                                    {'C', std::vector<int> {323,198}}, 
                                                    {'D', std::vector<int> {310,185}}, 
                                                    {'F', std::vector<int> {324,199}}, 
                                                    {'G', std::vector<int> {313,188,88,38}}, 
                                                    {'H', std::vector<int> {316,191}}, 
                                                    {'I', std::vector<int> {303,178,78}}, 
                                                    {'K', std::vector<int> {90,40}}, 
                                                    {'M', std::vector<int> {28}}, 
                                                    {'P', std::vector<int> {319,194,94,44}}, 
                                                    {'L', std::vector<int> {99,49,307,182,82,32}}, 
                                                    {'N', std::vector<int> {315,190}}, 
                                                    {'Q', std::vector<int> {91,41}}, 
                                                    {'R', std::vector<int> {318,193,93,43,89,39}}, 
                                                    {'S', std::vector<int> {317,192,92,42,314,189}}, 
                                                    {'T', std::vector<int> {308,183,83,33}}, 
                                                    {'V', std::vector<int> {301,176,76,26}}, 
                                                    {'W', std::vector<int> {48}}, 
                                                    {'Y', std::vector<int> {321,196}}, 

                                                    {'*', std::vector<int> {46,96,98}}, // STOP (UAG,UAA,UGA)
                                                    {'U', std::vector<int> {98}}, // Selenocysteine (UGA)
                                                    {'X', std::vector<int> {306,181,81,31,85,35,323,198,310,185,324,199,313,188,88,38,316,191,303,178,78,90,40,28,319,194,94,44,99,49,307,182,82,32,315,190,91,41,318,193,93,43,89,39,317,192,92,42,314,189,308,183,83,33,301,176,76,26,48,321,196,46,96,98}},  // Any (No constraint)
                                                    {'Z', std::vector<int> {91,41,85,35}},  // Q or E 
                                                    {'B', std::vector<int> {315,190,310,185}}, // N or D 

                                                    {'+', std::vector<int> {316,191,90,40,318,193,93,43,89,39}}, // Positively charged (H,K,R)
                                                    {'-', std::vector<int> {310,185,85,35}}, // Negatively charged (D,E)
                                                    {'=', std::vector<int> {313,188,88,38,306,181,81,31,319,194,94,44,301,176,76,26,99,49,307,182,82,32,303,178,78,28,48,324,199}}, // Hydrophobic (G,A,P,V,L,I,M,W,F)
                                                    {'&', std::vector<int> {317,192,92,42,314,189,308,183,83,33,321,196,315,190,91,41,323,198}}, // Polar (S,T,Y,N,Q,C)

                                                    {'@', std::vector<int> {324,199,48,321,196}}, // Aromatic (F,W,Y)
                                                    {'^', std::vector<int> {313,188,88,38,306,181,81,31,301,176,76,26,99,49,307,182,82,32,303,178,78}}, // Aliphatic (G,A,V,L,I)
                                                    {'$', std::vector<int> {317,192,92,42,314,189,323,198,308,183,83,33,28,98}}, // Hydroxyl or Sulfur/Selenium-containing (S,C,T,M,U)
                                                    {'%', std::vector<int> {310,185,85,35,315,190,91,41}},  //Acidic and their Amide (D,E,N,Q)

                                                    {'!', std::vector<int> {303,178,78,99,49,307,182,82,32,301,176,76,26,48,324,199,323,198}}, // Highly hydrophobic (I,L,V,W,F,C)
                                                    {'~', std::vector<int> {90,40,318,193,93,43,89,39}}, // Highly basic (K,R)
                                                    {'?', std::vector<int> {85,35,90,40,318,193,93,43,89,39,313,188,88,38,91,41,317,192,92,42,314,189,319,194,94,44,306,181,81,31}} // Enriched (E,K,R,G,Q,S,P,A
                                                    };  

const std::map<int,char>  codonAaDict = {{306,'A'},{181,'A'},{81,'A'},{31,'A'},
                                         {85,'E'},{35,'E'},
                                         {323,'C'},{198,'C'},
                                         {310,'D'},{185,'D'},
                                         {324,'F'},{199,'F'},
                                         {313,'G'},{188,'G'},{88,'G'},{38,'G'},
                                         {191,'H'},{316,'H'},
                                         {303,'I'},{178,'I'},{78,'I'},
                                         {90,'K'},{40,'K'},
                                         {28,'M'},
                                         {319,'P'},{194,'P'},{94,'P'},{44,'P'},
                                         {99,'L'},{49,'L'},{307,'L'},{182,'L'},{82,'L'},{32,'L'},
                                         {315,'N'},{190,'N'},
                                         {91,'Q'},{41,'Q'},
                                         {318,'R'},{193,'R'},{93,'R'},{43,'R'},{89,'R'},{39,'R'},
                                         {317,'S'},{192,'S'},{92,'S'},{42,'S'},{314,'S'},{189,'S'},
                                         {308,'T'},{183,'T'},{83,'T'},{33,'T'},
                                         {301,'V'},{176,'V'},{76,'V'},{26,'V'},
                                         {48,'W'},
                                         {321,'Y'},{196,'Y'},
                                         {98,'U'},
                                         {96,'*'},{46,'*'}};
                                         

class AAConstraint{
	public:			
		AAConstraint(int startPos, std::string  aaTarget, std::string  aaSeq, int maxBlosumScore, int aaSimilCstr, int strLen);
		~AAConstraint();

		std::string getTarget();
		char getTrgConstAt(int pos);
		std::string getAllowedAaAt(int pos);

		std::vector<int> getDomain(int pos);

		int getSimilarity();		
		bool maximizeScore();
		int getStartPos();
		int getLength();
		int isValid();		
		void toString();
		std::string getErrorMessage();
		
		int getMaxBlosumValue(int pos);
		int getBlosumIndex(int pos);		

	private:
		int startPos_;
		std::string aaTarget_;
		std::string aaSeq_;
		int maxBlosumScore_;
		int aaSimilCstr_;
		int strLen_;

		int size_;
		std::vector<std::string > allowedAAs;
		std::vector<std::vector<int> > codonDomains;

		int errCode_;
		Blosum* blosumHelper;

};

// Auxiliary functions 
std::string getAAConstraintsByClassification(std::string aaTarget, int aaSimilCstr);
void parseAAconstraints(std::vector<AAConstraint*> *aaConstraints, std::string aaTarget, std::string aaConst,std::string aaStartPos, int maxBlosumScore, int aaSimilCstr, int strLen);

#endif

