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
#ifdef _WIN32
	#include <windows>
#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	#include <unistd.h>
	#include <limits.h>
#endif

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <sstream>
#include"misc.h"

#define DEBUG 0 

int d2[13] = {0,13,0,14,0,0,0,15,0,0,0,0,16};
int simple[13] = {-1,2,-1,0,-1,-1,-1,1,-1,-1,-1,-1,3};

double R = 0.001987;
double T = 310.00;

/* 
  d1 breakdown

  00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 
  AU UA GC CG GU UG AC AG CA CU GA UC GG AA CC UU

  00 01 02 03 04 05 06 07
  0A 0U 0C 0G A0 U0 C0 G0
*/
                  // 0 1  2 3  4 5 6  7 8 9 10 11  12 13 14 15 16  17 18 19  20 21  22  23 24  25  26  27  28
int d1[29] = {1,4,-1,0,-1,2,2, 9,6,0,10, 3, -1, 7, 7, 4, 8, 11, 3, 6, -1, 1, -1,  5, 5, 12, 13, 14, 15};

//cannonical basepairs indexed by d1[getIndex(i,j)]
int cbp[29] = {1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

int guau[29] = {1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

// function takes structure and converts to 1 indes base pair array and store length at the first position
int* make_BasePair_Table(std::string structure, bool withUndet){
	int i, j,hx;
	int *stack, *table;

	hx=0;
	stack = (int *) malloc(sizeof(int)*(structure.length()+1));
	table = (int *) malloc(sizeof(int)*(structure.length()+1));
	table[0]=structure.length();
	for (i=0; i<structure.length(); i++){
		switch (structure[i]){
		case '.':
		     table[i+1]= -1;
		     break;
		case ',':
		     if(withUndet){
				table[i+1]= -2;
		     }
		     else{
				table[i+1]= -1;
		     }
		     break;

		case '(':
		     stack[hx++]=i;
		     break;
		case ')':
		     j = stack[--hx];
		     if (hx<0){
		     	//printf("%s \n",structure);
		     	printf("Unbalanced brackets in make_BasePair_Table\n");
		     	exit(0);
		     }
		     table[i+1]=j+1;
		     table[j+1]=i+1;
		     break;
		}
	}
	if (hx!=0){
		//printf("%s \n",structure);
		printf("Unbalanced brackets in make_BasePair_Table\n");
		exit(0);
	}
	free(stack);
	return table;
}

char ToNucl(int trans){

	char nucl = ' ';

	if (trans == 0){
		nucl = 'A';
	}
	else if(trans == 1){
		nucl = 'C';
	}
	else if(trans == 2){
		nucl = 'G';
	}
	else if(trans == 3){
		nucl = 'U';
	}

	return nucl;
}

int Nt2Int(char nucl){

	int trans = 0;
	switch(nucl){
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'U':
			return 3;
		default:
			printf("Invalid nucleotide in conversion %c", nucl);
			exit(0);
	}
	return trans;
}

// function takes sequence and converts to integer array through dictionary A->0 C->1, G->3, U->12
int * IntSeq(char * sequence)
{
  int * integer_seq;
  int i;
  char * start;
  integer_seq = (int *)malloc(strlen(sequence)*sizeof(int));
  i = 0; 
  start = sequence;
  while(*start != '\0'){
 
    if (*start == '0'){
      integer_seq[i] = strlen(sequence)-1;
    }
    else if(*start == 'A'){
      integer_seq[i] = 0;
    }
    else if(*start == 'C'){
      integer_seq[i] = 1;
    }
    else if(*start == 'G'){
      integer_seq[i] = 2;
    }
    else if(*start == 'U'){
      integer_seq[i] = 3;
    }
   
    start = start+ sizeof(char);
    i++;
  }
  return integer_seq;
}


// function takes sequence and converts to integer array through dictionary G->1, A->3, C->7, U->12
int * TransformSeq(char * sequence)
{
  int * integer_seq;
  int i;
  char * start;
  integer_seq = (int *)malloc(strlen(sequence)*sizeof(int));
  i = 0; 
  start = sequence;
  while(*start != '\0'){
 
    if (*start == '0'){
      integer_seq[i] = strlen(sequence)-1;
    }
    else if(*start == 'G'){
      integer_seq[i] = 1;
    }
    else if(*start == 'A'){
      integer_seq[i] = 3;
    }
    else if(*start == 'C'){
      integer_seq[i] = 7;
    }
    else if(*start == 'U'){
      integer_seq[i] = 12;
    }
   
    start = start+ sizeof(char);
    i++;
  }
  return integer_seq;
}


// reason we have d2 is to account for AA,GG,CC,UU, +12 because integer value for U is 12 
int GetIndex(int i, int j)
{
  int index;
  index = i-j+12+((i==j)*d2[i]);
  return index;
}

void PrintStr(char *ptr)
{
  char *start;
  start = ptr;
  start = start+ sizeof(char);
  while(*start != '\0'){
    printf("%c",*start);
    start = start + sizeof(char);
  }
  printf("\n");
}

void PrintInt(int *ptr, int len)
{
  int i;
  for(i = 0; i < len+1; i++){
    printf("%d", ptr[i]);
  }
  printf("\n");
}



char GetNucl(int trans){

  char nucl;

  if (trans == 1){
    nucl = 'G';
  }

  else if(trans == 3){
    nucl = 'A';
  }

  else if(trans == 7){
    nucl = 'C';
  }

  else if(trans == 12){
    nucl = 'U';
  }

  return nucl;
}

int CountChar(char * str, char c){

  int count;
  count = 0;

  while(*str != '\0'){
    if(*str == c){
      count++;
    }
    str++;
  }
  
  return count;
}

int **Allocate2Dmatrix(int a, int b){
  int i = 0;
  int j = 0;
  int **Matrix;
  

  Matrix = (int **) calloc(a,sizeof(int *));
  if(Matrix == NULL){
    printf("out of memory\n");
    exit(1);
  }
  for (i = 0; i <a; i++){
    Matrix[i] = (int *) calloc(b,sizeof(int));
    if(Matrix[i] == NULL){
      printf("out of memory\n");
      exit(1);
    }
  }

  return Matrix;
}

double **Allocate2DMatrixDouble(int a, int b){
  int i = 0;
  int j = 0;
  double **Matrix;

  Matrix = (double **) calloc(a,sizeof(double *));
  if(Matrix == NULL){
    printf("out of memory\n");
    exit(1);
  }
  for (i = 0; i <a; i++){
    Matrix[i] = (double *) calloc(b,sizeof(double));
    if(Matrix[i] == NULL){
      printf("out of memory\n");
      exit(1);
    }
  }
  return Matrix;
}

int **BPBetweenBool(int** Matrix, int * bp, int len){

  int pos,t,x,i,j,k,l,i_index,k_index;

  for(i_index = 0; i_index < len; i_index++){

    if(bp[i_index] != -1){
   
      i = i_index;
      j = bp[i];

      for(k_index=i_index+1;k_index < j; k_index++){
        if (bp[k_index] != -1){
          k = k_index;
          l = bp[k];
          break;
        }
      }
  
      if(k > i && j > l){
        Matrix[i][j] = 1;
      }
    
      for(t=l+1;t<j;t++){
        x = bp[t];
        if(x!=-1){
          Matrix[l][j] = 1;
          break;
        }
      }
    }
  }
  
  if(DEBUG){
    for(i = 0; i < len; i++){
      for(j = 0; j < len; j++){
        printf("%d ",Matrix[i][j]);
      }
      printf("\n");
    }
  }  
  
  return Matrix;
}

int * GetBPList(char * str, int len){

  int i = 0;
  int pointer = 0;
  int * left;
  int * bp;
  int leftbp;

  left = (int *)malloc(len*sizeof(int));
  bp = (int *)malloc(len*sizeof(int));

  for (i = 0; i < len; i++){
    bp[i] = -1;
  }
  
  for (i = 0; i < len; i++){

    if (str[i] == '('){
      left[pointer] = i;
      pointer++;
    }
  
    else if(str[i] == ')'){
      if (pointer == 0){
        printf("Unbalanced Secondary Structure\n");
        exit(1);
      }
      leftbp = left[pointer-1];
      bp[leftbp] = i;
      bp[i] = leftbp;
      pointer--;
    }
  }

  if(DEBUG){
   
    for(i = 0; i < len; i++){
      printf("%d %d\n",i,bp[i]);
    }
  
  }

  free(left);

  return bp;
}

int linearize2D(int row, int column){

  int index;
  index = row*16+column;

  return index;
}


bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
int parseTriple(std::string strTuple, int maxValue, std::string name, std::vector<std::tuple<int,int,int>>* tTuple){
	if(strTuple.length()>0){
		std::vector<std::string> vTuple;
		split(strTuple,',', vTuple);
		if(vTuple.size()>0){
			for(int i=0;i<vTuple.size();i++){
				std::vector<std::string> vTupleElems;
				split(vTuple[i],' ', vTupleElems);
				if(vTupleElems.size()==1){
					if(!is_number(vTupleElems[0])){
						printf("Invalid syntax in %s: Number of elements (%s) is not a positive number\n",name.c_str(),vTupleElems[0].c_str());
						return 0;
					}

					int nElems=atoi(vTupleElems[0].c_str());
					if(nElems<0 || nElems>maxValue){
						printf("Invalid syntax in %s: Number of elements must be between 0 and %d\n",name.c_str(),maxValue);
						return 0;
					}
					else{
						tTuple->push_back(std::make_tuple(nElems, 1, maxValue));
					}			

				}
				else if(vTupleElems.size()==3){
					if(!is_number(vTupleElems[0])){
						printf("Invalid syntax in %s: Number of elements (%s) is not a positive  number\n",name.c_str(),vTupleElems[0].c_str());
						return 0;
					}
					if(!is_number(vTupleElems[1])){
						printf("Invalid syntax in %s: First position (%s) is not a positive number\n",name.c_str(),vTupleElems[1].c_str());
						return 0;
					}
					if(!is_number(vTupleElems[2])){
						printf("Invalid syntax in %s: Last position (%s) is not a positive number\n",name.c_str(),vTupleElems[2].c_str());
						return 0;
					}
					int nElems=atoi(vTupleElems[0].c_str());
					int firstPos=atoi(vTupleElems[1].c_str());
					int lastPos=atoi(vTupleElems[2].c_str());
					if(nElems<0 || nElems>maxValue){
						printf("Invalid syntax in %s: Number of elements must be between 0 and %d\n",name.c_str(),maxValue);
						return 0;
					}
					if(firstPos < 1 || firstPos > maxValue || lastPos < 1 || lastPos > maxValue ){
						printf("Invalid syntax in %s: Starting and end positions must be between must be between 1 and %d\n",name.c_str(),maxValue);
						return 0;
					}
					if(firstPos > lastPos){
						printf("Invalid syntax in %s: Starting positions must be lower than ending position\n",name.c_str());
						return 0;
					}
					if(firstPos+nElems-1 > lastPos){
						printf("Invalid syntax in %s: Number of nucleotides is greater than the range specified\n",name.c_str());
						return 0;
					}
					
					tTuple->push_back(std::make_tuple(nElems, firstPos, lastPos));
				}
				else{
					printf("Invalid syntax in %s: Each constraint must contain 1 or 3 numbers separated by spaces\n",name.c_str());
					return 0;
				}
			}
		}		
	}
//	for(int i=0;i< tTuple->size();i++){
//		printf("%s: %d %d %d\n",name.c_str(),std::get<0>(tTuple->at(i)),std::get<1>(tTuple->at(i)),std::get<2>(tTuple->at(i)));
//	}
	return 1;
}

int validSecondaryStructure(std::string secStr, bool acceptCommas){
	
	int balance=0;
	int numberOfBPs=0;
	int nCutPoints=0;
	int nCommas=0;
	for(int i=0; i<secStr.length();i++){
		switch (secStr[i]){
			case '(':
				numberOfBPs++;
				balance++;
				break;
			case ')':
				balance--;
				break;
			case '.':
				break;
			case '&':
				nCutPoints++;
				if(nCutPoints>1){
					return -1;
				}
				break;
			case ',':
				if(!acceptCommas){
					return -1;
				}
				nCommas++;
				break;
			default:
				return -1;
				break;
				
		}
		if(balance<0){
			return -1;
		}
	}
	
	if(balance!=0){
		return -1;
	}
	if(acceptCommas){
		return numberOfBPs +nCommas/2;
	}
	else{
		return numberOfBPs;	
	}
}

int findCutPoint(std::string* secStr){
		int cutPoint = -1;
		for(int i=0; i<secStr->length();i++){
			if(secStr->at(i)=='&'){
				cutPoint=i;
			}
		}
		
		if(cutPoint!=-1){
			secStr->erase(cutPoint,1);
		}
		
		return cutPoint;
}

std::string getExecPath(std::string argv0){
	std::string exec_path=".";
	#if defined(_WIN32)
	// GetModuleFileName
	//_pgmptr
	#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	char buff[PATH_MAX];
	int bufsize = PATH_MAX-1;
	if(FILE *file = fopen("/proc/self/exe", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/exe", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			exec_path = std::string(buff);
		}		
	}
	else if(FILE *file = fopen("/proc/curproc/file", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/curproc/file", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			exec_path = std::string(buff);
		}		
	}
	else if(FILE *file = fopen("/proc/self/path/a.out", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/path/a.out", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			exec_path = std::string(buff);
		}		
	}
	else{
		exec_path= argv0;
	}
	size_t slash_position= exec_path.find_last_of('/');
	if(slash_position != std::string::npos){			
		exec_path=exec_path.substr(0,slash_position);
	}
	else{
		exec_path=".";
	}

	#endif
	
	return exec_path;
}
