

#ifndef ARRAYCLASS_H
#define ARRAYCLASS_H

#include "defines.h"

// arrayclass encapsulates the large 2-d arrays of w and v, used by
// the dynamic programming algorithm 
// shaped like this:
//      -----
//     |     |
//      |     |
//  i    |     |
//        |     |
//         -----
//            j

class arrayclass {
  private:
	int Size;

  public:
	int k;
	integersize **dg;
	integersize infinite;

	// the constructor allocates the space needed by the arrays
	arrayclass(int size, integersize energy = INFINITE_ENERGY);

	// the destructor deallocates the space used
	~arrayclass();

	// f is an integer function that references the correct element of
	// the array
	integersize &f(int i, int j);
};

inline integersize &arrayclass::f(int i, int j) {
	if (i > Size) {
		i -= Size;
		j -= Size;
	}

	if (i > j) {
		return infinite;
	}

	return dg[i][j];
}

//transpose of array class, for efficient traversal of columns
class arrayclassT {
  private:
	int Size;

  public:
	int k;
	integersize **dg;
	integersize infinite;

	// the constructor allocates the space needed by the arrays
	arrayclassT(int size, integersize energy = INFINITE_ENERGY);

	// the destructor deallocates the space used
	~arrayclassT();

	// f is an integer function that references the correct element of
	// the array
	integersize &f(int i, int j);
};

inline integersize &arrayclassT::f(int i, int j) {
	if (i > Size) {
		i -= Size;
		j -= Size;
	}

	if (i > j) {
		return infinite;
	}

	return dg[j][i];
}
#endif
