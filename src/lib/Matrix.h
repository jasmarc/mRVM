/*
 * Matrix.h
 *
 *  Created on: Feb 10, 2011
 *      Author: jason
 */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#ifndef MATRIX_H_
#define MATRIX_H_

namespace jason {

class Matrix {
public:
	Matrix();
	Matrix(char* filename);
	virtual ~Matrix();
	void Print();
private:
	int NumberOfRows(FILE *f);
	int NumberOfColumns(FILE *f);
	gsl_matrix* m;
};

}

#endif /* MATRIX_H_ */
