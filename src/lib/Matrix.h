/*
 * Matrix.h
 *
 *  Created on: Feb 10, 2011
 *      Author: jason
 */
#include <gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#ifndef MATRIX_H_
#define MATRIX_H_

namespace jason {

class Matrix {
public:
	Matrix();
	Matrix(double *data, size_t height, size_t width);
	Matrix(gsl_matrix* data);
	Matrix(const char* filename);
	virtual ~Matrix();
	size_t Height();
	size_t Width();
	double *Data();
	Matrix *Clone();
	Matrix *Invert();
	double Get(int row, int col);
	void Set(int row, int col, double val);
	void Print();
	Matrix *Multiply(Matrix *other);
private:
	int NumberOfRows(FILE *f);
	int NumberOfColumns(FILE *f);
	gsl_matrix *CloneGSLMatrix();
	gsl_matrix* m;
};

}

#endif /* MATRIX_H_ */
