#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <gsl_blas.h>
#include <gsl_linalg.h>
#include <gsl_matrix.h>
#include "Vector.h"

namespace jason {

class Vector;

class Matrix {
public:
	Matrix(Vector *vec);
	Matrix(double *data, size_t height, size_t width);
	Matrix(const char* filename);
	virtual ~Matrix();
	size_t Height();
	size_t Width();
	Matrix *Clone();
	Matrix *Invert();
	double Get(int row, int col);
	void Set(int row, int col, double val);
	Vector *Row(int row);
	Vector *Column(int col);
	void Print();
	Matrix* Multiply(Matrix *other);
	friend class Vector;
private:
	Matrix(gsl_matrix *mat);
	int NumberOfRows(FILE *f);
	int NumberOfColumns(FILE *f);
	gsl_matrix *CloneGSLMatrix();
	gsl_matrix* m;
};

}

#endif /* MATRIX_H_ */
