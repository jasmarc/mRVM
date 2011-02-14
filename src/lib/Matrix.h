#ifndef MATRIX_H_
#define MATRIX_H_

#include "Vector.h"
#include <iostream>
#include <gsl_matrix.h>

namespace jason {

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
	int NumberOfRows(FILE *f);
	int NumberOfColumns(FILE *f);
	gsl_matrix *CloneGSLMatrix();
	gsl_matrix* m;
};

}

#endif /* MATRIX_H_ */
