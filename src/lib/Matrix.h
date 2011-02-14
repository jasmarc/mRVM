// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_MATRIX_H_
#define SRC_LIB_MATRIX_H_

#include <stddef.h>

#include <gsl_blas.h>
#include <gsl_linalg.h>
#include <gsl_matrix.h>

#include "lib/Vector.h"

namespace jason {

class Vector;

class Matrix {
  public:
    explicit Matrix(Vector *vec);
    explicit Matrix(double *data, size_t height, size_t width);
    explicit Matrix(const char* filename);
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
    explicit Matrix(gsl_matrix *mat);
    int NumberOfRows(FILE *f);
    int NumberOfColumns(FILE *f);
    gsl_matrix *CloneGSLMatrix();
    gsl_matrix* m;
};
}

#endif  // SRC_LIB_MATRIX_H_
