// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_MATRIX_H_
#define SRC_LIB_MATRIX_H_

#include <stddef.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "lib/Vector.h"

namespace jason {

class Vector;

class Matrix {
  public:
    explicit Matrix(double *data, size_t height, size_t width);
    explicit Matrix(size_t height, size_t width);
    explicit Matrix(Vector *vec);
    explicit Matrix(const char* filename);
    virtual ~Matrix();
    void Write(const char* filename);
    size_t Height();
    size_t Width();
    void Invert();
    double Get(int row, int col);
    void Set(int row, int col, double val);
    void Add(Matrix *other);
    Vector *Row(size_t row);
    Vector *Column(size_t col);
    void SetColumn(size_t col, Vector *vec);
    void SetRow(size_t row, Vector *vec);
    void Sphere();
    void Sphere(Matrix *other);
    void NormalizeResults();
    void CacheMeansAndStdevs();
    Vector* GetMeans();
    Vector* GetStdevs();
    char *ToString();
    Matrix* Multiply(Matrix *other);
    Vector* Multiply(Vector *vec);
    friend class Vector;
    friend class Kernel;
  private:
    void Init();
    explicit Matrix(gsl_matrix *mat);
    size_t NumberOfRows(FILE *f);
    size_t NumberOfColumns(FILE *f);
    gsl_matrix* m;
    Vector *means;
    Vector *stdevs;
    char *to_str;
};
}

#endif  // SRC_LIB_MATRIX_H_
