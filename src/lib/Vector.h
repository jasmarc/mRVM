// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_VECTOR_H_
#define SRC_LIB_VECTOR_H_

#include <gsl_matrix.h>

#include "lib/Matrix.h"

namespace jason {

class Matrix;

class Vector {
  public:
    explicit Vector(size_t size);
    explicit Vector(double *data, size_t size);
    virtual ~Vector();
    size_t Size();
    double Get(size_t elem);
    void Set(size_t elem, double value);
    Matrix *RepmatVert(size_t k);
    Matrix *RepmatHoriz(size_t k);
    double Multiply(Vector *other);
    Vector *Multiply(Matrix *m);
    Vector *Subtract(Vector *other);
    Vector *Add(Vector *other);
    void Print();
    friend class Matrix;
  private:
    explicit Vector(gsl_vector *v);
    gsl_vector *v;
};
}

#endif  // SRC_LIB_VECTOR_H_
