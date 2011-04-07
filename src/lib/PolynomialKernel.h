// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_POLYNOMIALKERNEL_H_
#define SRC_LIB_POLYNOMIALKERNEL_H_

#include "lib/Vector.h"
#include "lib/Matrix.h"
#include "lib/Kernel.h"

namespace jason {

class Matrix;
class Vector;
class Kernel;

class PolynomialKernel: public jason::Kernel {
  public:
    PolynomialKernel(Matrix *m1, Matrix *m2, int n);
    virtual ~PolynomialKernel();
    double KernelElementFunction(Vector *vec1, Vector *vec2);
  private:
    int n;
};
}

#endif  // SRC_LIB_POLYNOMIALKERNEL_H_
