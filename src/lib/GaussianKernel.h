// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_GAUSSIANKERNEL_H_
#define SRC_LIB_GAUSSIANKERNEL_H_

#include "lib/Vector.h"
#include "lib/Matrix.h"
#include "lib/Kernel.h"

namespace jason {

class Matrix;
class Vector;
class Kernel;

class GaussianKernel: public jason::Kernel {
  public:
    explicit GaussianKernel(Vector *theta);
    GaussianKernel(Matrix *m1, Matrix *m2, int param);
    virtual ~GaussianKernel();
    double KernelElementFunction(Vector *vec1, Vector *vec2);
  private:
    Matrix *theta;
};
}

#endif  // SRC_LIB_GAUSSIANKERNEL_H_
