// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_LINEARKERNEL_H_
#define SRC_LIB_LINEARKERNEL_H_

#include "lib/Vector.h"
#include "lib/Matrix.h"
#include "lib/Kernel.h"

namespace jason {

class Matrix;
class Vector;
class Kernel;

class LinearKernel: public jason::Kernel {
  public:
    LinearKernel(Matrix *m1, Matrix *m2);
    virtual ~LinearKernel();
    double KernelElementFunction(Vector *vec1, Vector *vec2);
};
}

#endif  // SRC_LIB_LINEARKERNEL_H_
