// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_KERNEL_H_
#define SRC_LIB_KERNEL_H_

#include "lib/Matrix.h"

namespace jason {

class Matrix;

class Kernel: public jason::Matrix {
  public:
    Kernel();
    Kernel(Matrix *m1, Matrix *m2);
    virtual ~Kernel();
    virtual double KernelElementFunction(Vector *vec1, Vector *vec2) = 0;
};
}

#endif  // SRC_LIB_KERNEL_H_
