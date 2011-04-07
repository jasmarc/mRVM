// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_KERNEL_H_
#define SRC_LIB_KERNEL_H_

#include "lib/Matrix.h"

namespace jason {

class Matrix;

enum KernelType { LINEAR, POLYNOMIAL, GAUSSIAN };

class Kernel: public jason::Matrix {
  public:
    Kernel(Matrix *m1, Matrix *m2);
    virtual ~Kernel();
    void Init();
    virtual double KernelElementFunction(Vector *vec1, Vector *vec2) = 0;
  protected:
    Matrix *m1;
    Matrix *m2;
};
}

#endif  // SRC_LIB_KERNEL_H_
