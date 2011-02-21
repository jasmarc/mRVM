// Copyright 2011 Jason Marcell

#ifndef KERNEL_H_
#define KERNEL_H_

#include "lib/Matrix.h"

namespace jason {

class Kernel: public jason::Matrix {
  public:
    Kernel(Matrix *m1, Matrix *m2);
    virtual ~Kernel();
};
}

#endif  // KERNEL_H_
