// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_PREDICTOR_H_
#define SRC_LIB_PREDICTOR_H_

#include "lib/Matrix.h"
#include "lib/Kernel.h"

namespace jason {

class Matrix;
class Kernel;

class Predictor {
  public:
    Predictor(Matrix *w, Matrix *x_train, Matrix* x_predict, Kernel *kernel);
    virtual ~Predictor();
    Matrix* Predict();
  private:
    Matrix *QuadratureApproximation();
    Kernel *k;
    Matrix *w;
};
}

#endif  // SRC_LIB_PREDICTOR_H_
