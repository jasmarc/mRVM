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
    Predictor(Matrix *w, Vector *t, Matrix *x_train, Matrix* x_predict);
    virtual ~Predictor();
    void Predict();
  private:
    void QuadratureApproximation();
    Kernel *k;
    Matrix *w;
    Vector *t;
};
}

#endif  // SRC_LIB_PREDICTOR_H_