// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_PREDICTOR_H_
#define SRC_LIB_PREDICTOR_H_

#include "lib/Matrix.h"

namespace jason {

class Matrix;

class Predictor {
  public:
    Predictor(Matrix *w, Matrix *x_train, Matrix* x_predict);
    virtual ~Predictor();
};
}

#endif  // SRC_LIB_PREDICTOR_H_
