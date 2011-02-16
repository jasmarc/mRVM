// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_TRAINER_H_
#define SRC_LIB_TRAINER_H_

#include "lib/Matrix.h"

namespace jason {

class Matrix;

class Trainer {
  public:
    explicit Trainer(Matrix *matrix, Vector *labels, size_t classes);
    virtual ~Trainer();
    void Process();

  private:
    Matrix *x;  // Data Points
    Vector *t;  // Labels
    size_t samples, features, classes;

    Matrix *w;
    Matrix *k;
    Matrix *a;
    Matrix *y;

    Matrix *BuildKernel();
    void InitializeYAW();
    void UpdateA(double tau, double v);
    void UpdateW();
    void UpdateY();
};
}

#endif  // SRC_LIB_TRAINER_H_
