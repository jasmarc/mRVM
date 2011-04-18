// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_TRAINER_H_
#define SRC_LIB_TRAINER_H_

#include "lib/Matrix.h"
#include "lib/Kernel.h"

namespace jason {

class Matrix;
class Kernel;

class Trainer {
  public:
    explicit Trainer(Matrix *matrix, Vector *labels, size_t classes,
        Kernel *kernel);
    virtual ~Trainer();
    void Process(double tau, double upsilon);
    Matrix *GetW();

  private:
    Matrix *x;  // Data Points
    Vector *t;  // Labels
    size_t samples, features, classes;

    bool converged;

    Matrix *w;
    Kernel *k;
    Matrix *a;
    Matrix *y;

    void InitializeYAW();
    void UpdateA(double tau, double upsilon);
    void UpdateW();
    void UpdateY();
};
}

#endif  // SRC_LIB_TRAINER_H_
