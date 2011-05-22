// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_TRAINER2_H_
#define SRC_LIB_TRAINER2_H_

#include "lib/Matrix.h"
#include "lib/Kernel.h"

namespace jason {

class Matrix;
class Kernel;

class Trainer2 {
  public:
    explicit Trainer2(Matrix *matrix, Vector *labels, size_t classes,
        Kernel *kernel);
    virtual ~Trainer2();
    void Process(double tau, double upsilon);
    double Theta(double i, Matrix *q, Vector *s);
    Matrix *GetW();

  private:
    Matrix *x;  // Data Points
    Vector *t;  // Labels
    size_t samples, features, classes;

    bool converged;

    Matrix *w;
    Kernel *k;
    Vector *a;
    Matrix *y;
    Vector *active_samples;

    void InitializeYAW();
    size_t GetFirstSampleIndex();
    void UpdateA(double tau, double upsilon);
    void UpdateW();
    void UpdateY();
};
}

#endif  // SRC_LIB_TRAINER2_H_
