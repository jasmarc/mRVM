// Copyright 2011 Jason Marcell

#include <stdio.h>

#include "lib/Matrix.h"
#include "lib/Kernel.h"
#include "lib/Trainer.h"
#include "lib/Predictor.h"
#include "lib/GaussHermiteQuadrature.h"

using jason::Matrix;
using jason::Vector;
using jason::Trainer;
using jason::Kernel;
using jason::Predictor;
using jason::GaussHermiteQuadrature;

int main() {
  double m1_arr[] = { 1, 2, 9,
                      4, 6, 5,
                      3, 7, 4,
                      6, 2, 0 };

  double m2_arr[] = { 1, 2, 9,
                      4, 6, 5,
                      6, 2, 0 };

  double v_arr[] = { 0, 1, 0, 3 };

  printf("Starting...\n");

  Matrix *m1 = new Matrix(m1_arr, 4, 3);
  Vector *v = new Vector(v_arr, 4);
  Trainer *t = new Trainer(m1, v, 4);
  t->Process();

  Matrix *m2 = new Matrix(m2_arr, 3, 3);
  Predictor *p = new Predictor(t->GetW(), v, m1, m2);
  p->Predict();

  printf("End.\n");
  return 0;
}
