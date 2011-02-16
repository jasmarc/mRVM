// Copyright 2011 Jason Marcell

#include <stdio.h>

#include "lib/Matrix.h"
#include "lib/Trainer.h"

using jason::Matrix;
using jason::Vector;
using jason::Trainer;

int main() {
  double m1_arr[] = { 1, 2, 9,
                      4, 6, 5,
                      3, 7, 4,
                      6, 2, 0 };

  double v_arr[] = { 0, 1, 0, 3 };

  double m2_arr[] = { 20.25, -2.39, -1.92, -14.93,
      -2.39, 2.33, 2.03, -0.97,
      -1.92, 2.03, 4.38, -3.50,
      -14.93, -0.97, -3.50, 20.41 };

  double m3_arr[] = {1, 2, 3, 4};

  printf("Starting...\n");

  Matrix *m1 = new Matrix(m1_arr, 4, 3);
  Vector *v = new Vector(v_arr, 4);
  Trainer *t = new Trainer(m1, v, 4);
  t->Process();

  printf("End.\n");
  return 0;
}
