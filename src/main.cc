// Copyright 2011 Jason Marcell

#include <stdio.h>

#include "lib/Matrix.h"

using jason::Matrix;
using jason::Vector;

int main() {
  double m1_arr[] = { 1, 2, 9, 4, 6, 5, 3, 7, 4, 6, 2, 0 };

  double m2_arr[] = { 1, 2, 9, 4, 6, 5, 3, 7, 4 };

  double v_arr[] = { 1, 2, 3 };

  Matrix *m1 = new Matrix(m1_arr, 4, 3);
  printf("Starting...\n");
  m1->Sphere();
  m1->Print();
  printf("End.\n");
  return 0;
}
