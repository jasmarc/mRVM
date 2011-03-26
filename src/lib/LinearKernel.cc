// Copyright 2011 Jason Marcell

#include <math.h>

#include <gsl/gsl_matrix.h>

#include "lib/LinearKernel.h"

namespace jason {

LinearKernel::LinearKernel(Matrix *m1, Matrix *m2) : Kernel(m1, m2) {
}

LinearKernel::~LinearKernel() {
}

double LinearKernel::KernelElementFunction(Vector *vec1, Vector *vec2) {
  return vec1->Multiply(vec2);
}
}
