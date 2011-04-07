// Copyright 2011 Jason Marcell

#include <math.h>

#include <gsl/gsl_matrix.h>

#include "lib/LinearKernel.h"
#include "lib/Log.h"

namespace jason {

LinearKernel::LinearKernel(Matrix *m1, Matrix *m2) : Kernel(m1, m2) {
  LOG(DEBUG, "Linear Kernel constructor with params.\n");
}

LinearKernel::~LinearKernel() {
  LOG(DEBUG, "Linear Kernel constructor with no params.\n");
}

double LinearKernel::KernelElementFunction(Vector *vec1, Vector *vec2) {
  LOG(DEBUG, "Linear KernelElementFunction.\n");
  return vec1->Multiply(vec2);
}
}
