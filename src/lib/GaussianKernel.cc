// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/GaussianKernel.h"

namespace jason {

GaussianKernel::GaussianKernel(Vector *theta) {
  // TODO(jrm) danger, the m variable isn't getting set
  this->theta = new Matrix(theta);
}

GaussianKernel::GaussianKernel(Matrix *m1, Matrix *m2, Vector *theta)
  : Kernel(m1, m2) {
  this->theta = new Matrix(theta);
}

GaussianKernel::~GaussianKernel() {
}

double GaussianKernel::KernelElementFunction(Vector *vec1, Vector *vec2) {
  double ret;
  Vector *v1 = vec1->Subtract(vec2);
  Vector *v2 = v1->Multiply(theta);
  ret = v2->Multiply(v1);
  delete v1;
  delete v2;
  return exp(-0.5 * ret);
}
}
