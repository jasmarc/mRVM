// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/GaussianKernel.h"
#include "lib/Log.h"

namespace jason {

GaussianKernel::GaussianKernel(Matrix *m1, Matrix *m2, int param)
  : Kernel(m1, m2) {
  LOG(DEBUG, "Gaussian Kernel constructor with params.\n");
  Vector *theta = new Vector(m1->Width());
  for (size_t i = 0; i < theta->Size(); ++i) {
    theta->Set(i, static_cast<double>(param));
  }
  this->theta = new Matrix(theta);
  delete theta;
}

GaussianKernel::~GaussianKernel() {
}

double GaussianKernel::KernelElementFunction(Vector *vec1, Vector *vec2) {
  double ret = 1;
  Vector *v1 = vec1->Subtract(vec2);
  Vector *v2 = v1->Multiply(theta);
  ret = v2->Multiply(v1);
  delete v1;
  delete v2;
  return exp(-0.5 * ret);
}
}
