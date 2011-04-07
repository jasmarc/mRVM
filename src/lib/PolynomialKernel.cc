// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/PolynomialKernel.h"
#include "lib/Log.h"

namespace jason {

PolynomialKernel::PolynomialKernel(int n) {
  LOG(DEBUG, "Polynomial Kernel constructor with one param.\n");
  // TODO(jrm) danger, the m variable isn't getting set
  this->n = n;
}

PolynomialKernel::PolynomialKernel(Matrix *m1, Matrix *m2, int n)
  : Kernel(m1, m2) {
  LOG(DEBUG, "Polynomial Kernel constructor with params.\n");
  this->n = n;
}

PolynomialKernel::~PolynomialKernel() {
}

double PolynomialKernel::KernelElementFunction(Vector *vec1, Vector *vec2) {
  return pow(vec1->Multiply(vec2) + 1, n);
}
}
