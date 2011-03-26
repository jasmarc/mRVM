// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/PolynomialKernel.h"

namespace jason {

PolynomialKernel::PolynomialKernel(int n) {
  // TODO(jrm) danger, the m variable isn't getting set
  this->n = n;
}

PolynomialKernel::PolynomialKernel(Matrix *m1, Matrix *m2, int n)
  : Kernel(m1, m2) {
  this->n = n;
}

PolynomialKernel::~PolynomialKernel() {
}

double PolynomialKernel::KernelElementFunction(Vector *vec1, Vector *vec2) {
  return pow(vec1->Multiply(vec2) + 1, n);
}
}
