// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/Kernel.h"
#include "lib/Matrix.h"
#include "lib/Log.h"

namespace jason {

Kernel::Kernel(Matrix *m1, Matrix *m2) : Matrix(m1->Height(), m2->Height()) {
  LOG(DEBUG, "Base Kernel constructor with params.\n");
  this->m1 = m1;
  this->m2 = m2;
}

Kernel::~Kernel() {
}

void Kernel::Init() {
  LOG(DEBUG, "= Begin Base Kernel Init. =\n");
  this->m = gsl_matrix_alloc(m1->Height(), m2->Height());
  for (size_t row = 0; row < this->Height(); ++row) {
    for (size_t col = 0; col < this->Width(); ++col) {
      Vector* vec1 = m1->Row(row);
      Vector* vec2 = m2->Row(col);
      double elem = this->KernelElementFunction(vec1, vec2);
      double s1 = this->KernelElementFunction(vec1, vec1);
      double s2 = this->KernelElementFunction(vec2, vec2);
      elem = elem / sqrt(s1 * s2);
      this->Set(row, col, elem);
      delete vec2;
      delete vec1;
    }
  }
  LOG(DEBUG, "= End Base Kernel Init. =\n");
}
}
