// Copyright 2011 Jason Marcell

#include "math.h"

#include "lib/Kernel.h"

namespace jason {

Kernel::Kernel(Matrix *m1, Matrix *m2) : Matrix(m1->Height(), m2->Height()) {
  this->m = (m1->Multiply(m2))->m;  // TODO(jrm) don't do this
  for (size_t row = 0; row < this->Height(); ++row) {
    for (size_t col = 0; col < this->Width(); ++col) {
      Vector* vec1 = m1->Row(row);
      Vector* vec2 = m2->Row(col);
      double elem = this->Get(row, col);
      elem = elem / sqrt(vec1->Multiply(vec1) * vec2->Multiply(vec2));
      this->Set(row, col, elem);
      delete vec2;
      delete vec1;
    }  // for col
  }  // for row
}

Kernel::~Kernel() {
}
}
