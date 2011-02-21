// Copyright 2011 Jason Marcell

#include "lib/Predictor.h"
#include "lib/Kernel.h"

namespace jason {

Predictor::Predictor(Matrix *w, Matrix *x_train, Matrix *x_predict) {
  Vector *v;     // Labels
  Matrix *prob;  // Class probabilities
  Kernel *k = new Kernel(x_train, x_predict);
  k->Print();
}

Predictor::~Predictor() {
}
}
