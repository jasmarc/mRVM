// Copyright 2011 Jason Marcell

#include "lib/Predictor.h"
#include "lib/Kernel.h"
#include "lib/LinearKernel.h"
#include "lib/GaussHermiteQuadrature.h"
#include "lib/RandomNumberGenerator.h"

namespace jason {

Predictor::Predictor(Matrix *w, Vector *t, Matrix *x_train, Matrix *x_predict) {
//  Vector *v;     // Labels
//  Matrix *prob;  // Class probabilities
  this->k = new LinearKernel(x_train, x_predict);
  this->w = w;
  this->t = t;
}

Predictor::~Predictor() {
}

void Predictor::Predict() {
  this->QuadratureApproximation();
}

void Predictor::QuadratureApproximation() {
  RandomNumberGenerator *r = new RandomNumberGenerator();
  GaussHermiteQuadrature *g = new GaussHermiteQuadrature();
  double *points;
  double *weights;
  g->Process(3, &points, &weights);
  delete g;

  for (size_t n = 0; n < w->Height(); ++n) {
    size_t i = (size_t)t->Get(n);
    Vector *kn = k->Row(n);
    Vector *wi = w->Column(i);
    double wikn = wi->Multiply(kn);
    double sum = 0;
    for (size_t k = 0; k < 3; ++k) {
      double prod = 1;
      for (size_t j = 0; j < w->Width(); ++j) {
        if (j != i) {
          Vector *wj = w->Column(j);
          double wjkn = wj->Multiply(kn);
          prod *= r->GaussianCDF(points[k] + wikn - wjkn);
        }  // if
      }  // for j
      sum += weights[k]*prod;
    }  // for k
  }  // for n
  delete r;
}
}
