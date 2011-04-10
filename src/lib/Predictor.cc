// Copyright 2011 Jason Marcell

#include "lib/Predictor.h"
#include "lib/Kernel.h"
#include "lib/LinearKernel.h"
#include "lib/GaussHermiteQuadrature.h"
#include "lib/RandomNumberGenerator.h"
#include "lib/Log.h"

namespace jason {

Predictor::Predictor(Matrix *w, Matrix *x_train, Matrix *x_predict,
    Kernel *kernel) {
  this->k = kernel;
  this->w = w;
}

Predictor::~Predictor() {
}

Matrix* Predictor::Predict() {
  LOG(VERBOSE, "= Initializing Predictor Kernel. =\n");

  k->Init();

  LOG(VERBOSE, "= Printing Predictor Kernel: =\n");
  LOG(VERBOSE, "%s\n", k->ToString());

  Matrix *predictions = QuadratureApproximation();
  predictions->NormalizeResults();
  return predictions;
}

Matrix* Predictor::QuadratureApproximation() {
  RandomNumberGenerator *r = new RandomNumberGenerator();
  GaussHermiteQuadrature *g = new GaussHermiteQuadrature();
  double *points;
  double *weights;
  g->Process(3, &points, &weights);
  delete g;

  Matrix *result = new Matrix(k->Width(), w->Width());
  for (size_t n = 0; n < k->Width(); ++n) {
    for (size_t i = 0; i < w->Width(); ++i) {
      Vector *kn = k->Column(n);
      Vector *wi = w->Column(i);
      double wikn = wi->Multiply(kn);
      double sum = 0;
      for (size_t k = 0; k < 3; ++k) {
        double prod = 1;
        for (size_t j = 0; j < w->Width(); ++j) {
          if (j != i) {
            Vector *wj = w->Column(j);
            double wjkn = wj->Multiply(kn);
            delete wj;
            prod *= r->GaussianCDF(points[k] + wikn - wjkn);
          }  // if
        }  // for j
        sum += weights[k]*prod;
      }  // for k
      delete wi;
      delete kn;
      LOG(DEBUG, "sample n=%zu, class i=%zu, value=%f\n", n, i, sum);
      result->Set(n, i, sum);
    }  // for i
  }  // for n
  delete r;
  delete[] points;
  delete[] weights;
  return result;
}
}
