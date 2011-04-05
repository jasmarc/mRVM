// Copyright 2011 Jason Marcell

#include "lib/Predictor.h"
#include "lib/Kernel.h"
#include "lib/LinearKernel.h"
#include "lib/GaussHermiteQuadrature.h"
#include "lib/RandomNumberGenerator.h"
#include "lib/Log.h"

namespace jason {

Predictor::Predictor(Matrix *w, Matrix *x_train, Matrix *x_predict) {
  LOG(DEBUG, "== Beginning Predictor Constructor. ==\n\n");

  LOG(VERBOSE, "= printing x_predict before sphere: =\n");
  LOG(VERBOSE, "%s", x_predict->ToString());
  LOG(VERBOSE, "= end x_predict before sphere. =\n\n");

  x_predict->Sphere(x_train);

  LOG(DEBUG, "= printing x_predict after sphere: =\n");
  LOG(DEBUG, "%s", x_predict->ToString());
  LOG(DEBUG, "= end x_predict after sphere. =\n\n");

  this->k = new LinearKernel(x_train, x_predict);
  k->Init();

  LOG(VERBOSE, "= printing k: =\n");
  LOG(VERBOSE, "%s", k->ToString());
  LOG(VERBOSE, "= end k. =\n\n");

  this->w = w;

  LOG(DEBUG, "= printing w: =\n");
  LOG(DEBUG, "%s", w->ToString());
  LOG(DEBUG, "= end w. =\n\n");

  LOG(DEBUG, "== End Predictor Constructor. ==\n\n");
}

Predictor::~Predictor() {
}

void Predictor::Predict() {
  Matrix *predictions = QuadratureApproximation();
  predictions->NormalizeResults();
  LOG(NORMAL, "= predictions: =\n");
  LOG(NORMAL, "%s", predictions->ToString());
  LOG(NORMAL, "= end predictions. =\n\n");
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
            prod *= r->GaussianCDF(points[k] + wikn - wjkn);
          }  // if
        }  // for j
        sum += weights[k]*prod;
      }  // for k
      LOG(DEBUG, "sample n=%zu, class i=%zu, value=%f\n", n, i, sum);
      result->Set(n, i, sum);
    }  // for i
  }  // for n
  delete r;
  return result;
}
}
