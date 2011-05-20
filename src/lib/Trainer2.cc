// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/RandomNumberGenerator.h"
#include "lib/Trainer2.h"
#include "lib/LinearKernel.h"
#include "lib/Log.h"

#define INF 9999
#define EPSILON 0.001
#define MAX_ITER 100
namespace jason {

Trainer2::Trainer2(Matrix *matrix, Vector *labels, size_t classes,
    Kernel *kernel) {
  this->x = matrix;
  this->t = labels;
  this->samples = matrix->Height();
  this->features = matrix->Width();
  this->classes = classes;
  this->k = kernel;
  this->converged = false;
  this->active_samples = new Vector(samples);
}

Trainer2::~Trainer2() {
  delete y;
  delete a;
  delete w;
  delete active_samples;
}

void Trainer2::Process(double tau, double upsilon) {
  LOG(DEBUG, "== Beginning Trainer. ==\n\n");
  LOG(DEBUG, "= Initializing Train Kernel. =\n")

  k->Init();

  LOG(DEBUG, "= Printing Train Kernel: =\n");
  LOG(DEBUG, "%s\n", k->ToString());

  InitializeYAW();
  size_t current_sample_index = GetFirstSampleIndex();
  active_samples->Set(current_sample_index, 1);
  for (size_t i = 0; i < MAX_ITER && !converged; ++i) {
    LOG(DEBUG, "Iteration: %zu\n", i);
    Matrix *kstar = k->RemoveRowsReturnMatrix(active_samples);
    Matrix *astar = new Matrix(a->RemoveElementsReturnVector(active_samples));
    Matrix *kka_inv = kstar->Multiply(kstar);
    kka_inv->Add(astar);
    kka_inv->Invert();
    Vector *ki = k->Row(current_sample_index);
    double ai = a->Get(current_sample_index);
    bool is_already_in = (active_samples->Get(current_sample_index) == 1);
    bool does_contribute;
    double si;
    Vector *qci;
    if (i == 0) {
      double kiki = ki->Multiply(ki);
      Vector *kikstar = ki->Multiply(kstar);
      Vector *kstarki = kstar->Multiply(ki);
      Vector *kiy = ki->Multiply(y);
      Matrix *kstary = kstar->Multiply(y);
      Vector *temp1 = kikstar->Multiply(kka_inv);
      double temp2 = temp1->Multiply(kstarki);
      Vector *temp3 = temp1->Multiply(kstary);
      double sm = kiki - temp2;
      Vector *qcm = kiy->Subtract(temp3);
      if (!is_already_in) {
        si = sm;
        qci = qcm;
      } else {
        si = ai*sm/(ai - sm);
        qci = qcm;
        for (size_t i; i < qci->Size(); ++i) {
          qci->Set(i, ai*qci->Get(i) / (ai - sm));
        }
      }
      delete qcm;  // TODO(jrm): probably don't do this
      delete temp3;
      delete temp1;
      delete kstary;
      delete kiy;
      delete kstarki;
      delete kikstar;
    } else {
      // TODO(jrm): missing code
    }
    if (does_contribute && is_already_in) {
      // TODO(jrm): missing code
    } else if (does_contribute && !is_already_in) {
      // TODO(jrm): missing code
    } else if (!does_contribute && is_already_in) {
      // if (length(active_samples) != 1) {
      // TODO(jrm): missing code
      // }
    }
    delete qci;  // TODO(jrm): do this?
    delete ki;
    delete kka_inv;
    delete astar;
    delete kstar;
  }

  LOG(DEBUG, "= Printing w: =\n");
  LOG(DEBUG, "%s\n", w->ToString());

  LOG(DEBUG, "== End Trainer. ==\n");
}

Matrix *Trainer2::GetW() {
  return this->w;
}

void Trainer2::InitializeYAW() {
  LOG(DEBUG, "= InitializeYAW. =\n");
  y = new Matrix(samples, classes);
  a = new Vector(samples);
  w = new Matrix(samples, classes);
  RandomNumberGenerator *r = new RandomNumberGenerator();
  for (size_t row = 0; row < samples; ++row) {
    a->Set(row, INF);
    for (size_t col = 0; col < classes; ++col) {
      double y_val, a_val, w_val;
      if (t->Get(row) == col)
        y_val = r->SampleUniform(0, 10);
      else
        y_val = r->SampleUniform(0, 1);
      a_val = a->Get(row);
      w_val = r->SampleGaussian(sqrt(1/a_val));
      y->Set(row, col, y_val);
      w->Set(row, col, w_val);
    }
  }
  delete r;
}

size_t Trainer2::GetFirstSampleIndex() {
  Vector *v = new Vector(k->Height());
  Matrix *result = k->Multiply(y);
  for (size_t row = 0; row < result->Height(); ++row) {
    double num = 0.0;
    double den = 0.0;
    for (size_t col = 0; col < result->Width(); ++col) {
      num += pow(result->Get(row, col), 2.0);
      den += pow(k->Get(row, col), 2.0);
    }
    v->Set(row, (num / (classes * den)));
  }
  delete result;
  size_t max_index = 0;
  double max = -1.0;
  for (size_t i; i < v->Size(); ++i) {
    if (v->Get(i) > max) {
      max = v->Get(i);
      max_index = i;
    }
  }
  delete v;
  return max_index;
}

void Trainer2::UpdateA(double tau, double upsilon) {
  LOG(DEBUG, "= UpdateA. =\n");
  this->converged = true;
  // TODO(jrm): Serious re-work
  // result = (C*si^2) / (sum(qci.^2)-C*si);
  // if result<1e-5
  //   result = 1e-6;
  // end
}

void Trainer2::UpdateW() {
  LOG(DEBUG, "= UpdateW. =\n");
  LOG(DEBUG, "k is %zux%zu\n", k->Height(), k->Width());
  LOG(DEBUG, "y is %zux%zu\n", y->Height(), y->Width());
  LOG(DEBUG, "a is 1x%zu\n", a->Size());
  LOG(DEBUG, "w is %zux%zu\n", w->Height(), w->Width());
  // TODO(jrm): Serious re-work
  // W(active_samples,:) = KKA_inv * Kstar * Y';
}

void Trainer2::UpdateY() {
  LOG(DEBUG, "= UpdateY. =\n");
  LOG(DEBUG, "k is %zux%zu\n", k->Height(), k->Width());
  LOG(DEBUG, "y is %zux%zu\n", y->Height(), y->Width());
  LOG(DEBUG, "a is 1x%zu\n", a->Size());
  LOG(DEBUG, "w is %zux%zu\n", w->Height(), w->Width());
  RandomNumberGenerator *r = new RandomNumberGenerator();
  for (size_t n = 0; n < samples; ++n) {
    LOG(DEBUG, "n = %zu.\n", n);
    size_t i = (size_t)t->Get(n);
    Vector *k_n = k->Column(n);
    for (size_t c = 0; c < classes; ++c) {
      LOG(DEBUG, "c = %zu.\n", c);
      Vector *w_c = w->Column(c);
      Vector *w_i = w->Column(i);
      double wckn = w_c->Multiply(k_n);
      double wikn = w_i->Multiply(k_n);
      delete w_c;
      delete w_i;
      if (c != i) {
        double numerator = 0;
        double denominator = 0;
        for (int monte = 0; monte < 1000; ++monte) {
          double u = r->SampleGaussian(1.0);
          double num = r->GaussianPDF(wckn - wikn);
          double den = r->GaussianCDF(u + wikn - wckn);
          for (size_t j = 0; j < classes; ++j) {
            if (j != i && j != c) {
              LOG(DEBUG, "j = %zu.\n", j);
              Vector *w_j = w->Column(j);
              double wjkn = w_j->Multiply(k_n);
              num *= r->GaussianCDF(u + wikn - wjkn);
              den *= r->GaussianCDF(u + wikn - wjkn);
              delete w_j;
            }  // if
          }  // for j
          numerator   += num;
          denominator += den;
        }  // for monte
        if (denominator != 0) {
          double val = wckn - numerator / denominator;
          y->Set(n, c, val);
        } else {
          perror("Error! denominator equal to zero");
        }  // if
      } else {
        double val = wckn;
        for (size_t j = 0; j < classes; ++j) {
          if (j != i) {
            Vector *w_j = w->Column(j);
            double wjkn = w_j->Multiply(k_n);
            double y_nj = y->Get(n, j);
            val = val - (y_nj - wjkn);
            delete w_j;
          }  // if
        }  // for j
        y->Set(n, c, val);
      }  // if
    }  // for c
    delete k_n;
  }  // for n
  delete r;
}
}
