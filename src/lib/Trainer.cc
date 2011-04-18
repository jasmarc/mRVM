// Copyright 2011 Jason Marcell

#include <math.h>

#include "lib/RandomNumberGenerator.h"
#include "lib/Trainer.h"
#include "lib/LinearKernel.h"
#include "lib/Log.h"

#define EPSILON 0.001
#define MAX_ITER 100
namespace jason {

Trainer::Trainer(Matrix *matrix, Vector *labels, size_t classes,
    Kernel *kernel) {
  this->x = matrix;
  this->t = labels;
  this->samples = matrix->Height();
  this->features = matrix->Width();
  this->classes = classes;
  this->k = kernel;
  this->converged = false;
}

Trainer::~Trainer() {
  delete y;
  delete a;
  delete w;
}

void Trainer::Process(double tau, double upsilon) {
  LOG(DEBUG, "== Beginning Trainer. ==\n\n");
  LOG(DEBUG, "= Initializing Train Kernel. =\n")

  k->Init();

  LOG(DEBUG, "= Printing Train Kernel: =\n");
  LOG(DEBUG, "%s\n", k->ToString());

  InitializeYAW();
  for (size_t i = 0; i < MAX_ITER && !converged; ++i) {
    LOG(DEBUG, "Iteration: %zu\n", i);
    UpdateW();
    UpdateA(tau, upsilon);
    UpdateY();
  }

  LOG(DEBUG, "= Printing w: =\n");
  LOG(DEBUG, "%s\n", w->ToString());

  LOG(DEBUG, "== End Trainer. ==\n");
}

Matrix *Trainer::GetW() {
  return this->w;
}

void Trainer::InitializeYAW() {
  LOG(DEBUG, "= InitializeYAW. =\n");
  y = new Matrix(samples, classes);
  a = new Matrix(samples, classes);
  w = new Matrix(samples, classes);
  RandomNumberGenerator *r = new RandomNumberGenerator();
  for (size_t row = 0; row < samples; ++row) {
    for (size_t col = 0; col < classes; ++col) {
      double y_val, a_val, w_val;
      if (t->Get(row) == col)
        y_val = r->SampleUniform(0, 10);
      else
        y_val = r->SampleUniform(0, 1);
      a_val = 1;
      w_val = r->SampleGaussian(sqrt(1/a_val));
      y->Set(row, col, y_val);
      a->Set(row, col, a_val);
      w->Set(row, col, w_val);
    }
  }
  delete r;
}

void Trainer::UpdateA(double tau, double upsilon) {
  LOG(DEBUG, "= UpdateA. =\n");
  this->converged = true;
  for (size_t row = 0; row < samples; ++row) {
    for (size_t col = 0; col < classes; ++col) {
      double wval = w->Get(row, col);
      double oldval = a->Get(row, col);
      double newval = (2*tau + 1)/(wval*wval + 2*upsilon);
      a->Set(row, col, newval);
      LOG(DEBUG, "UpdateA: %.3f\t%.3f\t%.3f\n", oldval, newval, fabs(oldval - newval));
      if (fabs(oldval - newval) > EPSILON) {
        this->converged = false;
      }
    }
  }
}

void Trainer::UpdateW() {
  LOG(DEBUG, "= UpdateW. =\n");
  for (size_t col = 0; col < classes; ++col) {
    Vector *A_c = a->Column(col);
    Matrix *A = new Matrix(A_c);
    Vector *Y_c = y->Column(col);
    Matrix *w_temp1 = k->Multiply(k);
    w_temp1->Add(A);
    w_temp1->Invert();
    Matrix *w_temp2 = w_temp1->Multiply(k);
    Vector *W_c = w_temp2->Multiply(Y_c);
    w->SetColumn(col, W_c);
    delete w_temp2;
    delete w_temp1;
    delete W_c;
    delete Y_c;
    delete A;
    delete A_c;
  }
}

void Trainer::UpdateY() {
  LOG(DEBUG, "= UpdateY. =\n");
  RandomNumberGenerator *r = new RandomNumberGenerator();
  for (size_t n = 0; n < samples; ++n) {
    size_t i = (size_t)t->Get(n);
    Vector *k_n = k->Row(n);
    for (size_t c = 0; c < classes; ++c) {
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
