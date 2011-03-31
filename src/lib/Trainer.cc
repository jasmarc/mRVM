// Copyright 2011 Jason Marcell

#include "lib/RandomNumberGenerator.h"
#include "lib/Trainer.h"
#include "lib/LinearKernel.h"

namespace jason {

Trainer::Trainer(Matrix *matrix, Vector *labels, size_t classes) {
  this->x = matrix;
  this->t = labels;
  this->samples = matrix->Height();
  this->features = matrix->Width();
  this->classes = classes;
}

Trainer::~Trainer() {
}

void Trainer::Process() {
//  TODO(jrm): Add proper debugging statements
  printf("printing x before sphere:\n");
  x->Print();
  printf("end x before sphere.\n");

  x->Sphere();

  printf("printing x after sphere:\n");
  x->Print();
  printf("end x after sphere.\n");

  k = BuildKernel(x);
  k->Init();

  printf("printing kernel:\n");
  k->Print();
  printf("end kernel.\n");

  InitializeYAW();
  for (int i = 0; i < 20; ++i) {
    UpdateW();
    UpdateA(1.0, 1.0);
    UpdateY();
  }
  w->Print();
}

Matrix *Trainer::GetW() {
  return this->w;
}

Kernel *Trainer::BuildKernel(Matrix *m) {
  return new LinearKernel(m, m);
}

void Trainer::InitializeYAW() {
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
      w_val = r->SampleGaussian(1/a_val);
      y->Set(row, col, y_val);
      a->Set(row, col, a_val);
      w->Set(row, col, w_val);
    }
  }
  delete r;
}

void Trainer::UpdateA(double tau, double v) {
  for (size_t row = 0; row < samples; ++row) {
    for (size_t col = 0; col < classes; ++col) {
      double val = w->Get(row, col);
      a->Set(row, col, (2*tau + 1)/(val*val + 2*v));
    }
  }
}

void Trainer::UpdateW() {
  Matrix *w_temp;
  for (size_t col = 0; col < classes; ++col) {
    Vector *A_c = a->Column(col);
    Matrix *A = new Matrix(A_c);
    Vector *Y_c = y->Column(col);
    w_temp = k->Multiply(k);
    w_temp->Add(A);
    w_temp->Invert();
    w_temp = w_temp->Multiply(k);
    Vector *W_c = w_temp->Multiply(Y_c);
    w->SetColumn(col, W_c);
    delete W_c;
    delete Y_c;
    delete A;
    delete A_c;
  }
}

void Trainer::UpdateY() {
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
  }  // for n
  delete r;
}
}
