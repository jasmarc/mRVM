// Copyright 2011 Jason Marcell

#include "lib/Vector.h"

namespace jason {

Vector::Vector(double *data, size_t size) {
  this->v = gsl_vector_alloc(size);
  for (size_t i = 0; i < size; ++i) {
    gsl_vector_set(this->v, i, data[i]);
  }
}

void Vector::Print() {
  for (size_t j = 0; j < v->size; j++) {
    printf("%.2f ", gsl_vector_get(v, j));
  }
  printf("\n");
}

Vector::~Vector() {
  gsl_vector_free(this->v);
}

size_t Vector::Size() {
  return this->v->size;
}

double Vector::Get(size_t elem) {
  return gsl_vector_get(this->v, elem);
}

Matrix *Vector::RepmatVert(size_t k) {
  gsl_matrix *mat = gsl_matrix_alloc(v->size, k);
  for (size_t i = 0; i < k; ++i) {
    gsl_matrix_set_col(mat, i, v);
  }
  return new Matrix(mat);
}

Matrix *Vector::RepmatHoriz(size_t k) {
  gsl_matrix *mat = gsl_matrix_alloc(k, v->size);
  for (size_t i = 0; i < k; ++i) {
    gsl_matrix_set_row(mat, i, v);
  }
  return new Matrix(mat);
}

double Vector::Multiply(Vector *other) {
  double result;
  gsl_blas_ddot(this->v, other->v, &result);
  return result;
}

Vector::Vector(gsl_vector *v) {
  this->v = v;
}
}
