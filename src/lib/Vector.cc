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
  return new Matrix(mat);  // TODO(jrm) warning! newing up
}

Matrix *Vector::RepmatHoriz(size_t k) {
  gsl_matrix *mat = gsl_matrix_alloc(k, v->size);
  for (size_t i = 0; i < k; ++i) {
    gsl_matrix_set_row(mat, i, v);
  }
  return new Matrix(mat);  // TODO(jrm) warning! newing up
}

double Vector::Multiply(Vector *other) {
  double result;
  gsl_blas_ddot(this->v, other->v, &result);
  return result;
}

Vector *Vector::Multiply(Matrix *m) {
  gsl_vector *result = gsl_vector_alloc(m->Width());
  cblas_dgemm(CblasRowMajor,  // const enum CBLAS_ORDER Order
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransA
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransB
      1,                      // const int M
      m->Width(),             // const int N
      this->Size(),           // const int K
      1.0f,                   // const double alpha
      this->v->data,          // const double * A
      this->Size(),           // const int lda
      m->m->data,             // const double * B
      m->Width(),             // const int ldb
      0.0f,                   // const double beta
      result->data,           // double * C
      result->size);          // const int ldc
  return new Vector(result->data, m->Width());
  // TODO(jrm) warning! newing up
}

Vector *Vector::Subtract(Vector *other) {
  Vector *ret = new Vector(v);  // TODO(jrm) warning! newing up
  gsl_vector_sub(ret->v, other->v);
  return ret;
}

Vector *Vector::Add(Vector *other) {
  Vector *ret = new Vector(v);  // TODO(jrm) warning! newing up
  gsl_vector_add(ret->v, other->v);
  return ret;
}

Vector::Vector(gsl_vector *v) {
  this->v = v;
}
}
