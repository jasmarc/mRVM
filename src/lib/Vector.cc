// Copyright 2011 Jason Marcell

#include <cstring>

#include "lib/Vector.h"
#include "lib/Log.h"

namespace jason {

Vector::Vector(size_t size) {
  to_str = reinterpret_cast<char*> (malloc(256 * sizeof(*to_str)));
  // TODO(jrm): check if this is the same as gsl_vector_alloc
  this->v = gsl_vector_calloc(size);
}

Vector::Vector(gsl_vector *vec) {
  to_str = reinterpret_cast<char*> (malloc(256 * sizeof(*to_str)));
  this->v = vec;
}

Vector::Vector(double *data, size_t size) {
  to_str = reinterpret_cast<char*> (malloc(256 * sizeof(*to_str)));
  this->v = gsl_vector_alloc(size);
  for (size_t i = 0; i < size; ++i) {
    gsl_vector_set(this->v, i, data[i]);
  }
}

Vector::Vector(const char* filename) {
  to_str = reinterpret_cast<char*> (malloc(256 * sizeof(*to_str)));
  size_t rows;
  FILE *f;
  f = fopen(filename, "r");
  if (f) {
    rows = NumberOfElements(f);
    this->v = gsl_vector_alloc(rows);
    gsl_vector_fscanf(f, this->v);
    fclose(f);
  } else {
    perror("Error");
    throw("File read error.");
  }
}

char* Vector::ToString() {
  LOG(DEBUG, "Vector ToString\n");
  #define kElementSize 12
  this->to_str[0] = NULL;
  size_t total = 0;
  for (size_t j = 0; j < v->size; j++) {
    char temp[kElementSize];
    snprintf(temp, kElementSize, "%.3f\t", this->Get(j));
    strncat(this->to_str, temp, kElementSize);
    total += strlen(temp);
    if (total + kElementSize + 1 > 255) break;
  }
  char temp[1];
  snprintf(temp, kElementSize, "\n");
  strncat(this->to_str, temp, 1);
  return this->to_str;
}

Vector::~Vector() {
  free(to_str);
  gsl_vector_free(this->v);
}

size_t Vector::Size() {
  return this->v->size;
}

double Vector::Get(size_t elem) {
  return gsl_vector_get(this->v, elem);
}

void Vector::Set(size_t elem, double value) {
  gsl_vector_set(this->v, elem, value);
}

Matrix *Vector::RepmatVert(size_t k) {
  gsl_matrix *mat = gsl_matrix_alloc(v->size, k);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  for (size_t i = 0; i < k; ++i) {
    gsl_matrix_set_col(mat, i, v);
  }
  return new Matrix(mat);
}

Matrix *Vector::RepmatHoriz(size_t k) {
  gsl_matrix *mat = gsl_matrix_alloc(k, v->size);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  for (size_t i = 0; i < k; ++i) {
    gsl_matrix_set_row(mat, i, v);
  }
  return new Matrix(mat);
}

Vector *Vector::RemoveElementsReturnVector(Vector *rows) {
  LOG(DEBUG, "RemoveElementsReturnVector.\n");
  size_t new_height = 0;
  for (size_t row = 0; row < rows->Size(); ++row) {
    new_height += rows->Get(row);
  }
  LOG(DEBUG, "New height = %zu.\n", new_height);
  gsl_vector *new_v = gsl_vector_alloc(new_height);
  for (size_t ret_row = 0, row = 0; row < this->Size(); ++row) {
    if (rows->Get(row) == 1) {
      double val = this->Get(row);
      gsl_vector_set(new_v, ret_row++, val);
    }
  }
  return new Vector(new_v);
}

double Vector::Multiply(Vector *other) {
  LOG(DEBUG, "Multiplying two vectors: %zu %zu\n",
    this->Size(), other->Size());
  if (this->Size() != other->Size()) {
    fprintf(stderr, "Dimension Error.\n");
    exit(1);
  }
  double result;
  gsl_blas_ddot(this->v, other->v, &result);
  return result;
}

Vector *Vector::Multiply(Matrix *m) {
  LOG(DEBUG, "Multiplying a vector by a matrix: %zu %zux%zu\n",
    this->Size(), m->Height(), m->Width());
  if (this->Size() != m->Height()) {
    fprintf(stderr, "Dimension Error.\n");
    exit(1);
  }
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
  return new Vector(result);
  // TODO(jrm) warning! newing up
}

Vector *Vector::Subtract(Vector *other) {
  Vector *ret = new Vector(v->data, v->size);  // TODO(jrm) warning! newing up
  gsl_vector_sub(ret->v, other->v);
  return ret;
}

Vector *Vector::Add(Vector *other) {
  Vector *ret = new Vector(v->data, v->size);  // TODO(jrm) warning! newing up
  gsl_vector_add(ret->v, other->v);
  return ret;
}

size_t Vector::GetNumberOfClasses() {
  size_t max = 0;
  LOG(DEBUG, "Size = %zu\n", this->Size());
  for (double i = 0; i < this->Size(); ++i) {
    if (this->Get(i) > max)
      max = static_cast<int>(this->Get(i));
  }
  LOG(DEBUG, "Number of classes = %zu\n", max);
  return max + 1;
}

size_t Vector::NumberOfElements(FILE *f) {  // TODO(jrm): move to another class
  char lastChar = '\n';
  char currentChar = NULL;
  size_t count = 0;
  while ((currentChar = fgetc(f)) != EOF) {
    if (lastChar == '\n' && currentChar != '\n') {
      ++count;
    }
    lastChar = currentChar;
  }
  fseek(f, 0, SEEK_SET);
  return count;
}
}
