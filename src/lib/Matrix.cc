// Copyright 2011 Jason Marcell

#include <ctype.h>

#include <gsl/gsl_statistics.h>

#include "lib/Matrix.h"

namespace jason {

Matrix::Matrix() {
  // TODO(jrm) danger, the m variable isn't getting set
}

Matrix::Matrix(size_t height, size_t width) {
  this->m = gsl_matrix_alloc(height, width);
}

Matrix::Matrix(gsl_matrix *mat) {
  this->m = mat;
}

Matrix::Matrix(Vector *vec) {
  gsl_matrix *mat = gsl_matrix_alloc(vec->v->size, vec->v->size);
  gsl_vector_view diag = gsl_matrix_diagonal(mat);
  gsl_matrix_set_all(mat, 0.0);
  gsl_vector_memcpy(&diag.vector, vec->v);
  this->m = mat;
}

Matrix::Matrix(double *data, size_t height, size_t width) {
  this->m = gsl_matrix_alloc(height, width);
  for (size_t row = 0; row < height; ++row) {
    for (size_t col = 0; col < width; ++col) {
      gsl_matrix_set(m, row, col, *(data + row * width + col));
    }
  }
}

Matrix::Matrix(const char* filename) {
  int rows, cols;
  FILE *f;
  f = fopen(filename, "r");
  if (f) {
    rows = NumberOfRows(f);
    cols = NumberOfColumns(f);
    this->m = gsl_matrix_alloc(rows, cols);
    gsl_matrix_fscanf(f, this->m);
    fclose(f);
  } else {
    perror("Error");
    throw("File read error.");
  }
}

Matrix::~Matrix() {
  gsl_matrix_free(this->m);
}

Matrix *Matrix::Clone() {  // TODO(jrm) change to copy constructor
  return new Matrix(this->m->data, this->Width(), this->Height());
}

void Matrix::Invert() {
  int n = this->Width();
  gsl_matrix *inverse = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);
  int s = 0;
  gsl_linalg_LU_decomp(m, perm, &s);
  gsl_linalg_LU_invert(m, perm, inverse);
  gsl_permutation_free(perm);
  gsl_matrix_free(m);
  m = inverse;
}

double Matrix::Get(int row, int col) {
  return gsl_matrix_get(this->m, row, col);
}

void Matrix::Set(int row, int col, double val) {
  gsl_matrix_set(this->m, row, col, val);
}

void Matrix::Add(Matrix *other) {
  gsl_matrix_add(this->m, other->m);
}

Vector* Matrix::Row(size_t row) {
  gsl_vector *v = gsl_vector_alloc(this->Width());
  gsl_matrix_get_row(v, m, row);
  return new Vector(v);  // TODO(jrm) warning! newing up
}

Vector* Matrix::Column(size_t col) {
  gsl_vector *v = gsl_vector_alloc(this->Height());
  gsl_matrix_get_col(v, m, col);
  return new Vector(v);  // TODO(jrm) warning! newing up
}

void Matrix::SetRow(size_t row, Vector *vec) {
  gsl_matrix_set_row(m, row, vec->v);
}

void Matrix::SetColumn(size_t col, Vector *vec) {
  gsl_matrix_set_col(m, col, vec->v);
}

//  TODO(jrm): fix for general purpose sphering
void Matrix::Sphere() {
  Sphere(this);
}

void Matrix::Sphere(Matrix *other) {
  size_t height = this->Height();
  size_t width = this->Width();
  Vector *vec;
  for (size_t col = 0; col < width; ++col) {
    vec = other->Column(col);
    double mean = gsl_stats_mean(vec->v->data, 1, height);
    double stdev = gsl_stats_sd(vec->v->data, 1, height);
    gsl_vector_add_constant(vec->v, -mean);
    gsl_vector_scale(vec->v, 1.0 / stdev);
    gsl_matrix_set_col(m, col, vec->v);
    delete vec;
  }
}

void Matrix::Print() {  // TODO(jrm): turn into ToString
  for (size_t row = 0; row < this->Height(); ++row) {
    for (size_t col = 0; col < this->Width(); ++col) {
      printf("%.2f\t", this->Get(row, col));
    }
    printf("\n");
  }
}

char * Matrix::ToString() {
  char *ret = (char*) malloc(256 * sizeof(char));  // TODO(jrm): this should get freed
  ret[0] = NULL;
  for (size_t row = 0; row < this->Height(); ++row) {
    for (size_t col = 0; col < this->Width(); ++col) {
      snprintf(ret, 256 * sizeof(char), "%s%.2f\t", ret, this->Get(row, col));
    }
    snprintf(ret, 256 * sizeof(char), "%s\n", ret);
  }
  return ret;
}

size_t Matrix::Height() {
  return this->m->size1;
}

size_t Matrix::Width() {
  return this->m->size2;
}

Matrix* Matrix::Multiply(Matrix *other) {
  gsl_matrix *result = gsl_matrix_alloc(this->Height(), other->Height());
  cblas_dgemm(CblasRowMajor,  // const enum CBLAS_ORDER Order
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransA
      CblasTrans,             // const enum CBLAS_TRANSPOSE TransB
      this->Height(),         // const int M
      other->Height(),        // const int N
      this->Width(),          // const int K
      1.0f,                   // const double alpha
      this->m->data,          // const double * A
      this->Width(),          // const int lda
      other->m->data,         // const double * B
      other->Width(),         // const int ldb
      0.0f,                   // const double beta
      result->data,           // double * C
      result->size2);         // const int ldc
  return new Matrix(result->data, this->Height(), other->Height());
  // TODO(jrm) warning! newing up
}

Vector* Matrix::Multiply(Vector *vec) {
  gsl_vector *result = gsl_vector_calloc(vec->Size());
  cblas_dgemm(CblasRowMajor,  // const enum CBLAS_ORDER Order
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransA
      CblasTrans,             // const enum CBLAS_TRANSPOSE TransB
      this->Height(),         // const int M (height of A)
      vec->Size(),            // const int N (width of B)
      this->Width(),          // const int K (width of A)
      1.0f,                   // const double alpha
      this->m->data,          // const double * A
      m->size2,               // const int lda
      vec->v->data,           // const double * B
      vec->Size(),            // const int ldb
      0.0f,                   // const double beta
      result->data,           // double * C
      1);                     // const int ldc
  return new Vector(result->data, result->size);
  // TODO(jrm) warning! newing up
}

int Matrix::NumberOfRows(FILE *f) {
  char lastChar = '\n';
  char currentChar = NULL;
  int count = 0;
  while ((currentChar = fgetc(f)) != EOF) {
    if (lastChar == '\n' && currentChar != '\n') {
      ++count;
    }
    lastChar = currentChar;
  }
  fseek(f, 0, SEEK_SET);
  return count;
}

int Matrix::NumberOfColumns(FILE *f) {
  char lastChar = ' ';
  char currentChar = NULL;
  int count = 0;
  while ((currentChar = fgetc(f)) != '\n') {
    if (isspace(lastChar) && !isspace(currentChar)) {
      ++count;
    }
    lastChar = currentChar;
  }
  fseek(f, 0, SEEK_SET);
  return count;
}

gsl_matrix *Matrix::CloneGSLMatrix() {
  size_t height = this->Height();
  size_t width = this->Width();
  gsl_matrix *ret = gsl_matrix_alloc(height, width);
  for (size_t row = 0; row < height; ++row) {
    for (size_t col = 0; col < width; ++col) {
      gsl_matrix_set(ret, row, col, *(this->m->data + row * width + col));
    }
  }
  return ret;
}
}
