// Copyright 2011 Jason Marcell

#include <ctype.h>

#include <gsl/gsl_statistics.h>

#include "lib/Matrix.h"
#include "lib/Log.h"

namespace jason {

Matrix::Matrix(size_t height, size_t width) {
  to_str = reinterpret_cast<char*>(malloc(256 * sizeof(*to_str)));
  this->m = gsl_matrix_alloc(height, width);
}

Matrix::Matrix(gsl_matrix *mat) {
  to_str = reinterpret_cast<char*>(malloc(256 * sizeof(*to_str)));
  this->m = mat;
}

Matrix::Matrix(Vector *vec) {
  to_str = reinterpret_cast<char*>(malloc(256 * sizeof(*to_str)));
  gsl_matrix *mat = gsl_matrix_alloc(vec->v->size, vec->v->size);
  gsl_vector_view diag = gsl_matrix_diagonal(mat);
  gsl_matrix_set_all(mat, 0.0);
  gsl_vector_memcpy(&diag.vector, vec->v);
  this->m = mat;
}

Matrix::Matrix(double *data, size_t height, size_t width) {
  to_str = reinterpret_cast<char*>(malloc(256 * sizeof(*to_str)));
  this->m = gsl_matrix_alloc(height, width);
  for (size_t row = 0; row < height; ++row) {
    for (size_t col = 0; col < width; ++col) {
      gsl_matrix_set(m, row, col, *(data + row * width + col));
    }
  }
}

Matrix::Matrix(const char* filename) {  // TODO(jrm): move to another class
  to_str = reinterpret_cast<char*>(malloc(256 * sizeof(*to_str)));
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
  free(to_str);
  gsl_matrix_free(this->m);
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
  return new Vector(v->data, v->size);
}

Vector* Matrix::Column(size_t col) {
  gsl_vector *v = gsl_vector_alloc(this->Height());
  gsl_matrix_get_col(v, m, col);
  return new Vector(v->data, v->size);
}

void Matrix::SetRow(size_t row, Vector *vec) {
  gsl_matrix_set_row(m, row, vec->v);
}

void Matrix::SetColumn(size_t col, Vector *vec) {
  gsl_matrix_set_col(m, col, vec->v);
}

void Matrix::Sphere() {
  LOG(DEBUG, "= Calling sphere with no params =\n");
  Sphere(this);
}

void Matrix::Sphere(Matrix *other) {
  size_t width = this->Width();
  Vector *vec_other;
  Vector *vec_self;
  LOG(DEBUG, "= Start sphere =\n");
  LOG(DEBUG, "= this =\n");
  LOG(DEBUG, "%s", this->ToString());
  LOG(DEBUG, "= other =\n");
  LOG(DEBUG, "%s", other->ToString());
  LOG(DEBUG, "= starting loop =\n");
  for (size_t col = 0; col < width; ++col) {
    vec_other = other->Column(col);
    vec_self = this->Column(col);
    double mean = other->GetMeans()->Get(col);
    double stdev = other->GetStdevs()->Get(col);
    LOG(DEBUG, "col = %zu, mean = %f, stdev = %f\n", col, mean, stdev);
    gsl_vector_add_constant(vec_self->v, -mean);
    gsl_vector_scale(vec_self->v, 1.0 / stdev);
    gsl_matrix_set_col(m, col, vec_self->v);
    delete vec_other;
    delete vec_self;
  }
  LOG(DEBUG, "= this =\n");
  LOG(DEBUG, "%s", this->ToString());
  LOG(DEBUG, "= End sphere =\n");
}

void Matrix::NormalizeResults() {
  LOG(DEBUG, "= Start normalize results. =\n");
  for (size_t row = 0; row < this->Height(); ++row) {
    float sum = 0;
    for (size_t col = 0; col < this->Width(); ++col) {
      sum += this->Get(row, col);
    }
    for (size_t col = 0; col < this->Width(); ++col) {
      float val = this->Get(row, col);
      this->Set(row, col, val / sum);
    }
  }
  LOG(DEBUG, "= End normalize results. =\n");
}

void Matrix::CacheMeansAndStdevs() {
  size_t width = this->Width();
  size_t height = this->Height();
  means = new Vector(width);
  stdevs = new Vector(width);
  for (size_t col = 0; col < width; ++col) {
    Vector *vec = this->Column(col);
    double mean = gsl_stats_mean(vec->v->data, 1, height);
    // gsl_vector_add_constant(vec->v, -mean);
    double stdev = gsl_stats_sd(vec->v->data, 1, height);
    means->Set(col, mean);
    stdevs->Set(col, stdev);
    delete vec;
  }
}

Vector* Matrix::GetMeans() {
  return means;
}

Vector* Matrix::GetStdevs() {
  return stdevs;
}

char * Matrix::ToString() {
  to_str[0] = NULL;
  for (size_t row = 0; row < this->Height(); ++row) {
    for (size_t col = 0; col < this->Width(); ++col) {
      snprintf(to_str, 256 * sizeof(*to_str), "%s%.3f\t", to_str, this->Get(row, col));
    }
    snprintf(to_str, 256 * sizeof(*to_str), "%s\n", to_str);
  }
  return to_str;
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

size_t Matrix::NumberOfRows(FILE *f) {  // TODO(jrm): move to another class
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

size_t Matrix::NumberOfColumns(FILE *f) {  // TODO(jrm): move to another class
  char lastChar = ' ';
  char currentChar = NULL;
  size_t count = 0;
  while ((currentChar = fgetc(f)) != '\n') {
    if (isspace(lastChar) && !isspace(currentChar)) {
      ++count;
    }
    lastChar = currentChar;
  }
  fseek(f, 0, SEEK_SET);
  return count;
}
}
