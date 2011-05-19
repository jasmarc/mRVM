// Copyright 2011 Jason Marcell

#include <ctype.h>

#include <gsl/gsl_statistics.h>

#include <cstring>

#include "lib/Matrix.h"
#include "lib/Log.h"

namespace jason {

Matrix::Matrix(size_t height, size_t width) {
  LOG(DEBUG, "Matrix Constructor with height %zu and width %zu.\n",
    height, width);
  Init();
  this->m = gsl_matrix_alloc(height, width);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
}

Matrix::Matrix(gsl_matrix *mat) {
  LOG(DEBUG, "Matrix Constructor with mat.\n");
  Init();
  this->m = mat;
}

Matrix::Matrix(Vector *vec) {
  LOG(DEBUG, "Matrix Constructor with vector.\n");
  Init();
  gsl_matrix *mat = gsl_matrix_alloc(vec->v->size, vec->v->size);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  gsl_vector_view diag = gsl_matrix_diagonal(mat);
  gsl_matrix_set_all(mat, 0.0);
  gsl_vector_memcpy(&diag.vector, vec->v);
  this->m = mat;
}

Matrix::Matrix(double *data, size_t height, size_t width) {
  LOG(DEBUG, "Matrix Constructor with data, height, width.\n");
  Init();
  this->m = gsl_matrix_alloc(height, width);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  for (size_t row = 0; row < height; ++row) {
    for (size_t col = 0; col < width; ++col) {
      gsl_matrix_set(m, row, col, *(data + row * width + col));
    }
  }
}

Matrix::Matrix(const char* filename) {  // TODO(jrm): move to another class
  LOG(DEBUG, "Matrix Constructor with filename %s.\n", filename);
  Init();
  int rows, cols;
  FILE *f;
  f = fopen(filename, "r");
  if (f) {
    rows = NumberOfRows(f);
    cols = NumberOfColumns(f);
    LOG(DEBUG, "rows: %d, cols: %d\n", rows, cols);
    this->m = gsl_matrix_alloc(rows, cols);
    LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
    if (gsl_matrix_fscanf(f, this->m) == GSL_EFAILED) {
      perror("Error");
      throw("File read error.");
    }
    fclose(f);
  } else {
    perror("Error");
    throw("File read error.");
  }
}

void Matrix::Init() {
  LOG(DEBUG, "Matrix Init\n");
  this->to_str = reinterpret_cast<char*>(malloc(256 * sizeof(*to_str)));
  means = NULL;
  stdevs = NULL;
}

Matrix::~Matrix() {
  LOG(DEBUG, "Matrix Destructor.\n");
  free(this->to_str);
  gsl_matrix_free(this->m);
  LOG(DEBUG, "\t\t\tgsl_matrix_free\n");
  if (means != NULL) {
    delete means;
    delete stdevs;
  }
}

void Matrix::Write(const char* filename) {  // TODO(jrm): move to another class
  LOG(DEBUG, "Matrix Write.\n");
  FILE *f;
  f = fopen(filename, "w");
  if (f) {
    for (size_t row = 0; row < this->Height(); ++row) {
      for (size_t col = 0; col < this->Width(); ++col) {
        fprintf(f, "%.3f ", this->Get(row, col));
      }
      fprintf(f, "\n");
    }
    fclose(f);
  } else {
    perror("Error");
    throw("File write error.");
  }
}

Matrix* Matrix::RemoveRowsReturnMatrix(Vector *rows) {
  LOG(DEBUG, "RemoveRows (returns new Matrix).\n");
  size_t new_height = 0;
  for (size_t row = 0; row < rows->Size(); ++row) {
    new_height += rows->Get(row);
  }
  LOG(DEBUG, "New height = %zu.\n", new_height);
  gsl_matrix *new_m = gsl_matrix_alloc(new_height, this->Width());
  for (size_t ret_row = 0, row = 0; row < this->Height(); ++row) {
    if (rows->Get(row) == 1) {
      gsl_vector *v = gsl_vector_alloc(this->Width());
      gsl_matrix_get_row(v, m, row);
      gsl_matrix_set_row(new_m, ret_row++, v);
      gsl_vector_free(v);
    }
  }
  return new Matrix(new_m);
}

void Matrix::RemoveRows(Vector *rows) {
  LOG(DEBUG, "RemoveRows.\n");
  size_t new_height = 0;
  for (size_t row = 0; row < rows->Size(); ++row) {
    new_height += rows->Get(row);
  }
  LOG(DEBUG, "New height = %zu.\n", new_height);
  gsl_matrix *new_m = gsl_matrix_alloc(new_height, this->Width());
  for (size_t ret_row = 0, row = 0; row < this->Height(); ++row) {
    if (rows->Get(row) == 1) {
      gsl_vector *v = gsl_vector_alloc(this->Width());
      gsl_matrix_get_row(v, m, row);
      gsl_matrix_set_row(new_m, ret_row++, v);
      gsl_vector_free(v);
    }
  }
  gsl_matrix_free(this->m);
  this->m = new_m;
}

void Matrix::RemoveColumns(Vector *columns) {
  LOG(DEBUG, "RemoveColumns.\n");
  size_t new_width = 0;
  for (size_t column = 0; column < columns->Size(); ++column) {
    new_width += columns->Get(column);
  }
  LOG(DEBUG, "New width = %zu.\n", new_width);
  gsl_matrix *new_m = gsl_matrix_alloc(this->Height(), new_width);
  for (size_t ret_column = 0, column = 0; column < this->Width(); ++column) {
    if (columns->Get(column) == 1) {
      gsl_vector *v = gsl_vector_alloc(this->Height());
      gsl_matrix_get_col(v, m, column);
      gsl_matrix_set_col(new_m, ret_column++, v);
      gsl_vector_free(v);
    }
  }
  gsl_matrix_free(this->m);
  this->m = new_m;
}

void Matrix::Invert() {
  int n = this->Width();
  gsl_matrix *inverse = gsl_matrix_alloc(n, n);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  gsl_permutation *perm = gsl_permutation_alloc(n);
  int s = 0;
  gsl_linalg_LU_decomp(m, perm, &s);
  gsl_linalg_LU_invert(m, perm, inverse);
  gsl_permutation_free(perm);
  gsl_matrix_free(m);
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  m = inverse;
}

double Matrix::Get(int row, int col) {
  return gsl_matrix_get(this->m, row, col);
}

void Matrix::Set(int row, int col, double val) {
  gsl_matrix_set(this->m, row, col, val);
}

void Matrix::Add(Matrix *other) {
  LOG(DEBUG, "Adding a %zux%zu to a %zux%zu.\n",
    this->Height(), this->Width(), other->Height(), other->Width());
  if ((this->Width() != other->Width())
    || (this->Height() != other->Height())) {
    fprintf(stderr, "Dimension Error.\n");
    exit(1);
  }
  gsl_matrix_add(this->m, other->m);
}

Vector* Matrix::Row(size_t row) {
  gsl_vector *v = gsl_vector_alloc(this->Width());
  gsl_matrix_get_row(v, m, row);
  return new Vector(v);
}

Vector* Matrix::Column(size_t col) {
  gsl_vector *v = gsl_vector_alloc(this->Height());
  gsl_matrix_get_col(v, m, col);
  return new Vector(v);
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

char* Matrix::ToString() {
  LOG(DEBUG, "Matrix ToString\n");
  #define kElementSize 12
  this->to_str[0] = NULL;
  size_t total = 0;
  for (size_t row = 0; row < this->Height(); ++row) {
    for (size_t col = 0; col < this->Width(); ++col) {
      char temp[kElementSize];
      snprintf(temp, kElementSize, "%.3f\t", this->Get(row, col));
      strncat(this->to_str, temp, kElementSize);
      total += strlen(temp);
      if (total + kElementSize + 1 > 255) break;
    }
    char temp[1];
    snprintf(temp, kElementSize, "\n");
    strncat(this->to_str, temp, 1);
    total += 1;
    if (total + kElementSize > 255) break;
  }
  return this->to_str;
}

size_t Matrix::Height() {
  return this->m->size1;
}

size_t Matrix::Width() {
  return this->m->size2;
}

Matrix* Matrix::Multiply(Matrix *other) {
  LOG(DEBUG, "Multiplying a %zux%zu by a %zux%zu (transposed).\n",
    this->Height(), this->Width(), other->Height(), other->Width());
  if (this->Width() != other->Width()) {
    fprintf(stderr, "Dimension Error.\n");
    exit(1);
  }
  gsl_matrix *result = gsl_matrix_alloc(this->Height(), other->Height());
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  cblas_dgemm(CblasRowMajor,  // const enum CBLAS_ORDER Order
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransA
      CblasTrans,             // const enum CBLAS_TRANSPOSE TransB
      this->Height(),         // const int M
      other->Height(),        // const int N
      other->Width(),         // const int K
      1.0f,                   // const double alpha
      this->m->data,          // const double * A
      this->Width(),          // const int lda
      other->m->data,         // const double * B
      other->Width(),         // const int ldb
      0.0f,                   // const double beta
      result->data,           // double * C
      other->Height());       // const int ldc
  return new Matrix(result);
}

Matrix* Matrix::MultiplyNoTrans(Matrix *other) {
  LOG(DEBUG, "Multiplying a %zux%zu by a %zux%zu.\n",
    this->Height(), this->Width(), other->Height(), other->Width());
  if (this->Width() != other->Height()) {
    fprintf(stderr, "Dimension Error.\n");
    exit(1);
  }
  gsl_matrix *result = gsl_matrix_alloc(this->Height(), other->Width());
  LOG(DEBUG, "\t\t\tgsl_matrix_alloc\n");
  cblas_dgemm(CblasRowMajor,  // const enum CBLAS_ORDER Order
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransA
      CblasNoTrans,             // const enum CBLAS_TRANSPOSE TransB
      this->Height(),         // const int M
      other->Width(),         // const int N
      other->Height(),        // const int K
      1.0f,                   // const double alpha
      this->m->data,          // const double * A
      this->Width(),          // const int lda
      other->m->data,         // const double * B
      other->Width(),         // const int ldb
      0.0f,                   // const double beta
      result->data,           // double * C
      other->Width());        // const int ldc
  return new Matrix(result);
}

Vector* Matrix::Multiply(Vector *vec) {
  LOG(DEBUG, "Multiplying a %zux%zu by a vector of length %zu.\n",
    this->Height(), this->Width(), vec->Size());
  if (this->Width() != vec->Size()) {
    fprintf(stderr, "Dimension Error.\n");
    exit(1);
  }
  gsl_vector *result = gsl_vector_alloc(this->Height());
  cblas_dgemv(CblasRowMajor,  // const enum CBLAS_ORDER Order
      CblasNoTrans,           // const enum CBLAS_TRANSPOSE TransA
      this->Height(),         // const int M (height of A)
      this->Width(),          // const int N (width of A)
      1.0f,                   // const double alpha
      this->m->data,          // const double * A
      this->Width(),          // const int lda
      vec->v->data,           // const double * x
      1,                      // const int incx
      0.0f,                   // const double beta
      result->data,           // double * y
      1);                     // const int incy
  return new Vector(result);
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
