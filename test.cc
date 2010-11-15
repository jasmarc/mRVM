// Copyright 2010 Jason Marcell

#include <gtest/gtest.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl_blas.h>
#include <stdio.h>
#include <stdarg.h>
#include <list>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <vector>

#define FOREACH(m_itname, m_container) \
        for (typeof(m_container.begin()) m_itname = m_container.begin(); \
        m_itname != m_container.end(); \
        m_itname++ )

using std::set;
using std::vector;

int NumberOfRows(FILE *f) {
  char lastChar = '\n';
  char currentChar = NULL;
  int count = 0;
  while ((currentChar=fgetc(f)) != EOF) {
    if (lastChar == '\n' && currentChar != '\n') {
      ++count;
    }
    lastChar = currentChar;
  }
  fseek(f, 0, SEEK_SET);
  return count;
}

int NumberOfColumns(FILE *f) {
  char lastChar = ' ';
  char currentChar = NULL;
  int count = 0;
  while ((currentChar=fgetc(f)) != '\n') {
    if (isspace(lastChar) && !isspace(currentChar)) {
      ++count;
    }
    lastChar = currentChar;
  }
  fseek(f, 0, SEEK_SET);
  return count;
}

void ReadMatrix(const char *filename, gsl_matrix ** m) {
  int rows, cols;
  FILE *f;
  f = fopen(filename, "r");
  rows = NumberOfRows(f);
  cols = NumberOfColumns(f);
  *m = gsl_matrix_alloc(rows, cols);
  gsl_matrix_fscanf(f, *m);
  fclose(f);
}

void PrintSubMatrix(gsl_matrix * m,
    size_t k1,
    size_t k2,
    size_t n1,
    size_t n2) {
  printf("\nsubmatrix: %zu, %zu, %zu, %zu\n", k1, k2, n1, n2);
  gsl_matrix_view sub = gsl_matrix_submatrix(m, k1, k2, n1, n2);
  for (size_t i = 0; i < n1; i++) {
    for (size_t j = 0; j < n2; j++) {
      printf("%.2f ", gsl_matrix_get(&sub.matrix, i, j));
    }
    printf("\n");
  }
}

void PrintMatrix(gsl_matrix * m) {
  for (size_t i = 0; i < m->size1; ++i) {
    for (size_t  j = 0; j < m->size2; ++j) {
      printf("%.2f ", gsl_matrix_get(m, i, j));
    }
    printf("\n");
  }
}

void PrintVector(gsl_vector * v) {
  for (size_t j = 0; j < v->size; j++) {
    printf("%.2f ", gsl_vector_get(v, j));
  }
  printf("\n");
}

void SphereMatrix(gsl_matrix * m) {
  gsl_vector *v;
  v = gsl_vector_alloc(m->size2);
  for (size_t i = 0; i < m->size1; ++i) {
    float sum = 0.0;
    for (size_t j = 0; j < m->size2; ++j) {
      sum += gsl_matrix_get(m, i, j);
    }
    gsl_matrix_get_row(v, m, i);
    double mean = gsl_stats_mean(v->data, 1, m->size2);
    double stdev = gsl_stats_sd(v->data, 1, m->size2);
    gsl_vector_add_constant(v, -mean);
    gsl_vector_scale(v, 1.0/stdev);
    gsl_matrix_set_row(m, i, v);
  }
  gsl_vector_free(v);
}

void CrossValidation(gsl_matrix * m, size_t splits) {
  printf("initial matrix:\n");
  PrintMatrix(m);

  vector<gsl_matrix> matrices;
  size_t rows, cols;
  gsl_matrix *m_temp;
  // We create a vector of matrices to dump our rows into
  for (size_t i = 0; i < splits; ++i) {
    if (i < m->size1 % splits)
      // Some matrices have an extra row, e.g. If our original matrix
      // has 10 rows but we want 3 splits, 10/3 = 3 remainder 1
      // So two matrices will have three rows, but one matrix will
      // have an extra row
      rows = m->size1/splits + 1;
    else
      rows = m->size1/splits;
    cols = m->size2;
    m_temp = gsl_matrix_alloc(rows, cols);  // allocate
    gsl_matrix_set_zero(m_temp);            // zero it out
    matrices.push_back(*m_temp);            // add it to the vector
  }

  // Now let's loop through our matrices and populate them
  // ((j*splits) < m->size1) is to handle the case of extra rows
  for (size_t j = 0; (j <= rows) && ((j*splits) < m->size1); ++j) {
    vector<gsl_vector> row_vector;
    // Scoop up some rows from the original matrix
    // ((i + j*splits) < m->size1) is to handle the case of extra rows
    for (size_t i = 0; (i < splits) && ((i + j*splits) < m->size1); ++i) {
      gsl_vector *a_single_row = gsl_vector_alloc(cols);
      gsl_matrix_get_row(a_single_row, m, i + j*splits);
      row_vector.push_back(*a_single_row);
    }
    // Shuffle up those rows that we just scooped up
    random_shuffle(row_vector.begin(), row_vector.end());
    // Now spread those rows evenly across our new matrices
    // ((i + j*splits) < m->size1) is to handle the case of extra rows
    for (size_t i = 0; (i < splits) && ((i + j*splits) < m->size1); ++i) {
      gsl_matrix_set_row(&matrices[i], j, &row_vector[i]);
    }
  }

  // Now let's print our new matrices out
  vector<gsl_matrix>::iterator it = matrices.begin();
  FOREACH(it, matrices) {
    printf("\nMatrix\n");
    PrintMatrix(&(*it));
  }
}

gsl_matrix * DiagAlloc(gsl_vector * v) {
    gsl_matrix *mat = gsl_matrix_alloc(v->size, v->size);
    gsl_vector_view diag = gsl_matrix_diagonal(mat);
    gsl_matrix_set_all(mat, 0.0);
    gsl_vector_memcpy(&diag.vector, v);
    return mat;
}

gsl_matrix * RepMatHorizAlloc(gsl_vector * v, size_t k) {
    gsl_matrix *mat = gsl_matrix_alloc(k, v->size);
    for (size_t i = 0; i < k; ++i) {
      gsl_matrix_set_row(mat, i, v);
    }
    return mat;
}

gsl_matrix * RepMatVertAlloc(gsl_vector * v, size_t k) {
    gsl_matrix *mat = gsl_matrix_alloc(v->size, k);
    for (size_t i = 0; i < k; ++i) {
      gsl_matrix_set_col(mat, i, v);
    }
    return mat;
}

gsl_matrix * CreateMatrix(double * arr, size_t size1, size_t size2) {
  gsl_matrix *m = gsl_matrix_alloc(size1, size2);
  for (size_t i = 0; i < size1; ++i) {
    for (size_t j = 0; j < size2; ++j) {
      gsl_matrix_set(m, i, j, *(arr + i*size2 + j));
    }
  }
  return m;
}

gsl_vector * CreateVector(double * arr, size_t size) {
  gsl_vector *v = gsl_vector_alloc(size);
  for (size_t i = 0; i < size; ++i) {
    gsl_vector_set(v, i, arr[i]);
  }
  return v;
}

gsl_vector * MultiplyVecMat(gsl_vector *v, gsl_matrix *m) {
  gsl_vector *result = gsl_vector_calloc(v->size);
  cblas_dgemm(CblasRowMajor,   // const enum CBLAS_ORDER Order
              CblasNoTrans ,   // const enum CBLAS_TRANSPOSE TransA
              CblasNoTrans ,   // const enum CBLAS_TRANSPOSE TransB
              1,               // const int M
              m->size2,        // const int N
              v->size,         // const int K
              1.0f,            // const double alpha
              v->data,         // const double * A
              v->size,         // const int lda
              m->data,         // const double * B
              m->size2,        // const int ldb
              0.0f,            // const double beta
              result->data,    // double * C
              result->size);   // const int ldc
  return result;
}

gsl_vector * MultiplyMatVec(gsl_matrix *m, gsl_vector *v) {
  gsl_vector *result = gsl_vector_calloc(v->size);
  cblas_dgemm(CblasRowMajor,   // const enum CBLAS_ORDER Order
              CblasNoTrans ,   // const enum CBLAS_TRANSPOSE TransA
              CblasTrans ,     // const enum CBLAS_TRANSPOSE TransB
              m->size1,        // const int M (height of A)
              v->size,         // const int N (width of B)
              m->size2,        // const int K (width of A)
              1.0f,            // const double alpha
              m->data,         // const double * A
              m->size2,        // const int lda
              v->data,         // const double * B
              v->size,         // const int ldb
              0.0f,            // const double beta
              result->data,    // double * C
              1);              // const int ldc
  return result;
}

double MultiplyVecVec(gsl_vector *v1, gsl_vector *v2) {
  double result;
  gsl_blas_ddot(v1, v2, &result);
  return result;
}

gsl_matrix * MultiplyMatMat(gsl_matrix *m1, gsl_matrix *m2) {
  gsl_matrix *result = gsl_matrix_alloc(m1->size1, m2->size1);
  cblas_dgemm(CblasRowMajor,   // const enum CBLAS_ORDER Order
              CblasNoTrans ,   // const enum CBLAS_TRANSPOSE TransA
              CblasTrans ,     // const enum CBLAS_TRANSPOSE TransB
              m1->size1,       // const int M
              m2->size1,       // const int N
              m1->size2,       // const int K
              1.0f,            // const double alpha
              m1->data,        // const double * A
              m1->size2,       // const int lda
              m2->data,        // const double * B
              m2->size2,       // const int ldb
              0.0f,            // const double beta
              result->data,    // double * C
              result->size2);  // const int ldc
  return result;
}

double LinearKernel(gsl_vector *v1, gsl_vector *v2) {
  return MultiplyVecVec(v1, v2);
}

double StandardizedLinearKernel(gsl_vector *v1, gsl_vector *v2) {
  float kernel = LinearKernel(v1, v2);
  float root1 = sqrt(LinearKernel(v1, v1));
  float root2 = sqrt(LinearKernel(v2, v2));
  return kernel/(root1*root2);
}

double PolynomialKernel(gsl_vector *v1, gsl_vector *v2, double d) {
  return pow((MultiplyVecVec(v1, v2) + 1), d);
}

double StandardizedPolynomialKernel(gsl_vector *v1, gsl_vector *v2, double d) {
  float kernel = PolynomialKernel(v1, v2, d);
  float root1 = sqrt(PolynomialKernel(v1, v1, d));
  float root2 = sqrt(PolynomialKernel(v2, v2, d));
  return kernel/(root1*root2);
}

double GaussianKernel(gsl_vector *v1, gsl_vector *v2, gsl_vector *t) {
  double result;
  gsl_matrix *m = DiagAlloc(t);
  gsl_vector_sub(v1, v2);
  gsl_vector *v1m = MultiplyVecMat(v1, m);
  result = MultiplyVecVec(v1m, v1);
  result = gsl_sf_exp(result);
  return result;
}

TEST(MRVMTest, submatrix) {
  gsl_matrix *mm;

  ReadMatrix("test.dat", &mm);

  // y, x, height, width
  gsl_matrix_view sub = gsl_matrix_submatrix(mm, 8, 1, 2, 2);
  printf("matrix\n");
  PrintMatrix(mm);
  PrintSubMatrix(mm, 8, 1, 2, 2);

  gsl_matrix_free(mm);
}

TEST(MRVMTest, sphering) {
  gsl_matrix *mm;

  ReadMatrix("test.dat", &mm);

  // print
  printf("rows: %zu cols: %zu\n", mm->size1, mm->size2);
  PrintMatrix(mm);

  // sphere and print
  SphereMatrix(mm);
  printf("\nSphered:\n");
  PrintMatrix(mm);

  gsl_matrix_free(mm);
}

TEST(MRVMTest, CrossValidation) {
  gsl_matrix *mm;
  ReadMatrix("test2.dat", &mm);
  CrossValidation(mm, 3);
  gsl_matrix_free(mm);
}

TEST(MRVMTest, shuffle) {
  vector<int> vec;
  int values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  for (size_t i = 0; i < 10; ++i) {
    vec.push_back(values[i]);
  }
  for (vector<int>::iterator iter = vec.begin();
      iter != vec.end();
      ++iter) {
    printf("%d\n", *iter);
  }
  printf("\nNow Random:\n");
  random_shuffle(vec.begin(), vec.end());
  for (vector<int>::iterator iter = vec.begin();
      iter != vec.end();
      ++iter) {
    printf("%d\n", *iter);
  }
}

TEST(MRVMTest, inner_product) {
  gsl_vector *v1 = gsl_vector_alloc(5);
  gsl_vector *v2 = gsl_vector_alloc(5);
  for (size_t i = 0; i < 5; ++i) {
    gsl_vector_set(v1, i, 1.23 + i);
    gsl_vector_set(v2, i, 2.23 + i*0.3);
  }
  printf("vector1:\n");
  gsl_vector_fprintf(stdout, v1, "%.5g");
  printf("vector2:\n");
  gsl_vector_fprintf(stdout, v2, "%.5g");
  printf("vector product\n");
  double result;
  gsl_blas_ddot(v1, v2, &result);
  printf("%f\n", result);
  gsl_vector_free(v1);
  gsl_vector_free(v2);
}

TEST(MRVMTest, linear_kernel) {
  double v1_arr[] = { 1, 2 };
  double v2_arr[] = { 2, 3 };
  gsl_vector *v1 = CreateVector(v1_arr, 2);
  gsl_vector *v2 = CreateVector(v2_arr, 2);
  double answer = LinearKernel(v1, v2);
  ASSERT_EQ(8, answer);
}

TEST(MRVMTest, linear_kernel_from_file) {
  gsl_matrix *m;
  ReadMatrix("test2.dat", &m);
  gsl_vector *v1 = gsl_vector_alloc(m->size2);
  gsl_vector *v2 = gsl_vector_alloc(m->size2);
  gsl_matrix_get_row(v1, m, 0);
  gsl_matrix_get_row(v2, m, 1);
  PrintVector(v1);
  PrintVector(v2);
  double answer = LinearKernel(v1, v2);
  ASSERT_EQ(6, answer);
}

TEST(MRVMTest, standardized_linear_kernel) {
  double v1_arr[] = { 1, 2 };
  double v2_arr[] = { 2, 3 };
  gsl_vector *v1 = CreateVector(v1_arr, 2);
  gsl_vector *v2 = CreateVector(v2_arr, 2);
  double answer = StandardizedLinearKernel(v1, v2);
  ASSERT_EQ(0.99227787671366774, answer);
}

TEST(MRVMTest, polynomial_kernel) {
  double v1_arr[] = { 1, 2 };
  double v2_arr[] = { 4, 3 };
  gsl_vector *v1 = CreateVector(v1_arr, 2);
  gsl_vector *v2 = CreateVector(v2_arr, 2);
  double answer = PolynomialKernel(v1, v2, 2);
  ASSERT_EQ(121, answer);
}

TEST(MRVMTest, standardized_polynomial_kernel) {
  double v1_arr[] = { 1, 2 };
  double v2_arr[] = { 2, 3 };
  gsl_vector *v1 = CreateVector(v1_arr, 2);
  gsl_vector *v2 = CreateVector(v2_arr, 2);
  double answer = StandardizedPolynomialKernel(v1, v2, 2);
  ASSERT_EQ(0.9642857142857143, answer);
}

TEST(MRVMTest, gaussian_kernel) {
  double v1_arr[] = { 1, 2 };
  double v2_arr[] = { 4, 3 };
  double v3_arr[] = { 1, 2 };
  gsl_vector *v1 = CreateVector(v1_arr, 2);
  gsl_vector *v2 = CreateVector(v2_arr, 2);
  gsl_vector *v3 = CreateVector(v3_arr, 2);
  double answer = GaussianKernel(v1, v2, v3);
  ASSERT_EQ(59874.141715197817, answer);
}

TEST(MRVMTest, diagonal) {
  double v[] = {1, 2, 3};
  gsl_vector *v1 = CreateVector(v, 3);
  printf("Original Vector\n");
  PrintVector(v1);
  printf("Diagonal Matrix\n");
  PrintMatrix(DiagAlloc(v1));
}

TEST(MRVMTest, repmat) {
  double v[] = {1, 2, 3};
  gsl_vector *v1 = CreateVector(v, 3);
  printf("Original Vector\n");
  PrintVector(v1);
  printf("Horiz Repmat Matrix\n");
  PrintMatrix(RepMatHorizAlloc(v1, 7));
  printf("Vert Repmat Matrix\n");
  PrintMatrix(RepMatVertAlloc(v1, 7));
}

TEST(MRMVTest, mul_vec_mat) {
  double v_arr[] = { 1, 2 };
  double m_arr[] = { 1, 2,
                     3, 4 };
  gsl_vector *v = CreateVector(v_arr, 2);
  gsl_matrix *m = CreateMatrix(m_arr, 2, 2);
  printf("vector:\n");
  PrintVector(v);
  printf("matrix:\n");
  PrintMatrix(m);
  gsl_vector *result = MultiplyVecMat(v, m);
  printf("result1:\n");
  PrintVector(result);
  result = MultiplyMatVec(m, v);
  printf("result2:\n");
  PrintVector(result);
}

TEST(MRMVTest, full_linear) {
  gsl_matrix *m1, *m2;
  ReadMatrix("test3.dat", &m1);
  m2 = CreateMatrix(m1->data, m1->size1, m1->size2);
  printf("Matrix1:\n");
  PrintMatrix(m1);
  printf("Matrix2:\n");
  PrintMatrix(m2);
  printf("Product:\n");
  PrintMatrix(MultiplyMatMat(m1, m2));
}

TEST(MRMVTest, pairwise_division) {
  double m1_arr[] = { 1, 1,
                      1, 1 };
  gsl_matrix *m1 = CreateMatrix(m1_arr, 2, 2);
  double m2_arr[] = { 2, 2,
                      2, 2 };
  gsl_matrix *m2 = CreateMatrix(m2_arr, 2, 2);
  printf("Matrix1:\n");
  PrintMatrix(m1);
  printf("Matrix2:\n");
  PrintMatrix(m2);
  gsl_matrix_div_elements(m1, m2);
  printf("Quotient:\n");
  PrintMatrix(m1);
}
