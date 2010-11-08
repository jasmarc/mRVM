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
  for (size_t j = 0; j < m->size2; ++j) {
    for (size_t  i = 0; i < m->size1; ++i) {
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

gsl_matrix * CreateMatrix(double * arr, size_t size1, size_t size2) {
  gsl_matrix *m = gsl_matrix_alloc(size1, size2);
  for (size_t j = 0; j < size2; ++j) {
    for (size_t i = 0; i < size1; ++i) {
      gsl_matrix_set(m, i, j, *(arr + i + j*size2));
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

void MultiplyMatrix(double *a, unsigned int aheight, unsigned int awidth,
                    double *b, unsigned int bwidth, double *c) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, aheight, bwidth,
                awidth, 1.0f, a, aheight, b, awidth, 0.0f, c, aheight);
}

gsl_vector * MultiplyVecMat(gsl_vector *v, gsl_matrix *m) {
  gsl_vector *result = gsl_vector_calloc(v->size);
  MultiplyMatrix(v->data, 1, v->size, m->data, m->size2, result->data);
  return result;
}

gsl_vector * MultiplyMatVec(gsl_matrix *m, gsl_vector *v) {
  gsl_vector *result = gsl_vector_calloc(v->size);
  MultiplyMatrix(m->data, m->size1, m->size2, v->data, 1, result->data);
  return result;
}

double MultiplyVecVec(gsl_vector *v1, gsl_vector *v2) {
  double result;
  gsl_blas_ddot(v1, v2, &result);
  return result;
}

double LinearKernel(gsl_vector *v1, gsl_vector *v2) {
  return MultiplyVecVec(v1, v2);
}

double PolynomialKernel(gsl_vector *v1, gsl_vector *v2, double d) {
  return pow((MultiplyVecVec(v1, v2) + 1), d);
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

  // read matrix from file
  int rows, cols;
  FILE *f;
  f = fopen("test.dat", "r");
  fscanf(f, "%d %d", &rows, &cols);
  mm = gsl_matrix_alloc(rows, cols);
  gsl_matrix_fscanf(f, mm);
  fclose(f);

  // y, x, height, width
  gsl_matrix_view sub = gsl_matrix_submatrix(mm, 8, 1, 2, 2);
  printf("matrix\n");
  PrintMatrix(mm);
  PrintSubMatrix(mm, 8, 1, 2, 2);

  gsl_matrix_free(mm);
}

TEST(MRVMTest, sphering) {
  gsl_matrix *mm;

  // read matrix from file
  int rows, cols;
  FILE *f;
  f = fopen("test.dat", "r");
  fscanf(f, "%d %d", &rows, &cols);
  mm = gsl_matrix_alloc(rows, cols);
  gsl_matrix_fscanf(f, mm);
  fclose(f);

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

  // read matrix from file
  int rows, cols;
  FILE *f;
  f = fopen("test2.dat", "r");
  fscanf(f, "%d %d", &rows, &cols);
  mm = gsl_matrix_alloc(rows, cols);
  gsl_matrix_fscanf(f, mm);
  fclose(f);

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

TEST(MRVMTest, polynomial_kernel) {
  double v1_arr[] = { 1, 2 };
  double v2_arr[] = { 4, 3 };
  gsl_vector *v1 = CreateVector(v1_arr, 2);
  gsl_vector *v2 = CreateVector(v2_arr, 2);
  double answer = PolynomialKernel(v1, v2, 2);
  ASSERT_EQ(121, answer);
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
  printf("result2:\n");
  result = MultiplyMatVec(m, v);
  PrintVector(result);
}
