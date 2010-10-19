// Copyright 2010 Jason Marcell

#include <gtest/gtest.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
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

void print_sub_matrix(gsl_matrix * m,
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

void print_matrix(gsl_matrix * m) {
  for (size_t  i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      printf("%.2f ", gsl_matrix_get(m, i, j));
    }
    printf("\n");
  }
}

void sphere_matrix(gsl_matrix * m) {
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

void cross_validation(gsl_matrix * m, size_t splits) {
  printf("initial matrix:\n");
  print_matrix(m);

  vector<gsl_matrix> matrices;
  int rows, cols;
  gsl_matrix *m_temp;
  for (size_t i = 0; i < splits; ++i) {
    if (i < m->size1 % splits)
      rows = m->size1/splits + 1;
    else
      rows = m->size1/splits;
    cols = m->size2;
    m_temp = gsl_matrix_alloc(rows, cols);
    gsl_matrix_set_zero(m_temp);
    matrices.push_back(*m_temp);
  }

  for (size_t j = 0; (j <= rows) && ((j*splits) < m->size1); ++j) {
    vector<gsl_vector> row_vector;
    for (size_t i = 0; (i < splits) && ((i + j*splits) < m->size1); ++i) {
      gsl_vector *a_single_row = gsl_vector_alloc(cols);
      gsl_matrix_get_row(a_single_row, m, i + j*splits);
      row_vector.push_back(*a_single_row);
    }
    random_shuffle(row_vector.begin(), row_vector.end());
    for (size_t i = 0; (i < splits) && ((i + j*splits) < m->size1); ++i) {
      gsl_matrix_set_row(&matrices[i], j, &row_vector[i]);
    }
  }

  vector<gsl_matrix>::iterator it = matrices.begin();
  FOREACH(it, matrices) {
    printf("\nMatrix\n");
    print_matrix(&(*it));
  }
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
  print_matrix(mm);
  print_sub_matrix(mm, 8, 1, 2, 2);

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
  print_matrix(mm);

  // sphere and print
  sphere_matrix(mm);
  printf("\nSphered:\n");
  print_matrix(mm);

  gsl_matrix_free(mm);
}

TEST(MRVMTest, cross_validation) {
  gsl_matrix *mm;

  // read matrix from file
  int rows, cols;
  FILE *f;
  f = fopen("test2.dat", "r");
  fscanf(f, "%d %d", &rows, &cols);
  mm = gsl_matrix_alloc(rows, cols);
  gsl_matrix_fscanf(f, mm);
  fclose(f);

  cross_validation(mm, 6);

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
