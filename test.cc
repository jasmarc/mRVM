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

using std::set;

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
  const gsl_rng_type * T;
  gsl_rng * r;
  set<size_t> myset;

  // Setup random number generator
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  while (myset.size() < splits) {
    float rand = gsl_ran_flat(r, 0 + 1, m->size1 - 1);
    int split = static_cast<int>(rand + 0.5);
    myset.insert(split);
  }

  // y, x, height, width
  int previous = 0;
  for (set<size_t>::const_iterator iter = myset.begin();
      iter != myset.end();
      ++iter) {
    print_sub_matrix(m, previous, 0, *iter - previous, m->size2 - 1);
    previous = *iter;
  }
  print_sub_matrix(m, previous, 0, m->size1 - previous, m->size2 - 1);

  gsl_rng_free(r);
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
  f = fopen("test.dat", "r");
  fscanf(f, "%d %d", &rows, &cols);
  mm = gsl_matrix_alloc(rows, cols);
  gsl_matrix_fscanf(f, mm);
  fclose(f);

  cross_validation(mm, 2);

  gsl_matrix_free(mm);
}
