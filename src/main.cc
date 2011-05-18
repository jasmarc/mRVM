// Copyright 2011 Jason Marcell

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>

#include "lib/Matrix.h"
#include "lib/Vector.h"
#include "lib/Kernel.h"
#include "lib/LinearKernel.h"
#include "lib/PolynomialKernel.h"
#include "lib/GaussianKernel.h"
#include "lib/Trainer.h"
#include "lib/Predictor.h"
#include "lib/GaussHermiteQuadrature.h"
#include "lib/Log.h"
#include "./main.h"

int main(int argc, char **argv) {
  jason::main(argc, argv);
}

namespace jason {

int main(int argc, char **argv) {
  int opt = 0;
  int long_opt_index = 0;

  char *train_filename = NULL;
  char *labels_filename = NULL;
  char *test_filename = NULL;
  char *answers_filename = NULL;
  char *out_filename = NULL;
  KernelType kernel = LINEAR;
  char *str_kernel = NULL;
  int kernel_param = -1;
  double tau = 0;
  double upsilon = 0;

  // no arguments given
  if (argc == 1) {
    print_help(1);
  }

  struct option long_options[] = {
      { "help",     0, NULL,      'h' },
      { "version",  0, NULL,      'V' },
      { "verbose",  1, NULL,      'v' },
      { "train",    1, NULL,      'r' },
      { "labels",   1, NULL,      'l' },
      { "test",     1, NULL,      't' },
      { "answers",  1, NULL,      'a' },
      { "out",      1, NULL,      'o' },
      { "kernel",   1, NULL,      'k' },
      { "param",    1, NULL,      'p' },
      { "tau",      1, NULL,      'T' },
      { "upsilon",  1, NULL,      'u' },
      { 0,          0, 0,         0  }
  };

  while ((opt = getopt_long(argc, argv, "hVv:r:l:t:a:k:p:T:u:",
      long_options, &long_opt_index)) != -1) {
    switch (opt) {
    case 'h':
      print_help(0);
      break;
    case 'V':
      printf("%s %s\n\n", PACKAGE, VERSION);
      exit(0);
      break;
    case 'v':
      verbosity = atoi(optarg);
      break;
    case 'r':
      train_filename = optarg;
      break;
    case 'l':
      labels_filename = optarg;
      break;
    case 't':
      test_filename = optarg;
      break;
    case 'a':
      answers_filename = optarg;
      break;
    case 'o':
      out_filename = optarg;
      break;
    case 'k':
      handleKernelOption(&kernel, &str_kernel);
      break;
    case 'p':
      kernel_param = atoi(optarg);
      break;
    case 'T':
      tau = atof(optarg);
      break;
    case 'u':
      upsilon = atof(optarg);
      break;
    case ':':
      fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE,
        optopt);
      print_help(1);
      break;
    case '?':
      fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
      print_help(1);
      break;
    }
  }

  // print all remaining options
  for (; optind < argc; optind++)
    printf("argument: %s\n", argv[optind]);

  LOG(VERBOSE, "Verbosity level = %d\n", verbosity)
  LOG(VERBOSE, "Kernel          = %s\n", str_kernel);
  LOG(VERBOSE, "Training file   = %s\n", train_filename);
  LOG(VERBOSE, "Labels file     = %s\n", labels_filename);
  LOG(VERBOSE, "Test file       = %s\n", test_filename);
  LOG(VERBOSE, "Answers file    = %s\n", answers_filename);
  LOG(VERBOSE, "Out file        = %s\n", out_filename);
  LOG(VERBOSE, "Kernel param    = %d\n", kernel_param);
  LOG(VERBOSE, "Tau param       = %.3f\n", tau);
  LOG(VERBOSE, "Upsilon param   = %.3f\n", upsilon);

  if (train_filename == NULL) {
    fprintf(stderr, "%s: Error - Training file must be specified.\n\n",
        PACKAGE);
    print_help(1);
  } else if (labels_filename == NULL) {
    fprintf(stderr, "%s: Error - Labels file must be specified.\n\n", PACKAGE);
    print_help(1);
  } else if (test_filename == NULL) {
    fprintf(stderr, "%s: Error - Test file must be specified.\n\n", PACKAGE);
    print_help(1);
  } else if (kernel_param == -1 && kernel != LINEAR) {
    fprintf(stderr, "%s: Error - Must specify param for non-linear kernel.\n\n",
        PACKAGE);
    print_help(1);
  }

  run(train_filename, labels_filename, test_filename, answers_filename,
      out_filename, kernel, kernel_param, tau, upsilon);

  return 0;
}

void handleKernelOption(KernelType *kernel, char **kernel_str) {
  *kernel_str = optarg;
  if (strcmp(optarg, "LINEAR") == 0) {
    *kernel = LINEAR;
  } else if (strcmp(optarg, "POLYNOMIAL") == 0) {
    *kernel = POLYNOMIAL;
  } else if (strcmp(optarg, "GAUSSIAN") == 0) {
    *kernel = GAUSSIAN;
  } else {
    fprintf(stderr, "%s: Error - Unknown Kernel Specified.\n\n", PACKAGE);
    print_help(1);
  }
}

void print_help(int exval) {
  printf("%s, %s multi-class multi-kernel Relevance Vector Machines (mRVM)\n",
    PACKAGE, VERSION);
  printf("%s [options]\n\n", PACKAGE);

  printf("  -h, --help         print this help and exit\n");
  printf("  -V, --version      print version and exit\n\n");

  printf("  -r, --train   FILE set input file (required)\n");
  printf("  -l, --labels  FILE set input file (required)\n");
  printf("  -t, --test    FILE set input file (required)\n");
  printf("  -a, --answers FILE set input file (required)\n");
  printf("  -o, --out     FILE set output file\n");
  printf("  -k, --kernel       specify the kernel:\n");
  printf("                       LINEAR (default)\n");
  printf("                       POLYNOMIAL\n");
  printf("                       GAUSSIAN\n");
  printf("  -v, --verbose      set verbosity level:\n");
  printf("                       0 = No output\n");
  printf("                       1 = Normal (default)\n");
  printf("                       2 = Verbose\n");
  printf("                       3 = Debug\n");
  printf("  -p, --param n      set param for poly or gauss\n");
  printf("                     kernel to n.\n");
  printf("  -T, --tau n        set tau parameter\n");
  printf("  -u, --upsilon n    set upsilon parameter\n\n");

  printf("Based upon work by Psorakis, Damoulas, Girolami.\n");
  printf("Implementation by Marcell, jasonmarcell@gmail.com\n\n");

  exit(exval);
}

void run(char *train_filename, char *labels_filename, char *test_filename,
    char *answers_filename, char *out_filename, KernelType kernel_type,
    int kernel_param, double tau, double upsilon) {
  Matrix *train  = new Matrix(train_filename);
  Vector *labels = new Vector(labels_filename);
  Matrix *test   = new Matrix(test_filename);
  size_t classes = labels->GetNumberOfClasses();

  LOG(VERBOSE, "=== Starting... ===\n");

  train->CacheMeansAndStdevs();
  train->Sphere();
  test->Sphere(train);

  Kernel *train_kernel;
  Kernel *test_kernel;
  switch (kernel_type) {
  case LINEAR:
    LOG(DEBUG, "Creating Linear Kernels.\n");
    train_kernel = new LinearKernel(train, train);
    break;
  case POLYNOMIAL:
    LOG(DEBUG, "Creating Polynomial Kernels.\n");
    train_kernel = new PolynomialKernel(train, train, kernel_param);
    break;
  case GAUSSIAN:
    LOG(DEBUG, "Creating Gaussian Kernels.\n");
    train_kernel = new GaussianKernel(train, train, kernel_param);
    break;
  default:
    fprintf(stderr, "%s: Error - No such kernel.\n", PACKAGE);
    print_help(1);
  }

  // Pass in training points, labels, and number of classes
  Trainer *trainer = new Trainer(train, labels, classes, train_kernel);
  trainer->Process(tau, upsilon);

  switch (kernel_type) {
  case LINEAR:
    LOG(DEBUG, "Creating Linear Kernels.\n");
    test_kernel = new LinearKernel(train, test);
    break;
  case POLYNOMIAL:
    LOG(DEBUG, "Creating Polynomial Kernels.\n");
    test_kernel = new PolynomialKernel(train, test, kernel_param);
    break;
  case GAUSSIAN:
    LOG(DEBUG, "Creating Gaussian Kernels.\n");
    test_kernel = new GaussianKernel(train, test, kernel_param);
    break;
  default:
    fprintf(stderr, "%s: Error - No such kernel.\n", PACKAGE);
    print_help(1);
  }

  // Pass in the w matrix, the training points, and the testing points
  Predictor *predictor = new Predictor(trainer->GetW(), train, test,
      test_kernel);
  Matrix *predictions = predictor->Predict();

  LOG(VERBOSE, "= Predictions: =\n");
  LOG(VERBOSE, "%s\n", predictions->ToString());

  if (answers_filename) {
    Vector *answers = new Vector(answers_filename);
    PerformEvaluation(predictions, answers);
    delete answers;
  }

  if (out_filename) {
    LOG(VERBOSE, "Writing to file %s.\n", out_filename);
    predictions->Write(out_filename);
  }

  delete train;
  delete labels;
  delete test;
  delete predictions;
  delete predictor;
  delete trainer;
  delete train_kernel;
  delete test_kernel;

  LOG(VERBOSE, "=== End. ===\n");
}

void PerformEvaluation(Matrix *predictions, Vector *answers) {
  LOG(DEBUG, "= Evaluation =\n");
  size_t total_correct = 0;
  for (size_t row = 0; row < predictions->Height(); ++row) {
    double max = 0.0;
    size_t max_index = 0;
    for (size_t col = 0; col < predictions->Width(); ++col) {
      double current = predictions->Get(row, col);
      if (current > max) {
        max = current;
        max_index = col;
      }
    }
    size_t answer = answers->Get(row);
    LOG(VERBOSE, "Prediction: %zu, \tAnswer %zu\t", max_index, answer);
    if (max_index == answer) {
      LOG(VERBOSE, "\tCORRECT\n");
      total_correct += 1;
    } else {
      LOG(VERBOSE, "\tINCORRECT\n");
    }
  }
  LOG(NORMAL, "Percent correct: %.3f\n",
      static_cast<double>(total_correct) / predictions->Height());
}
}
