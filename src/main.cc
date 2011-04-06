// Copyright 2011 Jason Marcell

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>

#include "lib/Matrix.h"
#include "lib/Vector.h"
#include "lib/Kernel.h"
#include "lib/Trainer.h"
#include "lib/Predictor.h"
#include "lib/GaussHermiteQuadrature.h"
#include "lib/Log.h"

#define PACKAGE "mRVM"
#define VERSION "0.0.1"

using jason::Matrix;
using jason::Vector;
using jason::Trainer;
using jason::Kernel;
using jason::Predictor;
using jason::GaussHermiteQuadrature;
using jason::KernelType;
using jason::LINEAR;
using jason::POLYNOMIAL;
using jason::GAUSSIAN;

extern int verbosity;

void print_help(int exval);
void run(char *train_filename, char *labels_filename, char *test_filename);
void handleFile(char **filename);
void handleVerbosity();
void handleKernelOption(KernelType *kernel);

int main(int argc, char **argv) {
  int opt;
  int long_opt_index = 0;
  int longval;
  
  char *train_filename;
  char *labels_filename;
  char *test_filename;
  KernelType kernel = LINEAR;

  // no arguments given
  if (argc == 1) {
    print_help(1);
  }

  struct option long_options[] = {
      { "help",     0, NULL,      'h' },
      { "version",  0, NULL,      'V' },
      { "verbose",  1, &longval,  'v' },
      { "train",    1, &longval,  'r' },
      { "labels",   1, &longval,  'l' },
      { "test",     1, &longval,  't' },
      { "kernel",   1, &longval,  'k' },
      { 0,          0, 0,         0  }
  };

  while ((opt = getopt_long(argc, argv, "hVv:r:l:t:k:", long_options,
    &long_opt_index)) != -1) {
    switch (opt) {
    case 'h':
      print_help(0);
      break;
    case 'V':
      printf("%s %s\n\n", PACKAGE, VERSION);
      exit(0);
      break;
    case 'v':
      handleVerbosity();
      break;
    case 'r':
      handleFile(&train_filename);
      break;
    case 'l':
      handleFile(&labels_filename);
      break;
    case 't':
      handleFile(&test_filename);
      break;
    case 'k':
      handleKernelOption(&kernel);
      break;
    case ':':
      fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE,
        optopt);
      print_help(1);
      break;
    case '?':
      fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
      print_help(1);
    case 0:
      switch (longval) {
      case 'v':
        handleVerbosity();
        break;
      case 'r':
        handleFile(&train_filename);
        break;
      case 'l':
        handleFile(&labels_filename);
        break;
      case 't':
        handleFile(&test_filename);
        break;
      case 'k':
        handleKernelOption(&kernel);
        break;
      }
    }
  }

  // print all remaining options
  for (; optind < argc; optind++)
    printf("argument: %s\n", argv[optind]);

  const char *str_kernel = NULL;
  switch (kernel) {
  case LINEAR:
    str_kernel = "LINEAR";
    break;
  case POLYNOMIAL:
    str_kernel = "POLYNOMIAL";
    break;
  case GAUSSIAN:
    str_kernel = "GAUSSIAN";
    break;
  default:
    str_kernel = "UNKNOWN";
  }
  LOG(VERBOSE, "Verbosity level = %d\n", verbosity)
  LOG(VERBOSE, "Kernel          = %s\n", str_kernel);
  LOG(VERBOSE, "Training file   = %s\n", train_filename);
  LOG(VERBOSE, "Labels file     = %s\n", labels_filename);
  LOG(VERBOSE, "Test file       = %s\n", test_filename);

  run(train_filename, labels_filename, test_filename);

  return 0;
}

void handleFile(char **filename) {
  *filename = optarg;
}

void handleVerbosity() {
  verbosity = atoi(optarg);
}

void handleKernelOption(KernelType *kernel) {
  if (strcmp(optarg, "LINEAR") == 0) {
    *kernel = LINEAR;
  } else if (strcmp(optarg, "POLYNOMIAL") == 0) {
    *kernel = POLYNOMIAL;
  } else if (strcmp(optarg, "GAUSSIAN") == 0) {
    *kernel = GAUSSIAN;
  }
}

void print_help(int exval) {
  printf("%s, %s multi-class multi-kernel Relevance Vector Machines (mRVM)\n",
    PACKAGE, VERSION);
  printf("%s [-h] [-V] [-f FILE] [-o FILE]\n\n", PACKAGE);

  printf("  -h, --help         print this help and exit\n");
  printf("  -V, --version      print version and exit\n\n");

  printf("  -r, --train  FILE  set input file\n");
  printf("  -l, --labels FILE  set input file\n");
  printf("  -t, --test   FILE  set input file\n");
  printf("  -k, --kernel       specify the kernel:\n");
  printf("                       LINEAR\n");
  printf("                       POLYNOMIAL\n");
  printf("                       GAUSSIAN\n");
  printf("  -v, --verbose      set verbosity level:\n");
  printf("                       0 = No output\n");
  printf("                       1 = Normal (default)\n");
  printf("                       2 = Verbose\n");
  printf("                       3 = Debug\n");

  printf("Based upon work by Psorakis, Damoulas, Girolami.\n");
  printf("Implementation by Marcell, jasonmarcell@gmail.com\n\n");

  exit(exval);
}

void run(char *train_filename, char *labels_filename, char *test_filename) {
  Matrix *train  = new Matrix(train_filename);
  Vector *labels = new Vector(labels_filename);
  Matrix *test   = new Matrix(test_filename);
  size_t classes = labels->GetNumberOfClasses();

  LOG(VERBOSE, "=== Starting... ===\n");

  // Pass in training points, labels, and number of classes
  Trainer *trainer = new Trainer(train, labels, classes);
  trainer->Process();

  // Pass in the w matrix, the training points, and the testing points
  Predictor *predictor = new Predictor(trainer->GetW(), train, test);
  predictor->Predict();

  LOG(VERBOSE, "=== End. ===\n");
}
