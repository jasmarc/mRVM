// Copyright 2011 Jason Marcell

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>

#include "lib/Matrix.h"
#include "lib/Kernel.h"
#include "lib/Trainer.h"
#include "lib/Predictor.h"
#include "lib/GaussHermiteQuadrature.h"
#include "lib/Log.h"

#define PACKAGE "mRVM"
#define VERSION "0.0.1"

#define PACKAGE    "mRVM"
#define VERSION    "0.0.1"

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

void print_help(int exval);
void run();
void handleKernelOption(KernelType & kernel);

int main(int argc, char **argv) {
  int opt;
  int long_opt_index = 0;
  KernelType kernel = LINEAR;
  int longval;

  // no arguments given
  if(argc == 1) {
    print_help(1);
  }

  struct option long_options[] = {
      { "help",     0, NULL,      'h' },
      { "version",  0, NULL,      'V' },
      { "verbose",  0, NULL,      'v' },
      { "kernel",   1, &longval,  'k' },
      { 0,          0, 0,         0  }
  };

  while((opt = getopt_long(argc, argv, "hVvf:k:", long_options, &long_opt_index)) != -1) {
    switch(opt) {
    case 'h':
      print_help(0);
      break;
    case 'V':
      printf("%s %s\n\n", PACKAGE, VERSION);
      exit(0);
      break;
    case 'v':
      printf("%s: Verbose option is set `%c'\n", PACKAGE, optopt);
      break;
    case 'k':
      printf("%s: Kernel %s\n", PACKAGE, optarg);
      handleKernelOption(kernel);
    break;
    case 'f':
      printf("%s: Filename %s\n", PACKAGE, optarg);
      break;
    case ':':
      fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE, optopt);
      print_help(1);
      break;
    case '?':
      fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
      print_help(1);
    case 0:
      switch (longval) {
      case 'k':
        printf("%s: Kernel %s\n", PACKAGE, optarg);
        handleKernelOption(kernel);
        break;
      }
    }
  }

  // print all remaining options
  for(; optind < argc; optind++)
    printf("argument: %s\n", argv[optind]);

  run();

  return 0;
}

void handleKernelOption(KernelType &kernel) {
  if (strcmp(optarg, "LINEAR") == 0) {
    kernel = LINEAR;
  } else if (strcmp(optarg, "POLYNOMIAL") == 0) {
    kernel = POLYNOMIAL;
  } else if (strcmp(optarg, "GAUSSIAN") == 0) {
    kernel = GAUSSIAN;
  }
}

void print_help(int exval) {
  printf("%s, %s multi-class multi-kernel Relevance Vector Machines (mRVM)\n", PACKAGE, VERSION);
  printf("%s [-h] [-V] [-f FILE] [-o FILE]\n\n", PACKAGE);

  printf("  -h, --help      print this help and exit\n");
  printf("  -V, --version   print version and exit\n\n");

  printf("  -v, --verbose   set verbose flag\n");
  printf("  -k, --kernel    specify the kernel: LINEAR, POLYNOMIAL, GAUSSIAN\n");
  printf("  -f FILE         set input file\n");

  printf("Based upon work by Psorakis, Damoulas, Girolami.\n");
  printf("Implementation by Marcell, jasonmarcell@gmail.com\n\n");

  exit(exval);
}

void run()
{
  // The number of classes
  const size_t CLASSES = 2;

  // Training points
  double train_arr[] =  \
      {1, 0.9,
    0.8, 0.9,
    0.5, 0.6,
    0.1, 0.2,
    0.2, 0.3};
  Matrix *train = new Matrix(train_arr, 5, 2);

  // Labels
  double labels_arr[] = {0, 0, 0, 1, 1};
  Vector *labels = new Vector(labels_arr, 5);

  // Testing Points
  double test_arr[] =
      {1, 0.9,
    0.8, 0.9,
    0.1, 0.2,
    0.2, 0.3,
    0.7, 0.7};
  Matrix *test = new Matrix(test_arr, 5, 2);

  LOG(VERBOSE, "=== Starting... ===\n");

  // Pass in training points, labels, and number of classes
  Trainer *trainer = new Trainer(train, labels, CLASSES);
  trainer->Process();

  // Pass in the w matrix, the training points, and the testing points
  Predictor *predictor = new Predictor(trainer->GetW(), train, test);
  predictor->Predict();

  LOG(VERBOSE, "=== End. ===\n");
}
