// Copyright 2011 Jason Marcell

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>

#include "lib/Matrix.h"
#include "lib/Kernel.h"
#include "lib/Trainer.h"
#include "lib/Predictor.h"
#include "lib/GaussHermiteQuadrature.h"

#define PACKAGE    "mRVM"
#define VERSION    "0.0.1"

using jason::Matrix;
using jason::Vector;
using jason::Trainer;
using jason::Kernel;
using jason::Predictor;
using jason::GaussHermiteQuadrature;

void print_help(int exval);
void run();

int main(int argc, char **argv) {


  int opt;

  // no arguments given
  if(argc == 1) {
    print_help(1);
  }

  while((opt = getopt(argc, argv, "hVvf:o:")) != -1) {
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
    case 'f':
      printf("%s: Filename %s\n", PACKAGE, optarg);
      break;
    case 'o':
      printf("Output: %s\n", optarg);
      break;
    case ':':
      fprintf(stderr, "%s: Error - Option `%c' needs a value\n\n", PACKAGE, optopt);
      print_help(1);
      break;
    case '?':
      fprintf(stderr, "%s: Error - No such option: `%c'\n\n", PACKAGE, optopt);
      print_help(1);
    }
  }

  // print all remaining options
  for(; optind < argc; optind++)
    printf("argument: %s\n", argv[optind]);

  run();

  return 0;
}

void print_help(int exval) {
  printf("%s, %s multi-class multi-kernel Relevance Vector Machines (mRVM)\n", PACKAGE, VERSION);
  printf("%s [-h] [-V] [-f FILE] [-o FILE]\n\n", PACKAGE);

  printf("  -h              print this help and exit\n");
  printf("  -V              print version and exit\n\n");

  printf("  -v              set verbose flag\n");
  printf("  -f FILE         set input file\n");
  printf("  -o FILE         set output file\n\n");

  printf("Based upon work by Psorakis, Damoulas, Girolami.\n");
  printf("Implementation by Damoulas and Marcell, {damoulas, jrm425}@cs.cornell.edu\n\n");

  exit(exval);
}

void run()
{
    double m1_arr[] = {1, 2, 9, 4, 6, 5, 3, 7, 4, 6, 2, 0};
    double m2_arr[] = {1, 2, 9, 4, 6, 5, 6, 2, 0};
    double v_arr[] = {0, 1, 0, 3};
    printf("Starting...\n");
    Matrix *m1 = new Matrix(m1_arr, 4, 3);
    Vector *v = new Vector(v_arr, 4);
    Trainer *t = new Trainer(m1, v, 4);
    t->Process();
    Matrix *m2 = new Matrix(m2_arr, 3, 3);
    Predictor *p = new Predictor(t->GetW(), v, m1, m2);
    p->Predict();
    printf("End.\n");
}
