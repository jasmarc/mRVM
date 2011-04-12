// Copyright 2011 Jason Marcell

#ifndef SRC_MAIN_H_
#define SRC_MAIN_H_

#define PACKAGE "mRVM"
#define VERSION "0.0.2"

extern int verbosity;

namespace jason {

int main(int argc, char **argv);
void print_help(int exval);
void run(char *train_filename, char *labels_filename, char *test_filename,
    char *answers_filename, char *out_filename, KernelType kernel_type,
    int kernel_param, double tau, double upsilon);
void handleKernelOption(KernelType *kernel, char **kernel_str);
void PerformEvaluation(Matrix *predictions, Vector *answers);
}

#endif  // SRC_MAIN_H_
