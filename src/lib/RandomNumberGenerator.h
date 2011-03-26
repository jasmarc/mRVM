// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_RANDOMNUMBERGENERATOR_H_
#define SRC_LIB_RANDOMNUMBERGENERATOR_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace jason {

class RandomNumberGenerator {
  public:
    RandomNumberGenerator();
    virtual ~RandomNumberGenerator();
    double SampleGaussian(double sigma);
    double GaussianPDF(double val);
    double GaussianCDF(double val);
    double SampleUniform(double lower, double upper);
  private:
    gsl_rng * r;
};
}

#endif  // SRC_LIB_RANDOMNUMBERGENERATOR_H_
