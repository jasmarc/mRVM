#!/bin/bash
VERBOSITY=${1-1}

valgrind \
  --track-origins=yes \
  --leak-check=full \
  --show-reachable=yes \
./bin/mRVM \
  -k GAUSSIAN \
  -v $VERBOSITY \
  --train   ./data/dev.train.dat \
  --labels  ./data/dev.labels.dat \
  --test    ./data/dev.test.dat \
  --answers ./data/dev.answers.dat \
  --out     ./data/dev.out.dat \
  --param   2 \
  --tau     0.000001 \
  --upsilon 0.000001

# 1st Probit Kernel Regression, Kernel-based with multi-nomial probit likilhood
# 2nd Relevence Vector Machine (mRVM2), multi-class, top-down approach
# 3rd Relevence Vector Machine (mRVM1), multi-class, bottom-up approach