#!/bin/bash
VERBOSITY=${1-1}

./bin/mRVM \
  -k LINEAR \
  -v $VERBOSITY \
  --train   ./data/train.dat \
  --labels  ./data/labels.dat \
  --test    ./data/test.dat \
  --param   1 \
  --tau     0.75 \
  --upsilon 0.34
