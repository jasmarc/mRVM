#!/bin/bash
VERBOSITY=${1-1}

./bin/mRVM \
  -k GAUSSIAN \
  -v $VERBOSITY \
  -r ./data/train.dat \
  -l ./data/labels.dat \
  -t ./data/test.dat \
  -p 1
