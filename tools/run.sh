#!/bin/bash
VERBOSITY=${1-1}

# valgrind --leak-check=yes \
# --track-origins=yes \
./bin/mRVM \
  -k GAUSSIAN \
  -v $VERBOSITY \
  --train   ./data/train.dat \
  --labels  ./data/labels.dat \
  --test    ./data/test.dat \
  --answers ./data/answers.dat \
  --param   1 \
  --tau     0.75 \
  --upsilon 0.34
