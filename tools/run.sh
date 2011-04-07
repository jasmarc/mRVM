#!/bin/bash
VERBOSITY=${1-1}

# valgrind \
#   --track-origins=yes \
#   --leak-check=full \
#   --show-reachable=yes \
#   --suppressions=./tools/darwin9.supp \
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
