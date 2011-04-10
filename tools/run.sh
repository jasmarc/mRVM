#!/bin/bash
VERBOSITY=${1-1}

# valgrind \
#   --track-origins=yes \
#   --leak-check=full \
#   --show-reachable=yes \
./bin/mRVM \
  -k GAUSSIAN \
  -v $VERBOSITY \
  --train   ./data/iris.train.dat \
  --labels  ./data/iris.labels.dat \
  --test    ./data/iris.test.dat \
  --answers ./data/iris.answers.dat \
  --param   2 \
  --tau     0.0001 \
  --upsilon 0.0001
