# This little snippet can randomize a file.
# from http://www.commandlinefu.com/commands/view/2320/randomize-lines-in-a-file
# awk 'BEGIN{srand()}{print rand(),$0}' SOMEFILE | sort -n | cut -d ' ' -f2-

PREFIX=iris
LINES=`wc -l ./data/${PREFIX}.train.dat | awk '{print int(($1 + 1) / 10)}'`

cd ./data

# Clean up, just in case
rm train.a*
rm test.a*
rm labels.a*
rm answers.a*

# Create the labels and answers
split -l $LINES $PREFIX.labels.dat del.labels.
for i in del.labels.*; do
  cat `(ls del.labels.* | grep -v $i)` > ${i//del./}.dat;
  cp $i ${i//del.labels/answers}.dat;
done;
rm del.labels.*

# Create the train and test data
split -l $LINES $PREFIX.train.dat del.train.
for i in del.train.*; do
  cat `(ls del.train.* | grep -v $i)` > ${i//del./}.dat;
  cp $i ${i//del.train/test}.dat;
done;
rm del.train.*

for i in train.*.dat; do
  ../bin/mRVM \
    -k GAUSSIAN \
    -v 1 \
    --train   $i\
    --labels  ${i//train/labels}\
    --test    ${i//train/test}\
    --answers ${i//train/answers}\
    --param   2 \
    --tau     0.0001 \
    --upsilon 0.0001
done;

# Clean up
rm train.a*
rm test.a*
rm labels.a*
rm answers.a*

# Go back to the previous directory
cd -