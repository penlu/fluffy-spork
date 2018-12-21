#!/usr/bin/env bash

n=$2
N=$1
step=$(($N / 10))
MM=$(($N * 5))
SS=1000
SW=$(($N * 100))
p=30

# take lots of data...
DATADIR=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}
rm -rf $DATADIR
mkdir -p $DATADIR
RESULT=$DATADIR/outcomes.txt

echo "M sat unsat walkunsat unconverged unknown" > $RESULT

for M in `seq 0 $step $MM`; do
  echo $M

  ./test_survey.sh -n ${n} -N ${N} -M $M -S ${SS} -W ${SW} -p ${p} > /dev/null;

  # get result counts
  COUNT_SAT=$(grep "^survey: sat$" out/* | wc -l);
  COUNT_UNS=$(grep "^survey: unsat$" out/* | wc -l);
  COUNT_WUS=$(grep "^survey: walk unsat$" out/* | wc -l);
  COUNT_UNC=$(grep "^survey: unconverged$" out/* | wc -l);
  COUNT_UNK=$(grep "^survey: unknown$" out/* | wc -l);
  echo "$M $COUNT_SAT $COUNT_UNS $COUNT_WUS $COUNT_UNC $COUNT_UNK" >> $RESULT;

  mkdir -p $DATADIR/out_$M
  mv out/* $DATADIR/out_$M
done
