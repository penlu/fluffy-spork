#!/usr/bin/env bash

n=100
N=$1
step=$(($1 / 10))
MM=$(($1 * 5))
S=$(($1 * 100))
p=30

DATADIR=data/walksat_n${n}_N${N}_M${MM}_S${S}_p${p}
rm -rf $DATADIR
mkdir -p $DATADIR
RESULT=$DATADIR/outcomes.txt

echo "M sat unsat unknown" > $RESULT

for M in `seq 0 ${step} ${MM}`; do
  echo $M

  ./test_walk.sh -n ${n} -N ${N} -M $M -S ${S} -p ${p} > /dev/null;

  # get result counts
  COUNT_SAT=$(grep "^walk: sat$" out/* | wc -l);
  COUNT_UNS=$(grep "^walk: unsat$" out/* | wc -l);
  COUNT_UNK=$(grep "^walk: unknown$" out/* | wc -l);
  echo "$M $COUNT_SAT $COUNT_UNS $COUNT_UNK" >> $RESULT;

  mkdir -p $DATADIR/out_$M
  mv out/* $DATADIR/out_$M
done
