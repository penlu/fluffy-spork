#!/usr/bin/env bash

n=100
N=$1
step=$(($1 / 10))
MM=$(($1 * 5))
S=$(($1 * 100))
p=30

RESULT1=data/walksat_n${n}_N${N}_M${MM}_S${S}_p${p}.txt
RESULT2=data/walksat_n${n}_N${N}_M${MM}_S${S}_p${p}_steps.txt

echo "M sat unsat unknown" > $RESULT1
echo "M steps" > $RESULT2

run_sed () {
  sed -n -e 's/^walk: \([0-9]\+\) steps$/\1/p' out/$2.txt | tr -s "\n" " " >> $1;
}

export -f run_sed

for M in `seq 0 ${step} ${MM}`; do
  echo $M

  ./test_walk.sh -n ${n} -N ${N} -M $M -S ${S} -p ${p} > /dev/null;

  # get result counts
  COUNT_SAT=$(grep "^walk: sat$" out/* | wc -l);
  COUNT_UNS=$(grep "^walk: unsat$" out/* | wc -l);
  COUNT_UNK=$(grep "^walk: unknown$" out/* | wc -l);
  echo "$M $COUNT_SAT $COUNT_UNS $COUNT_UNK" >> $RESULT1;

  # get number of steps
  echo -n "$M " >> $RESULT2
  seq 1 ${n} | parallel -j32 run_sed $RESULT2
  echo >> $RESULT2
done
