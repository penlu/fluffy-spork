#!/usr/bin/env bash

n=100
N=1000
step=200
MM=6000
SS=10000
SW=10000
p=30

# take lots of data...
RESULT1=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}.txt
RESULT2=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}_walksteps.txt
RESULT3=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}_walkstart.txt
RESULT4=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}_unsat.txt
RESULT5=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}_unconverged.txt
RESULT6=data/survey_n${n}_N${N}_M${MM}_SS${SS}_SW${SW}_p${p}_timedout.txt

echo "M sat unsat walkunsat unconverged unknown" > $RESULT1
echo "M walksteps" > $RESULT2
echo "M walkstart" > $RESULT3
echo "M unsat" > $RESULT4
echo "M unconverged" > $RESULT5
echo "M timedout" > $RESULT5

run_sed () {
  sed -n -e 's/^walk: \([0-9]\+\) steps$/\1/p' out/$6.txt | tr -s "\n" " " >> $1;
  sed -n -e 's/^survey: starting walk at \([0-9]\+\) steps$/\1/p' out/$6.txt | tr -s "\n" " " >> $2;
  sed -n -e 's/^survey: unsat after \([0-9]\+\) steps$/\1/p' out/$6.txt | tr -s "\n" " " >> $3;
  sed -n -e 's/^survey: unconverged after \([0-9]\+\) steps$/\1/p' out/$6.txt | tr -s "\n" " " >> $4;
  sed -n -e 's/^survey: \([0-9]\+\) steps: timed out$/\1/p' out/$6.txt | tr -s "\n" " " >> $5;
}

export -f run_sed

for M in `seq 0 ${step} ${MM}`; do
  echo $M

  ./test_survey.sh -n ${n} -N ${N} -M $M -SS ${SS} -SW ${SW} -p ${p} > /dev/null;

  # get result counts
  COUNT_SAT=$(grep "^survey: sat$" out/* | wc -l);
  COUNT_UNS=$(grep "^survey: unsat$" out/* | wc -l);
  COUNT_WUS=$(grep "^survey: walk unsat$" out/* | wc -l);
  COUNT_UNC=$(grep "^survey: unconverged$" out/* | wc -l);
  COUNT_UNK=$(grep "^walk: unknown$" out/* | wc -l);
  echo "$M $COUNT_SAT $COUNT_UNS $COUNT_WUS $COUNT_UNC $COUNT_UNK" >> $RESULT1;

  # get number of steps
  echo -n "$M " >> $RESULT2
  echo -n "$M " >> $RESULT3
  echo -n "$M " >> $RESULT4
  echo -n "$M " >> $RESULT5
  echo -n "$M " >> $RESULT6
  seq 1 ${n} | parallel -j32 run_sed $RESULT2 $RESULT3 $RESULT4 $RESULT5 $RESULT6
  echo >> $RESULT2
  echo >> $RESULT3
  echo >> $RESULT4
  echo >> $RESULT5
  echo >> $RESULT6
done
