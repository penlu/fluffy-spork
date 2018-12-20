#!/usr/bin/env bash

usage() { echo "Usage: $0 [-n [how many tests]] [-N [how many vars]] [-M [how many clauses]] [-SS [how many survey iters]] [-SW [how many walksat steps]] [-p [walksat noise level 0-100]]"; exit 0; }

n=1000
N=1000
M=3800
SS=10000
SW=10000
p=30

while getopts ":n:N:M:S:p:" o; do
  case "${o}" in
    n)
      n=${OPTARG}
    ;;

    N)
      N=${OPTARG}
    ;;

    M)
      M=${OPTARG}
    ;;

    SS)
      SS=${OPTARG}
    ;;

    SW)
      SW=${OPTARG}
    ;;

    p)
      p=${OPTARG}
    ;;

    *)
      usage
    ;;
  esac
done

mkdir -p out
rm -f out/*.txt

run_test () {
  N=$1
  M=$2
  SS=$3
  SW=$4
  p=$5
  i=$6
  ./randsat 3 $N $M | ./survey $SS $SW $p > out/$i.txt;
  if [ $? -ne 0 ]; then
    echo $?;
  fi;
}

export -f run_test

seq 1 $n | parallel -j8 run_test $N $M $SS $SW $p

grep "unknown" out/* | wc -l
