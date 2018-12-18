usage() { echo "Usage: $0 [-n [how many tests]] [-N [how many vars]] [-M [how many clauses]] [-I [how many warning prop iters]]"; exit 0; }

n=1000
N=1000
M=3800
I=10000

while getopts ":n:N:M:I:" o; do
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

    I)
      I=${OPTARG}
    ;;

    *)
      usage
    ;;
  esac
done

mkdir -p out
rm -f out/*.txt

for i in `seq 1 $n`; do
  ./randsat 3 $N $M | ./warn $I > out/$i.txt;
  if [ $? -ne 0 ]; then
    echo $?;
  fi;
done

grep "unconverged" out/* | wc -l | tr -s "\n" " "
grep "unsat" out/* | wc -l
grep "UNSAT" out/* | wc -l
