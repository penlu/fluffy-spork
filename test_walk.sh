usage() { echo "Usage: $0 [-n [how many tests]] [-N [how many vars]] [-M [how many clauses]] [-S [how many walksat steps]] [-p [walksat noise level 0-100]]"; exit 0; }

n=1000
N=1000
M=3800
S=10000
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

    S)
      S=${OPTARG}
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

for i in `seq 1 $n`; do
  ./randsat 3 $N $M | ./walk $S $p > out/$i.txt;
  if [ $? -ne 0 ]; then
    echo $?;
  fi;
done

grep "unsat after" out/* | wc -l
