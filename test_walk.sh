mkdir -p out
rm -f out/*.txt
for i in `seq 1 1000`; do
  ./randsat 3 ${1} ${2} | ./walk 1000 30 > out/$i.txt;
  if [ $? -ne 0 ]; then
    echo $?;
  fi;
done
grep "unsat after" out/* | wc -l
