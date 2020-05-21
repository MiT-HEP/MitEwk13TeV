#!/bin/sh
s=$1
flav=$2
DIR=$3
INPUT=$4

echo "Computing QCD uncertainties for ${s}TeV ${flav} channel"
if [[ ${flav} == "ele" ]]; then
  FILES=("wep" "wem" "we" "wexpm" "zee" "wpez" "wmez" "wez")
elif [[ ${flav} == "mu" ]]; then
  FILES=("wmp" "wmm" "wm" "wmxpm" "zmm" "wpmz" "wmmz" "wmz")
fi

OFILE=${DIR}/QCD_${s}.txt
NQCD=6

for ((i=0;i<${#FILES[@]};++i)); do 
  echo ".... ${FILES[i]}" 
  FILE1=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}_amcPythia_dressed.txt
  sum=0.0
  echo -n "${FILES[i]} "  >> ${OFILE}
  
  a=$(grep acc ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
  # echo "acceptance is "  ${a}
  max=0
  for ((k=0;k<${NQCD};++k)); do
    a1=$(grep -w qcd${k} ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
    diff=$(echo "scale=16 ; ${a}-${a1} " | bc)
    # echo diff ${diff}
    # echo "try max  with diff = " ${diff} " and max = " ${max} 
    max=$(echo | awk -v x=${diff} -v m=${max} '{if(x>m) m=x; if(-x>m) m=-x;} { print m ;}')
    # echo "max is ${max}"
  done
  
  # echo "a="${a} " max = "${max} 
  rel=$(echo | awk -v x=${max} -v y=${a} '{printf ("%.5f", 100*(x/y))}')
  # echo "pct diff is " ${rel}
  echo ${rel} >> ${OFILE}

done