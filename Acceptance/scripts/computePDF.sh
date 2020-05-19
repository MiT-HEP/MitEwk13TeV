#!/bin/sh
s=$1
flav=$2
DIR=$3
INPUT=$4


echo "Computing PDF uncertainties for ${s}TeV ${flav} channel"
if [[ ${flav} == "ele" ]]; then
  FILES=("wep" "wem" "we" "wexpm" "zee" "wpez" "wmez" "wez")
elif [[ ${flav} == "mu" ]]; then
  FILES=("wmp" "wmm" "wm" "wmxpm" "zmm" "wpmz" "wmmz" "wmz")
fi
# only need to compare num_qcc and dnm_qcd for files

NPDFS=100
OFILE=${DIR}/PDF_${s}.txt

for ((i=0;i<${#FILES[@]};++i)); do 
  echo ".... ${FILES[i]}"
  FILE1=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}_amcPythia_dressed.txt
  sum=0.0
  echo -n "${FILES[i]} "  >> ${OFILE}
  
  a=$(grep acc ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
  
  for ((k=0;k<${NPDFS};++k)); do
    a1=$(grep -w pdf${k} ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
    diff=$(echo "scale=16 ; ${a}-${a1} " | bc)
    # echo "central value: " ${a}  " current value: " ${a1}  "  diff: " ${diff} 
    sum=$(echo "scale=16 ; ${sum} + ${diff}*${diff}" | bc)
    # echo "running total : " ${sum}
  done
  
  stddev=$(echo | awk -v x=${sum} -v n=${NPDFS} '{printf ("%.5f", sqrt(x))}')
  # echo "stddev " ${stddev}
  pct=$(echo | awk -v x=${a} -v s=${stddev} '{printf ("%.5f", 100*s/x)}')
  # echo "final " ${pct}
  echo ${pct} >> ${OFILE}
done