#!/bin/sh

s=$1
TYPE=$2
flav=$3
DIR=$4
INPUT=$5

if [[ ${flav} == "ele" ]]; then
  FILES=("wep" "wem" "we" "wexpm" "zee" "wpez" "wmez" "wez")
elif [[ ${flav} == "mu" ]]; then
  FILES=("wmp" "wmm" "wm" "wmxpm" "zmm" "wpmz" "wmmz" "wmz")
fi
# only need to compare num_qcc and dnm_qcd for files
if [[ ${TYPE} == "ewk" ]]; then 
  echo "Doing electroweak uncertainty for ${s}TeV ${flav} channel"
  CAT="EWK"
  mc1="powPhotos"
  mc2="powPythia"
  fsr="undressed"
elif [[ ${TYPE} == "res" ]]; then
  echo "Doing resummation uncertainty for ${s}TeV ${flav} channel"
  CAT="RES"
  mc1="amcPythia"
  mc2="ptWeight"
  fsr="dressed"
fi


OFILE=${DIR}/${CAT}_${s}.txt

for ((i=0;i<${#FILES[@]};++i)); do 
  echo ".... ${FILES[i]}"  
  FILE1=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}_${mc1}_${fsr}.txt
  # FILE2=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}_${mc2}_${fsr}.txt
  FILE2=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}_${mc1}_${mc2}.txt
  
  
  echo -n "${FILES[i]} " >> ${OFILE}
  
  a1=$(grep -w acc ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
  a2=$(grep -w acc ${FILE2} | awk '{x=$2; printf ("%.6f",x)} ')
  
  
  ans1=$(echo "scale=8 ; ${a1}/${a2} - 1" | bc)
  ans=$(echo | awk -v x=${ans1} 'function abs(v) {return v < 0 ? -v : v} {printf ("%.4f",100*abs(x))} ')
  echo ${ans} >> ${OFILE} 
done