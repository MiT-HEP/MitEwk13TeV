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
if [[ ${TYPE} == "bkg" ]]; then 
  echo "Doing BKG uncertainty for ${s}TeV ${flav} channel"
  CAT="BKG"
elif [[ ${TYPE} == "tagpt" ]]; then
  echo "Doing Tag pT uncertainty for ${s}TeV ${flav} channel"
  CAT="TAG"
elif [[ ${TYPE} == "fsr" ]]; then
  echo "Doing FSR uncertainty for ${s}TeV ${flav} channel"
  CAT="FSR"
elif [[ ${TYPE} == "stat" ]]; then
  echo "Doing stat uncertainty for ${s}TeV ${flav} channel"
  CAT="STAT"
elif [[ ${TYPE} == "mc" ]]; then
  echo "Doing MC uncertainty for ${s}TeV ${flav} channel"
  CAT="MC"
fi


OFILE=${DIR}/${CAT}_${s}.txt

for ((i=0;i<${#FILES[@]};++i)); do 
  echo ".... ${FILES[i]}"  
  FILE1=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}.txt
  # FILE2=${DIR}/${INPUT}/ACC_${FILES[i]}_${s}.txt
  
  # echo "Using file " ${FILE1}
  # echo ${TYPE}
  echo -n "${FILES[i]} " >> ${OFILE}
  
  a1=$(grep -w acc ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
  a2=$(grep -w sel_${TYPE} ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
  if [[ ${TYPE} == "stat" ]]; then
    a2=$(grep -w acc_stat ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
  fi
  avg=${a2}
  a3=0
  if [[ ${flav} == "mu" && ${TYPE} != "stat" ]]; then
    # echo "hello" 
    a3=$(grep -w sta_${TYPE} ${FILE1} | awk '{x=$2; printf ("%.6f",x)} ')
    avg=$(echo "scale=8 ; (${a2}+${a3})/2" | bc)
    # echo ${avg}
  fi
  # echo ${a1} ${a2} ${a3} ${avg}
  ans1=$(echo "scale=8 ; ${a1}/${avg} - 1" | bc)
  ans=$(echo | awk -v x=${ans1} 'function abs(v) {return v < 0 ? -v : v} {printf ("%.4f",100*abs(x))} ')
  echo ${ans} >> ${OFILE} 
done