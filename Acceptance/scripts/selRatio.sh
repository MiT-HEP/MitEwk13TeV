#!/bin/sh

s=$1
DIR=$2
IN=$3

echo "Making Ratios for ${s} TeV"

NUMS1=("wep"   "wep"  "wem"  "we"   "wmp"   "wmp"  "wmm"  "wm"  )
DNMS1=("wem"   "zee"  "zee"  "zee"  "wmm"   "zmm"  "zmm"  "zmm" )
RATIO=("wexpm" "wpez" "wmez" "wez"  "wmxpm" "wpmz" "wmmz" "wmz" )

TYPE=("fsr" "mc" "bkg" "tagpt")

for ((i=0;i<${#RATIO[@]};++i)); do 
  R=${RATIO[i]}
  FILE1=${DIR}/${IN}/ACC_${NUMS1[i]}_${s}.txt
  FILE2=${DIR}/${IN}/ACC_${DNMS1[i]}_${s}.txt
  OFILE=${DIR}/ACC_${RATIO[i]}_${s}.txt
  
  echo "....  ${NUMS1[i]} / ${DNMS1[i]} " 

  echo -n "acc "  >> ${OFILE}
  n1=$(grep -w acc ${FILE1} | awk '{num=$2; printf ("%.5f",num)} ')
  n2=$(grep -w acc ${FILE2} | awk '{num=$2; printf ("%.5f",num)} ')
  echo "scale=8 ; ${n1}/${n2} " | bc >> ${OFILE}
  
  echo -n "acc_stat "  >> ${OFILE}
  n1=$(grep -w acc_stat ${FILE1} | awk '{num=$2; printf ("%.5f",num)} ')
  n2=$(grep -w acc_stat ${FILE2} | awk '{num=$2; printf ("%.5f",num)} ')
  echo "scale=8 ; ${n1}/${n2} " | bc >> ${OFILE}
  
  for ((k=0;k<${#TYPE[@]};++k)); do
    typ=${TYPE[k]}
    echo -n "sel_${typ} "  >> ${OFILE}
    a1=$(grep -w sel_${typ} ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
    a2=$(grep -w sel_${typ} ${FILE2} | awk '{num=$2; printf ("%.6f",num)} ')
    echo "scale=8 ; ${a1}/${a2} " | bc >> ${OFILE}
  done
  
  if [[ ${NUMS1[i]} == "wmm" ||  ${NUMS1[i]} == "wmp" ||  ${NUMS1[i]} == "wm" ]]; then
    for ((k=0;k<${#TYPE[@]};++k)); do
      # echo "standalone"
      typ=${TYPE[k]}
      echo -n "sta_${typ} "  >> ${OFILE}
      a1=$(grep -w sta_${typ} ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
      a2=$(grep -w sta_${typ} ${FILE2} | awk '{num=$2; printf ("%.6f",num)} ')
      echo "scale=8 ; ${a1}/${a2} " | bc >> ${OFILE}
    done
  fi
done