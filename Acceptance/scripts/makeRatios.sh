#!/bin/sh

s=$1
fsr=$2
mc=$3
DIR=$4
IN=$5

echo "Making Ratios for ${TeV}"

NUMS1=("wep"   "wep"  "wem"  "we"   "wmp"   "wmp"  "wmm"  "wm"  )
DNMS1=("wem"   "zee"  "zee"  "zee"  "wmm"   "zmm"  "zmm"  "zmm" )
RATIO=("wexpm" "wpez" "wmez" "wez"  "wmxpm" "wpmz" "wmmz" "wmz" )


for ((i=0;i<${#RATIO[@]};++i)); do 
  R=${RATIO[i]}
  FILE1=${DIR}/${IN}/ACC_${NUMS1[i]}_${s}_${mc}_${fsr}.txt
  FILE2=${DIR}/${IN}/ACC_${DNMS1[i]}_${s}_${mc}_${fsr}.txt
  OFILE=${DIR}/ACC_${RATIO[i]}_${s}_${mc}_${fsr}.txt
  
  echo "....  ${NUMS1[i]} / ${DNMS1[i]} " 

  echo -n "acc "  >> ${OFILE}
  n1=$(grep -e acc ${FILE1} | awk '{num=$2; printf ("%.5f",num)} ')
  n2=$(grep -e acc ${FILE2} | awk '{num=$2; printf ("%.5f",num)} ')
  echo "scale=8 ; ${n1}/${n2} " | bc >> ${OFILE}
  
  for ((k=0;k<100;++k)); do
    echo -n "pdf${k} "  >> ${OFILE}
    n1=$(grep -w pdf${k} ${FILE1} | awk '{num=$2; printf ("%.8f",num)} ')
    n2=$(grep -w pdf${k} ${FILE2} | awk '{num=$2; printf ("%.8f",num)} ')
    echo "scale=8 ; ${n1}/${n2} " | bc >> ${OFILE}
  done
  
  for ((k=0;k<6;++k)); do
    echo -n "qcd${k} "  >> ${OFILE}
    n1=$(grep -w qcd${k} ${FILE1} | awk '{num=$2; printf ("%.8f",num)} ')
    n2=$(grep -w qcd${k} ${FILE2} | awk '{num=$2; printf ("%.8f",num)} ')
    echo "scale=8 ; ${n1}/${n2} " | bc >> ${OFILE}
  done  
done