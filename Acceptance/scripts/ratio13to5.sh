#!/bin/sh

mc=$1
fsr=$2
s=$3

# NUMS1=("wep" "wep" "wem" "we")
# DNMS1=("wem" "zee" "zee" "zee")
# OUTFS=("wexpm" "wpez" "wmez" "wez")
NUMS1=("wep")
DNMS1=("wem")
RATIO=("wexpm")
DIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/FINAL_GEN"

NQCD=6
NPDF=100
# ACC_wem_w0_13_amcPythia_dressed
for ((i=0;i<${#RATIO[@]};++i)); do 
  R=${RATIO[i]}
  FILE1=${DIR}/ACC_${NUMS1[i]}_${s}_${mc}_${fsr}.txt
  FILE2=${DIR}/ACC_${DNMS1[i]}_${s}_${mc}_${fsr}.txt
  OFILE=${DIR}/ACC_${RATIO[i]}_${s}_${mc}_${fsr}.txt
  echo ${FILE1}
  # go line by line & do division
  
  echo -n "acc "  >> ${OFILE}
  a1=$(grep -e acc ${FILE1} | awk '{a=$2; printf ("%.6f",a)} ')
  a2=$(grep -e acc ${FILE2} | awk '{a=$2; printf ("%.6f",a)} ')
  echo "scale=8 ; ${a1}/${a2} " | bc >> ${OFILE}
  
  for ((k=0;k<${NPDF};++k)); do
    echo -n "${TYPE[j]}_pdf${k} "  >> ${OFILE}
    a1=$(grep -e pdf${k} ${FILE1} | awk '{a=$2; printf ("%.6f",a)} ')
    a2=$(grep -e pdf${k} ${FILE2} | awk '{a=$2; printf ("%.6f",a)} ')
    echo "scale=8 ; ${a1}/${a2} " | bc >> ${OFILE}
  done
  
  for ((k=0;k<${NQCD};++k)); do
    echo -n "${TYPE[j]}_qcd${k} "  >> ${OFILE}
    a1=$(grep -e qcd${k} ${FILE1} | awk '{a=$2; printf ("%.6f",a)} ')
    a2=$(grep -e qcd${k} ${FILE2} | awk '{a=$2; printf ("%.6f",a)} ')
    echo "scale=8 ; ${a1}/${a2} " | bc >> ${OFILE}
  done
  
done