#!/bin/sh

s=$1
fsr=$2
TYPE=$3
DIR=$4


if [[ ${TYPE} == "powPhotos" || ${TYPE} == "powPythia" ]]; then
  XSECS=("11811" "8677")
elif [[ ${TYPE} == "amcPythia" &&  ${s} == "13" ]]; then
  echo "Set X-sections for W+jets"
  XSECS=("51410" "8417" "3306")
elif [[ ${TYPE} == "amcPythia" && ${s} == "5" ]]; then
    XSECS=("1")
fi

FILES=("wmp" "wmm" "wm" "wep" "wem" "we")

echo "Combining type ${TYPE}_${fsr}  at  ${s}TeV "

echo ${TOT}
NQCD=6
NPDF=100

## go through different combos
for ((j=0;j<${#FILES[@]};++j)); do
  echo ".... ${FILES[j]}"  
  TOT=0
  acc=0

  SMALLDIR="GEN_${FILES[j]}_${s}TeV_${TYPE}_${fsr}"
  OFILE=${DIR}/ACC_${FILES[j]}_${s}_${TYPE}_${fsr}.txt
  
  for ((i=0;i<${#XSECS[@]};++i)); do 
    XS=${XSECS[i]}
    FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${TYPE}_${fsr}_${i}.txt
    a1=$(grep -w acc ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
    if [[ 1 -eq "$(echo "${a1} > 0" | bc)" ]]; then 
      TOT=$(( ${TOT} + ${XS} ))
    fi
    acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
  done
  echo "acc" $(echo "scale=8 ; ${acc}/${TOT}" | bc) >> ${OFILE}

  
  for ((k=0;k<${NPDF};++k)); do
    acc=0
    for ((i=0;i<${#XSECS[@]};++i)); do 
      XS=${XSECS[i]}
      FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${TYPE}_${fsr}_${i}.txt
      a1=$(grep -w pdf${k} ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
      acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
    done
    echo "pdf${k}" $(echo "scale=8 ; ${acc}/${TOT}" | bc) >> ${OFILE}
  done
    
  for ((k=0;k<${NQCD};++k)); do
    acc=0
    for ((i=0;i<${#XSECS[@]};++i)); do 
      XS=${XSECS[i]}
      FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${TYPE}_${fsr}_${i}.txt
      a1=$(grep -w qcd${k} ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
      acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
    done
    echo "qcd${k}" $(echo "scale=8 ; ${acc}/${TOT}" | bc)>> ${OFILE}
  done
done