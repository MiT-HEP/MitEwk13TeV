#!/bin/sh

s=$1
TYPE=$2
DIR=$3


if [[ ${TYPE} == "powPhotos" || ${TYPE} == "powPythia" ]]; then
  XSECS=("11811" "8677")
elif [[ ${TYPE} == "" &&  ${s} == "13" ]]; then
  echo "Set X-sections for W+jets"
  XSECS=("51410" "8417" "3306")
elif [[ ${TYPE} == "" && ${s} == "5" ]]; then
    XSECS=("1")
fi

FILES=("wmp" "wmm" "wm" "wep" "wem" "we")
# FILES=("wmp")
TYPE=("fsr" "mc" "bkg" "tagpt")

echo "Combining  W+Jets  at  ${s}TeV "

echo ${TOT}
NQCD=6
NPDF=100

## go through different combos
for ((j=0;j<${#FILES[@]};++j)); do
  echo ".... ${FILES[j]}"  
  TOT=0
  acc=0

  SMALLDIR="SEL_${FILES[j]}_${s}TeV"
  OFILE=${DIR}/ACC_${FILES[j]}_${s}.txt
  
  for ((i=0;i<${#XSECS[@]};++i)); do 
    XS=${XSECS[i]}
    FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${i}.txt
    a1=$(grep -w acc ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
    if [[ 1 -eq "$(echo "${a1} > 0" | bc)" ]]; then 
      TOT=$(( ${TOT} + ${XS} ))
    fi
    acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
  done
  echo "acc" $(echo "scale=8 ; ${acc}/${TOT}" | bc) >> ${OFILE}
  
  acc=0
  a1=0
  for ((i=0;i<${#XSECS[@]};++i)); do 
    XS=${XSECS[i]}
    FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${i}.txt
    a1=$(grep -w acc_stat ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
    acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
  done
  echo "acc_stat" $(echo "scale=8 ; ${acc}/${TOT}" | bc) >> ${OFILE}
  
  for ((k=0;k<${#TYPE[@]};++k)); do
    typ=${TYPE[k]}
    acc=0
    for ((i=0;i<${#XSECS[@]};++i)); do 
      XS=${XSECS[i]}
      FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${i}.txt
      a1=$(grep -w sel_${typ} ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
      acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
    done
    echo "sel_${typ}" $(echo "scale=8 ; ${acc}/${TOT}" | bc) >> ${OFILE}
  done
  
  if [[ ${FILES[j]} == "wmm" ||  ${FILES[j]} == "wmp" ||  ${FILES[j]} == "wm" ]]; then
    for ((k=0;k<${#TYPE[@]};++k)); do
      typ=${TYPE[k]}
      acc=0
      for ((i=0;i<${#XSECS[@]};++i)); do 
        XS=${XSECS[i]}
        FILE1=${DIR}/${SMALLDIR}/ACC_${FILES[j]}_${s}_${i}.txt
        a1=$(grep -w sta_${typ} ${FILE1} | awk '{num=$2; printf ("%.6f",num)} ')
        acc=$(echo "scale=8 ; ${acc} + ${XS}*${a1}" | bc)
      done
      echo "sta_${typ}" $(echo "scale=8 ; ${acc}/${TOT}" | bc) >> ${OFILE}
    done
  fi

done