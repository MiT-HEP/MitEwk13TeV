#!/bin/sh

## ( ./makeTable.sh  outputName )

#ele or mu
s=$1
flav=$2
fname=$3
DIR=$4
IN=$5

if [[ ${flav} == "ele" ]]; then
  CATS=("wep" "wem" "we" "wexpm" "zee" "wpez" "wmez" "wez")
elif [[ ${flav} == "mu" ]]; then
  CATS=("wmp" "wmm" "wm" "wmxpm" "zmm" "wpmz" "wmmz" "wmz")
fi

CAPS=("FSR" "MC" "Background" "Tag pT" "Statistical" ) # add "Muon momentum scale "
TXTS=("FSR" "MC" "BKG" "TAG" "STAT")

## empty array to track total
SUM=(0 0 0 0 0 0 0 0) 
SF="13_out"


OUT=${DIR}/${flav}${fname}
## Fill this part in by hand, the script will take care of the math and add them into total
# THY=(0 0 0 0 0 0 0 0)

# ssww
for ((i=0;i<${#CAPS[@]};++i)); do 
  CAP=${CAPS[i]}
  echo -n ${CAP} >> ${OUT}
  for ((j=0;j<${#CATS[@]};++j)); do 
    FILE=${DIR}/${IN}/${TXTS[i]}_${s}.txt
    CAT=${CATS[j]}
    # cat ${FILE}
    # echo ${FILE} ${TXT}
    # set column delimiter
    echo -n " & " >> ${OUT}
    # grep for the identifier, get the #, convert to percent and round to 2 decimals
    UNC=$(grep -w ${CAT} ${FILE} | awk '{pct=$2; printf ("%.3f",pct)} ')
    # echo ${UNC}
    echo -n ${UNC} >> ${OUT}
    # echo ${SUM[j]} ${UNC}
    SUM[j]=$(echo ${SUM[j]} + ${UNC} \* ${UNC} | bc)
    # echo ${SUM[j]}
  done
  ## endline for latex table
  echo " \\\\ " >> ${OUT}
done

echo  " \\\\" >> ${OUT}
echo "\\hline \\hline" >> ${OUT}

echo -n "Total [\%]" >> ${OUT}
for   ((i=0;i<${#SUM[@]};++i)); do 
  echo -n " & " >> ${OUT}
  sqr=$(echo | awk -v x=${SUM[i]} '{printf ("%.2f", sqrt(x))}')
  echo -n ${sqr} >> ${OUT}
done
echo  " \\\\" >> ${OUT}


# THY=(0 0 0 0 0 0 0 0)
echo " " >> ${OUT}
echo "%%%% for the full table:" >> ${OUT}
echo -n "EFFS=(">> ${OUT}
for   ((i=0;i<${#SUM[@]};++i)); do 
  echo -n " \"" >> ${OUT}
  sqr=$(echo | awk -v x=${SUM[i]} '{printf ("%.2f", sqrt(x))}')
  echo -n ${sqr} >> ${OUT}
  echo -n "\"" >> ${OUT}
done
echo " )" >> ${OUT}