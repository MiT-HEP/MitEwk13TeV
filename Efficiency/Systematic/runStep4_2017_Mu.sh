#!/bin/bash
#./runstep4.sh Ele GsfSel etapt
# lepton=$1
# efftype=$2
# charge=$3
binvar=etapt
#binnum=0
lepton=Mu
efftype=SIT
charge=Negative
dirname=Zmm
workdir="/afs/cern.ch/work/s/sabrandt/public/SM/LowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/TagAndProbe"
filedir="/afs/cern.ch/work/s/sabrandt/public/SM/LowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2"
plotdir="/afs/cern.ch/work/s/sabrandt/public/SM/LowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/Systematic"
POSTFIX=_POWxPythia_v1
POSTFIX_ALT=_POWxPhotos_v1
# POSTFIX_ALT=_POWxPythia_v1
INPUT_MAIN=${filedir}/results/${dirname}/Data/${lepton}${efftype}Eff${POSTFIX}/${charge}
INPUT_ALT=${filedir}/results/${dirname}/Data/${lepton}${efftype}Eff${POSTFIX}/${charge}
INPUT_TOY=${filedir}/results/${lepton}${efftype}Eff${POSTFIX}${POSTFIX_ALT}_origtest/${charge}/
outputdir=${plotdir}/Results/${lepton}${efftype}Eff${POSTFIX}${POSTFIX_ALT}_testclosure/${charge}
# for ((binnum=0; binnum<64; binnum++))
for ((binnum=1; binnum<2; binnum++))
#for binnum in 0 2 3 5 6 8
do

mkdir -p ${outputdir}


effcb=`grep eff ${INPUT_ALT}/plots/fitres${binvar}_${binnum}.txt | awk '{print $3}'`
# effpl=`grep eff ${filedir}/${dirname}/Data/${lepton}${efftype}Eff${POSTFIX_ALT}/${charge}/plots/fitres${binvar}_${binnum}.txt | awk '{print $3}'`

echo ${effcb}
# echo ${effpl}


grep EFFICIENCY ${INPUT_TOY}/_${binvar}_${binnum} | awk '{diff = $2 - "'"$effcb"'"; pull = diff/$3; print pull }' > ${INPUT_TOY}/Sig_${binvar}_${binnum}.txt
# grep EFFICIENCY ${INPUT_TOY}/_${binvar}_${binnum} | awk '{diff = $2; pull = diff/$3; print pull }' > ${INPUT_TOY}/Sig_${binvar}_${binnum}.txt
# grep EFFICIENCY ${filedir}/${lepton}${efftype}Eff${POSTFIX_ALT}/_pl${binvar}_${binnum} | awk '{diff = $2 - "'"$effpl"'"; pull = diff/$3; print pull }' > ${filedir}/${lepton}${efftype}Eff/Bkg_${binvar}_${binnum}.txt

# #grep EFFICIENCY ${filedir}/${lepton}${efftype}Eff/${binvar}_${binnum} | awk '$2 ~ /^[0-9]\.[0-9]*$/' | awk '$3 ~ /^[0-9]\.[0-9]*$/' | awk '{diff = $2 - "'"$effcb"'"; pull = diff/$3; print pull }' > ${filedir}/${lepton}${efftype}Eff/Sig_${binvar}_${binnum}.txt
# #grep EFFICIENCY ${filedir}/${lepton}${efftype}Eff/PL/pl${binvar}_${binnum} | awk '$2 ~ /^[0-9]\.[0-9]*$/' | awk '$3 ~ /^[0-9]\.[0-9]*$/' | awk '{diff = $2 - "'"$effpl"'"; pull = diff/$3; print pull }' > ${filedir}/${lepton}${efftype}Eff/Bkg_${binvar}_${binnum}.txt

root -l -b << EOF
gSystem->Load("${plotdir}/SDrawSystematic_C.so")
SDrawSystematic("${INPUT_TOY}/Sig_${binvar}_${binnum}.txt", "${outputdir}", "Sig_${lepton}_${efftype}_${charge}_${binvar}_${binnum}")
.q
EOF
#SDrawSystematic("${filedir}/${lepton}${efftype}Eff/Bkg_${binvar}_${binnum}.txt", "${plotdir}", "Bkg_${lepton}_${efftype}_${binvar}_${binnum}")
echo "done with root" 
effmg=`grep eff ${INPUT_MAIN}/plots/fitres${binvar}_${binnum}.txt | awk '{print $3}'`
sigmamg=`grep eff ${INPUT_MAIN}/plots/fitres${binvar}_${binnum}.txt | awk '{print $5}'`

if [ "$sigmamg" == "<none>" ];then
        sigmamg=`grep eff ${INPUT_MAIN}/plots/fitres${binvar}_${binnum}.txt | sed 's/(+/ /g' | sed 's/,-/ /g' | sed 's/)/ /g' | awk '{print $5}'`
fi
echo "hello" 
head -n 1 ${outputdir}/Sig_${lepton}_${efftype}_${charge}_${binvar}_${binnum}.txt | awk '{print $1 * "'"$sigmamg"'" / "'"$effmg"'"}' > ${outputdir}/Sig_Final_${lepton}_${efftype}_${charge}_${binvar}_${binnum}.txt
# head -n 1 ${outputdir}/Bkg_${lepton}_${efftype}_${charge}_${binvar}_${binnum}.txt | awk '{print $1 * "'"$sigmamg"'" / "'"$effmg"'"}' > ${plotdir}/Results/Bkg_${lepton}_${efftype}_${charge}_${binvar}_${binnum}.txt


done
