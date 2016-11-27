#!/bin/bash
lepton=$1
efftype=$2
lepcharge=$3
binvar=$4
binnum=$5
toynum=$6

workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

CMSSW_BASE="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/"
TOP="$PWD"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP

#TOP="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/makePseudoData_C.so")
makePseudoData("${filedir}/${lepton}${efftype}Eff/MG${lepcharge}/plots/", "${filedir}/${lepton}${efftype}Eff/PL${lepcharge}/plots/", "${binvar}_${binnum}", "${TOP}/${lepton}${efftype}Eff/Step2Output/PL${lepcharge}/",${binnum},${binnum},${toynum})
.q
EOF

for ((psenum=0; psenum<${toynum};psenum++)); do
root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/doStep3_C.so")
doStep3("${TOP}/${lepton}${efftype}Eff/Step2Output/PL${lepcharge}","${binvar}_${binnum}_${psenum}.dat","${filedir}/${lepton}${efftype}Eff/MG${lepcharge}/plots/${binvar}_${binnum}.root","${filedir}/${lepton}${efftype}Eff/","${lepcharge}_pl${binvar}_${binnum}")
.q
EOF
rm ${TOP}/${lepton}${efftype}Eff/Step2Output/PL${lepcharge}/${binvar}_${binnum}_${psenum}.dat
done
