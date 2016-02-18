#!/bin/bash
lepton=$1
efftype=$2
binvar=$3
binnum=$4
toynum=$5

workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

CMSSW_BASE="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/"
TOP="$PWD"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP


root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/makePseudoData_C.so")
makePseudoData("${filedir}/${lepton}${efftype}Eff/CB/plots/", "${filedir}/${lepton}${efftype}Eff/MG/plots/", "${binvar}_${binnum}", "${filedir}/${lepton}${efftype}Eff/Step2Output/CB/",-1,${binnum},${toynum})
.q
EOF

for ((psenum=0; psenum<${toynum};psenum++)); do
root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/TagAndProbe/doStep3_C.so")
doStep3("${filedir}/${lepton}${efftype}Eff/Step2Output/CB","${binvar}_${binnum}_${psenum}.dat","${filedir}/${lepton}${efftype}Eff/MG/plots/${binvar}_${binnum}.root", "${filedir}/${lepton}${efftype}Eff","${binvar}_${binnum}")
.q
EOF
rm ${filedir}/${lepton}${efftype}Eff/Step2Output/CB/${binvar}_${binnum}_${psenum}.dat
done
