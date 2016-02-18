#!/bin/bash
#Mu HLT 0
#Mu Sta 5
#Mu SelIsoTrk 8

lepton=$1
efftype=$2
effnum=$3
#bkgmodel=$3
xaxisname=$4
yaxisname=$5

workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b << EOF
gSystem->Load("${workdir}/Selection/selectZmm_C.so")
selectZmm("${workdir}/Selection/zmm.conf", "${filedir}/${lepton}", 0)

.q
EOF
