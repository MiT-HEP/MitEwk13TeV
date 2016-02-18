#!/bin/bash
#Mu HLT 0
#Mu Sta 5
#Mu SelIsoTrk 8

lepton=$1
efftype=$2
#effnum=$3
bkgmodel=$3
#xaxisname=$4
#yaxisname=$5

workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b << EOF

gSystem->Load("${workdir}/Efficiency/FitBkg/RooVoigtianShape_cc.so")
gSystem->Load("${workdir}/Efficiency/FitBkg/RooCMSShape_cc.so")
gSystem->Load("${workdir}/Efficiency/FitBkg/CPlot_cc.so")
gSystem->Load("${workdir}/Efficiency/FitBkg/MitStyleRemix_cc.so")
gSystem->Load("${workdir}/Efficiency/FitBkg/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Efficiency/FitBkg/CEffUser2D_cc.so")
gSystem->Load("${workdir}/Efficiency/FitBkg/plotEff_C.so")
plotEff("${workdir}/Efficiency/FitBkg/${lepton}pteta.bins",2,7,2,7,"${filedir}/${lepton}GsfSelEff/Data/probes.root","${filedir}/${lepton}GsfSelEff/MGFitBkg","png",0,0,0,"${lepton}","GsfSel",0.7,1.02,2263,"${filedir}/${lepton}GsfSelEff/MC/probes.root")
.q
EOF
