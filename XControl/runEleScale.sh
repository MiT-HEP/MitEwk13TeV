#!/bin/bash
workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b << EOF
gSystem->Load("${workdir}/Utils/CPlot_cc.so")
gSystem->Load("${workdir}/Utils/MitStyleRemix_cc.so")
gSystem->Load("${workdir}/EleScale/EleScale_C.so")
EleScale()
.q
EOF
