#!/bin/bash
workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelZmmBinned_C.so")
computeAccSelZmmBinned("${workdir}/Acceptance/zmm.conf", "${filedir}/Results/Zmm/Ori/woPU", 0)
computeAccSelZmmBinned("${workdir}/Acceptance/zmm.conf", "${filedir}/Results/Zmm/Ori/wPU", 1)
.q
EOF
