#!/bin/bash
workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_C.so")
computeAccSelZeeBinned("${workdir}/Acceptance/zee.conf", "${filedir}/Results/Zee/Ori/", 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_Sys_C.so")
computeAccSelZeeBinned_Sys("${workdir}/Acceptance/zee.conf", "${filedir}/Results/Zee/SigUp", 1, "${filedir}/Results/EleGsfSelSigSys.root")
computeAccSelZeeBinned_Sys("${workdir}/Acceptance/zee.conf", "${filedir}/Results/Zee/BkgUp", 1, "${filedir}/Results/EleGsfSelBkgSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_Bin_C.so")
computeAccSelZeeBinned_Bin("${workdir}/Acceptance/zee.conf", "${filedir}/Results/Zee/Bin", 1)

.q
EOF
