#!/bin/bash
workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelZmmBinned_Charge_C.so")
computeAccSelZmmBinned_Charge("${workdir}/Acceptance/zmm.conf", "${filedir}/", "${filedir}/Results/Zmm/Central/", 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelZmmBinned_Sys_C.so")
computeAccSelZmmBinned_Sys("${workdir}/Acceptance/zmm.conf", "${filedir}/", "${filedir}/Results/Zmm/SigUp", 1, "${filedir}/Results/EffShapeSysToys/MuSITSigPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuSITSigNegaSys.root", "${filedir}/Results/EffShapeSysToys/MuStaSigPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuStaSigNegaSys.root")
computeAccSelZmmBinned_Sys("${workdir}/Acceptance/zmm.conf", "${filedir}/", "${filedir}/Results/Zmm/BkgUp", 1, "${filedir}/Results/EffShapeSysToys/MuSITBkgPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuSITBkgNegaSys.root", "${filedir}/Results/EffShapeSysToys/MuStaBkgPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuStaBkgNegaSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelZmmBinned_Bin_C.so")
computeAccSelZmmBinned_Bin("${workdir}/Acceptance/zmm.conf", "${filedir}/Results/Zmm/Bin", 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelZmmBinned_PileupSys_C.so")
computeAccSelZmmBinned_PileupSys("${workdir}/Acceptance/zmm.conf", "${filedir}/", "${filedir}/Results/Zmm/PileupUp", 1, "puWeightsUp")
computeAccSelZmmBinned_PileupSys("${workdir}/Acceptance/zmm.conf", "${filedir}/", "${filedir}/Results/Zmm/PileupDown", 1, "puWeightsDown")

.q
EOF
