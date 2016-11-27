#!/bin/bash
workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Inclusive_C.so")
computeAccSelWm_Inclusive("${workdir}/Acceptance/we_local.conf", "${filedir}/", "${filedir}/Results/Wm/Central", 0, 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Inclusive_Sys_C.so")
computeAccSelWm_Inclusive_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/", "${filedir}/Results/Wm/SigUp", 0, 1, "${filedir}/Results/EffShapeSysToys/MuSITSigPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuSITSigNegaSys.root", "${filedir}/Results/EffShapeSysToys/MuStaSigPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuStaSigNegaSys.root")
computeAccSelWm_Inclusive_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/", "${filedir}/Results/Wm/BkgUp", 0, 1, "${filedir}/Results/EffShapeSysToys/MuSITBkgPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuSITBkgNegaSys.root", "${filedir}/Results/EffShapeSysToys/MuStaBkgPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuStaBkgNegaSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Inclusive_Bin_C.so")
computeAccSelWm_Inclusive_Bin("${workdir}/Acceptance/we_local.conf", "${filedir}/", "${filedir}/Results/Wm/Bin", 0, 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Inclusive_PileupSys_C.so")
computeAccSelWm_Inclusive_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/", "${filedir}/Results/Wm/PileupUp", 0, 1, "puWeightsUp")
computeAccSelWm_Inclusive_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/", "${filedir}/Results/Wm/PileupDown", 0, 1, "puWeightsDown")

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Charge_C.so")
computeAccSelWm_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmp/Central", +1, 1)
computeAccSelWm_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmm/Central", -1, 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Sys_C.so")
computeAccSelWm_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmp/SigUp", +1, 1, "${filedir}/Results/EffShapeSysToys/MuSITSigPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuStaSigPosiSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmm/SigUp", -1, 1, "${filedir}/Results/EffShapeSysToys/MuSITSigNegaSys.root", "${filedir}/Results/EffShapeSysToys/MuStaSigNegaSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmp/BkgUp", +1, 1, "${filedir}/Results/EffShapeSysToys/MuSITBkgPosiSys.root", "${filedir}/Results/EffShapeSysToys/MuStaBkgPosiSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmm/BkgUp", -1, 1, "${filedir}/Results/EffShapeSysToys/MuSITBkgNegaSys.root", "${filedir}/Results/EffShapeSysToys/MuStaBkgNegaSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Bin_C.so")
computeAccSelWm_Bin("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmp/Bin", +1, 1)
computeAccSelWm_Bin("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmm/Bin", -1, 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_PileupSys_C.so")
computeAccSelWm_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmp/PileupUp", +1, 1, "puWeightsUp")
computeAccSelWm_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmm/PileupUp", -1, 1, "puWeightsUp")
computeAccSelWm_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmp/PileupDown", +1, 1, "puWeightsDown")
computeAccSelWm_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wmm/PileupDown", -1, 1, "puWeightsDown")

.q
EOF
