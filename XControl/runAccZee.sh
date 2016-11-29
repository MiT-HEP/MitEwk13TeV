#!/bin/bash
workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")
gSystem->Load("${workdir}/EleScale/EnergyScaleCorrection_class_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_Charge_C.so")
computeAccSelZeeBinned_Charge("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/Central/", 1, 1, 0)

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_Sys_C.so")
computeAccSelZeeBinned_Sys("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/SigUp", 1, 1, 0, "${filedir}/Results/EffShapeSysToys/EleGsfSelSigPosiSys.root", "${filedir}/Results/EffShapeSysToys/EleGsfSelSigNegaSys.root")
computeAccSelZeeBinned_Sys("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/BkgUp", 1, 1, 0, "${filedir}/Results/EffShapeSysToys/EleGsfSelBkgPosiSys.root", "${filedir}/Results/EffShapeSysToys/EleGsfSelBkgNegaSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_Bin_C.so")
computeAccSelZeeBinned_Bin("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/Bin", 1, 1, 0)

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_PileupSys_C.so")
computeAccSelZeeBinned_PileupSys("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/PileupUp", 1, 1, 0, "puWeightsUp")
computeAccSelZeeBinned_PileupSys("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/PileupDown", 1, 1, 0, "puWeightsDown")

gSystem->Load("${workdir}/Acceptance/computeAccSelZeeBinned_Charge_C.so")
computeAccSelZeeBinned_Charge("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/ScaleUp", 1, 1, 1
computeAccSelZeeBinned_Charge("${workdir}/Acceptance/zee.conf", "${filedir}/", "${filedir}/Results/Zee/ScaleDown", 1, 1, -1)
.q
EOF
