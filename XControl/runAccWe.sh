#!/bin/bash
workdir="/afs/cern.ch/work/x/xniu/public/Lumi/Ele/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")
gSystem->Load("${workdir}/EleScale/EnergyScaleCorrection_class_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_Charge_C.so")
computeAccSelWe_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/Central", +1, 1, 1, 0)
computeAccSelWe_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/Central", -1, 1, 1, 0)

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_Sys_C.so")
computeAccSelWe_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/SigUp", +1, 1, 1, 0, "${filedir}/Results/EffShapeSysToys/EleGsfSelSigPosiSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/SigUp", -1, 1, 1, 0, "${filedir}/Results/EffShapeSysToys/EleGsfSelSigNegaSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/BkgUp", +1, 1, 1, 0, "${filedir}/Results/EffShapeSysToys/EleGsfSelBkgPosiSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/BkgUp", -1, 1, 1, 0, "${filedir}/Results/EffShapeSysToys/EleGsfSelBkgNegaSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_Bin_C.so")
computeAccSelWe_Bin("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/Bin", +1, 1, 1, 0)
computeAccSelWe_Bin("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/Bin", -1, 1, 1, 0)

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_PileupSys_C.so")
computeAccSelWe_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/PileupUp", +1, 1, 1, 0, "puWeightsUp")
computeAccSelWe_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/PileupUp", -1, 1, 1, 0, "puWeightsUp")
computeAccSelWe_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/PileupDown", +1, 1, 1, 0, "puWeightsDown")
computeAccSelWe_PileupSys("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/PileupDown", -1, 1, 1, 0, "puWeightsDown")

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_Charge_C.so")
computeAccSelWe_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/ScaleUp", +1, 1, 1, 1)
computeAccSelWe_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/ScaleUp", -1, 1, 1, 1)
computeAccSelWe_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/ScaleDown", +1, 1, 1, -1)
computeAccSelWe_Charge("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/ScaleDown", -1, 1, 1, -1)
.q
EOF

#charge mis id
#gSystem->Load("${workdir}/Acceptance/computeAccSelWe_CMI_C.so")
#computeAccSelWe_CMI("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wep/ChargeMisID", +1, 1, 1, 0)
#computeAccSelWe_CMI("${workdir}/Acceptance/we_local.conf", "${filedir}/Results/Wem/ChargeMisID", -1, 1, 1, 0)
