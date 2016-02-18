#!/bin/bash
workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_C.so")
computeAccSelWm("${workdir}/Acceptance/wmp.conf", "${filedir}/Result/Wm/Ori", 0, 1)
computeAccSelWm("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wmp/Ori", +1, 1)
computeAccSelWm("${workdir}/Acceptance/wmm.conf", "${filedir}/Results/Wmm/Ori", -1, 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Sys_C.so")
computeAccSelWm_Sys("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wm/SigUp", 0, 1, "${filedir}/Results/MuSITSigSys.root", "${filedir}/Results/MuStaSigSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wmp/SigUp", +1, 1, "${filedir}/Results/MuSITSigSys.root", "${filedir}/Results/MuStaSigSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/wmm.conf", "${filedir}/Results/Wmm/SigUp", -1, 1, "${filedir}/Results/MuSITSigSys.root", "${filedir}/Results/MuStaSigSys.root")

computeAccSelWm_Sys("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wm/BkgUp", 0, 1, "${filedir}/Results/MuSITBkgSys.root", "${filedir}/Results/MuStaBkgSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wmp/BkgUp", +1, 1, "${filedir}/Results/MuSITBkgSys.root", "${filedir}/Results/MuStaBkgSys.root")
computeAccSelWm_Sys("${workdir}/Acceptance/wmm.conf", "${filedir}/Results/Wmm/BkgUp", -1, 1, "${filedir}/Results/MuSITBkgSys.root", "${filedir}/Results/MuStaBkgSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelWm_Bin_C.so")
computeAccSelWm_Bin("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wm/Bin", 0, 1)
computeAccSelWm_Bin("${workdir}/Acceptance/wmp.conf", "${filedir}/Results/Wmp/Bin", +1, 1)
computeAccSelWm_Bin("${workdir}/Acceptance/wmm.conf", "${filedir}/Results/Wmm/Bin", -1, 1)
.q
EOF
