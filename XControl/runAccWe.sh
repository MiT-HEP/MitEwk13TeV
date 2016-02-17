#!/bin/bash
workdir="/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

root -l -b <<EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_C.so")
computeAccSelWe("${workdir}/Acceptance/wep.conf", "${filedir}/Result/We/Ori", 0, 1)
computeAccSelWe("${workdir}/Acceptance/wep.conf", "${filedir}/Results/Wep/Ori", +1, 1)
computeAccSelWe("${workdir}/Acceptance/wem.conf", "${filedir}/Results/Wem/Ori", -1, 1)

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_Sys_C.so")
computeAccSelWe_Sys("${workdir}/Acceptance/wep.conf", "${filedir}/Results/We/SigUp", 0, 1, "${filedir}/Results/EleGsfSelSigSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/wep.conf", "${filedir}/Results/Wep/SigUp", +1, 1, "${filedir}/Results/EleGsfSelSigSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/wem.conf", "${filedir}/Results/Wem/SigUp", -1, 1, "${filedir}/Results/EleGsfSelSigSys.root")

computeAccSelWe_Sys("${workdir}/Acceptance/wep.conf", "${filedir}/Results/We/BkgUp", 0, 1, "${filedir}/Results/EleGsfSelBkgSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/wep.conf", "${filedir}/Results/Wep/BkgUp", +1, 1, "${filedir}/Results/EleGsfSelBkgSys.root")
computeAccSelWe_Sys("${workdir}/Acceptance/wem.conf", "${filedir}/Results/Wem/BkgUp", -1, 1, "${filedir}/Results/EleGsfSelBkgSys.root")

gSystem->Load("${workdir}/Acceptance/computeAccSelWe_Bin_C.so")
computeAccSelWe_Bin("${workdir}/Acceptance/wep.conf", "${filedir}/Results/We/Bin", 0, 1)
computeAccSelWe_Bin("${workdir}/Acceptance/wep.conf", "${filedir}/Results/Wep/Bin", +1, 1)
computeAccSelWe_Bin("${workdir}/Acceptance/wem.conf", "${filedir}/Results/Wem/Bin", -1, 1)
.q
EOF
