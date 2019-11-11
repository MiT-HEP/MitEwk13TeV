#! /bin/bash


# SYSDIR=("MuSIT" "MuSta" "MuStaDirect")
SYSDIR=("MuSIT")

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance
EFFDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results"
EFFSYSDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics/"

#
# W->munu
# #
# for SYS in ${SYSDIR}; do
  # Sta=${EFFSYSDIR}/${SYS}/SysUnc_MuStaEff.root
  # SIT=${EFFSYSDIR}/${SYS}/SysUnc_MuSITEff.root
  # root -l -q computeAccSelWm_Sys.C+\(\"w0j_13.conf\",\"${EFFDIR}/Zmm/\",\"Wmp0Sel13TeV_${SYS}\",1,0,\"${SIT}\",\"${Sta}\",1\)
  # root -l -q computeAccSelWm_Sys.C+\(\"w1j_13.conf\",\"${EFFDIR}/Zmm/\",\"Wmp1Sel13TeV_${SYS}\",1,0,\"${SIT}\",\"${Sta}\",1\)
  # root -l -q computeAccSelWm_Sys.C+\(\"w2j_13.conf\",\"${EFFDIR}/Zmm/\",\"Wmp2Sel13TeV_${SYS}\",1,0,\"${SIT}\",\"${Sta}\",1\)
  
  # # root -l -q computeAccSelWm_Sys.C+\(\"w0j_13.conf\",\"${EFFDIR}/Zmm/\",\"Wmm0Sel13TeV_${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",1\)
  # # root -l -q computeAccSelWm_Sys.C+\(\"w1j_13.conf\",\"${EFFDIR}/Zmm/\",\"Wmm1Sel13TeV_${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",1\)
  # # root -l -q computeAccSelWm_Sys.C+\(\"w2j_13.conf\",\"${EFFDIR}/Zmm/\",\"Wmm2Sel13TeV_${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",1\)
# done

root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wep13TeV_GenAcc_Dressed\",1,1\)
root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wem13TeV_GenAcc_Dressed\",1,-1\)
root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/We13TeV_GenAcc_Dressed\",1,0\)

# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wmp13TeV_GenAcc_Dressed\",1,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wmm13TeV_GenAcc_Dressed\",1,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wm13TeV_GenAcc_Dressed\",1,0\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wmp13TeV_GenAcc_Undressed\",0,1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wmm13TeV_GenAcc_Undressed\",0,-1\)
# root -l -q computeAccGenWm_Sys.C+\(\"wj_13.conf\",\"${OUTDIR}/Wm13TeV_GenAcc_Undressed\",0,0\)

# root -l -q computeAccGenWe_Sys.C+\(\"wj_13.conf\",\"WeGen13TeV_Dressed_025\",1,0\)
