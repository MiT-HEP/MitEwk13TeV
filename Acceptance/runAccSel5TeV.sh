#! /bin/bash


SYS=""
OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/FINAL_SEL/TEST
EFFDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_5TeV/results"

S13="13TeV"
S5="5TeV"

#### Efficiency Systematics Files
EFFSYSDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics"
Sta=${EFFSYSDIR}/SysUnc_MuStaEff.root
SIT=${EFFSYSDIR}/SysUnc_MuSITEff.root
GSF=${EFFSYSDIR}/SysUnc_EleGSFSelEff.root

# #
# #    W- > mu nu
# # #

  root -l -q computeAccSelWm.C+\(\"w5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_wmp_5TeV\",\"ACC_wmp_5\",1,0,\"${SIT}\",\"${Sta}\",0\)
  
  root -l -q computeAccSelWm.C+\(\"w5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_wmm_5TeV\",\"ACC_wmm_5\",-1,0,\"${SIT}\",\"${Sta}\",0\)
  
  root -l -q computeAccSelWm.C+\(\"w5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_wm_5TeV\",\"ACC_wm_5\",0,0,\"${SIT}\",\"${Sta}\",0\)

# ####
# ####     Z -> mu mu

  root -l -q computeAccSelZmm.C+\(\"z5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_zmm_5TeV${SYS}\",\"ACC_zmm_5\",0,\"${SIT}\",\"${Sta}\",0\)


# # #
# # # ##    W -> e nu
# # #


root -l -q computeAccSelWe.C+\(\"w5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_wep_5TeV${SYS}\",\"ACC_wep_5\",1,0,1,0,\"${GSF}\",0\)

root -l -q computeAccSelWe.C+\(\"w5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_wem_5TeV${SYS}\",\"ACC_wem_5\",-1,0,1,0,\"${GSF}\",0\)

root -l -q computeAccSelWe.C+\(\"w5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_we_5TeV${SYS}\",\"ACC_we_5\",0,0,1,0,\"${GSF}\",0\)


# # # #
# # # # Z->ee
# # # #
root -l -q computeAccSelZee.C+\(\"z5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_zee_5TeV${SYS}\",\"ACC_zee_5\",\"${S5}\",0,1,0,\"${GSF}\",0\)


# # rm *.so *.d
