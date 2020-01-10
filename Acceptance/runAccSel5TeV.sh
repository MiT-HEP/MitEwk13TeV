#! /bin/bash

SYSDIR=""

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/AccSel_5TeV_EleMedID
EFFDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_5TeV/results"
EFFSYSDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics"

S5="5TeV"

# #
# # W->munu
# # #
# # for SYS in SYSDIR; do
  Sta=${EFFSYSDIR}/SysUnc_MuStaEff.root
  SIT=${EFFSYSDIR}/SysUnc_MuSITEff.root
  # root -l -q computeAccSelWm_Sys.C+\(\"w_5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/PRELIM_Wmp_${S5}_wEff${SYS}\",1,0,\"${SIT}\",\"${Sta}\",0\)
  
  # root -l -q computeAccSelWm_Sys.C+\(\"w_5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/PRELIM_Wmm_${S5}_wEff${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",0\)
# done

#
# Z->mumu

  root -l -q computeAccSelZmmBinned_Sys.C+\(\"z_5.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/PRELIM_Zmm_${S5}_full\",0,\"${SIT}\",\"${Sta}\",0\)


# #
# # W->enu
# #
GSF=${EFFSYSDIR}/SysUnc_EleGSFSelEff.root
# root -l -q computeAccSelWe_Sys.C+\(\"w_5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wep_${S5}_wEff\",1,0,1,0,\"${GSF}\",1\)
# root -l -q computeAccSelWe_Sys.C+\(\"w_5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wem_${S5}_wEff\",-1,0,1,0,\"${GSF}\",1\)

# # #
# # # Z->ee
# # #
# root -l -q computeAccSelZeeBinned_Sys.C+\(\"z_5.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/PRELIM_Zee_${S5}_full\",\"${S5}\",0,1,0,\"${GSF}\",0\)


# rm *.so *.d
