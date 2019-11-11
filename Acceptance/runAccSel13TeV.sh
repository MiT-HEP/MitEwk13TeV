#! /bin/bash


# declare -a SYSDIR=("MuSIT" "MuSta" "MuStaDirect")
# declare -a SYSDIR=("MuSta" "MuStaDirect")
# SYSDIR=("MuStaDirect")
SYSDIR=""

OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance
EFFDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results"
# EFFDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV_v2/results"
EFFSYSDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics/"

S13="13TeV"
S5="5TeV"

# #
# # W->munu
# # #
# # for SYS in SYSDIR; do
  Sta=${EFFSYSDIR}/${SYS}/SysUnc_MuStaEff.root
  SIT=${EFFSYSDIR}/${SYS}/SysUnc_MuSITEff.root
  # root -l -q computeAccSelWm_Sys.C+\(\"w0j_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Wmp_0j_13TeV_wEff${SYS}\",1,0,\"${SIT}\",\"${Sta}\",1\)
  # root -l -q computeAccSelWm_Sys.C+\(\"w1j_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Wmp_1j_13TeV_wEff${SYS}\",1,0,\"${SIT}\",\"${Sta}\",1\)
  # root -l -q computeAccSelWm_Sys.C+\(\"w2j_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Wmp_2j_13TeV_wEff${SYS}\",1,0,\"${SIT}\",\"${Sta}\",1\)
  
  # root -l -q computeAccSelWm_Sys.C+\(\"w0j_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Wmm_0j_13TeV_wEff${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",1\)
  # root -l -q computeAccSelWm_Sys.C+\(\"w1j_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Wmm_1j_13TeV_wEff${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",1\)
  # root -l -q computeAccSelWm_Sys.C+\(\"w2j_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Wmm_2j_13TeV_wEff${SYS}\",-1,0,\"${SIT}\",\"${Sta}\",1\)
# done

#
# Z->mumu

  # root -l -q computeAccSelZmmBinned_Sys.C+\(\"z_13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/Zmm_13TeV_wEff\",0,\"${SIT}\",\"${Sta}\",1\)


# #
# # W->enu
# #
GSF=${EFFSYSDIR}/SysUnc_EleGSFSelEff.root
# root -l -q computeAccSelWe_Sys.C+\(\"w0j_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wep_0j_13TeV_wEff\",1,0,1,0,\"${GSF}\",1\)
root -l -q computeAccSelWe_Sys.C+\(\"w1j_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wep_1j_13TeV_wEff\",1,0,1,0,\"${GSF}\",1\)
root -l -q computeAccSelWe_Sys.C+\(\"w2j_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wep_2j_13TeV_wEff\",1,0,1,0,\"${GSF}\",1\)

# root -l -q computeAccSelWe_Sys.C+\(\"w0j_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wem_0j_13TeV_wEff\",-1,0,1,0,\"${GSF}\",1\)
# root -l -q computeAccSelWe_Sys.C+\(\"w1j_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wem_1j_13TeV_wEff\",-1,0,1,0,\"${GSF}\",1\)
# root -l -q computeAccSelWe_Sys.C+\(\"w2j_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Wem_2j_13TeV_wEff\",-1,0,1,0,\"${GSF}\",1\)

# # #
# # # Z->ee
# # #
# root -l -q computeAccSelZeeBinned_Sys.C+\(\"z_13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/Zee_13TeV_wEff\",\"${S13}\",0,1,0,\"${GSF}\",1\)


# rm *.so *.d
