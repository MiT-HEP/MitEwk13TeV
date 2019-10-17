#! /bin/bash


LUMI13=199.270 # Low PU 13 TeV 2017

# # # LUMI2=199.270 # not used in Low PU atm
# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/Wmunu/ntuples/
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/Wmunu/ntuples/
# # S="5TeV"
S="13TeV"

EFFDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics
FSIT=${EFFDIR}/SysUnc_MuSITEff.root
FSTA=${EFFDIR}/SysUnc_MuStaEff.root
FGSF=${EFFDIR}/SysUnc_EleGSFSelEff.root
# FSTA2=${EFFDIR}/_v2_MuSta_Direct/SysUnc_MuStaEff.root

# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"data_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm0_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm1_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm2_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wx_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zxx_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zz_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wz_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"ww_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"top_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/AntiWmunu/ntuples
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWmunu/ntuples

INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/AntiWmunu/ntuples/
OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWmunu/ntuples

# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"data_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm0_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm1_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm2_select.root\",\"${FSIT}\",\"${FSTA}\"\);
root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wx_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zxx_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zz_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wz_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"ww_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"top_select.root\",\"${FSIT}\",\"${FSTA}\"\);

#rm *.so *.d
