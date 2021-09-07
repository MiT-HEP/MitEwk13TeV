#! /bin/bash

# integrated luminosity for data
LUMI13=199.270 # Low PU 13 TeV 2017
LUMI5=294.4 # Low PU 5 TeV 2017 ele trigger

INPUTDIR=/tmp/arapyan/Wmunu/ntuples
OUTPUTDIR=/afs/cern.ch/work/a/arapyan/FilesSM2017/ntuples_new_mu/
# S="5TeV"
S="13TeV"

EFFDIR=/afs/cern.ch/work/a/arapyan/Run2/test/CMSSW_10_2_13/src/MitEwk13TeV/data
FSIT=${EFFDIR}/SysUnc_MuSITEff_FINAL.root
FSTA=${EFFDIR}/SysUnc_MuStaEff_FINAL.root
# FSTA2=${EFFDIR}/_v2_MuSta_Direct/SysUnc_MuStaEff.root

# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"data_select.root\",\"${FSIT}\",\"${FSTA}\"\);
root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm0_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm1_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm2_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wx_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zxx_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zz_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wz_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"ww_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"top_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/AntiWmunu/ntuples
# OUTPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_13TeV_wRecoil_wStat_2G_eff/AntiWmunu/ntuples

# INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV/AntiWmunu/ntuples
# # INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/Wmunu/ntuples/wparts
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_wRecoil/AntiWmunu/ntuples

# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"data_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm0_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm1_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wm2_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wx_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zxx_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zz_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wz_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"ww_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"top_select.root\",\"${FSIT}\",\"${FSTA}\"\);

#rm *.so *.d
