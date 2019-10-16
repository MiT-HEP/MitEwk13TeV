#! /bin/bash


# integrated luminosity for data
# LUMI=213.1 # Low PU 13 TeV 2017
LUMI13=199.270 # Low PU 13 TeV 2017
LUMI5=294.4 # Low PU 5 TeV 2017 ele trigger
# LUMI5=291.1 # Low PU 5 TeV 2017 muon trigger

# # # LUMI2=199.270 # not used in Low PU atm
# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/Wenu/ntuples/
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/Wenu/ntuples/
# S="5TeV"
S="13TeV"

EFFDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics
FSIT=${EFFDIR}/SysUnc_MuSITEff.root
FSTA=${EFFDIR}/SysUnc_MuStaEff.root
FGSF=${EFFDIR}/SysUnc_EleGSFSelEff.root

# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"data_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"we0_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"we1_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"we2_select.root\",\"${FGSF}\"\);
 # root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wx_select.root\",\"${FGSF}\"\);
 # root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zxx_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zz_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wz_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"ww_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"top_select.root\",\"${FGSF}\"\);

INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/AntiWenu/ntuples
OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWenu/ntuples


# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"data_select.root\",\"${FGSF}\"\);
root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"we_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wx_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zxx_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"zz_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"wz_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"ww_select.root\",\"${FGSF}\"\);
# root -l -q eleNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"${S}\",\"top_select.root\",\"${FGSF}\"\);

#rm *.so *.d
