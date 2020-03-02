#! /bin/bash


LUMI13=199.270 # Low PU 13 TeV 2017

# # # # # LUMI2=199.270 # not used in Low PU atm
### INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/Wmunu/ntuples/
# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_v3/Wmunu
 # ### INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v3/Wenu
 
INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v4_fixHist/Wmunu
OUTPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v4_fixHist_PP_full/Wmunu
# # S="5TeV"
S="13TeV"

EFFDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics
FSIT=${EFFDIR}/SysUnc_MuSITEff.root
FSTA=${EFFDIR}/SysUnc_MuStaEff.root
FGSF=${EFFDIR}/SysUnc_EleGSFSelEff.root
# FSTA2=${EFFDIR}/_v2_MuSta_Direct/SysUnc_MuStaEff.root

NSEC=1
ITH=0

NTUP=ntuples_${ITH}_${NSEC}

# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"data_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wm0_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wm1_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wm2_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx0_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx1_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
 # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx2_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
  # root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zxx_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zz_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wz_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"ww_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top1_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top2_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top3_select.raw.root\",\"${FSIT}\",\"${FSTA}\"\);

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_v3/AntiWmunu
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWmunu

INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v4_fixHist/AntiWmunu
OUTPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v4_fixHist_PP_full/AntiWmunu

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV/AntiWmunu/ntuples/
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWmunu/ntuples

NSEC=1
ITH=0

NTUP=ntuples_${ITH}_${NSEC}

root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"data_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wm0_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wm1_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wm2_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx0_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx1_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx2_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zxx_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zz_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wz_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"ww_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top1_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top2_select.root\",\"${FSIT}\",\"${FSTA}\"\);
# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top3_select.root\",\"${FSIT}\",\"${FSTA}\"\);

#rm *.so *.d
