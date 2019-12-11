#! /bin/bash

# # # LUMI2=199.270 # not used in Low PU atm
 # INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_v4/Wenu
 # OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/Wenu
INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017/Wenu
OUTPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017_PP_full/Wenu
# S="5TeV"
S="13TeV"

EFFDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics
FSIT=${EFFDIR}/SysUnc_MuSITEff.root
FSTA=${EFFDIR}/SysUnc_MuStaEff.root
FGSF=${EFFDIR}/SysUnc_EleGSFSelEff.root

NSEC=5
ITH=3

NTUP=ntuples_${ITH}_${NSEC}
# NTUP=ntuples

# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"data_select.root\",\"${FGSF}\"\);
 root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we0_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we1_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we2_select.root\",\"${FGSF}\"\);
 # root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx0_select.root\",\"${FGSF}\"\);
 # root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx1_select.root\",\"${FGSF}\"\);
 # root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx2_select.root\",\"${FGSF}\"\);
 # root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zxx_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zz_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wz_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"ww_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top1_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top2_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top3_select.root\",\"${FGSF}\"\);

# INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWenu
# INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5/AntiWenu_sieie
# OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP/AntiWenu_sieie
INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017/AntiWenu
OUTPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017_PP_full/AntiWenu

# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"data_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we0_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we1_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we2_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx0_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx1_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx2_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zxx_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zz_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wz_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"ww_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top1_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top2_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top3_select.root\",\"${FGSF}\"\);

#rm *.so *.d
