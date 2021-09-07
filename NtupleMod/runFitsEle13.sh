#! /bin/bash

# # # LUMI2=199.270 # not used in Low PU atm
INPUTDIR=/tmp/arapyan/Wenu
#OUTPUTDIR=/eos/cms/store/group/phys_smp/arapyan/ntuples_ele_update/
OUTPUTDIR=/tmp/arapyan/ntuples_ele_update/
# S="5TeV"
S="13TeV"

EFFDIR=/afs/cern.ch/work/a/arapyan/Run2/test/CMSSW_10_2_13/src/MitEwk13TeV/data/Efficiency/lowpu_13TeV/Systematics/
FSIT=${EFFDIR}/SysUnc_MuSITEff.root
FSTA=${EFFDIR}/SysUnc_MuStaEff.root
FGSF=${EFFDIR}/SysUnc_EleGSFSelEff.root

NSEC=1
ITH=0

NTUP=ntuples_${ITH}_${NSEC}
# NTUP=ntuples

#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"data_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we0_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we1_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we2_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx0_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx1_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wx2_select.root\",\"${FGSF}\"\);
#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zxx_select.root\",\"${FGSF}\"\);
root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"zz_select.root\",\"${FGSF}\"\);
root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"wz_select.root\",\"${FGSF}\"\);
root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"ww_select.root\",\"${FGSF}\"\);
root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top1_select.root\",\"${FGSF}\"\);
root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top2_select.root\",\"${FGSF}\"\);
root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"top3_select.root\",\"${FGSF}\"\);

#INPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_preapproval/AntiWenu
#OUTPUTDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_preapproval_PP/AntiWenu

#root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"data_select.root\",\"${FGSF}\"\);
# root -l -q -b eleNtupleMod.C+\(\"${OUTPUTDIR}/${NTUP}\",\"${INPUTDIR}/${NTUP}\",\"${S}\",\"we_select.root\",\"${FGSF}\"\);
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
