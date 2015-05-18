#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=flatWrapper.sh
#INPUT_DIR=/store/user/jlawhorn/NNPDF30_PWHG2_RW
INPUT_DIR=/store/user/arapyan/Run2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_Wjets_genbacon_7_4_t4/150513_161203/0000
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/CT10nlo-PYTHIA6-13-v2
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/aMCNLO/DYJetsToLL_500000
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/aMCNLO/WJetsToLNu_50001
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/flat-05-14

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp makeFlat.C $WORK_DIR
cp makeFlat_C* $WORK_DIR

#SCRAM_DIR=$1
#IN_FILE=$2
#OUT_FILE=$3
#RUN_MACRO=$4
#SO_FILE=$5

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${INPUT_DIR} | grep bacon`
for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${INPUT_DIR}`
do
    if [[ -e ${OUTPUT_DIR}/${file%-*}-flat.root ]]
    then
        continue
    else
	echo bsub -q 1nh ${SCRIPT} ${WORK_DIR} root://eoscms/${INPUT_DIR}/${file} ${OUTPUT_DIR}/${file%-*}-flat.root makeFlat.C makeFlat_C.so
	bsub -q 1nh ${SCRIPT} ${WORK_DIR} root://eoscms/${INPUT_DIR}/${file} ${OUTPUT_DIR}/${file%-*}-flat.root makeFlat.C makeFlat_C.so
   fi
done