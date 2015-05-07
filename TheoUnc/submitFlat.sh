#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=flatWrapper.sh
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/NNPDF30-GENSIM
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/CT10nlo-PYTHIA6-13-v2
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/aMCNLO/DYJetsToLL_500000
INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/aMCNLO/WJetsToLNu_50001
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/baseComp-05-06

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp makeFlat* $WORK_DIR

#SCRAM_DIR=$1
#IN_FILE=$2
#OUT_FILE=$3
#RUN_MACRO=$4
#SO_FILE=$5



#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/user/jlawhorn/NNPDF30-GENSIM/ |grep bacon`
#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/user/jlawhorn/CT10nlo-PYTHIA6-13-v2/ |grep bacon`
#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/user/jlawhorn/aMCNLO/DYJetsToLL_500000/`
for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/user/jlawhorn/aMCNLO/WJetsToLNu_50001/`
do
    if [[ -e ${OUTPUT_DIR}/DYJets-${file} ]]
    then
        continue
    else
	echo bsub -q 1nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${file} ${OUTPUT_DIR}/WJets-${file} makeFlat.C makeFlat_C.so
	bsub -q 1nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${file} ${OUTPUT_DIR}/WJets-${file} makeFlat.C makeFlat_C.so
    fi
done