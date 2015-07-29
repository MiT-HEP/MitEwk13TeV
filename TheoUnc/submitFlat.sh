#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=flatWrapper.sh
INPUT_DIR=/store/user/jlawhorn/amcnlo_w
OUTPUT_DIR=`pwd`

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp makeFlat.C $WORK_DIR
cp makeFlat_C* $WORK_DIR

#SCRAM_DIR=$1
#IN_FILE=$2
#OUT_FILE=$3
#BOS_ID=$4
#RUN_MACRO=$5
#SO_FILE=$6

for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${INPUT_DIR}`
#for file in `ls ${INPUT_DIR}/*_born_v_right/*root`
do
    if [[ -e ${file%.*}-flat.root ]]
    then
        continue
    else
	echo ${file%.*}-flat.root
	echo bsub -q 8nh ${SCRIPT} ${WORK_DIR} root://eoscms/${INPUT_DIR}/${file} ${OUTPUT_DIR}/${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
	bsub -q 8nh ${SCRIPT} ${WORK_DIR} root://eoscms/${INPUT_DIR}/${file} ${OUTPUT_DIR}/${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
	#echo bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${file} ${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
	#bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${file} ${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
   fi
done