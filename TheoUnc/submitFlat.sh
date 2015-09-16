#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=flatWrapper.sh
#INPUT_DIR=/store/user/jlawhorn/amcnlo_dy
INPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/HORACE_NOT_FUCKED/z_*_born_final/
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

#for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${INPUT_DIR} |grep p6.root`
for file in `ls ${INPUT_DIR}/*photos.root`
do
    if [[ -e ${OUTPUT_DIR}/${file%.*}-flat.root ]]
    then
        continue
    else
	#echo bsub -q 1nh ${SCRIPT} ${WORK_DIR} root://eoscms/${INPUT_DIR}/${file} ${OUTPUT_DIR}/${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
	#bsub -q 1nh ${SCRIPT} ${WORK_DIR} root://eoscms/${INPUT_DIR}/${file} ${OUTPUT_DIR}/${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
	echo bsub -q 1nh ${SCRIPT} ${WORK_DIR} ${file} ${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so
	bsub -q 1nh ${SCRIPT} ${WORK_DIR} ${file} ${file%.*}-flat.root 23 makeFlat.C makeFlat_C.so

   fi
done