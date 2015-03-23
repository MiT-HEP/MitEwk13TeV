#!/bin/bash

#WORK_DIR=$CMSSW_BASE
WORK_DIR=.
SCRIPT=accWrapper.sh
#INPUT_DIR=/afs/cern.ch/work/j/jlawhorn
#INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/Ewk8TeV
INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/PYTHIA-CT10nlo
#OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-8tev-acc
OUTPUT_DIR=.

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp acc*.C $WORK_DIR
cp acc*_C.so $WORK_DIR

for CHAN in wme #wpe wmm wpm zee zmm
do
    #echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/phil_wme_select_2.root ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo 0 52
    #./${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/phil_wme_select_2.root ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}_test.C acc${CHAN}_test_C.so CT10nlo 0 52

    #for AS in `seq 115 123`
    #do 
#	echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/phil_wme_select_2.root ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo_as_0${AS} 0 0
#	./${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/phil_wme_select_2.root ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}_test.C acc${CHAN}_test_C.so CT10nlo_as_0${AS} 0 0
    #done

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo 0 52
    ./${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo 0 52

    for AS in `seq 115 123`
    do 
	echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo_as_0${AS} 0 0
	./${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo_as_0${AS} 0 0
    done
done