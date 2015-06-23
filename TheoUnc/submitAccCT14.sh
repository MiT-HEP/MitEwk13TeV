#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=accWrapper.sh
INPUT_DIR=/afs/cern.ch/work/j/jlawhorn/theo-unc-06-10
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/CT14_amc
ENV_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/CT14_amc
LOG_DIR=/afs/cern.ch/work/j/jlawhorn/public/log-files

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp acceptGenW.C $WORK_DIR
cp acceptGenW_C.so $WORK_DIR
cp acceptGenZ.C $WORK_DIR
cp acceptGenZ_C.so $WORK_DIR

for CHAN in wpm wme wpe wmm
do
    if [[ ${CHAN} = "wme" ]] 
    then 
	CHANNUM=0
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct14nlo.out -e ${LOG_DIR}/${CHAN}_ct14nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo.txt acceptGenW.C acceptGenW_C.so CT14nlo 0 56 ${CHANNUM}

	for AS in 116 117 119 120
	do 
	    bsub -q 1nh -o ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so CT14nlo_as_0${AS} 0 0 ${CHANNUM}
	done

    elif [[ ${CHAN} = "wpe" ]]
    then 
	CHANNUM=1
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct14nlo.out -e ${LOG_DIR}/${CHAN}_ct14nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo.txt acceptGenW.C acceptGenW_C.so CT14nlo 0 56 ${CHANNUM}
	
	for AS in 116 117 119 120
	do 
	    bsub -q 1nh -o ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so CT14nlo_as_0${AS} 0 0 ${CHANNUM}
	done
	
    elif [[ ${CHAN} = "wmm" ]]
    then
	CHANNUM=2
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct14nlo.out -e ${LOG_DIR}/${CHAN}_ct14nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo.txt acceptGenW.C acceptGenW_C.so CT14nlo 0 56 ${CHANNUM}
	
	for AS in 116 117 119 120
	do 
	    bsub -q 1nh -o ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so CT14nlo_as_0${AS} 0 0 ${CHANNUM}
	done
	
    elif [[ ${CHAN} = "wpm" ]]
    then 
	CHANNUM=3
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct14nlo.out -e ${LOG_DIR}/${CHAN}_ct14nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo.txt acceptGenW.C acceptGenW_C.so CT14nlo 0 56 ${CHANNUM}
	
	for AS in 116 117 119 120
	do 
	bsub -q 1nh -o ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so CT14nlo_as_0${AS} 0 0 ${CHANNUM}
	done

    fi

done

for CHAN in zmm zee
do

    if [[ ${CHAN} = "zee" ]]
    then 
	CHANNUM=0
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct14nlo.out -e ${LOG_DIR}/${CHAN}_ct14nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo.txt acceptGenZ.C acceptGenZ_C.so CT14nlo 0 56 ${CHANNUM}

	for AS in 116 117 119 120
	do 
	    bsub -q 1nh -o ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so CT14nlo_as_0${AS} 0 0 ${CHANNUM}
	done
	
    elif [[ ${CHAN} = "zmm" ]]
    then 
	CHANNUM=1
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct14nlo.out -e ${LOG_DIR}/${CHAN}_ct14nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo.txt acceptGenZ.C acceptGenZ_C.so CT14nlo 0 56 ${CHANNUM}
	
	for AS in 116 117 119 120
	do 
	    bsub -q 1nh -o ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct14nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct14nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so CT14nlo_as_0${AS} 0 0 ${CHANNUM}
	done
	
    fi
    
done