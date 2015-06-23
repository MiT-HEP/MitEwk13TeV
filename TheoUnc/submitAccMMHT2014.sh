#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=accWrapper.sh
INPUT_DIR=/afs/cern.ch/work/j/jlawhorn/theo-unc-06-10
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/MMHT2014_amc
ENV_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/MMHT2014_amc
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
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}    

    elif [[ ${CHAN} = "wpe" ]]
    then 
	CHANNUM=1
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}    

    elif [[ ${CHAN} = "wmm" ]]
    then 
	CHANNUM=2
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}    

    elif [[ ${CHAN} = "wpm" ]]
    then 
	CHANNUM=3
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}    
    fi
done

for CHAN in zee zmm
do

    if [[ ${CHAN} = "zee" ]]
    then 
	CHANNUM=0
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}
    
    elif [[ ${CHAN} = "zmm" ]]
    then 
	CHANNUM=1
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
	bsub -q 8nh  -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}    
    fi
done
