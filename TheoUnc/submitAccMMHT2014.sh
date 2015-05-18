#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=accWrapper.sh
INPUT_DIR=/afs/cern.ch/work/j/jlawhorn/flat-05-14
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes-new
ENV_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes-new
LOG_DIR=/afs/cern.ch/work/j/jlawhorn/public/log-files

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp acceptGenW.C $WORK_DIR
cp acceptGenW_C.so $WORK_DIR
cp acceptGenZ.C $WORK_DIR
cp acceptGenZ_C.so $WORK_DIR

for CHAN in wpmNLOEWK #wmeNLOEWK wpeNLOEWK wmmNLOEWK wpmNLOEWK zeeNLOEWK zmmNLOEWK
do
    if [[ ${CHAN} = "wmeNLOEWK" ]]; then CHANNUM=0; 
    elif [[ ${CHAN} = "wpeNLOEWK" ]]; then CHANNUM=1;
    elif [[ ${CHAN} = "wmmNLOEWK" ]]; then CHANNUM=2;
    elif [[ ${CHAN} = "wpmNLOEWK" ]]; then CHANNUM=3; fi;

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
    #bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo_asmzsmallrange 0 4 ${CHANNUM}
    #bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenW.C acceptGenW_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}
    
done

#for CHAN in zeeNLOEWK zmmNLOEWK
#do

#    if [[ ${CHAN} = "zeeNLOEWK" ]]; then CHANNUM=0; 
#    elif [[ ${CHAN} = "zmmNLOEWK" ]]; then CHANNUM=1; fi;

#    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenW.C acceptGenZ_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}
#    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo68cl 0 50 ${CHANNUM}

#    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}
#    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acceptGenZ.C acceptGenZ_C.so MMHT2014nlo_asmzsmallrange 0 50 ${CHANNUM}
    
#done