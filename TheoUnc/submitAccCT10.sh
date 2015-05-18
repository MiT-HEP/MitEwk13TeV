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

#    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acceptGenW.C acceptGenW_C.so CT10nlo 0 52 ${CHANNUM}
#    bsub -q 1nd -o ${LOG_DIR}/${CHAN}_ct10nlo.out -e ${LOG_DIR}/${CHAN}_ct10nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acceptGenW.C acceptGenW_C.so CT10nlo 0 52 ${CHANNUM}

    for AS in 115 116 117 119 120 121
    do 
	echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so CT10nlo_as_0${AS} 0 0 ${CHANNUM}
#	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct10nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct10nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so CT10nlo_as_0${AS} 0 0 ${CHANNUM}
    done

done

#for CHAN in zeeNLOEWK zmmNLOEWK
#do

#    if [[ ${CHAN} = "zeeNLOEWK" ]]; then CHANNUM=0; 
#    elif [[ ${CHAN} = "zmmNLOEWK" ]]; then CHANNUM=1; fi;

#    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acceptGenZ.C acceptGenZ_C.so CT10nlo 0 52 ${CHANNUM}
#    bsub -q 1nd -o ${LOG_DIR}/${CHAN}_ct10nlo.out -e ${LOG_DIR}/${CHAN}_ct10nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acceptGenZ.C acceptGenZ_C.so CT10nlo 0 52 ${CHANNUM}

#    for AS in `seq 115 123`
#    do 
#	echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so CT10nlo_as_0${AS} 0 0 ${CHANNUM}
#	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct10nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct10nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so CT10nlo_as_0${AS} 0 0 ${CHANNUM}
#    done

#done