#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=accWrapper.sh
INPUT_DIR=root://eoscms.cern.ch//store/user/jlawhorn/PYTHIA-CT10nlo-RW-13
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-acc
ENV_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes
LOG_DIR=/afs/cern.ch/work/j/jlawhorn/public/log-files

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp acc*deta.C $WORK_DIR
cp acc*deta_C.so $WORK_DIR

for CHAN in wme wpe wmm wpm zee zmm
do
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so CT10nlo 0 52
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct10nlo.out -e ${LOG_DIR}/${CHAN}_ct10nlo.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so CT10nlo 0 52

    for AS in `seq 115 123`
    do 
	echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so CT10nlo_as_0${AS} 0 0
	bsub -q 8nh -o ${LOG_DIR}/${CHAN}_ct10nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_ct10nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so CT10nlo_as_0${AS} 0 0
    done

    AS=119
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 100
    bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 100
    AS=118
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 72
    bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 72
    AS=120
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 72
    bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 72
    AS=117
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 27
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 27
    AS=121
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 27
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 27
    AS=116
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 5
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 5
    AS=122
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 5
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF23_nlo_as_0${AS} 0 5

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl 0 40
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mstw2008nlo68cl.out -e ${LOG_DIR}/${CHAN}_mstw2008nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz+68cl 0 40
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.out -e ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz+68cl 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz-68cl 0 40
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.out -e ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz-68cl 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz+68clhalf 0 40
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.out -e ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz+68clhalf 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz-68clhalf 0 40
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.out -e ${LOG_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MSTW2008nlo68cl_asmz-68clhalf 0 40

    AS=119
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF30_nlo_as_0${AS} 0 100
    bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF30_nlo_as_0${AS} 0 100

    AS=117
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF30_nlo_as_0${AS} 0 72
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF30_nlo_as_0${AS} 0 72

    AS=121
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF30_nlo_as_0${AS} 0 72
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so NNPDF30_nlo_as_0${AS} 0 72

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MMHT2014nlo68cl 0 50
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo68cl.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MMHT2014nlo68cl 0 50

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MMHT2014nlo_asmzsmallrange 0 50
    bsub -q 8nh -o ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.out -e ${LOG_DIR}/${CHAN}_mmht2014nlo68cl_asmzsmallrange.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}_bacon.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_mmht2014nlo_asmzsmallrange.txt acc${CHAN}_deta.C acc${CHAN}_deta_C.so MMHT2014nlo_asmzsmallrange 0 50
    
done