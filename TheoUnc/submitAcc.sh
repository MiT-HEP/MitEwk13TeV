#!/bin/bash

WORK_DIR=$CMSSW_BASE
SCRIPT=accWrapper.sh
INPUT_DIR=/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN/CT10nlo_rw_13
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-acc

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp acc*.C $WORK_DIR
cp acc*_C.so $WORK_DIR

for CHAN in wme wpe wmm wpm zee zmm
do
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo 0 52
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo 0 52

    for AS in `seq 115 123`
    do 
	echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo_as_0${AS} 0 0
	bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_ct10nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so CT10nlo_as_0${AS} 0 0
    done

    AS=119
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 100
    bsub -q 1nd ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 100
    AS=118
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 72
    bsub -q 1nd ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 72
    AS=120
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 72
    bsub -q 1nd ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 72
    AS=117
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 27
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 27
    AS=121
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 27
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 27
    AS=116
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 5
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 5
    AS=122
    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 5
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_nnpdf23_nlo_as_0${AS}.txt acc${CHAN}.C acc${CHAN}_C.so NNPDF23_nlo_as_0${AS} 0 5

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl 0 40
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz+68cl 0 40
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68cl.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz+68cl 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz-68cl 0 40
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68cl.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz-68cl 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz+68clhalf 0 40
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz+68clhalf.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz+68clhalf 0 40

    echo ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz-68clhalf 0 40
    bsub -q 8nh ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/${CHAN}-bacon.root ${OUTPUT_DIR}/${CHAN}_mstw2008nlo68cl_asmz-68clhalf.txt acc${CHAN}.C acc${CHAN}_C.so MSTW2008nlo68cl_asmz-68clhalf 0 40
    
done