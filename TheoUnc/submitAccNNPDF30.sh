#!/bin/bash

#void acceptGenW(TString input="/afs/cern.ch/work/j/jlawhorn/flat-05-14/WJets.root",
#                TString outputDir="./",
#                TString pdfName="NNPDF30_nlo_as_0118",
#                Int_t setMin=0,
#                Int_t setMax=4,
#                Int_t proc=0) {

#SCRAM_DIR=$1
#  IN_FILE=$2
#  ENV_DIR=$3
# OUT_FILE=$4
#RUN_MACRO=$5
#  SO_FILE=$6
#      PDF=$7
#    F_REP=$8
#    L_REP=$9
#     CHAN=${10}

#enum PROC {wme=0, wpe, wmm, wpm};

WORK_DIR=$CMSSW_BASE
SCRIPT=accWrapper.sh
INPUT_DIR=/afs/cern.ch/work/j/jlawhorn/theo-unc-06-10
OUTPUT_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc
ENV_DIR=/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc
LOG_DIR=/afs/cern.ch/work/j/jlawhorn/public/log-files

cp rootlogon.C $WORK_DIR
cp lhapdf_init.sh $WORK_DIR
cp acceptGenW.C $WORK_DIR
cp acceptGenW_C.so $WORK_DIR
cp acceptGenZ.C $WORK_DIR
cp acceptGenZ_C.so $WORK_DIR

############
# CENTRAL 118
############

#CHAN=wme
#CHANNUM=0
#AS=118
#bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 100 ${CHANNUM}
#
#CHAN=wpe
#CHANNUM=1
#bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 100 ${CHANNUM}
#
#CHAN=wmm
#CHANNUM=2
#bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 100 ${CHANNUM}
#
#CHAN=wpm
#CHANNUM=3
#bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 100 ${CHANNUM}
#
#CHAN=zee
#CHANNUM=0
#bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 100 ${CHANNUM}
#
#CHAN=zmm
#CHANNUM=1
#bsub -q 1nd -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 100 ${CHANNUM}


############
# +1 119
############

#CHAN=wme
#CHANNUM=0
AS=119
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=wpe
#CHANNUM=1
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=wmm
#CHANNUM=2
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=wpm
#CHANNUM=3
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=zee
#CHANNUM=0
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}

CHAN=zmm
CHANNUM=1
bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}


############
# -1 117
############

#CHAN=wme
#CHANNUM=0
#AS=117
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=wpe
#CHANNUM=1
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=wmm
#CHANNUM=2
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=wpm
#CHANNUM=3
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=zee
#CHANNUM=0
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}
#
#CHAN=zmm
#CHANNUM=1
#bsub -q 8nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 72 ${CHANNUM}

############
# +3 121
############

#CHAN=wme
#CHANNUM=0
AS=121
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=wpe
#CHANNUM=1
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=wmm
#CHANNUM=2
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=wpm
#CHANNUM=3
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=zee
#CHANNUM=0
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}

CHAN=zmm
CHANNUM=1
bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}

############
# -3 115
############

#CHAN=wme
#CHANNUM=0
#AS=115
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=wpe
#CHANNUM=1
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=wmm
#CHANNUM=2
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WmJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=wpm
#CHANNUM=3
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/WpJetsToLNu.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenW.C acceptGenW_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=zee
#CHANNUM=0
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}
#
#CHAN=zmm
#CHANNUM=1
#bsub -q 1nh -o ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.out -e ${LOG_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.err ${SCRIPT} ${WORK_DIR} ${INPUT_DIR}/DYJetsToLL.root ${ENV_DIR} ${OUTPUT_DIR}/${CHAN}_nnpdf30_nlo_as_0${AS}.txt acceptGenZ.C acceptGenZ_C.so NNPDF30_nlo_as_0${AS} 0 5 ${CHANNUM}