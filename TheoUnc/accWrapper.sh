#!/bin/bash

#void acceptGenW(TString input="/afs/cern.ch/work/j/jlawhorn/flat-05-14/WJets.root",
#                TString outputDir="./",
#                TString pdfName="NNPDF30_nlo_as_0118",
#                Int_t setMin=0,
#                Int_t setMax=4,
#                Int_t proc=0) {

SCRAM_DIR=$1
  IN_FILE=$2
  ENV_DIR=$3
 OUT_FILE=$4
RUN_MACRO=$5
  SO_FILE=$6
      PDF=$7
    F_REP=$8
    L_REP=$9
     CHAN=${10}

WORK_DIR=`pwd`
echo `hostname`
echo "args:  $*"

cd ${SCRAM_DIR}/src
eval `scramv1 runtime -sh`
cd ${WORK_DIR}

cp ${SCRAM_DIR}/rootlogon.C .
cp ${SCRAM_DIR}/lhapdf_init.sh .
cp ${SCRAM_DIR}/${RUN_MACRO} .
cp ${SCRAM_DIR}/${SO_FILE} .

source lhapdf_init.sh

echo ${IN_FILE} ${PDF} > ${OUT_FILE}

echo root -l -q ${RUN_MACRO}+\(\"${IN_FILE}\",\"${ENV_DIR}/\",\"${PDF}\",${F_REP},${L_REP},${CHAN}\) >> ${OUT_FILE}
root -l -q ${RUN_MACRO}+\(\"${IN_FILE}\",\"${ENV_DIR}/\",\"${PDF}\",${F_REP},${L_REP},${CHAN}\) >> ${OUT_FILE}

status=`echo $?`
echo "Status - $status"

exit $status
