#!/bin/bash

SCRAM_DIR=$1
IN_FILE=$2
OUT_FILE=$3
RUN_MACRO=$4
SO_FILE=$5
PDF=$6
F_REP=$7
L_REP=$8

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

for i in `seq $F_REP $L_REP`;
do
    echo root -l -q ${RUN_MACRO}+\(\"${IN_FILE}\",\"${PDF}\",${i}\) >> ${OUT_FILE}
    root -l -q ${RUN_MACRO}+\(\"${IN_FILE}\",\"${PDF}\",${i}\) >> ${OUT_FILE}
done

status=`echo $?`
echo "Status - $status"

exit $status
