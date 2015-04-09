#!/bin/bash

SCRAM_DIR=$1
IN_FILE=$2
ENV_DIR=$3
OUT_FILE=$4
RUN_MACRO=$5
SO_FILE=$6
PDF=$7
F_REP=$8
L_REP=$9

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
    echo root -l -q ${RUN_MACRO}+\(\"${IN_FILE}\",\"${ENV_DIR}/\",\"${PDF}\",${i}\) >> ${OUT_FILE}
    root -l -q ${RUN_MACRO}+\(\"${IN_FILE}\",\"${ENV_DIR}/\",\"${PDF}\",${i}\) >> ${OUT_FILE}
done

status=`echo $?`
echo "Status - $status"

exit $status
