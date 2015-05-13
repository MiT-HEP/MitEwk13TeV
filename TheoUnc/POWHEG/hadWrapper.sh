#!/bin/bash

SCRAM_DIR=$1
HAD_SCRIPT=$2
IN_FILE=$3
EOS_OUT=$4
OUT_FILE=$5
MAX_EVTS=$6
SKIP_EVTS=$7

WORK_DIR=`pwd`
echo `hostname`
echo "args:  $*"

cd ${SCRAM_DIR}/src
eval `scramv1 runtime -sh`
cd ${WORK_DIR}

cp ${SCRAM_DIR}/${HAD_SCRIPT} .

cmsRun ${HAD_SCRIPT} ${IN_FILE} ${OUT_FILE} ${MAX_EVTS} ${SKIP_EVTS}

xrdcp ${OUT_FILE} ${EOS_OUT}

status=`echo $?`
echo "Status - $status"

exit $status
