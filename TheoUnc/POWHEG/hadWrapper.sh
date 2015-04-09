#!/bin/bash
#-----------------------------------------------------------------
# wrapper script for hadronizing events in LHE format to
# CMSSW GEN format
#-----------------------------------------------------------------
SCRAM_DIR=$1
IN_FILE=$2
OUT_FILE=$3

WORK_DIR=`pwd`
echo `hostname`
echo "args:  $*"

cd ${SCRAM_DIR}/src
eval `scramv1 runtime -sh`
cd ${WORK_DIR}

cp ${SCRAM_DIR}/hadronizer.py .

cmsRun hadronizer.py ${IN_FILE} ${OUT_FILE}

status=`echo $?`
echo "Status - $status"

exit $status
