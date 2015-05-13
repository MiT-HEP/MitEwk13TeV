#!/bin/bash

SCRAM_DIR=$1
PWHG_DIR=$2
RUN_PRGM=$3
CONF_FILE=$4

#----------------------------

WORK_DIR=`pwd`
echo `hostname`
echo "args:  $*"

cd ${SCRAM_DIR}/src
eval `scramv1 runtime -sh`
cd ${WORK_DIR}

cp ${PWHG_DIR}/lhapdf_init.sh .
cp ${PWHG_DIR}/${CONF_FILE} .

source lhapdf_init.sh

#----------------------------

mkdir ${CONF_FILE%-*}
mv ${CONF_FILE} ${CONF_FILE%-*}/${CONF_FILE}
cd ${CONF_FILE%-*}
echo ${CONF_FILE%-*} | ${RUN_PRGM}

cd ..
xrdcp ${CONF_FILE%-*} root://eoscms.cern.ch//store/user/jlawhorn/LHE/${CONF_FILE%-*}

status=`echo $?`
echo "Status - $status"

exit $status
