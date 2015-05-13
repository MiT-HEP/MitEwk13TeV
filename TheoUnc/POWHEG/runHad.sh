#!/bin/bash

#SCRAM_DIR=$1
#HAD_SCRIPT=$2
#IN_FILE=$3
#EOS_OUT=$4
#OUT_FILE=$5
#MAX_EVTS=$6
#SKIP_EVTS=$7

HAD_SCRIPT=pythia_with_powheg_veto.py

cp ${HAD_SCRIPT} ${CMSSW_BASE}/${HAD_SCRIPT}

for CHAN in zee zmm wpm wmm wme wpe
do
    for i in 0 100000 200000 300000 400000 500000 600000 700000 800000 900000
    do 
	bsub -q 2nd hadWrapper.sh ${CMSSW_BASE} ${HAD_SCRIPT} root://eoscms//eos/cms/store/user/jlawhorn/LHE/NNPDF30_PWHG2_RW/${CHAN}NLOEWK-events.lhe root://eoscms.cern.ch//store/user/jlawhorn/NNPDF30_PWHG2_RW/ ${CHAN}NLOEWK-pythia8-powheg2-${i}.root 100000 ${i}
    done 
done

