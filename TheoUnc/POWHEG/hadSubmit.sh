#!/bin/bash
#-----------------------------------------------------------------
# submission script for hadronizing events in LHE format to 
# CMSSW GEN format
#-----------------------------------------------------------------

INPUT_FOLDER=/afs/cern.ch/work/j/jlawhorn/POWHEG-LHE/CT10nlo_rw_13
OUTPUT_FOLDER=/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN/CT10nlo_rw_13

cp hadronizer.py ${CMSSW_BASE}

for CHAN in wpe wmm wpm zee zmm wme
do
    bsub -q 8nh hadWrapper.sh ${CMSSW_BASE} file:${INPUT_FOLDER}/${CHAN}-events.lhe file:${OUTPUT_FOLDER}/${CHAN}.root
done