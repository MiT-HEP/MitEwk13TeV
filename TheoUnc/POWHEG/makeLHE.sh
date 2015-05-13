#!/bin/bash

#SCRAM_DIR=$1
#PWHG_DIR=$2
#RUN_PRGM=$3
#CONF_FILE=$4

for file in `ls z*CT10*powheg.input`
do
    echo bsub -q 2nd powhegWrapper.sh /afs/cern.ch/user/j/jlawhorn/CMSSW_7_2_2_patch1/ `pwd` /afs/cern.ch/work/j/jlawhorn/POWHEG-BOX-V2/W_ew-BMNNPV/pwhg_main ${file}
    bsub -q 2nd powhegWrapper.sh /afs/cern.ch/user/j/jlawhorn/CMSSW_7_2_2_patch1/ `pwd` /afs/cern.ch/work/j/jlawhorn/POWHEG-BOX-V2/Z_ew-BMNNPV/pwhg_main ${file}
done
