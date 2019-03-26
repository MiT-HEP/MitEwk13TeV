#!/bin/bash

exe=./runStep2Step3.sh
PU=`echo "${2}"`
ID=`echo "${3}"`
NBINS=64 #64 bins for muon channels
NTOYS=1000 # should be default 1000
# EFFTYPE=MuHLTEff
EFFTYPE=MuSITEff
# EFFTYPE=MuStaEff
# declare -a EFFTYPES=("MuHLTEff" "MuSelEff" "MuStaEff") #for muons
# declare -a CHARGES=("Positive" "Negative") #for muons
FOLDER=Zmm # or Zee
CHARGE=Negative #or Positive
# # CHARGE=Positive
# for EFFTYPE in "${EFFTYPES[@]}"
# do
    # for CHARGE in "${CHARGES[@]}"
    # do
        echo "submitting jobs to calculate efficiencies for ${EFFTYPE} ${CHARGE}, total bins ${NBINS} with ${NTOYS} toys" 
        echo "executable              = "${exe} > tmp.sub
        echo "arguments               = \$(ClusterId) \$(ProcId)" ${NTOYS} ${FOLDER} ${EFFTYPE} ${CHARGE}>> tmp.sub
        echo "output                  = output/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).out" >> tmp.sub
        echo "error                   = error/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).err"  >> tmp.sub
        echo "log                     = log/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).log"    >> tmp.sub
        # echo "+JobFlavour = \"workday\"  " >> tmp.sub
        echo "+JobFlavour = \"longlunch\"  " >> tmp.sub
        # echo "+JobFlavour = \"espresso\"  " >> tmp.sub
        echo "queue ${NBINS}" >> tmp.sub
        condor_submit tmp.sub
    # done
# done
