#!/bin/bash

exe=./runStep2Step3.sh
PU=`echo "${2}"`
ID=`echo "${3}"`
NBINS=64 #64 bins for muon channels
# NBINS=96 #96 bins for electron channels
NTOYS=1000 # should be default 1000
# EFFTYPE=MuSITEff
EFFTYPE=MuStaEff
# EFFTYPE=GSFSelEff

# declare -a EFFTYPES=("MuHLTEff" "MuSelEff" "MuStaEff") #for muons
# declare -a CHARGES=("Positive" "Negative") #for muons
FOLDER=Zmm # or Zee
# FOLDER=Zee # or Zee
# CHARGE=Negative #or Positive
# # CHARGE=Positive
CHARGE=Combined # combine pos & neg for the muon standalone nightmare

## FSR alternative shape evaluation
POSTFIX=_POWxPythia_v0
POSTFIX_alt=_POWxPhotos_v0

## MC alternative shape eval
# POSTFIX=_aMCxPythia_v0
# POSTFIX_alt=_minloxPythia_v0

# POSTFIX=_aMCxPythia_staFit7
# POSTFIX_alt=_minloxPythia_staFit7

## BKG alternative shape eval
# POSTFIX=_aMCxPythia_v0
# POSTFIX_alt=_POWBKG_v0

# POSTFIX=_aMCxPythia_staFit7
# POSTFIX_alt=_POWBKG_staFit7

# for EFFTYPE in "${EFFTYPES[@]}"
# do
    # for CHARGE in "${CHARGES[@]}"
    # do
        echo "submitting jobs to calculate efficiencies for ${EFFTYPE} ${CHARGE}, total bins ${NBINS} with ${NTOYS} toys" 
        echo "executable              = "${exe} > tmp.sub
        echo "arguments               = \$(ClusterId) \$(ProcId)" ${NTOYS} ${FOLDER} ${EFFTYPE} ${CHARGE} ${POSTFIX} ${POSTFIX_alt}>> tmp.sub
        echo "output                  = output/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).out" >> tmp.sub
        echo "error                   = error/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).err"  >> tmp.sub
        echo "log                     = log/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).log"    >> tmp.sub
        # echo "+JobFlavour = \"tomorrow\"  " >> tmp.sub
        # echo "+JobFlavour = \"longlunch\"  " >> tmp.sub
        # echo "+JobFlavour = \"microcentury\"  " >> tmp.sub
        # echo "+JobFlavour = \"espresso\"  " >> tmp.sub
        echo "queue ${NBINS}" >> tmp.sub
        condor_submit tmp.sub
    # done
# done
