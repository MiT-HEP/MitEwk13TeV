#!/bin/bash

exe=./runStep2Step3.sh
PU=`echo "${2}"`
ID=`echo "${3}"`
# NBINS=96 #96 bins for muon channels
# NBINS=24 #24 bins for STANDALONE MUON channels
NBINS=36 #36 for Zee and Zmm
NTOYS=1000 # should be default 1000
# EFFTYPE=MuSITEff
# EFFTYPE=MuStaEff
EFFTYPE=EleGSFSelEff

# declare -a EFFTYPES=("MuHLTEff" "MuSelEff" "MuStaEff") #for muons
declare -a CHARGES=("Combined") #for muons
# declare -a CHARGES=("Positive" "Negative") #for muons
# declare -a CHARGES=("Negative") #for muons
# FOLDER=Zmm # or Zee
FOLDER=Zee # or Zee
# CHARGE=Negative #or Positive
# CHARGE=Positive
# CHARGE=Combined # combine pos & neg for the muon standalone nightmare
 
VERS=

# ## FSR alternative shape evaluation
POSTFIX=_POWxPythia${VERS}
POSTFIX_alt=_POWxPhotos${VERS}


# MC alternative shape eval
# POSTFIX=_aMCxPythia${VERS}
# POSTFIX_alt=_minloxPythia${VERS}

## BKG alternative shape eval
# POSTFIX=_aMCxPythia${VERS}
# POSTFIX_alt=_POWBKG${VERS}

# POSTFIX=_aMCxPythia_staFit7
# POSTFIX_alt=_POWBKG_staFit7
# bin = 0
# declare -a BINLIST=()
# for ((bin=12; bin<${NBINS};bin++))
# # for ((bin=16; bin<${NBINS};bin++))
# do
    for CHARGE in "${CHARGES[@]}"
    do
        echo "submitting jobs to calculate efficiencies for ${EFFTYPE} ${CHARGE}, total bins ${NBINS} with ${NTOYS} toys" 
        echo "executable              = "${exe} > tmp.sub
        echo "arguments               = \$(ClusterId) \$(ProcId)" ${NTOYS} ${FOLDER} ${EFFTYPE} ${CHARGE} ${POSTFIX} ${POSTFIX_alt}>> tmp.sub
        echo "output                  = output/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).out" >> tmp.sub
        echo "error                   = error/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).err"  >> tmp.sub
        echo "log                     = log/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).log"    >> tmp.sub
        echo "+JobFlavour = \"tomorrow\"  " >> tmp.sub
        # echo "+JobFlavour = \"longlunch\"  " >> tmp.sub
        # echo "+JobFlavour = \"workday\"  " >> tmp.sub
        # echo "+JobFlavour = \"espresso\"  " >> tmp.sub
        # echo "queue" >> tmp.sub
        echo "queue ${NBINS}" >> tmp.sub
        condor_submit tmp.sub
    done
# done
