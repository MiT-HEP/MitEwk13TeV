#!/bin/bash

exe=./runStep2Step3.sh
PU=`echo "${2}"`
ID=`echo "${3}"`
# NBINS=96 #96 bins for muon channels
NBINS=1 #96 bins for STANDALONE MUON channels
# NBINS=96 #96 bins for electron channels
NTOYS=100 # should be default 1000
# EFFTYPE=MuSITEff
EFFTYPE=MuStaEff
# EFFTYPE=EleGSFSelEff

# declare -a EFFTYPES=("MuHLTEff" "MuSelEff" "MuStaEff") #for muons
declare -a CHARGES=("Combined") #for muons
# declare -a CHARGES=("Positive" "Negative") #for muons
FOLDER=Zmm # or Zee
# FOLDER=Zee # or Zee
# CHARGE=Negative #or Positive
# CHARGE=Positive
# CHARGE=Combined # combine pos & neg for the muon standalone nightmare
 
VERS=_TEST_CHECKCLEANING
# #CLOSURE 
# POSTFIX=_aMCxPythia${VERS}
# POSTFIX_alt=_aMCxPythia${VERS}

# ## FSR alternative shape evaluation
# POSTFIX=_POWxPythia${VERS}
# POSTFIX_alt=_POWxPhotos${VERS}

# reversed? Test
# POSTFIX_alt=_POWxPythia${VERS}
# POSTFIX=_POWxPhotos${VERS}


# MC alternative shape eval
# POSTFIX=_aMCxPythia${VERS}
# POSTFIX_alt=_minloxPythia${VERS}

## BKG alternative shape eval
POSTFIX_alt=_aMCxPythia${VERS}
POSTFIX=_POWBKG${VERS}

# POSTFIX=_aMCxPythia_staFit7
# POSTFIX_alt=_POWBKG_staFit7

# declare -a BINLIST=()
for ((bin=0; bin<${NBINS};bin++))
# for ((bin=16; bin<${NBINS};bin++))
do
    for CHARGE in "${CHARGES[@]}"
    do
        echo "submitting jobs to calculate efficiencies for ${EFFTYPE} ${CHARGE}, total bins ${NBINS} with ${NTOYS} toys" 
        echo "executable              = "${exe} > tmp.sub
        echo "arguments               = \$(ClusterId) \$(ProcId)" ${bin} ${NTOYS} ${FOLDER} ${EFFTYPE} ${CHARGE} ${POSTFIX} ${POSTFIX_alt}>> tmp.sub
        echo "output                  = output/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).out" >> tmp.sub
        echo "error                   = error/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).err"  >> tmp.sub
        echo "log                     = log/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).log"    >> tmp.sub
        echo "+JobFlavour = \"tomorrow\"  " >> tmp.sub
        # echo "+JobFlavour = \"longlunch\"  " >> tmp.sub
        # echo "+JobFlavour = \"workday\"  " >> tmp.sub
        # echo "+JobFlavour = \"espresso\"  " >> tmp.sub
        echo "queue" >> tmp.sub
        # echo "queue ${NBINS}" >> tmp.sub
        condor_submit tmp.sub
    done
done
