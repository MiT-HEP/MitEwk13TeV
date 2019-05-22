#! /bin/bash


# integrated luminosity for data
# LUMI=213.1 # Low PU 13 TeV 2017
LUMI13=199.270 # Low PU 13 TeV 2017
LUMI5=294.4 # Low PU 5 TeV 2017 ele trigger
# LUMI5=291.1 # Low PU 5 TeV 2017 muon trigger

# LUMI2=199.270 # not used in Low PU atm

INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/
OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV_wRecoil/Wmunu/ntuples/

# root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"data_select.root\"\);
root -l -q muonNtupleMod.C+\(\"${OUTPUTDIR}\",\"${INPUTDIR}\",\"wm_select.raw.root\"\);

#rm *.so *.d
