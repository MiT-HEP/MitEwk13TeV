#!/bin/bash

# output ntuple directory
#NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Fixed_noR9
NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_5TeV
 
# integrated luminosity for data
LUMI=291.1 # with the trigger requirement in Brilcalc
# LUMI=303.7 # without the trigger requirement in Brilcalc

# root -l -q selectZmm.C+g\(\"zmm_eos.conf\",\"${NTUPDIR}/Zmumu\",0\)
# root -l -q selectZmm.C+\(\"zmm_5.conf\",\"${NTUPDIR}/Zmumu\",0,0,0\)
root -l -q selectWm.C+\(\"wm_5.conf\",\"${NTUPDIR}/Wmunu\",0,0,0\)
# root -l -q selectAntiWm.C+\(\"wm_5.conf\",\"${NTUPDIR}/AntiWmunu\",0,0,0\)

# root -l -q selectZee.C+\(\"zee_5.conf\",\"${NTUPDIR}/Zee\",1,0,0\)
# root -l -q selectWe.C+\(\"zee.conf\",\"${NTUPDIR}/Wenu\",0\)
#root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu\",0\)
#root -l -q rootlogon.plot.C plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

#rm *.so *.d
