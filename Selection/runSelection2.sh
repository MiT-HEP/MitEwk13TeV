#!/bin/bash

# output ntuple directory
#NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_wPrefire
NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV

# integrated luminosity for data
# LUMI=213.1 #doesn't matter here

# root -l -q selectWm.C+\(\"wm_13.conf\",\"${NTUPDIR}/Wmunu\",0,0,1\)
root -l -q selectAntiWm.C+\(\"wm_13.conf\",\"${NTUPDIR}/AntiWmunu\",0,0,1\)

# root -l -q selectWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu\",1,0,0,1\)
# root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu/AntiIso\",1,0,0,1\)
# root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu/AntiDEtaDPhi\",1,0,0,1\)

#rm *.so *.d
