#!/bin/bash

# output ntuple directory
#NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec
NTUPDIR=/data/t3home000/sabrandt/2018_12_02_13TeVlowPU_v1/

# integrated luminosity for data
LUMI=212

#root -l -q selectZmm.C+g\(\"zmm_eos.conf\",\"${NTUPDIR}/Zmumu\",0\)
root -l -q selectZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu\",0\)
#root -l -q selectZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen\",0\)
#root -l -q selectWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu\",0\)
#root -l -q selectAntiWm.C+\(\"wm.conf\",\"${NTUPDIR}/AntiWmunu\"\)
#root -l -q rootlogon.plot.C plotZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q rootlogon.plot.C plotWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu/ntuples\",\"Wmunu\",${LUMI}\)

#root -l -q selectZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee\",0\)
#root -l -q selectWe.C+\(\"zee.conf\",\"${NTUPDIR}/Wenu\",0\)
#root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu\",0\)
#root -l -q rootlogon.plot.C plotZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\)
#root -l -q rootlogon.plot.C plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

#rm *.so *.d
