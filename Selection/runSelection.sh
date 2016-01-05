#!/bin/bash

# output ntuple directory
NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec

# integrated luminosity for data
LUMI=2215

root -l -q selectZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu\",0\)
#root -l -q selectWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu\"\)
#root -l -q selectAntiWm.C+\(\"wm.conf\",\"${NTUPDIR}/AntiWmunu\"\)
#root -l -q rootlogon.plot.C plotZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q rootlogon.plot.C plotWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu/ntuples\",\"Wmunu\",${LUMI}\)

#root -l -q selectZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee\",1\)
#root -l -q selectWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu\",0\)
#root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu\",0\)
#root -l -q rootlogon.plot.C plotZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\)
#root -l -q rootlogon.plot.C plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

rm *.so *.d
