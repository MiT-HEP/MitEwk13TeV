#! /bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/work/j/jlawhorn/public/wz-ntuples/

# integrated luminosity for data
HFLUMI=0.001
PIXLUMI=0.001
LUMI=${HFLUMI}

#root -l -q selectZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu\"\)
#root -l -q selectWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu\"\)
#root -l -q selectAntiWm.C+\(\"wm.conf\",\"${NTUPDIR}/AntiWmunu\"\)
#root -l -q rootlogon.plot.C plotZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q rootlogon.plot.C plotWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu/ntuples\",\"Wmunu\",${LUMI}\)

root -l -q selectZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee\",0\)
#root -l -q selectWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu\",0\)
#root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu\",0\)
#root -l -q rootlogon.plot.C plotZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\)
#root -l -q rootlogon.plot.C plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

rm *.so *.d
