#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency

# integrated luminosity for data
LUMI=2318.3

root -l -q  plotZee.C+\(\"${NTUPDIR}/Ele/ntuples\",\"Zee\",${LUMI}\,0\)
root -l -q  plotZee.C+\(\"${NTUPDIR}/Ele/ntuples\",\"Zee\",${LUMI}\,1\)

