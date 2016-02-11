#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X

# integrated luminosity for data
LUMI=2263

root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\,0\)
root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\,1\)

