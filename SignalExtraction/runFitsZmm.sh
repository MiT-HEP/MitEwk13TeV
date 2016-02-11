#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X

# integrated luminosity for data
LUMI=2263

root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmm\",${LUMI}\,0\)
root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmm\",${LUMI}\,1\)

