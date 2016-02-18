#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X_nominal

# integrated luminosity for data
LUMI=2263

root -l -q  plotZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmm\",${LUMI}\)
root -l -q  plotZmmGenResScaleUncert.C+\(\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmm\",${LUMI}\)

