#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/user/a/amarini/work/Bacon/Run2/wz_flat_76X_newpu

# integrated luminosity for data
LUMI=2318.3

root -l -q  plotZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmm\",${LUMI}\)
root -l -q  plotZmmGenResScaleUncert.C+\(\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmm\",${LUMI}\)

