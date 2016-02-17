#! /bin/bash

#ntuple directory
#NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec
NTUPDIR=/afs/cern.ch/user/a/amarini/work/Bacon/Run2/wz_flat_76X

# integrated luminosity for data
LUMI=2263

root -l -q  plotZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmumu\",${LUMI}\)
root -l -q  plotZmmGenResScaleUncert.C+\(\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmumu\",${LUMI}\)

