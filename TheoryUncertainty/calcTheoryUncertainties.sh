#! /bin/bash

#ntuple directory
NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec

# integrated luminosity for data
LUMI=2215

root -l -q  plotZmmTheoryUnc.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmumu\",${LUMI}\)

