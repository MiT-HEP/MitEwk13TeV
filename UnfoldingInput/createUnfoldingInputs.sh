#! /bin/bash

# integrated luminosity for data
LUMI=2215

root -l -q  plotZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmumu\",${LUMI}\)
root -l -q  plotZmmGenResScaleUncert.C+\(\"${NTUPDIR}/ZmumuGen/ntuples\",\"Zmumu\",${LUMI}\)

