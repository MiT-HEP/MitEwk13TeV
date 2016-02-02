#!/bin/bash

# output ntuple directory
NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec

# integrated luminosity for data
LUMI=2215

root -l -q selectZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu\",0\)
root -l -q selectZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen\",0\)

rm *.so *.d
