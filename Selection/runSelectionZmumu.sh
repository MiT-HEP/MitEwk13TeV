#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X

# integrated luminosity for data
LUMI=2263

echo "->Running Zmm"
root -l -q selectZmm.C+\(\"zmm_eos.conf\",\"${NTUPDIR}/Zmumu\",0\)
echo "->Running ZmmGen"
root -l -q selectZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen\",0\)

rm *.so *.d
