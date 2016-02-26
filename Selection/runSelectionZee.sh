#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X

# integrated luminosity for data
LUMI=2263

echo "-> Running Zee"
root -l -q selectZee.C+\(\"zee_eos.conf\",\"${NTUPDIR}/Zee\",0\)
echo "-> Running ZeeGen"
root -l -q selectZeeGen.C+\(\"zeegen_eos.conf\",\"${NTUPDIR}/ZeeGen\")

rm *.so *.d
