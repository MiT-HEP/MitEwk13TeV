#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X_v2

# integrated luminosity for data
LUMI=2263

root -l -q selectZee.C+\(\"zee_eos.conf\",\"${NTUPDIR}/Zee\",0\)
#root -l -q selectZeeGen.C+\(\"zeegen_eos.conf\",\"${NTUPDIR}/ZeeGen\",0\)

rm *.so *.d
