#!/bin/bash

# output ntuple directory
NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec_76x

# integrated luminosity for data
LUMI=2215

root -l -q selectZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee\",0\)
#root -l -q selectZeeGen.C+\(\"zeegen.conf\",\"${NTUPDIR}/ZeeGen\",0\)

rm *.so *.d
