#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/user/a/amarini/work/Bacon/Run2/wz_flat_76X_dressed

# integrated luminosity for data
LUMI=2318.3

root -l -q  plotZeeTheoryUnc.C+\(\"zeegen.conf\",\"${NTUPDIR}/ZeeGen/ntuples\",\"Zee\",${LUMI}\)

