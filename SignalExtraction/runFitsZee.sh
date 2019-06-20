#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_2017ID

# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV

# integrated luminosity for data
LUMI=199.2

#root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\,0\)
root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"TEST_Zee_ReRun_CheckYield_PrefireDown\",${LUMI}\,1\)

