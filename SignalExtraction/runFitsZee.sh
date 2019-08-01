#! /bin/bash

#ntuple directory
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/TEST_Zee_ScaleUp5sig
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV_ele
NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV_looseID_raw_rmTrig
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV_ele_ecalEnergy
# NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV
OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction

# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV

# integrated luminosity for data
LUMI=199.2

#root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\,0\)
root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_Pt32_raw_loose_noTrig_noCat\",${LUMI}\,1\)

