#!/bin/bash

# output ntuple directory

# NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_v5
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017

NSEC=5
ITH=4

root -l -q selectZee.C+\(\"zee_13.conf\",\"${NTUPDIR}/Zee\",1,0,0,1,${NSEC},${ITH}\)
# root -l -q selectZee.C+\(\"zee_5.conf\",\"${NTUPDIR}/Zee\",1,0,0,0\)
# root -l -q selectZee.C+\(\"zee_13_minlo.conf\",\"${NTUPDIR}/Zee_minlo\",1,0,0,1\)
