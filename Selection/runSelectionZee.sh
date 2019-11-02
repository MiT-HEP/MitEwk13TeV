#!/bin/bash

# output ntuple directory

NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV
# NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV_Raw

root -l -q selectZee.C+\(\"zee_13.conf\",\"${NTUPDIR}/Zee\",1,0,0,1\)
# root -l -q selectZee.C+\(\"zee_5.conf\",\"${NTUPDIR}/Zee\",0,0,0,0\)
# root -l -q selectZee.C+\(\"zee_13_minlo.conf\",\"${NTUPDIR}/Zee_minlo\",1,0,0,1\)
