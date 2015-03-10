#! /bin/bash

# input ntuple directory
NTUPDIR=/scratch/klawhorn/EWKAnaStore/8TeV/Selection

# integrated luminosity for data
HFLUMI=0.018729
PIXLUMI=0.018479
LUMI=${HFLUMI}

root -l -q selectZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu\"\)
root -l -q selectWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu\"\)
root -l -q selectAntiWm.C+\(\"wm.conf\",\"${NTUPDIR}/AntiWmunu\"\)
root -l -q plotZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
root -l -q plotWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu/ntuples\",\"Wmunu\",${LUMI}\)

root -l -q selectZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee\",1\)
root -l -q selectWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu\",1\)
root -l -q selectAntiWe.C+\(\"we.conf\",\"${NTUPDIR}/AntiWenu\",1\)
root -l -q plotZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\)
root -l -q plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

rm *.so *.d
