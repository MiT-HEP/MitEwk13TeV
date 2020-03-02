#!/bin/bash

## Argument format is selectZee({conf},{output directory},{do scale correction},{nsigma},{do pu weight},{is13TeV},{number segments},{ith segment})

######################################################################################
##         13 TeV
## -----------------------------------------------------------------------------------
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV
NSEC=1
ITH=0

root -l -q selectZee.C+\(\"zee_13.conf\",\"${NTUPDIR}/Zee\",1,0,0,1,${NSEC},${ITH}\)

######################################################################################

######################################################################################
##         5 TeV
## -----------------------------------------------------------------------------------
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_5TeV
NSEC=1
ITH=0

root -l -q selectZee.C+\(\"zee_5.conf\",\"${NTUPDIR}/Zee\",1,0,0,0,${NSEC},${ITH}\)