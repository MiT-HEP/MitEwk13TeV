#!/bin/bash
## Argument format is selectZee({conf},{output directory},{do scale correction},{do pu weight},{is13TeV},{number segments},{ith segment})
#########################################################################
##         13 TeV
## ----------------------------------------------------------------------
# # output ntuple directory
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV
NSEC=1
ITH=0

root -l -q selectZmm.C+\(\"zmm_13.conf\",\"${NTUPDIR}/Zmumu\",0,0,1,${NSEC},${ITH}\)


#########################################################################
##         5 TeV
## ----------------------------------------------------------------------

NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_5TeV
NSEC=5
ITH=1

root -l -q selectZmm.C+\(\"zmm_5.conf\",\"${NTUPDIR}/Zmumu\",0,0,0,${NSEC},${ITH}\)
