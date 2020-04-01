#!/bin/bash

## Electron format is selectWe({conf},{output directory},{do scale correction},{nsigma},{do pu weight},{is13TeV},{number segments},{ith segment})

## Muon format is selectWm({conf},{output directory},{do scale correction},{do pu weight},{is13TeV},{number segments},{ith segment})

##############################################################################
##                    13 TeV
##############################################################################
# output ntuple directory
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_NEWPROD_TEST
NSEC=100
ITH=0

# root -l -q selectZmm.C+\(\"zmm_13.conf\",\"${NTUPDIR}/Zmumu\",0,0,1,${NSEC},${ITH}\)
# root -l -q -b selectWm.C+\(\"wm_13.conf\",\"${NTUPDIR}/Wmunu\",0,0,1,${NSEC},${ITH}\)
# root -l -q -b selectAntiWm.C+\(\"wm_13.conf\",\"${NTUPDIR}/AntiWmunu\",0,0,1,${NSEC},${ITH}\)

# root -l -q selectZee.C+\(\"zee_13.conf\",\"${NTUPDIR}/Zee\",1,0,0,1,${NSEC},${ITH}\)
# root -l -q -b selectWe.C+\(\"we_13.conf\",\"${NTUPDIR}/Wenu\",1,0,0,1,${NSEC},${ITH}\)
root -l -q -b selectAntiWe.C+\(\"we_13.conf\",\"${NTUPDIR}/AntiWenu\",1,0,0,1,${NSEC},${ITH}\)



# root -l -q -b selectEMu.C+\(\"wm_13.conf\",\"${NTUPDIR}/SelEmu_SingleMuon\",1,0,1,${NSEC},${ITH}\)

# ##############################################################################
# ##                    5 TeV
# ##############################################################################
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_5TeV_NEWPROD

NSEC=3
ITH=2

# root -l -q selectZmm.C+\(\"zmm_5.conf\",\"${NTUPDIR}/Zmumu\",0,0,0,${NSEC},${ITH}\)
# root -l -q -b selectWm.C+\(\"wm_5.conf\",\"${NTUPDIR}/Wmunu\",0,0,0,${NSEC},${ITH}\)
# root -l -q -b selectAntiWm.C+\(\"wm_5.conf\",\"${NTUPDIR}/AntiWmunu\",0,0,0,${NSEC},${ITH}\)

# root -l -q selectZee.C+\(\"zee_5.conf\",\"${NTUPDIR}/Zee\",1,0,0,0,${NSEC},${ITH}\)
# root -l -q -b selectWe.C+\(\"we_5.conf\",\"${NTUPDIR}/Wenu\",1,0,0,0,${NSEC},${ITH}\)
# root -l -q -b selectAntiWe.C+\(\"we_5.conf\",\"${NTUPDIR}/AntiWenu\",1,0,0,0,${NSEC},${ITH}\)


# #rm *.so *.d
