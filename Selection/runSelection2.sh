#!/bin/bash

# output ntuple directory
#NTUPDIR=/data/blue/Bacon/Run2/wz_flat_diffxsec
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_13TeV_v4
NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_5TeV_JSONv1
# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV
NSEC=1
ITH=0
# integrated luminosity for data
# LUMI=213.1 # doesn't matter here

# root -l -q -b selectZmm.C+g\(\"zmm_eos.conf\",\"${NTUPDIR}/Zmumu\",0\) 
# root -l -q -b selectZmm.C+\(\"zmm_13.conf\",\"${NTUPDIR}/Zmumu\",0,0,1\)
# root -l -q -b selectZmm.C+\(\"zmm_5.conf\",\"${NTUPDIR}/Zmumu\",0,0,0\)
# root -l -q -b selectZmm.C+\(\"zmm_13_minlo.conf\",\"${NTUPDIR}/Zmumu_minlo\",0,0,1\)
#root -l -q -b selectZmmGen.C+\(\"zmmgen.conf\",\"${NTUPDIR}/ZmumuGen\",0\)
# root -l -q -b selectWm.C+\(\"wm_13.conf\",\"${NTUPDIR}/Wmunu\",0,0,1,${NSEC},${ITH}\)
# root -l -q -b selectAntiWm.C+\(\"wm_13.conf\",\"${NTUPDIR}/AntiWmunu\",0,0,1,${NSEC},${ITH}\)
#root -l -q -b rootlogon.plot.C plotZmm.C+\(\"zmm.conf\",\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q -b rootlogon.plot.C plotWm.C+\(\"wm.conf\",\"${NTUPDIR}/Wmunu/ntuples\",\"Wmunu\",${LUMI}\)

# root -l -q -b selectZee.C+\(\"zee_13.conf\",\"${NTUPDIR}/Zee\",1,1,0,1\)
# root -l -q -b selectZee.C+\(\"zee_13_minlo.conf\",\"${NTUPDIR}/Zee_minlo\",1,0,0,1\)
# root -l -q -b selectWe.C+\(\"we_13.conf\",\"${NTUPDIR}/Wenu\",1,0,0,1,${NSEC},${ITH}\)
root -l -q -b selectAntiWe.C+\(\"we_13.conf\",\"${NTUPDIR}/AntiWenu\",1,0,0,1,${NSEC},${ITH}\)
# root -l -q -b selectAntiWe.C+\(\"we_13.conf\",\"${NTUPDIR}/AntiWenu_sieie\",1,0,0,1\)

# root -l -q -b selectAntiWe.C+\(\"we0_13.conf\",\"${NTUPDIR}/AntiWenu_dEta_dPhi\",1,0,0,1,${NSEC},${ITH}\)
# root -l -q -b selectAntiWe.C+\(\"wx0_13.conf\",\"${NTUPDIR}/AntiWenu_dEta_dPhi\",1,0,0,1,${NSEC},${ITH}\)
#root -l -q -b rootlogon.plot.C plotZee.C+\(\"zee.conf\",\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\)
#root -l -q -b rootlogon.plot.C plotWe.C+\(\"we.conf\",\"${NTUPDIR}/Wenu/ntuples\",\"Wenu\",${LUMI}\)

#rm *.so *.d
