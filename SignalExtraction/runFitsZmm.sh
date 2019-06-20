#! /bin/bash

#ntuple directory
# NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X_newpu
NTUPDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV

# integrated luminosity for data
LUMI13=199.270

root -b -l -q fitZm.C+\(\"TEST_Zmm_StatFixed_PlotCheck_3Gaus\",${LUMI13}\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"TEST_Zmm_ReRun_CheckYield_noPrefire\",${LUMI13}\,0\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"checkTOP\",${LUMI13}\,0\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"Zmm\",${LUMI13}\,1\)
# root -l -q  plotZmmResScaleUncert.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"Zmm\",${LUMI13}\)
