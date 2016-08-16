#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2308
LUMI2=96

#  root -l -q fitWm.C+\(\"Wmunu_pileup_simultaneous_check\",${LUMI},${LUMI2},0\)
# root -l -q fitWm.C+\(\"Wmunu_pileup_free_fixPtCut\",${LUMI},${LUMI2},0\)
# root -l -q fitWm.C+\(\"AntiWmunu_WshapeDiff2\",${LUMI},${LUMI2},0\)
# root -l -q fitBinsWm.C+\(\"Wmunu_etaBins_test\",${LUMI},${LUMI2},0\)
# root -l -q fitZm.C+\(\"Zmumu_WLike\",${LUMI},0\)
#root -l -q fitWm_mc.C+\(\"Wmunu_mc\",${LUMI},0\)
# root -l -q fitZm.C+\(\"Zmumu_PF_testToys\",${LUMI},0\)

root -l -q fitWe.C+\(\"Wenu_test_pepe1_free\",${LUMI},${LUMI},0\)
#root -l -q fitZee.C+\(\"Zee\",${LUMI},0\)
#root -l -q fitZee2.C+\(\"May23/Zee2\",${LUMI},0\)

# root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\,0\)
# root -l -q  plotZmm.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\,1\)
# root -l -q  plotZmmResScaleUncert.C+\(\"${NTUPDIR}/Zmumu/ntuples\",\"Zmumu\",${LUMI}\)
#root -l -q   plotZee.C+\(\"Zee\",${LUMI}\)

# root -l -q postFitWe.C+\(\"Wenu_post_simultaneous_fixLumi_fixPlots_fixCombine\",${LUMI},0\)
# root -l -q postFitWm.C+\(\"Wmunu_post_free_fixPtCut\",${LUMI},0\)
# root -l -q postFitTestCombine.C+\(\"Wmunu_test\",${LUMI},0\)
# root -l -q postFitWm.C+\(\"WmunuPF_post_WShapeDiff\",${LUMI},0\)
#rm *.so *.d
