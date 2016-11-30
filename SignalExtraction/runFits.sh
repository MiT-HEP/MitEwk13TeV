#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2305
LUMI2=96

# root -l -q fitWm.C+\(\"testWmunu\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_Central_Charge_35GeV\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_Recoil_Inclusive\",${LUMI},${LUMI2},0\) #recoil inclusive
# root -l -q fitWm.C+\(\"Wmunu_Recoil_RooKeys\",${LUMI},${LUMI2},0\) #recoil RooKeys
# root -l -q fitWm.C+\(\"Wmunu_Pileup_Up\",${LUMI},${LUMI2},0\) #pileup Up
# root -l -q fitWm.C+\(\"Wmunu_Pileup_Down\",${LUMI},${LUMI2},0\) #pileup Down
# root -l -q fitWm.C+\(\"Wmunu_Ewk_Fix\",${LUMI},${LUMI2},0\) #pileup Down
# root -l -q fitWm.C+\(\"Wmunu_QCD_Free\",${LUMI},${LUMI2},0\) #pileup Down

root -l -q fitWe.C+\(\"Wenu_TestGaussEwk_35GeV\",${LUMI},${LUMI2},0\) #central

# root -l -q fitWe.C+\(\"Wenu_Central_Charge_35GeV\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"testWenu\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"Wenu_Central_Charge\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"Wenu_Recoil_Inclusive\",${LUMI},${LUMI2},0\) #recoil inclusive
# root -l -q fitWe.C+\(\"Wenu_Recoil_RooKeys\",${LUMI},${LUMI2},0\) #recoil RooKeys
# root -l -q fitWe.C+\(\"Wenu_Pileup_Up\",${LUMI},${LUMI2},0\) #pileup Up
# root -l -q fitWe.C+\(\"Wenu_Pileup_Down\",${LUMI},${LUMI2},0\) #pileup Down
 #################################################################3
 
#  root -l -q fitWm.C+\(\"Wmunu_pileup_simultaneous_check\",${LUMI},${LUMI2},0\)
# root -l -q fitWm.C+\(\"Wmunu_pileup_free_fixPtCut\",${LUMI},${LUMI2},0\)
# root -l -q fitWm.C+\(\"AntiWmunu_WshapeDiff2\",${LUMI},${LUMI2},0\)
# root -l -q fitBinsWm.C+\(\"Wmunu_etaBins_test\",${LUMI},${LUMI2},0\)
# root -l -q fitZm.C+\(\"Zmumu_WLike\",${LUMI},0\)
#root -l -q fitWm_mc.C+\(\"Wmunu_mc\",${LUMI},0\)
# root -l -q fitZm.C+\(\"Zmumu_PF_testToys\",${LUMI},0\)

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
