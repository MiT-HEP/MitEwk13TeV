#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2305
LUMI2=96

# root -l -q fitWm.C+\(\"print_Bins\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_QCDFree_Print\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"Wenu_PF_bkg\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_Puppi\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_Puppi_MakeWorkspace\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWmPt.C+\(\"Wmunu_Puppi_lepPt_QCDData\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe2.C+\(\"Wenu_Part2\",${LUMI},${LUMI2},0\) #central
root -l -q fitWm.C+\(\"Wmunu_Default\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"Wenu_plot_eta3\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_Central_print_sf\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWm.C+\(\"Wmunu_Recoil_Inclusive\",${LUMI},${LUMI2},0\) #recoil inclusive
# root -l -q fitWm.C+\(\"Wmunu_Recoil_RooKeys\",${LUMI},${LUMI2},0\) #recoil RooKeys
# root -l -q fitWm.C+\(\"Wmunu_Pileup_Up\",${LUMI},${LUMI2},0\) #pileup Up
# root -l -q fitWm.C+\(\"Wmunu_Pileup_Down\",${LUMI},${LUMI2},0\) #pileup Down
# root -l -q fitWm.C+\(\"Wmunu_Ewk_Fix\",${LUMI},${LUMI2},0\) #pileup Down
# root -l -q fitWm.C+\(\"Wmunu_QCD_Free\",${LUMI},${LUMI2},0\) #pileup Down

# root -l -q fitWe.C+\(\"Wenu_EleID_QCD\",${LUMI},${LUMI2},0\) #central

#  root -l -q fitWlikeZe.C+\(\"Zee_wLike_FixSel\",${LUMI},${LUMI2},0\) #central

# root -l -q fitWe.C+\(\"Wenu_Central_Charge_35GeV\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"testWenu\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"Wenu_Central_Charge\",${LUMI},${LUMI2},0\) #central
# root -l -q fitWe.C+\(\"Wenu_Recoil_Inclusive\",${LUMI},${LUMI2},0\) #recoil inclusive
# root -l -q fitWe.C+\(\"Wenu_Recoil_RooKeys\",${LUMI},${LUMI2},0\) #recoil RooKeys
# root -l -q fitWe.C+\(\"Wenu_Pileup_Up\",${LUMI},${LUMI2},0\) #pileup Up
# root -l -q fitWe.C+\(\"Wenu_Pileup_Down\",${LUMI},${LUMI2},0\) #pileup Down
 #################################################################3

# root -l -q checkFP.C+\(\"checkFP_test\",${LUMI},0\) #central

#rm *.so *.d
