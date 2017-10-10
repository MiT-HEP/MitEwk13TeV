#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2305
LUMI2=96

 # root -l -q fitWe.C+\(\"Wenu_newBacon_JUL5recoil_Zee_workspace\",${LUMI},${LUMI},0\)
 # mv ./Wenu_pdfTemplates.root Wenu_newBacon_JUL5recoil_Zee_workspace
 
# root -l -q fitWe_toys_dualinput.C+\(\"toys_Zmm_Vs_Zmm_test_a3free\",${LUMI},${LUMI},0\)

root -l -q fitWe_test_parfit.C+\(\"fitWe_test_RooCat\",${LUMI},${LUMI},0\)

# root -l -q drawTemplates.C+\(\"Wenu_toys_2_QCDPdf_Bin1.0GeV_pdfs\",${LUMI},${LUMI},0\)
#rm *.so *.d
