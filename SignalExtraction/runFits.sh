#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=212 # Low PU 13 TeV 2017
LUMI2=90.5 # not used in Low PU atm

LUMI_0=222.904
LUMI_1=223.971
LUMI_2=224.02
LUMI_3=215.887
LUMI_4=224.085
LUMI_5=223.632
LUMI_6=145.584
LUMI_7=224.058
LUMI_8=224.008
LUMI_9=224.172

# root -l -q plotResiduals.C+\(\"Zmm_lumi0_Wlike_ErrorBand_all\",${LUMI_0},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi1_Wlike_ErrorBand_all\",${LUMI_1},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi2_Wlike_ErrorBand_all\",${LUMI_2},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi3_Wlike_ErrorBand_all\",${LUMI_3},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi4_Wlike_ErrorBand_all\",${LUMI_4},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi5_Wlike_ErrorBand_all\",${LUMI_5},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi6_Wlike_ErrorBand_all\",${LUMI_6},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi7_Wlike_ErrorBand_all\",${LUMI_7},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi8_Wlike_ErrorBand_all\",${LUMI_8},0\) 
# root -l -q plotResiduals.C+\(\"Zmm_lumi9_Wlike_ErrorBand_all\",${LUMI_9},0\) 
#root -l -q plotResidualsEle.C+\(\"Zee_ErrorBand_Wlike_all_zoom\",${LUMI},0\) &
#sleep 100
#root -l -q fitWe.C+\(\"Wenu_default\",${LUMI},${LUMI},0\) #central

# root -l -q fitZm.C+\(\"Zmm_pf_lowPU13_pt25_rebin_fits_plotZeeMCwithZmmMC_lowPU_scEleCorr\",$LUMI,0\);
root -l -q fitZm.C+\(\"Zee_pf_lowPU13_pt25_zee_relval_fixfilter_wDAta\",$LUMI,0\);
# root -l -q fitZe.C+\(\"Zee_pf_lowPU13_pt25_rebin_fits_Rel\",$LUMI,0\);

#rm *.so *.d
