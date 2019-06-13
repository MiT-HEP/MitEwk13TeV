#! /bin/bash


# integrated luminosity for data
# LUMI=213.1 # Low PU 13 TeV 2017
LUMI13=199.270 # Low PU 13 TeV 2017
LUMI5=294.4 # Low PU 5 TeV 2017 ele trigger
# LUMI5=291.1 # Low PU 5 TeV 2017 muon trigger

# LUMI2=199.270 # not used in Low PU atm

# root -l -q plotResiduals.C+\(\"Zmm_lumi9_Wlike_ErrorBand_all\",${LUMI_9},0\) 
#root -l -q plotResidualsEle.C+\(\"Zee_ErrorBand_Wlike_all_zoom\",${LUMI},0\) &
#sleep 100
#root -l -q fitWe.C+\(\"Wenu_default\",${LUMI},${LUMI},0\) #central

# root -l -q fitZm.C+\(\"Zmm_pf_lowPU13_pt25_rebin_fits_plotZeeMCwithZmmMC_lowPU_scEleCorr\",$LUMI,0\);
# root -l -q fitZm.C+\(\"Zmm_newNtuple_LepEffSFs\",$LUMI,0\);
# root -l -q fitZm.C+\(\"Zmm_13TeVMC_5TeVData_raw\",$LUMI,0\);
# root -l -q fitZe.C+\(\"Zee_13TeV_2017EleID\",$LUMI13,0\);

# root -l -q fitZm.C+\(\"TEST_Zmm_13TeV_incl_bigBin_v1\",$LUMI13,0\);
# root -l -q fitZm.C+\(\"TEST_Zmm_13TeV_2017ID_Eff_v2\",$LUMI13,0\);
# root -l -q fitWlikeZm.C+\(\"TEST_Zmm_Wlike_incl_v0\",$LUMI13,0\);

# root -l -q fitZm.C+\(\"Zmm_5TeV_Try2_NoPrefire\",$LUMI5,0\);
# root -l -q fitZe.C+\(\"Zee_13TeV_AllCorr_noPrefire_fixLUMI\",$LUMI13,0\);
# root -l -q fitZe.C+\(\"Zee_5TeV_Try2_WithPrefire_EE\",$LUMI5,0\);

# root -l -q fitWm.C+\(\"TEST_Wm13_Independent_Pepe1_constEWK_v2\",$LUMI13,$LUMI13,0\);

root -l -q fitWm_lumis_2d.C+\(\"TEST_Wm13_MET_IndepPepe2_gap02\",$LUMI13,$LUMI13,0\);
# root -l -q fitWm_lumis_2d.C+\(\"TEST_Wm13_mT_25_v2_full\",$LUMI13,$LUMI13,0\);
# root -l -q fitWe.C+\(\"Wenu_checkBkgShape_Fits_wSig\",$LUMI13,$LUMI13,0\);
# root -l -q plotResiduals.C+\(\"TEST_Zmm_bigBin_Stat_v_Centr\",$LUMI13,0\);

#rm *.so *.d
