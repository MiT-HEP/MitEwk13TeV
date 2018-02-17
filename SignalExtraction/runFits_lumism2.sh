#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2305
LUMI2=96

# root -l -q drawTemplates.C+\(\"Wmunu_compareEfficiencies_sec1_8\",${LUMI},${LUMI},0\)

 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P0_0.5gev\",${LUMI},${LUMI2},0,\"0\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P0_0.5gev
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P1_0.5gev\",${LUMI},${LUMI2},0,\"1\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P1_0.5gev
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P2_0.5gev\",${LUMI},${LUMI2},0,\"2\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P2_0.5gev
# root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P3_0.5gev\",${LUMI},${LUMI2},0,\"3\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P3_0.5gev
 
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P4_0.5gev\",${LUMI},${LUMI2},0,\"4\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P4_0.5gev
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P5_0.5gev\",${LUMI},${LUMI2},0,\"5\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P5_0.5gev
 
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P6_0.5gev\",${LUMI},${LUMI},0,\"6\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P6_0.5gev
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P7_0.5gev\",${LUMI},${LUMI},0,\"7\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P7_0.5gev
 
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P8_0.5gev\",${LUMI},${LUMI},0,\"8\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P8_0.5gev
 # root -l -q fitWm_lumis.C+\(\"Wmunu_MedID_newBacon_10Sec_P9_0.5gev\",${LUMI},${LUMI},0,\"9\"\) #central
 # mv ./Wmunu_pdfTemplates.root Wmunu_MedID_newBacon_10Sec_P9_0.5gev
 
#  root -l -q fitWe_lumis.C+\(\"Wenu_fixEff_smallBin_lumi0\",${LUMI},${LUMI},0,\"0\"\) #central
#  root -l -q fitWe_lumis.C+\(\"Wenu_fixEff_smallBin_lumi1\",${LUMI},${LUMI},0,\"1\"\) #central
#  root -l -q fitWe_lumis.C+\(\"Wenu_fixEff_smallBin_lumi2\",${LUMI},${LUMI},0,\"2\"\) #central
 
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi0_13b\",${LUMI},${LUMI},0,\"0\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi1_13b\",${LUMI},${LUMI},0,\"1\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi2_13b\",${LUMI},${LUMI},0,\"2\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi3_13b\",${LUMI},${LUMI},0,\"3\"\) #central
  # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi4_13b\",${LUMI},${LUMI},0,\"4\"\) #central
 #  root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi5_13b\",${LUMI},${LUMI},0,\"5\"\) #central
   #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi6_13b\",${LUMI},${LUMI},0,\"6\"\) #central
 #  root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi7_13b\",${LUMI},${LUMI},0,\"7\"\) #central
   #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi8_13b\",${LUMI},${LUMI},0,\"8\"\) #central
  # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi9_13b\",${LUMI},${LUMI},0,\"9\"\) #central
  
  root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi0_16_extraBin_v3\",222.904,11.775,0,\"0\"\) #central
  root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi1_16_extraBin_v3\",223.971,12.249,0,\"1\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi2_16_extraBin\",224.02,11.201,0,\"2\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi3_16_extraBin\",215.887,10.794,0,\"3\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi4_16_extraBin\",224.085,9.319,0,\"4\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi5_16_extraBin\",223.632,7.901,0,\"5\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi6_16_extraBin\",145.584,6.432,0,\"6\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi7_16_extraBin\",224.058,7.424,0,\"7\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi8_16_extraBin\",224.008,8.785,0,\"8\"\) #central
  #root -l -q fitWm_lumis_2d.C+\(\"Wmunu_2d_lumi9_16_extraBin\",224.172,5.262,0,\"9\"\) #central
  
  # root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi1\",${LUMI},${LUMI},0,\"1\"\) #central
 # mv ./Wenu_pdfTemplates.root Wenu_2d_lumi1
 
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_MedID_newBacon_10Sec_P7_0.5gev_ZeeCorr_noconst\",${LUMI},${LUMI},0,\"7\"\) #central
 # mv ./Wenu_pdfTemplates.root Wenu_MedID_newBacon_10Sec_P7_0.5gev_ZeeCorr_noconst
 
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_MedID_newBacon_10Sec_P8_0.5gev_ZeeCorr_noconst\",${LUMI},${LUMI},0,\"8\"\) #central
 # mv ./Wenu_pdfTemplates.root Wenu_MedID_newBacon_10Sec_P8_0.5gev_ZeeCorr_noconst
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_MedID_newBacon_10Sec_P9_0.5gev_ZeeCorr_noconst\",${LUMI},${LUMI},0,\"9\"\) #central
 # mv ./Wenu_pdfTemplates.root Wenu_MedID_newBacon_10Sec_P9_0.5gev_ZeeCorr_noconst

#loop
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P0_0.5gev_constEWK\",${LUMI},${LUMI},0,\"0\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P1_0.5gev_constEWK\",${LUMI},${LUMI},0,\"1\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P2_0.5gev_constEWK\",${LUMI},${LUMI},0,\"2\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P3_0.5gev_constEWK\",${LUMI},${LUMI},0,\"3\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P4_0.5gev_constEWK\",${LUMI},${LUMI},0,\"4\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P5_0.5gev_consta1\",${LUMI},${LUMI},0,\"5\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P6_0.5gev_consta1\",${LUMI},${LUMI},0,\"6\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P7_0.5gev_consta1\",${LUMI},${LUMI},0,\"7\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P8_0.5gev_consta1\",${LUMI},${LUMI},0,\"8\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P9_0.5gev_consta1\",${LUMI},${LUMI},0,\"9\"\) #central
# root -l -q fitWe_lumis.C+\(\"Wenu_MedID_newBacon_10Sec_P10_0.5gev_allfree\",${LUMI},${LUMI},0,\"10\"\) #central

#rm *.so *.d
