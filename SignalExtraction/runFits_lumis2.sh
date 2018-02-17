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
 
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_incl\",228.821,228.821,0,\"0\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi0_12_newEff_m1_v2\",228.821,228.821,0,\"0\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi1_12_newEff_m1_v2\",229.95,229.95,0,\"1\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi2_12_newEff_m1_v2\",230.001,230.001,0,\"2\"\) #central
   #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi3_12_newEff_m1_v2\",221.651,221.651,0,\"3\"\) #central
 #  root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi4_12_newEff_m1_v2\",230.067,230.067,0,\"4\"\) #central
 #  root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi5_12_newEff_m1_v2\",229.602,229.602,0,\"5\"\) #central
   #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi6_12_newEff_m1_v2\",145.584,145.584,0,\"6\"\) #central
 #  root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi7_12_newEff_m1_v2\",230.04,230.04,0,\"7\"\) #central
   #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi8_12_newEff_m1_v2\",229.989,229.989,0,\"8\"\) #central
  # root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi9_12_newEff_m1_v2\",230.156,230.156,0,\"9\"\) #central
 
 root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi0_16_1ewk_extraBin_v4_pQ\",222.904,222.904,0,\"0\"\) #central
 root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi1_16_1ewk_extraBin_v4_pQ\",223.971,223.971,0,\"1\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi2_16_1ewk_extraBin_v2\",224.02,224.02,0,\"2\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi3_16_1ewk_extraBin_v2\",215.887,215.887,0,\"3\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi4_16_1ewk_extraBin_v2\",224.085,224.085,0,\"4\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi5_16_1ewk_extraBin_v2\",223.632,223.632,0,\"5\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi6_16_1ewk_extraBin_v2\",145.584,145.584,0,\"6\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi7_16_1ewk_extraBin_v2\",224.058,224.058,0,\"7\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi8_16_1ewk_extraBin_v2\",224.008,224.008,0,\"8\"\) #central
 #root -l -q fitWe_lumis_2d.C+\(\"Wenu_2d_lumi9_16_1ewk_extraBin_v2\",224.172,224.172,0,\"9\"\) #central
  
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
