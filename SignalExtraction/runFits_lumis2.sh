#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2152
LUMI2=96

# EXTENSION=met_pt30to45GeV_central
EXTENSION=bar_35
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


# root -l -q compare10Segments.C+\(\"Zee_Wlike_sectionCompare\",$LUMI_0,$LUMI_0,0\) #central


# root -l -q fitWlikeZe.C+\(\"Zee_Wlike_full_WepRecoil\",$LUMI,$LUMI,0\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi0_Wlike_${EXTENSION}\",$LUMI_0,$LUMI_0,0,\"0\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi1_Wlike_${EXTENSION}\",$LUMI_1,$LUMI_1,0,\"1\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi2_Wlike_${EXTENSION}\",$LUMI_2,$LUMI_2,0,\"2\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi3_Wlike_${EXTENSION}\",$LUMI_3,$LUMI_3,0,\"3\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi4_Wlike_${EXTENSION}\",$LUMI_4,$LUMI_4,0,\"4\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi5_Wlike_${EXTENSION}\",$LUMI_5,$LUMI_5,0,\"5\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi6_Wlike_${EXTENSION}\",$LUMI_6,$LUMI_6,0,\"6\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi7_Wlike_${EXTENSION}\",$LUMI_7,$LUMI_7,0,\"7\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi8_Wlike_${EXTENSION}\",$LUMI_8,$LUMI_8,0,\"8\"\) #central
# root -l -q fitWlikeZe.C+\(\"Zee_lumi9_Wlike_${EXTENSION}\",$LUMI_9,$LUMI_9,0,\"9\"\) #central
 # root -l -q fitZe.C+\(\"Zee_plot_diag2\",$LUMI,0\) #central
 #root -l -q fitZe.C+\(\"Zee_lumi0_met_full\",222.904,0\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi0_mt_incl_mt50_sieie\",$LUMI_0,$LUMI_0,0,\"0\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi0_mt_50GeV_08035_inclBkg\",222.904,${LUMI},0,\"0\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_full_test_iso_smalbkg_pep_eff0_v2\",$LUMI,$LUMI,0,\"0\"\) #central 
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi0_${EXTENSION}\",$LUMI_0,$LUMI_0,0,\"0\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi1_${EXTENSION}\",$LUMI_1,$LUMI_1,0,\"1\"\)  #central
 root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi2_${EXTENSION}\",$LUMI_2,$LUMI_2,0,\"2\"\) #central
 root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi3_${EXTENSION}\",$LUMI_3,$LUMI_3,0,\"3\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi4_${EXTENSION}\",$LUMI_4,$LUMI_4,0,\"4\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi5_${EXTENSION}\",223.632,223.632,0,\"5\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi6_${EXTENSION}\",145.584,145.584,0,\"6\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi7_${EXTENSION}\",224.058,224.058,0,\"7\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi8_${EXTENSION}\",224.008,224.008,0,\"8\"\) #central
 # root -l -q fitWe_lumis_2d.C+\(\"Wenu_lumi9_${EXTENSION}\",224.172,224.172,0,\"9\"\) #central


#rm *.so *.d
