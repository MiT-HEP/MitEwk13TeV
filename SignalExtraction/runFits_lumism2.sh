#! /bin/bash

#ntuple directory
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/SM/newBacon/

# integrated luminosity for data
LUMI=2152
LUMI2=96



 EXTENSION=bar
 
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

LUMI2_0=11.775
LUMI2_1=12.249
LUMI2_2=11.201
LUMI2_3=10.794
LUMI2_4=9.319
LUMI2_5=7.901
LUMI2_6=6.432
LUMI2_7=7.424
LUMI2_8=8.785
LUMI2_9=5.262


 # root -l -q fitZm.C+\(\"Zmm_pf_lowPU13_pt25\",240,0\) #central
 # root -l -q fitZm.C+\(\"Zmm_full_ptandshitty\",$LUMI,0\) #central
 # root -l -q fitZm.C+\(\"Zmm_plot_lumi0_${EXTENSION}\",$LUMI_0,0,\"0\"\) #central
 # root -l -q fitZm.C+\(\"Zmm_plot_lumi1_${EXTENSION}\",$LUMI_1,0,\"1\"\) #central
 # root -l -q fitZm.C+\(\"Zmm_plot_lumi2_${EXTENSION}\",$LUMI_2,0,\"2\"\) #central
 # root -l -q fitZm.C+\(\"Zmm_plot_lumi3_${EXTENSION}\",$LUMI_3,0,\"3\"\) #central
 # root -l -q fitZm.C+\(\"Zmm_plot_lumi4_${EXTENSION}\",$LUMI_4,0,\"4\"\) #central
 # # root -l -q fitZm.C+\(\"Zmm_plot_lumi5_${EXTENSION}\",$LUMI_5,0,\"5\"\) #central
 # root -l -q fitZm.C+\(\"Zmm_plot_lumi6_${EXTENSION}\",$LUMI_6,0,\"6\"\) #central
 #root -l -q fitZm.C+\(\"Zmm_plot_g25_375\",$LUMI,0\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_Wlike_craprecoilRW\",$LUMI,0\) #central

# root -l -q compare10SegmentsMuons.C+\(\"Zmm_Wlike_sectionCompare\",$LUMI_0,$LUMI_0,0\) #central

 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi0_Wlike_${EXTENSION}\",${LUMI_0},${LUMI2_0},0,\"0\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi1_Wlike_${EXTENSION}\",${LUMI_1},${LUMI2_1},0,\"1\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi2_Wlike_${EXTENSION}\",${LUMI_2},${LUMI2_2},0,\"2\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi3_Wlike_${EXTENSION}\",${LUMI_3},${LUMI2_3},0,\"3\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi4_Wlike_${EXTENSION}\",${LUMI_4},${LUMI2_4},0,\"4\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi5_Wlike_${EXTENSION}\",${LUMI_5},${LUMI2_5},0,\"5\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi6_Wlike_${EXTENSION}\",${LUMI_6},${LUMI2_6},0,\"6\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi7_Wlike_${EXTENSION}\",${LUMI_7},${LUMI2_7},0,\"7\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi8_Wlike_${EXTENSION}\",${LUMI_8},${LUMI2_8},0,\"8\"\) #central
 # root -l -q fitWlikeZm.C+\(\"Zmm_lumi9_Wlike_${EXTENSION}\",${LUMI_9},${LUMI2_9},0,\"9\"\) #central

 root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi0_${EXTENSION}\",222.904,${LUMI2},0,\"0\"\) #central
 root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi1_${EXTENSION}\",223.971,${LUMI2},0,\"1\"\) #central
 root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi2_${EXTENSION}\",224.02,${LUMI2},0,\"2\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi3_${EXTENSION}\",215.887,${LUMI2},0,\"3\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi4_${EXTENSION}\",224.085,${LUMI2},0,\"4\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi5_${EXTENSION}\",223.632,7.901,0,\"5\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi6_${EXTENSION}\",145.584,6.432,0,\"6\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi7_${EXTENSION}\",224.058,7.424,0,\"7\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi8_${EXTENSION}\",224.008,8.785,0,\"8\"\) #central
 # root -l -q fitWm_lumis_2d.C+\(\"Wmunu_lumi9_${EXTENSION}\",224.172,5.262,0,\"9\"\) #central

#rm *.so *.d
