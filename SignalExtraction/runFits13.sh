#! /bin/bash


# integrated luminosity for data
LUMI13=199.270 # Low PU 13 TeV 2017


# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_wRecoil_wStat_2G_eff
NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP
# NTUPLEDIR=/afs/cern.ch/user/s/sabrandt/work/public/TEST_Wenu_Corrected_v3
# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV_wRecoil
OUT=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/SignalExtraction

# root -l -q fitWlnu.C+\(\"${OUT}/TEST_Wmunu13_fixAll\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q fitWlnu.C+\(\"${OUT}/PRELIM_Wmunu_noRecBkg\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q fitWlnu.C+\(\"${OUT}/PRELIM_Wm\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q fitWlnu.C+\(\"${OUT}/PRELIM_Wenu_thyFix\",\"${NTUPLEDIR}\",\"Wenu\",${LUMI13},${LUMI13}\);
# root -l -q fitWlnu.C+\(\"${OUT}/TEST_Wmunu5_fixRecoil\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI5},${LUMI5}\);
# root -l -q fitWlikeZm.C+\(\"TEST_WlikeZm_Inclusive_3G_ptRW\",$LUMI13,$LUMI13\);

# # root -l -q plotResiduals.C+\(\"TEST_2019June18_CheckFixedStaUnc_2GausStat_Keys_Eta_unzoom\",${LUMI13}\); 
# root -l -q plotWshapes.C+\(\"${OUT}/TEST_PlotSWs\",${LUMI13}\); 

#rm *.so *.d
 