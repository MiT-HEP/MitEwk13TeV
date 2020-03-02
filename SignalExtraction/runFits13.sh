#! /bin/bash


# integrated luminosity for data
LUMI13=199.270 # Low PU 13 TeV 2017


# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_wRecoil_wStat_2G_eff
# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP
NTUPLEDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017_PP_noStat
# NTUPLEDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v4_fixHist_PP_noStat
# NTUPLEDIR=/afs/cern.ch/user/s/sabrandt/work/public/TEST_Wenu_Corrected_v3
# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV_wRecoil
OUT=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/SignalExtraction

# root -l -q -b fitWlnu.C+\(\"${OUT}/TEST_Wmunu13_fixAll\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu.C+\(\"${OUT}/PRELIM_Wmunu_noRecBkg\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu.C+\(\"${OUT}/PRELIM_Wm\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu.C+\(\"${OUT}/PRELIM_Wmunu_pt30\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu.C+\(\"${OUT}/TEST_Wmunu_mt_mt50_origIso\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu_test.C+\(\"${OUT}/TEST_Wenu13_IsoControl_Barrel_pt30\",\"${NTUPLEDIR}\",\"Wenu\",${LUMI13},${LUMI13}\);
root -l -q -b fitWlnu.C+\(\"${OUT}/TEST_NewEleNtuples_fullMT\",\"${NTUPLEDIR}\",\"Wenu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu_Aram.C+\(\"${OUT}/PRELIM_Wenu_newTrigger_testAram\",\"${NTUPLEDIR}\",\"Wenu\",${LUMI13},${LUMI13}\);
# root -l -q -b fitWlnu.C+\(\"${OUT}/TEST_Wmunu5_fixRecoil\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI5},${LUMI5}\);
# root -l -q -b fitWlikeZm.C+\(\"TEST_WlikeZm_Inclusive_3G_ptRW\",$LUMI13,$LUMI13\);

# # root -l -q -b plotResiduals.C+\(\"TEST_2019June18_CheckFixedStaUnc_2GausStat_Keys_Eta_unzoom\",${LUMI13}\); 
# root -l -q -b plotWshapes.C+\(\"${OUT}/TEST_PlotSWs\",${LUMI13}\); 

#rm *.so *.d
 