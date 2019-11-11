#! /bin/bash


# integrated luminosity for data
LUMI5=291.1 # Low PU 5 TeV 2017 muon trigger


NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV_PP
OUT=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/SignalExtraction

# root -l -q fitWlnu.C+\(\"${OUT}/TEST_Wmunu5_fixRecoil_10th_mT\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI5},${LUMI5}\);
# root -l -q fitWlnu.C+\(\"${OUT}/TEST_Wmunu5_test_fixedRoch_fixHist\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI5},${LUMI5}\);

# root -l -q fitWlnu.C+\(\"${OUT}/PRELIM_Wmunu_5TeV_noRecBkg\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI5},${LUMI5}\);
root -l -q fitWlnu.C+\(\"${OUT}/PRELIM_Wenu_5TeV_iso_025-045\",\"${NTUPLEDIR}\",\"Wenu\",${LUMI5},${LUMI5}\);
# root -l -q fitWlikeZm.C+\(\"TEST_WlikeZm_Inclusive_3G_ptRW\",$LUMI5,$LUMI5\);

# root -l -q plotResiduals.C+\(\"TEST_2019June18_CheckFixedStaUnc_2GausStat_Keys_Eta_unzoom\",${LUMI5}\); 

#rm *.so *.d
