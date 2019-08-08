#! /bin/bash


# integrated luminosity for data
# LUMI=213.1 # Low PU 13 TeV 2017
# LUMI13=199.270 # Low PU 13 TeV 2017
LUMI13=199.270 # Low PU 13 TeV 2017
LUMI5=294.4 # Low PU 5 TeV 2017 ele trigger
# LUMI5=291.1 # Low PU 5 TeV 2017 muon trigger

# LUMI2=199.270 # not used in Low PU atm

# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_wRecoil_wStat_2G_eff
# NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_wRecoil
NTUPLEDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV_wGM
OUT=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/SignalExtraction

# root -l -q fitWlnu.C+\(\"${OUT}/TEST_Wmunu13_fixRecoil_10th_mT\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI13},${LUMI13}\);
root -l -q fitWlnu.C+\(\"${OUT}/TEST_Wmunu5_test_fixedRoch_fixHist\",\"${NTUPLEDIR}\",\"Wmunu\",${LUMI5},${LUMI5}\);
# root -l -q fitWlikeZm.C+\(\"TEST_WlikeZm_Inclusive_3G_ptRW\",$LUMI13,$LUMI13\);

# root -l -q plotResiduals.C+\(\"TEST_2019June18_CheckFixedStaUnc_2GausStat_Keys_Eta_unzoom\",${LUMI13}\); 

#rm *.so *.d
