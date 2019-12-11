#! /bin/bash

#ntuple directory

# NTUPDIR_13=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV_PP
# NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV
NTUPDIR_13=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV_v5_EleMedID2017
OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction

# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV

# integrated luminosity for data
LUMI13=199.2
LUMI5=291.1
S13="13TeV"
S5="5TeV"

# root -b -l -q fitWlikeZe.C+\(\"${OUTDIR}/PLOT_Ze_Wlike_We_pos_13TeV_noEleCorr\",\"${NTUPDIR_13}\",\"${S13}\",${LUMI13},0\)
#root -b -l-q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\,0\)
# root -b -l-q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_LooseNoTrig_RawCorr_Compare_txt_SCALE_SMEAR_MINE\",${LUMI}\,1\)

# # root -b -l -q fitZe.C+\(\"${NTUPDIR_13}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_MET_MedID_newEff_${S13}\",\"${S13}\",${LUMI13},0\) 
# root -b -l-q  plotZee.C+\(\"${NTUPDIR_13}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_MET_MedID_newEff_${S13}\",${LUMI}\,1\)
# root -b -l-q  plotZee.C+\(\"${NTUPDIR_13}/Zee/ntuples\",\"${OUTDIR}/TEST_NewTrigger\",\"${S13}\",${LUMI}\,1\)
root -b -l-q  fitZee.C+\(\"${NTUPDIR_13}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_MET_MedID_newEff_${S13}\",\"${S13}\",${LUMI13}\,1\)
# root -b -l-q  testZee.C+\(\"${NTUPDIR_13}/Zee/ntuples\",\"${OUTDIR}/PRELIM_Zee_5TeV_newHLT\",\"${S5}\",${LUMI5}\,1\)
# root -b -l-q  fitZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_5TeV_wEffs13_tryAgain\",\"${S5}\",${LUMI5}\,1\)

# root -b -l-q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_5_check\",\"${S5}\",${LUMI5}\,1\)
