#! /bin/bash

#ntuple directory

# NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV
NTUPDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_5TeV
# NTUPDIR=/eos/cms/store/user/sabrandt/StandardModel/Ntuples2017GH/LowPU2017ID_13TeV
OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction

# NTUPDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU2017ID_13TeV

# integrated luminosity for data
LUMI=199.2
LUMI5=291.1
S13="13TeV"
S5="5TeV"

#root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"Zee\",${LUMI}\,0\)
# root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_LooseNoTrig_RawCorr_Compare_txt_SCALE_SMEAR_MINE\",${LUMI}\,1\)
# root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_corrEcalE_all_SCCorr_startRaw\",${LUMI}\,1\)
# root -l -q  plotZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_13TeVZee_RootOutput\",\"${S13}\",${LUMI}\,1\)
# root -l -q  fitZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/FINAL_Zee_13TeV_noUncBand\",\"${S13}\",${LUMI}\,1\)
# root -l -q  testZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_zee_5TeV_checkPlots\",\"${S5}\",${LUMI5}\,1\)
root -l -q  fitZee.C+\(\"${NTUPDIR}/Zee/ntuples\",\"${OUTDIR}/TEST_Zee_5TeV_wEffs13\",\"${S5}\",${LUMI5}\,1\)

