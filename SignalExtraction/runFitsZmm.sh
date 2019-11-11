#! /bin/bash

#ntuple directory
# NTUPDIR=/afs/cern.ch/work/a/amarini/Bacon/Run2/wz_flat_76X_newpu
# NTUPDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_13TeV
NTUPDIR_13=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_13TeV
NTUPDIR_5=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/LowPU2017ID_5TeV
OUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction

# integrated luminosity for data
LUMI13=199.270
LUMI5=291.107
S13="13TeV"
S5="5TeV"

# root -b -l -q fitWlikeZm.C+\(\"${OUTDIR}/PLOT_ZmWlike_Wm_${S13}_RecoilClosure\",\"${NTUPDIR_13}\",\"${S13}\",${LUMI13},0\)
# root -b -l -q fitZm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"${OUTDIR}/FINAL_ZmmMET_rmEEjet_keys_${S13}\",\"${S13}\",${LUMI13},0\)
# root -b -l -q fitZm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"${OUTDIR}/TEST_print_${S13}\",\"${S13}\",${LUMI13},0\)
root -b -l -q  fitZmm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"${OUTDIR}/FINAL_Zmm_13TeV_v3_adjSta\",\"${S13}\",${LUMI13}\,0\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_5}/Zmumu/ntuples\",\"${OUTDIR}/FINAL_Zmm_5TeV\",\"${S5}\",${LUMI5}\,0\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_5}/Zmumu/ntuples\",\"${OUTDIR}/plotZmm5TeV_2017Roch_noBKG\",${LUMI5}\,0\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"checkTOP\",${LUMI13}\,0\)
# root -b -l -q  plotZmm.C+\(\"${NTUPDIR_13}/Zmumu/ntuples\",\"Zmm\",${LUMI13}\,1\)
# root -l -q  plotResiduals.C+\(\"${OUTDIR}/Zmm_METenv_all\",${LUMI13}\)

### 
# 5TeV 

# root -b -l -q fitZm.C+\(\"${NTUPDIR_5}/Zmumu/ntuples\",\"${OUTDIR}/PRELIM_Zmm_5TeV_Recoil\",\"${S5}\",${LUMI5},0\)

# root -b -l -q  fitZmm.C+\(\"${NTUPDIR_5}/Zmumu/ntuples\",\"${OUTDIR}/FINAL_Zmm_5TeV_noRoch\",\"${S5}\",${LUMI5}\,0\)