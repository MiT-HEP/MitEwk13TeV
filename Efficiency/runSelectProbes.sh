#! /bin/bash

# INPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV
# INPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_matchDataTrig
# INPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2 # this one is ok, has prefire weight, ele scale is not right
INPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_wPrefire # this one is ok, has prefire weight, ele scale is not right
# INPUTDIR=/afs/cern.ch/work/a/arapyan/public/flat_fixed
# OUTPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_Efficiency_5TeV
# OUTPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_5TeV_Try2_Efficiency
OUTPUTDIR=/afs/cern.ch/user/s/sabrandt/work/public/LowPU_13TeV_Efficiency_v1/probes #_v1 made from ntuples which include prefire weight

#
# Select probes for muon efficiencies
#

# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuHLTEff_Test\",0,0,1\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm_MuL1Eff_Test\",1,0,1\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuSelEff_Test\",2,0,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuTrkEff_Test\",3,0,1\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_Test\",4,0,1\)

# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuHLTEff_Test\",0,0,0\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuL1Eff_Test\",1,0,0\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuSelEff_Test\",2,0,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuTrkEff_Test\",3,0,0\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuStaEff_Test\",4,0,0\)

# These are the ones for low PU 13 TeV
# mc
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_HLTEff_Test\",0,0,1,0\)
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_GSFSelEff_Test\",4,0,1,0\)
# # # # data
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_HLTEff_Test\",0,0,0,0\)
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_GSFSelEff_Test\",4,0,0,0\)

# iso efficiency for specific runs (to study excess in data of leptons in the positive endcap of the ECAL)
# mc
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleIsoEff\",7,1,1,0\)
# data
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleIsoEff_run251562\",7,0,0,251562\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleIsoEff_run251883\",7,0,0,251883\)

# trigger efficiency for specific runs (to study excess in data of leptons in the positive endcap of the ECAL)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251244\",0,0,0,251244\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251251\",0,0,0,251251\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251252\",0,0,0,251252\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251561\",0,0,0,251561\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251562\",0,0,0,251562\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251643\",0,0,0,251643\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251721\",0,0,0,251721\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff_run251883\",0,0,0,251883\)

# reconstruction efficiency for specific runs (to study excess in data of leptons in the positive endcap of the ECAL)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251244\",3,0,0,251244\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251251\",3,0,0,251251\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251252\",3,0,0,251252\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251561\",3,0,0,251561\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251562\",3,0,0,251562\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251643\",3,0,0,251643\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251721\",3,0,0,251721\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff_run251883\",3,0,0,251883\)

# id+iso efficiency for specific runs (to study excess in data of leptons in the positive endcap of the ECAL)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251244\",2,0,0,251244\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251251\",2,0,0,251251\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251252\",2,0,0,251252\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251561\",2,0,0,251561\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251562\",2,0,0,251562\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251643\",2,0,0,251643\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251721\",2,0,0,251721\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff_run251883\",2,0,0,251883\)

# rm *.so *.d
