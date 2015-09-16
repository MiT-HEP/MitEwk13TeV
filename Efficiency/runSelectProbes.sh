#! /bin/bash

INPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-ntuples
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency

#
# Select probes for muon efficiencies
#

#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select_08_06.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff\",0,1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuL1Eff\",1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select_08_06.root\",\"${OUTPUTDIR}/Zmm_MuSelEff\",2,1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select_08_06.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff\",3,1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_trkEff_v2/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff_v2\",9,1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuStaEff\",4,1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select_08_06.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso\",5,1,1,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select_08_06.root\",\"${OUTPUTDIR}/Zmm_MuSelStaEff_iso\",8,1,1,0\)

#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select_08_04.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff\",0,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuL1Eff\",1,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select_08_04.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff\",2,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select_08_04.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff\",3,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_trkEff_v2/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff_v2\",9,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff\",4,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select_08_04.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso\",5,0,0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select_08_04.root\",\"${OUTPUTDIR}/DataZmm_MuSelStaEff_iso\",8,0,0,0\)

#
# Select probes for electron efficiencies
#
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",0,1,1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleL1Eff\",1,1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",2,1,1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleGsfEff\",3,1,1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",4,1,1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"../Selection/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleSCEff\",5,1,1,0\)

#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff\",0,0,0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleL1Eff\",1,0,0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff\",2,0,0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff\",3,0,0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff\",4,0,0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"../Selection/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSCEff\",5,0,0,0\)

# id efficiency for specific runs (to study excess in data of leptons in the positive endcap of the ECAL)
# mc
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleIDEff\",6,1,1,0\)
# data
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleIDEff_run251562\",6,0,0,251562\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleIDEff_run251883\",6,0,0,251883\)

# iso efficiency for specific runs (to study excess in data of leptons in the positive endcap of the ECAL)
# mc
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_factorize_id_iso/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleIsoEff\",7,1,1,0\)
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

rm *.so *.d
