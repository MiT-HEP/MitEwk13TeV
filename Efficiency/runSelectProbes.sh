#! /bin/bash

INPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-ntuples
OUTPUTDIR=/afs/cern.ch/work/c/cmedlock/public/wz-efficiency

#
# Select probes for muon efficiencies
#

root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuHLTEff\",0,1,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuL1Eff\",1,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuSelEff\",2,1,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuTrkEff\",3,1,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuStaEff\",4,1,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.root\",\"${OUTPUTDIR}/Zmm_MuStaEff_iso\",5,1,1\)

root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuHLTEff\",0,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuL1Eff\",1,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuSelEff\",2,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuTrkEff\",3,0\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff\",4,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZmm_MuStaEff_iso\",5,0\)

#
# Select probes for electron efficiencies
#
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleHLTEff\",0,1,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleL1Eff\",1,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleSelEff\",2,1,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleGsfEff\",3,1,1\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee_EleGsfSelEff\",4,1,1\)

#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleHLTEff\",0,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleL1Eff\",1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleSelEff\",2,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfEff\",3,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/DataZee_EleGsfSelEff\",4,0\)

rm *.so *.d
