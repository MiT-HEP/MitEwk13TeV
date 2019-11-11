#! /bin/bash

INPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/LowPU2017ID_13TeV #this one for ele at the moment
v=
OUTPUTDIR=/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/probes
#
# Select probes for muon efficiencies
#

# ## Zmm aMC@NLO with weights for the pythia and powheg w/ photos reweightings
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuHLTEff_tagPt\",0,0,1\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuSITEff_tagPt\",8,0,1\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_tagPt\",4,0,1\)

## Zmm Minlo sample for one of the uncertainties
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_minlo/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuHLTEff_minlo\",0,0,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_minlo/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuSITEff_minlo\",8,0,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_minlo/ntuples/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_minlo\",4,0,1\)

## Zmm Data
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuHLTEff_tagPt\",0,0,0\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuSITEff_tagPt\",8,0,0\)
# root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuStaEff_tagPt\",4,0,0\)

# These are the ones for low PU 13 TeV
# mc
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleHLTEff\",0,0,1,0\)
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleGSFSelEff\",4,0,1,0\)

# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleHLTEff_tagPt\",0,0,1,0\)
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleGSFSelEff_tagPt\",4,0,1,0\)

root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_minlo/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleHLTEff_minlo\",0,0,1,0\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_minlo/ntuples/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleGSFSelEff_minlo\",4,0,1,0\)
# # # # data
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleHLTEff\",0,0,0,0\)
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleGSFSelEff\",4,0,0,0\)

# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleHLTEff_tagPt\",0,0,0,0\)
# root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleGSFSelEff_tagPt\",4,0,0,0\)


# rm *.so *.d
