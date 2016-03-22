#! /bin/bash
#./runofficial ZPt "p_{T}^{#mu^{+}#mu^{-}}     [GeV]"
#variable=$1
#xaxislabel=$2

workdir="/afs/cern.ch/user/x/xniu/WZXSection/Copy/CMSSW_7_6_3_patch2/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"

LUMI=2318.3
LUMIUNCERT=2.7

ITERATIONSZPT=4
ITERATIONSPHISTAR=3
ITERATIONSZRAP=3
ITERATIONSLEP1PT=5
ITERATIONSLEP2PT=5
ITERATIONSLEPNEGPT=5
ITERATIONSLEPPOSPT=5
ITERATIONSLEP1ETA=3
ITERATIONSLEP2ETA=3

PATH=${PATH}:./bin/
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/afs/cern.ch/user/x/xniu/WZXSection/Copy/CMSSW_7_6_3_patch2/src/BootStrap/bin

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt "p_{T}^{#mu^{+}#mu^{-}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZPt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap "|y^{#mu^{+}#mu^{-}}|" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm ZRap EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar "#phi_{#eta}*" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm PhiStar EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt "p_{T} (leading muon) [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Pt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt "p_{T} (2nd leading muon) [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Pt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta "|#eta| (leading muon)" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep1Eta EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta "p_{T} (leading muon) [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm Lep2Eta EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt "p_{T}^{#mu^{+}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepPosPt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${filedir}/MuUnfoldingInput/zmmph_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt "p_{T}^{#mu^{-}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/MuUnfoldingInput/zmm_UnfoldInputs.root ${workdir}/SignalExtraction/Zmm/Zmm_DataBkg.root Zmm LepNegPt EffBin ${LUMI} ${ITERATIONSZPT}
