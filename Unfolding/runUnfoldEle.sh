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

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt "p_{T}^{#e^{+}#e^{-}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZPt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap "|y^{#e^{+}#e^{-}}|" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee ZRap EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar "#phi_{#eta}*" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee PhiStar EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt "p_{T} (leading ele) [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Pt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt "p_{T} (2nd leading ele) [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Pt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta "|#eta| (leading ele)" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep1Eta EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta "p_{T} (leading ele) [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee Lep2Eta EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt "p_{T}^{#e^{+}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepPosPt EffBin ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${filedir}/EleUnfoldingInput/zeeph_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt "p_{T}^{#e^{-}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/EleUnfoldingInput/zee_UnfoldInputs.root ${workdir}/SignalExtraction/Zee/Zee_DataBkg.root Zee LepNegPt EffBin ${LUMI} ${ITERATIONSZPT}
