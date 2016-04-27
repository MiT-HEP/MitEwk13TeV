#! /bin/bash

workdir="../SignalExtraction/Zee"
filedir="../UnfoldingInput/Zee"

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
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/afs/cern.ch/work/k/kfiekas/Analysis/WZCrossSection/CMSSW_7_4_14/src/BootStrap/bin

RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt "p_{T}^{e^{+}e^{-}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZPt EffBin ${LUMI} ${ITERATIONSZPT}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZRAP}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZRAP}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZRAP}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap "|y^{e^{+}e^{-}}|" ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap EWK 0.30 ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap EWK 0.30 ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap Top 0.10 ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap Top 0.10 ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap EffSigShape ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap EffBkgShape ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee ZRap EffBin ${LUMI} ${ITERATIONSZRAP}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSPHISTAR}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSPHISTAR}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar "#phi_{#eta}*" ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar EWK 0.30 ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar EWK 0.30 ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar Top 0.10 ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar Top 0.10 ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar EffSigShape ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar EffBkgShape ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee PhiStar EffBin ${LUMI} ${ITERATIONSPHISTAR}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1PT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1PT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt "p_{T} (leading electron) [GeV]" ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt Top 0.10 ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt Top 0.10 ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt EffSigShape ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt EffBkgShape ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Pt EffBin ${LUMI} ${ITERATIONSLEP1PT}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2PT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2PT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt "p_{T} (2nd leading electron) [GeV]" ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt Top 0.10 ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt Top 0.10 ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt EffSigShape ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt EffBkgShape ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Pt EffBin ${LUMI} ${ITERATIONSLEP2PT}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1ETA}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta "|#eta| (leading electron)" ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta Top 0.10 ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta Top 0.10 ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta EffSigShape ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta EffBkgShape ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep1Eta EffBin ${LUMI} ${ITERATIONSLEP1ETA}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2ETA}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta "p_{T} (leading electron) [GeV]" ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta Top 0.10 ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta Top 0.10 ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta EffSigShape ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta EffBkgShape ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee Lep2Eta EffBin ${LUMI} ${ITERATIONSLEP2ETA}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPPOSPT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt "p_{T}^{e^{+}} [GeV]" ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt Top 0.10 ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt Top 0.10 ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt EffSigShape ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt EffBkgShape ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepPosPt EffBin ${LUMI} ${ITERATIONSLEPPOSPT}


RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPNEGPT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_UnfoldModelSys ${filedir}/zee_UnfoldInputs.root ${filedir}/zeeph_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt "p_{T}^{e^{-}} [GeV]" ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt Top 0.10 ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt Top 0.10 ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_StatEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt EffSigShape ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt EffBkgShape ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_SystEffSys ${filedir}/zee_UnfoldInputs.root ${workdir}/Zee_DataBkg.root Zee LepNegPt EffBin ${LUMI} ${ITERATIONSLEPNEGPT}