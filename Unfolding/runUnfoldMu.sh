#! /bin/bash

workdir="../SignalExtraction/Zmm"
filedir="../UnfoldingInput/Zmm"

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

RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt "p_{T}^{#mu^{+}#mu^{-}} [GeV]" ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt EWK 0.30 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt EWK 0.30 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt Top 0.10 ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt Top 0.10 ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt EffSigShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt EffBkgShape ${LUMI} ${ITERATIONSZPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZPt EffBin ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_ResScaleSys Zmm ${LUMI} ${ITERATIONSZPT}

RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZRAP}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZRAP}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZRAP}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap "|y^{#mu^{+}#mu^{-}}|" ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap EWK 0.30 ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap EWK 0.30 ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap Top 0.10 ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap Top 0.10 ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap EffSigShape ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap EffBkgShape ${LUMI} ${ITERATIONSZRAP}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm ZRap EffBin ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_ResScaleSys Zmm ${LUMI} ${ITERATIONSZRAP}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSPHISTAR}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSPHISTAR}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar "#phi_{#eta}*" ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar EWK 0.30 ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar EWK 0.30 ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar Top 0.10 ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar Top 0.10 ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar EffSigShape ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar EffBkgShape ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm PhiStar EffBin ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_ResScaleSys Zmm ${LUMI} ${ITERATIONSPHISTAR}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1PT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1PT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt "p_{T} (leading muon) [GeV]" ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt Top 0.10 ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt Top 0.10 ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt EffSigShape ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt EffBkgShape ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Pt EffBin ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP1PT}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2PT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2PT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt "p_{T} (2nd leading muon) [GeV]" ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt Top 0.10 ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt Top 0.10 ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt EffSigShape ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt EffBkgShape ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Pt EffBin ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP2PT}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1ETA}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta "|#eta| (leading muon)" ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta Top 0.10 ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta Top 0.10 ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta EffSigShape ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta EffBkgShape ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep1Eta EffBin ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2ETA}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta "p_{T} (leading muon) [GeV]" ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta EWK 0.30 ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta EWK 0.30 ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta Top 0.10 ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta Top 0.10 ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta EffSigShape ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta EffBkgShape ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm Lep2Eta EffBin ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPPOSPT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt "p_{T}^{#mu^{+}} [GeV]" ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt Top 0.10 ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt Top 0.10 ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt EffSigShape ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt EffBkgShape ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepPosPt EffBin ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}


RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPNEGPT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_UnfoldModelSys ${filedir}/zmm_UnfoldInputs.root ${filedir}/zmmph_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt "p_{T}^{#mu^{-}} [GeV]" ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_UnfoldMatrixSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt EWK 0.30 ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt EWK 0.30 ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt Top 0.10 ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_BkgSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt Top 0.10 ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_StatEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt EffSigShape ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt EffBkgShape ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldData_SystEffSys ${filedir}/zmm_UnfoldInputs.root ${workdir}/Zmm_DataBkg.root Zmm LepNegPt EffBin ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
