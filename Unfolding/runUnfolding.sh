#! /bin/bash

LUMI=2215
LUMIUNCERT=4.6

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
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/kfiekas/CMSSW_7_4_14/src/BootStrap/bin

RooUnfoldDataZPt Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldDataZPt Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_BinEffSys Zmumu ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_ResScaleSys Zmumu ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataPhiStar Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_BinEffSys Zmumu ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_ResScaleSys Zmumu ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataZRapidity Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_BinEffSys Zmumu ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_ResScaleSys Zmumu ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataLep1Pt Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_BinEffSys Zmumu ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_ResScaleSys Zmumu ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep2Pt Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_BinEffSys Zmumu ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_ResScaleSys Zmumu ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLepNegPt Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_BinEffSys Zmumu ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_ResScaleSys Zmumu ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepPosPt Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_BinEffSys Zmumu ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_ResScaleSys Zmumu ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLep1Eta Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_BinEffSys Zmumu ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_ResScaleSys Zmumu ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep2Eta Zmumu ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta Zmumu ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta Zmumu ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_UnfoldModelSys Zmumu ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_UnfoldMatrixSys Zmumu ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_EWKBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_EWKBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_TopBkgSys Zmumu ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_TopBkgSys Zmumu ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_StatEffSys Zmumu ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_StatEffSys Zmumu ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_BinEffSys Zmumu ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_BkgShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_SigShapeEffSys Zmumu ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_ResScaleSys Zmumu ${LUMI} ${ITERATIONSLEP2ETA}



