#! /bin/bash

LUMI=2263
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

RooUnfoldDataZPt Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZPT}
RooUnfoldDataZPt Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSZPT}
RooUnfoldDataZPt_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSZPT}
RooUnfoldDataZPt_BinEffSys Zmm ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataZPt_ResScaleSys Zmm ${LUMI} ${ITERATIONSZPT}
RooUnfoldDataPhiStar Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_BinEffSys Zmm ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataPhiStar_ResScaleSys Zmm ${LUMI} ${ITERATIONSPHISTAR}
RooUnfoldDataZRapidity Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_BinEffSys Zmm ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataZRapidity_ResScaleSys Zmm ${LUMI} ${ITERATIONSZRAP}
RooUnfoldDataLep1Pt Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_BinEffSys Zmm ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep1Pt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP1PT}
RooUnfoldDataLep2Pt Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_BinEffSys Zmm ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLep2Pt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP2PT}
RooUnfoldDataLepNegPt Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_BinEffSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepNegPt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEPNEGPT}
RooUnfoldDataLepPosPt Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_BinEffSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLepPosPt_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEPPOSPT}
RooUnfoldDataLep1Eta Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_BinEffSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep1Eta_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP1ETA}
RooUnfoldDataLep2Eta Zmm ${LUMI} ${LUMIUNCERT} 0 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta Zmm ${LUMI} ${LUMIUNCERT} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta Zmm ${LUMI} ${LUMIUNCERT} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_UnfoldModelSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_UnfoldMatrixSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_EWKBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_EWKBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_TopBkgSys Zmm ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_TopBkgSys Zmm ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_StatEffSys Zmm ${LUMI} 1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_StatEffSys Zmm ${LUMI} -1 ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_BinEffSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_BkgShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_SigShapeEffSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}
RooUnfoldDataLep2Eta_ResScaleSys Zmm ${LUMI} ${ITERATIONSLEP2ETA}



