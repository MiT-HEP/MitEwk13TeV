//================================================================================================
// 
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

//helper class to handle rochester corrections
#include <rochcor2015.h>
#include <muresolution_run2.h>

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

//=== MAIN MACRO ================================================================================================= 

void plotZmm(const TString  inputDir,    // input directory
	     const TString  outputDir,   // output directory
             const Double_t lumi,        // integrated luminosity (/fb)
	     const Bool_t   normToData=0 //draw MC normalized to data
) {
  gBenchmark->Start("plotZmm");
  gStyle->SetTitleOffset(1.100,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  const Double_t mu_MASS  = 0.1057;
  //
  // input ntuple file names
  //
  enum { eData, eZmm, eEWK, eTop };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  fnamev.push_back(inputDir + TString("/") + TString("data_select.root")); typev.push_back(eData);
  fnamev.push_back(inputDir + TString("/") + TString("zmm_select.raw.root"));   typev.push_back(eZmm);
  fnamev.push_back(inputDir + TString("/") + TString("ewk_select.raw.root"));  typev.push_back(eEWK);
  fnamev.push_back(inputDir + TString("/") + TString("top_select.raw.root"));  typev.push_back(eTop);

  //
  // Fit options
  //
  const Int_t    NBINS     = 60;
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;

  // efficiency files
  
  const TString dataHLTEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/MG/eff.root";
  const TString dataHLTEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/MG/eff.root";
  const TString zmmHLTEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/CT/eff.root";
  const TString zmmHLTEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/CT/eff.root";

  const TString dataSelEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/MG/eff.root";
  const TString dataSelEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/MG/eff.root";
  const TString zmmSelEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/CT/eff.root";
  const TString zmmSelEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/CT/eff.root";

  const TString dataTrkEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/MG/eff.root";
  const TString dataTrkEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/MG/eff.root";
  const TString zmmTrkEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/CT/eff.root";
  const TString zmmTrkEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/CT/eff.root";

  const TString dataStaEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/MG/eff.root";
  const TString dataStaEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/MG/eff.root";
  const TString zmmStaEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/CT/eff.root";
  const TString zmmStaEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/CT/eff.root";

  // efficiency files 2Bins

  const TString dataHLTEff2BinName_pos = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuHLTEff/MG/eff.root";
  const TString dataHLTEff2BinName_neg = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuHLTEff/MG/eff.root";
  const TString zmmHLTEff2BinName_pos  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuHLTEff/CT/eff.root";
  const TString zmmHLTEff2BinName_neg  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuHLTEff/CT/eff.root";

  const TString dataSelEff2BinName_pos = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/MG/eff.root";
  const TString dataSelEff2BinName_neg = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/MG/eff.root";
  const TString zmmSelEff2BinName_pos  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/CT/eff.root";
  const TString zmmSelEff2BinName_neg  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/CT/eff.root";

  const TString dataTrkEff2BinName_pos = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/MG/eff.root";
  const TString dataTrkEff2BinName_neg = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/MG/eff.root";
  const TString zmmTrkEff2BinName_pos  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/CT/eff.root";
  const TString zmmTrkEff2BinName_neg  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuSITEff/CT/eff.root";

  const TString dataStaEff2BinName_pos = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuStaEff/MG/eff.root";
  const TString dataStaEff2BinName_neg = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuStaEff/MG/eff.root";
  const TString zmmStaEff2BinName_pos  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuStaEff/CT/eff.root";
  const TString zmmStaEff2BinName_neg  = "/data/blue/xniu/WZXSection/NewMu/Uncertainty/MuStaEff/CT/eff.root";

  TString StaEffSignalShapeSys = "/data/blue/xniu/WZXSection/NewMu/MuStaSigSys.root";
  TString StaEffBackgroundShapeSys = "/data/blue/xniu/WZXSection/NewMu/MuStaBkgSys.root";
  TString SelEffSignalShapeSys = "/data/blue/xniu/WZXSection/NewMu/MuSITSigSys.root";
  TString SelEffBackgroundShapeSys = "/data/blue/xniu/WZXSection/NewMu/MuSITBkgSys.root";

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");


  // plot output file format
  const TString format("png");

  // setup efficiency shape systematics
  TFile *StaSigSysFile = new TFile(StaEffSignalShapeSys);
  TH2D *hStaSigSys = (TH2D*)StaSigSysFile->Get("h");
  TFile *StaBkgSysFile = new TFile(StaEffBackgroundShapeSys);
  TH2D *hStaBkgSys = (TH2D*)StaBkgSysFile->Get("h");
  TFile *SelSigSysFile = new TFile(SelEffSignalShapeSys);
  TH2D *hSelSigSys = (TH2D*)SelSigSysFile->Get("h");
  TFile *SelBkgSysFile = new TFile(SelEffBackgroundShapeSys);
  TH2D *hSelBkgSys = (TH2D*)SelBkgSysFile->Get("h");

  Int_t yield = 0;
  Double_t yield_zmm = 0, yield_zmm_unc=0;
  Double_t yield_ewk = 0, yield_ewk_unc=0;
   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // event category enumeration
  enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk }; // event category enum

 
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  // histograms for full selection
  double ZPtBins[35]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,375,500,1000};
  double PhiStarBins[28]={0,0.01,0.012,0.014,0.017,0.021,0.025,0.030,0.036,0.043,0.052,0.062,0.074,0.089,0.11,0.13,0.15,0.18,0.22,0.27,0.32,0.38,0.46,0.55,0.66,0.79,0.95,1.1};
  double Lep1PtBins[26]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,154,171,190,211,234,265,300};
  double Lep2PtBins[21]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,157,200};
  double LepNegPtBins[26]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,154,171,190,211,234,265,300};
  double LepPosPtBins[26]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,154,171,190,211,234,265,300};


  TH1D *hData = new TH1D("hData","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  TH1D *hZmm  = new TH1D("hZmm", "",NBINS,MASS_LOW,MASS_HIGH); hZmm->Sumw2();
  TH1D *hEWK  = new TH1D("hEWK", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hTop  = new TH1D("hTop", "",NBINS,MASS_LOW,MASS_HIGH); hTop->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();

  TH1D *hDataNPV = new TH1D("hDataNPV","",50,0,50); hDataNPV->Sumw2();
  TH1D *hZmmNPV  = new TH1D("hZmmNPV", "",50,0,50); hZmmNPV->Sumw2();
  TH1D *hEWKNPV  = new TH1D("hEWKNPV", "",50,0,50); hEWKNPV->Sumw2();
  TH1D *hTopNPV  = new TH1D("hTopNPV", "",50,0,50); hTopNPV->Sumw2();
  TH1D *hMCNPV   = new TH1D("hMCNPV",  "",50,0,50); hMCNPV->Sumw2();

  TH1D *hDataZPt = new TH1D("hDataZPt","",34,ZPtBins); hDataZPt->Sumw2();
  TH1D *hZmmZPt  = new TH1D("hZmmZPt", "",34,ZPtBins); hZmmZPt->Sumw2();
  TH1D *hEWKZPt  = new TH1D("hEWKZPt", "",34,ZPtBins); hEWKZPt->Sumw2();
  TH1D *hTopZPt  = new TH1D("hTopZPt", "",34,ZPtBins); hTopZPt->Sumw2();
  TH1D *hMCZPt   = new TH1D("hMCZPt",  "",34,ZPtBins); hMCZPt->Sumw2();

  TH1D *hEWKZPt_EffBin  = new TH1D("hEWKZPt_EffBin", "",34,ZPtBins); hEWKZPt_EffBin->Sumw2();
  TH1D *hTopZPt_EffBin  = new TH1D("hTopZPt_EffBin", "",34,ZPtBins); hTopZPt_EffBin->Sumw2();
  TH1D *hEWKZPt_EffStatUp  = new TH1D("hEWKZPt_EffStatUp", "",34,ZPtBins); hEWKZPt_EffStatUp->Sumw2();
  TH1D *hTopZPt_EffStatUp  = new TH1D("hTopZPt_EffStatUp", "",34,ZPtBins); hTopZPt_EffStatUp->Sumw2();
  TH1D *hEWKZPt_EffStatDown  = new TH1D("hEWKZPt_EffStatDown", "",34,ZPtBins); hEWKZPt_EffStatDown->Sumw2();
  TH1D *hTopZPt_EffStatDown  = new TH1D("hTopZPt_EffStatDown", "",34,ZPtBins); hTopZPt_EffStatDown->Sumw2();
  TH1D *hEWKZPt_EffSigShape  = new TH1D("hEWKZPt_EffSigShape", "",34,ZPtBins); hEWKZPt_EffSigShape->Sumw2();
  TH1D *hTopZPt_EffSigShape  = new TH1D("hTopZPt_EffSigShape", "",34,ZPtBins); hTopZPt_EffSigShape->Sumw2();
  TH1D *hEWKZPt_EffBkgShape  = new TH1D("hEWKZPt_EffBkgShape", "",34,ZPtBins); hEWKZPt_EffBkgShape->Sumw2();
  TH1D *hTopZPt_EffBkgShape  = new TH1D("hTopZPt_EffBkgShape", "",34,ZPtBins); hTopZPt_EffBkgShape->Sumw2();

  TH1D *hDataPhiStar = new TH1D("hDataPhiStar","",27,PhiStarBins); hDataPhiStar->Sumw2();
  TH1D *hZmmPhiStar  = new TH1D("hZmmPhiStar", "",27,PhiStarBins); hZmmPhiStar->Sumw2();
  TH1D *hEWKPhiStar  = new TH1D("hEWKPhiStar", "",27,PhiStarBins); hEWKPhiStar->Sumw2();
  TH1D *hTopPhiStar  = new TH1D("hTopPhiStar", "",27,PhiStarBins); hTopPhiStar->Sumw2();
  TH1D *hMCPhiStar   = new TH1D("hMCPhiStar",  "",27,PhiStarBins); hMCPhiStar->Sumw2();
  TH1D *hEWKPhiStar_EffBin  = new TH1D("hEWKPhiStar_EffBin", "",27,PhiStarBins); hEWKPhiStar_EffBin->Sumw2();
  TH1D *hTopPhiStar_EffBin  = new TH1D("hTopPhiStar_EffBin", "",27,PhiStarBins); hTopPhiStar_EffBin->Sumw2();
  TH1D *hEWKPhiStar_EffStatUp  = new TH1D("hEWKPhiStar_EffStatUp", "",27,PhiStarBins); hEWKPhiStar_EffStatUp->Sumw2();
  TH1D *hTopPhiStar_EffStatUp  = new TH1D("hTopPhiStar_EffStatUp", "",27,PhiStarBins); hTopPhiStar_EffStatUp->Sumw2();
  TH1D *hEWKPhiStar_EffStatDown  = new TH1D("hEWKPhiStar_EffStatDown", "",27,PhiStarBins); hEWKPhiStar_EffStatDown->Sumw2();
  TH1D *hTopPhiStar_EffStatDown  = new TH1D("hTopPhiStar_EffStatDown", "",27,PhiStarBins); hTopPhiStar_EffStatDown->Sumw2();
  TH1D *hEWKPhiStar_EffSigShape  = new TH1D("hEWKPhiStar_EffSigShape", "",27,PhiStarBins); hEWKPhiStar_EffSigShape->Sumw2();
  TH1D *hTopPhiStar_EffSigShape  = new TH1D("hTopPhiStar_EffSigShape", "",27,PhiStarBins); hTopPhiStar_EffSigShape->Sumw2();
  TH1D *hEWKPhiStar_EffBkgShape  = new TH1D("hEWKPhiStar_EffBkgShape", "",27,PhiStarBins); hEWKPhiStar_EffBkgShape->Sumw2();
  TH1D *hTopPhiStar_EffBkgShape  = new TH1D("hTopPhiStar_EffBkgShape", "",27,PhiStarBins); hTopPhiStar_EffBkgShape->Sumw2();

  
  TH1D *hDataZRap = new TH1D("hDataZRap","",24,0,2.4); hDataZRap->Sumw2();
  TH1D *hZmmZRap  = new TH1D("hZmmZRap", "",24,0,2.4); hZmmZRap->Sumw2();
  TH1D *hEWKZRap  = new TH1D("hEWKZRap", "",24,0,2.4); hEWKZRap->Sumw2();
  TH1D *hTopZRap  = new TH1D("hTopZRap", "",24,0,2.4); hTopZRap->Sumw2();
  TH1D *hMCZRap   = new TH1D("hMCZRap",  "",24,0,2.4); hMCZRap->Sumw2();

  TH1D *hEWKZRap_EffBin  = new TH1D("hEWKZRap_EffBin", "",24,0,2.4); hEWKZRap_EffBin->Sumw2();
  TH1D *hTopZRap_EffBin  = new TH1D("hTopZRap_EffBin", "",24,0,2.4); hTopZRap_EffBin->Sumw2();
  TH1D *hEWKZRap_EffStatUp  = new TH1D("hEWKZRap_EffStatUp", "",24,0,2.4); hEWKZRap_EffStatUp->Sumw2();
  TH1D *hTopZRap_EffStatUp  = new TH1D("hTopZRap_EffStatUp", "",24,0,2.4); hTopZRap_EffStatUp->Sumw2();
  TH1D *hEWKZRap_EffStatDown  = new TH1D("hEWKZRap_EffStatDown", "",24,0,2.4); hEWKZRap_EffStatDown->Sumw2();
  TH1D *hTopZRap_EffStatDown  = new TH1D("hTopZRap_EffStatDown", "",24,0,2.4); hTopZRap_EffStatDown->Sumw2();
  TH1D *hEWKZRap_EffSigShape  = new TH1D("hEWKZRap_EffSigShape", "",24,0,2.4); hEWKZRap_EffSigShape->Sumw2();
  TH1D *hTopZRap_EffSigShape  = new TH1D("hTopZRap_EffSigShape", "",24,0,2.4); hTopZRap_EffSigShape->Sumw2();
  TH1D *hEWKZRap_EffBkgShape  = new TH1D("hEWKZRap_EffBkgShape", "",24,0,2.4); hEWKZRap_EffBkgShape->Sumw2();
  TH1D *hTopZRap_EffBkgShape  = new TH1D("hTopZRap_EffBkgShape", "",24,0,2.4); hTopZRap_EffBkgShape->Sumw2();


  TH1D *hDataLep1Pt = new TH1D("hDataLep1Pt","",25,Lep1PtBins); hDataLep1Pt->Sumw2();
  TH1D *hZmmLep1Pt  = new TH1D("hZmmLep1Pt", "",25,Lep1PtBins); hZmmLep1Pt->Sumw2();
  TH1D *hEWKLep1Pt  = new TH1D("hEWKLep1Pt", "",25,Lep1PtBins); hEWKLep1Pt->Sumw2();
  TH1D *hTopLep1Pt  = new TH1D("hTopLep1Pt", "",25,Lep1PtBins); hTopLep1Pt->Sumw2();
  TH1D *hMCLep1Pt   = new TH1D("hMCLep1Pt",  "",25,Lep1PtBins); hMCLep1Pt->Sumw2();

  TH1D *hEWKLep1Pt_EffBin  = new TH1D("hEWKLep1Pt_EffBin", "",25,Lep1PtBins); hEWKLep1Pt_EffBin->Sumw2();
  TH1D *hTopLep1Pt_EffBin  = new TH1D("hTopLep1Pt_EffBin", "",25,Lep1PtBins); hTopLep1Pt_EffBin->Sumw2();
  TH1D *hEWKLep1Pt_EffStatUp  = new TH1D("hEWKLep1Pt_EffStatUp", "",25,Lep1PtBins); hEWKLep1Pt_EffStatUp->Sumw2();
  TH1D *hTopLep1Pt_EffStatUp  = new TH1D("hTopLep1Pt_EffStatUp", "",25,Lep1PtBins); hTopLep1Pt_EffStatUp->Sumw2();
  TH1D *hEWKLep1Pt_EffStatDown  = new TH1D("hEWKLep1Pt_EffStatDown", "",25,Lep1PtBins); hEWKLep1Pt_EffStatDown->Sumw2();
  TH1D *hTopLep1Pt_EffStatDown  = new TH1D("hTopLep1Pt_EffStatDown", "",25,Lep1PtBins); hTopLep1Pt_EffStatDown->Sumw2();
  TH1D *hEWKLep1Pt_EffSigShape  = new TH1D("hEWKLep1Pt_EffSigShape", "",25,Lep1PtBins); hEWKLep1Pt_EffSigShape->Sumw2();
  TH1D *hTopLep1Pt_EffSigShape  = new TH1D("hTopLep1Pt_EffSigShape", "",25,Lep1PtBins); hTopLep1Pt_EffSigShape->Sumw2();
  TH1D *hEWKLep1Pt_EffBkgShape  = new TH1D("hEWKLep1Pt_EffBkgShape", "",25,Lep1PtBins); hEWKLep1Pt_EffBkgShape->Sumw2();
  TH1D *hTopLep1Pt_EffBkgShape  = new TH1D("hTopLep1Pt_EffBkgShape", "",25,Lep1PtBins); hTopLep1Pt_EffBkgShape->Sumw2();


  TH1D *hDataLep2Pt = new TH1D("hDataLep2Pt","",20,Lep2PtBins); hDataLep2Pt->Sumw2();
  TH1D *hZmmLep2Pt  = new TH1D("hZmmLep2Pt", "",20,Lep2PtBins); hZmmLep2Pt->Sumw2();
  TH1D *hEWKLep2Pt  = new TH1D("hEWKLep2Pt", "",20,Lep2PtBins); hEWKLep2Pt->Sumw2();
  TH1D *hTopLep2Pt  = new TH1D("hTopLep2Pt", "",20,Lep2PtBins); hTopLep2Pt->Sumw2();
  TH1D *hMCLep2Pt   = new TH1D("hMCLep2Pt",  "",20,Lep2PtBins); hMCLep2Pt->Sumw2();

  TH1D *hEWKLep2Pt_EffBin  = new TH1D("hEWKLep2Pt_EffBin", "",20,Lep2PtBins); hEWKLep2Pt_EffBin->Sumw2();
  TH1D *hTopLep2Pt_EffBin  = new TH1D("hTopLep2Pt_EffBin", "",20,Lep2PtBins); hTopLep2Pt_EffBin->Sumw2();
  TH1D *hEWKLep2Pt_EffStatUp  = new TH1D("hEWKLep2Pt_EffStatUp", "",20,Lep2PtBins); hEWKLep2Pt_EffStatUp->Sumw2();
  TH1D *hTopLep2Pt_EffStatUp  = new TH1D("hTopLep2Pt_EffStatUp", "",20,Lep2PtBins); hTopLep2Pt_EffStatUp->Sumw2();
  TH1D *hEWKLep2Pt_EffStatDown  = new TH1D("hEWKLep2Pt_EffStatDown", "",20,Lep2PtBins); hEWKLep2Pt_EffStatDown->Sumw2();
  TH1D *hTopLep2Pt_EffStatDown  = new TH1D("hTopLep2Pt_EffStatDown", "",20,Lep2PtBins); hTopLep2Pt_EffStatDown->Sumw2();
  TH1D *hEWKLep2Pt_EffSigShape  = new TH1D("hEWKLep2Pt_EffSigShape", "",20,Lep2PtBins); hEWKLep2Pt_EffSigShape->Sumw2();
  TH1D *hTopLep2Pt_EffSigShape  = new TH1D("hTopLep2Pt_EffSigShape", "",20,Lep2PtBins); hTopLep2Pt_EffSigShape->Sumw2();
  TH1D *hEWKLep2Pt_EffBkgShape  = new TH1D("hEWKLep2Pt_EffBkgShape", "",20,Lep2PtBins); hEWKLep2Pt_EffBkgShape->Sumw2();
  TH1D *hTopLep2Pt_EffBkgShape  = new TH1D("hTopLep2Pt_EffBkgShape", "",20,Lep2PtBins); hTopLep2Pt_EffBkgShape->Sumw2();

  TH1D *hDataLepNegPt = new TH1D("hDataLepNegPt","",25,LepNegPtBins); hDataLepNegPt->Sumw2();
  TH1D *hZmmLepNegPt  = new TH1D("hZmmLepNegPt", "",25,LepNegPtBins); hZmmLepNegPt->Sumw2();
  TH1D *hEWKLepNegPt  = new TH1D("hEWKLepNegPt", "",25,LepNegPtBins); hEWKLepNegPt->Sumw2();
  TH1D *hTopLepNegPt  = new TH1D("hTopLepNegPt", "",25,LepNegPtBins); hTopLepNegPt->Sumw2();
  TH1D *hMCLepNegPt   = new TH1D("hMCLepNegPt",  "",25,LepNegPtBins); hMCLepNegPt->Sumw2();

  TH1D *hEWKLepNegPt_EffBin  = new TH1D("hEWKLepNegPt_EffBin", "",25,LepNegPtBins); hEWKLepNegPt_EffBin->Sumw2();
  TH1D *hTopLepNegPt_EffBin  = new TH1D("hTopLepNegPt_EffBin", "",25,LepNegPtBins); hTopLepNegPt_EffBin->Sumw2();
  TH1D *hEWKLepNegPt_EffStatUp  = new TH1D("hEWKLepNegPt_EffStatUp", "",25,LepNegPtBins); hEWKLepNegPt_EffStatUp->Sumw2();
  TH1D *hTopLepNegPt_EffStatUp  = new TH1D("hTopLepNegPt_EffStatUp", "",25,LepNegPtBins); hTopLepNegPt_EffStatUp->Sumw2();
  TH1D *hEWKLepNegPt_EffStatDown  = new TH1D("hEWKLepNegPt_EffStatDown", "",25,LepNegPtBins); hEWKLepNegPt_EffStatDown->Sumw2();
  TH1D *hTopLepNegPt_EffStatDown  = new TH1D("hTopLepNegPt_EffStatDown", "",25,LepNegPtBins); hTopLepNegPt_EffStatDown->Sumw2();
  TH1D *hEWKLepNegPt_EffSigShape  = new TH1D("hEWKLepNegPt_EffSigShape", "",25,LepNegPtBins); hEWKLepNegPt_EffSigShape->Sumw2();
  TH1D *hTopLepNegPt_EffSigShape  = new TH1D("hTopLepNegPt_EffSigShape", "",25,LepNegPtBins); hTopLepNegPt_EffSigShape->Sumw2();
  TH1D *hEWKLepNegPt_EffBkgShape  = new TH1D("hEWKLepNegPt_EffBkgShape", "",25,LepNegPtBins); hEWKLepNegPt_EffBkgShape->Sumw2();
  TH1D *hTopLepNegPt_EffBkgShape  = new TH1D("hTopLepNegPt_EffBkgShape", "",25,LepNegPtBins); hTopLepNegPt_EffBkgShape->Sumw2();


  TH1D *hDataLepPosPt = new TH1D("hDataLepPosPt","",25,LepPosPtBins); hDataLepPosPt->Sumw2();
  TH1D *hZmmLepPosPt  = new TH1D("hZmmLepPosPt", "",25,LepPosPtBins); hZmmLepPosPt->Sumw2();
  TH1D *hEWKLepPosPt  = new TH1D("hEWKLepPosPt", "",25,LepPosPtBins); hEWKLepPosPt->Sumw2();
  TH1D *hTopLepPosPt  = new TH1D("hTopLepPosPt", "",25,LepPosPtBins); hTopLepPosPt->Sumw2();
  TH1D *hMCLepPosPt   = new TH1D("hMCLepPosPt",  "",25,LepPosPtBins); hMCLepPosPt->Sumw2();

  TH1D *hEWKLepPosPt_EffBin  = new TH1D("hEWKLepPosPt_EffBin", "",25,LepPosPtBins); hEWKLepPosPt_EffBin->Sumw2();
  TH1D *hTopLepPosPt_EffBin  = new TH1D("hTopLepPosPt_EffBin", "",25,LepPosPtBins); hTopLepPosPt_EffBin->Sumw2();
  TH1D *hEWKLepPosPt_EffStatUp  = new TH1D("hEWKLepPosPt_EffStatUp", "",25,LepPosPtBins); hEWKLepPosPt_EffStatUp->Sumw2();
  TH1D *hTopLepPosPt_EffStatUp  = new TH1D("hTopLepPosPt_EffStatUp", "",25,LepPosPtBins); hTopLepPosPt_EffStatUp->Sumw2();
  TH1D *hEWKLepPosPt_EffStatDown  = new TH1D("hEWKLepPosPt_EffStatDown", "",25,LepPosPtBins); hEWKLepPosPt_EffStatDown->Sumw2();
  TH1D *hTopLepPosPt_EffStatDown  = new TH1D("hTopLepPosPt_EffStatDown", "",25,LepPosPtBins); hTopLepPosPt_EffStatDown->Sumw2();
  TH1D *hEWKLepPosPt_EffSigShape  = new TH1D("hEWKLepPosPt_EffSigShape", "",25,LepPosPtBins); hEWKLepPosPt_EffSigShape->Sumw2();
  TH1D *hTopLepPosPt_EffSigShape  = new TH1D("hTopLepPosPt_EffSigShape", "",25,LepPosPtBins); hTopLepPosPt_EffSigShape->Sumw2();
  TH1D *hEWKLepPosPt_EffBkgShape  = new TH1D("hEWKLepPosPt_EffBkgShape", "",25,LepPosPtBins); hEWKLepPosPt_EffBkgShape->Sumw2();
  TH1D *hTopLepPosPt_EffBkgShape  = new TH1D("hTopLepPosPt_EffBkgShape", "",25,LepPosPtBins); hTopLepPosPt_EffBkgShape->Sumw2();



  TH1D *hDataLep1Eta = new TH1D("hDataLep1Eta","",24,0,2.4); hDataLep1Eta->Sumw2();
  TH1D *hZmmLep1Eta  = new TH1D("hZmmLep1Eta", "",24,0,2.4); hZmmLep1Eta->Sumw2();
  TH1D *hEWKLep1Eta  = new TH1D("hEWKLep1Eta", "",24,0,2.4); hEWKLep1Eta->Sumw2();
  TH1D *hTopLep1Eta  = new TH1D("hTopLep1Eta", "",24,0,2.4); hTopLep1Eta->Sumw2();
  TH1D *hMCLep1Eta   = new TH1D("hMCLep1Eta",  "",24,0,2.4); hMCLep1Eta->Sumw2();

  TH1D *hEWKLep1Eta_EffBin  = new TH1D("hEWKLep1Eta_EffBin", "",24,0,2.4); hEWKLep1Eta_EffBin->Sumw2();
  TH1D *hTopLep1Eta_EffBin  = new TH1D("hTopLep1Eta_EffBin", "",24,0,2.4); hTopLep1Eta_EffBin->Sumw2();
  TH1D *hEWKLep1Eta_EffStatUp  = new TH1D("hEWKLep1Eta_EffStatUp", "",24,0,2.4); hEWKLep1Eta_EffStatUp->Sumw2();
  TH1D *hTopLep1Eta_EffStatUp  = new TH1D("hTopLep1Eta_EffStatUp", "",24,0,2.4); hTopLep1Eta_EffStatUp->Sumw2();
  TH1D *hEWKLep1Eta_EffStatDown  = new TH1D("hEWKLep1Eta_EffStatDown", "",24,0,2.4); hEWKLep1Eta_EffStatDown->Sumw2();
  TH1D *hTopLep1Eta_EffStatDown  = new TH1D("hTopLep1Eta_EffStatDown", "",24,0,2.4); hTopLep1Eta_EffStatDown->Sumw2();
  TH1D *hEWKLep1Eta_EffSigShape  = new TH1D("hEWKLep1Eta_EffSigShape", "",24,0,2.4); hEWKLep1Eta_EffSigShape->Sumw2();
  TH1D *hTopLep1Eta_EffSigShape  = new TH1D("hTopLep1Eta_EffSigShape", "",24,0,2.4); hTopLep1Eta_EffSigShape->Sumw2();
  TH1D *hEWKLep1Eta_EffBkgShape  = new TH1D("hEWKLep1Eta_EffBkgShape", "",24,0,2.4); hEWKLep1Eta_EffBkgShape->Sumw2();
  TH1D *hTopLep1Eta_EffBkgShape  = new TH1D("hTopLep1Eta_EffBkgShape", "",24,0,2.4); hTopLep1Eta_EffBkgShape->Sumw2();


  TH1D *hDataLep2Eta = new TH1D("hDataLep2Eta","",24,0,2.4); hDataLep2Eta->Sumw2();
  TH1D *hZmmLep2Eta  = new TH1D("hZmmLep2Eta", "",24,0,2.4); hZmmLep2Eta->Sumw2();
  TH1D *hEWKLep2Eta  = new TH1D("hEWKLep2Eta", "",24,0,2.4); hEWKLep2Eta->Sumw2();
  TH1D *hTopLep2Eta  = new TH1D("hTopLep2Eta", "",24,0,2.4); hTopLep2Eta->Sumw2();
  TH1D *hMCLep2Eta   = new TH1D("hMCLep2Eta",  "",24,0,2.4); hMCLep2Eta->Sumw2();

  TH1D *hEWKLep2Eta_EffBin  = new TH1D("hEWKLep2Eta_EffBin", "",24,0,2.4); hEWKLep2Eta_EffBin->Sumw2();
  TH1D *hTopLep2Eta_EffBin  = new TH1D("hTopLep2Eta_EffBin", "",24,0,2.4); hTopLep2Eta_EffBin->Sumw2();
  TH1D *hEWKLep2Eta_EffStatUp  = new TH1D("hEWKLep2Eta_EffStatUp", "",24,0,2.4); hEWKLep2Eta_EffStatUp->Sumw2();
  TH1D *hTopLep2Eta_EffStatUp  = new TH1D("hTopLep2Eta_EffStatUp", "",24,0,2.4); hTopLep2Eta_EffStatUp->Sumw2();
  TH1D *hEWKLep2Eta_EffStatDown  = new TH1D("hEWKLep2Eta_EffStatDown", "",24,0,2.4); hEWKLep2Eta_EffStatDown->Sumw2();
  TH1D *hTopLep2Eta_EffStatDown  = new TH1D("hTopLep2Eta_EffStatDown", "",24,0,2.4); hTopLep2Eta_EffStatDown->Sumw2();
  TH1D *hEWKLep2Eta_EffSigShape  = new TH1D("hEWKLep2Eta_EffSigShape", "",24,0,2.4); hEWKLep2Eta_EffSigShape->Sumw2();
  TH1D *hTopLep2Eta_EffSigShape  = new TH1D("hTopLep2Eta_EffSigShape", "",24,0,2.4); hTopLep2Eta_EffSigShape->Sumw2();
  TH1D *hEWKLep2Eta_EffBkgShape  = new TH1D("hEWKLep2Eta_EffBkgShape", "",24,0,2.4); hEWKLep2Eta_EffBkgShape->Sumw2();
  TH1D *hTopLep2Eta_EffBkgShape  = new TH1D("hTopLep2Eta_EffBkgShape", "",24,0,2.4); hTopLep2Eta_EffBkgShape->Sumw2();

    
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Int_t   q1, q2;
  TLorentzVector *lep1=0, *lep2=0;

  //
  // HLT efficiency
  //
  cout << "Loading trigger efficiencies..." << endl;
  
  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
    
  TFile *zmmHLTEffFile_pos = new TFile(zmmHLTEffName_pos);
  CEffUser2D zmmHLTEff_pos;
  zmmHLTEff_pos.loadEff((TH2D*)zmmHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEffFile_neg = new TFile(zmmHLTEffName_neg);
  CEffUser2D zmmHLTEff_neg;
  zmmHLTEff_neg.loadEff((TH2D*)zmmHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataHLTEff2BinFile_pos = new TFile(dataHLTEff2BinName_pos);
  CEffUser2D dataHLTEff2Bin_pos;
  dataHLTEff2Bin_pos.loadEff((TH2D*)dataHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEff2BinFile_neg = new TFile(dataHLTEff2BinName_neg);
  CEffUser2D dataHLTEff2Bin_neg;
  dataHLTEff2Bin_neg.loadEff((TH2D*)dataHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrhEtaPt"));
    
  TFile *zmmHLTEff2BinFile_pos = new TFile(zmmHLTEff2BinName_pos);
  CEffUser2D zmmHLTEff2Bin_pos;
  zmmHLTEff2Bin_pos.loadEff((TH2D*)zmmHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEff2BinFile_neg = new TFile(zmmHLTEff2BinName_neg);
  CEffUser2D zmmHLTEff2Bin_neg;
  zmmHLTEff2Bin_neg.loadEff((TH2D*)zmmHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEff2BinFile_neg->Get("hErrhEtaPt"));

 
  
  //
  // Selection efficiency
  //
  cout << "Loading selection efficiencies..." << endl;
  
  TFile *dataSelEffFile_pos = new TFile(dataSelEffName_pos);
  CEffUser2D dataSelEff_pos;
  dataSelEff_pos.loadEff((TH2D*)dataSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEffFile_neg = new TFile(dataSelEffName_neg);
  CEffUser2D dataSelEff_neg;
  dataSelEff_neg.loadEff((TH2D*)dataSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEffFile_pos = new TFile(zmmSelEffName_pos);
  CEffUser2D zmmSelEff_pos;
  zmmSelEff_pos.loadEff((TH2D*)zmmSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEffFile_neg = new TFile(zmmSelEffName_neg);
  CEffUser2D zmmSelEff_neg;
  zmmSelEff_neg.loadEff((TH2D*)zmmSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataSelEff2BinFile_pos = new TFile(dataSelEff2BinName_pos);
  CEffUser2D dataSelEff2Bin_pos;
  dataSelEff2Bin_pos.loadEff((TH2D*)dataSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEff2BinFile_neg = new TFile(dataSelEff2BinName_neg);
  CEffUser2D dataSelEff2Bin_neg;
  dataSelEff2Bin_neg.loadEff((TH2D*)dataSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEff2BinFile_pos = new TFile(zmmSelEff2BinName_pos);
  CEffUser2D zmmSelEff2Bin_pos;
  zmmSelEff2Bin_pos.loadEff((TH2D*)zmmSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEff2BinFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEff2BinFile_neg = new TFile(zmmSelEff2BinName_neg);
  CEffUser2D zmmSelEff2Bin_neg;
  zmmSelEff2Bin_neg.loadEff((TH2D*)zmmSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEff2BinFile_neg->Get("hErrhEtaPt"));

 
  //
  // Standalone efficiency
  //
  cout << "Loading standalone efficiencies..." << endl;
  
  TFile *dataStaEffFile_pos = new TFile(dataStaEffName_pos);
  CEffUser2D dataStaEff_pos;
  dataStaEff_pos.loadEff((TH2D*)dataStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEffFile_neg = new TFile(dataStaEffName_neg);
  CEffUser2D dataStaEff_neg;
  dataStaEff_neg.loadEff((TH2D*)dataStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_pos = new TFile(zmmStaEffName_pos);
  CEffUser2D zmmStaEff_pos;
  zmmStaEff_pos.loadEff((TH2D*)zmmStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_neg = new TFile(zmmStaEffName_neg);
  CEffUser2D zmmStaEff_neg;
  zmmStaEff_neg.loadEff((TH2D*)zmmStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataStaEff2BinFile_pos = new TFile(dataStaEff2BinName_pos);
  CEffUser2D dataStaEff2Bin_pos;
  dataStaEff2Bin_pos.loadEff((TH2D*)dataStaEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEff2BinFile_neg = new TFile(dataStaEff2BinName_neg);
  CEffUser2D dataStaEff2Bin_neg;
  dataStaEff2Bin_neg.loadEff((TH2D*)dataStaEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEff2BinFile_pos = new TFile(zmmStaEff2BinName_pos);
  CEffUser2D zmmStaEff2Bin_pos;
  zmmStaEff2Bin_pos.loadEff((TH2D*)zmmStaEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEff2BinFile_neg = new TFile(zmmStaEff2BinName_neg);
  CEffUser2D zmmStaEff2Bin_neg;
  zmmStaEff2Bin_neg.loadEff((TH2D*)zmmStaEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEff2BinFile_neg->Get("hErrhEtaPt"));
 
  //
  // Tracker efficiency
  //
  cout << "Loading track efficiencies..." << endl;
  
  TFile *dataTrkEffFile_pos = new TFile(dataTrkEffName_pos);
  CEffUser2D dataTrkEff_pos;
  dataTrkEff_pos.loadEff((TH2D*)dataTrkEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataTrkEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataTrkEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataTrkEffFile_neg = new TFile(dataTrkEffName_neg);
  CEffUser2D dataTrkEff_neg;
  dataTrkEff_neg.loadEff((TH2D*)dataTrkEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataTrkEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataTrkEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEffFile_pos = new TFile(zmmTrkEffName_pos);
  CEffUser2D zmmTrkEff_pos;
  zmmTrkEff_pos.loadEff((TH2D*)zmmTrkEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmTrkEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmTrkEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEffFile_neg = new TFile(zmmTrkEffName_neg);
  CEffUser2D zmmTrkEff_neg;
  zmmTrkEff_neg.loadEff((TH2D*)zmmTrkEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmTrkEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmTrkEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataTrkEff2BinFile_pos = new TFile(dataTrkEff2BinName_pos);
  CEffUser2D dataTrkEff2Bin_pos;
  dataTrkEff2Bin_pos.loadEff((TH2D*)dataTrkEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataTrkEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataTrkEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataTrkEff2BinFile_neg = new TFile(dataTrkEff2BinName_neg);
  CEffUser2D dataTrkEff2Bin_neg;
  dataTrkEff2Bin_neg.loadEff((TH2D*)dataTrkEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataTrkEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataTrkEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEff2BinFile_pos = new TFile(zmmTrkEff2BinName_pos);
  CEffUser2D zmmTrkEff2Bin_pos;
  zmmTrkEff2Bin_pos.loadEff((TH2D*)zmmTrkEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmTrkEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmTrkEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmTrkEff2BinFile_neg = new TFile(zmmTrkEff2BinName_neg);
  CEffUser2D zmmTrkEff2Bin_neg;
  zmmTrkEff2Bin_neg.loadEff((TH2D*)zmmTrkEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmTrkEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmTrkEff2BinFile_neg->Get("hErrhEtaPt"));

  //Setting up rochester corrections
  rochcor2015 *rmcor = new rochcor2015();

  TFile *infile=0;
  TTree *intree=0;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree -> SetBranchStatus("*",0);
    intree -> SetBranchStatus("runNum",1);
    intree -> SetBranchStatus("lumiSec",1);
    intree -> SetBranchStatus("evtNum",1);
    intree -> SetBranchStatus("category",1);
    intree -> SetBranchStatus("npv",1);
    intree -> SetBranchStatus("npu",1);
    intree -> SetBranchStatus("scale1fb",1);
    intree -> SetBranchStatus("scale1fbUp",1);
    intree -> SetBranchStatus("scale1fbDown",1);
    intree -> SetBranchStatus("q1",1);
    intree -> SetBranchStatus("q2",1);
    intree -> SetBranchStatus("lep1",1);
    intree -> SetBranchStatus("lep2",1);

    intree->SetBranchAddress("runNum",     &runNum);      // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);      // event number
    intree->SetBranchAddress("category",   &category);    // dilepton category
    intree->SetBranchAddress("npv",        &npv);	  // number of primary vertices
    intree->SetBranchAddress("npu",        &npu);	  // number of in-time PU events (MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",   &scale1fbDown);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("q1",         &q1);	  // charge of tag lepton
    intree->SetBranchAddress("q2",         &q2);	  // charge of probe lepton
    intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
   
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;
      if(q1*q2>0) continue;
      
      float mass = 0;
      float pt = 0;
      float rapidity = 0;
      float phiacop=0;
      float costhetastar=0;
      float phistar=0;
     
      Double_t weight=1;
      if(typev[ifile]!=eData) {
	weight *= scale1fb*lumi;
      }
      
      // fill Z events passing selection
      if((category==eMuMu2HLT) || (category==eMuMu1HLT) || (category==eMuMu1HLT1L1)) {
        if(typev[ifile]==eData) { 

	  TLorentzVector mu1;
	  TLorentzVector mu2;
	  mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
	  mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
	  float qter1=1.0;
	  float qter2=1.0;

	  rmcor->momcor_data(mu1,q1,0,qter1);
	  rmcor->momcor_data(mu2,q2,0,qter2);

	  Double_t lp1 = mu1.Pt();
	  Double_t lp2 = mu2.Pt();
	  Double_t lq1 = q1;
	  Double_t lq2 = q2;

	  TLorentzVector l1, l2;
	  if(lp1>lp2)
	    {
	      l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),mu_MASS);
	      l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),mu_MASS);
	    }
	  else
	    {
	      l1.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),mu_MASS);
	      l2.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),mu_MASS);
	      lq1=q2;
	      lq2=q1;
	    }

	  mass=(l1+l2).M();
	  pt =(l1+l2).Pt();
	  rapidity = (l1+l2).Rapidity();

	  phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
	  if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
	  else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
	  phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));
	  
	  if(mass        < MASS_LOW)  continue;
	  if(mass        > MASS_HIGH) continue;
	  if(l1.Pt()        < PT_CUT)    continue;
	  if(l2.Pt()        < PT_CUT)    continue;

	  hData->Fill(mass); 
	  hDataNPV->Fill(npv);
	  hDataZPt->Fill(pt); 
	  hDataPhiStar->Fill(phistar); 
	  hDataLep1Pt->Fill(l1.Pt()); 
	  hDataLep2Pt->Fill(l2.Pt()); 
	  if(lq1<0)
	    {
	      hDataLepNegPt->Fill(l1.Pt()); 
	      hDataLepPosPt->Fill(l2.Pt());
	    }
	  else 
	    {
	      hDataLepNegPt->Fill(l2.Pt()); 
	      hDataLepPosPt->Fill(l1.Pt());
	    }
	  hDataLep1Eta->Fill(fabs(l1.Eta())); 
	  hDataLep2Eta->Fill(fabs(l2.Eta())); 
	  hDataZRap->Fill(fabs(rapidity));
	  
	  yield++;
	
	} else {

	  TLorentzVector mu1;
	  TLorentzVector mu2;
	  mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
	  mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
	  float qter1=1.0;
	  float qter2=1.0;

	  rmcor->momcor_mc(mu1,q1,0,qter1);
	  rmcor->momcor_mc(mu2,q2,0,qter2);

	  Double_t lp1 = mu1.Pt();
	  Double_t lp2 = mu2.Pt();
	  Double_t lq1 = q1;
	  Double_t lq2 = q2;

	  TLorentzVector l1, l2;
	  if(lp1>lp2)
	    {
	      l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),mu_MASS);
	      l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),mu_MASS);
	    }
	  else
	    {
	      l1.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),mu_MASS);
	      l2.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),mu_MASS);
	      lq1=q2;
	      lq2=q1;
	    }

	  double mll=(l1+l2).M();
	  Double_t effdata, effmc;
	  Double_t corr=1;
	  Double_t eff2Bindata, eff2Binmc;
	  Double_t corr2Bin=1;
	  Double_t corrUp=1;
	  Double_t corrDown=1;
	  Double_t effSigShapedata;
	  Double_t corrSigShape=1;
	  Double_t effBkgShapedata;
	  Double_t corrBkgShape=1;


	  if(mll       < MASS_LOW)  continue;
	  if(mll       > MASS_HIGH) continue;
	  if(lp1        < PT_CUT)    continue;
	  if(lp2        < PT_CUT)    continue;
	  effdata=1; effmc=1;    
          if(lq1>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff((l1.Eta()), l1.Pt())); 
            effmc   *= (1.-zmmHLTEff_pos.getEff((l1.Eta()), l1.Pt())); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff((l1.Eta()), l1.Pt())); 
            effmc   *= (1.-zmmHLTEff_neg.getEff((l1.Eta()), l1.Pt())); 
          }
          if(lq2>0) {
            effdata *= (1.-dataHLTEff_pos.getEff((l2.Eta()), l2.Pt())); 
            effmc   *= (1.-zmmHLTEff_pos.getEff((l2.Eta()), l2.Pt()));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff((l2.Eta()), l2.Pt())); 
            effmc   *= (1.-zmmHLTEff_neg.getEff((l2.Eta()), l2.Pt()));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;
	  corrSigShape *= effdata/effmc;
	  corrBkgShape *= effdata/effmc;
    
          effdata=1; effmc=1;
	  effSigShapedata=1;
	  effBkgShapedata=1;
	  if(lq1>0) { 
	    effdata *= dataSelEff_pos.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmSelEff_pos.getEff((l1.Eta()), l1.Pt()); 
	    effSigShapedata *= dataSelEff_pos.getEff((l1.Eta()), l1.Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(l1.Eta()), hSelSigSys->GetYaxis()->FindBin(l1.Pt())); 
	    effBkgShapedata *= dataSelEff_pos.getEff((l1.Eta()), l1.Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(l1.Eta()), hSelBkgSys->GetYaxis()->FindBin(l1.Pt())); 
	  } else {
	    effdata *= dataSelEff_neg.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmSelEff_neg.getEff((l1.Eta()), l1.Pt()); 
	    effSigShapedata *= dataSelEff_neg.getEff((l1.Eta()), l1.Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(l1.Eta()), hSelSigSys->GetYaxis()->FindBin(l1.Pt()));
	    effBkgShapedata *= dataSelEff_neg.getEff((l1.Eta()), l1.Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(l1.Eta()), hSelBkgSys->GetYaxis()->FindBin(l1.Pt())); 
	  }
	  if(lq2>0) {
	    effdata *= dataSelEff_pos.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmSelEff_pos.getEff((l2.Eta()), l2.Pt());
	    effSigShapedata *= dataSelEff_pos.getEff((l2.Eta()), l2.Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(l2.Eta()), hSelSigSys->GetYaxis()->FindBin(l2.Pt())); 
	    effBkgShapedata *= dataSelEff_pos.getEff((l2.Eta()), l2.Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(l2.Eta()), hSelBkgSys->GetYaxis()->FindBin(l2.Pt())); 
	  } else {
	    effdata *= dataSelEff_neg.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmSelEff_neg.getEff((l2.Eta()), l2.Pt());
	    effSigShapedata *= dataSelEff_neg.getEff((l2.Eta()), l2.Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(l2.Eta()), hSelSigSys->GetYaxis()->FindBin(l2.Pt()));
	    effBkgShapedata *= dataSelEff_neg.getEff((l2.Eta()), l2.Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(l2.Eta()), hSelBkgSys->GetYaxis()->FindBin(l2.Pt()));
	  }
	  corr *= effdata/effmc;
	  corrSigShape *= effSigShapedata/effmc;
	  corrBkgShape *= effBkgShapedata/effmc;
	  
	  effdata=1; effmc=1;
	  effSigShapedata=1;
	  effBkgShapedata=1;
	  if(lq1>0) { 
	    effdata *= dataStaEff_pos.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmStaEff_pos.getEff((l1.Eta()), l1.Pt()); 
	    effSigShapedata *= dataStaEff_neg.getEff((l1.Eta()), l1.Pt())*hStaSigSys->GetBinContent(hStaSigSys->GetXaxis()->FindBin(l1.Eta()), hStaSigSys->GetYaxis()->FindBin(l1.Pt()));
	    effBkgShapedata *= dataStaEff_neg.getEff((l1.Eta()), l1.Pt())*hStaBkgSys->GetBinContent(hStaBkgSys->GetXaxis()->FindBin(l1.Eta()), hStaBkgSys->GetYaxis()->FindBin(l1.Pt())); 
	  } else {
	    effdata *= dataStaEff_neg.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmStaEff_neg.getEff((l1.Eta()), l1.Pt()); 
	    effSigShapedata *= dataStaEff_neg.getEff((l1.Eta()), l1.Pt())*hStaSigSys->GetBinContent(hStaSigSys->GetXaxis()->FindBin(l1.Eta()), hStaSigSys->GetYaxis()->FindBin(l1.Pt()));
	    effBkgShapedata *= dataStaEff_neg.getEff((l1.Eta()), l1.Pt())*hStaBkgSys->GetBinContent(hStaBkgSys->GetXaxis()->FindBin(l1.Eta()), hStaBkgSys->GetYaxis()->FindBin(l1.Pt())); 
	  }
	  if(lq2>0) {
	    effdata *= dataStaEff_pos.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmStaEff_pos.getEff((l2.Eta()), l2.Pt());
	    effSigShapedata *= dataStaEff_pos.getEff((l2.Eta()), l2.Pt())*hStaSigSys->GetBinContent(hStaSigSys->GetXaxis()->FindBin(l2.Eta()), hStaSigSys->GetYaxis()->FindBin(l2.Pt())); 
	    effBkgShapedata *= dataStaEff_pos.getEff((l2.Eta()), l2.Pt())*hStaBkgSys->GetBinContent(hStaBkgSys->GetXaxis()->FindBin(l2.Eta()), hStaBkgSys->GetYaxis()->FindBin(l2.Pt())); 
	  } else {
	    effdata *= dataStaEff_neg.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmStaEff_neg.getEff((l2.Eta()), l2.Pt());
	    effSigShapedata *= dataStaEff_neg.getEff((l2.Eta()), l2.Pt())*hStaSigSys->GetBinContent(hStaSigSys->GetXaxis()->FindBin(l2.Eta()), hStaSigSys->GetYaxis()->FindBin(l2.Pt()));
	    effBkgShapedata *= dataStaEff_neg.getEff((l2.Eta()), l2.Pt())*hStaBkgSys->GetBinContent(hStaBkgSys->GetXaxis()->FindBin(l2.Eta()), hStaBkgSys->GetYaxis()->FindBin(l2.Pt()));
	  }
	  corr *= effdata/effmc; 
	  corrSigShape *= effSigShapedata/effmc;
	  corrBkgShape *= effBkgShapedata/effmc;
	  
          effdata=1; effmc=1;
          if(lq1>0) { 
            effdata *= dataTrkEff_pos.getEff((l1.Eta()), l1.Pt()); 
            effmc   *= zmmTrkEff_pos.getEff((l1.Eta()), l1.Pt()); 
          } else {
            effdata *= dataTrkEff_neg.getEff((l1.Eta()), l1.Pt()); 
            effmc   *= zmmTrkEff_neg.getEff((l1.Eta()), l1.Pt()); 
          }
          if(lq2>0) {
            effdata *= dataTrkEff_pos.getEff((l2.Eta()), l2.Pt()); 
            effmc   *= zmmTrkEff_pos.getEff((l2.Eta()), l2.Pt());
          } else {
            effdata *= dataTrkEff_neg.getEff((l2.Eta()), l2.Pt()); 
            effmc   *= zmmTrkEff_neg.getEff((l2.Eta()), l2.Pt());
          }
          //corr *= effdata/effmc;
	  //corr=1;


	  // scale factor uncertainties   
	  
	  double var=0.;      
	  
	  // TRACKER
	  if(lq1>0) {
	    Double_t effdata = dataTrkEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataTrkEff_pos.getErrLow(l1.Eta(), l1.Pt()), dataTrkEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmTrkEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmTrkEff_pos.getErrLow(l1.Eta(), l1.Pt()), zmmTrkEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    //var+=errTrk*errTrk;
	  } else {
	    Double_t effdata = dataTrkEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataTrkEff_neg.getErrLow(l1.Eta(), l1.Pt()), dataTrkEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmTrkEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmTrkEff_neg.getErrLow(l1.Eta(), l1.Pt()), zmmTrkEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    //var+=errTrk*errTrk;
	  }
	  
	  if(lq2>0) {
	    Double_t effdata = dataTrkEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataTrkEff_pos.getErrLow(l2.Eta(), l2.Pt()), dataTrkEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmTrkEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmTrkEff_pos.getErrLow(l2.Eta(), l2.Pt()), zmmTrkEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    //var+=errTrk*errTrk;
	  } else {
	    Double_t effdata = dataTrkEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataTrkEff_neg.getErrLow(l2.Eta(), l2.Pt()), dataTrkEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmTrkEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmTrkEff_neg.getErrLow(l2.Eta(), l2.Pt()), zmmTrkEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    //var+=errTrk*errTrk;
	  }
	  
	  // STANDALONE
	  if(lq1>0) {
	    Double_t effdata = dataStaEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(l1.Eta(), l1.Pt()), dataStaEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmStaEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(l1.Eta(), l1.Pt()), zmmStaEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSta*errSta;
	  } else {
	    Double_t effdata = dataStaEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(l1.Eta(), l1.Pt()), dataStaEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmStaEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(l1.Eta(), l1.Pt()), zmmStaEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSta*errSta;
	  }
	  
	  if(lq2>0) {
	    Double_t effdata = dataStaEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(l2.Eta(), l2.Pt()), dataStaEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmStaEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(l2.Eta(), l2.Pt()), zmmStaEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errSta = ((effdata/effmc))*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSta*errSta;
	  } else {
	    Double_t effdata = dataStaEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(l2.Eta(), l2.Pt()), dataStaEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmStaEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(l2.Eta(), l2.Pt()), zmmStaEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSta*errSta;
	  }
	  
	  // SELECTION
	  if(lq1>0) {
	    Double_t effdata = dataSelEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(l1.Eta(), l1.Pt()), dataSelEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmSelEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(l1.Eta(), l1.Pt()), zmmSelEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSel*errSel;
	  } else {
	    Double_t effdata = dataSelEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(l1.Eta(), l1.Pt()), dataSelEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmSelEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(l1.Eta(), l1.Pt()), zmmSelEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSel*errSel;
	  }
	  
	  if(lq2>0) {
	    Double_t effdata = dataSelEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(l2.Eta(), l2.Pt()), dataSelEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmSelEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(l2.Eta(), l2.Pt()), zmmSelEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSel*errSel;
	  } else {
	    Double_t effdata = dataSelEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(l2.Eta(), l2.Pt()), dataSelEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmSelEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(l2.Eta(), l2.Pt()), zmmSelEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errSel*errSel;
	  }

	  
	  //HLT
	  if(lq1>0) {
	    Double_t effdata = dataHLTEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(l1.Eta(), l1.Pt()), dataHLTEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmHLTEff_pos.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(l1.Eta(), l1.Pt()), zmmHLTEff_pos.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errHLT*errHLT;
	  } else {
	    Double_t effdata = dataHLTEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(l1.Eta(), l1.Pt()), dataHLTEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t effmc   = zmmHLTEff_neg.getEff(l1.Eta(), l1.Pt());
	    Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(l1.Eta(), l1.Pt()), zmmHLTEff_neg.getErrHigh(l1.Eta(), l1.Pt()));
	    Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errHLT*errHLT;
	  }
	  
	  if(lq2>0) {
	    Double_t effdata = dataHLTEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(l2.Eta(), l2.Pt()), dataHLTEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmHLTEff_pos.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(l2.Eta(), l2.Pt()), zmmHLTEff_pos.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errHLT*errHLT;
	  } else {
	    Double_t effdata = dataHLTEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(l2.Eta(), l2.Pt()), dataHLTEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t effmc   = zmmHLTEff_neg.getEff(l2.Eta(), l2.Pt());
	    Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(l2.Eta(), l2.Pt()), zmmHLTEff_neg.getErrHigh(l2.Eta(), l2.Pt()));
	    Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
	    var+=errHLT*errHLT;
	  }
	  
	  corrUp=corr+sqrt(var);
	  corrDown=corr-sqrt(var);  

	  eff2Bindata=1; eff2Binmc=1;    
          if(lq1>0) { 
            eff2Bindata *= (1.-dataHLTEff2Bin_pos.getEff((l1.Eta()), l1.Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_pos.getEff((l1.Eta()), l1.Pt())); 
          } else {
            eff2Bindata *= (1.-dataHLTEff2Bin_neg.getEff((l1.Eta()), l1.Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_neg.getEff((l1.Eta()), l1.Pt())); 
          }
          if(lq2>0) {
            eff2Bindata *= (1.-dataHLTEff2Bin_pos.getEff((l2.Eta()), l2.Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_pos.getEff((l2.Eta()), l2.Pt()));
          } else {
            eff2Bindata *= (1.-dataHLTEff2Bin_neg.getEff((l2.Eta()), l2.Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_neg.getEff((l2.Eta()), l2.Pt()));
          }
          eff2Bindata = 1.-eff2Bindata;
          eff2Binmc   = 1.-eff2Binmc;
          corr2Bin *= eff2Bindata/eff2Binmc;
    
          eff2Bindata=1; eff2Binmc=1;
          if(lq1>0) { 
            eff2Bindata *= dataSelEff2Bin_pos.getEff((l1.Eta()), l1.Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_pos.getEff((l1.Eta()), l1.Pt()); 
          } else {
            eff2Bindata *= dataSelEff2Bin_neg.getEff((l1.Eta()), l1.Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_neg.getEff((l1.Eta()), l1.Pt()); 
          }
          if(lq2>0) {
            eff2Bindata *= dataSelEff2Bin_pos.getEff((l2.Eta()), l2.Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_pos.getEff((l2.Eta()), l2.Pt());
          } else {
            eff2Bindata *= dataSelEff2Bin_neg.getEff((l2.Eta()), l2.Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_neg.getEff((l2.Eta()), l2.Pt());
          }
          corr2Bin *= eff2Bindata/eff2Binmc;
    
	  eff2Bindata=1; eff2Binmc=1;
          if(lq1>0) { 
            eff2Bindata *= dataStaEff2Bin_pos.getEff((l1.Eta()), l1.Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_pos.getEff((l1.Eta()), l1.Pt()); 
          } else {
            eff2Bindata *= dataStaEff2Bin_neg.getEff((l1.Eta()), l1.Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_neg.getEff((l1.Eta()), l1.Pt()); 
          }
          if(lq2>0) {
            eff2Bindata *= dataStaEff2Bin_pos.getEff((l2.Eta()), l2.Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_pos.getEff((l2.Eta()), l2.Pt());
          } else {
            eff2Bindata *= dataStaEff2Bin_neg.getEff((l2.Eta()), l2.Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_neg.getEff((l2.Eta()), l2.Pt());
          }
	  corr2Bin *= eff2Bindata/eff2Binmc; 
	  
          eff2Bindata=1; eff2Binmc=1;
          if(lq1>0) { 
            eff2Bindata *= dataTrkEff2Bin_pos.getEff((l1.Eta()), l1.Pt()); 
            eff2Binmc   *= zmmTrkEff2Bin_pos.getEff((l1.Eta()), l1.Pt()); 
          } else {
            eff2Bindata *= dataTrkEff2Bin_neg.getEff((l1.Eta()), l1.Pt()); 
            eff2Binmc   *= zmmTrkEff2Bin_neg.getEff((l1.Eta()), l1.Pt()); 
          }
          if(lq2>0) {
            eff2Bindata *= dataTrkEff2Bin_pos.getEff((l2.Eta()), l2.Pt()); 
            eff2Binmc   *= zmmTrkEff2Bin_pos.getEff((l2.Eta()), l2.Pt());
          } else {
            eff2Bindata *= dataTrkEff2Bin_neg.getEff((l2.Eta()), l2.Pt()); 
            eff2Binmc   *= zmmTrkEff2Bin_neg.getEff((l2.Eta()), l2.Pt());
          }
          //corr2Bin *= eff2Bindata/eff2Binmc;
	  //corr2Bin=1;
	
	  mass = (l1+l2).M();
	  pt = (l1+l2).Pt();
	  rapidity = (l1+l2).Rapidity();

	  phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
	  if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
	  else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
	  phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));

	  if(typev[ifile]==eZmm) 
	    {
	      yield_zmm += weight*corr;
	      yield_zmm_unc += weight*weight*corr*corr;
	      hZmm->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
	      hZmmNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);
	      hZmmZPt->Fill(pt,weight*corr); 
	      hMCZPt->Fill(pt,weight*corr);
	      hZmmPhiStar->Fill(phistar,weight*corr); 
	      hMCPhiStar->Fill(phistar,weight*corr);
	      hZmmZRap->Fill(fabs(rapidity),weight*corr); 
	      hMCZRap->Fill(fabs(rapidity),weight*corr);
	      hZmmLep1Pt->Fill(l1.Pt(),weight*corr); 
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);
	      hZmmLep2Pt->Fill(l2.Pt(),weight*corr); 
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);
	      if(lq1<0)
		{
		  hZmmLepNegPt->Fill(l1.Pt(),weight*corr); 
		  hMCLepNegPt->Fill(l1.Pt(),weight*corr);
		  hZmmLepPosPt->Fill(l2.Pt(),weight*corr); 
		  hMCLepPosPt->Fill(l2.Pt(),weight*corr);
		}
	      else 
		{
		  hZmmLepNegPt->Fill(l2.Pt(),weight*corr); 
		  hMCLepNegPt->Fill(l2.Pt(),weight*corr);
		  hZmmLepPosPt->Fill(l1.Pt(),weight*corr); 
		  hMCLepPosPt->Fill(l1.Pt(),weight*corr);
		}
	      hZmmLep1Eta->Fill(fabs(l1.Eta()),weight*corr); 
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
	      hZmmLep2Eta->Fill(fabs(l2.Eta()),weight*corr); 
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
	  if(typev[ifile]==eEWK) 
	    {
	      yield_ewk += weight*corr;
	      yield_ewk_unc += weight*weight*corr*corr;
	      hEWK->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);

	      hEWKNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);

	      hEWKZPt->Fill(pt,weight*corr); 
	      hEWKZPt_EffBin->Fill(pt,weight*corr2Bin); 
	      hEWKZPt_EffStatUp->Fill(pt,weight*corrUp);
	      hEWKZPt_EffStatDown->Fill(pt,weight*corrDown);
	      hEWKZPt_EffSigShape->Fill(pt,weight*corrSigShape);
	      hEWKZPt_EffBkgShape->Fill(pt,weight*corrBkgShape);
	      hMCZPt->Fill(pt,weight*corr);

	      hEWKPhiStar->Fill(phistar,weight*corr);
	      hEWKPhiStar_EffBin->Fill(pt,weight*corr2Bin); 
	      hEWKPhiStar_EffStatUp->Fill(pt,weight*corrUp);
	      hEWKPhiStar_EffStatDown->Fill(pt,weight*corrDown);
	      hEWKPhiStar_EffSigShape->Fill(pt,weight*corrSigShape);
	      hEWKPhiStar_EffBkgShape->Fill(pt,weight*corrBkgShape);
	      hMCPhiStar->Fill(phistar,weight*corr);

	      hEWKZRap->Fill(fabs(rapidity),weight*corr); 
	      hEWKZRap_EffBin->Fill(fabs(rapidity),weight*corr2Bin);
	      hEWKZRap_EffStatUp->Fill(fabs(rapidity),weight*corrUp);
	      hEWKZRap_EffStatDown->Fill(fabs(rapidity),weight*corrDown);
	      hEWKZRap_EffSigShape->Fill(fabs(rapidity),weight*corrSigShape);
	      hEWKZRap_EffBkgShape->Fill(fabs(rapidity),weight*corrBkgShape);
	      hMCZRap->Fill(fabs(rapidity),weight*corr);

	      hEWKLep1Pt->Fill(l1.Pt(),weight*corr);
	      hEWKLep1Pt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
	      hEWKLep1Pt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
	      hEWKLep1Pt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
	      hEWKLep1Pt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
	      hEWKLep1Pt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);

	      hEWKLep2Pt->Fill(l2.Pt(),weight*corr);
	      hEWKLep2Pt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
	      hEWKLep2Pt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
	      hEWKLep2Pt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
	      hEWKLep2Pt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
	      hEWKLep2Pt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);

	      if(lq1<0)
		{
		  hEWKLepNegPt->Fill(l1.Pt(),weight*corr);
		  hEWKLepNegPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
		  hEWKLepNegPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
		  hEWKLepNegPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
		  hEWKLepNegPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
		  hEWKLepNegPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
		  hMCLepNegPt->Fill(l1.Pt(),weight*corr);
		  
		  hEWKLepPosPt->Fill(l2.Pt(),weight*corr);
		  hEWKLepPosPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
		  hEWKLepPosPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
		  hEWKLepPosPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
		  hEWKLepPosPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
		  hEWKLepPosPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
		  hMCLepPosPt->Fill(l2.Pt(),weight*corr);
		}
	      else
		{
		  hEWKLepNegPt->Fill(l2.Pt(),weight*corr);
		  hEWKLepNegPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
		  hEWKLepNegPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
		  hEWKLepNegPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
		  hEWKLepNegPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
		  hEWKLepNegPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
		  hMCLepNegPt->Fill(l2.Pt(),weight*corr);
		  
		  hEWKLepPosPt->Fill(l1.Pt(),weight*corr);
		  hEWKLepPosPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
		  hEWKLepPosPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
		  hEWKLepPosPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
		  hEWKLepPosPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
		  hEWKLepPosPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
		  hMCLepPosPt->Fill(l1.Pt(),weight*corr);
		}

	      hEWKLep1Eta->Fill(fabs(l1.Eta()),weight*corr); 
	      hEWKLep1Eta_EffBin->Fill(fabs(l1.Eta()),weight*corr2Bin);
	      hEWKLep1Eta_EffStatUp->Fill(fabs(l1.Eta()),weight*corrUp);
	      hEWKLep1Eta_EffStatDown->Fill(fabs(l1.Eta()),weight*corrDown);
	      hEWKLep1Eta_EffSigShape->Fill(fabs(l1.Eta()),weight*corrSigShape);
	      hEWKLep1Eta_EffBkgShape->Fill(fabs(l1.Eta()),weight*corrBkgShape);
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);

	      hEWKLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	      hEWKLep2Eta_EffBin->Fill(fabs(l2.Eta()),weight*corr2Bin);
	      hEWKLep2Eta_EffStatUp->Fill(fabs(l2.Eta()),weight*corrUp);
	      hEWKLep2Eta_EffStatDown->Fill(fabs(l2.Eta()),weight*corrDown);
	      hEWKLep2Eta_EffSigShape->Fill(fabs(l2.Eta()),weight*corrSigShape);
	      hEWKLep2Eta_EffBkgShape->Fill(fabs(l2.Eta()),weight*corrBkgShape);
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
	  if(typev[ifile]==eTop) 
	    {
	      yield_ewk += weight*corr;
	      yield_ewk_unc += weight*weight*corr*corr;
	      hTop->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);

	      hTopNPV->Fill(npv,weight*corr); 
	      hMCNPV->Fill(npv,weight*corr);

	      hTopZPt->Fill(pt,weight*corr); 
	      hTopZPt_EffBin->Fill(pt,weight*corr2Bin); 
	      hTopZPt_EffStatUp->Fill(pt,weight*corrUp);
	      hTopZPt_EffStatDown->Fill(pt,weight*corrDown);
	      hTopZPt_EffSigShape->Fill(pt,weight*corrSigShape); 
	      hTopZPt_EffBkgShape->Fill(pt,weight*corrBkgShape); 
	      hMCZPt->Fill(pt,weight*corr);

	      hTopPhiStar->Fill(phistar,weight*corr);
	      hTopPhiStar_EffBin->Fill(pt,weight*corr2Bin); 
	      hTopPhiStar_EffStatUp->Fill(pt,weight*corrUp);
	      hTopPhiStar_EffStatDown->Fill(pt,weight*corrDown);
	      hTopPhiStar_EffSigShape->Fill(pt,weight*corrSigShape); 
	      hTopPhiStar_EffBkgShape->Fill(pt,weight*corrBkgShape); 
	      hMCPhiStar->Fill(phistar,weight*corr);

	      hTopZRap->Fill(fabs(rapidity),weight*corr); 
	      hTopZRap_EffBin->Fill(fabs(rapidity),weight*corr2Bin);
	      hTopZRap_EffStatUp->Fill(fabs(rapidity),weight*corrUp);
	      hTopZRap_EffStatDown->Fill(fabs(rapidity),weight*corrDown);
	      hTopZRap_EffSigShape->Fill(fabs(rapidity),weight*corrSigShape); 
	      hTopZRap_EffBkgShape->Fill(fabs(rapidity),weight*corrBkgShape); 
	      hMCZRap->Fill(fabs(rapidity),weight*corr);

	      hTopLep1Pt->Fill(l1.Pt(),weight*corr);
	      hTopLep1Pt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
	      hTopLep1Pt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
	      hTopLep1Pt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
	      hTopLep1Pt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
	      hTopLep1Pt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);

	      hTopLep2Pt->Fill(l2.Pt(),weight*corr);
	      hTopLep2Pt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
	      hTopLep2Pt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
	      hTopLep2Pt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
	      hTopLep2Pt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
	      hTopLep2Pt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);

	      if(lq1<0)
		{
		  hTopLepNegPt->Fill(l1.Pt(),weight*corr);
		  hTopLepNegPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
		  hTopLepNegPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
		  hTopLepNegPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
		  hTopLepNegPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
		  hTopLepNegPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
		  hMCLepNegPt->Fill(l1.Pt(),weight*corr);
		  
		  hTopLepPosPt->Fill(l2.Pt(),weight*corr);
		  hTopLepPosPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
		  hTopLepPosPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
		  hTopLepPosPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
		  hTopLepPosPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
		  hTopLepPosPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
		  hMCLepPosPt->Fill(l2.Pt(),weight*corr);
		}
	      else
		{
		  hTopLepNegPt->Fill(l2.Pt(),weight*corr);
		  hTopLepNegPt_EffBin->Fill(l2.Pt(),weight*corr2Bin);
		  hTopLepNegPt_EffStatUp->Fill(l2.Pt(),weight*corrUp);
		  hTopLepNegPt_EffStatDown->Fill(l2.Pt(),weight*corrDown);
		  hTopLepNegPt_EffSigShape->Fill(l2.Pt(),weight*corrSigShape);
		  hTopLepNegPt_EffBkgShape->Fill(l2.Pt(),weight*corrBkgShape);
		  hMCLepNegPt->Fill(l2.Pt(),weight*corr);
		  
		  hTopLepPosPt->Fill(l1.Pt(),weight*corr);
		  hTopLepPosPt_EffBin->Fill(l1.Pt(),weight*corr2Bin);
		  hTopLepPosPt_EffStatUp->Fill(l1.Pt(),weight*corrUp);
		  hTopLepPosPt_EffStatDown->Fill(l1.Pt(),weight*corrDown);
		  hTopLepPosPt_EffSigShape->Fill(l1.Pt(),weight*corrSigShape);
		  hTopLepPosPt_EffBkgShape->Fill(l1.Pt(),weight*corrBkgShape);
		  hMCLepPosPt->Fill(l1.Pt(),weight*corr);
		}

	      hTopLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
	      hTopLep1Eta_EffBin->Fill(fabs(l1.Eta()),weight*corr2Bin);
	      hTopLep1Eta_EffStatUp->Fill(fabs(l1.Eta()),weight*corrUp);
	      hTopLep1Eta_EffStatDown->Fill(fabs(l1.Eta()),weight*corrDown);
	      hTopLep1Eta_EffSigShape->Fill(fabs(l1.Eta()),weight*corrSigShape);
	      hTopLep1Eta_EffBkgShape->Fill(fabs(l1.Eta()),weight*corrBkgShape);
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);

	      hTopLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	      hTopLep2Eta_EffBin->Fill(fabs(l2.Eta()),weight*corr2Bin);
	      hTopLep2Eta_EffStatUp->Fill(fabs(l2.Eta()),weight*corrUp);
	      hTopLep2Eta_EffStatDown->Fill(fabs(l2.Eta()),weight*corrDown);
	      hTopLep2Eta_EffSigShape->Fill(fabs(l2.Eta()),weight*corrSigShape);
	      hTopLep2Eta_EffBkgShape->Fill(fabs(l2.Eta()),weight*corrBkgShape);
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
	}
      }
    }
    
    delete infile;
    infile=0, intree=0;
  } 

  outFile->cd();

  hDataZPt->Write();
  hEWKZPt->Write();
  hEWKZPt_EffBin->Write();
  hEWKZPt_EffStatUp->Write();
  hEWKZPt_EffStatDown->Write();
  hEWKZPt_EffSigShape->Write();
  hEWKZPt_EffBkgShape->Write();
  hTopZPt->Write();
  hTopZPt_EffBin->Write();
  hTopZPt_EffStatUp->Write();
  hTopZPt_EffStatDown->Write();
  hTopZPt_EffSigShape->Write();
  hTopZPt_EffBkgShape->Write();

  hDataPhiStar->Write();
  hEWKPhiStar->Write();
  hEWKPhiStar_EffBin->Write();
  hEWKPhiStar_EffStatUp->Write();
  hEWKPhiStar_EffStatDown->Write();
  hEWKPhiStar_EffSigShape->Write();
  hEWKPhiStar_EffBkgShape->Write();
  hTopPhiStar->Write();
  hTopPhiStar_EffBin->Write();
  hTopPhiStar_EffStatUp->Write();
  hTopPhiStar_EffStatDown->Write();
  hTopPhiStar_EffSigShape->Write();
  hTopPhiStar_EffBkgShape->Write();

 
  hDataZRap->Write();
  hEWKZRap->Write();
  hEWKZRap_EffBin->Write();
  hEWKZRap_EffStatUp->Write();
  hEWKZRap_EffStatDown->Write();
  hEWKZRap_EffSigShape->Write();
  hEWKZRap_EffBkgShape->Write();
  hTopZRap->Write();
  hTopZRap_EffBin->Write();
  hTopZRap_EffStatUp->Write();
  hTopZRap_EffStatDown->Write();
  hTopZRap_EffSigShape->Write();
  hTopZRap_EffBkgShape->Write();

  hDataLep1Pt->Write();
  hEWKLep1Pt->Write();
  hEWKLep1Pt_EffBin->Write();
  hEWKLep1Pt_EffStatUp->Write();
  hEWKLep1Pt_EffStatDown->Write();
  hEWKLep1Pt_EffSigShape->Write();
  hEWKLep1Pt_EffBkgShape->Write();
  hTopLep1Pt->Write();
  hTopLep1Pt_EffBin->Write();
  hTopLep1Pt_EffStatUp->Write();
  hTopLep1Pt_EffStatDown->Write();
  hTopLep1Pt_EffSigShape->Write();
  hTopLep1Pt_EffBkgShape->Write();

  hDataLep2Pt->Write();
  hEWKLep2Pt->Write();
  hEWKLep2Pt_EffBin->Write();
  hEWKLep2Pt_EffStatUp->Write();
  hEWKLep2Pt_EffStatDown->Write();
  hEWKLep2Pt_EffSigShape->Write();
  hEWKLep2Pt_EffBkgShape->Write();
  hTopLep2Pt->Write();
  hTopLep2Pt_EffBin->Write();
  hTopLep2Pt_EffStatUp->Write();
  hTopLep2Pt_EffStatDown->Write();
  hTopLep2Pt_EffSigShape->Write();
  hTopLep2Pt_EffBkgShape->Write();

  hDataLepNegPt->Write();
  hEWKLepNegPt->Write();
  hEWKLepNegPt_EffBin->Write();
  hEWKLepNegPt_EffStatUp->Write();
  hEWKLepNegPt_EffStatDown->Write();
  hEWKLepNegPt_EffSigShape->Write();
  hEWKLepNegPt_EffBkgShape->Write();
  hTopLepNegPt->Write();
  hTopLepNegPt_EffBin->Write();
  hTopLepNegPt_EffStatUp->Write();
  hTopLepNegPt_EffStatDown->Write();
  hTopLepNegPt_EffSigShape->Write();
  hTopLepNegPt_EffBkgShape->Write();

  hDataLepPosPt->Write();
  hEWKLepPosPt->Write();
  hEWKLepPosPt_EffBin->Write();
  hEWKLepPosPt_EffStatUp->Write();
  hEWKLepPosPt_EffStatDown->Write();
  hEWKLepPosPt_EffSigShape->Write();
  hEWKLepPosPt_EffBkgShape->Write();
  hTopLepPosPt->Write();
  hTopLepPosPt_EffBin->Write();
  hTopLepPosPt_EffStatUp->Write();
  hTopLepPosPt_EffStatDown->Write();
  hTopLepPosPt_EffSigShape->Write();
  hTopLepPosPt_EffBkgShape->Write();

  hDataLep1Eta->Write();
  hEWKLep1Eta->Write();
  hEWKLep1Eta_EffBin->Write();
  hEWKLep1Eta_EffStatUp->Write();
  hEWKLep1Eta_EffStatDown->Write();
  hEWKLep1Eta_EffSigShape->Write();
  hEWKLep1Eta_EffBkgShape->Write();
  hTopLep1Eta->Write();
  hTopLep1Eta_EffBin->Write();
  hTopLep1Eta_EffStatUp->Write();
  hTopLep1Eta_EffStatDown->Write();
  hTopLep1Eta_EffSigShape->Write();
  hTopLep1Eta_EffBkgShape->Write();

  hDataLep2Eta->Write();
  hEWKLep2Eta->Write();
  hEWKLep2Eta_EffBin->Write();
  hEWKLep2Eta_EffStatUp->Write();
  hEWKLep2Eta_EffStatDown->Write();
  hEWKLep2Eta_EffSigShape->Write();
  hEWKLep2Eta_EffBkgShape->Write();
  hTopLep2Eta->Write();
  hTopLep2Eta_EffBin->Write();
  hTopLep2Eta_EffStatUp->Write();
  hTopLep2Eta_EffStatDown->Write();
  hTopLep2Eta_EffSigShape->Write();
  hTopLep2Eta_EffBkgShape->Write();

  outFile->Write();
  outFile->Close(); 

  double MCscale=hData->Integral()/hMC->Integral();

  if(normToData)
    {

      cout<<"Normalized to data: "<<MCscale<<endl;
      
      hZmm->Scale(MCscale);
      hMC->Scale(MCscale);
      hEWK->Scale(MCscale);
      hTop->Scale(MCscale);
      hZmmNPV->Scale(MCscale);
      hMCNPV->Scale(MCscale);
      hEWKNPV->Scale(MCscale);
      hTopNPV->Scale(MCscale);
      hZmmZPt->Scale(MCscale);
      hMCZPt->Scale(MCscale);
      hEWKZPt->Scale(MCscale);
      hTopZPt->Scale(MCscale);
      hZmmPhiStar->Scale(MCscale);
      hMCPhiStar->Scale(MCscale);
      hEWKPhiStar->Scale(MCscale);
      hTopPhiStar->Scale(MCscale);
      hZmmZRap->Scale(MCscale);
      hMCZRap->Scale(MCscale);
      hEWKZRap->Scale(MCscale);
      hTopZRap->Scale(MCscale);
      hZmmLep1Pt->Scale(MCscale);
      hMCLep1Pt->Scale(MCscale);
      hEWKLep1Pt->Scale(MCscale);
      hTopLep1Pt->Scale(MCscale);
      hZmmLep2Pt->Scale(MCscale);
      hMCLep2Pt->Scale(MCscale);
      hEWKLep2Pt->Scale(MCscale);
      hTopLep2Pt->Scale(MCscale);
      hZmmLepNegPt->Scale(MCscale);
      hMCLepNegPt->Scale(MCscale);
      hEWKLepNegPt->Scale(MCscale);
      hTopLepNegPt->Scale(MCscale);
      hZmmLepPosPt->Scale(MCscale);
      hMCLepPosPt->Scale(MCscale);
      hEWKLepPosPt->Scale(MCscale);
      hTopLepPosPt->Scale(MCscale);
      hZmmLep1Eta->Scale(MCscale);
      hMCLep1Eta->Scale(MCscale);
      hEWKLep1Eta->Scale(MCscale);
      hTopLep1Eta->Scale(MCscale);
      hZmmLep2Eta->Scale(MCscale);
      hMCLep2Eta->Scale(MCscale);
      hEWKLep2Eta->Scale(MCscale);
      hTopLep2Eta->Scale(MCscale);
    }

  std::cout << hData->Integral() << std::endl;
  std::cout << hMC->Integral() << std::endl;
  std::cout << hZmm->Integral() << std::endl;
  std::cout << hEWK->Integral() << std::endl;
  std::cout << hTop->Integral() << std::endl;

  for(int j=0;j!=hDataZPt->GetNbinsX();++j)
    {
      hDataZPt->SetBinContent(j+1,hDataZPt->GetBinContent(j+1)/hDataZPt->GetBinWidth(j+1));
      hDataZPt->SetBinError(j+1,hDataZPt->GetBinError(j+1)/hDataZPt->GetBinWidth(j+1));
      hMCZPt->SetBinContent(j+1,hMCZPt->GetBinContent(j+1)/hMCZPt->GetBinWidth(j+1));
      hMCZPt->SetBinError(j+1,hMCZPt->GetBinError(j+1)/hMCZPt->GetBinWidth(j+1));
      hZmmZPt->SetBinContent(j+1,hZmmZPt->GetBinContent(j+1)/hZmmZPt->GetBinWidth(j+1));
      hZmmZPt->SetBinError(j+1,hZmmZPt->GetBinError(j+1)/hZmmZPt->GetBinWidth(j+1));
      hEWKZPt->SetBinContent(j+1,hEWKZPt->GetBinContent(j+1)/hEWKZPt->GetBinWidth(j+1));
      hEWKZPt->SetBinError(j+1,hEWKZPt->GetBinError(j+1)/hEWKZPt->GetBinWidth(j+1));
      hTopZPt->SetBinContent(j+1,hTopZPt->GetBinContent(j+1)/hTopZPt->GetBinWidth(j+1));
      hTopZPt->SetBinError(j+1,hTopZPt->GetBinError(j+1)/hTopZPt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataPhiStar->GetNbinsX();++j)
    {
      hDataPhiStar->SetBinContent(j+1,hDataPhiStar->GetBinContent(j+1)/hDataPhiStar->GetBinWidth(j+1));
      hDataPhiStar->SetBinError(j+1,hDataPhiStar->GetBinError(j+1)/hDataPhiStar->GetBinWidth(j+1));
      hMCPhiStar->SetBinContent(j+1,hMCPhiStar->GetBinContent(j+1)/hMCPhiStar->GetBinWidth(j+1));
      hMCPhiStar->SetBinError(j+1,hMCPhiStar->GetBinError(j+1)/hMCPhiStar->GetBinWidth(j+1));
      hZmmPhiStar->SetBinContent(j+1,hZmmPhiStar->GetBinContent(j+1)/hZmmPhiStar->GetBinWidth(j+1));
      hZmmPhiStar->SetBinError(j+1,hZmmPhiStar->GetBinError(j+1)/hZmmPhiStar->GetBinWidth(j+1));
      hEWKPhiStar->SetBinContent(j+1,hEWKPhiStar->GetBinContent(j+1)/hEWKPhiStar->GetBinWidth(j+1));
      hEWKPhiStar->SetBinError(j+1,hEWKPhiStar->GetBinError(j+1)/hEWKPhiStar->GetBinWidth(j+1));
      hTopPhiStar->SetBinContent(j+1,hTopPhiStar->GetBinContent(j+1)/hTopPhiStar->GetBinWidth(j+1));
      hTopPhiStar->SetBinError(j+1,hTopPhiStar->GetBinError(j+1)/hTopPhiStar->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLep1Pt->GetNbinsX();++j)
    {
      hDataLep1Pt->SetBinContent(j+1,hDataLep1Pt->GetBinContent(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hDataLep1Pt->SetBinError(j+1,hDataLep1Pt->GetBinError(j+1)/hDataLep1Pt->GetBinWidth(j+1));
      hMCLep1Pt->SetBinContent(j+1,hMCLep1Pt->GetBinContent(j+1)/hMCLep1Pt->GetBinWidth(j+1));
      hMCLep1Pt->SetBinError(j+1,hMCLep1Pt->GetBinError(j+1)/hMCLep1Pt->GetBinWidth(j+1));
      hZmmLep1Pt->SetBinContent(j+1,hZmmLep1Pt->GetBinContent(j+1)/hZmmLep1Pt->GetBinWidth(j+1));
      hZmmLep1Pt->SetBinError(j+1,hZmmLep1Pt->GetBinError(j+1)/hZmmLep1Pt->GetBinWidth(j+1));
      hEWKLep1Pt->SetBinContent(j+1,hEWKLep1Pt->GetBinContent(j+1)/hEWKLep1Pt->GetBinWidth(j+1));
      hEWKLep1Pt->SetBinError(j+1,hEWKLep1Pt->GetBinError(j+1)/hEWKLep1Pt->GetBinWidth(j+1));
      hTopLep1Pt->SetBinContent(j+1,hTopLep1Pt->GetBinContent(j+1)/hTopLep1Pt->GetBinWidth(j+1));
      hTopLep1Pt->SetBinError(j+1,hTopLep1Pt->GetBinError(j+1)/hTopLep1Pt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLep2Pt->GetNbinsX();++j)
    {
      hDataLep2Pt->SetBinContent(j+1,hDataLep2Pt->GetBinContent(j+1)/hDataLep2Pt->GetBinWidth(j+1));
      hDataLep2Pt->SetBinError(j+1,hDataLep2Pt->GetBinError(j+1)/hDataLep2Pt->GetBinWidth(j+1));
      hMCLep2Pt->SetBinContent(j+1,hMCLep2Pt->GetBinContent(j+1)/hMCLep2Pt->GetBinWidth(j+1));
      hMCLep2Pt->SetBinError(j+1,hMCLep2Pt->GetBinError(j+1)/hMCLep2Pt->GetBinWidth(j+1));
      hZmmLep2Pt->SetBinContent(j+1,hZmmLep2Pt->GetBinContent(j+1)/hZmmLep2Pt->GetBinWidth(j+1));
      hZmmLep2Pt->SetBinError(j+1,hZmmLep2Pt->GetBinError(j+1)/hZmmLep2Pt->GetBinWidth(j+1));
      hEWKLep2Pt->SetBinContent(j+1,hEWKLep2Pt->GetBinContent(j+1)/hEWKLep2Pt->GetBinWidth(j+1));
      hEWKLep2Pt->SetBinError(j+1,hEWKLep2Pt->GetBinError(j+1)/hEWKLep2Pt->GetBinWidth(j+1));
      hTopLep2Pt->SetBinContent(j+1,hTopLep2Pt->GetBinContent(j+1)/hTopLep2Pt->GetBinWidth(j+1));
      hTopLep2Pt->SetBinError(j+1,hTopLep2Pt->GetBinError(j+1)/hTopLep2Pt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLepPosPt->GetNbinsX();++j)
    {
      hDataLepPosPt->SetBinContent(j+1,hDataLepPosPt->GetBinContent(j+1)/hDataLepPosPt->GetBinWidth(j+1));
      hDataLepPosPt->SetBinError(j+1,hDataLepPosPt->GetBinError(j+1)/hDataLepPosPt->GetBinWidth(j+1));
      hMCLepPosPt->SetBinContent(j+1,hMCLepPosPt->GetBinContent(j+1)/hMCLepPosPt->GetBinWidth(j+1));
      hMCLepPosPt->SetBinError(j+1,hMCLepPosPt->GetBinError(j+1)/hMCLepPosPt->GetBinWidth(j+1));
      hZmmLepPosPt->SetBinContent(j+1,hZmmLepPosPt->GetBinContent(j+1)/hZmmLepPosPt->GetBinWidth(j+1));
      hZmmLepPosPt->SetBinError(j+1,hZmmLepPosPt->GetBinError(j+1)/hZmmLepPosPt->GetBinWidth(j+1));
      hEWKLepPosPt->SetBinContent(j+1,hEWKLepPosPt->GetBinContent(j+1)/hEWKLepPosPt->GetBinWidth(j+1));
      hEWKLepPosPt->SetBinError(j+1,hEWKLepPosPt->GetBinError(j+1)/hEWKLepPosPt->GetBinWidth(j+1));
      hTopLepPosPt->SetBinContent(j+1,hTopLepPosPt->GetBinContent(j+1)/hTopLepPosPt->GetBinWidth(j+1));
      hTopLepPosPt->SetBinError(j+1,hTopLepPosPt->GetBinError(j+1)/hTopLepPosPt->GetBinWidth(j+1));
    }

  for(int j=0;j!=hDataLepPosPt->GetNbinsX();++j)
    {
      hDataLepPosPt->SetBinContent(j+1,hDataLepPosPt->GetBinContent(j+1)/hDataLepPosPt->GetBinWidth(j+1));
      hDataLepPosPt->SetBinError(j+1,hDataLepPosPt->GetBinError(j+1)/hDataLepPosPt->GetBinWidth(j+1));
      hMCLepPosPt->SetBinContent(j+1,hMCLepPosPt->GetBinContent(j+1)/hMCLepPosPt->GetBinWidth(j+1));
      hMCLepPosPt->SetBinError(j+1,hMCLepPosPt->GetBinError(j+1)/hMCLepPosPt->GetBinWidth(j+1));
      hZmmLepPosPt->SetBinContent(j+1,hZmmLepPosPt->GetBinContent(j+1)/hZmmLepPosPt->GetBinWidth(j+1));
      hZmmLepPosPt->SetBinError(j+1,hZmmLepPosPt->GetBinError(j+1)/hZmmLepPosPt->GetBinWidth(j+1));
      hEWKLepPosPt->SetBinContent(j+1,hEWKLepPosPt->GetBinContent(j+1)/hEWKLepPosPt->GetBinWidth(j+1));
      hEWKLepPosPt->SetBinError(j+1,hEWKLepPosPt->GetBinError(j+1)/hEWKLepPosPt->GetBinWidth(j+1));
      hTopLepPosPt->SetBinContent(j+1,hTopLepPosPt->GetBinContent(j+1)/hTopLepPosPt->GetBinWidth(j+1));
      hTopLepPosPt->SetBinError(j+1,hTopLepPosPt->GetBinError(j+1)/hTopLepPosPt->GetBinWidth(j+1));
    }
    
  TH1D *hZmumuDiff = makeDiffHist(hData,hMC,"hZmumuDiff");
  hZmumuDiff->SetMarkerStyle(kFullCircle); 
  hZmumuDiff->SetMarkerSize(0.9);

  TH1D *hZmumuNPVDiff = makeDiffHist(hDataNPV,hMCNPV,"hZmumuNPVDiff");
  hZmumuNPVDiff->SetMarkerStyle(kFullCircle); 
  hZmumuNPVDiff->SetMarkerSize(0.9);

  TH1D *hZmumuZPtDiff = makeDiffHist(hDataZPt,hMCZPt,"hZmumuZPtDiff");
  hZmumuZPtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuZPtDiff->SetMarkerSize(0.9);

  TH1D *hZmumuPhiStarDiff = makeDiffHist(hDataPhiStar,hMCPhiStar,"hZmumuPhiStarDiff");
  hZmumuPhiStarDiff->SetMarkerStyle(kFullCircle); 
  hZmumuPhiStarDiff->SetMarkerSize(0.9);

  TH1D *hZmumuZRapDiff = makeDiffHist(hDataZRap,hMCZRap,"hZmumuZRapDiff");
  hZmumuZRapDiff->SetMarkerStyle(kFullCircle); 
  hZmumuZRapDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLep1PtDiff = makeDiffHist(hDataLep1Pt,hMCLep1Pt,"hZmumuLep1PtDiff");
  hZmumuLep1PtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLep1PtDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLep2PtDiff = makeDiffHist(hDataLep2Pt,hMCLep2Pt,"hZmumuLep2PtDiff");
  hZmumuLep2PtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLep2PtDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLepNegPtDiff = makeDiffHist(hDataLepNegPt,hMCLepNegPt,"hZmumuLepNegPtDiff");
  hZmumuLepNegPtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLepNegPtDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLepPosPtDiff = makeDiffHist(hDataLepPosPt,hMCLepPosPt,"hZmumuLepPosPtDiff");
  hZmumuLepPosPtDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLepPosPtDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLep1EtaDiff = makeDiffHist(hDataLep1Eta,hMCLep1Eta,"hZmumuLep1EtaDiff");
  hZmumuLep1EtaDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLep1EtaDiff->SetMarkerSize(0.9);

  TH1D *hZmumuLep2EtaDiff = makeDiffHist(hDataLep2Eta,hMCLep2Eta,"hZmumuLep2EtaDiff");
  hZmumuLep2EtaDiff->SetMarkerStyle(kFullCircle); 
  hZmumuLep2EtaDiff->SetMarkerSize(0.9);
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  char normtext[100];
  sprintf(normtext,"MC normalized to data (#times %.2f)",MCscale);  

  string norm="";
  if(normToData)norm="_norm";
  
  // plot colors
  Int_t linecolorZ   = kOrange-3;
  Int_t fillcolorZ   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorTop = kGreen+2;
  Int_t fillcolorTop = kGreen-5;
  Int_t ratioColor   = kGray+2;

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1); 
  TGaxis::SetMaxDigits(3);
  
  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plotZmumu("zmm"+norm,"","",ylabel);
  plotZmumu.AddHist1D(hData,"data","E");
  plotZmumu.AddToStack(hZmm,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumu.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumu.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumu.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumu.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZmumu.TransLegend(0.1,-0.05);
  plotZmumu.Draw(c,kFALSE,format,1);

  CPlot plotZmumuDiff("zmm"+norm,"","M(#mu^{+}#mu^{-}) [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuDiff.AddHist1D(hZmumuDiff,"EX0",ratioColor);
  plotZmumuDiff.SetYRange(-0.2,0.2);
  plotZmumuDiff.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZmumuDiff.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZmumuDiff.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZmumuDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumu2("zmmlog"+norm,"","",ylabel);
  plotZmumu2.AddHist1D(hData,"data","E");
  plotZmumu2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumu2.AddToStack(hTop,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumu2.AddToStack(hZmm,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumu2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumu2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumu2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumu2.SetLogy();
  plotZmumu2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZmumu2.TransLegend(0.1,-0.05);
  plotZmumu2.Draw(c,kTRUE,format,1);

  //
  // NPV
  // 

  sprintf(ylabel,"Events");
  CPlot plotZmumuNPV("zmmNPV"+norm,"","",ylabel);
  plotZmumuNPV.AddHist1D(hDataNPV,"data","E");
  plotZmumuNPV.AddToStack(hZmmNPV,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuNPV.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuNPV.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuNPV.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuNPV.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZmumuNPV.TransLegend(0.1,-0.05);
  plotZmumuNPV.Draw(c,kFALSE,format,1);

  CPlot plotZmumuNPVDiff("zmmNPV"+norm,"","npv","#frac{Data-Pred}{Data}");
  plotZmumuNPVDiff.AddHist1D(hZmumuNPVDiff,"EX0",ratioColor);
  plotZmumuNPVDiff.SetYRange(-0.2,0.2);
  plotZmumuNPVDiff.AddLine(0, 0,50, 0,kBlack,1);
  plotZmumuNPVDiff.AddLine(0, 0.1,50, 0.1,kBlack,3);
  plotZmumuNPVDiff.AddLine(0,-0.1,50,-0.1,kBlack,3);
  plotZmumuNPVDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuNPV2("zmmNPVlog"+norm,"","",ylabel);
  plotZmumuNPV2.AddHist1D(hDataNPV,"data","E");
  plotZmumuNPV2.AddToStack(hEWKNPV,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuNPV2.AddToStack(hTopNPV,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuNPV2.AddToStack(hZmmNPV,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuNPV2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuNPV2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuNPV2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuNPV2.SetLogy();
  plotZmumuNPV2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZmumuNPV2.TransLegend(0.1,-0.05);
  plotZmumuNPV2.Draw(c,kTRUE,format,1);

  //
  // Z Pt
  //   
  sprintf(ylabel,"Events / 1 GeV");
  CPlot plotZmumuZPt("zmmZPt"+norm,"","",ylabel);
  plotZmumuZPt.AddHist1D(hDataZPt,"data","E");
  plotZmumuZPt.AddToStack(hZmmZPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuZPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuZPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuZPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuZPt.SetLogx();
  plotZmumuZPt.SetLogy(0);
  plotZmumuZPt.SetYRange(0.01,1.2*(hDataZPt->GetMaximum() + sqrt(hDataZPt->GetMaximum())));
  plotZmumuZPt.TransLegend(0.1,-0.05);
  plotZmumuZPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuZPtDiff("zmmZPt"+norm,"","p_{T}^{#mu^{+}#mu^{-}} [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuZPtDiff.AddHist1D(hZmumuZPtDiff,"EX0",ratioColor);
  plotZmumuZPtDiff.SetLogx();
  plotZmumuZPtDiff.SetYRange(-0.2,0.2);
  plotZmumuZPtDiff.AddLine(0, 0,1000, 0,kBlack,1);
  plotZmumuZPtDiff.AddLine(0, 0.1,1000, 0.1,kBlack,3);
  plotZmumuZPtDiff.AddLine(0,-0.1,1000,-0.1,kBlack,3);
  plotZmumuZPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuZPt2("zmmZPtlog"+norm,"","",ylabel);
  plotZmumuZPt2.AddHist1D(hDataZPt,"data","E");
  plotZmumuZPt2.AddToStack(hEWKZPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuZPt2.AddToStack(hTopZPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuZPt2.AddToStack(hZmmZPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);\
  plotZmumuZPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuZPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuZPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuZPt2.SetLogx();
  plotZmumuZPt2.SetLogy();
  plotZmumuZPt2.SetYRange(1e-6*(hDataZPt->GetMaximum()),10*(hDataZPt->GetMaximum()));
  plotZmumuZPt2.TransLegend(0.1,-0.05);
  plotZmumuZPt2.Draw(c,kTRUE,format,1);

  //
  // Phi*
  //   
  sprintf(ylabel,"Events / 1.0");
  CPlot plotZmumuPhiStar("zmmPhiStar"+norm,"","",ylabel);
  plotZmumuPhiStar.AddHist1D(hDataPhiStar,"data","E");
  plotZmumuPhiStar.AddToStack(hZmmPhiStar,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuPhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuPhiStar.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuPhiStar.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuPhiStar.SetLogx();
  plotZmumuPhiStar.SetLogy(0);
  plotZmumuPhiStar.SetYRange(0.01,1.2*(hDataPhiStar->GetMaximum() + sqrt(hDataPhiStar->GetMaximum())));
  plotZmumuPhiStar.TransLegend(0.1,-0.05);
  plotZmumuPhiStar.Draw(c,kFALSE,format,1);

  CPlot plotZmumuPhiStarDiff("zmmPhiStar"+norm,"","#phi_{#eta}*","#frac{Data-Pred}{Data}");
  plotZmumuPhiStarDiff.AddHist1D(hZmumuPhiStarDiff,"EX0",ratioColor);
  plotZmumuPhiStarDiff.SetLogx();
  plotZmumuPhiStarDiff.SetYRange(-0.2,0.2);
  plotZmumuPhiStarDiff.AddLine(0, 0,1.1, 0,kBlack,1);
  plotZmumuPhiStarDiff.AddLine(0, 0.1,1.1, 0.1,kBlack,3);
  plotZmumuPhiStarDiff.AddLine(0,-0.1,1.1,-0.1,kBlack,3);
  plotZmumuPhiStarDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuPhiStar2("zmmPhiStarlog"+norm,"","",ylabel);
  plotZmumuPhiStar2.AddHist1D(hDataPhiStar,"data","E");
  plotZmumuPhiStar2.AddToStack(hEWKPhiStar,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuPhiStar2.AddToStack(hTopPhiStar,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuPhiStar2.AddToStack(hZmmPhiStar,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuPhiStar2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuPhiStar2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuPhiStar2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZmumuPhiStar2.SetLogx();
  plotZmumuPhiStar2.SetLogy();
  plotZmumuPhiStar2.SetYRange(1e-5*(hDataPhiStar->GetMaximum()),10*(hDataPhiStar->GetMaximum()));
  plotZmumuPhiStar2.TransLegend(0.1,-0.05);
  plotZmumuPhiStar2.Draw(c,kTRUE,format,1);

  //
  // Z Rapidity
  //   
  sprintf(ylabel,"Events / %.1f ",hDataZRap->GetBinWidth(1));
  CPlot plotZmumuZRap("zmmZRap"+norm,"","",ylabel);
  plotZmumuZRap.AddHist1D(hDataZRap,"data","E");
  plotZmumuZRap.AddToStack(hZmmZRap,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuZRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuZRap.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuZRap.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZmumuZRap.SetLogy(0);
  plotZmumuZRap.SetYRange(0.01,1.25*(hDataZRap->GetMaximum() + sqrt(hDataZRap->GetMaximum())));
  plotZmumuZRap.TransLegend(0.1,-0.05);
  plotZmumuZRap.Draw(c,kFALSE,format,1);

  CPlot plotZmumuZRapDiff("zmmZRap"+norm,"","|y^{#mu^{+}#mu^{-}}|","#frac{Data-Pred}{Data}");
  plotZmumuZRapDiff.AddHist1D(hZmumuZRapDiff,"EX0",ratioColor);
  plotZmumuZRapDiff.SetYRange(-0.2,0.2);
  plotZmumuZRapDiff.AddLine(0, 0,2.4, 0,kBlack,1);
  plotZmumuZRapDiff.AddLine(0, 0.1,2.4, 0.1,kBlack,3);
  plotZmumuZRapDiff.AddLine(0,-0.1,2.4,-0.1,kBlack,3);
  plotZmumuZRapDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuZRap2("zmmZRaplog"+norm,"","",ylabel);
  plotZmumuZRap2.AddHist1D(hDataZRap,"data","E");
  plotZmumuZRap2.AddToStack(hEWKZRap,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuZRap2.AddToStack(hTopZRap,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuZRap2.AddToStack(hZmmZRap,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuZRap2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuZRap2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuZRap2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZmumuZRap2.SetLogy();
  plotZmumuZRap2.SetYRange(1e-4*(hDataZRap->GetMaximum()),500*(hDataZRap->GetMaximum()));
  plotZmumuZRap2.TransLegend(0.1,-0.05);
  plotZmumuZRap2.Draw(c,kTRUE,format,1);

 //
  // Lep1 Pt
  //   
  sprintf(ylabel,"Events / %.1f GeV",hDataLep1Pt->GetBinWidth(1));
  CPlot plotZmumuLep1Pt("zmmLep1Pt"+norm,"","",ylabel);
  plotZmumuLep1Pt.AddHist1D(hDataLep1Pt,"data","E");
  plotZmumuLep1Pt.AddToStack(hZmmLep1Pt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep1Pt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLep1Pt.SetLogx();
  plotZmumuLep1Pt.SetLogy(0);
  plotZmumuLep1Pt.SetYRange(0.01,1.2*(hDataLep1Pt->GetMaximum() + sqrt(hDataLep1Pt->GetMaximum())));
  plotZmumuLep1Pt.TransLegend(0.1,-0.05);
  plotZmumuLep1Pt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep1PtDiff("zmmLep1Pt"+norm,"","p_{T}(leading muon) [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuLep1PtDiff.AddHist1D(hZmumuLep1PtDiff,"EX0",ratioColor);
  plotZmumuLep1PtDiff.SetLogx();
  plotZmumuLep1PtDiff.SetYRange(-0.2,0.2);
  plotZmumuLep1PtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZmumuLep1PtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZmumuLep1PtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZmumuLep1PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep1Pt2("zmmLep1Ptlog"+norm,"","",ylabel);
  plotZmumuLep1Pt2.AddHist1D(hDataLep1Pt,"data","E");
  plotZmumuLep1Pt2.AddToStack(hEWKLep1Pt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLep1Pt2.AddToStack(hTopLep1Pt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLep1Pt2.AddToStack(hZmmLep1Pt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep1Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep1Pt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuLep1Pt2.SetLogx();
  plotZmumuLep1Pt2.SetLogy();
  plotZmumuLep1Pt2.SetYRange(1e-5*(hDataLep1Pt->GetMaximum()),10*(hDataLep1Pt->GetMaximum()));
  plotZmumuLep1Pt2.TransLegend(0.1,-0.05);
  plotZmumuLep1Pt2.Draw(c,kTRUE,format,1);

  //
  // Lep2 Pt
  //   
  sprintf(ylabel,"Events / %.1f GeV",hDataLep2Pt->GetBinWidth(1));
  CPlot plotZmumuLep2Pt("zmmLep2Pt"+norm,"","",ylabel);
  plotZmumuLep2Pt.AddHist1D(hDataLep2Pt,"data","E");
  plotZmumuLep2Pt.AddToStack(hZmmLep2Pt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep2Pt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLep2Pt.SetLogx();
  plotZmumuLep2Pt.SetLogy(0);
  plotZmumuLep2Pt.SetYRange(0.01,1.2*(hDataLep2Pt->GetMaximum() + sqrt(hDataLep2Pt->GetMaximum())));
  plotZmumuLep2Pt.TransLegend(0.1,-0.05);
  plotZmumuLep2Pt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep2PtDiff("zmmLep2Pt"+norm,"","p_{T}(2nd leading muon) [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuLep2PtDiff.AddHist1D(hZmumuLep2PtDiff,"EX0",ratioColor);
  plotZmumuLep2PtDiff.SetLogx();
  plotZmumuLep2PtDiff.SetYRange(-0.2,0.2);
  plotZmumuLep2PtDiff.AddLine(25, 0,200, 0,kBlack,1);
  plotZmumuLep2PtDiff.AddLine(25, 0.1,200, 0.1,kBlack,3);
  plotZmumuLep2PtDiff.AddLine(25,-0.1,200,-0.1,kBlack,3);
  plotZmumuLep2PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep2Pt2("zmmLep2Ptlog"+norm,"","",ylabel);
  plotZmumuLep2Pt2.AddHist1D(hDataLep2Pt,"data","E");
  plotZmumuLep2Pt2.AddToStack(hEWKLep2Pt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLep2Pt2.AddToStack(hTopLep2Pt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLep2Pt2.AddToStack(hZmmLep2Pt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);\
  plotZmumuLep2Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep2Pt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuLep2Pt2.SetLogx();
  plotZmumuLep2Pt2.SetLogy();
  plotZmumuLep2Pt2.SetYRange(1e-5*(hDataLep2Pt->GetMaximum()),10*(hDataLep2Pt->GetMaximum()));
  plotZmumuLep2Pt2.TransLegend(0.1,-0.05);
  plotZmumuLep2Pt2.Draw(c,kTRUE,format,1);

  //
  // LepNeg Pt
  //   
  sprintf(ylabel,"Events / %.1f GeV",hDataLepNegPt->GetBinWidth(1));
  CPlot plotZmumuLepNegPt("zmmLepNegPt"+norm,"","",ylabel);
  plotZmumuLepNegPt.AddHist1D(hDataLepNegPt,"data","E");
  plotZmumuLepNegPt.AddToStack(hZmmLepNegPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepNegPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLepNegPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLepNegPt.SetLogx();
  plotZmumuLepNegPt.SetLogy(0);
  plotZmumuLepNegPt.SetYRange(0.01,1.2*(hDataLepNegPt->GetMaximum() + sqrt(hDataLepNegPt->GetMaximum())));
  plotZmumuLepNegPt.TransLegend(0.1,-0.05);
  plotZmumuLepNegPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLepNegPtDiff("zmmLepNegPt"+norm,"","p_{T}^{#mu^{-}} [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuLepNegPtDiff.AddHist1D(hZmumuLepNegPtDiff,"EX0",ratioColor);
  plotZmumuLepNegPtDiff.SetLogx();
  plotZmumuLepNegPtDiff.SetYRange(-0.2,0.2);
  plotZmumuLepNegPtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZmumuLepNegPtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZmumuLepNegPtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZmumuLepNegPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLepNegPt2("zmmLepNegPtlog"+norm,"","",ylabel);
  plotZmumuLepNegPt2.AddHist1D(hDataLepNegPt,"data","E");
  plotZmumuLepNegPt2.AddToStack(hEWKLepNegPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLepNegPt2.AddToStack(hTopLepNegPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLepNegPt2.AddToStack(hZmmLepNegPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLepNegPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepNegPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLepNegPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuLepNegPt2.SetLogx();
  plotZmumuLepNegPt2.SetLogy();
  plotZmumuLepNegPt2.SetYRange(1e-5*(hDataLepNegPt->GetMaximum()),10*(hDataLepNegPt->GetMaximum()));
  plotZmumuLepNegPt2.TransLegend(0.1,-0.05);
  plotZmumuLepNegPt2.Draw(c,kTRUE,format,1);

  //
  // LepPos Pt
  //   
  sprintf(ylabel,"Events / %.1f GeV",hDataLepPosPt->GetBinWidth(1));
  CPlot plotZmumuLepPosPt("zmmLepPosPt"+norm,"","",ylabel);
  plotZmumuLepPosPt.AddHist1D(hDataLepPosPt,"data","E");
  plotZmumuLepPosPt.AddToStack(hZmmLepPosPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepPosPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLepPosPt.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLepPosPt.SetLogx();
  plotZmumuLepPosPt.SetLogy(0);
  plotZmumuLepPosPt.SetYRange(0.01,1.2*(hDataLepPosPt->GetMaximum() + sqrt(hDataLepPosPt->GetMaximum())));
  plotZmumuLepPosPt.TransLegend(0.1,-0.05);
  plotZmumuLepPosPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLepPosPtDiff("zmmLepPosPt"+norm,"","p_{T}^{#mu^{+}} [GeV]","#frac{Data-Pred}{Data}");
  plotZmumuLepPosPtDiff.AddHist1D(hZmumuLepPosPtDiff,"EX0",ratioColor);
  plotZmumuLepPosPtDiff.SetLogx();
  plotZmumuLepPosPtDiff.SetYRange(-0.2,0.2);
  plotZmumuLepPosPtDiff.AddLine(25, 0,300, 0,kBlack,1);
  plotZmumuLepPosPtDiff.AddLine(25, 0.1,300, 0.1,kBlack,3);
  plotZmumuLepPosPtDiff.AddLine(25,-0.1,300,-0.1,kBlack,3);
  plotZmumuLepPosPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLepPosPt2("zmmLepPosPtlog"+norm,"","",ylabel);
  plotZmumuLepPosPt2.AddHist1D(hDataLepPosPt,"data","E");
  plotZmumuLepPosPt2.AddToStack(hEWKLepPosPt,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLepPosPt2.AddToStack(hTopLepPosPt,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLepPosPt2.AddToStack(hZmmLepPosPt,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);\
  plotZmumuLepPosPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepPosPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLepPosPt2.AddTextBox(normtext,0.6,0.55,0.95,0.60,0,13,0.03,-1);
  plotZmumuLepPosPt2.SetLogx();
  plotZmumuLepPosPt2.SetLogy();
  plotZmumuLepPosPt2.SetYRange(1e-5*(hDataLepPosPt->GetMaximum()),10*(hDataLepPosPt->GetMaximum()));
  plotZmumuLepPosPt2.TransLegend(0.1,-0.05);
  plotZmumuLepPosPt2.Draw(c,kTRUE,format,1);

  //
  // Lep1 Eta
  //   
  sprintf(ylabel,"Events / %.1f ",hDataLep1Eta->GetBinWidth(1));
  CPlot plotZmumuLep1Eta("zmmLep1Eta"+norm,"","",ylabel);
  plotZmumuLep1Eta.AddHist1D(hDataLep1Eta,"data","E");
  plotZmumuLep1Eta.AddToStack(hZmmLep1Eta,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Eta.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep1Eta.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLep1Eta.SetLogy(0);
  plotZmumuLep1Eta.SetYRange(0.01,1.2*(hDataLep1Eta->GetMaximum() + sqrt(hDataLep1Eta->GetMaximum())));
  plotZmumuLep1Eta.TransLegend(0.1,-0.05);
  plotZmumuLep1Eta.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep1EtaDiff("zmmLep1Eta"+norm,"","|#eta| (leading muon)","#frac{Data-Pred}{Data}");
  plotZmumuLep1EtaDiff.AddHist1D(hZmumuLep1EtaDiff,"EX0",ratioColor);
  plotZmumuLep1EtaDiff.SetYRange(-0.2,0.2);
  plotZmumuLep1EtaDiff.AddLine(0, 0,2.4, 0,kBlack,1);
  plotZmumuLep1EtaDiff.AddLine(0, 0.1,2.4, 0.1,kBlack,3);
  plotZmumuLep1EtaDiff.AddLine(0,-0.1,2.4,-0.1,kBlack,3);
  plotZmumuLep1EtaDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep1Eta2("zmmLep1Etalog"+norm,"","",ylabel);
  plotZmumuLep1Eta2.AddHist1D(hDataLep1Eta,"data","E");
  plotZmumuLep1Eta2.AddToStack(hEWKLep1Eta,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLep1Eta2.AddToStack(hTopLep1Eta,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLep1Eta2.AddToStack(hZmmLep1Eta,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep1Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Eta2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep1Eta2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZmumuLep1Eta2.SetLogy();
  plotZmumuLep1Eta2.SetYRange(1e-4*(hDataLep1Eta->GetMaximum()),500*(hDataLep1Eta->GetMaximum()));
  plotZmumuLep1Eta2.TransLegend(0.1,-0.05);
  plotZmumuLep1Eta2.Draw(c,kTRUE,format,1);

  //
  // Lep2 Eta
  //   
  sprintf(ylabel,"Events / %.1f ",hDataLep2Eta->GetBinWidth(1));
  CPlot plotZmumuLep2Eta("zmmLep2Eta"+norm,"","",ylabel);
  plotZmumuLep2Eta.AddHist1D(hDataLep2Eta,"data","E");
  plotZmumuLep2Eta.AddToStack(hZmmLep2Eta,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Eta.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep2Eta.AddTextBox(normtext,0.6,0.65,0.95,0.70,0,13,0.03,-1);
  plotZmumuLep2Eta.SetLogy(0);
  plotZmumuLep2Eta.SetYRange(0.01,1.2*(hDataLep2Eta->GetMaximum() + sqrt(hDataLep2Eta->GetMaximum())));
  plotZmumuLep2Eta.TransLegend(0.1,-0.05);
  plotZmumuLep2Eta.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep2EtaDiff("zmmLep2Eta"+norm,"","|#eta| (2nd leading muon)","#frac{Data-Pred}{Data}");
  plotZmumuLep2EtaDiff.AddHist1D(hZmumuLep2EtaDiff,"EX0",ratioColor);
  plotZmumuLep2EtaDiff.SetYRange(-0.2,0.2);
  plotZmumuLep2EtaDiff.AddLine(0, 0,2.4, 0,kBlack,1);
  plotZmumuLep2EtaDiff.AddLine(0, 0.1,2.4, 0.1,kBlack,3);
  plotZmumuLep2EtaDiff.AddLine(0,-0.1,2.4,-0.1,kBlack,3);
  plotZmumuLep2EtaDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep2Eta2("zmmLep2Etalog"+norm,"","",ylabel);
  plotZmumuLep2Eta2.AddHist1D(hDataLep2Eta,"data","E");
  plotZmumuLep2Eta2.AddToStack(hEWKLep2Eta,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumuLep2Eta2.AddToStack(hTopLep2Eta,"t#bar{t}",fillcolorTop,linecolorTop);
  plotZmumuLep2Eta2.AddToStack(hZmmLep2Eta,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumuLep2Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Eta2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  if(normToData)plotZmumuLep2Eta2.AddTextBox(normtext,0.21,0.75,0.45,0.80,0,13,0.03,-1);
  plotZmumuLep2Eta2.SetLogy();
  plotZmumuLep2Eta2.SetYRange(1e-4*(hDataLep2Eta->GetMaximum()),500*(hDataLep2Eta->GetMaximum()));
  plotZmumuLep2Eta2.TransLegend(0.1,-0.05);
  plotZmumuLep2Eta2.Draw(c,kTRUE,format,1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << " The Zmm event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  cout << " The Zmm expected event yield is " << yield_zmm << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  cout << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     

  gBenchmark->Show("plotZmm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = (TH1D*)hData->Clone("hDiff");
  hDiff->SetName(name);
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff=0;
    Double_t err=0;
    if(hData->GetBinContent(ibin)!=0)
      {
	diff = diff0/hData->GetBinContent(ibin);
	err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
      }
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  
  return hDiff;
}

