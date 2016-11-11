//================================================================================================
//
// Perform fit to extract W->enu signal
//
//  * outputs plots and fit results summary
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
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_asym2.hh"    // class to handle recoil corrections for MET
//#include "../Utils/RecoilCorrector.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
// #include "ZBackgrounds.hh"
#include "RooCategory.h"

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe);

// make webpage
void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void fitWe(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
           const Double_t lumi2,
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("fitWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // some flags to handle Recoil corrections
  bool doKeys = false;
  bool doInclusive = false;
  // some flags to handle the pileup Up/Down systematics
  bool pileupUp = false;
  bool pileupDown = false;
  
  
  vector<Float_t> upper; 
  upper.push_back(2.5); upper.push_back(2.3); upper.push_back(2.1); 
  upper.push_back(1.9); upper.push_back(1.7); upper.push_back(1.5);
  upper.push_back(1.2); upper.push_back(0.8); upper.push_back(0.4);
  //upper.push_back(2.5); upper.push_back(2.1); upper.push_back(1.5);upper.push_back(0);upper.push_back(-1.5); upper.push_back(-2.1);
  vector<Float_t> lower; 
  lower.push_back(2.3);  lower.push_back(2.1);  lower.push_back(1.9);
  lower.push_back(1.7);  lower.push_back(1.5);  lower.push_back(1.2);
  lower.push_back(0.8); lower.push_back(0.4);   lower.push_back(0);
  //lower.push_back(2.1); lower.push_back(1.5);lower.push_back(0);lower.push_back(-1.5); lower.push_back(-2.1); lower.push_back(-2.5); 
  vector<Float_t> fpcorr;
  //fpcorr.push_back(0.00502415); fpcorr.push_back(0.0141985); fpcorr.push_back(0.0136337); fpcorr.push_back(0.0155871); fpcorr.push_back(0.0144434); fpcorr.push_back(0.000100203);
  fpcorr.push_back(-0.000826797); fpcorr.push_back(0.00520995); fpcorr.push_back(0.0120186); fpcorr.push_back(0.0121537); fpcorr.push_back(0.0215729); fpcorr.push_back(0.0192083); fpcorr.push_back(0.014545); fpcorr.push_back(0.0155789);  fpcorr.push_back(0.0115706);


  // MET histogram binning and range
  const Int_t    NBINS  = 75;
  const Double_t METMAX = 150;

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.5;
  
  // efficiency files

  const TString baseDir = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/";
  const TString dataHLTEffName_pos = baseDir + "EleHLTEff/MGpositive/eff.root";
  const TString dataHLTEffName_neg = baseDir + "EleHLTEff/MGnegative/eff.root";
  const TString zeeHLTEffName_pos  = baseDir + "EleHLTEff/CTpositive/eff.root";
  const TString zeeHLTEffName_neg  = baseDir + "EleHLTEff/CTnegative/eff.root";
  
  const TString dataGsfSelEffName_pos = baseDir + "EleGsfSelEff/MGpositive/eff.root";
  const TString dataGsfSelEffName_neg = baseDir + "EleGsfSelEff/MGnegative/eff.root";
  const TString zeeGsfSelEffName_pos  = baseDir + "EleGsfSelEff/CTpositive/eff.root";
  const TString zeeGsfSelEffName_neg  = baseDir + "EleGsfSelEff/CTnegative/eff.root";

  //efficiency files 2Bins

  const TString dataHLTEff2BinName_pos = baseDir + "EleHLTEff/MGpositive/eff.root";
  const TString dataHLTEff2BinName_neg = baseDir + "EleHLTEff/MGnegative/eff.root";
  const TString zeeHLTEff2BinName_pos  = baseDir + "EleHLTEff/CTpositive/eff.root";
  const TString zeeHLTEff2BinName_neg  = baseDir + "EleHLTEff/CTnegative/eff.root";
  
  const TString dataGsfSelEff2BinName_pos = baseDir + "EleGsfSelEff/MGpositive/eff.root";
  const TString dataGsfSelEff2BinName_neg = baseDir + "EleGsfSelEff/MGnegative/eff.root";
  const TString zeeGsfSelEff2BinName_pos  = baseDir + "EleGsfSelEff/CTpositive/eff.root";
  const TString zeeGsfSelEff2BinName_neg  = baseDir + "EleGsfSelEff/CTnegative/eff.root";

  TString GsfSelEffSignalShapeSys = baseDir + "Results/EleGsfSelSigSys.root";
  TString GsfSelEffBackgroundShapeSys = baseDir + "Results/EleGsfSelBkgSys.root";
  
    // setup efficiency shape systematics
  TFile *GsfSelSigSysFile = new TFile(GsfSelEffSignalShapeSys);
  TH2D *hGsfSelSigSys = (TH2D*)GsfSelSigSysFile->Get("h");
  TFile *GsfSelBkgSysFile = new TFile(GsfSelEffBackgroundShapeSys);
  TH2D *hGsfSelBkgSys = (TH2D*)GsfSelBkgSysFile->Get("h");

  // file format for output plots
  const TString format("png"); 

  // recoil correction
  //, (!) uncomment to perform corrections to recoil from W-MC/Z-MC



  // ======================= Recoil Corrections ================================
  const TString directory("/afs/cern.ch/user/d/dalfonso/public/WZ/oct7");
  // for Puppi, inclusive
  RecoilCorrector *recoilCorr = new  RecoilCorrector("","");
  recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi/",directory.Data()));
  recoilCorr->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg/",directory.Data()));
  recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("","");
  recoilCorrm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi/",directory.Data()));
  recoilCorrm->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg/",directory.Data()));
  recoilCorrm->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  
  // placeholders until recoil files are fixed
  RecoilCorrector *recoilCorrPuUp = new  RecoilCorrector("","");
  recoilCorrPuUp->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_PileupUp/",directory.Data()));
  recoilCorrPuUp->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_PileupUp/",directory.Data()));
  recoilCorrPuUp->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_PileupUp/",directory.Data()));
  
  RecoilCorrector *recoilCorrPuUpm = new  RecoilCorrector("","");
  recoilCorrPuUpm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_PileupUp/",directory.Data()));
  recoilCorrPuUpm->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_PileupUp/",directory.Data()));
  recoilCorrPuUpm->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_PileupUp/",directory.Data()));
  
  RecoilCorrector *recoilCorrPuDown = new  RecoilCorrector("","");
  recoilCorrPuDown->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_PileupDown/",directory.Data()));
  recoilCorrPuDown->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_PileupDown/",directory.Data()));
  recoilCorrPuDown->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_PileupDown/",directory.Data()));
  
  RecoilCorrector *recoilCorrPuDownm = new  RecoilCorrector("","");
  recoilCorrPuDownm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_PileupDown/",directory.Data()));
  recoilCorrPuDownm->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_PileupDown/",directory.Data()));
  recoilCorrPuDownm->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_PileupDown/",directory.Data()));

  // --------------------- Eta-binned recoil corrections -----------------------
  RecoilCorrector *recoilCorr05 = new  RecoilCorrector("","");
  recoilCorr05->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_rap05/",directory.Data()));
  recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  recoilCorr05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));
  
  RecoilCorrector *recoilCorrm05 = new  RecoilCorrector("","");
  recoilCorrm05->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_rap05/",directory.Data()));
  recoilCorrm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  recoilCorrm05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));

  RecoilCorrector *recoilCorr051 = new  RecoilCorrector("","");
  recoilCorr051->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_rap05-1/",directory.Data()));
  recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  recoilCorr051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));

  RecoilCorrector *recoilCorrm051 = new  RecoilCorrector("","");
  recoilCorrm051->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_rap05-1/",directory.Data()));
  recoilCorrm051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  recoilCorrm051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));

  RecoilCorrector *recoilCorr1 = new  RecoilCorrector("","");
  recoilCorr1->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_rap1/",directory.Data()));
  recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  recoilCorr1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data()));

  RecoilCorrector *recoilCorrm1 = new  RecoilCorrector("","");
  recoilCorrm1->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_rap1/",directory.Data()));
  recoilCorrm1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  recoilCorrm1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data())); 
  
  // ---------------------- KEYS -------------
  RecoilCorrector *recoilCorrKeys05 = new  RecoilCorrector("","");
  recoilCorrKeys05->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys_rap05/",directory.Data()));
  recoilCorrKeys05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  recoilCorrKeys05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));
  
  RecoilCorrector *recoilCorrKeysm05 = new  RecoilCorrector("","");
  recoilCorrKeysm05->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys_rap05/",directory.Data()));
  recoilCorrKeysm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  recoilCorrKeysm05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));

  RecoilCorrector *recoilCorrKeys051 = new  RecoilCorrector("","");
  recoilCorrKeys051->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys_rap05-1/",directory.Data()));
  recoilCorrKeys051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  recoilCorrKeys051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));

  RecoilCorrector *recoilCorrKeysm051 = new  RecoilCorrector("","");
  recoilCorrKeysm051->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys_rap05-1/",directory.Data()));
  recoilCorrKeysm051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  recoilCorrKeysm051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));

  RecoilCorrector *recoilCorrKeys1 = new  RecoilCorrector("","");
  recoilCorrKeys1->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys_rap1/",directory.Data()));
  recoilCorrKeys1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  recoilCorrKeys1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data()));

  RecoilCorrector *recoilCorrKeysm1 = new  RecoilCorrector("","");
  recoilCorrKeysm1->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys_rap1/",directory.Data()));
  recoilCorrKeysm1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  recoilCorrKeysm1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data())); 
  
  // ===========================================================================
  
  TFile *_rdWmp = new TFile("shapeDiff/Wep_relDiff.root");
  TFile *_rdWmm = new TFile("shapeDiff/Wem_relDiff.root");
  
  TH1D *hh_diffm = new TH1D("hh_diffm","hh_diffm",75,0,150);
  TH1D *hh_diffp = new TH1D("hh_diffp","hh_diffp",75,0,150);
  
  hh_diffm = (TH1D*)_rdWmm->Get("hh_diff");
  hh_diffp = (TH1D*)_rdWmp->Get("hh_diff");
  
  TFile *_rat1 = new TFile("shapeDiff/zmm_PDFUnc.root");
  TH1D *hh_mc;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_mc = (TH1D*)_rat1->Get("hZPtTruthNominal");
  hh_mc->Scale(1/hh_mc->Integral());
  
  TFile *_rat2 = new TFile("shapeDiff/UnfoldingOutputZPt.root");
  TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_diff = (TH1D*)_rat2->Get("hUnfold");
  hh_diff->Scale(1/hh_diff->Integral());
  hh_diff->Divide(hh_mc);
//   TCanvas *b = new TCanvas("b","b",800,800);
  //hh_diff->Smooth(100);
//   hh_diff->Draw("");
//   b->SaveAs("smoothing.png");
//   return;

  enum { eData, eWenu, eEWK , eBKG, eAntiData, eAntiWenu, eAntiEWK};  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  

  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wenu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wenu/ntuples/we_select.root");   typev.push_back(eWenu);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wenu/ntuples/ewk_select1.root");  typev.push_back(eEWK);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wenu/ntuples/boson_select.root");  typev.push_back(eBKG);
  //fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wenu/ntuples/ewk_select.root");  typev.push_back(eEWK);
  //fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wenu/ntuples/top_select.root");  typev.push_back(eEWK);

 
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWenu/ntuples/data_select.root"); typev.push_back(eAntiData);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWenu/ntuples/we_select.root");   typev.push_back(eAntiWenu);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWenu/ntuples/ewk_select.root");  typev.push_back(eAntiEWK);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWenu/ntuples/top_select.root");  typev.push_back(eAntiEWK);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet  = new TH1D("hDataMet", "",NBINS,0,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm = new TH1D("hDataMetm","",NBINS,0,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp = new TH1D("hDataMetp","",NBINS,0,METMAX); hDataMetp->Sumw2();
  TH1D *hWenuMet  = new TH1D("hWenuMet", "",NBINS,0,METMAX); hWenuMet->Sumw2();
  TH1D *hWenuMetp = new TH1D("hWenuMetp","",NBINS,0,METMAX); hWenuMetp->Sumw2();
  TH1D *hWenuMetm = new TH1D("hWenuMetm","",NBINS,0,METMAX); hWenuMetm->Sumw2();
  TH1D *hEWKMet   = new TH1D("hEWKMet",  "",NBINS,0,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp  = new TH1D("hEWKMetp", "",NBINS,0,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm  = new TH1D("hEWKMetm", "",NBINS,0,METMAX); hEWKMetm->Sumw2();
  
  TH1D *hEWKMet_PileupUp   = new TH1D("hEWKMet_PileupUp",  "",NBINS,0,METMAX); hEWKMet_PileupUp->Sumw2();
  TH1D *hEWKMetp_PileupUp  = new TH1D("hEWKMetp_PileupUp", "",NBINS,0,METMAX); hEWKMetp_PileupUp->Sumw2();
  TH1D *hEWKMetm_PileupUp  = new TH1D("hEWKMetm_PileupUp", "",NBINS,0,METMAX); hEWKMetm_PileupUp->Sumw2();
  
  TH1D *hEWKMet_PileupDown   = new TH1D("hEWKMet_PileupDown",  "",NBINS,0,METMAX); hEWKMet_PileupDown->Sumw2();
  TH1D *hEWKMetp_PileupDown  = new TH1D("hEWKMetp_PileupDown", "",NBINS,0,METMAX); hEWKMetp_PileupDown->Sumw2();
  TH1D *hEWKMetm_PileupDown  = new TH1D("hEWKMetm_PileupDown", "",NBINS,0,METMAX); hEWKMetm_PileupDown->Sumw2();
  
  // Recoil shape uncertainties
  TH1D *hWenuMet_RecoilUp  = new TH1D("hWenuMet_RecoilUp", "",NBINS,0,METMAX); hWenuMet_RecoilUp->Sumw2();
  TH1D *hWenuMetp_RecoilUp = new TH1D("hWenuMetp_RecoilUp","",NBINS,0,METMAX); hWenuMetp_RecoilUp->Sumw2();
  TH1D *hWenuMetm_RecoilUp = new TH1D("hWenuMetm_RecoilUp","",NBINS,0,METMAX); hWenuMetm_RecoilUp->Sumw2();
  TH1D *hWenuMet_RecoilDown  = new TH1D("hWenuMet_RecoilDown", "",NBINS,0,METMAX); hWenuMet_RecoilDown->Sumw2();
  TH1D *hWenuMetp_RecoilDown = new TH1D("hWenuMetp_RecoilDown","",NBINS,0,METMAX); hWenuMetp_RecoilDown->Sumw2();
  TH1D *hWenuMetm_RecoilDown = new TH1D("hWenuMetm_RecoilDown","",NBINS,0,METMAX); hWenuMetm_RecoilDown->Sumw2();
  
  // 
  TH1D *hWenuMet_ScaleUp  = new TH1D("hWenuMet_ScaleUp", "",NBINS,0,METMAX); hWenuMet_ScaleUp->Sumw2();
  TH1D *hWenuMetp_ScaleUp = new TH1D("hWenuMetp_ScaleUp","",NBINS,0,METMAX); hWenuMetp_ScaleUp->Sumw2();
  TH1D *hWenuMetm_ScaleUp = new TH1D("hWenuMetm_ScaleUp","",NBINS,0,METMAX); hWenuMetm_ScaleUp->Sumw2();
  TH1D *hWenuMet_ScaleDown  = new TH1D("hWenuMet_ScaleDown", "",NBINS,0,METMAX); hWenuMet_ScaleDown->Sumw2();
  TH1D *hWenuMetp_ScaleDown = new TH1D("hWenuMetp_ScaleDown","",NBINS,0,METMAX); hWenuMetp_ScaleDown->Sumw2();
  TH1D *hWenuMetm_ScaleDown = new TH1D("hWenuMetm_ScaleDown","",NBINS,0,METMAX); hWenuMetm_ScaleDown->Sumw2();
  
  TH1D *hWenuMet_PileupUp  = new TH1D("hWenuMet_PileupUp", "",NBINS,0,METMAX); hWenuMet_PileupUp->Sumw2();
  TH1D *hWenuMetp_PileupUp = new TH1D("hWenuMetp_PileupUp","",NBINS,0,METMAX); hWenuMetp_PileupUp->Sumw2();
  TH1D *hWenuMetm_PileupUp = new TH1D("hWenuMetm_PileupUp","",NBINS,0,METMAX); hWenuMetm_PileupUp->Sumw2();
  TH1D *hWenuMet_PileupDown  = new TH1D("hWenuMet_PileupDown", "",NBINS,0,METMAX); hWenuMet_PileupDown->Sumw2();
  TH1D *hWenuMetp_PileupDown = new TH1D("hWenuMetp_PileupDown","",NBINS,0,METMAX); hWenuMetp_PileupDown->Sumw2();
  TH1D *hWenuMetm_PileupDown = new TH1D("hWenuMetm_PileupDown","",NBINS,0,METMAX); hWenuMetm_PileupDown->Sumw2();

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet","",  NBINS,0,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,0,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,0,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWenuMet  = new TH1D("hAntiWenuMet","", NBINS,0,METMAX); hAntiWenuMet->Sumw2();
  TH1D *hAntiWenuMetp = new TH1D("hAntiWenuMetp","",NBINS,0,METMAX); hAntiWenuMetp->Sumw2();
  TH1D *hAntiWenuMetm = new TH1D("hAntiWenuMetm","",NBINS,0,METMAX); hAntiWenuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet", "",  NBINS,0,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,0,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,0,METMAX); hAntiEWKMetm->Sumw2();

  cout << "Loading trigger efficiencies..." << endl;

  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_pos = new TFile(zeeHLTEffName_pos);
  CEffUser2D zeeHLTEff_pos;
  zeeHLTEff_pos.loadEff((TH2D*)zeeHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_neg = new TFile(zeeHLTEffName_neg);
  CEffUser2D zeeHLTEff_neg;
  zeeHLTEff_neg.loadEff((TH2D*)zeeHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrhEtaPt"));
   
  TFile *dataHLTEff2BinFile_pos = new TFile(dataHLTEff2BinName_pos);
  CEffUser2D dataHLTEff2Bin_pos;
  dataHLTEff2Bin_pos.loadEff((TH2D*)dataHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEff2BinFile_neg = new TFile(dataHLTEff2BinName_neg);
  CEffUser2D dataHLTEff2Bin_neg;
  dataHLTEff2Bin_neg.loadEff((TH2D*)dataHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrhEtaPt"));
    
  TFile *zeeHLTEff2BinFile_pos = new TFile(zeeHLTEff2BinName_pos);
  CEffUser2D zeeHLTEff2Bin_pos;
  zeeHLTEff2Bin_pos.loadEff((TH2D*)zeeHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zeeHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEff2BinFile_neg = new TFile(zeeHLTEff2BinName_neg);
  CEffUser2D zeeHLTEff2Bin_neg;
  zeeHLTEff2Bin_neg.loadEff((TH2D*)zeeHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zeeHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeHLTEff2BinFile_neg->Get("hErrhEtaPt"));

  //
  // Selection efficiency
  //
  cout << "Loading GSF+selection efficiencies..." << endl;
  
  TFile *dataGsfSelEffFile_pos = new TFile(dataGsfSelEffName_pos);
  CEffUser2D dataGsfSelEff_pos;
  dataGsfSelEff_pos.loadEff((TH2D*)dataGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataGsfSelEffFile_neg = new TFile(dataGsfSelEffName_neg);
  CEffUser2D dataGsfSelEff_neg;
  dataGsfSelEff_neg.loadEff((TH2D*)dataGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_pos = new TFile(zeeGsfSelEffName_pos);
  CEffUser2D zeeGsfSelEff_pos;
  zeeGsfSelEff_pos.loadEff((TH2D*)zeeGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_neg = new TFile(zeeGsfSelEffName_neg);
  CEffUser2D zeeGsfSelEff_neg;
  zeeGsfSelEff_neg.loadEff((TH2D*)zeeGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataGsfSelEff2BinFile_pos = new TFile(dataGsfSelEff2BinName_pos);
  CEffUser2D dataGsfSelEff2Bin_pos;
  dataGsfSelEff2Bin_pos.loadEff((TH2D*)dataGsfSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataGsfSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataGsfSelEff2BinFile_neg = new TFile(dataGsfSelEff2BinName_neg);
  CEffUser2D dataGsfSelEff2Bin_neg;
  dataGsfSelEff2Bin_neg.loadEff((TH2D*)dataGsfSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataGsfSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zeeGsfSelEff2BinFile_pos = new TFile(zeeGsfSelEff2BinName_pos);
  CEffUser2D zeeGsfSelEff2Bin_pos;
  zeeGsfSelEff2Bin_pos.loadEff((TH2D*)zeeGsfSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_pos->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEff2BinFile_neg = new TFile(zeeGsfSelEff2BinName_neg);
  CEffUser2D zeeGsfSelEff2Bin_neg;
  zeeGsfSelEff2Bin_neg.loadEff((TH2D*)zeeGsfSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_neg->Get("hErrhEtaPt"));
 
    
  
    //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0, *lep_raw=0;
//   TLorentzVector *sc=0;
  
  
  TFile *infile=0;
  TTree *intree=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	  assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",  &genVy);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("puppiMet",      &met);       // MET
    intree->SetBranchAddress("puppiMetPhi",   &metPhi);    // phi(MET)
//     intree->SetBranchAddress("mvaMet",      &met);       // MET
//     intree->SetBranchAddress("mvaMetPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",       &mt);        // transverse mass
    intree->SetBranchAddress("puppiU1",       &u1);        // parallel component of recoil
    intree->SetBranchAddress("puppiU2",       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("lep_raw",      &lep_raw);       // lepton 4-vector
//     intree->SetBranchAddress("sc",       &sc);        // electron Supercluster 4-vector
  

  
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) std::cout << "On Entry.... " << ientry << std::endl;
      
      // 2d Vectors to correct MET for lepton scaling/smearing
      TVector2 vLepRaw((lep_raw->Pt())*cos(lep_raw->Phi()),(lep_raw->Pt())*sin(lep_raw->Phi()));
      TVector2 vLepCor((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
    
      double pU1         = 0;  //--
      double pU2         = 0;  //--

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

//       if(lep->Pt()        < PT_CUT)  continue;	
      if(fabs(lep->Eta()) > ETA_CUT) continue;
  
      mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );

      
      
      effdata=1; effmc=1;
      if(q>0) { 
        effdata *= (1.-dataHLTEff_pos.getEff(lep->Eta(), lep->Pt())); 
        effmc   *= (1.-zeeHLTEff_pos.getEff(lep->Eta(), lep->Pt())); 
      } else {
        effdata *= (1.-dataHLTEff_neg.getEff(lep->Eta(), lep->Pt())); 
        effmc   *= (1.-zeeHLTEff_neg.getEff(lep->Eta(), lep->Pt())); 
      }
      effdata = 1.-effdata;
      effmc   = 1.-effmc;
      corr *= effdata/effmc;
      corrSigShape *= effdata/effmc;
      corrBkgShape *= effdata/effmc;
  
      effdata=1; effmc=1;
      effSigShapedata=1;
      effBkgShapedata=1;
      if(q>0) {
        effdata *= dataGsfSelEff_pos.getEff(lep->Eta(), lep->Pt()); 
        effmc   *= zeeGsfSelEff_pos.getEff(lep->Eta(), lep->Pt());
        effSigShapedata *= dataGsfSelEff_pos.getEff((lep->Eta()), lep->Pt())*hGsfSelSigSys->GetBinContent(hGsfSelSigSys->GetXaxis()->FindBin(lep->Eta()), hGsfSelSigSys->GetYaxis()->FindBin(lep->Pt())); 
        effBkgShapedata *= dataGsfSelEff_pos.getEff((lep->Eta()), lep->Pt())*hGsfSelBkgSys->GetBinContent(hGsfSelBkgSys->GetXaxis()->FindBin(lep->Eta()), hGsfSelBkgSys->GetYaxis()->FindBin(lep->Pt()));
      } else {
        effdata *= dataGsfSelEff_neg.getEff(lep->Eta(), lep->Pt()); 
        effmc   *= zeeGsfSelEff_neg.getEff(lep->Eta(), lep->Pt()); 
        effSigShapedata *= dataGsfSelEff_neg.getEff((lep->Eta()), lep->Pt())*hGsfSelSigSys->GetBinContent(hGsfSelSigSys->GetXaxis()->FindBin(lep->Eta()), hGsfSelSigSys->GetYaxis()->FindBin(lep->Pt())); 
        effBkgShapedata *= dataGsfSelEff_neg.getEff((lep->Eta()), lep->Pt())*hGsfSelBkgSys->GetBinContent(hGsfSelBkgSys->GetXaxis()->FindBin(lep->Eta()), hGsfSelBkgSys->GetYaxis()->FindBin(lep->Pt()));
      }
      corr *= effdata/effmc;
      corrSigShape *= effSigShapedata/effmc;
      corrBkgShape *= effBkgShapedata/effmc;
      //corr=1;
      
      if(typev[ifile]==eData) {
        if(lep->Pt()        < PT_CUT)  continue;
        // Correct MET to use raw lepton
        TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
        Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
        hDataMet->Fill(corrMetWithLepton); 
        if(q>0) { hDataMetp->Fill(corrMetWithLepton); } 
        else    { hDataMetm->Fill(corrMetWithLepton); }
      }
      else if(typev[ifile]==eAntiData) {
        if(lep->Pt()        < PT_CUT)  continue;
        // Correct MET to use raw lepton
//         TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
//         Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
        Double_t corrMetWithLepton = met;
        hAntiDataMet->Fill(corrMetWithLepton);
        if(q>0) { hAntiDataMetp->Fill(corrMetWithLepton); } 
        else    { hAntiDataMetm->Fill(corrMetWithLepton); }
      }
      else {
        Double_t weight = 1; Double_t weightUp = 1; Double_t weightDown = 1;
        Double_t weight2 =1;
        if(pileupUp) {
          weight *= scale1fbUp*lumi*corr;
          weight2 *= scale1fbUp*lumi2*corr;
        } else if(pileupDown){
          weight *= scale1fbDown*lumi*corr;
          weight2 *= scale1fbDown*lumi2*corr;
        } else {
          weight2*=scale1fb*lumi2*corr;
          weight *= scale1fb*lumi*corr;
        }

        Double_t eleCorr=0.0;
        for (Int_t i=0; i<9; i++) {
          if( fabs(lep->Eta())>lower[i] && fabs(lep->Eta())<upper[i] ) {
            eleCorr=fpcorr[i];
          }
        }//footprint?

        // apply recoil corrections to W MC
        Double_t lepPt = lep->Pt();
        //lepPt = lepPt + eleCorr*lepPt;
//         Double_t lepPtup = lep->Pt();
//         Double_t lepPtdown = lep->Pt();

        //std::cout << "hey hey " << typev[ifile] << " " << eAntiWenu << " " << eWenu << std::endl;
        if(typev[ifile]==eWenu || typev[ifile]==eBKG) {
        //if(typev[ifile]==eWenu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          if(lepPt       > PT_CUT) {
            double bin = 0;
            for(int i = 1; i <= hh_diff->GetNbinsX();++i){
              if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ 
                bin = i;
                break;
              }
            }
            double w2 = 1.0;//hh_diff->GetBinContent(bin);
            
	    corrMet=met, corrMetPhi=metPhi;
	    hWenuMet_PileupUp->Fill(corrMet,weightUp);
	    hWenuMet_PileupDown->Fill(corrMet,weightDown);
	    if(q>0) {
	      if(doInclusive) recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
          else if(pileupUp)    recoilCorrPuUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
          else if(pileupDown)  recoilCorrPuDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
          else if(doKeys){
            if(fabs(genVy)<0.5)
              recoilCorrKeys05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              recoilCorrKeys051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else
              recoilCorrKeys1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
          } else {
            if(fabs(genVy)<0.5)
              recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else
              recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
          }
          TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
	      if(typev[ifile]==eWenu)
          {// FIX to include vector sum
            hWenuMet->Fill(corrMetWithLepton,weight);
            //corrMet = sqrt( 2.0 * (lep->Pt()) * (corrMet) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),corrMetPhi))) );
            hWenuMetp->Fill(corrMetWithLepton,weight*w2);
            hWenuMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
            hWenuMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
          }
	      else
          {
            hEWKMet->Fill(corrMetWithLepton,weight);
            hEWKMetp->Fill(corrMetWithLepton,weight);
            hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
            hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
          }
	    } else { 
	      if(doInclusive)recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
          else if(pileupUp)    recoilCorrPuUpm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
          else if(pileupDown)  recoilCorrPuDownm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
	      else if(doKeys){
            if(fabs(genVy)<0.5)
              recoilCorrKeysm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              recoilCorrKeysm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else
              recoilCorrKeysm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
           } else {
            if(fabs(genVy)<0.5)
              recoilCorrm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              recoilCorrm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
            else
              recoilCorrm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
           }            
          TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
	      if(typev[ifile]==eWenu)
          {             
            hWenuMet->Fill(corrMetWithLepton,weight);
            hWenuMetm->Fill(corrMetWithLepton,weight*w2); corrMet=met, corrMetPhi=metPhi;
            hWenuMetm_PileupUp->Fill(corrMetWithLepton,weightUp); corrMet=met, corrMetPhi=metPhi;
            hWenuMetm_PileupDown->Fill(corrMetWithLepton,weightDown); corrMet=met, corrMetPhi=metPhi;
          }
	      else
          {
            hEWKMet->Fill(corrMetWithLepton,weight);
            hEWKMetm->Fill(corrMetWithLepton,weight); corrMet=met, corrMetPhi=metPhi; 
            hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); corrMet=met, corrMetPhi=metPhi; 
            hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); corrMet=met, corrMetPhi=metPhi; 
          }
	    }
	    // unused ??
	    corrMet=met, corrMetPhi=metPhi;
	    hWenuMet_RecoilUp->Fill(corrMet,weight);
	    if(q>0) {
	      hWenuMetp_RecoilUp->Fill(corrMet,weight); 
	      corrMet=met, corrMetPhi=metPhi;
	    } else { 
	      hWenuMetm_RecoilUp->Fill(corrMet,weight); 
	      corrMet=met, corrMetPhi=metPhi;
	    }
	    corrMet=met, corrMetPhi=metPhi;
	    hWenuMet_RecoilDown->Fill(corrMet,weight);
	    corrMet=met, corrMetPhi=metPhi;
	    if(q>0){ 
          hWenuMetp_RecoilDown->Fill(corrMet,weight); 
          corrMet=met, corrMetPhi=metPhi;
        } else { 
          hWenuMetm_RecoilDown->Fill(corrMet,weight); 
          corrMet=met, corrMetPhi=metPhi;
        }
	  }
	}
	else if(typev[ifile]==eAntiWenu)
	  {
	    Double_t corrMet=met, corrMetPhi=metPhi;
	    if(lepPt        > PT_CUT) {
	      hAntiWenuMet->Fill(corrMet,weight2);
	      if(q>0) {
            if(doInclusive)recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupUp)    recoilCorrPuUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupDown)  recoilCorrPuDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
            else if(doKeys){
              if(fabs(genVy)<0.5)
                recoilCorrKeys05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorrKeys051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorrKeys1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            } else {
              if(fabs(genVy)<0.5)
                recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            }  
//             TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
//             Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            Double_t corrMetWithLepton = met;
            hAntiWenuMetp->Fill(corrMetWithLepton,weight2); // *w2 
            corrMet=met, corrMetPhi=metPhi;
	      } else {
            if(doInclusive)recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupUp)    recoilCorrPuUpm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupDown)  recoilCorrPuDownm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0);
            else if(doKeys){
              if(fabs(genVy)<0.5)
                recoilCorrKeysm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorrKeysm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorrKeysm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            } else {
              if(fabs(genVy)<0.5)
                recoilCorrm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorrm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorrm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep_raw->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            }  
//             TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
//             Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            Double_t corrMetWithLepton = met;
            hAntiWenuMetm->Fill(corrMetWithLepton,weight2); //*w2);
            corrMet=met, corrMetPhi=metPhi;
	      }
	    }
	  }
        if(typev[ifile]==eEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          hEWKMet->Fill(corrMetWithLepton,weight);
          hEWKMet_PileupUp->Fill(corrMetWithLepton,weightUp);
          hEWKMet_PileupDown->Fill(corrMetWithLepton,weightDown);
          if(q>0) {
            hEWKMetp->Fill(corrMetWithLepton,weight); 
            hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp); 
            hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown); 
          }
          else { 
            hEWKMetm->Fill(corrMetWithLepton,weight); 
            hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); 
            hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); 
          }
        }
        if(typev[ifile]==eAntiEWK) { 
          if(lep->Pt()        < PT_CUT)  continue;
          // don't have raw lepton info in ntuple at the moment
//           TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
//           Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          Double_t corrMetWithLepton = met;
          hAntiEWKMet->Fill(corrMetWithLepton,weight2);
          if(q>0) { hAntiEWKMetp->Fill(corrMetWithLepton,weight2); }
          else    { hAntiEWKMetm->Fill(corrMetWithLepton,weight2); }
        }
      }
    }
  }  
  delete infile;
  infile=0, intree=0;   
  
  
  //   Calculate the shapes for W+
  TH1D *corrP = (TH1D*) hWenuMetp->Clone("up");
  *corrP = (*hh_diffp)*(*hWenuMetp);
  hWenuMetp_RecoilUp->Add(corrP,hWenuMetp,-1);
  hWenuMetp_RecoilDown->Add(corrP,hWenuMetp,1);
//   hWenuMetp_ScaleUp->Add(corrP,hWenuMetp,-1);
//   hWenuMetp_ScaleDown->Add(corrP,hWenuMetp,1);
  // Calculate the shapes for W-
  TH1D *corrM = (TH1D*) hWenuMetm->Clone("up");
  *corrM = (*hh_diffm)*(*hWenuMetm);
  std::cout << hh_diffm->GetBinContent(2)<< std::endl;
  hWenuMetm_RecoilUp->Add(corrM,hWenuMetm,-1);
  std::cout << corrM->GetBinContent(2)<< std::endl;
  std::cout << hWenuMetm->GetBinContent(2)<< std::endl;
  std::cout << hWenuMetm_RecoilUp->GetBinContent(2)<< std::endl;
  hWenuMetm_RecoilDown->Add(corrM,hWenuMetm,1);
//   hWenuMetm_ScaleUp->Add(corrM,hWenuMetm,-1);
//   hWenuMetm_ScaleDown->Add(corrM,hWenuMetm,1);
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  //RooRealVar nSig("nSig","nSig",(hWenuMet->Integral()),0,hWenuMet->Integral());
  //RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  //RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWenuMet->Integral());
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWenuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hAntiDataMet->Integral()),0,hAntiDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWenuMet->Integral());
  dewk.setConstant(kTRUE);
  //   nAntiSig.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));

  //hWenuMetp->Scale(6.7446e+06/hWenuMetp->Integral());
  //RooRealVar nSigp("nSigp","nSigp",hDataMetp->Integral()*0.8,0,hDataMetp->Integral());
  RooRealVar nSigp("nSigp","nSigp",1.0*hWenuMetp->Integral(),0,hDataMetp->Integral());
  //nSigp.setConstant(kTRUE);
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",0.2*hDataMetp->Integral(),0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5);
  cewkp.setVal(hEWKMetp->Integral()/hWenuMetp->Integral());
  cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWenuMetp->Integral(),0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.95*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWenuMetp->Integral());
  //std::cout << "TTTT " << hAntiEWKMetp->Integral() << " "<< hAntiWenuMetp->Integral() << std::endl;
  dewkp.setConstant(kTRUE);
//   nAntiSigp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));

  //hWenuMetm->Scale(5.2551e+06/hWenuMetm->Integral());
  RooRealVar nSigm("nSigm","nSigm",1.0*hWenuMetm->Integral(),0,hDataMetm->Integral());
  //RooRealVar nSigm("nSigm","nSigm",hDataMetm->Integral()*0.8,0,hDataMetm->Integral());
  //nSigm.setConstant(kTRUE);
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm->Integral()*0.2,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWenuMetm->Integral());
  cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
   RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWenuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.95*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWenuMetm->Integral());
  dewkm.setConstant(kTRUE);
//   nAntiSigm.setConstant(kTRUE);
  RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));

  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wenuMet ("wenuMET", "wenuMET", RooArgSet(pfmet),hWenuMet);  RooHistPdf pdfWe ("we", "we", pfmet,wenuMet, 1);
  RooDataHist wenuMetp("wenuMETp","wenuMETp",RooArgSet(pfmet),hWenuMetp); RooHistPdf pdfWep("wep","wep",pfmet,wenuMetp,1);
  RooDataHist wenuMetm("wenuMETm","wenuMETm",RooArgSet(pfmet),hWenuMetm); RooHistPdf pdfWem("wem","wem",pfmet,wenuMetm,1); 
  
  RooDataHist wenuMet_RecoilUp("wenuMET_RecoilUp", "wenuMET_RecoilUp", RooArgSet(pfmet),hWenuMet_RecoilUp);  RooHistPdf pdfWe_RecoilUp("we_RecoilUp", "we_RecoilUp", pfmet,wenuMet_RecoilUp, 1);
  RooDataHist wenuMetp_RecoilUp("wenuMETp_RecoilUp","wenuMETp_RecoilUp",RooArgSet(pfmet),hWenuMetp_RecoilUp); RooHistPdf pdfWep_RecoilUp("wep_RecoilUp","wep_RecoilUp",pfmet,wenuMetp_RecoilUp,1);
  
  RooDataHist wenuMetm_RecoilUp("wenuMETm_RecoilUp","wenuMETm_RecoilUp",RooArgSet(pfmet),hWenuMetm_RecoilUp); RooHistPdf pdfWem_RecoilUp("wem_RecoilUp","wem_RecoilUp",pfmet,wenuMetm_RecoilUp,1); 
  
  RooDataHist wenuMet_RecoilDown("wenuMET_RecoilDown", "wenuMET_RecoilDown", RooArgSet(pfmet),hWenuMet_RecoilDown);  RooHistPdf pdfWe_RecoilDown("we_RecoilDown", "we_RecoilDown", pfmet,wenuMet_RecoilDown, 1);
  RooDataHist wenuMetp_RecoilDown("wenuMETp_RecoilDown","wenuMETp_RecoilDown",RooArgSet(pfmet),hWenuMetp_RecoilDown); RooHistPdf pdfWep_RecoilDown("wep_RecoilDown","wep_RecoilDown",pfmet,wenuMetp_RecoilDown,1);
  RooDataHist wenuMetm_RecoilDown("wenuMETm_RecoilDown","wenuMETm_RecoilDown",RooArgSet(pfmet),hWenuMetm_RecoilDown); RooHistPdf pdfWem_RecoilDown("wem_RecoilDown","wem_RecoilDown",pfmet,wenuMetm_RecoilDown,1);

  
  RooDataHist wenuMet_ScaleUp("wenuMET_ScaleUp", "wenuMET_ScaleUp", RooArgSet(pfmet),hWenuMet_ScaleUp);  RooHistPdf pdfWe_ScaleUp("we_ScaleUp", "we_ScaleUp", pfmet,wenuMet_ScaleUp, 1);
  RooDataHist wenuMetp_ScaleUp("wenuMETp_ScaleUp","wenuMETp_ScaleUp",RooArgSet(pfmet),hWenuMetp_ScaleUp); RooHistPdf pdfWep_ScaleUp("wep_ScaleUp","wep_ScaleUp",pfmet,wenuMetp_ScaleUp,1);
  RooDataHist wenuMetm_ScaleUp("wenuMETm_ScaleUp","wenuMETm_ScaleUp",RooArgSet(pfmet),hWenuMetm_ScaleUp); RooHistPdf pdfWem_ScaleUp("wem_ScaleUp","wem_ScaleUp",pfmet,wenuMetm_ScaleUp,1); 
  RooDataHist wenuMet_ScaleDown("wenuMET_ScaleDown", "wenuMET_ScaleDown", RooArgSet(pfmet),hWenuMet_ScaleDown);  RooHistPdf pdfWe_ScaleDown("we_ScaleDown", "we_ScaleDown", pfmet,wenuMet_ScaleDown, 1);
  RooDataHist wenuMetp_ScaleDown("wenuMETp_ScaleDown","wenuMETp_ScaleDown",RooArgSet(pfmet),hWenuMetp_ScaleDown); RooHistPdf pdfWep_ScaleDown("wep_ScaleDown","wep_ScaleDown",pfmet,wenuMetp_ScaleDown,1);
  RooDataHist wenuMetm_ScaleDown("wenuMETm_ScaleDown","wenuMETm_ScaleDown",RooArgSet(pfmet),hWenuMetm_ScaleDown); RooHistPdf pdfWem_ScaleDown("wem_ScaleDown","wem_ScaleDown",pfmet,wenuMetm_ScaleDown,1); 
  
  RooDataHist wenuMet_PileupUp("wenuMET_PileupUp", "wenuMET_PileupUp", RooArgSet(pfmet),hWenuMet_PileupUp);  RooHistPdf pdfWe_PileupUp("we_PileupUp", "we_PileupUp", pfmet,wenuMet_PileupUp, 1);
  RooDataHist wenuMetp_PileupUp("wenuMETp_PileupUp","wenuMETp_PileupUp",RooArgSet(pfmet),hWenuMetp_PileupUp); RooHistPdf pdfWep_PileupUp("wep_PileupUp","wep_PileupUp",pfmet,wenuMetp_PileupUp,1);
  RooDataHist wenuMetm_PileupUp("wenuMETm_PileupUp","wenuMETm_PileupUp",RooArgSet(pfmet),hWenuMetm_PileupUp); RooHistPdf pdfWem_PileupUp("wem_PileupUp","wem_PileupUp",pfmet,wenuMetm_PileupUp,1); 
  RooDataHist wenuMet_PileupDown("wenuMET_PileupDown", "wenuMET_PileupDown", RooArgSet(pfmet),hWenuMet_PileupDown);  RooHistPdf pdfWe_PileupDown("we_PileupDown", "we_PileupDown", pfmet,wenuMet_PileupDown, 1);
  RooDataHist wenuMetp_PileupDown("wenuMETp_PileupDown","wenuMETp_PileupDown",RooArgSet(pfmet),hWenuMetp_PileupDown); RooHistPdf pdfWep_PileupDown("wep_PileupDown","wep_PileupDown",pfmet,wenuMetp_PileupDown,1);
  RooDataHist wenuMetm_PileupDown("wenuMETm_PileupDown","wenuMETm_PileupDown",RooArgSet(pfmet),hWenuMetm_PileupDown); RooHistPdf pdfWem_PileupDown("wem_PileupDown","wem_PileupDown",pfmet,wenuMetm_PileupDown,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
  RooDataHist ewkMet_PileupUp("ewkMET_PileupUp", "ewkMET_PileupUp", RooArgSet(pfmet),hEWKMet_PileupUp);  RooHistPdf pdfEWK_PileupUp("ewk_PileupUp", "ewk_PileupUp", pfmet,ewkMet_PileupUp, 1);
  RooDataHist ewkMetp_PileupUp("ewkMETp_PileupUp","ewkMETp_PileupUp",RooArgSet(pfmet),hEWKMetp_PileupUp); RooHistPdf pdfEWKp_PileupUp("ewkp_PileupUp","ewkp_PileupUp",pfmet,ewkMetp_PileupUp,1);
  RooDataHist ewkMetm_PileupUp("ewkMETm_PileupUp","ewkMETm_PileupUp",RooArgSet(pfmet),hEWKMetm_PileupUp); RooHistPdf pdfEWKm_PileupUp("ewkm_PileupUp","ewkm_PileupUp",pfmet,ewkMetm_PileupUp,1); 
  RooDataHist ewkMet_PileupDown("ewkMET_PileupDown", "ewkMET_PileupDown", RooArgSet(pfmet),hEWKMet_PileupDown);  RooHistPdf pdfEWK_PileupDown("ewk_PileupDown", "ewk_PileupDown", pfmet,ewkMet_PileupDown, 1);
  RooDataHist ewkMetp_PileupDown("ewkMETp_PileupDown","ewkMETp_PileupDown",RooArgSet(pfmet),hEWKMetp_PileupDown); RooHistPdf pdfEWKp_PileupDown("ewkp_PileupDown","ewkp_PileupDown",pfmet,ewkMetp_PileupDown,1);
  RooDataHist ewkMetm_PileupDown("ewkMETm_PileupDown","ewkMETm_PileupDown",RooArgSet(pfmet),hEWKMetm_PileupDown); RooHistPdf pdfEWKm_PileupDown("ewkm_PileupDown","ewkm_PileupDown",pfmet,ewkMetm_PileupDown,1); 
  
  // Anti-Signal PDFs
  RooDataHist awenuMet ("awenuMET", "awenuMET", RooArgSet(pfmet),hAntiWenuMet);
  RooHistPdf apdfWe ("awe", "awe", pfmet,awenuMet, 1);
  RooDataHist awenuMetp("awenuMETp","awenuMETp",RooArgSet(pfmet),hAntiWenuMetp);
  RooHistPdf apdfWep("awep","awep",pfmet,awenuMetp,1);
  RooDataHist awenuMetm("awenuMETm","awenuMETm",RooArgSet(pfmet),hAntiWenuMetm);
  RooHistPdf apdfWem("awem","awem",pfmet,awenuMetm,1); 
  
  // Anti-EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 

//   RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
//   RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
//   RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  // QCD Pdfs
  //CExponential qcd(pfmet,kTRUE);
  //CExponential qcdp(pfmet,kTRUE);
  //CExponential qcdm(pfmet,kTRUE);
   CPepeModel2 qcd("qcd",pfmet);
   CPepeModel2 qcdp("qcdp",pfmet);
   CPepeModel2 qcdm("qcdm",pfmet);

   CPepeModel2 aqcd("aqcd",pfmet, qcd.a1);
   CPepeModel2 aqcdp("aqcdp",pfmet, qcdp.a1);
   CPepeModel2 aqcdm("aqcdm",pfmet, qcdm.a1);
//    CPepeModel2 aqcd("aqcd",pfmet);
//    CPepeModel2 aqcdp("aqcdp",pfmet);
//    CPepeModel2 aqcdm("aqcdm",pfmet);

   //CPepeModel1 aqcdp("aqcdp",pfmet,qcdp.a1,qcdp.sigma);
   //CPepeModel1 aqcdm("aqcdm",pfmet,qcdm.a1,qcdm.sigma);

   //qcdp.a1 = qcdm.a1;
   
//    CPepeModel2 aqcd("aqcd",pfmet);
//    CPepeModel2 aqcdp("aqcdp",pfmet);
//    CPepeModel2 aqcdm("aqcdm",pfmet);  
   
   //RooRealVar mm("mm","mm",0.0185,-1,1) ;
   //RooRealVar sm("sm","sm",0.0006,0.0,10) ;
   //RooRealVar mp("mp","mp",0.015826,-1,1) ;
   //RooRealVar sp("sp","sp",0.0003,0.0,10) ;
   //RooGaussian gaussm("gaussm","gauss(xm,mm,sm)",qcdm.a1,mm,sm) ;
   //RooGaussian gaussp("gaussp","gauss(xp,mp,sp)",qcdp.a1,mp,sp) ;
   //RooGaussian constm("constm","constm",*qcdm.a1,RooConst(0.0185),RooConst(0.0003)) ;
   //RooGaussian constp("constp","constp",*qcdp.a1,RooConst(0.015826),RooConst(0.0001)) ;
   RooGaussian constm("constm","constm",nAntiSigm,RooConst(hAntiWenuMetm->Integral()),RooConst(0.05*hAntiWenuMetm->Integral()));
   RooGaussian constp("constp","constp",nAntiSigp,RooConst(hAntiWenuMetp->Integral()),RooConst(0.05*hAntiWenuMetp->Integral()));
   //RooGaussian constp("constp","constp",*qcdp.a1,RooConst(0.015826),RooConst(0.0003)) ;
   

  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("Selectp");
  rooCat.defineType("Selectm");

  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWe,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWep,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWem,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm)); 

  // Anti-selection PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWe,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWep,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWem,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
//     RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWe,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWep,pdfEWKp,qcdpc),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWem,pdfEWKm,qcdmc),RooArgList(nSigm,nEWKm,nQCDm)); 
  
//     RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWe,pdfEWK,pdfQCD),RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWep,pdfEWKp,pdfQCD),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWem,pdfEWKm,pdfQCD),RooArgList(nSigm,nEWKm,nQCDm)); 
    
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(pdfMetp, "Selectp");
  pdfTotal.addPdf(pdfMetm,"Selectm");

  RooSimultaneous pdfTotalp("pdfTotalp","pdfTotalp",rooCat);
  pdfTotalp.addPdf(pdfMetp, "Select");
  pdfTotalp.addPdf(apdfMetp,"Anti");
  
  RooSimultaneous pdfTotalm("pdfTotalm","pdfTotalm",rooCat);
  pdfTotalm.addPdf(pdfMetm, "Select");
  pdfTotalm.addPdf(apdfMetm,"Anti");

  //
  // Perform fits
  //
  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);
 
  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
  RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
  RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);

  RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
			 Import("Select", dataMetp),
			 Import("Anti",   antiMetp));

  RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
			 Import("Select", dataMetm),
			 Import("Anti", antiMetm));
  
  cout << "Starting values for Wenu yields: " << endl;
  cout << "Selected: " << hDataMet->Integral() << endl;
  cout << "   sig: " << hWenuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;
  cout << "   qcd: " << hDataMet->Integral()-hWenuMet->Integral()-hEWKMet->Integral() << endl;

  cout << "Starting values for Wenu_p yields: " << endl;
  cout << "Selected: " << hDataMetp->Integral() << endl;
  cout << "   sig: " << hWenuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hDataMetp->Integral()-hWenuMetp->Integral()-hEWKMetp->Integral() << endl;

  cout << "Starting values for Wenu_m yields: " << endl;
  cout << "Selected: " << hDataMetm->Integral() << endl;
  cout << "   sig: " << hWenuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  cout << "   qcd: " << hDataMetm->Integral()-hWenuMetm->Integral()-hEWKMetm->Integral() << endl;

  cout << "Starting values for AntiWenu_p yields: " << endl;
  cout << "Selected: " << hAntiDataMetp->Integral() << endl;
  cout << "   sig: " << hAntiWenuMetp->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetp->Integral() << endl;

  cout << "Starting values for AntiWenu_m yields: " << endl;
  cout << "Selected: " << hAntiDataMetm->Integral() << endl;
  cout << "   sig: " << hAntiWenuMetm->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());

  RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  
  RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",0.3*(hAntiDataMet->Integral()),0,hAntiDataMet->Integral());
  RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",0.3*(hAntiDataMet->Integral()),0,hAntiDataMet->Integral());
  
  
//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",40000,0,100000);
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",40000,0,100000);

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);
  combine_workspace.import(pepe2Pdf_qcdp_norm);
  combine_workspace.import(pepe2Pdf_qcdm_norm);

  combine_workspace.import(pdfWe);
  combine_workspace.import(pdfWep);
  combine_workspace.import(pdfWem);
  combine_workspace.import(pdfWe_RecoilUp);
  combine_workspace.import(pdfWep_RecoilUp);
  combine_workspace.import(pdfWem_RecoilUp);
  combine_workspace.import(pdfWe_RecoilDown);
  combine_workspace.import(pdfWep_RecoilDown);
  combine_workspace.import(pdfWem_RecoilDown);
  combine_workspace.import(pdfWe_ScaleUp);
  combine_workspace.import(pdfWep_ScaleUp);
  combine_workspace.import(pdfWem_ScaleUp);
  combine_workspace.import(pdfWe_ScaleDown);
  combine_workspace.import(pdfWep_ScaleDown);
  combine_workspace.import(pdfWem_ScaleDown);
  combine_workspace.import(pdfWe_PileupUp);
  combine_workspace.import(pdfWep_PileupUp);
  combine_workspace.import(pdfWem_PileupUp);
  combine_workspace.import(pdfWe_PileupDown);
  combine_workspace.import(pdfWep_PileupDown);
  combine_workspace.import(pdfWem_PileupDown);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKm);
  combine_workspace.import(pdfEWK_PileupUp);
  combine_workspace.import(pdfEWKp_PileupUp);
  combine_workspace.import(pdfEWKm_PileupUp);
  combine_workspace.import(pdfEWK_PileupDown);
  combine_workspace.import(pdfEWKp_PileupDown);
  combine_workspace.import(pdfEWKm_PileupDown);
  combine_workspace.import(*(qcd.model));
  //combine_workspace.import(qcdpn);
  combine_workspace.import(*(qcdp.model));
  combine_workspace.import(*(qcdm.model));
  
  combine_workspace.import(antiMet);
  combine_workspace.import(antiMetp);
  combine_workspace.import(antiMetm);
  combine_workspace.import(pepe2Pdf_aqcdp_norm);
  combine_workspace.import(pepe2Pdf_aqcdm_norm);
  
  combine_workspace.import(apdfWe);
  combine_workspace.import(apdfWep);
  combine_workspace.import(apdfWem);
  combine_workspace.import(apdfEWK);
  combine_workspace.import(apdfEWKp);
  combine_workspace.import(apdfEWKm);
  combine_workspace.import(*(aqcd.model));
  //combine_workspace.import(qcdpn);
  combine_workspace.import(*(aqcdp.model));
  combine_workspace.import(*(aqcdm.model));

//   combine_workspace.import(pdfQCD);
//   //combine_workspace.import(qcdpn);
//   combine_workspace.import(pdfQCDp);
//   combine_workspace.import(pdfQCDm);

  combine_workspace.writeToFile("Wenu_pdfTemplates.root");

  RooFitResult *fitRes = pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE)); 
//   RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),ExternalConstraints(constp),Minos(kTRUE),Save(kTRUE));
  RooFitResult *fitResp = pdfTotalp.fitTo(dataTotalp,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntip = apdfMetp.fitTo(antiMetp,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntip = apdfMetp.fitTo(antiMetp,Extended(),ExternalConstraints(constm),Minos(kTRUE),Save(kTRUE));
  RooDataHist dataTotal("dataTotal,","dataTotal,", RooArgList(pfmet), Index(rooCat),Import("Selectm", dataMetm),Import("Selectp",   antiMetm));
  
  //RooFitResult *fitResm = pdfTotal.fitTo(dataTotal,Extended(),Minos(kTRUE),Save(kTRUE));
  RooFitResult *fitResm = pdfTotalm.fitTo(dataTotalm,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntim = apdfMetm.fitTo(antiMetm,Extended(),ExternalConstraints(constm),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntim = apdfMetm.fitTo(antiMetm,Extended(),Minos(kTRUE),Save(kTRUE));
    
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWenuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWenuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle);
  hMetmDiff->SetMarkerSize(0.9);
  
  // the diff hists for the anti-selection
  
  TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  
  hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetp->GetNbinsX(); ++ibin){hPdfAntiMetp->SetBinError(ibin, hAntiWenuMetp->GetBinError(ibin));}
  hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetm->GetNbinsX(); ++ibin){hPdfAntiMetm->SetBinError(ibin, hAntiWenuMetm->GetBinError(ibin));}
  hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  hAntiMetmDiff->SetMarkerStyle(kFullCircle);
  hAntiMetmDiff->SetMarkerSize(0.9);
    

   
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.02);
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
  gStyle->SetTitleOffset(1.100,"Y");
  TGaxis::SetMaxDigits(3);
  
  char ylabel[100];  // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  (8 TeV)",lumi*1000.);
  else         sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000.);
  
  // plot colors
  Int_t linecolorW   = kOrange-3;
  Int_t fillcolorW   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorQCD = kViolet+2;
  Int_t fillcolorQCD = kViolet-5;
  Int_t ratioColor   = kGray+2;
  
  //
  // Dummy histograms for TLegend
  // (I can't figure out how to properly pass RooFit objects...)
  //
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  
  TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  hDummyW->SetLineColor(linecolorW);
  hDummyW->SetFillColor(fillcolorW);
  hDummyW->SetFillStyle(1001);
  
  TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  hDummyEWK->SetLineColor(linecolorEWK);
  hDummyEWK->SetFillColor(fillcolorEWK);
  hDummyEWK->SetFillStyle(1001);
  
  TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  hDummyQCD->SetLineColor(linecolorQCD);
  hDummyQCD->SetFillColor(fillcolorQCD);
  hDummyQCD->SetFillStyle(1001);
   
  //
  // W MET plot
  //
  RooPlot *weframe = pfmet.frame(Bins(NBINS));
  weframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(weframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(weframe,LineColor(linecolorW));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
  pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfWe)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("wenu_fitmet",weframe,"","mT [GeV]",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrowe#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox("CMS",0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox(lumitext,0.63,0.92,0.95,0.99,0);
  plotMet.Draw(c,kTRUE,format,1);

  CPlot plotMetDiff("wenu_fitmet","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-0.20,0.20);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  plotMetDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMet.SetName("wenu_fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
    
  //
  // W+ MET plot
  //
  RooPlot *wepframe = pfmet.frame(Bins(NBINS));    
  wepframe->GetYaxis()->SetNdivisions(505);
  wepframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wepframe,LineColor(linecolorW));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWep)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("wenu_fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrowe^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("wenu_fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.20,0.20);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("wenu_fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  plotMetp.Draw(c,kTRUE,"pdf",1);
  
  // Anti-selection W+ background fits
  RooPlot *awepframe = pfmet.frame(Bins(NBINS));    
  awepframe->GetYaxis()->SetNdivisions(505);
  awepframe->GetXaxis()->SetLabelOffset(2.0);
  antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetp.plotOn(awepframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetp.plotOn(awepframe,LineColor(linecolorW));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),LineColor(linecolorEWK));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(*(aqcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(*(aqcdp.model))),LineColor(linecolorQCD));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfWep)),LineColor(linecolorW),LineStyle(2));
  antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
  CPlot plotAntiMetp("wenu_fitantimetp",awepframe,"","",ylabel);
  plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrowe^{+}#nu","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotAntiMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetp.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetpDiff("wenu_fitantimetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
  plotAntiMetpDiff.SetYRange(-0.20,0.20);
  plotAntiMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotAntiMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotAntiMetpDiff.Draw(c,kTRUE,format,2);
  plotAntiMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotAntiMetp.SetName("wenu_fitantimetplog");
  plotAntiMetp.SetLogy();
  plotAntiMetp.SetYRange(1e-3*(hAntiDataMetp->GetMaximum()),10*(hAntiDataMetp->GetMaximum()));
  plotAntiMetp.Draw(c,kTRUE,format,1);
  plotAntiMetp.Draw(c,kTRUE,"pdf",1);
  
  //
  // W- MET plot
  //
  RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
  wemframe->GetYaxis()->SetNdivisions(505);
  wemframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wemframe,LineColor(linecolorW));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfWem)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("wenu_fitmetm",wemframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrowe^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("wenu_fitmetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetmDiff.SetYRange(-0.2,0.2);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("wenu_fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
  plotMetm.Draw(c,kTRUE,"pdf",1);
  
  // Anti-selection W- background fits
  // fix these
  RooPlot *awemframe = pfmet.frame(Bins(NBINS)); 
  awemframe->GetYaxis()->SetNdivisions(505);
  awemframe->GetXaxis()->SetLabelOffset(2.0);
  antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetm.plotOn(awemframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetm.plotOn(awemframe,LineColor(linecolorW));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),LineColor(linecolorEWK));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(*(aqcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(*(aqcdm.model))),LineColor(linecolorQCD));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfWem)),LineColor(linecolorW),LineStyle(2));
  antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetm->GetBinWidth(1));
  CPlot plotAntiMetm("wenu_fitantimetm",awemframe,"","",ylabel);
  plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrowe^{-}#bar{#nu}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotAntiMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetm.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetmDiff("wenu_fitantimetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotAntiMetmDiff.SetYRange(-0.2,0.2);
  plotAntiMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotAntiMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotAntiMetmDiff.AddHist1D(hAntiMetmDiff,"EX0",ratioColor);
  plotAntiMetmDiff.Draw(c,kTRUE,format,2);
  plotAntiMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotAntiMetm.SetName("wenu_fitantimetmlog");
  plotAntiMetm.SetLogy();
  plotAntiMetm.SetYRange(1e-3*(hAntiDataMetm->GetMaximum()),10*(hAntiDataMetm->GetMaximum()));
  plotAntiMetm.Draw(c,kTRUE,format,1);
  plotAntiMetm.Draw(c,kTRUE,"pdf",1);
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  
  ios_base::fmtflags flags;
  
  Double_t chi2prob, chi2ndf;
  Double_t ksprob, ksprobpe;
  
  chi2prob = hDataMet->Chi2Test(hPdfMet,"PUW");
  chi2ndf  = hDataMet->Chi2Test(hPdfMet,"CHI2/NDFUW");
  ksprob   = hDataMet->KolmogorovTest(hPdfMet);
  ksprobpe = hDataMet->KolmogorovTest(hPdfMet,"DX");
  sprintf(txtfname,"%s/fitresWe.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMet->Integral() << endl;
  txtfile << "  Signal: " << nSig.getVal() << " +/- " << nSig.getPropagatedError(*fitRes) << endl;
  txtfile << "     QCD: " << nQCD.getVal() << " +/- " << nQCD.getPropagatedError(*fitRes) << endl;
  txtfile << "   Other: " << nEWK.getVal() << " +/- " << nEWK.getPropagatedError(*fitRes) << endl;
  txtfile << endl; 
  txtfile.flags(flags);
  
  fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitRes);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();
  
  chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
  ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
  sprintf(txtfname,"%s/fitresWep.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetp->Integral() << endl;
  txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << "AntiSelected: " << hAntiDataMetp->Integral() << endl;
  txtfile << "  AntiSignal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResp) << endl;
  txtfile << "     AntiQCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResp) << endl;
  txtfile << "   AntiOther: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << endl;  
  txtfile.flags(flags);
  
  fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResp);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  ksprob   = hDataMetm->KolmogorovTest(hPdfMetm);
  ksprobpe = hDataMetm->KolmogorovTest(hPdfMetm,"DX");  
  sprintf(txtfname,"%s/fitresWem.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetm->Integral() << endl;
  txtfile << "  Signal: " << nSigm.getVal() << " +/- " << nSigm.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nQCDm.getVal() << " +/- " << nQCDm.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nEWKm.getVal() << " +/- " << nEWKm.getPropagatedError(*fitResm) << endl;
  
  txtfile << "AntiSelected: " << hAntiDataMetm->Integral() << endl;
  txtfile << "  Signal: " << nAntiSigm.getVal() << " +/- " << nAntiSigm.getPropagatedError(*fitResm) << endl;
  txtfile << "     QCD: " << nAntiQCDm.getVal() << " +/- " << nAntiQCDm.getPropagatedError(*fitResm) << endl;
  txtfile << "   Other: " << nAntiEWKm.getVal() << " +/- " << nAntiEWKm.getPropagatedError(*fitResm) << endl;
  txtfile << endl;
  txtfile.flags(flags);
  
  fitResm->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResm);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitWe");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
/*TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
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
}*/
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    std::cout << "bin # " << ibin << std::endl;
    std::cout << "fit bin content " << hFit->GetBinContent(ibin) << std::endl;
    std::cout << "fit bin err " << hFit->GetBinError(ibin) << std::endl;
    std::cout << "data bin content " << hData->GetBinContent(ibin) << std::endl;
    std::cout << "data bin err " << hData->GetBinError(ibin) << std::endl;
    std::cout << "diff =  " << diff << std::endl;
//     Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    std::cout << "err = " << err << std::endl;
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
//--------------------------------------------------------------------------------------------------
void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList parlist = res->floatParsFinal();
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist.getSize(); i++) {
    for(Int_t j=0; j<parlist.getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe)
{
  ios_base::fmtflags flags = os.flags();
  
  os << "  Chi2 Test" << endl;
  os << " -----------" << endl;
  os << "       prob = " << chi2prob << endl;
  os << "   chi2/ndf = " << chi2ndf << endl;
  os << endl;
  os << "  KS Test" << endl;
  os << " ---------" << endl;
  os << "   prob = " << ksprob << endl;
  os << "   prob = " << ksprobpe << " with 1000 pseudo-experiments" << endl;
  os << endl;
 
  os.flags(flags);
}
			   
//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/WenuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wenu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmet.png\"><img src=\"fitmet.png\" alt=\"fitmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetp.png\"><img src=\"fitmetp.png\" alt=\"fitmetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetm.png\"><img src=\"fitmetm.png\" alt=\"fitmetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetlog.png\"><img src=\"fitmetlog.png\" alt=\"fitmetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetplog.png\"><img src=\"fitmetplog.png\" alt=\"fitmetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetmlog.png\"><img src=\"fitmetmlog.png\" alt=\"fitmetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}
