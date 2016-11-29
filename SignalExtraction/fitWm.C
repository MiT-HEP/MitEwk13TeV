//================================================================================================
//
// Perform fit to extract W->munu signal
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

#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_asym2.hh"
// #include "../Utils/RecoilCorrector_addJets.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// #include "ZBackgrounds.hh"

//helper class to handle rochester corrections
#include <rochcor2015r.h>
#include <muresolution_run2r.h>


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
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
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

void fitWm(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)'
           const Double_t lumi2, // lumi for the anti-isolation trigger
       const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // some flags to handle Recoil corrections
  bool doKeys = true;
  bool doInclusive = false;
  // some flags to handle the pileup Up/Down systematics
  bool pileupUp = false;
  bool pileupDown = false;
  
  
  // MET histogram binning and range
  const Int_t    NBINS   = 75;
  const Double_t METMAX  = 150;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t mu_MASS = 0.1057;
  const int NTOYS = 100;
  
    // efficiency files
 
  const TString baseDir = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/"; 
  const TString dataHLTEffName_pos = baseDir + "MuHLTEff/MGpositive/eff.root";
  const TString dataHLTEffName_neg = baseDir + "MuHLTEff/MGnegative/eff.root";
  const TString zmmHLTEffName_pos  = baseDir + "MuHLTEff/CTpositive/eff.root";
  const TString zmmHLTEffName_neg  = baseDir + "MuHLTEff/CTnegative/eff.root";

  const TString dataSelEffName_pos = baseDir + "MuSITEff/MGpositive/eff.root";
  const TString dataSelEffName_neg = baseDir + "MuSITEff/MGnegative/eff.root";
  const TString zmmSelEffName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmSelEffName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataTrkEffName_pos = baseDir + "MuSITEff/MGpositive/eff.root";
  const TString dataTrkEffName_neg = baseDir + "MuSITEff/MGnegative/eff.root";
  const TString zmmTrkEffName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmTrkEffName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataStaEffName_pos = baseDir + "MuStaEff/MGpositive/eff.root";
  const TString dataStaEffName_neg = baseDir + "MuStaEff/MGnegative/eff.root";
  const TString zmmStaEffName_pos  = baseDir + "MuStaEff/CTpositive/eff.root";
  const TString zmmStaEffName_neg  = baseDir + "MuStaEff/CTnegative/eff.root";

  // efficiency files 2Bins

  const TString dataHLTEff2BinName_pos = baseDir + "MuHLTEff/MGpositive/eff.root";
  const TString dataHLTEff2BinName_neg = baseDir + "MuHLTEff/MGnegative/eff.root";
  const TString zmmHLTEff2BinName_pos  = baseDir + "MuHLTEff/CTpositive/eff.root";
  const TString zmmHLTEff2BinName_neg  = baseDir + "MuHLTEff/CTnegative/eff.root";

  const TString dataSelEff2BinName_pos = baseDir + "MuSITEff/MGpositive/eff.root";
  const TString dataSelEff2BinName_neg = baseDir + "MuSITEff/MGnegative/eff.root";
  const TString zmmSelEff2BinName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmSelEff2BinName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataTrkEff2BinName_pos = baseDir + "MuSITEff/MGpositive/eff.root";
  const TString dataTrkEff2BinName_neg = baseDir + "MuSITEff/MGnegative/eff.root";
  const TString zmmTrkEff2BinName_pos  = baseDir + "MuSITEff/CTpositive/eff.root";
  const TString zmmTrkEff2BinName_neg  = baseDir + "MuSITEff/CTnegative/eff.root";

  const TString dataStaEff2BinName_pos = baseDir + "MuStaEff/MGpositive/eff.root";
  const TString dataStaEff2BinName_neg = baseDir + "MuStaEff/MGnegative/eff.root";
  const TString zmmStaEff2BinName_pos  = baseDir + "MuStaEff/CTpositive/eff.root";
  const TString zmmStaEff2BinName_neg  = baseDir + "MuStaEff/CTnegative/eff.root";

  TString StaEffSignalShapeSys     = baseDir + "Results/MuStaSigSys.root";
  TString StaEffBackgroundShapeSys = baseDir + "Results/MuStaBkgSys.root";
  TString SelEffSignalShapeSys     = baseDir + "Results/MuSITSigSys.root";
  TString SelEffBackgroundShapeSys = baseDir + "Results/MuSITBkgSys.root";

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);


  // plot output file format
//   const TString format("all");

  // setup efficiency shape systematics
  TFile *StaSigSysFile = new TFile(StaEffSignalShapeSys);
  TH2D *hStaSigSys = (TH2D*)StaSigSysFile->Get("h");
  TFile *StaBkgSysFile = new TFile(StaEffBackgroundShapeSys);
  TH2D *hStaBkgSys = (TH2D*)StaBkgSysFile->Get("h");
  TFile *SelSigSysFile = new TFile(SelEffSignalShapeSys);
  TH2D *hSelSigSys = (TH2D*)SelSigSysFile->Get("h");
  TFile *SelBkgSysFile = new TFile(SelEffBackgroundShapeSys);
  TH2D *hSelBkgSys = (TH2D*)SelBkgSysFile->Get("h");

  

  // file format for output plots
  const TString format("png"); 


  // ===================== Recoil correction files ============================
  const TString directory("/afs/cern.ch/user/d/dalfonso/public/WZ/nov26");
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
  
  // ==========================================================================
  
  // ----------------------------------------------------
  // Load the plots of relative difference, to create the up/down shapes
  TFile *_rdWmp = new TFile("shapeDiff/Wmp_relDiff.root");
  TFile *_rdWmm = new TFile("shapeDiff/Wmm_relDiff.root");
  
  TH1D *hh_diffm = new TH1D("hh_diffm","hh_diffm",75,0,150);
  TH1D *hh_diffp = new TH1D("hh_diffp","hh_diffp",75,0,150);
  
  hh_diffm = (TH1D*)_rdWmm->Get("hh_diff");
  hh_diffp = (TH1D*)_rdWmp->Get("hh_diff");
  // -----------------------------------------------------
  
  // Load the Z data and Z MC Pt spectra
  TFile *_rat1 = new TFile("shapeDiff/zmm_PDFUnc.root");
  TH1D *hh_mc;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_mc = (TH1D*)_rat1->Get("hZPtTruthNominal"); 
  hh_mc->Scale(1/hh_mc->Integral()); // normalize
  
  TFile *_rat2 = new TFile("shapeDiff/UnfoldingOutputZPt.root");
  TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_diff = (TH1D*)_rat2->Get("hUnfold");
  hh_diff->Scale(1/hh_diff->Integral()); // normalize
  hh_diff->Divide(hh_mc);
   
  //
  // input ntuple file names
  //
  enum { eData, eWmunu, eEWK, eBKG, eAntiData, eAntiWmunu, eAntiEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/wm_select.raw.root");   typev.push_back(eWmunu);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/boson_select.raw.root");  typev.push_back(eBKG);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/ewk_select1.root");  typev.push_back(eEWK);

  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWmunu/ntuples/data_select.root"); typev.push_back(eAntiData);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWmunu/ntuples/wm_select.root");   typev.push_back(eAntiWmunu);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWmunu/ntuples/ewk_select.root");  typev.push_back(eAntiEWK);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/AntiWmunu/ntuples/top_select.root");  typev.push_back(eAntiEWK);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet   = new TH1D("hDataMet","",  NBINS,0,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm  = new TH1D("hDataMetm","", NBINS,0,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp  = new TH1D("hDataMetp","", NBINS,0,METMAX); hDataMetp->Sumw2();
  TH1D *hWmunuMet  = new TH1D("hWmunuMet","", NBINS,0,METMAX); hWmunuMet->Sumw2();
  TH1D *hWmunuMetp = new TH1D("hWmunuMetp","",NBINS,0,METMAX); hWmunuMetp->Sumw2();
  TH1D *hWmunuMetm = new TH1D("hWmunuMetm","",NBINS,0,METMAX); hWmunuMetm->Sumw2();
  
  
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,0,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,0,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,0,METMAX); hEWKMetm->Sumw2();
  
  TH1D *hEWKMet_PileupUp   = new TH1D("hEWKMet_PileupUp",  "",NBINS,0,METMAX); hEWKMet_PileupUp->Sumw2();
  TH1D *hEWKMetp_PileupUp  = new TH1D("hEWKMetp_PileupUp", "",NBINS,0,METMAX); hEWKMetp_PileupUp->Sumw2();
  TH1D *hEWKMetm_PileupUp  = new TH1D("hEWKMetm_PileupUp", "",NBINS,0,METMAX); hEWKMetm_PileupUp->Sumw2();
  
  TH1D *hEWKMet_PileupDown   = new TH1D("hEWKMet_PileupDown",  "",NBINS,0,METMAX); hEWKMet_PileupDown->Sumw2();
  TH1D *hEWKMetp_PileupDown  = new TH1D("hEWKMetp_PileupDown", "",NBINS,0,METMAX); hEWKMetp_PileupDown->Sumw2();
  TH1D *hEWKMetm_PileupDown  = new TH1D("hEWKMetm_PileupDown", "",NBINS,0,METMAX); hEWKMetm_PileupDown->Sumw2();
  TH1D *hWmunuMet_RecoilUp  = new TH1D("hWmunuMet_RecoilUp", "",NBINS,0,METMAX); hWmunuMet_RecoilUp->Sumw2();
  TH1D *hWmunuMetp_RecoilUp = new TH1D("hWmunuMetp_RecoilUp","",NBINS,0,METMAX); hWmunuMetp_RecoilUp->Sumw2();
  TH1D *hWmunuMetm_RecoilUp = new TH1D("hWmunuMetm_RecoilUp","",NBINS,0,METMAX); hWmunuMetm_RecoilUp->Sumw2();
  TH1D *hWmunuMet_RecoilDown  = new TH1D("hWmunuMet_RecoilDown", "",NBINS,0,METMAX); hWmunuMet_RecoilDown->Sumw2();
  TH1D *hWmunuMetp_RecoilDown = new TH1D("hWmunuMetp_RecoilDown","",NBINS,0,METMAX); hWmunuMetp_RecoilDown->Sumw2();
  TH1D *hWmunuMetm_RecoilDown = new TH1D("hWmunuMetm_RecoilDown","",NBINS,0,METMAX); hWmunuMetm_RecoilDown->Sumw2();

  TH1D *hWmunuMet_RecoilCUp  = new TH1D("hWmunuMet_RecoilCUp", "",NBINS,0,METMAX); hWmunuMet_RecoilCUp->Sumw2();
  TH1D *hWmunuMetp_RecoilCUp = new TH1D("hWmunuMetp_RecoilCUp","",NBINS,0,METMAX); hWmunuMetp_RecoilCUp->Sumw2();
  TH1D *hWmunuMetm_RecoilCUp = new TH1D("hWmunuMetm_RecoilCUp","",NBINS,0,METMAX); hWmunuMetm_RecoilCUp->Sumw2();
  TH1D *hWmunuMet_RecoilCDown  = new TH1D("hWmunuMet_RecoilCDown", "",NBINS,0,METMAX); hWmunuMet_RecoilCDown->Sumw2();
  TH1D *hWmunuMetp_RecoilCDown = new TH1D("hWmunuMetp_RecoilCDown","",NBINS,0,METMAX); hWmunuMetp_RecoilCDown->Sumw2();
  TH1D *hWmunuMetm_RecoilCDown = new TH1D("hWmunuMetm_RecoilCDown","",NBINS,0,METMAX); hWmunuMetm_RecoilCDown->Sumw2();

  TH1D *hWmunuMet_ScaleUp  = new TH1D("hWmunuMet_ScaleUp", "",NBINS,0,METMAX); hWmunuMet_ScaleUp->Sumw2();
  TH1D *hWmunuMetp_ScaleUp = new TH1D("hWmunuMetp_ScaleUp","",NBINS,0,METMAX); hWmunuMetp_ScaleUp->Sumw2();
  TH1D *hWmunuMetm_ScaleUp = new TH1D("hWmunuMetm_ScaleUp","",NBINS,0,METMAX); hWmunuMetm_ScaleUp->Sumw2();
  TH1D *hWmunuMet_ScaleDown  = new TH1D("hWmunuMet_ScaleDown", "",NBINS,0,METMAX); hWmunuMet_ScaleDown->Sumw2();
  TH1D *hWmunuMetp_ScaleDown = new TH1D("hWmunuMetp_ScaleDown","",NBINS,0,METMAX); hWmunuMetp_ScaleDown->Sumw2();
  TH1D *hWmunuMetm_ScaleDown = new TH1D("hWmunuMetm_ScaleDown","",NBINS,0,METMAX); hWmunuMetm_ScaleDown->Sumw2();
  
  TH1D *hWmunuMet_PileupUp  = new TH1D("hWmunuMet_PileupUp", "",NBINS,0,METMAX); hWmunuMet_PileupUp->Sumw2();
  TH1D *hWmunuMetp_PileupUp = new TH1D("hWmunuMetp_PileupUp","",NBINS,0,METMAX); hWmunuMetp_PileupUp->Sumw2();
  TH1D *hWmunuMetm_PileupUp = new TH1D("hWmunuMetm_PileupUp","",NBINS,0,METMAX); hWmunuMetm_PileupUp->Sumw2();
  TH1D *hWmunuMet_PileupDown  = new TH1D("hWmunuMet_PileupDown", "",NBINS,0,METMAX); hWmunuMet_PileupDown->Sumw2();
  TH1D *hWmunuMetp_PileupDown = new TH1D("hWmunuMetp_PileupDown","",NBINS,0,METMAX); hWmunuMetp_PileupDown->Sumw2();
  TH1D *hWmunuMetm_PileupDown = new TH1D("hWmunuMetm_PileupDown","",NBINS,0,METMAX); hWmunuMetm_PileupDown->Sumw2();

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet","",  NBINS,0,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,0,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,0,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWmunuMet  = new TH1D("hAntiWmunuMet","", NBINS,0,METMAX); hAntiWmunuMet->Sumw2();
  TH1D *hAntiWmunuMetp = new TH1D("hAntiWmunuMetp","",NBINS,0,METMAX); hAntiWmunuMetp->Sumw2();
  TH1D *hAntiWmunuMetm = new TH1D("hAntiWmunuMetm","",NBINS,0,METMAX); hAntiWmunuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet", "",  NBINS,0,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,0,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,0,METMAX); hAntiEWKMetm->Sumw2();
  
//    //
//   // HLT efficiency
//   //
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
// 
//   
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
  TLorentzVector *lep=0, *genV=0;
  Float_t pfChIso, pfGamIso, pfNeuIso;
    
  rochcor2015 *rmcor = new rochcor2015();
  
  TFile *infile=0;
  TTree *intree=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);    assert(infile);
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
    intree->SetBranchAddress("genV",     &genV);       // lepton 4-vector
    intree->SetBranchAddress("pfChIso",  &pfChIso);
    intree->SetBranchAddress("pfGamIso", &pfGamIso);
    intree->SetBranchAddress("pfNeuIso", &pfNeuIso);
  
    Double_t mt=-999;

    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) std::cout << "On Entry.... " << ientry << std::endl;

      // vector containing raw lepton info for correcting MET
      TVector2 vLepRaw((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));

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

      if(fabs(lep->Eta()) > ETA_CUT) continue;
  
      mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );

      effdata=1; effmc=1;    
      if(q>0) {
        effdata *= (1.-dataHLTEff_pos.getEff((lep->Eta()), lep->Pt())); 
        effmc   *= (1.-zmmHLTEff_pos.getEff((lep->Eta()), lep->Pt())); 
      } else {
        effdata *= (1.-dataHLTEff_neg.getEff((lep->Eta()), lep->Pt())); 
        effmc   *= (1.-zmmHLTEff_neg.getEff((lep->Eta()), lep->Pt())); 
      }
      effdata = 1.-effdata;
      effmc   = 1.-effmc;
      corr *= effdata/effmc;
//       corrSigShape *= effdata/effmc;
//       corrBkgShape *= effdata/effmc;
    
      effdata=1; effmc=1;
      effSigShapedata=1;
      effBkgShapedata=1;
      if(q>0) {
        effdata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
//         effSigShapedata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(lep->Eta()), hSelSigSys->GetYaxis()->FindBin(lep->Pt())); 
//         effBkgShapedata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(lep->Eta()), hSelBkgSys->GetYaxis()->FindBin(lep->Pt())); 
      } else {
        effdata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
//         effSigShapedata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(lep->Eta()), hSelSigSys->GetYaxis()->FindBin(lep->Pt()));
//         effBkgShapedata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(lep->Eta()), hSelBkgSys->GetYaxis()->FindBin(lep->Pt())); 
      }
      corr *= effdata/effmc;
//       corrSigShape *= effSigShapedata/effmc;
//       corrBkgShape *= effBkgShapedata/effmc;
      
      effdata=1; effmc=1;
      effSigShapedata=1;
      effBkgShapedata=1;
      if(q>0) { 
        effdata *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
//         effSigShapedata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt())*hStaSigSys->GetBinContent(hStaSigSys->GetXaxis()->FindBin(lep->Eta()), hStaSigSys->GetYaxis()->FindBin(lep->Pt()));
//         effBkgShapedata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt())*hStaBkgSys->GetBinContent(hStaBkgSys->GetXaxis()->FindBin(lep->Eta()), hStaBkgSys->GetYaxis()->FindBin(lep->Pt())); 
      } else {
        effdata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
//         effSigShapedata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt())*hStaSigSys->GetBinContent(hStaSigSys->GetXaxis()->FindBin(lep->Eta()), hStaSigSys->GetYaxis()->FindBin(lep->Pt()));
//         effBkgShapedata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt())*hStaBkgSys->GetBinContent(hStaBkgSys->GetXaxis()->FindBin(lep->Eta()), hStaBkgSys->GetYaxis()->FindBin(lep->Pt())); 
      }
      corr *= effdata/effmc; 
//       corrSigShape *= effSigShapedata/effmc;
//       corrBkgShape *= effBkgShapedata/effmc;
      
      effdata=1; effmc=1;
      if(q>0) { 
        effdata *= dataTrkEff_pos.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmTrkEff_pos.getEff((lep->Eta()), lep->Pt()); 
      } else {
        effdata *= dataTrkEff_neg.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmTrkEff_neg.getEff((lep->Eta()), lep->Pt()); 
      }

      double var=0.;      
      
      // TRACKER
      if(q>0) {
        Double_t effdata = dataTrkEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataTrkEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataTrkEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmTrkEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmTrkEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmTrkEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        //var+=errTrk*errTrk;
      } else {
        Double_t effdata = dataTrkEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataTrkEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataTrkEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmTrkEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmTrkEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmTrkEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        //var+=errTrk*errTrk;
      }
//       std::cout << "made it standalone!" << std::endl;
      // STANDALONE
      if(q>0) {
        Double_t effdata = dataStaEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataStaEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmStaEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmStaEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSta*errSta;
      } else {
        Double_t effdata = dataStaEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataStaEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmStaEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmStaEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSta*errSta;
      }
//       std::cout << "made it sel!" << std::endl;
      // SELECTION
      if(q>0) {
        Double_t effdata = dataSelEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataSelEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmSelEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmSelEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSel*errSel;
      } else {
        Double_t effdata = dataSelEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataSelEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmSelEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmSelEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSel*errSel;
      }
// std::cout << "made it hlt!" << std::endl;
      //HLT
      if(q>0) {
        Double_t effdata = dataHLTEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataHLTEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmHLTEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmHLTEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errHLT*errHLT;
      } else {
        Double_t effdata = dataHLTEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataHLTEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmHLTEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmHLTEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errHLT*errHLT;
      }

      
      corrUp=corr+sqrt(var);
      corrDown=corr-sqrt(var);  
// std::cout << "made it 0!" << std::endl;
      eff2Bindata=1; eff2Binmc=1;    
          if(q>0) { 
            eff2Bindata *= (1.-dataHLTEff2Bin_pos.getEff((lep->Eta()), lep->Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_pos.getEff((lep->Eta()), lep->Pt())); 
          } else {
            eff2Bindata *= (1.-dataHLTEff2Bin_neg.getEff((lep->Eta()), lep->Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_neg.getEff((lep->Eta()), lep->Pt())); 
          }
          eff2Bindata = 1.-eff2Bindata;
          eff2Binmc   = 1.-eff2Binmc;
          corr2Bin *= eff2Bindata/eff2Binmc;
//     std::cout << "made it 1!" << std::endl;
          eff2Bindata=1; eff2Binmc=1;
          if(q>0) { 
            eff2Bindata *= dataSelEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
          } else {
            eff2Bindata *= dataSelEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
          }
          corr2Bin *= eff2Bindata/eff2Binmc;
//     std::cout << "made it 2!" << std::endl;
      eff2Bindata=1; eff2Binmc=1;
          if(q>0) { 
            eff2Bindata *= dataStaEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
          } else {
            eff2Bindata *= dataStaEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
          }
      corr2Bin *= eff2Bindata/eff2Binmc; 
//       std::cout << "made it 3!" << std::endl;
          eff2Bindata=1; eff2Binmc=1;
          if(q>0) { 
            eff2Bindata *= dataTrkEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmTrkEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
          } else {
            eff2Bindata *= dataTrkEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmTrkEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
          }
      if(typev[ifile]==eData || typev[ifile]==eAntiData){
        // Apply the Rochester Corrections to data (anti-isolation)
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;
        rmcor->momcor_data(mu1,q,0,qter1);
        if(mu1.Pt()        < PT_CUT)  continue;
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
        // calculate the corrected MET
        TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
        Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
        if(typev[ifile]==eData) {
          hDataMet->Fill(corrMetWithLepton);
          if(q>0) { hDataMetp->Fill(corrMetWithLepton); }
          else    { hDataMetm->Fill(corrMetWithLepton); }
        } else if(typev[ifile]==eAntiData) {
          hAntiDataMet->Fill(corrMetWithLepton);
          if(q>0) { hAntiDataMetp->Fill(corrMetWithLepton); } 
          else    { hAntiDataMetm->Fill(corrMetWithLepton); }   
        }
      } else {
        Double_t weight = 1;Double_t weightUp = 1;Double_t weightDown = 1;
        Double_t weight2 =1;
        //corr = 1.0;

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
        
        
        // Prepare 2-d Vector for raw Muon

        // Do some Rochester corrections for MC
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;
        rmcor->momcor_mc(mu1,q,0,qter1);
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
        Double_t lepPt = mu1.Pt();
        // change to have rochester corrected muon and raw lepton with MET corrected same way as electron channel 
        if(typev[ifile]==eWmunu || typev[ifile]==eBKG) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          if(lep->Pt()        > PT_CUT) {
            double bin = 0;
            for(int i = 1; i <= hh_diff->GetNbinsX();++i){
              if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            }
            double w2 = 1.0;//hh_diff->GetBinContent(bin);
            
            corrMet=met, corrMetPhi=metPhi;
            hWmunuMet->Fill(corrMet,weight);
//             hWmunuMet_PileupUp->Fill(corrMet,weightUp);
//             hWmunuMet_PileupDown->Fill(corrMet,weightDown);
            if(q>0) {
              if(doInclusive) recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
              else if(pileupUp)    recoilCorrPuUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
              else if(pileupDown)  recoilCorrPuDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
              else if(doKeys){
                if(fabs(genVy)<0.5)
                  recoilCorrKeys05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                  recoilCorrKeys051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else
                  recoilCorrKeys1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
              } else {
                if(fabs(genVy)<0.5)
                  recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                  recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else
                  recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
              }
              // Compute the corrected MET value
              TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
              Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
              if(typev[ifile]==eWmunu)
              {
                hWmunuMetp->Fill(corrMetWithLepton,weight);
//                 hWmunuMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
//                 hWmunuMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
              }
              else
              {
                hEWKMetp->Fill(corrMetWithLepton,weight);
              }

              corrMet=met, corrMetPhi=metPhi;
            } else {
              if(doInclusive) recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
              else if(pileupUp)    recoilCorrPuUpm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
              else if(pileupDown)  recoilCorrPuDownm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
              else if(doKeys){
                if(fabs(genVy)<0.5)
                  recoilCorrKeysm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                  recoilCorrKeysm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else
                  recoilCorrKeysm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
              } else {
                if(fabs(genVy)<0.5)
                  recoilCorrm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                  recoilCorrm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                else
                  recoilCorrm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
              }
              // Compute the corrected MET value
              TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
              Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
              if(typev[ifile]==eWmunu)
              {
                hWmunuMetm->Fill(corrMetWithLepton,weight);
//                 hWmunuMetm_PileupUp->Fill(corrMetWithLepton,weightUp);
//                 hWmunuMetm_PileupDown->Fill(corrMetWithLepton,weightDown);
              }
              else
              {
                hEWKMetm->Fill(corrMetWithLepton,weight);
              }
              corrMet=met, corrMetPhi=metPhi;
            }
            corrMet=met, corrMetPhi=metPhi;
          }
          //             recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),1,q);
          // unused?
          hWmunuMet_RecoilUp->Fill(corrMet,weight);
          if(q>0) {
            pU1 = 0; pU2 = 0; 
            hWmunuMetp_RecoilUp->Fill(corrMet,weight); 
            corrMet=met, corrMetPhi=metPhi;
          } else { 
            pU1 = 0; pU2 = 0; 
            hWmunuMetm_RecoilUp->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          }
          hWmunuMet_RecoilDown->Fill(corrMet,weight);
          if(q>0) {
            pU1 = 0; pU2 = 0; 
            hWmunuMetp_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          } else {
            pU1 = 0; pU2 = 0; 
            hWmunuMetm_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          }
        }
        if(typev[ifile]==eAntiWmunu) {
          if(lep->Pt()        < PT_CUT)  continue;
          Double_t corrMet=met, corrMetPhi=metPhi;
//           TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
//           Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
//        this histogram isn't actually used in any results / fits
          hAntiWmunuMet->Fill(corrMet,weight2);
          if(q>0) {              
            pU1 = 0; pU2 = 0; 
            if(doInclusive) recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupUp)    recoilCorrPuUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupDown)  recoilCorrPuDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
            else if(doKeys){
              if(fabs(genVy)<0.5)
                recoilCorrKeys05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorrKeys051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorrKeys1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            } else {
              if(fabs(genVy)<0.5)
                recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            }
            //recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
            Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            hAntiWmunuMetp->Fill(corrMetWithLepton,weight2); 
            corrMet = met; corrMetPhi = metPhi;
          } 
          else {
            pU1 = 0; pU2 = 0; 
            if(doInclusive) recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupUp)    recoilCorrPuUpm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
            else if(pileupDown)  recoilCorrPuDownm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0);
            else if(doKeys){
              if(fabs(genVy)<0.5)
                recoilCorrKeysm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorrKeysm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorrKeysm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            } else {
              if(fabs(genVy)<0.5)
                recoilCorrm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                recoilCorrm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              else
                recoilCorrm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            }
            //recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
            Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            hAntiWmunuMetm->Fill(corrMetWithLepton,weight2);
            corrMet = met; corrMetPhi = metPhi; 
          }
        }
        if(typev[ifile]==eEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          hEWKMet->Fill(corrMetWithLepton,weight);
//           hEWKMet_PileupUp->Fill(corrMetWithLepton,weightUp);
//           hEWKMet_PileupDown->Fill(corrMetWithLepton,weightDown);
          if(q>0) {
            hEWKMetp->Fill(corrMetWithLepton,weight); 
            // RooFit doesn't see these histograms, only for combine, I will remove completely later
//             hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp); 
//             hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown); 
          }
          else { 
            hEWKMetm->Fill(corrMetWithLepton,weight); 
//             hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); 
//             hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); 
          }
        }
        if(typev[ifile]==eAntiEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          hAntiEWKMet->Fill(corrMetWithLepton,weight2);
          if(q>0) { hAntiEWKMetp->Fill(corrMetWithLepton,weight2); }
          else    { hAntiEWKMetm->Fill(corrMetWithLepton,weight2); }
        }
      }
    }
  }
  delete infile;
  infile=0, intree=0;   
 
  // Here rescale the up/down whatever
//   Calculate the shapes for W+
  TH1D *corrP = (TH1D*) hWmunuMetp->Clone("up");
  *corrP = (*hh_diffp)*(*hWmunuMetp);
  hWmunuMetp_RecoilUp->Add(corrP,hWmunuMetp,-1);
  hWmunuMetp_RecoilDown->Add(corrP,hWmunuMetp,1);
//   hWmunuMetp_ScaleUp->Add(corrP,hWmunuMetp,-1);
//   hWmunuMetp_ScaleDown->Add(corrP,hWmunuMetp,1);
  // Calculate the shapes for W-
  TH1D *corrM = (TH1D*) hWmunuMetm->Clone("up");
  *corrM = (*hh_diffm)*(*hWmunuMetm);
  hWmunuMetm_RecoilUp->Add(corrM,hWmunuMetm,-1);
  hWmunuMetm_RecoilDown->Add(corrM,hWmunuMetm,1);
//   hWmunuMetm_ScaleUp->Add(corrM,hWmunuMetm,-1);
//   hWmunuMetm_ScaleDown->Add(corrM,hWmunuMetm,1);
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWmunuMet->Integral());
//   cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWmunuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
//   dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",1.0*(hWmunuMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral(),0,hAntiDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",hDataMetp->Integral()*0.3,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWmunuMetp->Integral());
//   cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  //RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral()*1.0,0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
//   dewkp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",1.0*(hWmunuMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm->Integral()*0.3,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWmunuMetm->Integral());
//   cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWmunuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
//   dewkm.setConstant(kTRUE);
  RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));
  
  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,wmunuMetm,1); 
  RooDataHist wmunuMet_RecoilUp("wmunuMET_RecoilUp", "wmunuMET_RecoilUp", RooArgSet(pfmet),hWmunuMet_RecoilUp);  RooHistPdf pdfWm_RecoilUp("wm_RecoilUp", "wm_RecoilUp", pfmet,wmunuMet_RecoilUp, 1);
  RooDataHist wmunuMetp_RecoilUp("wmunuMETp_RecoilUp","wmunuMETp_RecoilUp",RooArgSet(pfmet),hWmunuMetp_RecoilUp); RooHistPdf pdfWmp_RecoilUp("wmp_RecoilUp","wmp_RecoilUp",pfmet,wmunuMetp_RecoilUp,1);
  RooDataHist wmunuMetm_RecoilUp("wmunuMETm_RecoilUp","wmunuMETm_RecoilUp",RooArgSet(pfmet),hWmunuMetm_RecoilUp); RooHistPdf pdfWmm_RecoilUp("wmm_RecoilUp","wmm_RecoilUp",pfmet,wmunuMetm_RecoilUp,1); 
  RooDataHist wmunuMet_RecoilDown("wmunuMET_RecoilDown", "wmunuMET_RecoilDown", RooArgSet(pfmet),hWmunuMet_RecoilDown);  RooHistPdf pdfWm_RecoilDown("wm_RecoilDown", "wm_RecoilDown", pfmet,wmunuMet_RecoilDown, 1);
  RooDataHist wmunuMetp_RecoilDown("wmunuMETp_RecoilDown","wmunuMETp_RecoilDown",RooArgSet(pfmet),hWmunuMetp_RecoilDown); RooHistPdf pdfWmp_RecoilDown("wmp_RecoilDown","wmp_RecoilDown",pfmet,wmunuMetp_RecoilDown,1);
  RooDataHist wmunuMetm_RecoilDown("wmunuMETm_RecoilDown","wmunuMETm_RecoilDown",RooArgSet(pfmet),hWmunuMetm_RecoilDown); RooHistPdf pdfWmm_RecoilDown("wmm_RecoilDown","wmm_RecoilDown",pfmet,wmunuMetm_RecoilDown,1); 
   RooDataHist wmunuMet_ScaleUp("wmunuMET_ScaleUp", "wmunuMET_ScaleUp", RooArgSet(pfmet),hWmunuMet_ScaleUp);  RooHistPdf pdfWm_ScaleUp("wm_ScaleUp", "wm_ScaleUp", pfmet,wmunuMet_ScaleUp, 1);
  RooDataHist wmunuMetp_ScaleUp("wmunuMETp_ScaleUp","wmunuMETp_ScaleUp",RooArgSet(pfmet),hWmunuMetp_ScaleUp); RooHistPdf pdfWmp_ScaleUp("wmp_ScaleUp","wmp_ScaleUp",pfmet,wmunuMetp_ScaleUp,1);
  RooDataHist wmunuMetm_ScaleUp("wmunuMETm_ScaleUp","wmunuMETm_ScaleUp",RooArgSet(pfmet),hWmunuMetm_ScaleUp); RooHistPdf pdfWmm_ScaleUp("wmm_ScaleUp","wmm_ScaleUp",pfmet,wmunuMetm_ScaleUp,1); 
  RooDataHist wmunuMet_ScaleDown("wmunuMET_ScaleDown", "wmunuMET_ScaleDown", RooArgSet(pfmet),hWmunuMet_ScaleDown);  RooHistPdf pdfWm_ScaleDown("wm_ScaleDown", "wm_ScaleDown", pfmet,wmunuMet_ScaleDown, 1);
  RooDataHist wmunuMetp_ScaleDown("wmunuMETp_ScaleDown","wmunuMETp_ScaleDown",RooArgSet(pfmet),hWmunuMetp_ScaleDown); RooHistPdf pdfWmp_ScaleDown("wmp_ScaleDown","wmp_ScaleDown",pfmet,wmunuMetp_ScaleDown,1);
  RooDataHist wmunuMetm_ScaleDown("wmunuMETm_ScaleDown","wmunuMETm_ScaleDown",RooArgSet(pfmet),hWmunuMetm_ScaleDown); RooHistPdf pdfWmm_ScaleDown("wmm_ScaleDown","wmm_ScaleDown",pfmet,wmunuMetm_ScaleDown,1); 
  
    RooDataHist wmunuMet_PileupUp("wmunuMET_PileupUp", "wmunuMET_PileupUp", RooArgSet(pfmet),hWmunuMet_PileupUp);  RooHistPdf pdfWm_PileupUp("wm_PileupUp", "wm_PileupUp", pfmet,wmunuMet_PileupUp, 1);
  RooDataHist wmunuMetp_PileupUp("wmunuMETp_PileupUp","wmunuMETp_PileupUp",RooArgSet(pfmet),hWmunuMetp_PileupUp); RooHistPdf pdfWmp_PileupUp("wmp_PileupUp","wmp_PileupUp",pfmet,wmunuMetp_PileupUp,1);
  RooDataHist wmunuMetm_PileupUp("wmunuMETm_PileupUp","wmunuMETm_PileupUp",RooArgSet(pfmet),hWmunuMetm_PileupUp); RooHistPdf pdfWmm_PileupUp("wmm_PileupUp","wmm_PileupUp",pfmet,wmunuMetm_PileupUp,1); 
  RooDataHist wmunuMet_PileupDown("wmunuMET_PileupDown", "wmunuMET_PileupDown", RooArgSet(pfmet),hWmunuMet_PileupDown);  RooHistPdf pdfWm_PileupDown("wm_PileupDown", "wm_PileupDown", pfmet,wmunuMet_PileupDown, 1);
  RooDataHist wmunuMetp_PileupDown("wmunuMETp_PileupDown","wmunuMETp_PileupDown",RooArgSet(pfmet),hWmunuMetp_PileupDown); RooHistPdf pdfWmp_PileupDown("wmp_PileupDown","wmp_PileupDown",pfmet,wmunuMetp_PileupDown,1);
  RooDataHist wmunuMetm_PileupDown("wmunuMETm_PileupDown","wmunuMETm_PileupDown",RooArgSet(pfmet),hWmunuMetm_PileupDown); RooHistPdf pdfWmm_PileupDown("wmm_PileupDown","wmm_PileupDown",pfmet,wmunuMetm_PileupDown,1); 
  
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
  
//   // test using the reversed dEta, dPhi cuts as background
//   RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
//   RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
//   RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  // QCD Pdfs
  
  //CExponential qcd(pfmet,kTRUE);
  //CExponential qcdp(pfmet,kTRUE);
  //CExponential qcdm(pfmet,kTRUE);
  // comment back in for qcd functional form
  CPepeModel2 qcd("qcd",pfmet);
  CPepeModel2 qcdp("qcdp",pfmet);
  CPepeModel2 qcdm("qcdm",pfmet);
  //CPepeModel2 qcdm("qcdm",pfmet);
  //CPepeModel1 qcdm("qcdm",pfmet);
//  qcdp.a1->setConstant(kTRUE);
//   qcdp.a1->setVal(0.24);
//    qcdm.a1->setVal(0.24);

//   RooRealVar a1ConstMeanP("a1ConstMeanP","a1ConstMeanP",0.18);
//   RooRealVar a1ConstSigmaP("a1ConstSigmaP","a1ConstSigmaP",0.003);
// //   RooRealVar a1ConstMeanM("a1ConstMeanM","a1ConstMeanM",0.215);
// //   RooRealVar a1ConstSigmaM("a1ConstSigmaM","a1ConstSigmaM",0.007);
// 
//   RooRealVar sigConstMeanP("sigConstMeanP","sigConstMeanP",17.6);
//   RooRealVar sigConstSigmaP("sigConstSigmaP","sigConstSigmaP",0.2);
//   RooRealVar sigConstMeanM("sigConstMeanM","sigConstMeanM",16.9);
//   RooRealVar sigConstSigmaM("sigConstSigmaM","sigConstSigmaM",0.5);
//   
// // RooRealVar f("f","f",0.5,0.,1.) ;
//   RooGaussian fconsta1p("fconsta1p","fconsta1p",*(qcdp.a1),a1ConstMeanP,a1ConstSigmaP);
// //   RooGaussian fconsta1m("fconsta1m","fconsta1m",*(qcdm.a1),a1ConstMeanM,a1ConstSigmaM);
//   
//   RooGaussian fconstsigp("fconstsigp","fconstsigp",*(qcdp.sigma),sigConstMeanP,sigConstSigmaP);
//   RooGaussian fconstsigm("fconstsigm","fconstsigm",*(qcdm.sigma),sigConstMeanM,sigConstSigmaM);
//  
//   RooProdPdf qcdpc("qcdpc","qcdpc",RooArgSet(*(qcdp.model),fconsta1p,fconstsigp));
//   RooProdPdf qcdmc("qcdmc","qcdmc",RooArgSet(*(qcdm.model),fconstsigm));
//   RooProdPdf qcdmc("qcdmc","qcdmc",RooArgSet(qcdp,fconstraint));

    RooGaussian constm("constm","constm",nEWKm,RooConst(hEWKMetm->Integral()),RooConst(0.15*hEWKMetm->Integral()));
    RooGaussian constp("constp","constp",nEWKp,RooConst(hEWKMetp->Integral()),RooConst(0.15*hEWKMetp->Integral()));
 
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));
  
//   // constrained background?
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,qcdpc),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,qcdmc),RooArgList(nSigm,nEWKm,nQCDm));
  
  
    
  
  // Anti-Signal PDFs
  RooDataHist awmunuMet ("awmunuMET", "awmunuMET", RooArgSet(pfmet),hAntiWmunuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,awmunuMet, 1);
  RooDataHist awmunuMetp("awmunuMETp","awmunuMETp",RooArgSet(pfmet),hAntiWmunuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,awmunuMetp,1);
  RooDataHist awmunuMetm("awmunuMETm","awmunuMETm",RooArgSet(pfmet),hAntiWmunuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,awmunuMetm,1); 
  
  // Anti-EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  
//   // Anti-QCD Pdfs
  CPepeModel2 aqcd("aqcd",pfmet,qcd.a1);
  CPepeModel2 aqcdp("aqcdp",pfmet,qcdp.a1);
  CPepeModel2 aqcdm("aqcdm",pfmet,qcdm.a1);

//   CPepeModel2 aqcd("aqcd",pfmet);
//   CPepeModel2 aqcdp("aqcdp",pfmet);
//   CPepeModel2 aqcdm("aqcdm",pfmet);
  
  // Anti-selection PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
  // PDF for simultaneous fit
  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("Select");
  rooCat.defineType("Anti");
  
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(pdfMet, "Select");
  //pdfTotal.addPdf(apdfMet,"Anti");
  
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

//   RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
// 			 Import("Select", dataMetp),
// 			 Import("Anti",   antiMetp));
// 
//   RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
// 			 Import("Select", dataMetm),
// 			 Import("Anti", antiMetm));
  
  cout << "Starting values for Wmunu yields: " << endl;
  cout << "Selected: " << hDataMet->Integral() << endl;
  cout << "   sig: " << hWmunuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;
  cout << "   qcd: " << hDataMet->Integral()-hWmunuMet->Integral()-hEWKMet->Integral() << endl;

  cout << "Starting values for Wmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hDataMetp->Integral()-hWmunuMetp->Integral()-hEWKMetp->Integral() << endl;

  cout << "Starting values for Wmunu_m yields: " << endl;
  cout << "   sig: " << hWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  cout << "   qcd: " << hDataMetm->Integral()-hWmunuMetm->Integral()-hEWKMetm->Integral() << endl;
  
  
  cout << "Starting values for AntiWmunu yields: " << endl;
  cout << "Selected: " << hAntiDataMet->Integral() << endl;
  cout << "   sig: " << hAntiWmunuMet->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMet->Integral() << endl;
  cout << "   qcd: " << hAntiDataMet->Integral()-hAntiWmunuMet->Integral()-hAntiEWKMet->Integral() << endl;

  cout << "Starting values for AntiWmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetp->Integral()-hAntiWmunuMetp->Integral()-hAntiEWKMetp->Integral() << endl;

  cout << "Starting values for AntiWmunu_m yields: " << endl;
  cout << "   sig: " << hAntiWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetm->Integral()-hAntiWmunuMetm->Integral()-hAntiEWKMetm->Integral() << endl;

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());

  RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  
  RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  
  
//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",40000,0,100000);
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",40000,0,100000);

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);
  combine_workspace.import(pepe2Pdf_qcdp_norm);
  combine_workspace.import(pepe2Pdf_qcdm_norm);

  combine_workspace.import(pdfWm);
  combine_workspace.import(pdfWmp);
  combine_workspace.import(pdfWmm);
  combine_workspace.import(pdfWm_RecoilUp);
  combine_workspace.import(pdfWmp_RecoilUp);
  combine_workspace.import(pdfWmm_RecoilUp);
  combine_workspace.import(pdfWm_RecoilDown);
  combine_workspace.import(pdfWmp_RecoilDown);
  combine_workspace.import(pdfWmm_RecoilDown);
  combine_workspace.import(pdfWm_ScaleUp);
  combine_workspace.import(pdfWmp_ScaleUp);
  combine_workspace.import(pdfWmm_ScaleUp);
  combine_workspace.import(pdfWm_ScaleDown);
  combine_workspace.import(pdfWmp_ScaleDown);
  combine_workspace.import(pdfWmm_ScaleDown);
  combine_workspace.import(pdfWm_PileupUp);
  combine_workspace.import(pdfWmp_PileupUp);
  combine_workspace.import(pdfWmm_PileupUp);
  combine_workspace.import(pdfWm_PileupDown);
  combine_workspace.import(pdfWmp_PileupDown);
  combine_workspace.import(pdfWmm_PileupDown);
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
  
  combine_workspace.import(apdfWm);
  combine_workspace.import(apdfWmp);
  combine_workspace.import(apdfWmm);
  combine_workspace.import(apdfEWK);
  combine_workspace.import(apdfEWKp);
  combine_workspace.import(apdfEWKm);
  combine_workspace.import(*(aqcd.model));
  //combine_workspace.import(qcdpn);
  combine_workspace.import(*(aqcdp.model));
  combine_workspace.import(*(aqcdm.model));

  combine_workspace.writeToFile("Wmunu_pdfTemplates.root");

  RooDataHist dataTotal("dataTotal","dataTotal", RooArgList(pfmet), Index(rooCat),
            Import("Select", dataMet),
            Import("Anti",   antiMet));
  RooFitResult *fitRes = pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));//dataTotal.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));
  
 
  RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
             Import("Select", dataMetp),
             Import("Anti",   antiMetp));
             
//   RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntip = apdfMetp.fitTo(antiMetp,Extended(),Minos(kTRUE),Save(kTRUE));
  RooFitResult *fitResp = pdfTotalp.fitTo(dataTotalp,Extended(),ExternalConstraints(constp),RooFit::Strategy(2),Minos(kTRUE),Save(kTRUE));
  
//   RooDataHist dataMetp2("dataMetp2", "dataMetp2", RooArgSet(pfmet), hDataMetp);
  RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
             Import("Select", dataMetm),
             Import("Anti", antiMetm));
//   RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),ExternalConstraints(constm),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
//   RooFitResult *fitResAntim = apdfMetm.fitTo(antiMetm,Extended(),Minos(kTRUE),Save(kTRUE));
RooFitResult *fitResm = pdfTotalm.fitTo(dataTotalm,Extended(),ExternalConstraints(constm),Minimizer("MIGRAD","minimize"),RooFit::Strategy(2),Minos(kTRUE),Save(kTRUE));
  
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
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
  hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
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
  gStyle->SetTitleOffset(1.100,"Y");
  TGaxis::SetMaxDigits(3);
  
  char ylabel[100];  // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 8 TeV",lumi*1000.);
  else         sprintf(lumitext,"%.1f fb^{-1}  at  #sqrt{s} = 13 TeV",lumi/1000.);
  
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
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("wmunu_fitmet",weframe,"","mT [GeV]",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox("CMS",0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox(lumitext,0.63,0.92,0.95,0.99,0);
  plotMet.Draw(c,kTRUE,format,1);

  CPlot plotMetDiff("wmunu_fitmet","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-0.20,0.20);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  plotMetDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMet.SetName("wmunu_fitmetlog");
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
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("wmunu_fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("wmunu_fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.20,0.20);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("wmunu_fitmetplog");
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
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfWmp)),LineColor(linecolorW),LineStyle(2));
  antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
  CPlot plotAntiMetp("wmunu_fitantimetp",awepframe,"","",ylabel);
  plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotAntiMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetp.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetpDiff("wmunu_fitantimetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
  plotAntiMetpDiff.SetYRange(-0.20,0.20);
  plotAntiMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotAntiMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotAntiMetpDiff.Draw(c,kTRUE,format,2);
  plotAntiMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotAntiMetp.SetName("wmunu_fitantimetplog");
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
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("wmunu_fitmetm",wemframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("wmunu_fitmetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetmDiff.SetYRange(-0.2,0.2);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("wmunu_fitmetmlog");
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
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfWmm)),LineColor(linecolorW),LineStyle(2));
  antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetm->GetBinWidth(1));
  CPlot plotAntiMetm("wmunu_fitantimetm",awemframe,"","",ylabel);
  plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotAntiMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetm.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetmDiff("wmunu_fitantimetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotAntiMetmDiff.SetYRange(-0.2,0.2);
  plotAntiMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotAntiMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotAntiMetmDiff.AddHist1D(hAntiMetmDiff,"EX0",ratioColor);
  plotAntiMetmDiff.Draw(c,kTRUE,format,2);
  plotAntiMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotAntiMetm.SetName("wmunu_fitantimetmlog");
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
  sprintf(txtfname,"%s/fitresWm.txt",CPlot::sOutDir.Data());
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
  sprintf(txtfname,"%s/fitresWmp.txt",CPlot::sOutDir.Data());
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
  sprintf(txtfname,"%s/fitresWmm.txt",CPlot::sOutDir.Data());
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
  
  gBenchmark->Show("fitWm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
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
  sprintf(htmlfname,"%s/WmunuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wmunu</title></head>" << endl;
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
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimet.png\"><img src=\"fitantimet.png\" alt=\"fitantimet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetp.png\"><img src=\"fitantimetp.png\" alt=\"fitantimetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetm.png\"><img src=\"fitantimetm.png\" alt=\"fitantimetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetlog.png\"><img src=\"fitantimetlog.png\" alt=\"fitantimetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetplog.png\"><img src=\"fitantimetplog.png\" alt=\"fitantimetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetmlog.png\"><img src=\"fitantimetmlog.png\" alt=\"fitantimetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}
