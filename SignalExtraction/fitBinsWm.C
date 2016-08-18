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

void fitBinsWm(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)'
           const Double_t lumi2, // lumi for the anti-isolation trigger
       const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS   = 75;
  const Double_t METMAX  = 150;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t mu_MASS = 0.1057;
  const int NTOYS = 100;
  
    // efficiency files
 
  const TString baseDir = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/OldMuStore/25ns76X_SMP15011/"; 
  const TString dataHLTEffName_pos = baseDir + "MuHLTEff/MG/eff.root";
  const TString dataHLTEffName_neg = baseDir + "MuHLTEff/MG/eff.root";
  const TString zmmHLTEffName_pos  = baseDir + "MuHLTEff/CT/eff.root";
  const TString zmmHLTEffName_neg  = baseDir + "MuHLTEff/CT/eff.root";

  const TString dataSelEffName_pos = baseDir + "MuSITEff/MG/eff.root";
  const TString dataSelEffName_neg = baseDir + "MuSITEff/MG/eff.root";
  const TString zmmSelEffName_pos  = baseDir + "MuSITEff/CT/eff.root";
  const TString zmmSelEffName_neg  = baseDir + "MuSITEff/CT/eff.root";

  const TString dataTrkEffName_pos = baseDir + "MuSITEff/MG/eff.root";
  const TString dataTrkEffName_neg = baseDir + "MuSITEff/MG/eff.root";
  const TString zmmTrkEffName_pos  = baseDir + "MuSITEff/CT/eff.root";
  const TString zmmTrkEffName_neg  = baseDir + "MuSITEff/CT/eff.root";

  const TString dataStaEffName_pos = baseDir + "MuStaEff/MG/eff.root";
  const TString dataStaEffName_neg = baseDir + "MuStaEff/MG/eff.root";
  const TString zmmStaEffName_pos  = baseDir + "MuStaEff/CT/eff.root";
  const TString zmmStaEffName_neg  = baseDir + "MuStaEff/CT/eff.root";

  // efficiency files 2Bins

  const TString dataHLTEff2BinName_pos = baseDir + "MuHLTEff/1MG/eff.root";
  const TString dataHLTEff2BinName_neg = baseDir + "MuHLTEff/1MG/eff.root";
  const TString zmmHLTEff2BinName_pos  = baseDir + "MuHLTEff/1CT/eff.root";
  const TString zmmHLTEff2BinName_neg  = baseDir + "MuHLTEff/1CT/eff.root";

  const TString dataSelEff2BinName_pos = baseDir + "MuSITEff/1MG/eff.root";
  const TString dataSelEff2BinName_neg = baseDir + "MuSITEff/1MG/eff.root";
  const TString zmmSelEff2BinName_pos  = baseDir + "MuSITEff/1CT/eff.root";
  const TString zmmSelEff2BinName_neg  = baseDir + "MuSITEff/1CT/eff.root";

  const TString dataTrkEff2BinName_pos = baseDir + "MuSITEff/1MG/eff.root";
  const TString dataTrkEff2BinName_neg = baseDir + "MuSITEff/1MG/eff.root";
  const TString zmmTrkEff2BinName_pos  = baseDir + "MuSITEff/1CT/eff.root";
  const TString zmmTrkEff2BinName_neg  = baseDir + "MuSITEff/1CT/eff.root";

  const TString dataStaEff2BinName_pos = baseDir + "MuStaEff/1MG/eff.root";
  const TString dataStaEff2BinName_neg = baseDir + "MuStaEff/1MG/eff.root";
  const TString zmmStaEff2BinName_pos  = baseDir + "MuStaEff/1CT/eff.root";
  const TString zmmStaEff2BinName_neg  = baseDir + "MuStaEff/1CT/eff.root";

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

  
  Double_t vWPtBins[] = {0,0.20,0.40,0.60,0.80,1.0,1.2,1.4,1.6,1.85,2.10,2.40};
  Int_t nWBins = sizeof(vWPtBins)/sizeof(Double_t)-1;

  // file format for output plots
  const TString format("png"); 


  // for Puppi
  RecoilCorrector *recoilCorr = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorr->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorr->loadRooWorkspacesMC("../Recoil/WmpMC_newBacon/");
  
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorrm->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorrm->loadRooWorkspacesMC("../Recoil/WmmMC_newBacon/");
  
  // pileup up
  RecoilCorrector *recoilCorrUp = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorrUp->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorrUp->loadRooWorkspacesMC("../Recoil/WmpMCPuppi_PileupUp/");
  
  RecoilCorrector *recoilCorrmUp = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorrmUp->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorrmUp->loadRooWorkspacesMC("../Recoil/WmmMCPuppi_PileupUp/");
  
  // pileup down
  RecoilCorrector *recoilCorrDown = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorrDown->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorrDown->loadRooWorkspacesMC("../Recoil/WmpMCPuppi_PileupDown/");
  
  RecoilCorrector *recoilCorrmDown = new  RecoilCorrector("../Recoil/ZmmMCPuppi_newBacon/fits_puppi.root","fcnPF"); // get tgraph from here? what i guess its mean so it doesn't matter?
  recoilCorrmDown->loadRooWorkspacesData("../Recoil/ZmmDataPuppi_newBacon/");
  recoilCorrmDown->loadRooWorkspacesMC("../Recoil/WmmMCPuppi_PileupDown/");
  
 
  
  // ----------------------------------------------------
  // Load the plots of relative difference, to create the up/down shapes
  TFile *_rdWmp = new TFile("shapeDiff_etaBins/Wmp_relDiff.root");
  TFile *_rdWmm = new TFile("shapeDiff_etaBins/Wmm_relDiff.root");
  
  //TH1D *hh_diffm = new TH1D("hh_diffm","hh_diffm",75,0,150);
  //TH1D *hh_diffp = new TH1D("hh_diffp","hh_diffp",75,0,150);
  
 // hh_diffm = (TH1D*)_rdWmm->Get("hh_diff");
 // hh_diffp = (TH1D*)_rdWmp->Get("hh_diff");
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
  enum { eData, eWmunu, eEWK, eQCD, eAntiData, eAntiWmunu, eAntiEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/wm_select.raw.root");   typev.push_back(eWmunu);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/ewk_select.root");  typev.push_back(eEWK);
  fnamev.push_back("/afs/cern.ch/work/a/arapyan/public/flat_ntuples/Wmunu/ntuples/top_select.raw.root");  typev.push_back(eEWK);

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
  
  
  
  TH2D *hDataMet2D   = new TH2D("hDataMet2D","",  NBINS,0,METMAX,nWBins,vWPtBins); //hDataMet->Sumw2();
  TH2D *hDataMetm2D  = new TH2D("hDataMetm2D","", NBINS,0,METMAX,nWBins,vWPtBins); //hDataMetm->Sumw2();  
  TH2D *hDataMetp2D  = new TH2D("hDataMetp2D","", NBINS,0,METMAX,nWBins,vWPtBins); //hDataMetp->Sumw2();
  TH2D *hWmunuMet2D  = new TH2D("hWmunuMet2D","", NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMet->Sumw2();
  TH2D *hWmunuMetp2D = new TH2D("hWmunuMetp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hWmunuMetm2D = new TH2D("hWmunuMetm2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  TH2D *hWmunuMetp_RecoilUp2D = new TH2D("hWmunuMetp_RecoilUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hWmunuMetm_RecoilUp2D = new TH2D("hWmunuMetm_RecoilUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  TH2D *hWmunuMetp_RecoilDown2D = new TH2D("hWmunuMetp_RecoilDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hWmunuMetm_RecoilDown2D = new TH2D("hWmunuMetm_RecoilDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  
  TH2D *hWmunuMetp_PileupUp2D = new TH2D("hWmunuMetp_PileupUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hWmunuMetm_PileupUp2D = new TH2D("hWmunuMetm_PileupUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  TH2D *hWmunuMetp_PileupDown2D = new TH2D("hWmunuMetp_PileupDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hWmunuMetm_PileupDown2D = new TH2D("hWmunuMetm_PileupDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  
  TH2D *hEWKMet2D    = new TH2D("hEWKMet2D", "",  NBINS,0,METMAX,nWBins,vWPtBins); //hEWKMet->Sumw2();
  TH2D *hEWKMetp2D   = new TH2D("hEWKMetp2D", "", NBINS,0,METMAX,nWBins,vWPtBins); //hEWKMetp->Sumw2();
  TH2D *hEWKMetm2D   = new TH2D("hEWKMetm2D", "", NBINS,0,METMAX,nWBins,vWPtBins); //hEWKMetm->Sumw2();
  
  TH2D *hEWKMetp_PileupUp2D = new TH2D("hEWKMetp_PileupUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hEWKMetm_PileupUp2D = new TH2D("hEWKMetm_PileupUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  TH2D *hEWKMetp_PileupDown2D = new TH2D("hEWKMetp_PileupDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hEWKMetm_PileupDown2D = new TH2D("hEWKMetm_PileupDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  
  
  TH2D *hAntiDataMet2D   = new TH2D("hAntiDataMet2D","",  NBINS,0,METMAX,nWBins,vWPtBins); //hDataMet->Sumw2();
  TH2D *hAntiDataMetm2D  = new TH2D("hAntiDataMetm2D","", NBINS,0,METMAX,nWBins,vWPtBins); //hDataMetm->Sumw2();  
  TH2D *hAntiDataMetp2D  = new TH2D("hAntiDataMetp2D","", NBINS,0,METMAX,nWBins,vWPtBins); //hDataMetp->Sumw2();
  TH2D *hAntiWmunuMet2D  = new TH2D("hAntiWmunuMet2D","", NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMet->Sumw2();
  TH2D *hAntiWmunuMetp2D = new TH2D("hAntiWmunuMetp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hAntiWmunuMetm2D = new TH2D("hAntiWmunuMetm2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  
  TH2D *hAntiWmunuMetp_RecoilUp2D = new TH2D("hAntiWmunuMetp_RecoilUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hAntiWmunuMetm_RecoilUp2D = new TH2D("hAntiWmunuMetm_RecoilUp2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  TH2D *hAntiWmunuMetp_RecoilDown2D = new TH2D("hAntiWmunuMetp_RecoilDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetp->Sumw2();
  TH2D *hAntiWmunuMetm_RecoilDown2D = new TH2D("hAntiWmunuMetm_RecoilDown2D","",NBINS,0,METMAX,nWBins,vWPtBins); //hWmunuMetm->Sumw2();
  
  TH2D *hAntiEWKMet2D    = new TH2D("hAntiEWKMet2D", "",  NBINS,0,METMAX,nWBins,vWPtBins); //hEWKMet->Sumw2();
  TH2D *hAntiEWKMetp2D   = new TH2D("hAntiEWKMetp2D", "", NBINS,0,METMAX,nWBins,vWPtBins); //hEWKMetp->Sumw2();
  TH2D *hAntiEWKMetm2D   = new TH2D("hAntiEWKMetm2D", "", NBINS,0,METMAX,nWBins,vWPtBins); //hEWKMetm->Sumw2();
  
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
  Float_t genVPt, genVPhi;
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
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("puppiMet",      &met);       // MET
    intree->SetBranchAddress("puppiMetPhi",   &metPhi);    // phi(MET)
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
//     for(UInt_t ientry=0; ientry<10000; ientry++) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) std::cout << "On Entry.... " << ientry << std::endl;

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
//       std::cout << "made it!" << std::endl;
      if(typev[ifile]==eData) {
          
        // Apply the Rochester Corrections to data
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;
        rmcor->momcor_data(mu1,q,0,qter1);
        if(mu1.Pt()        < PT_CUT)  continue;
        hDataMet->Fill(met);
        hDataMet2D->Fill(met,fabs(lep->Eta()));
        if(q>0) { 
          hDataMetp2D->Fill(met,fabs(lep->Eta())); 
          hDataMetp->Fill(met); 
        } else { 
          hDataMetm2D->Fill(met,fabs(lep->Eta())); 
          hDataMetm->Fill(met); 
        }
      } else if(typev[ifile]==eAntiData) {
        // Apply the Rochester Corrections to data (anti-isolation)
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;
        rmcor->momcor_data(mu1,q,0,qter1);
        if(mu1.Pt()        < PT_CUT)  continue;
        
        hAntiDataMet->Fill(met);
        if(q>0) { 
          hAntiDataMetp2D->Fill(met,fabs(lep->Eta())); 
          hAntiDataMetp->Fill(met); 
        } else { 
          hAntiDataMetm2D->Fill(met,fabs(lep->Eta())); 
          hAntiDataMetm->Fill(met); 
        }
      } else {
        Double_t weight = 1;Double_t weightUp = 1;Double_t weightDown = 1;
        Double_t weight2 =1;
        weight2*=scale1fb*lumi2*corr;
        weight *= scale1fb*lumi*corr;
        weightUp *= scale1fbUp*lumi*corr;
        weightDown *= scale1fbDown*lumi*corr;
        
        // Do some Rochester corrections for MC
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        float qter1=1.0;

        rmcor->momcor_mc(mu1,q,0,qter1);
        Double_t lepPt = mu1.Pt();
        
        if(typev[ifile]==eWmunu) {
          Double_t corrMet=met, corrMetPhi=metPhi;
          
          if(lepPt        > PT_CUT) {
            double bin = 0;
            for(int i = 1; i <= hh_diff->GetNbinsX();++i){
              if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            }
            double w2 = hh_diff->GetBinContent(bin);
            
              corrMet=met, corrMetPhi=metPhi;
              hWmunuMet->Fill(corrMet,weight);
              hWmunuMet_PileupUp->Fill(corrMet,weightUp);
              hWmunuMet_PileupDown->Fill(corrMet,weightDown);
              if(q>0) {
                recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
                hWmunuMetp->Fill(corrMet,weight); 
                hWmunuMetp2D->Fill(corrMet,fabs(lep->Eta()),weight); 
                corrMet=met, corrMetPhi=metPhi;
                recoilCorrUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
                hWmunuMetp_PileupUp->Fill(corrMet,weightUp);
                hWmunuMetp_PileupUp2D->Fill(corrMet,fabs(lep->Eta()),weightUp);
                corrMet=met, corrMetPhi=metPhi;
                recoilCorrDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
                hWmunuMetp_PileupDown->Fill(corrMet,weightDown);
                hWmunuMetp_PileupDown2D->Fill(corrMet,fabs(lep->Eta()),weightDown);
                corrMet=met, corrMetPhi=metPhi;
              } else { 
                recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
                hWmunuMetm->Fill(corrMet,weight); 
                hWmunuMetm2D->Fill(corrMet,fabs(lep->Eta()),weight); 
                corrMet=met, corrMetPhi=metPhi;
                recoilCorrmUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
                hWmunuMetm_PileupUp->Fill(corrMet,weightUp);
                hWmunuMetm_PileupUp2D->Fill(corrMet,fabs(lep->Eta()),weightUp);
                corrMet=met, corrMetPhi=metPhi;
                recoilCorrmDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
                hWmunuMetm_PileupDown->Fill(corrMet,weightDown);
                hWmunuMetm_PileupDown2D->Fill(corrMet,fabs(lep->Eta()),weightDown);
                corrMet=met, corrMetPhi=metPhi;
              }
              corrMet=met, corrMetPhi=metPhi;
            }
//             //recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),1,q);
//             hWmunuMet_RecoilUp->Fill(corrMet,weight);
//             if(q>0) {
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetp_RecoilUp->Fill(corrMet,weight); 
//               corrMet=met, corrMetPhi=metPhi;
//             } else { 
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetm_RecoilUp->Fill(corrMet,weight);
//               corrMet=met, corrMetPhi=metPhi;
//             }
//             hWmunuMet_RecoilDown->Fill(corrMet,weight);
//             if(q>0) {
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetp_RecoilDown->Fill(corrMet,weight);
//               corrMet=met, corrMetPhi=metPhi;
//             } else {
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetm_RecoilDown->Fill(corrMet,weight);
//               corrMet=met, corrMetPhi=metPhi;
//             }
          }
//           Double_t lepPtup = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),1),getEleResCorr(lep->Eta(),1)));  // (!) uncomment to apply scale/res corrections to MC
//           if(lepPtup        > PT_CUT) {
//             corrMet=met, corrMetPhi=metPhi;
//             hWmunuMet_ScaleUp->Fill(corrMet,weight);
//             if(q>0){
//             pU1 = 0; pU2 = 0; 
//             hWmunuMetp_ScaleUp->Fill(corrMet,weight); 
//             corrMet=met, corrMetPhi=metPhi;
//             } else {
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetm_ScaleUp->Fill(corrMet,weight);
//               corrMet=met, corrMetPhi=metPhi;
//             }
//           }
//           Double_t lepPtdown = (gRandom->Gaus((lep->Pt())*getEleScaleCorr(lep->Eta(),-1),getEleResCorr(lep->Eta(),-1)));  // (!) uncomment to apply scale/res corrections to MC
//           if(lepPtdown        > PT_CUT) {
//             corrMet=met, corrMetPhi=metPhi;
//             hWmunuMet_ScaleDown->Fill(corrMet,weight);
//             if(q>0) {
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetp_ScaleDown->Fill(corrMet,weight);
//               corrMet=met, corrMetPhi=metPhi;
//             } else { 
//               pU1 = 0; pU2 = 0; 
//               hWmunuMetm_ScaleDown->Fill(corrMet,weight); 
//             }
//           }
//         }
        if(typev[ifile]==eAntiWmunu) {
          if(lep->Pt()        < PT_CUT)  continue;
          Double_t corrMet=met, corrMetPhi=metPhi;
          hAntiWmunuMet->Fill(corrMet,weight2);
          hAntiWmunuMet2D->Fill(corrMet,fabs(lep->Eta()),weight2);
          if(q>0) {               
            pU1 = 0; pU2 = 0; 
            recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            hAntiWmunuMetp2D->Fill(corrMet,fabs(lep->Eta()),weight2);
            hAntiWmunuMetp->Fill(corrMet,weight2); 
            corrMet = met; corrMetPhi = metPhi;
          } else { 
            pU1 = 0; pU2 = 0; 
            recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0);
            hAntiWmunuMetm2D->Fill(corrMet,fabs(lep->Eta()),weight2);
            hAntiWmunuMetm->Fill(corrMet,weight2);
            corrMet = met; corrMetPhi = metPhi; 
          }
        }
        if(typev[ifile]==eEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
          hEWKMet->Fill(met,weight);
          hEWKMet_PileupUp->Fill(met,weight);
          hEWKMet_PileupDown->Fill(met,weight);
          if(q>0) {
            hEWKMetp->Fill(met,weight); 
            hEWKMetp_PileupUp->Fill(met,weightUp);
            hEWKMetp_PileupDown->Fill(met,weightDown);
            hEWKMetp2D->Fill(met,fabs(lep->Eta()),weight); 
            hEWKMetp_PileupDown2D->Fill(met,fabs(lep->Eta()),weightDown);
          }
          else { 
            hEWKMetm->Fill(met,weight);
            hEWKMetm_PileupUp->Fill(met,weightUp);
            hEWKMetm_PileupDown->Fill(met,weightDown);
            hEWKMetm2D->Fill(met,fabs(lep->Eta()),weight);
            hEWKMetm_PileupDown2D->Fill(met,fabs(lep->Eta()),weightDown);
          }
        }
        if(typev[ifile]==eAntiEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
          hAntiEWKMet->Fill(met,weight2);
          hAntiEWKMet2D->Fill(met,fabs(lep->Eta()),weight2);
          if(q>0) { 
            hAntiEWKMetp2D->Fill(met,fabs(lep->Eta()),weight); 
            hAntiEWKMetp->Fill(met,weight2); 
          } else { 
            hAntiEWKMetm2D->Fill(met,fabs(lep->Eta()),weight); 
            hAntiEWKMetm->Fill(met,weight2); 
          }
        }
      }
    }
    std::cout << "blah blabh" << std::endl;
  }
  delete infile;
  infile=0, intree=0;   
 
//   // Here rescale the up/down whatever
  // Calculate the shapes for W+
//   TH1D *corrP = (TH1D*) hWmunuMetp->Clone("up");
//   *corrP = (*hh_diffp)*(*hWmunuMetp);
//   hWmunuMetp_RecoilUp->Add(corrP,hWmunuMetp,-1);
//   hWmunuMetp_RecoilDown->Add(corrP,hWmunuMetp,1);
//   hWmunuMetp_ScaleUp->Add(corrP,hWmunuMetp,-1);
//   hWmunuMetp_ScaleDown->Add(corrP,hWmunuMetp,1);
//   // Calculate the shapes for W-
//   TH1D *corrM = (TH1D*) hWmunuMetm->Clone("up");
//   *corrM = (*hh_diffm)*(*hWmunuMetm);
//   hWmunuMetm_RecoilUp->Add(corrM,hWmunuMetm,-1);
//   hWmunuMetm_RecoilDown->Add(corrM,hWmunuMetm,1);
//   hWmunuMetm_ScaleUp->Add(corrM,hWmunuMetm,-1);
//   hWmunuMetm_ScaleDown->Add(corrM,hWmunuMetm,1);
  
  
  
  for(int iBin = 1; iBin < nWBins+1; ++iBin){
    
    // get the up/down shape creating histograms
    TH1D *hh_diffm = new TH1D("hh_diffm","hh_diffm",75,0,150);
    TH1D *hh_diffp = new TH1D("hh_diffp","hh_diffp",75,0,150);
    
    char hname[100];
    sprintf(hname,"hh_diff_%i",iBin);
    hh_diffm = (TH1D*)_rdWmm->Get(hname);
    hh_diffp = (TH1D*)_rdWmp->Get(hname);
    
    std::cout << "wut" << std::endl; 
    TH1D *hWmunuMet_px = (TH1D*) hWmunuMet2D->ProjectionX("hWmunuMet_px",iBin,iBin);
    TH1D *hWmunuMetp_px = (TH1D*) hWmunuMetp2D->ProjectionX("hWmunuMetp_px",iBin,iBin);
    TH1D *hWmunuMetm_px = (TH1D*) hWmunuMetm2D->ProjectionX("hWmunuMetm_px",iBin,iBin);
    TH1D *hDataMet_px = (TH1D*) hDataMet2D->ProjectionX("hDataMet_px",iBin,iBin);
    TH1D *hDataMetp_px = (TH1D*) hDataMetp2D->ProjectionX("hDataMetp_px",iBin,iBin);
    TH1D *hDataMetm_px = (TH1D*) hDataMetm2D->ProjectionX("hDataMetm_px",iBin,iBin);
    TH1D *hEWKMet_px = (TH1D*) hEWKMet2D->ProjectionX("hEwkMet_px",iBin,iBin);
    TH1D *hEWKMetp_px = (TH1D*) hEWKMetp2D->ProjectionX("hEwkMetp_px",iBin,iBin);
    TH1D *hEWKMetm_px = (TH1D*) hEWKMetm2D->ProjectionX("hEwkMetm_px",iBin,iBin);
    
//     TH1D *hWmunuMet_PileupUp_px = (TH1D*) hWmunuMet_PileupUp2D->ProjectionX("hWmunuMet_PileupUp_px",iBin,iBin);
    TH1D *hWmunuMetp_PileupUp_px = (TH1D*) hWmunuMetp_PileupUp2D->ProjectionX("hWmunuMetp_PileupUp_px",iBin,iBin);
    TH1D *hWmunuMetm_PileupUp_px = (TH1D*) hWmunuMetm_PileupUp2D->ProjectionX("hWmunuMetm_PileupUp_px",iBin,iBin);
    
//     TH1D *hWmunuMet_PileupDown_px = (TH1D*) hWmunuMet_PileupDown2D->ProjectionX("hWmunuMet_PileupDown_px",iBin,iBin);
    TH1D *hWmunuMetp_PileupDown_px = (TH1D*) hWmunuMetp_PileupDown2D->ProjectionX("hWmunuMetp_PileupDown_px",iBin,iBin);
    TH1D *hWmunuMetm_PileupDown_px = (TH1D*) hWmunuMetm_PileupDown2D->ProjectionX("hWmunuMetm_PileupDown_px",iBin,iBin);
    
//     TH1D *hEWKMet_PileupUp_px = (TH1D*) hEWKMet_PileupUp2D->ProjectionX("hEwkMet_PileupUp_px",iBin,iBin);
    TH1D *hEWKMetp_PileupUp_px = (TH1D*) hEWKMetp_PileupUp2D->ProjectionX("hEwkMetp_PileupUp_px",iBin,iBin);
    TH1D *hEWKMetm_PileupUp_px = (TH1D*) hEWKMetm_PileupUp2D->ProjectionX("hEwkMetm_PileupUp_px",iBin,iBin);
    
//     TH1D *hEWKMet_PileupDown_px = (TH1D*) hEWKMet_PileupDown2D->ProjectionX("hEwkMet_PileupDown_px",iBin,iBin);
    TH1D *hEWKMetp_PileupDown_px = (TH1D*) hEWKMetp_PileupDown2D->ProjectionX("hEwkMetp_PileupDown_px",iBin,iBin);
    TH1D *hEWKMetm_PileupDown_px = (TH1D*) hEWKMetm_PileupDown2D->ProjectionX("hEwkMetm_PileupDown_px",iBin,iBin);
    
    std::cout << "hello" << std::endl;
    TH1D *hAntiWmunuMet_px  = (TH1D*) hAntiWmunuMet2D->ProjectionX("hWmunuMet_px",iBin,iBin);
    TH1D *hAntiWmunuMetp_px = (TH1D*) hAntiWmunuMetp2D->ProjectionX("hWmunuMetp_px",iBin,iBin);
    TH1D *hAntiWmunuMetm_px = (TH1D*) hAntiWmunuMetm2D->ProjectionX("hWmunuMetm_px",iBin,iBin);
    TH1D *hAntiDataMet_px  = (TH1D*) hAntiDataMet2D->ProjectionX("hDataMet_px",iBin,iBin);
    TH1D *hAntiDataMetp_px = (TH1D*) hAntiDataMetp2D->ProjectionX("hDataMetp_px",iBin,iBin);
    TH1D *hAntiDataMetm_px = (TH1D*) hAntiDataMetm2D->ProjectionX("hDataMetm_px",iBin,iBin);
    TH1D *hAntiEWKMet_px  = (TH1D*) hAntiEWKMet2D->ProjectionX("hEwkMet_px",iBin,iBin);
    TH1D *hAntiEWKMetp_px = (TH1D*) hAntiEWKMetp2D->ProjectionX("hEwkMetp_px",iBin,iBin);
    TH1D *hAntiEWKMetm_px = (TH1D*) hAntiEWKMetm2D->ProjectionX("hEwkMetm_px",iBin,iBin);
    
    //
    // Declare fit parameters for signal and background yields
    // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
    std::cout << "blah" << std::endl;
    //  Calculate the shapes for W+
    TH1D *corrP = (TH1D*) hWmunuMetp_px->Clone("up");
    *corrP = (*hh_diffp)*(*hWmunuMetp_px);
    hWmunuMetp_RecoilUp->Add(corrP,hWmunuMetp_px,-1);
    hWmunuMetp_RecoilDown->Add(corrP,hWmunuMetp_px,1);
    //hWmunuMetp_ScaleUp->Add(corrP,hWmunuMetp_px,-1);
    //hWmunuMetp_ScaleDown->Add(corrP,hWmunuMetp_px,1);
    // Calculate the shapes for W-
    TH1D *corrM = (TH1D*) hWmunuMetm_px->Clone("up");
    *corrM = (*hh_diffm)*(*hWmunuMetm_px);
    hWmunuMetm_RecoilUp->Add(corrM,hWmunuMetm_px,-1);
    hWmunuMetm_RecoilDown->Add(corrM,hWmunuMetm_px,1);
    //hWmunuMetm_ScaleUp->Add(corrM,hWmunuMetm_px,-1);
    //hWmunuMetm_ScaleDown->Add(corrM,hWmunuMetm_px,1);
    std::cout << "after recoil " << std::endl;    
      //
    RooRealVar nSig("nSig","nSig",0.7*(hDataMet_px->Integral()),0,hDataMet_px->Integral());
    RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet_px->Integral()),0,hDataMet_px->Integral());
    RooRealVar cewk("cewk","cewk",0.1,0,5) ;
    cewk.setVal(hEWKMet_px->Integral()/hWmunuMet_px->Integral());
    cewk.setConstant(kTRUE);
    RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
    RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWmunuMet->Integral()*0.9,0,hAntiDataMet->Integral());
    RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
    RooRealVar dewk("dewk","dewk",0.1,0,5) ;
    dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
    dewk.setConstant(kTRUE);
  //   nAntiSig.setConstant(kTRUE);
    RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
    
    RooRealVar nSigp("nSigp","nSigp",hDataMetp_px->Integral()*0.8,0,hDataMetp_px->Integral());
    RooRealVar nQCDp("nQCDp","nQCDp",hDataMetp_px->Integral()*0.1,0,hDataMetp_px->Integral()*0.5);
    RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
    cewkp.setVal(hEWKMetp_px->Integral()/hWmunuMetp_px->Integral());
    cewkp.setConstant(kTRUE);
    RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
    RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral()*0.9,0,hAntiDataMetp->Integral());
    RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
    RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
    dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
    dewkp.setConstant(kTRUE);
  //   nAntiSigp.setConstant(kTRUE);
    RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
    
    RooRealVar nSigm("nSigm","nSigm",hDataMetm_px->Integral()*0.9,0,hDataMetm_px->Integral());
    RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm_px->Integral()*0.2,0,hDataMetm_px->Integral()*0.5);
    RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
    cewkm.setVal(hEWKMetm_px->Integral()/hWmunuMetm_px->Integral());
    cewkm.setConstant(kTRUE);
    RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
    RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWmunuMetm->Integral()*0.9,0,hAntiDataMetm->Integral());
    RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
    RooRealVar dewkm("dewkm","dewkm",0.1,0,5);
    dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
    dewkm.setConstant(kTRUE);
  //   nAntiSigm.setConstant(kTRUE);
    RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));

   std::cout << "helllllo" << std::endl;
    
    
    //
    // Construct PDFs for fitting
    //
    RooRealVar pfmet("pfmet","pfmet",0,METMAX);
    pfmet.setBins(NBINS);
    
    // Signal PDFs
    RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(pfmet),hWmunuMet_px);  RooHistPdf pdfWm ("wm", "wm", pfmet,wmunuMet, 1);
    RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(pfmet),hWmunuMetp_px); RooHistPdf pdfWmp("wmp","wmp",pfmet,wmunuMetp,1);
    RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(pfmet),hWmunuMetm_px); RooHistPdf pdfWmm("wmm","wmm",pfmet,wmunuMetm,1); 
    RooDataHist wmunuMet_RecoilUp("wmunuMET_RecoilUp", "wmunuMET_RecoilUp", RooArgSet(pfmet),hWmunuMet_RecoilUp);  RooHistPdf pdfWm_RecoilUp("wm_RecoilUp", "wm_RecoilUp", pfmet,wmunuMet_RecoilUp, 1);
    RooDataHist wmunuMetp_RecoilUp("wmunuMETp_RecoilUp","wmunuMETp_RecoilUp",RooArgSet(pfmet),hWmunuMetp_RecoilUp); RooHistPdf pdfWmp_RecoilUp("wmp_RecoilUp","wmp_RecoilUp",pfmet,wmunuMetp_RecoilUp,1);
    RooDataHist wmunuMetm_RecoilUp("wmunuMETm_RecoilUp","wmunuMETm_RecoilUp",RooArgSet(pfmet),hWmunuMetm_RecoilUp); RooHistPdf pdfWmm_RecoilUp("wmm_RecoilUp","wmm_RecoilUp",pfmet,wmunuMetm_RecoilUp,1); 
    RooDataHist wmunuMet_RecoilDown("wmunuMET_RecoilDown", "wmunuMET_RecoilDown", RooArgSet(pfmet),hWmunuMet_RecoilDown);  RooHistPdf pdfWm_RecoilDown("wm_RecoilDown", "wm_RecoilDown", pfmet,wmunuMet_RecoilDown, 1);
    RooDataHist wmunuMetp_RecoilDown("wmunuMETp_RecoilDown","wmunuMETp_RecoilDown",RooArgSet(pfmet),hWmunuMetp_RecoilDown); RooHistPdf pdfWmp_RecoilDown("wmp_RecoilDown","wmp_RecoilDown",pfmet,wmunuMetp_RecoilDown,1);
    RooDataHist wmunuMetm_RecoilDown("wmunuMETm_RecoilDown","wmunuMETm_RecoilDown",RooArgSet(pfmet),hWmunuMetm_RecoilDown); RooHistPdf pdfWmm_RecoilDown("wmm_RecoilDown","wmm_RecoilDown",pfmet,wmunuMetm_RecoilDown,1); 
    
//     RooDataHist wmunuMet_PileupUp ("wmunuMET_PileupUp", "wmunuMET_PileupUp", RooArgSet(pfmet),hWmunuMet_PileupUp_px);  RooHistPdf  pdfWm_PileupUp("wm_PileupUp", "wm_PileupUp", pfmet,wmunuMet_PileupUp, 1);
    RooDataHist wmunuMetp_PileupUp("wmunuMETp_PileupUp","wmunuMETp_PileupUp",RooArgSet(pfmet),hWmunuMetp_PileupUp_px); RooHistPdf pdfWmp_PileupUp("wmp_PileupUp","wmp_PileupUp",pfmet,wmunuMetp_PileupUp,1);
    RooDataHist wmunuMetm_PileupUp("wmunuMETm_PileupUp","wmunuMETm_PileupUp",RooArgSet(pfmet),hWmunuMetm_PileupUp_px); RooHistPdf pdfWmm_PileupUp("wmm_PileupUp","wmm_PileupUp",pfmet,wmunuMetm_PileupUp,1); 
    
//     RooDataHist wmunuMet_PileupDown ("wmunuMET_PileupDown", "wmunuMET_PileupDown", RooArgSet(pfmet),hWmunuMet_PileupDown_px);  RooHistPdf pdfWm_PileupDown ("wm_PileupDown", "wm", pfmet,wmunuMet_PileupDown, 1);
    RooDataHist wmunuMetp_PileupDown("wmunuMETp_PileupDown","wmunuMETp_PileupDown",RooArgSet(pfmet),hWmunuMetp_PileupDown_px); RooHistPdf pdfWmp_PileupDown("wmp_PileupDown","wmp_PileupDown",pfmet,wmunuMetp_PileupDown,1);
    RooDataHist wmunuMetm_PileupDown("wmunuMETm_PileupDown","wmunuMETm_PileupDown",RooArgSet(pfmet),hWmunuMetm_PileupDown_px); RooHistPdf pdfWmm_PileupDown("wmm_PileupDown","wmm_PileupDown",pfmet,wmunuMetm_PileupDown,1); 
    
    RooDataHist wmunuMet_ScaleUp("wmunuMET_ScaleUp", "wmunuMET_ScaleUp", RooArgSet(pfmet),hWmunuMet_ScaleUp);  RooHistPdf pdfWm_ScaleUp("wm_ScaleUp", "wm_ScaleUp", pfmet,wmunuMet_ScaleUp, 1);
    RooDataHist wmunuMetp_ScaleUp("wmunuMETp_ScaleUp","wmunuMETp_ScaleUp",RooArgSet(pfmet),hWmunuMetp_ScaleUp); RooHistPdf pdfWmp_ScaleUp("wmp_ScaleUp","wmp_ScaleUp",pfmet,wmunuMetp_ScaleUp,1);
    RooDataHist wmunuMetm_ScaleUp("wmunuMETm_ScaleUp","wmunuMETm_ScaleUp",RooArgSet(pfmet),hWmunuMetm_ScaleUp); RooHistPdf pdfWmm_ScaleUp("wmm_ScaleUp","wmm_ScaleUp",pfmet,wmunuMetm_ScaleUp,1); 
    RooDataHist wmunuMet_ScaleDown("wmunuMET_ScaleDown", "wmunuMET_ScaleDown", RooArgSet(pfmet),hWmunuMet_ScaleDown);  RooHistPdf pdfWm_ScaleDown("wm_ScaleDown", "wm_ScaleDown", pfmet,wmunuMet_ScaleDown, 1);
    RooDataHist wmunuMetp_ScaleDown("wmunuMETp_ScaleDown","wmunuMETp_ScaleDown",RooArgSet(pfmet),hWmunuMetp_ScaleDown); RooHistPdf pdfWmp_ScaleDown("wmp_ScaleDown","wmp_ScaleDown",pfmet,wmunuMetp_ScaleDown,1);
    RooDataHist wmunuMetm_ScaleDown("wmunuMETm_ScaleDown","wmunuMETm_ScaleDown",RooArgSet(pfmet),hWmunuMetm_ScaleDown); RooHistPdf pdfWmm_ScaleDown("wmm_ScaleDown","wmm_ScaleDown",pfmet,wmunuMetm_ScaleDown,1); 
    
    // EWK+top PDFs
    RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet_px);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
    RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp_px); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
    RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm_px); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
    
//     RooDataHist ewkMet_PileupUp ("ewkMET_PileupUp", "ewkMET_PileupUp", RooArgSet(pfmet),hEWKMet_PileupUp_px);  RooHistPdf pdfEWK_PileupUp ("ewk_PileupUp", "ewk_PileupUp", pfmet,ewkMet_PileupUp, 1);
    RooDataHist ewkMetp_PileupUp("ewkMETp_PileupUp","ewkMETp_PileupUp",RooArgSet(pfmet),hEWKMetp_PileupUp_px); RooHistPdf pdfEWKp_PileupUp("ewkp_PileupUp","ewkp_PileupUp",pfmet,ewkMetp_PileupUp,1); 
    RooDataHist ewkMetm_PileupUp("ewkMETm_PileupUp","ewkMETm_PileupUp",RooArgSet(pfmet),hEWKMetm_PileupUp_px); RooHistPdf pdfEWKm_PileupUp("ewkm_PileupUp","ewkm_PileupUp",pfmet,ewkMetm_PileupUp,1); 
    
//     RooDataHist ewkMet_PileupDown ("ewkMET_PileupDown", "ewkMET_PileupDown", RooArgSet(pfmet),hEWKMet_PileupDown_px);  RooHistPdf pdfEWK_PileupDown ("ewk_PileupDown", "ewk_PileupDown", pfmet,ewkMet_PileupDown, 1);
    RooDataHist ewkMetp_PileupDown("ewkMETp_PileupDown","ewkMETp_PileupDown",RooArgSet(pfmet),hEWKMetp_PileupDown_px); RooHistPdf pdfEWKp_PileupDown("ewkp_PileupDown","ewkp_PileupDown",pfmet,ewkMetp_PileupDown,1); 
    RooDataHist ewkMetm_PileupDown("ewkMETm_PileupDown","ewkMETm_PileupDown",RooArgSet(pfmet),hEWKMetm_PileupDown_px); RooHistPdf pdfEWKm_PileupDown("ewkm_PileupDown","ewkm_PileupDown",pfmet,ewkMetm_PileupDown,1); 
    
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
  //   CPepeModel1 aqcd("aqcd",pfmet, qcd.a1);
  //   CPepeModel1 aqcdp("aqcdp",pfmet, qcdp.a1);
  //   CPepeModel1 aqcdm("aqcdm",pfmet, qcdm.a1);

    CPepeModel2 aqcd("aqcd",pfmet);
    CPepeModel2 aqcdp("aqcdp",pfmet);
    CPepeModel2 aqcdm("aqcdm",pfmet);
    
    // Anti-selection PDFs
    RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
    RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
    RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
    
    // PDF for simultaneous fit
    RooCategory rooCat("rooCat","rooCat");
    rooCat.defineType("Selectp");
    rooCat.defineType("Selectm");
    
    RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
    pdfTotal.addPdf(pdfMet, "Select");
    //pdfTotal.addPdf(apdfMet,"Anti");
    
    RooSimultaneous pdfTotalp("pdfTotalp","pdfTotalp",rooCat);
    pdfTotalp.addPdf(pdfMetp, "Selectp");
  //   pdfTotalp.addPdf(pdfMetm,"Selectm");
    //pdfTotalp.addPdf(apdfMetp,"Anti");
    
    RooSimultaneous pdfTotalm("pdfTotalm","pdfTotalm",rooCat);
    pdfTotalm.addPdf(pdfMetm, "Select");
    //pdfTotalm.addPdf(apdfMetm,"Anti");

    
    //
    // Perform fits
    //

    RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet), hDataMet_px);
    RooDataHist dataMetp("dataMetp", "dataMetp", RooArgSet(pfmet), hDataMetp_px);
    RooDataHist dataMetm("dataMetm", "dataMetm", RooArgSet(pfmet), hDataMetm_px);
    
    RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
    RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
    RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);

    cout << "Starting values for Wmunu yields: " << endl;
    cout << "   sig: " << hWmunuMet_px->Integral() << endl;
    cout << "   EWK: " << hEWKMet_px->Integral() << endl;
    cout << "   qcd: " << hDataMet_px->Integral()-hWmunuMet_px->Integral()-hEWKMet_px->Integral() << endl;

    cout << "Starting values for Wmunu_p yields: " << endl;
    cout << "   sig: " << hWmunuMetp_px->Integral() << endl;
    cout << "   EWK: " << hEWKMetp_px->Integral() << endl;
    cout << "   qcd: " << hDataMetp_px->Integral()-hWmunuMetp_px->Integral()-hEWKMetp_px->Integral() << endl;

    cout << "Starting values for Wmunu_m yields: " << endl;
    cout << "   sig: " << hWmunuMetm_px->Integral() << endl;
    cout << "   EWK: " << hEWKMetm_px->Integral() << endl;
    cout << "   qcd: " << hDataMetm_px->Integral()-hWmunuMetm_px->Integral()-hEWKMetm_px->Integral() << endl;

  //   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  //   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
    RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",hDataMetp_px->Integral()*0.2,0,hDataMetp_px->Integral()*0.5);
    RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",hDataMetm_px->Integral()*0.2,0,hDataMetm_px->Integral()*0.5);
    
    RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",hDataMetp_px->Integral()*0.2,0,hDataMetp_px->Integral());
    RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",hDataMetm_px->Integral()*0.2,0,hDataMetm_px->Integral());

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
//     combine_workspace.import(pdfWm_PileupUp);
    combine_workspace.import(pdfWmp_PileupUp);
    combine_workspace.import(pdfWmm_PileupUp);
//     combine_workspace.import(pdfWm_PileupDown);
    combine_workspace.import(pdfWmp_PileupDown);
    combine_workspace.import(pdfWmm_PileupDown);
    combine_workspace.import(pdfWm_ScaleUp);
    combine_workspace.import(pdfWmp_ScaleUp);
    combine_workspace.import(pdfWmm_ScaleUp);
    combine_workspace.import(pdfWm_ScaleDown);
    combine_workspace.import(pdfWmp_ScaleDown);
    combine_workspace.import(pdfWmm_ScaleDown);
    combine_workspace.import(pdfEWK);
    combine_workspace.import(pdfEWKp);
    combine_workspace.import(pdfEWKm);
//     combine_workspace.import(pdfEWK_PileupUp);
    combine_workspace.import(pdfEWKp_PileupUp);
    combine_workspace.import(pdfEWKm_PileupUp);
//     combine_workspace.import(pdfEWK_PileupDown);
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
    

  //   char outfile[100];
    stringstream outfile;
    outfile << outputDir << "/Wmunu_pdfTemplates_" << iBin << ".root";
  //   sprintf(outfile, "%s/Wmunu_pdfTemplates_%i.root",outputDir, iBin);
    combine_workspace.writeToFile(outfile.str().c_str());

    
    RooDataHist dataTotal("dataTotal","dataTotal", RooArgList(pfmet), Index(rooCat),
              Import("Select", dataMet),
              Import("Anti",   antiMet));
    RooFitResult *fitRes = 0;//pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));//dataTotal.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));
    
  
    
    RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
              Import("Select", dataMetp),
              Import("Anti",   antiMetp));
              
    RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
    RooFitResult *fitResAntip = apdfMetp.fitTo(antiMetp,Extended(),Minos(kTRUE),Save(kTRUE));
  // RooFitResult *fitResp = pdfMetp.fitTo(dataTotalp,Extended(),Minos(kTRUE),Save(kTRUE));
    

  //   RooDataHist dataMetp2("dataMetp2", "dataMetp2", RooArgSet(pfmet), hDataMetp);
    RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
              Import("Select", dataMetm),
              Import("Anti", antiMetm));
    RooFitResult *fitResm = pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
  //   RooFitResult *fitResm = pdfMetm.fitTo(dataTotalm,Extended(),Minos(kTRUE),Save(kTRUE));
    RooFitResult *fitResAntim = apdfMetm.fitTo(antiMetm,Extended(),Minos(kTRUE),Save(kTRUE));
  //   RooFitResult *fitResm = pdfTotalp.fitTo(dataTotalm,Extended(),Minos(kTRUE),Save(kTRUE));
    
    //
    // Use histogram version of fitted PDFs to make ratio plots
    // (Will also use PDF histograms later for Chi^2 and KS tests)
    //
    TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
    hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
    TH1D *hMetDiff = makeDiffHist(hDataMet_px,hPdfMet,"hMetDiff");
    hMetDiff->SetMarkerStyle(kFullCircle);
    hMetDiff->SetMarkerSize(0.9);
    
    TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
    hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
    TH1D *hMetpDiff = makeDiffHist(hDataMetp_px,hPdfMetp,"hMetpDiff");
    hMetpDiff->SetMarkerStyle(kFullCircle);
    hMetpDiff->SetMarkerSize(0.9);
      
    TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
    hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
    TH1D *hMetmDiff = makeDiffHist(hDataMetm_px,hPdfMetm,"hMetmDiff");
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
    else         sprintf(lumitext,"%.0f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
    
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
    RooPlot *wmframe = pfmet.frame(Bins(NBINS)); 
    wmframe->GetYaxis()->SetNdivisions(505);
    dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMet.plotOn(wmframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMet.plotOn(wmframe,LineColor(linecolorW));
    pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
    pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
    pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
    pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
    pdfMet.plotOn(wmframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
    dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
    
    
    stringstream plotName;
    plotName.str(""); plotName << "fitmet_" << iBin;
    sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
    CPlot plotMet(plotName.str().c_str(),wmframe,"","",ylabel);
    plotMet.SetLegend(0.68,0.57,0.93,0.77);
    plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
    plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
    plotMet.Draw(c,kFALSE,format,1);

    CPlot plotMetDiff(plotName.str().c_str(),"","#slash{E}_{T} [GeV]","MC/Data");
    //CPlot plotMetDiff("fitmet","","mT [GeV]","#chi");
    plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
    plotMetDiff.SetYRange(-1,1);
    plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotMetDiff.Draw(c,kTRUE,format,2);
    
    plotName.str(""); plotName << "fitmetlog_" << iBin;
    plotMet.SetName(plotName.str().c_str());
    plotMet.SetLogy();
    plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
    plotMet.Draw(c,kTRUE,format,1);
      
    RooPlot *awmframe = pfmet.frame(Bins(NBINS));    
    antiMet.plotOn(awmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    apdfMet.plotOn(awmframe,FillColor(fillcolorW),DrawOption("F"));
    apdfMet.plotOn(awmframe,LineColor(linecolorW));
    apdfMet.plotOn(awmframe,Components(RooArgSet(apdfEWK,*(aqcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
    apdfMet.plotOn(awmframe,Components(RooArgSet(apdfEWK,*(aqcd.model))),LineColor(linecolorEWK));
    apdfMet.plotOn(awmframe,Components(RooArgSet(*(aqcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
    apdfMet.plotOn(awmframe,Components(RooArgSet(*(aqcd.model))),LineColor(linecolorQCD));
    apdfMet.plotOn(awmframe,Components(RooArgSet(apdfWm)),LineColor(linecolorW),LineStyle(2));
    antiMet.plotOn(awmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
    
    plotName.str(""); plotName << "fitantimet_" << iBin;
    sprintf(ylabel,"Events / %.1f GeV",hAntiDataMet->GetBinWidth(1));
    CPlot plotAntiMet(plotName.str().c_str(),awmframe,"","",ylabel);
    plotAntiMet.SetLegend(0.68,0.57,0.93,0.77);
    plotAntiMet.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotAntiMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
    plotAntiMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotAntiMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotAntiMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotAntiMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotAntiMet.SetYRange(0.1,1.1*(hAntiDataMet->GetMaximum())); 
    plotAntiMet.Draw(c,kFALSE,format,1);

    CPlot plotAntiMetDiff(plotName.str().c_str(),"","#slash{E}_{T} [GeV]","#chi");
    plotAntiMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
    plotAntiMetDiff.SetYRange(-0.2,0.2);
    plotAntiMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotAntiMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotAntiMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotAntiMetDiff.Draw(c,kTRUE,format,2);
    
    plotName.str(""); plotName << "fitantimetlog_" << iBin;
    plotAntiMet.SetName(plotName.str().c_str());
    plotAntiMet.SetLogy();
    plotAntiMet.SetYRange(1e-3*(hAntiDataMet->GetMaximum()),10*(hAntiDataMet->GetMaximum()));
    plotAntiMet.Draw(c,kTRUE,format,1);
      
    //
    // W+ MET plot
    //
    RooPlot *wmpframe = pfmet.frame(Bins(NBINS));
    wmpframe->GetYaxis()->SetNdivisions(505);
    dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMetp.plotOn(wmpframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMetp.plotOn(wmpframe,LineColor(linecolorW));
    pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
    pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
    pdfMetp.plotOn(wmpframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
    pdfMetp.plotOn(wmpframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
    pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
    dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
    
  
    plotName.str(""); plotName << "fitmetp_" << iBin;
    sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
    CPlot plotMetp(plotName.str().c_str(),wmpframe,"","",ylabel);
    plotMetp.SetLegend(0.68,0.57,0.93,0.77);
    plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
    plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMetp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMetp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  //  plotMetp.SetYRange(0.1,1.1*(hDataMetp->GetMaximum()));
    //plotMetp.SetYRange(0.1,4100);
    plotMetp.Draw(c,kFALSE,format,1);

    CPlot plotMetpDiff(plotName.str().c_str(),"","#slash{E}_{T} [GeV]","#chi");
    //CPlot plotMetpDiff("fitmetp","","mT [GeV]","#chi");
    plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
    plotMetpDiff.SetYRange(-0.2,0.2);
    plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetpDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotMetpDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotMetpDiff.Draw(c,kTRUE,format,2);
    
    plotName.str(""); plotName << "fitmetplog_" << iBin;
    plotMetp.SetName(plotName.str().c_str());
    plotMetp.SetLogy();
    plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
    plotMetp.Draw(c,kTRUE,format,1);

    RooPlot *awmpframe = pfmet.frame(Bins(NBINS));    
    antiMetp.plotOn(awmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    apdfMetp.plotOn(awmpframe,FillColor(fillcolorW),DrawOption("F"));
    apdfMetp.plotOn(awmpframe,LineColor(linecolorW));
    apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
    apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),LineColor(linecolorEWK));
    apdfMetp.plotOn(awmpframe,Components(RooArgSet(*(aqcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
    apdfMetp.plotOn(awmpframe,Components(RooArgSet(*(aqcdp.model))),LineColor(linecolorQCD));
    apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfWmp)),LineColor(linecolorW),LineStyle(2));
    antiMetp.plotOn(awmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
    
    plotName.str(""); plotName << "fitantimetp_" << iBin;
    sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
    CPlot plotAntiMetp(plotName.str().c_str(),awmpframe,"","",ylabel);
    plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
    plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
    plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotAntiMetp.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotAntiMetp.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  //  plotAntiMetp.SetYRange(0.1,1.1*(hAntiDataMetp->GetMaximum()));
    plotAntiMetp.SetYRange(0.1,1500);
    plotAntiMetp.Draw(c,kFALSE,format,1);

    CPlot plotAntiMetpDiff(plotName.str().c_str(),"","#slash{E}_{T} [GeV]","#chi");
    plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
    plotAntiMetpDiff.SetYRange(-0.2,0.2);
    plotAntiMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotAntiMetpDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotAntiMetpDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotAntiMetpDiff.Draw(c,kTRUE,format,2);
    
    plotName.str(""); plotName << "fitantimetplog_" << iBin;
    plotAntiMetp.SetName(plotName.str().c_str());
    plotAntiMetp.SetLogy();
    plotAntiMetp.SetYRange(1e-3*(hAntiDataMetp->GetMaximum()),10*(hAntiDataMetp->GetMaximum()));
    plotAntiMetp.Draw(c,kTRUE,format,1);
    
    //
    // W- MET plot
    //
    RooPlot *wmmframe = pfmet.frame(Bins(NBINS)); 
    wmmframe->GetYaxis()->SetNdivisions(505);
    dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMetm.plotOn(wmmframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMetm.plotOn(wmmframe,LineColor(linecolorW));
    pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
    pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
    pdfMetm.plotOn(wmmframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
    pdfMetm.plotOn(wmmframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
    pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
    dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    
    plotName.str(""); plotName << "fitmetm_" << iBin;
    sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
    CPlot plotMetm(plotName.str().c_str(),wmmframe,"","",ylabel);
    plotMetm.SetLegend(0.68,0.57,0.93,0.77);
    plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
    plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMetm.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMetm.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  //  plotMetm.SetYRange(0.1,1.1*(hDataMetm->GetMaximum()));
  //plotMetm.SetYRange(0.1,4100);
    plotMetm.Draw(c,kFALSE,format,1);

    CPlot plotMetmDiff(plotName.str().c_str(),"","#slash{E}_{T} [GeV]","#chi");
    //CPlot plotMetmDiff("fitmetm","","mT [GeV]","#chi");
    plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
    plotMetmDiff.SetYRange(-0.2,0.2);
    plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetmDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotMetmDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotMetmDiff.Draw(c,kTRUE,format,2);
    
    plotName.str(""); plotName << "fitmetmlog_" << iBin;
    plotMetm.SetName(plotName.str().c_str());
    plotMetm.SetLogy();
    plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
    plotMetm.Draw(c,kTRUE,format,1);

    RooPlot *awmmframe = pfmet.frame(Bins(NBINS)); 
    antiMetm.plotOn(awmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    apdfMetm.plotOn(awmmframe,FillColor(fillcolorW),DrawOption("F"));
    apdfMetm.plotOn(awmmframe,LineColor(linecolorW));
    apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
    apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),LineColor(linecolorEWK));
    apdfMetm.plotOn(awmmframe,Components(RooArgSet(*(aqcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
    apdfMetm.plotOn(awmmframe,Components(RooArgSet(*(aqcdm.model))),LineColor(linecolorQCD));
    apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfWmm)),LineColor(linecolorW),LineStyle(2));
    antiMetm.plotOn(awmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    
    plotName.str(""); plotName << "fitantimetm_" << iBin;
    sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
    CPlot plotAntiMetm(plotName.str().c_str(),awmmframe,"","",ylabel);
    plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
    plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
    plotAntiMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotAntiMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotAntiMetm.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotAntiMetm.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  //  plotAntiMetm.SetYRange(0.1,1.1*(hAntiDataMetm->GetMaximum()));
    plotAntiMetm.SetYRange(0.1,1500);
    plotAntiMetm.Draw(c,kFALSE,format,1);

    CPlot plotAntiMetmDiff(plotName.str().c_str(),"","#slash{E}_{T} [GeV]","#chi");
    plotAntiMetmDiff.AddHist1D(hAntiMetmDiff,"EX0",ratioColor);
    plotAntiMetmDiff.SetYRange(-0.2,0.2);
    plotAntiMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotAntiMetmDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotAntiMetmDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotAntiMetmDiff.Draw(c,kTRUE,format,2);
    
    plotName.str(""); plotName << "fitantimetmlog_" << iBin;
    plotAntiMetm.SetName(plotName.str().c_str());
    plotAntiMetm.SetLogy();
    plotAntiMetm.SetYRange(1e-3*(hAntiDataMetm->GetMaximum()),10*(hAntiDataMetm->GetMaximum()));
    plotAntiMetm.Draw(c,kTRUE,format,1);

      
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
    sprintf(txtfname,"%s/fitresWm_%i.txt",CPlot::sOutDir.Data(),iBin);
    txtfile.open(txtfname);
    assert(txtfile.is_open());
    
    flags = txtfile.flags();
    txtfile << setprecision(10);
    txtfile << " *** Yields *** " << endl;
    txtfile << "Selected: " << hDataMet->Integral() << endl;
    //txtfile << "  Signal: " << nSig.getVal() << " +/- " << nSig.getPropagatedError(*fitRes) << endl;
    //txtfile << "     QCD: " << nQCD.getVal() << " +/- " << nQCD.getPropagatedError(*fitRes) << endl;
    //txtfile << "   Other: " << nEWK.getVal() << " +/- " << nEWK.getPropagatedError(*fitRes) << endl;
    txtfile << endl;
    txtfile.flags(flags);
    
    //fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
    txtfile << endl;
    //printCorrelations(txtfile, fitRes);
    txtfile << endl;
    printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
    txtfile.close();
    
    chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
    chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
    ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
    ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
    sprintf(txtfname,"%s/fitresWmp_%i.txt",CPlot::sOutDir.Data(),iBin);
    txtfile.open(txtfname);
    assert(txtfile.is_open());
    
    flags = txtfile.flags();
    txtfile << setprecision(10);
    txtfile << " *** Yields *** " << endl;
    txtfile << "Selected: " << hDataMetp->Integral() << endl;
    //txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
    //txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
    //txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
    txtfile << endl; 
    txtfile.flags(flags);
    
    //fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
    txtfile << endl;
    //printCorrelations(txtfile, fitResp);
    txtfile << endl;
    printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
    txtfile.close();

    chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
    chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
    ksprob   = hDataMetm->KolmogorovTest(hPdfMetm);
    ksprobpe = hDataMetm->KolmogorovTest(hPdfMetm,"DX");  
    sprintf(txtfname,"%s/fitresWmm_%i.txt",CPlot::sOutDir.Data(),iBin);
    txtfile.open(txtfname);
    assert(txtfile.is_open());
    
    flags = txtfile.flags();
    txtfile << setprecision(10);
    txtfile << " *** Yields *** " << endl;
    txtfile << "Selected: " << hDataMetm->Integral() << endl;
    txtfile << "  Signal: " << nSigm.getVal() << " +/- " << nSigm.getPropagatedError(*fitResm) << endl;
    txtfile << "     QCD: " << nQCDm.getVal() << " +/- " << nQCDm.getPropagatedError(*fitResm) << endl;
    txtfile << "   Other: " << nEWKm.getVal() << " +/- " << nEWKm.getPropagatedError(*fitResm) << endl;
    txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResm) << endl;
    txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResm) << endl;
    txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResm) << endl;
    
    txtfile << "reversed isolation sample" << endl;
    txtfile << "Selected: " << hDataMetm->Integral() << endl;
    txtfile << "  Signal: " << nAntiSigm.getVal() << " +/- " << nAntiSigm.getPropagatedError(*fitResm) << endl;
    txtfile << "     QCD: " << nAntiQCDm.getVal() << " +/- " << nAntiQCDm.getPropagatedError(*fitResm) << endl;
    txtfile << "   Other: " << nAntiEWKm.getVal() << " +/- " << nAntiEWKm.getPropagatedError(*fitResm) << endl;
    txtfile << "  Signal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResm) << endl;
    txtfile << "     QCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResm) << endl;
    txtfile << "   Other: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResm) << endl;
    txtfile << endl;
    txtfile.flags(flags);
    
    fitResm->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
    txtfile << endl;
    fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
    txtfile << endl;
    fitResAntim->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
    txtfile << endl;
    fitResAntip->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
    txtfile << endl;
    printCorrelations(txtfile, fitResm);
    txtfile << endl;
    printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
    txtfile.close();

  }
  
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
    std::cout << "data " << hData->GetBinContent(ibin) << std::endl;
    std::cout << "mc " << hFit->GetBinContent(ibin) << std::endl;
    std::cout << "diff " << diff << std::endl;
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
//     if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    hDiff->SetBinContent(ibin,diff);
//     else      hDiff->SetBinContent(ibin,0);
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
