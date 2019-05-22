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
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"

// #include "ZBackgrounds.hh"

//helper class to handle rochester corrections
// #include <rochcor2015r.h>
// #include <muresolution_run2r.h>
#include <../RochesterCorr/RoccoR.cc>

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
#include "RooMinuit.h"
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

void fitWm_lumis_2d(const TString  outputDir,   // output directory 
           const Double_t lumi,        // integrated luminosity (/fb)'
           const Double_t lumi2, // lumi for the anti-isolation trigger
       const Double_t nsigma=0,     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
       const TString input_section = "1"
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================     

  // file format for output plots
  const TString format("png"); 
  
  bool doMTCut = false;
  bool doMET = true;
  bool doTemplate = false;

  bool pileupUp = false;
  bool pileupDown = false;
  double yscale=0.5;
  
  // Double_t vIsoBins[] = {0.0,0.15,0.45,0.55,0.65,0.75};
  Double_t vIsoBins[] = {0.0,0.2,0.25,0.35,0.45,0.55,0.65,0.75};
  // Double_t vIsoBins[] = {0.0,0.15,0.65,1.5};
  // Double_t vIsoBins[] = {0.0,0.25,0.35,0.45,0.55,0.65,0.75};
  // Double_t vIsoBins[] = {0.0,0.25,0.35,0.45,0.55};
  int nIsoBins = sizeof(vIsoBins)/sizeof(vIsoBins[0])-1;
  std::cout << "size of isobin array is " << nIsoBins << std::endl;
  
  // MET histogram binning and range
  // const Int_t    NBINS   = 75*4;
  // const Int_t    NBINS   = 125;
  // const Int_t    NBINS   = 50;
  // const Int_t    NBINS   = 100;
  const Int_t    NBINS   = 75;
  const Double_t METMIN  = 0;
  const Double_t METMAX  = 150;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t MT_CUT = 50.0;
  
  const Double_t mu_MASS = 0.1057;

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);

  
  // not sure i'll be needing these?
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
  enum { eData, eWlnu, eEWK, eBKG, eZxx, eWx, eTtb, eDib, eQCD, eAntiData, eAntiWlnu, eAntiEWK, eAntiQCD, eAntiTtb, eAntiDib, eAntiWx, eAntiZxx };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/wm_select.raw.root");   typev.push_back(eWlnu);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/wx_select.raw.root");  typev.push_back(eWx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/zxx_select.raw.root");  typev.push_back(eZxx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/zz_select.raw.root");  typev.push_back(eDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/ww_select.raw.root");  typev.push_back(eDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/wz_select.raw.root");  typev.push_back(eDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wlnu/ntuples/top_select.raw.root");  typev.push_back(eTtb);

  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/data_select.root"); typev.push_back(eAntiData);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/wx_select.root"); typev.push_back(eAntiWx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/zxx_select.root"); typev.push_back(eAntiZxx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/ww_select.root"); typev.push_back(eAntiDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/wz_select.root"); typev.push_back(eAntiDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/zz_select.root"); typev.push_back(eAntiDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/wm_select.root"); typev.push_back(eAntiWlnu);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWlnu/ntuples/top_select.root");  typev.push_back(eAntiTtb);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet   = new TH1D("hDataMet","",  NBINS,METMIN,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm  = new TH1D("hDataMetm","", NBINS,METMIN,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp  = new TH1D("hDataMetp","", NBINS,METMIN,METMAX); hDataMetp->Sumw2();
  TH1D *hWlnuMet  = new TH1D("hWlnuMet","", NBINS,METMIN,METMAX); hWlnuMet->Sumw2();
  TH1D *hWlnuMetp = new TH1D("hWlnuMetp","",NBINS,METMIN,METMAX); hWlnuMetp->Sumw2();
  TH1D *hWlnuMetm = new TH1D("hWlnuMetm","",NBINS,METMIN,METMAX); hWlnuMetm->Sumw2();
  
  TH1D *hQCDMet    = new TH1D("hQCDMet", "",  NBINS,METMIN,METMAX); hQCDMet->Sumw2();
  TH1D *hQCDMetp   = new TH1D("hQCDMetp", "", NBINS,METMIN,METMAX); hQCDMetp->Sumw2();
  TH1D *hQCDMetm   = new TH1D("hQCDMetm", "", NBINS,METMIN,METMAX); hQCDMetm->Sumw2();
  
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,METMIN,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,METMIN,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,METMIN,METMAX); hEWKMetm->Sumw2();
  
  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet", "", NBINS,METMIN,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,METMIN,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,METMIN,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWlnuMet   = new TH1D("hAntiWlnuMet", "", NBINS,METMIN,METMAX); hAntiWlnuMet->Sumw2();
  TH1D *hAntiWlnuMetp  = new TH1D("hAntiWlnuMetp","", NBINS,METMIN,METMAX); hAntiWlnuMetp->Sumw2();
  TH1D *hAntiWlnuMetm  = new TH1D("hAntiWlnuMetm","", NBINS,METMIN,METMAX); hAntiWlnuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet",  "", NBINS,METMIN,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,METMIN,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,METMIN,METMAX); hAntiEWKMetm->Sumw2();
  
  TH1D *hAntiQCDMet    = new TH1D("hAntiQCDMet", "",  NBINS,METMIN,METMAX); hAntiQCDMet->Sumw2();
  TH1D *hAntiQCDMetp   = new TH1D("hAntiQCDMetp", "", NBINS,METMIN,METMAX); hAntiQCDMetp->Sumw2();
  TH1D *hAntiQCDMetm   = new TH1D("hAntiQCDMetm", "", NBINS,METMIN,METMAX); hAntiQCDMetm->Sumw2();
  
  TH1D *hMuonEtaDatap = new TH1D("hMuonEtaDatap","",25,0,2.5); hMuonEtaDatap->Sumw2();
  TH1D *hMuonEtaDatam = new TH1D("hMuonEtaDatam","",25,0,2.5); hMuonEtaDatam->Sumw2();
  TH1D *hMuonEtaMCp   = new TH1D("hMuonEtaMCp",  "",25,0,2.5); hMuonEtaMCp->Sumw2();
  TH1D *hMuonEtaMCm   = new TH1D("hMuonEtaMCm",  "",25,0,2.5); hMuonEtaMCm->Sumw2();
  TH1D *hMuonEtaAntiDatap = new TH1D("hMuonEtaAntiDatap","",25,0,2.5); hMuonEtaAntiDatap->Sumw2();
  TH1D *hMuonEtaAntiDatam = new TH1D("hMuonEtaAntiDatam","",25,0,2.5); hMuonEtaAntiDatam->Sumw2();
  
  TH1D *hDataMetpPhi   = new TH1D("hDataMetpPhi","",  100,-3.15, 6.30); hDataMetpPhi->Sumw2();
  TH1D *hDataMetmPhi   = new TH1D("hDataMetmPhi","",  100,-3.15, 6.30); hDataMetmPhi->Sumw2();
  TH1D *hWlnuMetpPhi   = new TH1D("hWlnuMetpPhi","",  100,-3.15, 6.30); hWlnuMetpPhi->Sumw2();
  TH1D *hWlnuMetmPhi   = new TH1D("hWlnuMetmPhi","",  100,-3.15, 6.30); hWlnuMetmPhi->Sumw2();
  
  
  
  // TH1D **hAntiDataMetIsoBins   = new TH1D*[5];// hAntiDataMet->Sumw2();
  TH1D **hDataMetm2d  = new TH1D*[nIsoBins];// hAntiDataMetm->Sumw2();  
  TH1D **hDataMetp2d  = new TH1D*[nIsoBins];// hAntiDataMetp->Sumw2();
  // TH1D **hAntiWlnuMetIsoBins   = new TH1D*[5];// hAntiWlnuMet->Sumw2();
  TH1D **hWlnuMetp2d  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D **hWlnuMetm2d  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // TH1D **hAntiEWKMetIsoBins    = new TH1D*[5];// hAntiEWKMet->Sumw2();
  TH1D **hEWKMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hEWKMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hDibMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hDibMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hTtbMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hTtbMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hWxMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hWxMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hZxxMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hZxxMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hMetpIsoValues = new TH1D*[nIsoBins];
  TH1D **hMetmIsoValues = new TH1D*[nIsoBins];
  // Create a histogram pointer in each space in the array
  for(int i = 0; i < nIsoBins; i++){
    hDataMetm2d[i]  = new TH1D(("hDataMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDataMetp2d[i]  = new TH1D(("hDataMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWlnuMetp2d[i]  = new TH1D(("hWlnuMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWlnuMetm2d[i]  = new TH1D(("hWlnuMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hEWKMetp2d[i]  = new TH1D(("hEwkMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hEWKMetm2d[i]  = new TH1D(("hEwkMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	  hDibMetp2d[i]  = new TH1D(("hDibMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDibMetm2d[i]  = new TH1D(("hDibMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	  hTtbMetp2d[i]  = new TH1D(("hTtbMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hTtbMetm2d[i]  = new TH1D(("hTtbMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
  	hWxMetp2d[i]  = new TH1D(("hWxMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWxMetm2d[i]  = new TH1D(("hWxMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	  hZxxMetp2d[i]  = new TH1D(("hZxxMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hZxxMetm2d[i]  = new TH1D(("hZxxMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 

    if(i==0){
     hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],0.15);
     hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],0.15);
    }else{
	   hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
     hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
    }
  }
  


  double tolerance = ROOT::Math::MinimizerOptions::DefaultTolerance();
  string algo = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
  string type = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
  // ROOT::Math::MinimizerOptions::SetStrategy(2);
  int strategy= ROOT::Math::MinimizerOptions::DefaultStrategy();
  // int strategy= ROOT::Math::MinimizerOptions::Strategy();

  int precision= ROOT::Math::MinimizerOptions::DefaultPrecision();
  int MaxFunctionCalls= ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls();
  int MaxIterations= ROOT::Math::MinimizerOptions::DefaultMaxIterations();

  cout << "DEFAULTS: algo " << algo.c_str() << " type " << type.c_str() << " tolerance " << tolerance << " strategy " << strategy << " precision " << precision << " MaxIterations " << MaxIterations << " MaxFunctionCalls " << MaxFunctionCalls << endl;
  

//   
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum, npv, npu;
  Float_t genVPt, genVPhi, genVy, genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0;//, *lep_raw=0, *genV=0, *genLep=0;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Double_t metCorrLep, metCorrMain, metCorrEta, metCorrStat, metCorrKeys, mtCorr;
  Double_t metCorrLepPhi, metCorrMainPhi, metCorrEtaPhi, metCorrStatPhi, metCorrKeysPhi;
  Double_t totalEvtWeight=1, effSFweight=1, relIso;
    
  
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
    intree->SetBranchAddress("genVy",    &genVy);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("prefireWeight", &prefireWeight);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb",      &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",    &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",  &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("metCorrLep",      &metCorrLep);     // MET including lepton scale/smear
    intree->SetBranchAddress("metCorrLepPhi",   &metCorrLepPhi);  // MET phi including lepton scale/smear
    intree->SetBranchAddress("metCorrMain",     &metCorrMain);     // MET including lepton scale/smear (w/ main recoil corrs)
    intree->SetBranchAddress("metCorrMainPhi",  &metCorrMainPhi);  // MET phi including lepton scale/smear (w/ main recoil corrs)
    intree->SetBranchAddress("metCorrEta",      &metCorrEta);      // MET including lepton scale/smear (w/ eta recoil corrs)
    intree->SetBranchAddress("metCorrEtaPhi",   &metCorrEtaPhi);   // MET phi including lepton scale/smear  (w/ eta recoil corrs)
    intree->SetBranchAddress("metCorrStat",     &metCorrStat);     // MET including lepton scale/smear (w/ stat unc recoil corrs)
    intree->SetBranchAddress("metCorrStatPhi",  &metCorrStatPhi);  // MET phi including lepton scale/smear (w/ stat unc recoil corrs)
    intree->SetBranchAddress("metCorrKeys",     &metCorrKeys);     // MET including lepton scale/smear (w/ keyspdf recoil corrs)
    intree->SetBranchAddress("metCorrKeysPhi",  &metCorrKeysPhi);  // MET phi including lepton scale/smear (w/ keyspdf recoil corrs)
    // intree->SetBranchAddress(met_name.c_str(),      &met);       // MET
    // intree->SetBranchAddress(metPhi_name.c_str(),   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",        &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",           &mt);        // transverse mass
    intree->SetBranchAddress("mtCorr",       &mtCorr);        // transverse mass
    intree->SetBranchAddress(u1_name.c_str(),       &u1);        // parallel component of recoil
    intree->SetBranchAddress(u2_name.c_str(),       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("genLep",      &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV",         &genV);       // lepton 4-vector
    intree->SetBranchAddress("pfChIso",      &pfChIso);
    intree->SetBranchAddress("pfGamIso",     &pfGamIso);
    intree->SetBranchAddress("pfNeuIso",     &pfNeuIso);
    intree->SetBranchAddress("pfCombIso",    &pfCombIso);       // lepton 4-vector
    intree->SetBranchAddress("relIso",       &relIso);       // relative isolation for the lepton
  
    UInt_t iterator=15;
    // UInt_t iterator=1;
    if(typev[ifile]==eData||typev[ifile]==eAntiData)iterator=1;
    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    // for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; ientry<((int)intree->GetEntries()); ientry+=iterator) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;

// figure out later what to do
        // if(typev[ifile]==eWlnu || typev[ifile]==eWx || typev[ifile]==eZxx) {
          // // what was this
            // double bin = 0;
            // for(int i = 1; i <= hh_diff->GetNbinsX();++i){
              // if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            // }
            // double w2 = 1.0;//hh_diff->GetBinContent(bin);

      if(lep->Pt() < PT_CUT) continue;
      if(fabs(lep->Eta()) > ETA_CUT) continue;
      if(domtCorrCut&&(mtCorr<mtCorr_CUT)) continue;
          
      // set up the event weights for the MC reweighting
      Double_t weight=totalEvtWeight*lumi;
      Double_t weight2=totalEvtWeight*lumi2;    
      
      if(typev[ifile]==eData) {
        hDataMet->Fill(metCorrLep);
        if(q>0) {
          doMET ? hDataMetp->Fill(metCorrLep) : hDataMetp->Fill(mtCorr);
          hDataMetpPhi->Fill(metCorrLepPhi);
          hMuonEtaDatap->Fill(fabs(lep->Eta()));
          doMET ? hDataMetp2d[0]->Fill(metCorrLep) : hDataMetp2d[0]->Fill(mtCorr);
          hMetpIsoValues[0]->Fill(relIso);
        } else {
          doMET ? hDataMetm->Fill(metCorrLep) : hDataMetm->Fill(mtCorr);
          hMuonEtaDatam->Fill(fabs(lep->Eta()));
          hDataMetmPhi->Fill(metCorrLepPhi);
          doMET ? hDataMetm2d[0]->Fill(metCorrLep) : hDataMetm2d[0]->Fill(mtCorr);
          hMetmIsoValues[0]->Fill(relIso);
        }
      } else if(typev[ifile]==eAntiData) {
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            hAntiDataMet->Fill(metCorrLep);
            if(q>0) { 
              doMET ? hAntiDataMetp->Fill(metCorrLep) : hAntiDataMetp->Fill(mtCorr);
              hMuonEtaAntiDatap->Fill(fabs(lep->Eta()));
              doMET ? hDataMetp2d[it]->Fill(metCorrLep) : hDataMetp2d[it]->Fill(mtCorr);
              hMetpIsoValues[it]->Fill(relIso);
            } else { 
              doMET ? hAntiDataMetm->Fill(metCorrLep) : hAntiDataMetm->Fill(mtCorr);
              hMuonEtaAntiDatap->Fill(fabs(lep->Eta()));
              doMET ? hDataMetm2d[it]->Fill(metCorrLep) : hDataMetm2d[it]->Fill(mtCorr);
              hMetmIsoValues[it]->Fill(relIso);
              break;
            }
          }
        }
      } else if(typev[ifile]==eWlnu ) {
        hWlnuMet->Fill(metCorrLep,weight);
        if(q>0){
          hMuonEtaMCp->Fill(fabs(lep->Eta()),weight);
          hWlnuMetpPhi->Fill(metCorrLepPhi);
          doMET ? hWlnuMetp->Fill(metCorrLep,weight) : hWlnuMetp->Fill(mtCorr,weight);
          doMET ? hWlnuMetp2d[0]->Fill(metCorrLep, weight) : hWlnuMetp2d[0]->Fill(mtCorr, weight);
        } else {
          hMuonEtaMCm->Fill(fabs(lep->Eta()),weight);
          hWlnuMetmPhi->Fill(metCorrLepPhi);
          doMET ? hWlnuMetm->Fill(metCorrLep,weight) : hWlnuMetm->Fill(mtCorr,weight);
          doMET ? hWlnuMetm2d[0]->Fill(metCorrLep,weight) : hWlnuMetm2d[0]->Fill(mtCorr,weight);
        }
      } else if() {
        doMET ? hEWKMet->Fill(metCorrLep,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          doMET ? hEWKMetp->Fill(metCorrLep,weight); : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hWxMetp2d[0] ->Fill(metCorrLep,weight) : hWxMetp2d[0] ->Fill(mtCorr,weight);
          doMET ? hEWKMetp2d[0]->Fill(metCorrLep,weight) : hEWKMetp2d[0]->Fill(mtCorr,weight);
        } else {
          doMET ? hEWKMetm->Fill(metCorrLep,weight): hEWKMetm->Fill(mtCorr,weight);
          doMET ? hWxMetm2d[0] ->Fill(metCorrLep,weight) : hWxMetm2d[0] ->Fill(mtCorr,weight);
          doMET ? hEWKMetm2d[0]->Fill(metCorrLep,weight) : hEWKMetm2d[0]->Fill(mtCorr,weight);
        }
      } else if(typev[ifile]==eZxx){
        doMET ? hEWKMet->Fill(metCorrLep,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
 				  doMET ? hEWKMetp->Fill(metCorrLep,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hZxxMetp2d[0]->Fill(metCorrLep,weight) : hZxxMetp2d[0]->Fill(mtCorr,weight);
          doMET ? hEWKMetp2d[0]->Fill(metCorrLep,weight) : hEWKMetp2d[0]->Fill(mtCorr,weight);
        } else {
          doMET ? hEWKMetm->Fill(metCorrLep,weight) :  hEWKMetm->Fill(mtCorr,weight);
          doMET ? hZxxMetm2d[0]->Fill(metCorrLep,weight) : hZxxMetm2d[0]->Fill(mtCorr,weight);
          doMet ? hEWKMetm2d[0]->Fill(metCorrLep,weight) : hEWKMetm2d[0]->Fill(mtCorr,weight);
        }
      } else if(typev[ifile]==eDib) {
        doMET ? hEWKMet->Fill(metCorrLep,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          doMET ? hEWKMetp->Fill(metCorrLep,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hDibMetp2d[it]->Fill(metCorrLep,weight) : hDibMetp2d[it]->Fill(mtCorr,weight);
          doMET ? hEWKMetp2d[it]->Fill(metCorrLep,weight) : hEWKMetp2d[it]->Fill(mtCorr,weight);
        } else {
          doMET ? hEWKMetm->Fill(metCorrLep,weight) : hEWKMetm->Fill(mtCorr,weight);
          doMET ? hDibMetm2d[it]->Fill(metCorrLep,weight) : hDibMetm2d[it]->Fill(mtCorr,weight);
          doMET ? hEWKMetm2d[it]->Fill(metCorrLep,weight) : hEWKMetm2d[it]->Fill(mtCorr,weight);
        }
      } else if(typev[ifile]==eTtb) {
        doMET ? hEWKMet->Fill(metCorrLep,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          doMET ? hEWKMetp->Fill(metCorrLep,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hTtbMetp2d[it]->Fill(metCorrLep,weight) : hTtbMetp2d[it]->Fill(mtCorr,weight);
          doMET ? hEWKMetp2d[it]->Fill(metCorrLep,weight) : hEWKMetp2d[it]->Fill(mtCorr,weight);
        } else {
          doMET ? hEWKMetm->Fill(metCorrLep,weight) : hEWKMetm->Fill(mtCorr,weight);
          doMET ? hTtbMetm2d[it]->Fill(metCorrLep,weight) : hTtbMetm2d[it]->Fill(mtCorr,weight);
          doMET ? hEWKMetm2d[it]->Fill(metCorrLep,weight) : hEWKMetm2d[it]->Fill(mtCorr,weight);
        }
      } else if(typev[ifile]==eAntiWlnu){
        hAntiWlnuMet->Fill(corrMet,weight2);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0) {              
              doMET ? hAntiWlnuMetp->Fill(metCorrLep,weight2) : hAntiWlnuMetp->Fill(mtCorr,weight2);
              doMET ? hWlnuMetp2d[it]->Fill(metCorrLep, weight2) : hWlnuMetp2d[it]->Fill(mtCorr, weight2);
            } else {
              doMET ? hAntiWlnuMetm->Fill(metCorrLep,weight2) : hAntiWlnuMetm->Fill(mtCorr,weight2);
              doMET ? hWlnuMetm2d[it]->Fill(metCorrLep, weight2) : hWlnuMetm2d[it]->Fill(mtCorr, weight2);
            }
          }
        }
      } else if(typev[ifile]==eAntiWx){
        doMET ? hAntiEWKMet->Fill(metCorrLep,weight2) : hAntiEWKMet->Fill(mtCorr,weight2);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              doMET ? hAntiEWKMetp->Fill(metCorrLep,weight2) :  hAntiEWKMetp->Fill(mtCorr,weight2);
              doMET ? hWxMetp2d[it]->Fill(metCorrLep, weight2) : hWxMetp2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetp2d[it]->Fill(metCorrLep, weight2) : hEWKMetp2d[it]->Fill(mtCorr, weight2);
            } else {
              doMET ? hAntiEWKMetm->Fill(metCorrLep,weight2) : hAntiEWKMetm->Fill(mtCorr,weight2);
              doMET ? hWxMetm2d[it]->Fill(metCorrLep, weight2) : hWxMetm2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetm2d[it]->Fill(metCorrLep, weight2) : hEWKMetm2d[it]->Fill(mtCorr, weight2);
            }
            break;
          }
        }
      } else if(typev[ifile]==eAntiZxx){
        doMET ? hAntiEWKMet->Fill(metCorrLep,weight2) : hAntiEWKMet->Fill(mtCorr,weight2);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill(metCorrLep,weight2); 
              doMET ? hZxxMetp2d[it]->Fill(metCorrLep, weight2) : hZxxMetp2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetp2d[it]->Fill(metCorrLep, weight2) : hEWKMetp2d[it]->Fill(mtCorr, weight2);
            } else {
              hAntiEWKMetm->Fill(metCorrLep,weight2); 
              doMET ? hZxxMetm2d[it]->Fill(metCorrLep, weight2) : hZxxMetm2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetm2d[it]->Fill(metCorrLep, weight2) : hEWKMetm2d[it]->Fill(mtCorr, weight2);
            }
          }
        }
      } else if(typev[ifile]==eAntiDib){
        doMET ? hAntiEWKMet->Fill(metCorrLep,weight2) : hAntiEWKMet->Fill(mtCorr,weight2);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill(metCorrLep,weight2); 
              doMET ? hDibMetp2d[it]->Fill(metCorrLep, weight2) : hDibMetp2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetp2d[it]->Fill(metCorrLep, weight2) : hEWKMetp2d[it]->Fill(mtCorr, weight2);
            } else {
              hAntiEWKMetm->Fill(metCorrLep,weight2); 
              doMET ? hDibMetm2d[it]->Fill(metCorrLep, weight2) : hDibMetm2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetm2d[it]->Fill(metCorrLep, weight2) : hEWKMetm2d[it]->Fill(mtCorr, weight2);
            }
          }
        }
      } else if(typev[ifile]==eAntiTtb){
        doMET ? hAntiEWKMet->Fill(metCorrLep,weight2) : hAntiEWKMet->Fill(mtCorr,weight2);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill(metCorrLep,weight2); 
              doMET ? hTtbMetp2d[it]->Fill(metCorrLep, weight2) : hTtbMetp2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetp2d[it]->Fill(metCorrLep, weight2) : hEWKMetp2d[it]->Fill(mtCorr, weight2);
            } else {
              hAntiEWKMetm->Fill(metCorrLep,weight2); 
              doMET ? hTtbMetm2d[it]->Fill(metCorrLep, weight2) : hTtbMetm2d[it]->Fill(mtCorr, weight2);
              doMET ? hEWKMetm2d[it]->Fill(metCorrLep, weight2) : hEWKMetm2d[it]->Fill(mtCorr, weight2);
            }
          }
        }
      } else if(typev[ifile]==eAntiQCD) {
        hAntiQCDMet->Fill(metCorrLep,weight2);
        q>0 ? hAntiQCDMetp->Fill(metCorrLep,weight2) : hAntiQCDMetm->Fill(metCorrLep,weight2); 
      }
    }
  }
  delete infile;
  infile=0, intree=0;   
   
    ofstream txtfile2;
    char txtfname2[100];
    std::cout << "Printing We+. " << std::endl;
    sprintf(txtfname2,"%s/isoAvg.txt",CPlot::sOutDir.Data());
    txtfile2.open(txtfname2);
    assert(txtfile2.is_open());
 
     // flags = txtfile.flags();
    // txtfile << setprecision(10);
    // txtfile << " *** Yields *** " << endl;
   char plotname[100];
  TCanvas *c2 = MakeCanvas("c2","c2",800,600);
  for(int j = 0; j < nIsoBins; ++j){

     int binmaxWp = hMetpIsoValues[j]->GetMaximumBin(); 
	 double medianWp = hMetpIsoValues[j]->GetXaxis()->GetBinCenter(binmaxWp);
	 
	 int binmaxWm = hMetmIsoValues[j]->GetMaximumBin(); 
	 double medianWm = hMetmIsoValues[j]->GetXaxis()->GetBinCenter(binmaxWm);
	 
     txtfile2 << "bin+ : " << j << " ; average : " << hMetpIsoValues[j]->GetMean() << " ; median: " << medianWp <<  "  entries: " << hMetpIsoValues[j]->GetEntries() << endl;//<< " low edge " << hMetpIsoValues[j]->GetBinLowEdge(0) << " low edge " << hMetpIsoValues[j]->GetBinLowEdge(1000) << endl;
     txtfile2 << "bin- : " << j << " ; average : " << hMetmIsoValues[j]->GetMean() << " ; median: " << medianWm << "  entries: " << hMetmIsoValues[j]->GetEntries() << endl;
	 sprintf(plotname,"%s/hMetpIso_%d.png",CPlot::sOutDir.Data(),j);
     c2->Clear();//hMetpIsoValues[j]->Draw();
     if(j==0){
		 c2->SetLogy(1);
		 }else{
	     c2->SetLogy(0);
		 }
	 hMetpIsoValues[j]->Draw();
	 c2->Update();c2->SaveAs(plotname);
	 sprintf(plotname,"%s/hMetmIso_%d.png",CPlot::sOutDir.Data(),j);
     c2->Clear();
	  // if(j==0){
		 // c2->SetLogy(kTrue);
		 // }else{
	     // hMetmIsoValues[j]->GetYaxis()->SetLogy(kFalse);
		 // }
	 hMetmIsoValues[j]->Draw();c2->Update();c2->SaveAs(plotname);
  }
  
    txtfile2.close();
    sprintf(plotname,"%s/hMetpPhi_SR.png",CPlot::sOutDir.Data());
    c2->Clear();
    double norm = hDataMetpPhi->Integral()/hWlnuMetpPhi->Integral();
    hWlnuMetpPhi->Scale(norm);
    hWlnuMetpPhi->Draw();
    hDataMetpPhi->Draw("same");
	c2->Update();c2->SaveAs(plotname);
    sprintf(plotname,"%s/hMetmPhi_SR.png",CPlot::sOutDir.Data());
    c2->Clear();
    norm = hDataMetmPhi->Integral()/hWlnuMetmPhi->Integral();
    hWlnuMetmPhi->Scale(norm);
    hWlnuMetmPhi->Draw();
    hDataMetmPhi->Draw("same");
	c2->Update();c2->SaveAs(plotname);
    
  
  ofstream txtfile3;
  char txtfname3[100];
  std::cout << "Printing We+. " << std::endl;
  sprintf(txtfname3,"%s/isoBins_Yields.txt",CPlot::sOutDir.Data());
  txtfile3.open(txtfname3);
  assert(txtfile3.is_open());
	txtfile3 << "W+ " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  Data# " << hDataMetp2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  Wsig# " << hWlnuMetp2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  ewk# "  << hEWKMetp2d[i]->Integral()  << std::endl;
	}
	
	txtfile3 << "W- " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  Data# " << hDataMetm2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  Wsig# " << hWlnuMetm2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  ewk# "  << hEWKMetm2d[i]->Integral()  << std::endl;
	}

	txtfile3.close();
  
  
  std::cout << " ==== nEvents === " << std::endl;
  std::cout << "sig single  = " << hWlnuMetm->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWlnuMetm2d[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetm->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetm2d[0]->Integral() << std::endl;
  
  std::cout << "sig single  = " << hWlnuMetp->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWlnuMetp2d[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetp->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetp2d[0]->Integral() << std::endl;
  
  
 
 std::cout << "blah" << std::endl;
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWlnuMet->Integral());
//   cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWlnuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWlnuMet->Integral());
  dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",0.7*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWlnuMetp->Integral(),0,hAntiDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",hDataMetp->Integral()*0.3,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWlnuMetp->Integral());
//   cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  //RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWlnuMetp->Integral()*1.0,0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWlnuMetp->Integral());
  dewkp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",0.7*(hDataMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm->Integral()*0.3,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWlnuMetm->Integral());
//   cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWlnuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWlnuMetm->Integral());
  dewkm.setConstant(kTRUE);
  RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));
  
   vector<RooRealVar*> nSigp_(nIsoBins), nQCDp_(nIsoBins), dewkp_(nIsoBins), cewkp_(nIsoBins);
   vector<RooRealVar*> nSigm_(nIsoBins), nQCDm_(nIsoBins), dewkm_(nIsoBins), cewkm_(nIsoBins);
   vector<RooRealVar*> nEWKp_(nIsoBins), nEWKm_(nIsoBins);
   vector<RooRealVar*> nTtbp_(nIsoBins), nTtbm_(nIsoBins);
   vector<RooRealVar*> nDibp_(nIsoBins), nDibm_(nIsoBins);
   vector<RooRealVar*> dwxp_(nIsoBins),  dwxm_(nIsoBins);
   vector<RooFormulaVar*> nWxp_(nIsoBins),  nWxm_(nIsoBins);
   vector<RooRealVar*> nZxxp_(nIsoBins), nZxxm_(nIsoBins);
   std::cout << "blah 2" << std::endl;
   // vector<RooFormulaVar*> nEWKp_(nIsoBins), nEWKm_(nIsoBins);
   // RooRealVar* nEWKpSR = new RooRealVar("nEWKpSR","nEWKpSR",hEWKMetp2d[0]->Integral());
   // RooRealVar* nEWKmSR = new RooRealVar("nEWKmSR","nEWKmSR",hEWKMetm2d[0]->Integral());
  char nname[50];
  char formula[50];
  for (int j = 0; j < nIsoBins; ++j){
	  // W+
	  double qcdFac = 0.9;
      if(j==0) qcdFac = 0.15;
      sprintf(nname, "nAntiSigp%d",j);
      nSigp_[j] = new RooRealVar(nname,nname,hWlnuMetp2d[j]->Integral(),0,hDataMetp2d[j]->Integral());      
      // sprintf(nname, "nAntiQCDp%d",j);
      // nQCDp_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetp2d[j]->Integral(),0,hDataMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiTtbp%d",j);
      nTtbp_[j] = new RooRealVar(nname,nname,hTtbMetp2d[j]->Integral(),0.8*hTtbMetp2d[j]->Integral(),hDataMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiDibp%d",j);
      nDibp_[j] = new RooRealVar(nname,nname,hDibMetp2d[j]->Integral(),0.8*hDibMetp2d[j]->Integral(),hDataMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiZxxp%d",j);
      nZxxp_[j] = new RooRealVar(nname,nname,hZxxMetp2d[j]->Integral(),0.8*hZxxMetp2d[j]->Integral(),hDataMetp2d[j]->Integral());
	  
	  double qcd_remainder = hDataMetp2d[j]->Integral() - hWlnuMetp2d[j]->Integral() - hEWKMetp2d[j]->Integral() ;
	  if(j==0){
         sprintf(nname, "nAntiQCDp%d",j);
		 // nQCDp_[j] = new RooRealVar(nname,nname,0.4*hDataMetp2d[j]->Integral(),0.0*hDataMetp2d[j]->Integral(),0.6*hDataMetp2d[j]->Integral());
         nQCDp_[j] = new RooRealVar(nname,nname,qcd_remainder,0,2.0*qcd_remainder);
         // nQCDp_[j]->setConstant(kTRUE);
	  } else {
         sprintf(nname, "nAntiQCDp%d",j);
         // nQCDp_[j] = new RooRealVar(nname,nname,hDataMetp2d[j]->Integral(),0.95*hDataMetp2d[j]->Integral(),1.05*hDataMetp2d[j]->Integral());
         nQCDp_[j] = new RooRealVar(nname,nname,hDataMetp2d[j]->Integral(),0.95*hDataMetp2d[j]->Integral(),1.05*hDataMetp2d[j]->Integral());
	     // nQCDp_[j]->setConstant(kTRUE);
      }
	  // nQCDp_[j]->setConstant(kTRUE);
	  
	  	  // the ratio between W->enu and W->tau nu
	  sprintf(nname, "dwxp%d",j);
      dwxp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dwxp_[j]->setVal(hWxMetp2d[j]->Integral()/hWlnuMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiWxp%d",j); sprintf(formula,"dwxp%d*nAntiSigp%d",j,j);
	  nWxp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*dwxp_[j],*nSigp_[j]));
	  
      // sprintf(nname, "dewkp%d",j);
      // dewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkp_[j]->setVal(hEWKMetp2d[j]->Integral()/hWlnuMetp2d[j]->Integral());
	  
      sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"dewkp%d*nAntiSigp%d",j,j);
      // // nEWKp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigp_[j],*dewkp_[j]));
	  nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetp2d[j]->Integral(),0.0,(2.0)*hEWKMetp2d[j]->Integral());
	
	  // sprintf(nname, "nAntiDibp%d",j); 
	  // nDibp_[j] = new RooRealVar(nname,nname,hDibMetp2d[j]->Integral(),0.0,(2.0)*hDibMetp2d[j]->Integral());	  
	  
	  // sprintf(nname, "nAntiTtbp%d",j); 
	  // nTtbp_[j] = new RooRealVar(nname,nname,hTtbMetp2d[j]->Integral(),0.0,(2.0)*hTtbMetp2d[j]->Integral());
	
	 
      // W-	 
	  sprintf(nname, "nAntiSigm%d",j);
      nSigm_[j] = new RooRealVar(nname,nname,hWlnuMetm2d[j]->Integral(),0,hDataMetm2d[j]->Integral());
	  
	  	  
	  sprintf(nname, "nAntiTtbm%d",j);
      nTtbm_[j] = new RooRealVar(nname,nname,hTtbMetm2d[j]->Integral(),0.8*hTtbMetm2d[j]->Integral(),hDataMetm2d[j]->Integral());
	  
	  sprintf(nname, "nAntiDibm%d",j);
      nDibm_[j] = new RooRealVar(nname,nname,hDibMetm2d[j]->Integral(),0.8*hDibMetm2d[j]->Integral(),hDataMetm2d[j]->Integral());
	  
	  sprintf(nname, "nAntiZxxm%d",j);
      nZxxm_[j] = new RooRealVar(nname,nname,hZxxMetm2d[j]->Integral(),0.8*hZxxMetm2d[j]->Integral(),hDataMetm2d[j]->Integral());

      // sprintf(nname, "nAntiQCDm%d",j);
      // nQCDm_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetp2d[j]->Integral(),0,hDataMetm2d[j]->Integral());
	  
	  qcd_remainder = hDataMetm2d[j]->Integral() - hWlnuMetm2d[j]->Integral() - hEWKMetm2d[j]->Integral() ;
	  if(j==0){
         sprintf(nname, "nAntiQCDm%d",j);
         // nQCDm_[j] = new RooRealVar(nname,nname,0.4*hDataMetm2d[j]->Integral(),0.0*hDataMetm2d[j]->Integral(),0.6*hDataMetm2d[j]->Integral());
         nQCDm_[j] = new RooRealVar(nname,nname,qcd_remainder,0,qcd_remainder*2.0);
         // nQCDm_[j]->setConstant(kTRUE);
	  } else {
		 sprintf(nname, "nAntiQCDm%d",j);
         // nQCDm_[j] = new RooRealVar(nname,nname,hDataMetm2d[j]->Integral(),0.95*hDataMetm2d[j]->Integral(),1.05*hDataMetm2d[j]->Integral());
         nQCDm_[j] = new RooRealVar(nname,nname,hDataMetm2d[j]->Integral(),0.95*hDataMetm2d[j]->Integral(),1.05*hDataMetm2d[j]->Integral());
	     // nQCDm_[j]->setConstant(kTRUE);
      }
	  
	  // nQCDm_[j]->setConstant(kTRUE);
	  
	  	  // the W->tau nu / W->e nu ratio
	  sprintf(nname, "dwxm%d",j);
      dwxm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dwxm_[j]->setVal(hWxMetm2d[j]->Integral()/hWlnuMetm2d[j]->Integral());
	  
	  sprintf(nname, "nAntiWxm%d",j); sprintf(formula,"dwxm%d*nAntiSigm%d",j,j);
	  nWxm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*dwxm_[j],*nSigm_[j]));
	  
      // sprintf(nname, "dewkm%d",j);
      // dewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkm_[j]->setVal(hEWKMetm2d[j]->Integral()/hWlnuMetm2d[j]->Integral());
	  
      sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      // // nEWKm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigm_[j],*dewkm_[j]));
	  nEWKm_[j] = new RooRealVar(nname,nname,hEWKMetm2d[j]->Integral(),0.0,(2.0)*hEWKMetm2d[j]->Integral());
	  
	  
	  std::cout << "blah 3." << j << std::endl;
	  
  }
  
   std::cout << "finished normalizations" << std::endl;
  
  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist WlnuMet ("WlnuMET", "WlnuMET", RooArgSet(pfmet),hWlnuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,WlnuMet, 1);
  RooDataHist WlnuMetp("WlnuMETp","WlnuMETp",RooArgSet(pfmet),hWlnuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,WlnuMetp,1);
  RooDataHist WlnuMetm("WlnuMETm","WlnuMETm",RooArgSet(pfmet),hWlnuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,WlnuMetm,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
//   // test using the reversed dEta, dPhi cuts as background
  RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
  RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
  RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  
  vector<RooDataHist*> WlnuMetp_(nIsoBins), WlnuMetm_(nIsoBins);
  vector<RooDataHist*> ewkMetp_(nIsoBins), ewkMetm_(nIsoBins);
  vector<RooDataHist*> qcdMetp_(nIsoBins), qcdMetm_(nIsoBins);
  vector<RooDataHist*> ttbMetp_(nIsoBins), dibMetp_(nIsoBins), wxMetp_(nIsoBins), zxxMetp_(nIsoBins);
  vector<RooDataHist*> ttbMetm_(nIsoBins), dibMetm_(nIsoBins), wxMetm_(nIsoBins), zxxMetm_(nIsoBins);
  
  vector<RooHistPdf*> pdfWep_(nIsoBins), pdfEWKp_(nIsoBins);
  vector<RooHistPdf*> pdfQCDp_(nIsoBins), pdfQCDm_(nIsoBins);
  vector<RooHistPdf*> pdfWem_(nIsoBins), pdfEWKm_(nIsoBins);
  vector<RooHistPdf*> pdfTtbp_(nIsoBins), pdfDibp_(nIsoBins), pdfWxp_(nIsoBins), pdfZxxp_(nIsoBins);
  vector<RooHistPdf*> pdfTtbm_(nIsoBins), pdfDibm_(nIsoBins), pdfWxm_(nIsoBins), pdfZxxm_(nIsoBins);
  
  
  TH1D *hDatapFirstBinSubtracted = (TH1D*) hDataMetp2d[1]->Clone("hDatapFirstBinSubtracted");
  TH1D *hDatamFirstBinSubtracted = (TH1D*) hDataMetm2d[1]->Clone("hDatamFirstBinSubtracted");
  hDatapFirstBinSubtracted->Add(hWlnuMetp2d[1],-1);
  hDatapFirstBinSubtracted->Add(hEWKMetp2d[1],-1);
  hDatamFirstBinSubtracted->Add(hWlnuMetm2d[1],-1);
  hDatamFirstBinSubtracted->Add(hEWKMetm2d[1],-1);
  
   std::cout << "making pdfs from histograms" << std::endl;
  // char hname[50];
  // char pname[50];
  for (int j = 0; j < nIsoBins; ++j){
      // signal pdfs
      sprintf(nname, "WlnuMETp%d",j); WlnuMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWlnuMetp2d[j]);
      sprintf(nname, "wep%d",j); pdfWep_[j] = new RooHistPdf(nname,nname,pfmet,*WlnuMetp_[j],1);      
      // ewk pdfs
      sprintf(nname, "ewkMETp%d",j); ewkMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hEWKMetp2d[j]);
      sprintf(nname, "ewkp%d",j); pdfEWKp_[j] = new RooHistPdf(nname,nname,pfmet,*ewkMetp_[j],1);
      
      if(j==0){
          sprintf(nname, "qcdMetp%d",j); qcdMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDatapFirstBinSubtracted);
      } else{
          sprintf(nname, "qcdMetp%d",j); qcdMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDataMetp2d[j]);
      }
      sprintf(nname, "qcdp%d",j); pdfQCDp_[j] = new RooHistPdf(nname,nname,pfmet,*qcdMetp_[j],1);
      
	  // Split EWK into W, Z, ttbar, and di-boson
	  sprintf(nname, "ttbMETp%d",j); ttbMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hTtbMetp2d[j]);
      sprintf(nname, "ttbp%d",j); pdfTtbp_[j] = new RooHistPdf(nname,nname,pfmet,*ttbMetp_[j],1);
	  
	  sprintf(nname, "dibMETp%d",j); dibMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDibMetp2d[j]);
      sprintf(nname, "dibp%d",j); pdfDibp_[j] = new RooHistPdf(nname,nname,pfmet,*dibMetp_[j],1);
	  
	  sprintf(nname, "wxMETp%d",j); wxMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWxMetp2d[j]);
      sprintf(nname, "wxp%d",j); pdfWxp_[j] = new RooHistPdf(nname,nname,pfmet,*wxMetp_[j],1);
	  
	  sprintf(nname, "zxxMETp%d",j); zxxMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hZxxMetp2d[j]);
      sprintf(nname, "zxxp%d",j); pdfZxxp_[j] = new RooHistPdf(nname,nname,pfmet,*zxxMetp_[j],1);
	  
	  // ----------------------------------------- W- ---------------------------------
	  
	  sprintf(nname, "WlnuMETm%d",j); WlnuMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWlnuMetm2d[j]);
      sprintf(nname, "wem%d",j); pdfWem_[j] = new RooHistPdf(nname,nname,pfmet,*WlnuMetm_[j],1);      
      // ewk pdfs
      sprintf(nname, "ewkMETm%d",j); ewkMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hEWKMetm2d[j]);
      sprintf(nname, "ewkm%d",j); pdfEWKm_[j] = new RooHistPdf(nname,nname,pfmet,*ewkMetm_[j],1);
	  
	  if(j==0){
          sprintf(nname, "qcdMetm%d",j); qcdMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDatamFirstBinSubtracted);
      } else{
          sprintf(nname, "qcdMetm%d",j); qcdMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDataMetm2d[j]);
      }
      sprintf(nname, "qcdm%d",j); pdfQCDm_[j] = new RooHistPdf(nname,nname,pfmet,*qcdMetm_[j],1);
      
	  // Split the EWK into W, Z, ttbar, and di-boson
	  sprintf(nname, "ttbMETm%d",j); ttbMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hTtbMetm2d[j]);
      sprintf(nname, "ttbm%d",j); pdfTtbm_[j] = new RooHistPdf(nname,nname,pfmet,*ttbMetm_[j],1);
	  
	  sprintf(nname, "dibMETm%d",j); dibMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDibMetm2d[j]);
      sprintf(nname, "dibm%d",j); pdfDibm_[j] = new RooHistPdf(nname,nname,pfmet,*dibMetm_[j],1);
	  
	  sprintf(nname, "wxMETm%d",j); wxMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWxMetm2d[j]);
      sprintf(nname, "wxm%d",j); pdfWxm_[j] = new RooHistPdf(nname,nname,pfmet,*wxMetm_[j],1);
	  
	  sprintf(nname, "zxxMETm%d",j); zxxMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hZxxMetm2d[j]);
      sprintf(nname, "zxxm%d",j); pdfZxxm_[j] = new RooHistPdf(nname,nname,pfmet,*zxxMetm_[j],1);
  }
  
  
  // QCD Pdfs
  
  //CExponential qcd(pfmet,kTRUE);
  //CExponential qcdp(pfmet,kTRUE);
  //CExponential qcdm(pfmet,kTRUE);
  // comment back in for qcd functional form
  CPepeModel2 qcd("qcd",pfmet);
  CPepeModel2 qcdp("qcdp",pfmet);
  CPepeModel2 qcdm("qcdm",pfmet);
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
    
    RooGaussian constantim("constantim","constantim",nAntiSigm,RooConst(hAntiWlnuMetm->Integral()),RooConst(0.15*hAntiWlnuMetm->Integral()));
    RooGaussian constantip("constantip","constantip",nAntiSigp,RooConst(hAntiWlnuMetp->Integral()),RooConst(0.15*hAntiWlnuMetp->Integral()));
	

 
//   // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));

//   // Signal + Background PDFs
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,pdfQCD),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,pdfQCDp),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,pdfQCDm),RooArgList(nSigm,nEWKm,nQCDm));
  
//   // constrained background?
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,qcdpc),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,qcdmc),RooArgList(nSigm,nEWKm,nQCDm));
  
  
    
  
  // Anti-Signal PDFs
  RooDataHist aWlnuMet ("aWlnuMET", "aWlnuMET", RooArgSet(pfmet),hAntiWlnuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,aWlnuMet, 1);
  RooDataHist aWlnuMetp("aWlnuMETp","aWlnuMETp",RooArgSet(pfmet),hAntiWlnuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,aWlnuMetp,1);
  RooDataHist aWlnuMetm("aWlnuMETm","aWlnuMETm",RooArgSet(pfmet),hAntiWlnuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,aWlnuMetm,1); 
  
  // Anti-EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  
//   // Anti-QCD Pdfs
  CPepeModel2 aqcd("aqcd",pfmet, qcd.a1);
  CPepeModel2 aqcdp("aqcdp",pfmet, qcdp.a1);
  CPepeModel2 aqcdm("aqcdm",pfmet, qcdm.a1);

//   CPepeModel2 aqcd("aqcd",pfmet);
//   CPepeModel2 aqcdp("aqcdp",pfmet);
//   CPepeModel2 aqcdm("aqcdm",pfmet);
  
  // Anti-selection PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
  
  // regular linear dependence
  // vector <CPepeModel2isobinsMuons*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  // qcdp_[0] = new CPepeModel2isobinsMuons("qcdp2d0",pfmet, 0.075);
  // qcdm_[0] = new CPepeModel2isobinsMuons("qcdm2d0",pfmet, 0.075);
  
  // quadratic dependence
  // vector <CPepeModel2isobinsQuad*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  // qcdp_[0] = new CPepeModel2isobinsQuad("qcdp2d0",pfmet, 0.075);
  // qcdm_[0] = new CPepeModel2isobinsQuad("qcdm2d0",pfmet, 0.075);
  
  // totally uncorrelated
  vector <CPepeModel2*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  qcdp_[0] = new CPepeModel2("qcdp2d0",pfmet);
  qcdm_[0] = new CPepeModel2("qcdm2d0",pfmet);
  
  for (int j = 1; j < nIsoBins; ++j){
	  
      // quadratic isolation dependence
	  // sprintf(nname, "qcdp2d%d",j);
	  // qcdp_[j] = new CPepeModel2isobinsQuad(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->b1, qcdp_[0]->b2, qcdp_[0]->b3, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  // printf(nname, "qcdm2d%d",j);
	  // qcdm_[j] = new CPepeModel2isobinsQuad(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2,  qcdm_[0]->b1, qcdm_[0]->b2, qcdm_[0]->b3, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
	  
	  // // Original isolation linear dependence
      // sprintf(nname, "qcdp2d%d",j);
      // qcdp_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  // sprintf(nname, "qcdm2d%d",j);
      // qcdm_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
  
      sprintf(nname, "qcdp2d%d",j);
      qcdp_[j] = new CPepeModel2(nname,pfmet);
	  
	  sprintf(nname, "qcdm2d%d",j);
      qcdm_[j] = new CPepeModel2(nname,pfmet);
  
  }
   std::cout << "just finished making the pepes" << std::endl;
   
   // create the appropriate gaussian constraint for the signal region electroweak+ttbar 
    RooGaussian constm_sr("constm_sr","constm_sr",*nEWKm_[0],RooConst(nEWKm_[0]->getVal()),RooConst(0.10*nEWKm_[0]->getVal()));
    RooGaussian constp_sr("constp_sr","constp_sr",*nEWKp_[0],RooConst(nEWKp_[0]->getVal()),RooConst(0.10*nEWKp_[0]->getVal()));
	
	RooGaussian constm_cr("constm_cr","constm_cr",*nEWKm_[1],RooConst(nEWKm_[1]->getVal()),RooConst(0.10*nEWKm_[1]->getVal()));
    RooGaussian constp_cr("constp_cr","constp_cr",*nEWKp_[1],RooConst(nEWKp_[1]->getVal()),RooConst(0.10*nEWKp_[1]->getVal()));
	
	// RooGaussian constm_sr("constm_sr","constm_sr",*dewkm_[0],RooConst(dewkm_[0]->getVal()),RooConst(0.15*dewkm_[0]->getVal()));
    // RooGaussian constp_sr("constp_sr","constp_sr",*dewkp_[0],RooConst(dewkp_[0]->getVal()),RooConst(0.15*dewkp_[0]->getVal()));
	
	
	     std::cout << "nsigp " << nSigp_[0]->getVal() << std::endl;
		 std::cout << "nttbp " << nTtbp_[0]->getVal() << std::endl;
		 std::cout << "nqcdp " << nQCDp_[0]->getVal() << std::endl;
		 std::cout << "nwxp " << nWxp_[0]->getVal() << std::endl;
		 std::cout << "nzxxp " << nZxxp_[0]->getVal() << std::endl;
	
		// Fill the actual values later
	// Wx sample - constrain dewk values w/ Gaussian
	RooGaussian const_wxm("const_wxm","const_wxm",*dwxm_[0],RooConst(dwxm_[0]->getVal()),RooConst(0.009*dwxm_[0]->getVal()));
    RooGaussian const_wxp("const_wxp","const_wxp",*dwxp_[0],RooConst(dwxp_[0]->getVal()),RooConst(0.009*dwxp_[0]->getVal()));
	
	RooGaussian const_ttbm("const_ttbm","const_ttbm",*nTtbm_[0],RooConst(nTtbm_[0]->getVal()*815/830),RooConst(0.0686*nTtbm_[0]->getVal()));
    RooGaussian const_ttbp("const_ttbp","const_ttbp",*nTtbp_[0],RooConst(nTtbp_[0]->getVal()*815/830),RooConst(0.0686*nTtbp_[0]->getVal()));
	
	RooGaussian const_dibm("const_dibm","const_dibm",*nDibm_[0],RooConst(nDibm_[0]->getVal()),RooConst(0.009*nDibm_[0]->getVal()));
    RooGaussian const_dibp("const_dibp","const_dibp",*nDibp_[0],RooConst(nDibp_[0]->getVal()),RooConst(0.009*nDibp_[0]->getVal()));
	
	RooGaussian const_zxxm("const_zxxm","const_zxxm",*nZxxm_[0],RooConst(nZxxm_[0]->getVal()*5688/5833),RooConst(0.0449*nZxxm_[0]->getVal()));
    RooGaussian const_zxxp("const_zxxp","const_zxxp",*nZxxp_[0],RooConst(nZxxp_[0]->getVal()*5688/5833),RooConst(0.0449*nZxxp_[0]->getVal()));
	
	
	
	std::cout << "make constraints" << std::endl;

	// RooGaussian constm_cr1("constm_cr1","constm_cr1",*nEWKm_[1],RooConst(nEWKm_[1]->getVal()),RooConst(0.05*nEWKm_[1]->getVal()));
    // RooGaussian constp_cr1("constp_cr1","constp_cr1",*nEWKp_[1],RooConst(nEWKp_[1]->getVal()),RooConst(0.05*nEWKp_[1]->getVal()));

	// RooGaussian constm_cr2("constm_cr2","constm_cr2",*nEWKm_[2],RooConst(nEWKm_[2]->getVal()),RooConst(0.05*nEWKm_[2]->getVal()));
    // RooGaussian constp_cr2("constp_cr2","constp_cr2",*nEWKp_[2],RooConst(nEWKp_[2]->getVal()),RooConst(0.05*nEWKp_[2]->getVal()));

	// RooGaussian constm_cr3("constm_cr3","constm_cr3",*nEWKm_[3],RooConst(nEWKm_[3]->getVal()),RooConst(0.05*nEWKm_[3]->getVal()));
    // RooGaussian constp_cr3("constp_cr3","constp_cr3",*nEWKp_[3],RooConst(nEWKp_[3]->getVal()),RooConst(0.05*nEWKp_[3]->getVal()));

	// RooGaussian constm_cr4("constm_cr4","constm_cr4",*nEWKm_[4],RooConst(nEWKm_[4]->getVal()),RooConst(0.05*nEWKm_[4]->getVal()));
    // RooGaussian constp_cr4("constp_cr4","constp_cr4",*nEWKp_[4],RooConst(nEWKp_[4]->getVal()),RooConst(0.05*nEWKp_[4]->getVal()));

	// RooGaussian constm_cr5("constm_cr5","constm_cr5",*nEWKm_[5],RooConst(nEWKm_[5]->getVal()),RooConst(0.05*nEWKm_[5]->getVal()));
    // RooGaussian constp_cr5("constp_cr5","constp_cr5",*nEWKp_[5],RooConst(nEWKp_[5]->getVal()),RooConst(0.05*nEWKp_[5]->getVal()));

	
     // contstuct the PDFs for the 2 test isolation bins
  // Test the vector map 
  // // map <string, RooAddPdf> pdfMetp_;
   // map <int, RooAddPdf*> pdfMetp_;
   vector < RooAddPdf*> pdfMetp_(nIsoBins), pdfMetm_(nIsoBins);
  // make the pdfs?
  for(int j = 0; j < nIsoBins; ++j){
	 // if(j==0 || j==1){
	 if(j==0||j==1||j==2||j==3||j==4){
		 sprintf(nname,"pdfWep%d",j);

		 if(doMET&&(!doTemplate)){
             pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j]));
         } else {
		     pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*pdfQCDp_[j]),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j]));
         }
		 // pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfWxp_[j],*pdfZxxp_[j],*pdfDibp_[j],*pdfTtbp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nWxp_[j],*nZxxp_[j],*nDibp_[j],*nTtbp_[j],*nQCDp_[j]));
		 std::cout << "blah" << std::endl;
		 
		 sprintf(nname,"pdfWem%d",j);
		 if(doMET&&(!doTemplate)){
             pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j]));
         }else {
             pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*pdfQCDm_[j]),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j]));
         }		
        // pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfWxm_[j],*pdfZxxm_[j],*pdfDibm_[j],*pdfTtbm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nWxm_[j],*nZxxm_[j],*nDibm_[j],*nTtbm_[j],*nQCDm_[j]));
		 // sprintf(nname,"pdfWep%d",j);
		 // pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfDibp_[j],*pdfTtbp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nDibp_[j],*nTtbp_[j],*nQCDp_[j]));
		 
		 // sprintf(nname,"pdfWem%d",j);
		 // pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfDibm_[j],*pdfTtbm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nDibm_[j],*nTtbm_[j]*nQCDm_[j]));
	 } else {
         sprintf(nname,"pdfWep%d",j);
		 if(doMET&&(!doTemplate)){
             pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*(qcdp_[j]->model)),RooArgList(*nQCDp_[j]));
		 } else {
             pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfQCDp_[j]),RooArgList(*nQCDp_[j]));
		 }
         
		 sprintf(nname,"pdfWem%d",j);
		 if (doMET&&(!doTemplate) ) {
             pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*(qcdm_[j]->model)),RooArgList(*nQCDm_[j]));
         } else {
		     pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfQCDm_[j]),RooArgList(*nQCDm_[j]));
         }
		 
	 }
  }
  
  std::cout << "made pdfs" << std::endl;
  
  
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


  
    // ------------------------------------------------------------------------
  // Make the data histograms for the RooCategory Fit
  
  // put them into a vector map, make sure the names match the RooCategory names
  map<string, RooDataHist*> dataMetp_, dataMetm_;
  std::vector<RooDataHist*> dataMetpHist_(nIsoBins),dataMetmHist_(nIsoBins);
  // map<string, RooDataHist*> dataMetm_;
  for(int i = 0; i < nIsoBins; ++i){
	  sprintf(nname, "isop%d",i);
	  dataMetp_[nname]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetp2d[i]);
      sprintf(nname, "isom%d",i);
	  dataMetm_[nname]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetm2d[i]);
      sprintf(nname, "isop%d",i);
      dataMetpHist_[i]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetp2d[i]);
	  sprintf(nname, "isom%d",i);
      dataMetmHist_[i]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetm2d[i]);
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  // set up RooCategory for data
  // test only 2 bins at first
  
  RooCategory rooCat2dTest("rooCat2dTest","rooCat2dTest");
  for (int i = 0; i < nIsoBins; ++ i){
	  sprintf(nname,"iso%d",i); rooCat2dTest.defineType(nname);
  }

  RooDataHist combDatap("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetp_));
  RooDataHist combDatam("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetm_));
  
  RooSimultaneous simPdfp("simPdfp","simultaneous pdf W+",rooCat2dTest) ;
  RooSimultaneous simPdfm("simPdfm","simultaneous pdf W-",rooCat2dTest) ;
  for(int i = 0; i < nIsoBins; ++i ){
	  sprintf(nname, "isop%d", i); simPdfp.addPdf(*pdfMetp_[i],nname);
	  sprintf(nname, "isom%d", i); simPdfm.addPdf(*pdfMetm_[i],nname);
  }

  
  cout << "Starting values for Wlnu yields: " << endl;
  cout << "Selected: " << hDataMet->Integral() << endl;
  cout << "   sig: " << hWlnuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;
  cout << "   qcd: " << hDataMet->Integral()-hWlnuMet->Integral()-hEWKMet->Integral() << endl;

  cout << "Starting values for Wlnu_p yields: " << endl;
  cout << "   sig: " << hWlnuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hDataMetp->Integral()-hWlnuMetp->Integral()-hEWKMetp->Integral() << endl;

  cout << "Starting values for Wlnu_m yields: " << endl;
  cout << "   sig: " << hWlnuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  cout << "   qcd: " << hDataMetm->Integral()-hWlnuMetm->Integral()-hEWKMetm->Integral() << endl;
  
  
  cout << "Starting values for AntiWlnu yields: " << endl;
  cout << "Selected: " << hAntiDataMet->Integral() << endl;
  cout << "   sig: " << hAntiWlnuMet->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMet->Integral() << endl;
  cout << "   qcd: " << hAntiDataMet->Integral()-hAntiWlnuMet->Integral()-hAntiEWKMet->Integral() << endl;

  cout << "Starting values for AntiWlnu_p yields: " << endl;
  cout << "   sig: " << hWlnuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetp->Integral()-hAntiWlnuMetp->Integral()-hAntiEWKMetp->Integral() << endl;

  cout << "Starting values for AntiWlnu_m yields: " << endl;
  cout << "   sig: " << hAntiWlnuMetm->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetm->Integral()-hAntiWlnuMetm->Integral()-hAntiEWKMetm->Integral() << endl;

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());

  RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  
  RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  
  
  
  std::vector<RooRealVar*> pepe2Pdf_qcdp_norms_(nIsoBins), pepe2Pdf_qcdm_norms_(nIsoBins);
  for(int j = 0; j < nIsoBins; ++j){
     // sprintf(nname,"pepe2Pdf_qcdp_norm");
     if(j==0){
       sprintf(nname,"pepe2Pdf_qcdp2d%i_norm",j);
       pepe2Pdf_qcdp_norms_[j] = new RooRealVar(nname,nname,0.4*hDataMetp2d[j]->Integral(),0,0.6*hDataMetp2d[j]->Integral());
       sprintf(nname,"pepe2Pdf_qcdm2d%i_norm",j);
       pepe2Pdf_qcdm_norms_[j] = new RooRealVar(nname,nname,0.4*hDataMetm2d[j]->Integral(),0,0.6*hDataMetm2d[j]->Integral());
     } else {
       sprintf(nname,"pepe2Pdf_qcdp2d%i_norm",j);
       pepe2Pdf_qcdp_norms_[j] = new RooRealVar(nname,nname,hDataMetp2d[j]->Integral(),0.95*hDataMetp2d[j]->Integral(),1.05*hDataMetp2d[j]->Integral());
       sprintf(nname,"pepe2Pdf_qcdm2d%i_norm",j);
       pepe2Pdf_qcdm_norms_[j] = new RooRealVar(nname,nname,hDataMetm2d[j]->Integral(),0.95*hDataMetm2d[j]->Integral(),1.05*hDataMetm2d[j]->Integral());
     }
  }
  std::cout << "mark 3 " << std::endl;
  
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
//   combine_workspace.import(*(qcd.model));
//   combine_workspace.import(*(qcdp.model));
//   combine_workspace.import(*(qcdm.model));

  combine_workspace.import(pdfQCD);
  combine_workspace.import(pdfQCDp);
  combine_workspace.import(pdfQCDm);
  
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

  sprintf(nname,"%s/Wlnu_pdfTemplates.root",CPlot::sOutDir.Data());
  combine_workspace.writeToFile(nname);
  
  
   // separate file for the binned workspace just so i don't fuck up the original one
  RooWorkspace binned_workspace("binned_workspace");
  // loop through the number of bins, import the appropriate pdf for each one
  for(int j = 0; j < nIsoBins; ++j){
      // binned_workspace.import(dataMet);
      binned_workspace.import(*dataMetpHist_[j]);
      binned_workspace.import(*dataMetmHist_[j]);
	  // QCD normalization RooRealVars
      binned_workspace.import(*pepe2Pdf_qcdp_norms_[j]);
      binned_workspace.import(*pepe2Pdf_qcdm_norms_[j]);
	  // QCD shapes
      // binned_workspace.import(*(qcd.model));
      binned_workspace.import(*(qcdp_[j]->model));
      binned_workspace.import(*(qcdm_[j]->model));
      std::cout << "importing   " << qcdm_[j]->model->GetName() << std::endl;
	  //The DataHist
	  // combine_workspace.import(*wenuMet_[j]);
	  binned_workspace.import(*WlnuMetp_[j]);
	  binned_workspace.import(*WlnuMetm_[j]);
      // binned_workspace.import(pdfWe);
      // binned_workspace.import(*pdfWep_[j]);
      // binned_workspace.import(*pdfWem_[j]);
	  // The DataHist
	  // combine_workspace.import(*ewkMet_[j]);
	  binned_workspace.import(*ewkMetp_[j]);
	  binned_workspace.import(*ewkMetm_[j]);
      
    binned_workspace.import(*wxMetp_[j]);
	  binned_workspace.import(*wxMetm_[j]);
      
    binned_workspace.import(*zxxMetp_[j]);
	  binned_workspace.import(*zxxMetm_[j]);
      
    binned_workspace.import(*ttbMetp_[j]);
	  binned_workspace.import(*ttbMetm_[j]);
      
    binned_workspace.import(*dibMetp_[j]);
	  binned_workspace.import(*dibMetm_[j]);
	  //The Pdfs
      // binned_workspace.import(pdfEWK);
      // binned_workspace.import(*pdfEWKp_[j]);
      // binned_workspace.import(*pdfEWKm_[j]);
	  
	  binned_workspace.import(*hDataMetp2d[j]);
	  binned_workspace.import(*hDataMetm2d[j]);
	  binned_workspace.import(*hWlnuMetp2d[j]);
	  binned_workspace.import(*hWlnuMetm2d[j]);
	  binned_workspace.import(*hEWKMetp2d[j]);
	  binned_workspace.import(*hEWKMetm2d[j]);
    binned_workspace.import(*hWxMetp2d[j]);
	  binned_workspace.import(*hWxMetm2d[j]);
    binned_workspace.import(*hZxxMetp2d[j]);
	  binned_workspace.import(*hZxxMetm2d[j]);
    binned_workspace.import(*hDibMetp2d[j]);
	  binned_workspace.import(*hDibMetm2d[j]);
    binned_workspace.import(*hTtbMetp2d[j]);
	  binned_workspace.import(*hTtbMetm2d[j]);
	  // binned_workspace.import(*hDataMetp2d[j]);
	  

  }// end of loop, now save the workspace
  sprintf(nname, "%s/Wlnu_pdfTemplates_binned.root",CPlot::sOutDir.Data());
  binned_workspace.writeToFile(nname);
  
  

  
  
     // had commented out the Min2/Strat2 when running some of the other options
  // RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),RooFit::Strategy(2)/*,Minimizer("Minuit2","minimize")*/,ExternalConstraints(RooArgList(const_wxp,const_zxxp,const_dibp,const_ttbp)),PrintEvalErrors(-1));
  // RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),RooFit::Strategy(2)/*,Minimizer("Minuit2","minimize")*/,ExternalConstraints(RooArgList(const_wxm,const_zxxm,const_dibm,const_ttbm)),PrintEvalErrors(-1));
  
  RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),ExternalConstraints(RooArgSet(constp_sr/*,constp_cr*/)),RooFit::Strategy(2),Minos(kTRUE),/*Minimizer("Minuit2","minimize"),*/PrintEvalErrors(-1));
  RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),ExternalConstraints(RooArgSet(constm_sr/*,constm_cr*/)),RooFit::Strategy(2),Minos(kTRUE),/*Minimizer("Minuit2","minimize"),*/PrintEvalErrors(-1));
  
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWlnuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWlnuMetm->GetBinError(ibin));}
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
   for(int ibin = 1; ibin < hPdfAntiMetp->GetNbinsX(); ++ibin){hPdfAntiMetp->SetBinError(ibin, hAntiWlnuMetp->GetBinError(ibin));}
  hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetm->GetNbinsX(); ++ibin){hPdfAntiMetm->SetBinError(ibin, hAntiWlnuMetm->GetBinError(ibin));}
  hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  hAntiMetmDiff->SetMarkerStyle(kFullCircle);
  hAntiMetmDiff->SetMarkerSize(0.9);
    
    char plotname2[100];
    //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
    // label for lumi
  // char lumitext[100];
  // if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  (8 TeV)",lumi*1000.);
  // else         sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000.);
  
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

  
    
  Double_t chi2probm, chi2ndfm, chi2probp, chi2ndfp;
  Double_t ksprobm, ksprobpem, ksprobp, ksprobpep;
  
    ofstream txtfile4;
    char txtfname4[100];
    std::cout << "Printing We+. " << std::endl;
    sprintf(txtfname4,"%s/chi2.txt",CPlot::sOutDir.Data());
    txtfile4.open(txtfname4);
    assert(txtfile4.is_open());
  
  for(int i = 0; i < nIsoBins; ++i){
      
    // set up the diff plot
    // turn this into its own function later?
    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetp = (TH1D*)(pdfMetp_[i]->createHistogram("hPdfMetp", pfmet));
    for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWlnuMetp2d[i]->GetBinError(ibin));}
    hPdfMetp->Scale((nSigp_[i]->getVal()+nEWKp_[i]->getVal()+nQCDp_[i]->getVal())/hPdfMetp->Integral());
    TH1D *hMetpDiff = makeDiffHist(hDataMetp2d[i],hPdfMetp,"hMetpDiff");
    hMetpDiff->SetMarkerStyle(kFullCircle);
    hMetpDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    // turn this part also into its own function later
    RooPlot *wepframe = pfmet.frame(Bins(NBINS)); 
    wepframe->GetYaxis()->SetNdivisions(505);
    wepframe->GetXaxis()->SetLabelOffset(2.0);
    sprintf(nname,"isop%d",i);
    sprintf(plotname2, "wep_fitmetp_bin%i",i);
    dataMetp_[nname]->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMetp_[i]->plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMetp_[i]->plotOn(wepframe,LineColor(linecolorW));
    if(!doMET||doTemplate){
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*pdfQCDp_[i])),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*pdfQCDp_[i])),LineColor(linecolorEWK));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfQCDp_[i])),LineColor(linecolorQCD));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfQCDp_[i])),FillColor(fillcolorQCD),DrawOption("F"));
    } else {
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*(qcdp_[i]->model))),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*(qcdp_[i]->model))),LineColor(linecolorEWK));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*(qcdp_[i]->model))),FillColor(fillcolorQCD),DrawOption("F"));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*(qcdp_[i]->model))),LineColor(linecolorQCD));
    } 
    pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfWep_[i])),LineColor(linecolorW),LineStyle(2));
    dataMetp_[nname]->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
   
    sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
    CPlot plotMetp(plotname2,wepframe,"","",ylabel);
    plotMetp.SetLegend(0.68,0.57,0.93,0.77);
    plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMetp.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
    plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
    plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
    plotMetp.Draw(c,kFALSE,format,1);
    
    std::cout << "Draw the W plot diff" << std::endl;
    CPlot plotMetpDiff(plotname2,"","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
    hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
    hMetpDiff->GetYaxis()->SetLabelSize(0.11);
    plotMetpDiff.SetYRange(-yscale,yscale);
    plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
    plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
    plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
    plotMetpDiff.Draw(c,kTRUE,format,2);
    plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
    std::cout << "Draw the W plot log" << std::endl;
    sprintf(plotname2,"wep_fitmetp_bin%i_log",i);
    plotMetp.SetName(plotname2);
    plotMetp.SetLogy();
    plotMetp.SetYRange(1e-5*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
    plotMetp.Draw(c,kTRUE,format,1);
    plotMetp.Draw(c,kTRUE,"pdf",1);
   
    chi2probp = hDataMetp2d[i]->Chi2Test(hPdfMetp,"PUW");
    chi2ndfp  = hDataMetp2d[i]->Chi2Test(hPdfMetp,"CHI2/NDFUW");
    ksprobp   = hDataMetp2d[i]->KolmogorovTest(hPdfMetp);
    ksprobpep = hDataMetp2d[i]->KolmogorovTest(hPdfMetp,"DX"); 
   
    // Do the W-  plots here
    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetm = (TH1D*)(pdfMetm_[i]->createHistogram("hPdfMetm", pfmet));
    for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWlnuMetm2d[i]->GetBinError(ibin));}
    hPdfMetm->Scale((nSigm_[i]->getVal()+nEWKm_[i]->getVal()+nQCDm_[i]->getVal())/hPdfMetm->Integral());
    TH1D *hMetmDiff = makeDiffHist(hDataMetm2d[i],hPdfMetm,"hMetmDiff");
    hMetmDiff->SetMarkerStyle(kFullCircle);
    hMetmDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    // turn this part also into its own function later
    RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
    wemframe->GetYaxis()->SetNdivisions(505);
    wemframe->GetXaxis()->SetLabelOffset(2.0);
    sprintf(nname,"isom%d",i);
    sprintf(plotname2, "wem_fitmetm_bin%i",i);
    dataMetm_[nname]->plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMetm_[i]->plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMetm_[i]->plotOn(wemframe,LineColor(linecolorW));
    if(!doMET||doTemplate){
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*pdfQCDm_[i])),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*pdfQCDm_[i])),LineColor(linecolorEWK));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfQCDm_[i])),FillColor(fillcolorQCD),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfQCDm_[i])),LineColor(linecolorQCD));
    } else {
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*(qcdm_[i]->model))),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*(qcdm_[i]->model))),LineColor(linecolorEWK));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*(qcdm_[i]->model))),FillColor(fillcolorQCD),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*(qcdm_[i]->model))),LineColor(linecolorQCD));
    }
    pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfWem_[i])),LineColor(linecolorW),LineStyle(2));
    dataMetm_[nname]->plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
   
    sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
    CPlot plotMetm(plotname2,wemframe,"","",ylabel);
    plotMetm.SetLegend(0.68,0.57,0.93,0.77);
    plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
    plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
    plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
    plotMetm.Draw(c,kFALSE,format,1);
    
    std::cout << "Draw the W plot diff" << std::endl;
    CPlot plotMetmDiff(plotname2,"","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
    hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
    hMetmDiff->GetYaxis()->SetLabelSize(0.11);
    plotMetmDiff.SetYRange(-yscale,yscale);
    plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
    plotMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
    plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
    plotMetmDiff.Draw(c,kTRUE,format,2);
    plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
    std::cout << "Draw the W plot log" << std::endl;
    sprintf(plotname2,"wem_fitmetm_bin%i_log",i);
    plotMetm.SetName(plotname2);
    plotMetm.SetLogy();
    plotMetm.SetYRange(1e-5*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
    plotMetm.Draw(c,kTRUE,format,1);
    plotMetm.Draw(c,kTRUE,"pdf",1);
	
	chi2probm = hDataMetm2d[i]->Chi2Test(hPdfMetm,"PUW");
    chi2ndfm  = hDataMetm2d[i]->Chi2Test(hPdfMetm,"CHI2/NDFUW");
    ksprobm   = hDataMetm2d[i]->KolmogorovTest(hPdfMetm);
    ksprobpem = hDataMetm2d[i]->KolmogorovTest(hPdfMetm,"DX"); 
	
    txtfile4 << "bin+ : " << i << " ; chi2 : " << chi2ndfp << endl;
	txtfile4 << "bin- : " << i << " ; chi2 : " << chi2ndfm << endl;
	
  
   std::cout << "blah " << i << std::endl;
  
  }
  
  txtfile4.close();
  
  std::cout << "hi" << std::endl;
  ofstream txtfile;
  char txtfname1[100];
  std::cout << "Printing We+. " << std::endl;
  sprintf(txtfname1,"%s/fitresWe2dCatTest.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname1);
  assert(txtfile.is_open());
  
  // flags = txtfile.flags();
  // txtfile << setprecision(10);
  // txtfile << " *** Yields *** " << endl;
  txtfile << "Data+: " << hDataMetp2d[0]->Integral() << endl;
  txtfile << "Data-: " << hDataMetm2d[0]->Integral() << endl;
  // txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  // txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  // txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  // txtfile << "AntiSelected: " << hAntiDataMetp->Integral() << endl;
  // txtfile << "  AntiSignal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResp) << endl;
  // txtfile << "     AntiQCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResp) << endl;
  // txtfile << "   AntiOther: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResp) << endl;
  // txtfile << endl;  
  // txtfile.flags(flags);
  
   
  
  fitResp2dCatTest->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose); 
  txtfile << endl;
  printCorrelations(txtfile, fitResp2dCatTest);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2probp, chi2ndfp, ksprobp, ksprobpep);
  txtfile<< endl;
  
  fitResm2dCatTest->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResm2dCatTest);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2probm, chi2ndfm, ksprobm, ksprobpem);
  txtfile <<endl;
  
  txtfile.close();
  
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
    if(hData->GetBinContent(ibin) == 0) diff = 0;
    
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    if(hData->GetBinContent(ibin) == 0) err = 0;
    std::cout << "ibin " << ibin << "  data " << hData->GetBinContent(ibin) << "  fit " << hFit->GetBinContent(ibin) << "  err " << err << std::endl;
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
  ios_base::fmtCorrflags flags = os.flags();
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
  ios_base::fmtCorrflags flags = os.flags();
  
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
  sprintf(htmlfname,"%s/WlnuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wlnu</title></head>" << endl;
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
