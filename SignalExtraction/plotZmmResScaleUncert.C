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

#include "TStopwatch.h" //PROFILE

#endif

TStopwatch sw_[10]; //PROFILE 

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

//=== MAIN MACRO ================================================================================================= 

void plotZmmResScaleUncert(const TString  inputDir,    // input directory
			   const TString  outputDir,   // output directory
			   const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmmResScaleUncert");
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

   
  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg_ResScale.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");


  // plot output file format
  const TString format("png");
   
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
  TH1D *hTop  = new TH1D("hTop", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();

  TH1D *hDataZPt = new TH1D("hDataZPt","",34,ZPtBins); hDataZPt->Sumw2();
  TH1D *hZmmZPt  = new TH1D("hZmmZPt", "",34,ZPtBins); hZmmZPt->Sumw2();
  TH1D *hEWKZPt  = new TH1D("hEWKZPt", "",34,ZPtBins); hEWKZPt->Sumw2();
  TH1D *hTopZPt  = new TH1D("hTopZPt", "",34,ZPtBins); hTopZPt->Sumw2();
  TH1D *hMCZPt   = new TH1D("hMCZPt",  "",34,ZPtBins); hMCZPt->Sumw2();

  TH1D *hDataPhiStar = new TH1D("hDataPhiStar","",27,PhiStarBins); hDataPhiStar->Sumw2();
  TH1D *hZmmPhiStar  = new TH1D("hZmmPhiStar", "",27,PhiStarBins); hZmmPhiStar->Sumw2();
  TH1D *hEWKPhiStar  = new TH1D("hEWKPhiStar", "",27,PhiStarBins); hEWKPhiStar->Sumw2();
  TH1D *hTopPhiStar  = new TH1D("hTopPhiStar", "",27,PhiStarBins); hTopPhiStar->Sumw2();
  TH1D *hMCPhiStar   = new TH1D("hMCPhiStar",  "",27,PhiStarBins); hMCPhiStar->Sumw2();

  TH1D *hDataZRap = new TH1D("hDataZRap","",24,0,2.4); hDataZRap->Sumw2();
  TH1D *hZmmZRap  = new TH1D("hZmmZRap", "",24,0,2.4); hZmmZRap->Sumw2();
  TH1D *hEWKZRap  = new TH1D("hEWKZRap", "",24,0,2.4); hEWKZRap->Sumw2();
  TH1D *hTopZRap  = new TH1D("hTopZRap", "",24,0,2.4); hTopZRap->Sumw2();
  TH1D *hMCZRap   = new TH1D("hMCZRap",  "",24,0,2.4); hMCZRap->Sumw2();

  TH1D *hDataLep1Pt = new TH1D("hDataLep1Pt","",25,Lep1PtBins); hDataLep1Pt->Sumw2();
  TH1D *hZmmLep1Pt  = new TH1D("hZmmLep1Pt", "",25,Lep1PtBins); hZmmLep1Pt->Sumw2();
  TH1D *hEWKLep1Pt  = new TH1D("hEWKLep1Pt", "",25,Lep1PtBins); hEWKLep1Pt->Sumw2();
  TH1D *hTopLep1Pt  = new TH1D("hTopLep1Pt", "",25,Lep1PtBins); hTopLep1Pt->Sumw2();
  TH1D *hMCLep1Pt   = new TH1D("hMCLep1Pt",  "",25,Lep1PtBins); hMCLep1Pt->Sumw2();

  TH1D *hDataLep2Pt = new TH1D("hDataLep2Pt","",20,Lep2PtBins); hDataLep2Pt->Sumw2();
  TH1D *hZmmLep2Pt  = new TH1D("hZmmLep2Pt", "",20,Lep2PtBins); hZmmLep2Pt->Sumw2();
  TH1D *hEWKLep2Pt  = new TH1D("hEWKLep2Pt", "",20,Lep2PtBins); hEWKLep2Pt->Sumw2();
  TH1D *hTopLep2Pt  = new TH1D("hTopLep2Pt", "",20,Lep2PtBins); hTopLep2Pt->Sumw2();
  TH1D *hMCLep2Pt   = new TH1D("hMCLep2Pt",  "",20,Lep2PtBins); hMCLep2Pt->Sumw2();

  TH1D *hDataLepNegPt = new TH1D("hDataLepNegPt","",25,LepNegPtBins); hDataLepNegPt->Sumw2();
  TH1D *hZmmLepNegPt  = new TH1D("hZmmLepNegPt", "",25,LepNegPtBins); hZmmLepNegPt->Sumw2();
  TH1D *hEWKLepNegPt  = new TH1D("hEWKLepNegPt", "",25,LepNegPtBins); hEWKLepNegPt->Sumw2();
  TH1D *hTopLepNegPt  = new TH1D("hTopLepNegPt", "",25,LepNegPtBins); hTopLepNegPt->Sumw2();
  TH1D *hMCLepNegPt   = new TH1D("hMCLepNegPt",  "",25,LepNegPtBins); hMCLepNegPt->Sumw2();

  TH1D *hDataLepPosPt = new TH1D("hDataLepPosPt","",25,LepPosPtBins); hDataLepPosPt->Sumw2();
  TH1D *hZmmLepPosPt  = new TH1D("hZmmLepPosPt", "",25,LepPosPtBins); hZmmLepPosPt->Sumw2();
  TH1D *hEWKLepPosPt  = new TH1D("hEWKLepPosPt", "",25,LepPosPtBins); hEWKLepPosPt->Sumw2();
  TH1D *hTopLepPosPt  = new TH1D("hTopLepPosPt", "",25,LepPosPtBins); hTopLepPosPt->Sumw2();
  TH1D *hMCLepPosPt   = new TH1D("hMCLepPosPt",  "",25,LepPosPtBins); hMCLepPosPt->Sumw2();

  TH1D *hDataLep1Eta = new TH1D("hDataLep1Eta","",24,0,2.4); hDataLep1Eta->Sumw2();
  TH1D *hZmmLep1Eta  = new TH1D("hZmmLep1Eta", "",24,0,2.4); hZmmLep1Eta->Sumw2();
  TH1D *hEWKLep1Eta  = new TH1D("hEWKLep1Eta", "",24,0,2.4); hEWKLep1Eta->Sumw2();
  TH1D *hTopLep1Eta  = new TH1D("hTopLep1Eta", "",24,0,2.4); hTopLep1Eta->Sumw2();
  TH1D *hMCLep1Eta   = new TH1D("hMCLep1Eta",  "",24,0,2.4); hMCLep1Eta->Sumw2();

  
  TH1D *hDataLep2Eta = new TH1D("hDataLep2Eta","",24,0,2.4); hDataLep2Eta->Sumw2();
  TH1D *hZmmLep2Eta  = new TH1D("hZmmLep2Eta", "",24,0,2.4); hZmmLep2Eta->Sumw2();
  TH1D *hEWKLep2Eta  = new TH1D("hEWKLep2Eta", "",24,0,2.4); hEWKLep2Eta->Sumw2();
  TH1D *hTopLep2Eta  = new TH1D("hTopLep2Eta", "",24,0,2.4); hTopLep2Eta->Sumw2();
  TH1D *hMCLep2Eta   = new TH1D("hMCLep2Eta",  "",24,0,2.4); hMCLep2Eta->Sumw2();

  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  Float_t pfCombIso1, pfCombIso2;

  
  TH2D *h=0;

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
  

  //Setting up rochester corrections
  rochcor2015 *rmcor = new rochcor2015(1234);

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
    intree -> SetBranchStatus("scale1fb",1);
    intree -> SetBranchStatus("q1",1);
    intree -> SetBranchStatus("q2",1);
    intree -> SetBranchStatus("lep1",1);
    intree -> SetBranchStatus("lep2",1);

    intree->SetBranchAddress("runNum",     &runNum);      // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);      // event number
    //intree->SetBranchAddress("matchGen",   &matchGen);    // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category",   &category);    // dilepton category
    //intree->SetBranchAddress("npv",        &npv);	  // number of primary vertices
    //intree->SetBranchAddress("npu",        &npu);	  // number of in-time PU events (MC)
    //intree->SetBranchAddress("genVPt",     &genVPt);      // GEN Z boson pT (signal MC)
    //intree->SetBranchAddress("genVPhi",    &genVPhi);     // GEN Z boson phi (signal MC)
    //intree->SetBranchAddress("genVy",      &genVy);       // GEN Z boson rapidity (signal MC)
    //intree->SetBranchAddress("genVMass",   &genVMass);    // GEN Z boson mass (signal MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    //intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);    // event weight per 1/fb (MC)
    //intree->SetBranchAddress("scale1fbDown",   &scale1fbDown);    // event weight per 1/fb (MC)
    //intree->SetBranchAddress("met",        &met);	  // MET
    //intree->SetBranchAddress("metPhi",     &metPhi);      // phi(MET)
    //intree->SetBranchAddress("sumEt",      &sumEt);       // Sum ET
    //intree->SetBranchAddress("u1",         &u1);	  // parallel component of recoil
    //intree->SetBranchAddress("u2",         &u2);	  // perpendicular component of recoil
    intree->SetBranchAddress("q1",         &q1);	  // charge of tag lepton
    intree->SetBranchAddress("q2",         &q2);	  // charge of probe lepton
    //intree->SetBranchAddress("dilep",      &dilep);       // dilepton 4-vector
    intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
    //intree->SetBranchAddress("pfCombIso1", &pfCombIso1);  // combined PF isolation of tag lepton
    //intree->SetBranchAddress("pfCombIso2", &pfCombIso2);  // combined PF isolation of probe lepton
   
    //
    // loop over events
    //
    cout<<"Beginning of Loop: nentries"<<intree->GetEntries() << endl;
    sw_[0].Start();

    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      sw_[6].Start(0);
      intree->GetEntry(ientry);
      sw_[6].Stop();
      if ( (ientry & 16383) == 0 )
	{
	  sw_[0].Stop();
	  cout<< "Took " <<  sw_[0] .RealTime() << "s CPU " << sw_[0] . CpuTime()<< endl;
	  sw_[0].Reset();
	  sw_[0].Start();
	  cout<< "Took 1 " <<  sw_[1] .RealTime() << "s CPU " << sw_[1] . CpuTime()<< endl;
	  sw_[1].Reset();
	  cout<< "Took 2 " <<  sw_[2] .RealTime() << "s CPU " << sw_[2] . CpuTime()<< endl;
	  sw_[2].Reset();
	  cout<< "Took 3 " <<  sw_[3] .RealTime() << "s CPU " << sw_[3] . CpuTime()<< endl;
	  sw_[3].Reset();
	  cout<< "Took 4 " <<  sw_[4] .RealTime() << "s CPU " << sw_[4] . CpuTime()<< endl;
	  sw_[4].Reset();
	  cout<< "Took 5 " <<  sw_[5] .RealTime() << "s CPU " << sw_[5] . CpuTime()<< endl;
	  sw_[5].Reset();
	  cout<< "Took 6 " <<  sw_[6] .RealTime() << "s CPU " << sw_[6] . CpuTime()<< endl;
	  sw_[6].Reset();
	  cout<<endl;
	}

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
 	     
      // fill Z events passing selection (MuMu2HLT + MuMu1HLT)
      if((category==eMuMu2HLT) || (category==eMuMu1HLT) || (category==eMuMu1HLT1L1)) {
        if(typev[ifile]==eData) { 
          sw_[5].Start(0); 
	  TLorentzVector mu1;
	  TLorentzVector mu2;
	  mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
	  mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
	  float qter1=1.0;
	  float qter2=1.0;
          sw_[5].Stop(); 

	  sw_[1].Start(0);
	  rmcor->momcor_data(mu1,q1,0,qter1);
	  rmcor->momcor_data(mu2,q2,0,qter2);
	  sw_[1].Stop();
	  sw_[2].Start(0);

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

	  sw_[2].Stop();
	  sw_[3].Start(0);

	  hData->Fill(mass); 
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
	  sw_[3].Stop();
	} else {
	  sw_[4].Start(0);
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
	     
          effdata=1; effmc=1;
	  if(lq1>0) { 
	    effdata *= dataSelEff_pos.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmSelEff_pos.getEff((l1.Eta()), l1.Pt()); 
	  } else {
	    effdata *= dataSelEff_neg.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmSelEff_neg.getEff((l1.Eta()), l1.Pt()); 
	  }
	  if(lq2>0) {
	    effdata *= dataSelEff_pos.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmSelEff_pos.getEff((l2.Eta()), l2.Pt());
	  } else {
	    effdata *= dataSelEff_neg.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmSelEff_neg.getEff((l2.Eta()), l2.Pt());
	  }
	  corr *= effdata/effmc;

	  effdata=1; effmc=1;
	  if(lq1>0) { 
	    effdata *= dataStaEff_pos.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmStaEff_pos.getEff((l1.Eta()), l1.Pt()); 
	  } else {
	    effdata *= dataStaEff_neg.getEff((l1.Eta()), l1.Pt()); 
	    effmc   *= zmmStaEff_neg.getEff((l1.Eta()), l1.Pt()); 
	  }
	  if(lq2>0) {
	    effdata *= dataStaEff_pos.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmStaEff_pos.getEff((l2.Eta()), l2.Pt());
	  } else {
	    effdata *= dataStaEff_neg.getEff((l2.Eta()), l2.Pt()); 
	    effmc   *= zmmStaEff_neg.getEff((l2.Eta()), l2.Pt());
	  }
	  corr *= effdata/effmc; 
	  
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
	
	  mass = (l1+l2).M();
	  pt = (l1+l2).Pt();
	  rapidity = (l1+l2).Rapidity();

	  phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
	  if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
	  else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
	  phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));

	  if(typev[ifile]==eZmm) 
	    {
	      hZmm->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
	      hZmmZPt->Fill(pt,weight*corr); 
	      hMCZPt->Fill(pt,weight*corr);
	      hZmmPhiStar->Fill(phistar,weight*corr); 
	      hMCPhiStar->Fill(phistar,weight*corr);
	      hZmmZRap->Fill(fabs(rapidity),weight*corr); 
	      hMCZRap->Fill(fabs(rapidity),weight*corr);
	      hZmmLep1Pt->Fill(l1.Pt(),weight*corr); 
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);
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
	      hZmmLep2Pt->Fill(l2.Pt(),weight*corr); 
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);
	      hZmmLep1Eta->Fill(fabs(l1.Eta()),weight*corr); 
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
	      hZmmLep2Eta->Fill(fabs(l2.Eta()),weight*corr); 
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
	  if(typev[ifile]==eEWK) 
	    {
	      hEWK->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);

	      hEWKZPt->Fill(pt,weight*corr); 
	      hMCZPt->Fill(pt,weight*corr);

	      hEWKPhiStar->Fill(phistar,weight*corr);
	      hMCPhiStar->Fill(phistar,weight*corr);

	      hEWKZRap->Fill(fabs(rapidity),weight*corr); 
	      hMCZRap->Fill(fabs(rapidity),weight*corr);

	      hEWKLep1Pt->Fill(l1.Pt(),weight*corr);
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);

	      hEWKLep2Pt->Fill(l2.Pt(),weight*corr);
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);

	      if(lq1<0)
		{
		  hEWKLepNegPt->Fill(l1.Pt(),weight*corr);
		  hMCLepNegPt->Fill(l1.Pt(),weight*corr);
		  
		  hEWKLepPosPt->Fill(l2.Pt(),weight*corr);
		  hMCLepPosPt->Fill(l2.Pt(),weight*corr);
		}
	      else
		{
		  hEWKLepNegPt->Fill(l2.Pt(),weight*corr);
		  hMCLepNegPt->Fill(l2.Pt(),weight*corr);
		  
		  hEWKLepPosPt->Fill(l1.Pt(),weight*corr);
		  hMCLepPosPt->Fill(l1.Pt(),weight*corr);
		}

	      hEWKLep1Eta->Fill(fabs(l1.Eta()),weight*corr); 
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);

	      hEWKLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
	  if(typev[ifile]==eTop) 
	    {
	      hTop->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);

	      hTopZPt->Fill(pt,weight*corr); 
	      hMCZPt->Fill(pt,weight*corr);

	      hTopPhiStar->Fill(phistar,weight*corr);
	      hMCPhiStar->Fill(phistar,weight*corr);

	      hTopZRap->Fill(fabs(rapidity),weight*corr); 
	      hMCZRap->Fill(fabs(rapidity),weight*corr);

	      hTopLep1Pt->Fill(l1.Pt(),weight*corr);
	      hMCLep1Pt->Fill(l1.Pt(),weight*corr);

	      hTopLep2Pt->Fill(l2.Pt(),weight*corr);
	      hMCLep2Pt->Fill(l2.Pt(),weight*corr);

	      if(lq1<0)
		{
		  hTopLepNegPt->Fill(l1.Pt(),weight*corr);
		  hMCLepNegPt->Fill(l1.Pt(),weight*corr);
		  
		  hTopLepPosPt->Fill(l2.Pt(),weight*corr);
		  hMCLepPosPt->Fill(l2.Pt(),weight*corr);
		}
	      else
		{
		  hTopLepNegPt->Fill(l2.Pt(),weight*corr);
		  hMCLepNegPt->Fill(l2.Pt(),weight*corr);
		  
		  hTopLepPosPt->Fill(l1.Pt(),weight*corr);
		  hMCLepPosPt->Fill(l1.Pt(),weight*corr);
		}

	      hTopLep1Eta->Fill(fabs(l1.Eta()),weight*corr);
	      hMCLep1Eta->Fill(fabs(l1.Eta()),weight*corr);

	      hTopLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	      hMCLep2Eta->Fill(fabs(l2.Eta()),weight*corr);
	    }
	  sw_[4].Stop();
	}//end MC
      } //end category
    }//end loop ientry
    
    delete infile;
    infile=0, intree=0;
  } // end loop ifile

  

  outFile->cd();

  hDataZPt->Write();
  hEWKZPt->Write();
  hTopZPt->Write();

  hDataPhiStar->Write();
  hEWKPhiStar->Write();
  hTopPhiStar->Write();
 
  hDataZRap->Write();
  hEWKZRap->Write();
  hTopZRap->Write();
 
  hDataLep1Pt->Write();
  hEWKLep1Pt->Write();
  hTopLep1Pt->Write();
  
  hDataLep2Pt->Write();
  hEWKLep2Pt->Write();
  hTopLep2Pt->Write();

  hDataLepNegPt->Write();
  hEWKLepNegPt->Write();
  hTopLepNegPt->Write();

  hDataLepPosPt->Write();
  hEWKLepPosPt->Write();
  hTopLepPosPt->Write();
 
  hDataLep1Eta->Write();
  hEWKLep1Eta->Write();
  hTopLep1Eta->Write();
  
  hDataLep2Eta->Write();
  hEWKLep2Eta->Write();
  hTopLep2Eta->Write();
  
  outFile->Write();
  outFile->Close(); 

  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     

  gBenchmark->Show("plotZmmResScaleUncert");
}


