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
void create( vector<TH1D> &v , const string name, int nbins, double xlow, double xhigh , const int ntoys) {for(int itoys=0;itoys!=ntoys;++itoys) {v.push_back(TH1D((name+Form("_%d",itoys)).c_str(),"",nbins,xlow,xhigh)); v[itoys].Sumw2();}}
void create( vector<TH1D> &v , const string name, int nbins, double* x , const int ntoys) {for(int itoys=0;itoys!=ntoys;++itoys) {v.push_back(TH1D((name+Form("_%d",itoys)).c_str(),"",nbins,x)); v[itoys].Sumw2();}}

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

  //----- 
  const int NTOYS = 100;

  // efficiency files
  const TString dataHLTEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/NewBin/MG/eff.root";
  const TString dataHLTEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/NewBin/MG/eff.root";
  const TString zmmHLTEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/NewBin/CT/eff.root";
  const TString zmmHLTEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuHLTEff/NewBin/CT/eff.root";

  const TString dataSelEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/MG/eff.root";
  const TString dataSelEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/MG/eff.root";
  const TString zmmSelEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/CT/eff.root";
  const TString zmmSelEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/CT/eff.root";

  const TString dataTrkEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/MG/eff.root";
  const TString dataTrkEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/MG/eff.root";
  const TString zmmTrkEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/CT/eff.root";
  const TString zmmTrkEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuSITEff/NewBin/CT/eff.root";

  const TString dataStaEffName_pos = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/NewBin/MG/eff.root";
  const TString dataStaEffName_neg = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/NewBin/MG/eff.root";
  const TString zmmStaEffName_pos  = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/NewBin/CT/eff.root";
  const TString zmmStaEffName_neg  = "/data/blue/xniu/WZXSection/NewMu/MuStaEff/NewBin/CT/eff.root";

   
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
  double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
  double PhiStarBins[]={0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.01,0.012,0.014,0.016,0.018,0.021,0.024,0.027,0.030,0.034,0.038,0.044,0.050,0.058,0.066,0.076,0.088,0.10,0.12,0.14,0.16,0.18,0.20,0.24,0.28,0.34,0.42,0.52,0.64,0.8,1.0,1.5,2,3};
  double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double Lep2PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
  double LepNegPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double LepPosPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  
  vector<TH1D> hData; create(hData,"hData",NBINS,MASS_LOW,MASS_HIGH,NTOYS);
  vector<TH1D> hZmm; create(hZmm,"hZmm",NBINS,MASS_LOW,MASS_HIGH,NTOYS);
  vector<TH1D> hEWK; create(hEWK,"hEWK",NBINS,MASS_LOW,MASS_HIGH,NTOYS);
  vector<TH1D> hTop; create(hTop,"hTop",NBINS,MASS_LOW,MASS_HIGH,NTOYS);
  vector<TH1D> hMC; create(hMC,"hMC",NBINS,MASS_LOW,MASS_HIGH,NTOYS);

  const int nBinsZPt= sizeof(ZPtBins)/sizeof(double)-1;
  vector<TH1D> hDataZPt ; create(hDataZPt,"hDataZPt",nBinsZPt,ZPtBins,NTOYS);
  vector<TH1D> hZmmZPt ; create(hZmmZPt,"hZmmZPt",nBinsZPt,ZPtBins,NTOYS);
  vector<TH1D> hEWKZPt ; create(hEWKZPt,"hEWKZPt",nBinsZPt,ZPtBins,NTOYS);
  vector<TH1D> hTopZPt ; create(hTopZPt,"hTopZPt",nBinsZPt,ZPtBins,NTOYS);
  vector<TH1D> hMCZPt ; create(hMCZPt,"hMCZPt",nBinsZPt,ZPtBins,NTOYS);

  const int nBinsPhiStar= sizeof(PhiStarBins)/sizeof(double)-1;
  vector<TH1D> hDataPhiStar ; create(hDataPhiStar,"hDataPhiStar",nBinsPhiStar,PhiStarBins,NTOYS);
  vector<TH1D> hZmmPhiStar ; create(hZmmPhiStar,"hZmmPhiStar",nBinsPhiStar,PhiStarBins,NTOYS);
  vector<TH1D> hEWKPhiStar ; create(hEWKPhiStar,"hEWKPhiStar",nBinsPhiStar,PhiStarBins,NTOYS);
  vector<TH1D> hTopPhiStar ; create(hTopPhiStar,"hTopPhiStar",nBinsPhiStar,PhiStarBins,NTOYS);
  vector<TH1D> hMCPhiStar ; create(hMCPhiStar,"hMCPhiStar",nBinsPhiStar,PhiStarBins,NTOYS);

  vector<TH1D> hDataZRap ; create(hDataZRap,"hDataZRap",24,0,2.4,NTOYS);
  vector<TH1D> hZmmZRap ; create(hZmmZRap,"hZmmZRap",24,0,2.4,NTOYS);
  vector<TH1D> hEWKZRap ; create(hEWKZRap,"hEWKZRap",24,0,2.4,NTOYS);
  vector<TH1D> hTopZRap ; create(hTopZRap,"hTopZRap",24,0,2.4,NTOYS);
  vector<TH1D> hMCZRap ; create(hMCZRap,"hMCZRap",24,0,2.4,NTOYS);

  const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
  vector<TH1D> hDataLep1Pt ; create(hDataLep1Pt,"hDataLep1Pt",nBinsLep1Pt,Lep1PtBins,NTOYS);
  vector<TH1D> hZmmLep1Pt ; create(hZmmLep1Pt,"hZmmLep1Pt",nBinsLep1Pt,Lep1PtBins,NTOYS);
  vector<TH1D> hEWKLep1Pt ; create(hEWKLep1Pt,"hEWKLep1Pt",nBinsLep1Pt,Lep1PtBins,NTOYS);
  vector<TH1D> hTopLep1Pt ; create(hTopLep1Pt,"hTopLep1Pt",nBinsLep1Pt,Lep1PtBins,NTOYS);
  vector<TH1D> hMCLep1Pt ; create(hMCLep1Pt,"hMCLep1Pt",nBinsLep1Pt,Lep1PtBins,NTOYS);

  const int nBinsLep2Pt= sizeof(Lep2PtBins)/sizeof(double)-1;
  vector<TH1D> hDataLep2Pt ; create(hDataLep2Pt,"hDataLep2Pt",nBinsLep2Pt,Lep2PtBins,NTOYS);
  vector<TH1D> hZmmLep2Pt ; create(hZmmLep2Pt,"hZmmLep2Pt",nBinsLep2Pt,Lep2PtBins,NTOYS);
  vector<TH1D> hEWKLep2Pt ; create(hEWKLep2Pt,"hEWKLep2Pt",nBinsLep2Pt,Lep2PtBins,NTOYS);
  vector<TH1D> hTopLep2Pt ; create(hTopLep2Pt,"hTopLep2Pt",nBinsLep2Pt,Lep2PtBins,NTOYS);
  vector<TH1D> hMCLep2Pt ; create(hMCLep2Pt,"hMCLep2Pt",nBinsLep2Pt,Lep2PtBins,NTOYS);

  const int nBinsLepNegPt= sizeof(LepNegPtBins)/sizeof(double)-1;
  vector<TH1D> hDataLepNegPt ; create(hDataLepNegPt,"hDataLepNegPt",nBinsLepNegPt,LepNegPtBins,NTOYS);
  vector<TH1D> hZmmLepNegPt ; create(hZmmLepNegPt,"hZmmLepNegPt",nBinsLepNegPt,LepNegPtBins,NTOYS);
  vector<TH1D> hEWKLepNegPt ; create(hEWKLepNegPt,"hEWKLepNegPt",nBinsLepNegPt,LepNegPtBins,NTOYS);
  vector<TH1D> hTopLepNegPt ; create(hTopLepNegPt,"hTopLepNegPt",nBinsLepNegPt,LepNegPtBins,NTOYS);
  vector<TH1D> hMCLepNegPt ; create(hMCLepNegPt,"hMCLepNegPt",nBinsLepNegPt,LepNegPtBins,NTOYS);

  const int nBinsLepPosPt= sizeof(LepPosPtBins)/sizeof(double)-1;
  vector<TH1D> hDataLepPosPt ; create(hDataLepPosPt,"hDataLepPosPt",nBinsLepPosPt,LepPosPtBins,NTOYS);
  vector<TH1D> hZmmLepPosPt ; create(hZmmLepPosPt,"hZmmLepPosPt",nBinsLepPosPt,LepPosPtBins,NTOYS);
  vector<TH1D> hEWKLepPosPt ; create(hEWKLepPosPt,"hEWKLepPosPt",nBinsLepPosPt,LepPosPtBins,NTOYS);
  vector<TH1D> hTopLepPosPt ; create(hTopLepPosPt,"hTopLepPosPt",nBinsLepPosPt,LepPosPtBins,NTOYS);
  vector<TH1D> hMCLepPosPt ; create(hMCLepPosPt,"hMCLepPosPt",nBinsLepPosPt,LepPosPtBins,NTOYS);

  vector<TH1D> hDataLep1Eta ; create(hDataLep1Eta,"hDataLep1Eta",24,0,2.4,NTOYS);
  vector<TH1D> hZmmLep1Eta ; create(hZmmLep1Eta,"hZmmLep1Eta",24,0,2.4,NTOYS);
  vector<TH1D> hEWKLep1Eta ; create(hEWKLep1Eta,"hEWKLep1Eta",24,0,2.4,NTOYS);
  vector<TH1D> hTopLep1Eta ; create(hTopLep1Eta,"hTopLep1Eta",24,0,2.4,NTOYS);
  vector<TH1D> hMCLep1Eta ; create(hMCLep1Eta,"hMCLep1Eta",24,0,2.4,NTOYS);

  
  vector<TH1D> hDataLep2Eta ; create(hDataLep2Eta,"hDataLep2Eta",24,0,2.4,NTOYS);
  vector<TH1D> hZmmLep2Eta ; create(hZmmLep2Eta,"hZmmLep2Eta",24,0,2.4,NTOYS);
  vector<TH1D> hEWKLep2Eta ; create(hEWKLep2Eta,"hEWKLep2Eta",24,0,2.4,NTOYS);
  vector<TH1D> hTopLep2Eta ; create(hTopLep2Eta,"hTopLep2Eta",24,0,2.4,NTOYS);
  vector<TH1D> hMCLep2Eta ; create(hMCLep2Eta,"hMCLep2Eta",24,0,2.4,NTOYS);

  
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Int_t   q1, q2;
  TLorentzVector *lep1=0, *lep2=0;
  
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
  vector<rochcor2015> vRocToys;
  for (int i=0 ; i<NTOYS; ++i) vRocToys.push_back(rochcor2015(1234+i*1000));
   

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
    intree->SetBranchAddress("category",   &category);    // dilepton category
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
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
 	     
      for(int itoys=0;itoys!=NTOYS;++itoys)
	{
      // fill Z events passing selection (MuMu2HLT + MuMu1HLT)
      if((category==eMuMu2HLT) || (category==eMuMu1HLT) || (category==eMuMu1HLT1L1)) {
        if(typev[ifile]==eData) { 
	  TLorentzVector mu1;
	  TLorentzVector mu2;
	  mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
	  mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
	  float qter1=1.0;
	  float qter2=1.0;
          
	  vRocToys[itoys].momcor_data(mu1,q1,0,qter1);
	  vRocToys[itoys].momcor_data(mu2,q2,0,qter2);
	  
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

	  hData[itoys].Fill(mass); 
	  hDataZPt[itoys].Fill(pt); 
	  hDataPhiStar[itoys].Fill(phistar); 
	  hDataLep1Pt[itoys].Fill(l1.Pt()); 
	  hDataLep2Pt[itoys].Fill(l2.Pt()); 
	  if(lq1<0)
	    {
	      hDataLepNegPt[itoys].Fill(l1.Pt()); 
	      hDataLepPosPt[itoys].Fill(l2.Pt());
	    }
	  else 
	    {
	      hDataLepNegPt[itoys].Fill(l2.Pt()); 
	      hDataLepPosPt[itoys].Fill(l1.Pt());
	    }
	  hDataLep1Eta[itoys].Fill(fabs(l1.Eta())); 
	  hDataLep2Eta[itoys].Fill(fabs(l2.Eta())); 
	  hDataZRap[itoys].Fill(fabs(rapidity));
	} else {
	  TLorentzVector mu1;
	  TLorentzVector mu2;
	  mu1.SetPtEtaPhiM(lep1->Pt(),lep1->Eta(),lep1->Phi(),mu_MASS);
	  mu2.SetPtEtaPhiM(lep2->Pt(),lep2->Eta(),lep2->Phi(),mu_MASS);
	  float qter1=1.0;
	  float qter2=1.0;

	  vRocToys[itoys].momcor_mc(mu1,q1,0,qter1);
	  vRocToys[itoys].momcor_mc(mu2,q2,0,qter2);

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
	      hZmm[itoys].Fill(mass,weight*corr); 
	      hMC[itoys].Fill(mass,weight*corr);
	      hZmmZPt[itoys].Fill(pt,weight*corr); 
	      hMCZPt[itoys].Fill(pt,weight*corr);
	      hZmmPhiStar[itoys].Fill(phistar,weight*corr); 
	      hMCPhiStar[itoys].Fill(phistar,weight*corr);
	      hZmmZRap[itoys].Fill(fabs(rapidity),weight*corr); 
	      hMCZRap[itoys].Fill(fabs(rapidity),weight*corr);
	      hZmmLep1Pt[itoys].Fill(l1.Pt(),weight*corr); 
	      hMCLep1Pt[itoys].Fill(l1.Pt(),weight*corr);
	      if(lq1<0)
		{
		  hZmmLepNegPt[itoys].Fill(l1.Pt(),weight*corr); 
		  hMCLepNegPt[itoys].Fill(l1.Pt(),weight*corr);
		  hZmmLepPosPt[itoys].Fill(l2.Pt(),weight*corr); 
		  hMCLepPosPt[itoys].Fill(l2.Pt(),weight*corr);
		}
	      else 
		{
		  hZmmLepNegPt[itoys].Fill(l2.Pt(),weight*corr); 
		  hMCLepNegPt[itoys].Fill(l2.Pt(),weight*corr);
		  hZmmLepPosPt[itoys].Fill(l1.Pt(),weight*corr); 
		  hMCLepPosPt[itoys].Fill(l1.Pt(),weight*corr);
		}
	      hZmmLep2Pt[itoys].Fill(l2.Pt(),weight*corr); 
	      hMCLep2Pt[itoys].Fill(l2.Pt(),weight*corr);
	      hZmmLep1Eta[itoys].Fill(fabs(l1.Eta()),weight*corr); 
	      hMCLep1Eta[itoys].Fill(fabs(l1.Eta()),weight*corr);
	      hZmmLep2Eta[itoys].Fill(fabs(l2.Eta()),weight*corr); 
	      hMCLep2Eta[itoys].Fill(fabs(l2.Eta()),weight*corr);
	    }
	  if(typev[ifile]==eEWK) 
	    {
	      hEWK[itoys].Fill(mass,weight*corr); 
	      hMC[itoys].Fill(mass,weight*corr);

	      hEWKZPt[itoys].Fill(pt,weight*corr); 
	      hMCZPt[itoys].Fill(pt,weight*corr);

	      hEWKPhiStar[itoys].Fill(phistar,weight*corr);
	      hMCPhiStar[itoys].Fill(phistar,weight*corr);

	      hEWKZRap[itoys].Fill(fabs(rapidity),weight*corr); 
	      hMCZRap[itoys].Fill(fabs(rapidity),weight*corr);

	      hEWKLep1Pt[itoys].Fill(l1.Pt(),weight*corr);
	      hMCLep1Pt[itoys].Fill(l1.Pt(),weight*corr);

	      hEWKLep2Pt[itoys].Fill(l2.Pt(),weight*corr);
	      hMCLep2Pt[itoys].Fill(l2.Pt(),weight*corr);

	      if(lq1<0)
		{
		  hEWKLepNegPt[itoys].Fill(l1.Pt(),weight*corr);
		  hMCLepNegPt[itoys].Fill(l1.Pt(),weight*corr);
		  
		  hEWKLepPosPt[itoys].Fill(l2.Pt(),weight*corr);
		  hMCLepPosPt[itoys].Fill(l2.Pt(),weight*corr);
		}
	      else
		{
		  hEWKLepNegPt[itoys].Fill(l2.Pt(),weight*corr);
		  hMCLepNegPt[itoys].Fill(l2.Pt(),weight*corr);
		  
		  hEWKLepPosPt[itoys].Fill(l1.Pt(),weight*corr);
		  hMCLepPosPt[itoys].Fill(l1.Pt(),weight*corr);
		}

	      hEWKLep1Eta[itoys].Fill(fabs(l1.Eta()),weight*corr); 
	      hMCLep1Eta[itoys].Fill(fabs(l1.Eta()),weight*corr);

	      hEWKLep2Eta[itoys].Fill(fabs(l2.Eta()),weight*corr);
	      hMCLep2Eta[itoys].Fill(fabs(l2.Eta()),weight*corr);
	    }
	  if(typev[ifile]==eTop) 
	    {
	      hTop[itoys].Fill(mass,weight*corr); 
	      hMC[itoys].Fill(mass,weight*corr);

	      hTopZPt[itoys].Fill(pt,weight*corr); 
	      hMCZPt[itoys].Fill(pt,weight*corr);

	      hTopPhiStar[itoys].Fill(phistar,weight*corr);
	      hMCPhiStar[itoys].Fill(phistar,weight*corr);

	      hTopZRap[itoys].Fill(fabs(rapidity),weight*corr); 
	      hMCZRap[itoys].Fill(fabs(rapidity),weight*corr);

	      hTopLep1Pt[itoys].Fill(l1.Pt(),weight*corr);
	      hMCLep1Pt[itoys].Fill(l1.Pt(),weight*corr);

	      hTopLep2Pt[itoys].Fill(l2.Pt(),weight*corr);
	      hMCLep2Pt[itoys].Fill(l2.Pt(),weight*corr);

	      if(lq1<0)
		{
		  hTopLepNegPt[itoys].Fill(l1.Pt(),weight*corr);
		  hMCLepNegPt[itoys].Fill(l1.Pt(),weight*corr);
		  
		  hTopLepPosPt[itoys].Fill(l2.Pt(),weight*corr);
		  hMCLepPosPt[itoys].Fill(l2.Pt(),weight*corr);
		}
	      else
		{
		  hTopLepNegPt[itoys].Fill(l2.Pt(),weight*corr);
		  hMCLepNegPt[itoys].Fill(l2.Pt(),weight*corr);
		  
		  hTopLepPosPt[itoys].Fill(l1.Pt(),weight*corr);
		  hMCLepPosPt[itoys].Fill(l1.Pt(),weight*corr);
		}

	      hTopLep1Eta[itoys].Fill(fabs(l1.Eta()),weight*corr);
	      hMCLep1Eta[itoys].Fill(fabs(l1.Eta()),weight*corr);

	      hTopLep2Eta[itoys].Fill(fabs(l2.Eta()),weight*corr);
	      hMCLep2Eta[itoys].Fill(fabs(l2.Eta()),weight*corr);
	    }
	}//end MC
      } //end category
	}//end toys
    }//end loop ientry
    
    delete infile;
    infile=0, intree=0;
  } // end loop ifile

  

  outFile->cd();

  for(int itoys=0;itoys!=NTOYS;++itoys)
    {
  hDataZPt[itoys].Write();
  hEWKZPt[itoys].Write();
  hTopZPt[itoys].Write();

  hDataPhiStar[itoys].Write();
  hEWKPhiStar[itoys].Write();
  hTopPhiStar[itoys].Write();
 
  hDataZRap[itoys].Write();
  hEWKZRap[itoys].Write();
  hTopZRap[itoys].Write();
 
  hDataLep1Pt[itoys].Write();
  hEWKLep1Pt[itoys].Write();
  hTopLep1Pt[itoys].Write();
  
  hDataLep2Pt[itoys].Write();
  hEWKLep2Pt[itoys].Write();
  hTopLep2Pt[itoys].Write();

  hDataLepNegPt[itoys].Write();
  hEWKLepNegPt[itoys].Write();
  hTopLepNegPt[itoys].Write();

  hDataLepPosPt[itoys].Write();
  hEWKLepPosPt[itoys].Write();
  hTopLepPosPt[itoys].Write();
 
  hDataLep1Eta[itoys].Write();
  hEWKLep1Eta[itoys].Write();
  hTopLep1Eta[itoys].Write();
  
  hDataLep2Eta[itoys].Write();
  hEWKLep2Eta[itoys].Write();
  hTopLep2Eta[itoys].Write();
    }
  
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


