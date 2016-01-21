//================================================================================================
//
// Make plots of various distributions after Zmm selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TLatex.h>
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include "TLorentzVector.h"               // 4-vector class

#include "ConfParse.hh"                   // input conf file parser
#include "../Utils/CSample.hh"            // helper class to handle samples
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/LeptonCorr.hh"

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

//helper class to handle rochester corrections
#include <rochcor2015.h>
#include <muresolution_run2.h>

#include "TStopwatch.h" //PROFILE

#endif


TStopwatch sw_[10]; //PROFILE 


//=== MAIN MACRO ================================================================================================= 

void plotZmmGen(const TString  conf,            // input file
             const TString  inputDir,        // input directory
	     const TString  outputDir,       // output directory
	     const Double_t lumi             // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmmGen");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  const Double_t mu_MASS  = 0.1057;
  
  const TString format("png");
  
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
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
 
  vector<TString>  snamev;      // sample name 
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  
  // setup efficiency shape systematics
  TFile *StaSigSysFile = new TFile(StaEffSignalShapeSys);
  TH2D *hStaSigSys = (TH2D*)StaSigSysFile->Get("h");
  TFile *StaBkgSysFile = new TFile(StaEffBackgroundShapeSys);
  TH2D *hStaBkgSys = (TH2D*)StaBkgSysFile->Get("h");
  TFile *SelSigSysFile = new TFile(SelEffSignalShapeSys);
  TH2D *hSelSigSys = (TH2D*)SelSigSysFile->Get("h");
  TFile *SelBkgSysFile = new TFile(SelEffBackgroundShapeSys);
  TH2D *hSelBkgSys = (TH2D*)SelBkgSysFile->Get("h");
  

  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  npv, npu;
  UInt_t  triggerDec;
  UInt_t  goodPV;
  UInt_t  matchTrigger;
  UInt_t  ngenlep;
  TLorentzVector *genlep1=0, *genlep2=0;
  Int_t   genq1, genq2;
  UInt_t nlep;
  TLorentzVector *lep1=0, *lep2=0;
  Int_t   q1, q2;
  Float_t scale1fbGen,scale1fb, scale1fbUp, scale1fbDown;

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

  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Read input file and get the TTrees
    TString infilename=inputDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    cout << "Processing " << infilename << "..." << endl;
    infile = TFile::Open(infilename);assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);
    
    intree->SetBranchAddress("runNum",   &runNum);     // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);     // event number
    intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("npv",      &npv);        // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);        // number of in-time PU events (MC)
    intree->SetBranchAddress("triggerDec",   &triggerDec);    // event pass the trigger
    intree->SetBranchAddress("goodPV",   &goodPV);    // event has a good PV
    intree->SetBranchAddress("matchTrigger",   &matchTrigger);    // event has at least one lepton matched to the trigger
    intree->SetBranchAddress("ngenlep",     &ngenlep);      // number of gen leptons
    intree->SetBranchAddress("genlep1",   &genlep1);     // gen lepton1 4-vector
    intree->SetBranchAddress("genlep2",   &genlep2);     // gen lepton2 4-vector
    intree->SetBranchAddress("genq1",     &genq1);     // charge gen lepton1
    intree->SetBranchAddress("genq2",     &genq2);     // charge gen lepton2
    intree->SetBranchAddress("nlep",     &nlep);      // number of leptons
    intree->SetBranchAddress("lep1",       &lep1);     // lepton1 4-vector
    intree->SetBranchAddress("lep2",       &lep2);     // lepton2 4-vector
    intree->SetBranchAddress("q1",       &q1);     // charge lepton1
    intree->SetBranchAddress("q2",       &q2);     // charge lepton2
    intree->SetBranchAddress("scale1fbGen",   &scale1fbGen);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",   &scale1fbDown);    // event weight per 1/fb (MC)
    
    //
    // Set up output file
    //
    TString outfilename = outputDir + TString("/") + snamev[isam] + TString("_UnfoldInputs.root");
    TFile *outFile = new TFile(outfilename,"RECREATE");

    const TString plotDir = outputDir + TString("/Plots_") + snamev[isam];
    gSystem->mkdir(plotDir,kTRUE);
    CPlot::sOutDir = plotDir;
    
    //
    // Create histograms
    //
    double ZPtBins[35]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,375,500,1000};
    
    double PhiStarBins[28]={0,0.01,0.012,0.014,0.017,0.021,0.025,0.030,0.036,0.043,0.052,0.062,0.074,0.089,0.11,0.13,0.15,0.18,0.22,0.27,0.32,0.38,0.46,0.55,0.66,0.79,0.95,1.1};
    
    double Lep1PtBins[26]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,154,171,190,211,234,265,300};
    double Lep2PtBins[21]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,157,200};
    
    double LepNegPtBins[26]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,154,171,190,211,234,265,300};
    double LepPosPtBins[26]={25,27.5,30.3,33.3,36.6,40.3,44.3,48.7,53.6,58.9,64.8,71.3,78.5,86.3,94.9,104,115,126,139,154,171,190,211,234,265,300};
     
    TH1D *hMassMC  = new TH1D("hMassMC","",30,60,120); hMassMC->Sumw2();
    TH1D *hZPtReco  = new TH1D("hZPtReco","",34,ZPtBins); hZPtReco->Sumw2();
    TH1D *hZPtTruth  = new TH1D("hZPtTruth","",34,ZPtBins); hZPtTruth->Sumw2();
    TH2D *hZPtMatrix  = new TH2D("hZPtMatrix","",34,ZPtBins,34,ZPtBins); hZPtMatrix->Sumw2();
    TH1D *hZPtReco_EffBin  = new TH1D("hZPtReco_EffBin","",34,ZPtBins); hZPtReco_EffBin->Sumw2();
    TH2D *hZPtMatrix_EffBin  = new TH2D("hZPtMatrix_EffBin","",34,ZPtBins,34,ZPtBins); hZPtMatrix_EffBin->Sumw2();
    TH1D *hZPtReco_EffStatUp  = new TH1D("hZPtReco_EffStatUp","",34,ZPtBins); hZPtReco_EffStatUp->Sumw2();
    TH2D *hZPtMatrix_EffStatUp  = new TH2D("hZPtMatrix_EffStatUp","",34,ZPtBins,34,ZPtBins); hZPtMatrix_EffStatUp->Sumw2();
    TH1D *hZPtReco_EffStatDown  = new TH1D("hZPtReco_EffStatDown","",34,ZPtBins); hZPtReco_EffStatDown->Sumw2();
    TH2D *hZPtMatrix_EffStatDown  = new TH2D("hZPtMatrix_EffStatDown","",34,ZPtBins,34,ZPtBins); hZPtMatrix_EffStatDown->Sumw2();
    TH1D *hZPtReco_EffSigShape  = new TH1D("hZPtReco_EffSigShape","",34,ZPtBins); hZPtReco_EffSigShape->Sumw2();
    TH2D *hZPtMatrix_EffSigShape  = new TH2D("hZPtMatrix_EffSigShape","",34,ZPtBins,34,ZPtBins); hZPtMatrix_EffSigShape->Sumw2();
    TH1D *hZPtReco_EffBkgShape  = new TH1D("hZPtReco_EffBkgShape","",34,ZPtBins); hZPtReco_EffBkgShape->Sumw2();
    TH2D *hZPtMatrix_EffBkgShape  = new TH2D("hZPtMatrix_EffBkgShape","",34,ZPtBins,34,ZPtBins); hZPtMatrix_EffBkgShape->Sumw2();
    
    TH1D *hPhiStarReco  = new TH1D("hPhiStarReco","",27,PhiStarBins); hPhiStarReco->Sumw2();
    TH1D *hPhiStarTruth  = new TH1D("hPhiStarTruth","",27,PhiStarBins); hPhiStarTruth->Sumw2();
    TH2D *hPhiStarMatrix  = new TH2D("hPhiStarMatrix","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix->Sumw2();
    TH1D *hPhiStarReco_EffBin  = new TH1D("hPhiStarReco_EffBin","",27,PhiStarBins); hPhiStarReco_EffBin->Sumw2();
    TH2D *hPhiStarMatrix_EffBin  = new TH2D("hPhiStarMatrix_EffBin","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix_EffBin->Sumw2();
    TH1D *hPhiStarReco_EffStatUp  = new TH1D("hPhiStarReco_EffStatUp","",27,PhiStarBins); hPhiStarReco_EffStatUp->Sumw2();
    TH2D *hPhiStarMatrix_EffStatUp  = new TH2D("hPhiStarMatrix_EffStatUp","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix_EffStatUp->Sumw2();
    TH1D *hPhiStarReco_EffStatDown  = new TH1D("hPhiStarReco_EffStatDown","",27,PhiStarBins); hPhiStarReco_EffStatDown->Sumw2();
    TH2D *hPhiStarMatrix_EffStatDown  = new TH2D("hPhiStarMatrix_EffStatDown","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix_EffStatDown->Sumw2();
    TH1D *hPhiStarReco_EffSigShape  = new TH1D("hPhiStarReco_EffSigShape","",27,PhiStarBins); hPhiStarReco_EffSigShape->Sumw2();
    TH2D *hPhiStarMatrix_EffSigShape  = new TH2D("hPhiStarMatrix_EffSigShape","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix_EffSigShape->Sumw2();
    TH1D *hPhiStarReco_EffBkgShape  = new TH1D("hPhiStarReco_EffBkgShape","",27,PhiStarBins); hPhiStarReco_EffBkgShape->Sumw2();
    TH2D *hPhiStarMatrix_EffBkgShape  = new TH2D("hPhiStarMatrix_EffBkgShape","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix_EffBkgShape->Sumw2();
    
    
    TH1D *hZRapReco  = new TH1D("hZRapReco","",24,0,2.4); hZRapReco->Sumw2();
    TH1D *hZRapTruth  = new TH1D("hZRapTruth","",24,0,2.4); hZRapTruth->Sumw2();
    TH2D *hZRapMatrix  = new TH2D("hZRapMatrix","",24,0,2.4,24,0,2.4); hZRapMatrix->Sumw2();
    TH1D *hZRapReco_EffBin  = new TH1D("hZRapReco_EffBin","",24,0,2.4); hZRapReco_EffBin->Sumw2();
    TH2D *hZRapMatrix_EffBin  = new TH2D("hZRapMatrix_EffBin","",24,0,2.4,24,0,2.4); hZRapMatrix_EffBin->Sumw2();
    TH1D *hZRapReco_EffStatUp  = new TH1D("hZRapReco_EffStatUp","",24,0,2.4); hZRapReco_EffStatUp->Sumw2();
    TH2D *hZRapMatrix_EffStatUp  = new TH2D("hZRapMatrix_EffStatUp","",24,0,2.4,24,0,2.4); hZRapMatrix_EffStatUp->Sumw2();
    TH1D *hZRapReco_EffStatDown  = new TH1D("hZRapReco_EffStatDown","",24,0,2.4); hZRapReco_EffStatDown->Sumw2();
    TH2D *hZRapMatrix_EffStatDown  = new TH2D("hZRapMatrix_EffStatDown","",24,0,2.4,24,0,2.4); hZRapMatrix_EffStatDown->Sumw2();
    TH1D *hZRapReco_EffSigShape  = new TH1D("hZRapReco_EffSigShape","",24,0,2.4); hZRapReco_EffSigShape->Sumw2();
    TH2D *hZRapMatrix_EffSigShape  = new TH2D("hZRapMatrix_EffSigShape","",24,0,2.4,24,0,2.4); hZRapMatrix_EffSigShape->Sumw2();
    TH1D *hZRapReco_EffBkgShape  = new TH1D("hZRapReco_EffBkgShape","",24,0,2.4); hZRapReco_EffBkgShape->Sumw2();
    TH2D *hZRapMatrix_EffBkgShape  = new TH2D("hZRapMatrix_EffBkgShape","",24,0,2.4,24,0,2.4); hZRapMatrix_EffBkgShape->Sumw2();
    
    
    TH1D *hLep1PtReco  = new TH1D("hLep1PtReco","",25,Lep1PtBins); hLep1PtReco->Sumw2();
    TH1D *hLep1PtTruth  = new TH1D("hLep1PtTruth","",25,Lep1PtBins); hLep1PtTruth->Sumw2();
    TH2D *hLep1PtMatrix  = new TH2D("hLep1PtMatrix","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix->Sumw2();
    TH1D *hLep1PtReco_EffBin  = new TH1D("hLep1PtReco_EffBin","",25,Lep1PtBins); hLep1PtReco_EffBin->Sumw2();
    TH2D *hLep1PtMatrix_EffBin  = new TH2D("hLep1PtMatrix_EffBin","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix_EffBin->Sumw2();
    TH1D *hLep1PtReco_EffStatUp  = new TH1D("hLep1PtReco_EffStatUp","",25,Lep1PtBins); hLep1PtReco_EffStatUp->Sumw2();
    TH2D *hLep1PtMatrix_EffStatUp  = new TH2D("hLep1PtMatrix_EffStatUp","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix_EffStatUp->Sumw2();
    TH1D *hLep1PtReco_EffStatDown  = new TH1D("hLep1PtReco_EffStatDown","",25,Lep1PtBins); hLep1PtReco_EffStatDown->Sumw2();
    TH2D *hLep1PtMatrix_EffStatDown  = new TH2D("hLep1PtMatrix_EffStatDown","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix_EffStatDown->Sumw2();
    TH1D *hLep1PtReco_EffSigShape  = new TH1D("hLep1PtReco_EffSigShape","",25,Lep1PtBins); hLep1PtReco_EffSigShape->Sumw2();
    TH2D *hLep1PtMatrix_EffSigShape  = new TH2D("hLep1PtMatrix_EffSigShape","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix_EffSigShape->Sumw2();
    TH1D *hLep1PtReco_EffBkgShape  = new TH1D("hLep1PtReco_EffBkgShape","",25,Lep1PtBins); hLep1PtReco_EffBkgShape->Sumw2();
    TH2D *hLep1PtMatrix_EffBkgShape  = new TH2D("hLep1PtMatrix_EffBkgShape","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix_EffBkgShape->Sumw2();
    
    TH1D *hLep2PtReco  = new TH1D("hLep2PtReco","",20,Lep2PtBins); hLep2PtReco->Sumw2();
    TH1D *hLep2PtTruth  = new TH1D("hLep2PtTruth","",20,Lep2PtBins); hLep2PtTruth->Sumw2();
    TH2D *hLep2PtMatrix  = new TH2D("hLep2PtMatrix","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix->Sumw2();
    TH1D *hLep2PtReco_EffBin  = new TH1D("hLep2PtReco_EffBin","",20,Lep2PtBins); hLep2PtReco_EffBin->Sumw2();
    TH2D *hLep2PtMatrix_EffBin  = new TH2D("hLep2PtMatrix_EffBin","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix_EffBin->Sumw2();
    TH1D *hLep2PtReco_EffStatUp  = new TH1D("hLep2PtReco_EffStatUp","",20,Lep2PtBins); hLep2PtReco_EffStatUp->Sumw2();
  TH2D *hLep2PtMatrix_EffStatUp  = new TH2D("hLep2PtMatrix_EffStatUp","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix_EffStatUp->Sumw2();
  TH1D *hLep2PtReco_EffStatDown  = new TH1D("hLep2PtReco_EffStatDown","",20,Lep2PtBins); hLep2PtReco_EffStatDown->Sumw2();
  TH2D *hLep2PtMatrix_EffStatDown  = new TH2D("hLep2PtMatrix_EffStatDown","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix_EffStatDown->Sumw2();
  TH1D *hLep2PtReco_EffSigShape  = new TH1D("hLep2PtReco_EffSigShape","",20,Lep2PtBins); hLep2PtReco_EffSigShape->Sumw2();
  TH2D *hLep2PtMatrix_EffSigShape  = new TH2D("hLep2PtMatrix_EffSigShape","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix_EffSigShape->Sumw2();
  TH1D *hLep2PtReco_EffBkgShape  = new TH1D("hLep2PtReco_EffBkgShape","",20,Lep2PtBins); hLep2PtReco_EffBkgShape->Sumw2();
  TH2D *hLep2PtMatrix_EffBkgShape  = new TH2D("hLep2PtMatrix_EffBkgShape","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix_EffBkgShape->Sumw2();

  TH1D *hLepNegPtReco  = new TH1D("hLepNegPtReco","",25,LepNegPtBins); hLepNegPtReco->Sumw2();
  TH1D *hLepNegPtTruth  = new TH1D("hLepNegPtTruth","",25,LepNegPtBins); hLepNegPtTruth->Sumw2();
  TH2D *hLepNegPtMatrix  = new TH2D("hLepNegPtMatrix","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix->Sumw2();
  TH1D *hLepNegPtReco_EffBin  = new TH1D("hLepNegPtReco_EffBin","",25,LepNegPtBins); hLepNegPtReco_EffBin->Sumw2();
  TH2D *hLepNegPtMatrix_EffBin  = new TH2D("hLepNegPtMatrix_EffBin","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix_EffBin->Sumw2();
  TH1D *hLepNegPtReco_EffStatUp  = new TH1D("hLepNegPtReco_EffStatUp","",25,LepNegPtBins); hLepNegPtReco_EffStatUp->Sumw2();
  TH2D *hLepNegPtMatrix_EffStatUp  = new TH2D("hLepNegPtMatrix_EffStatUp","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix_EffStatUp->Sumw2();
  TH1D *hLepNegPtReco_EffStatDown  = new TH1D("hLepNegPtReco_EffStatDown","",25,LepNegPtBins); hLepNegPtReco_EffStatDown->Sumw2();
  TH2D *hLepNegPtMatrix_EffStatDown  = new TH2D("hLepNegPtMatrix_EffStatDown","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix_EffStatDown->Sumw2();
  TH1D *hLepNegPtReco_EffSigShape  = new TH1D("hLepNegPtReco_EffSigShape","",25,LepNegPtBins); hLepNegPtReco_EffSigShape->Sumw2();
  TH2D *hLepNegPtMatrix_EffSigShape  = new TH2D("hLepNegPtMatrix_EffSigShape","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix_EffSigShape->Sumw2();
  TH1D *hLepNegPtReco_EffBkgShape  = new TH1D("hLepNegPtReco_EffBkgShape","",25,LepNegPtBins); hLepNegPtReco_EffBkgShape->Sumw2();
  TH2D *hLepNegPtMatrix_EffBkgShape  = new TH2D("hLepNegPtMatrix_EffBkgShape","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix_EffBkgShape->Sumw2();

  TH1D *hLepPosPtReco  = new TH1D("hLepPosPtReco","",25,LepPosPtBins); hLepPosPtReco->Sumw2();
  TH1D *hLepPosPtTruth  = new TH1D("hLepPosPtTruth","",25,LepPosPtBins); hLepPosPtTruth->Sumw2();
  TH2D *hLepPosPtMatrix  = new TH2D("hLepPosPtMatrix","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix->Sumw2();
  TH1D *hLepPosPtReco_EffBin  = new TH1D("hLepPosPtReco_EffBin","",25,LepPosPtBins); hLepPosPtReco_EffBin->Sumw2();
  TH2D *hLepPosPtMatrix_EffBin  = new TH2D("hLepPosPtMatrix_EffBin","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix_EffBin->Sumw2();
  TH1D *hLepPosPtReco_EffStatUp  = new TH1D("hLepPosPtReco_EffStatUp","",25,LepPosPtBins); hLepPosPtReco_EffStatUp->Sumw2();
  TH2D *hLepPosPtMatrix_EffStatUp  = new TH2D("hLepPosPtMatrix_EffStatUp","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix_EffStatUp->Sumw2();
  TH1D *hLepPosPtReco_EffStatDown  = new TH1D("hLepPosPtReco_EffStatDown","",25,LepPosPtBins); hLepPosPtReco_EffStatDown->Sumw2();
  TH2D *hLepPosPtMatrix_EffStatDown  = new TH2D("hLepPosPtMatrix_EffStatDown","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix_EffStatDown->Sumw2();
  TH1D *hLepPosPtReco_EffSigShape  = new TH1D("hLepPosPtReco_EffSigShape","",25,LepPosPtBins); hLepPosPtReco_EffSigShape->Sumw2();
  TH2D *hLepPosPtMatrix_EffSigShape  = new TH2D("hLepPosPtMatrix_EffSigShape","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix_EffSigShape->Sumw2();
  TH1D *hLepPosPtReco_EffBkgShape  = new TH1D("hLepPosPtReco_EffBkgShape","",25,LepPosPtBins); hLepPosPtReco_EffBkgShape->Sumw2();
  TH2D *hLepPosPtMatrix_EffBkgShape  = new TH2D("hLepPosPtMatrix_EffBkgShape","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix_EffBkgShape->Sumw2();

  TH1D *hLep1EtaReco  = new TH1D("hLep1EtaReco","",24,0,2.4); hLep1EtaReco->Sumw2();
  TH1D *hLep1EtaTruth  = new TH1D("hLep1EtaTruth","",24,0,2.4); hLep1EtaTruth->Sumw2();
  TH2D *hLep1EtaMatrix  = new TH2D("hLep1EtaMatrix","",24,0,2.4,24,0,2.4); hLep1EtaMatrix->Sumw2();
  TH1D *hLep1EtaReco_EffBin  = new TH1D("hLep1EtaReco_EffBin","",24,0,2.4); hLep1EtaReco_EffBin->Sumw2();
  TH2D *hLep1EtaMatrix_EffBin  = new TH2D("hLep1EtaMatrix_EffBin","",24,0,2.4,24,0,2.4); hLep1EtaMatrix_EffBin->Sumw2();
  TH1D *hLep1EtaReco_EffStatUp  = new TH1D("hLep1EtaReco_EffStatUp","",24,0,2.4); hLep1EtaReco_EffStatUp->Sumw2();
  TH2D *hLep1EtaMatrix_EffStatUp  = new TH2D("hLep1EtaMatrix_EffStatUp","",24,0,2.4,24,0,2.4); hLep1EtaMatrix_EffStatUp->Sumw2();
  TH1D *hLep1EtaReco_EffStatDown  = new TH1D("hLep1EtaReco_EffStatDown","",24,0,2.4); hLep1EtaReco_EffStatDown->Sumw2();
  TH2D *hLep1EtaMatrix_EffStatDown  = new TH2D("hLep1EtaMatrix_EffStatDown","",24,0,2.4,24,0,2.4); hLep1EtaMatrix_EffStatDown->Sumw2();
  TH1D *hLep1EtaReco_EffSigShape  = new TH1D("hLep1EtaReco_EffSigShape","",24,0,2.4); hLep1EtaReco_EffSigShape->Sumw2();
  TH2D *hLep1EtaMatrix_EffSigShape  = new TH2D("hLep1EtaMatrix_EffSigShape","",24,0,2.4,24,0,2.4); hLep1EtaMatrix_EffSigShape->Sumw2();
  TH1D *hLep1EtaReco_EffBkgShape  = new TH1D("hLep1EtaReco_EffBkgShape","",24,0,2.4); hLep1EtaReco_EffBkgShape->Sumw2();
  TH2D *hLep1EtaMatrix_EffBkgShape  = new TH2D("hLep1EtaMatrix_EffBkgShape","",24,0,2.4,24,0,2.4); hLep1EtaMatrix_EffBkgShape->Sumw2();


  TH1D *hLep2EtaReco  = new TH1D("hLep2EtaReco","",24,0,2.4); hLep2EtaReco->Sumw2();
  TH1D *hLep2EtaTruth  = new TH1D("hLep2EtaTruth","",24,0,2.4); hLep2EtaTruth->Sumw2();
  TH2D *hLep2EtaMatrix  = new TH2D("hLep2EtaMatrix","",24,0,2.4,24,0,2.4); hLep2EtaMatrix->Sumw2();
  TH1D *hLep2EtaReco_EffBin  = new TH1D("hLep2EtaReco_EffBin","",24,0,2.4); hLep2EtaReco_EffBin->Sumw2();
  TH2D *hLep2EtaMatrix_EffBin  = new TH2D("hLep2EtaMatrix_EffBin","",24,0,2.4,24,0,2.4); hLep2EtaMatrix_EffBin->Sumw2();
  TH1D *hLep2EtaReco_EffStatUp  = new TH1D("hLep2EtaReco_EffStatUp","",24,0,2.4); hLep2EtaReco_EffStatUp->Sumw2();
  TH2D *hLep2EtaMatrix_EffStatUp  = new TH2D("hLep2EtaMatrix_EffStatUp","",24,0,2.4,24,0,2.4); hLep2EtaMatrix_EffStatUp->Sumw2();
  TH1D *hLep2EtaReco_EffStatDown  = new TH1D("hLep2EtaReco_EffStatDown","",24,0,2.4); hLep2EtaReco_EffStatDown->Sumw2();
  TH2D *hLep2EtaMatrix_EffStatDown  = new TH2D("hLep2EtaMatrix_EffStatDown","",24,0,2.4,24,0,2.4); hLep2EtaMatrix_EffStatDown->Sumw2();
  TH1D *hLep2EtaReco_EffSigShape  = new TH1D("hLep2EtaReco_EffSigShape","",24,0,2.4); hLep2EtaReco_EffSigShape->Sumw2();
  TH2D *hLep2EtaMatrix_EffSigShape  = new TH2D("hLep2EtaMatrix_EffSigShape","",24,0,2.4,24,0,2.4); hLep2EtaMatrix_EffSigShape->Sumw2();
  TH1D *hLep2EtaReco_EffBkgShape  = new TH1D("hLep2EtaReco_EffBkgShape","",24,0,2.4); hLep2EtaReco_EffBkgShape->Sumw2();
  TH2D *hLep2EtaMatrix_EffBkgShape  = new TH2D("hLep2EtaMatrix_EffBkgShape","",24,0,2.4,24,0,2.4); hLep2EtaMatrix_EffBkgShape->Sumw2();

  
  //
  // loop over events
  //
  cout<<"Beginning of Loop: nentries"<<intree->GetEntries() << endl;
  sw_[0].Start();

  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    if ( ientry %10000 == 0 )
      {
	sw_[0].Stop();
	cout<< "Took " <<  sw_[0] .RealTime() << "s CPU " << sw_[0] . CpuTime()<< endl;
	sw_[0].Reset();
	sw_[0].Start();
      }
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

    TLorentzVector *dilep=new TLorentzVector(0,0,0,0);
    dilep->operator+=(l1);
    dilep->operator+=(l2);

    float phiacop=0;
    float costhetastar=0;
    float phistar=0;
    phiacop=TMath::Pi()-fabs(l1.DeltaPhi(l2));
    if(lq1<0) costhetastar=tanh(float((l1.Rapidity()-l2.Rapidity())/2));
    else costhetastar=tanh(float((l2.Rapidity()-l1.Rapidity())/2));
    phistar=tan(phiacop/2)*sqrt(1-pow(costhetastar,2));
    
    TLorentzVector *gendilep=new TLorentzVector(0,0,0,0);
    gendilep->operator+=(*genlep1);
    gendilep->operator+=(*genlep2);

    float genphiacop=0;
    float gencosthetastar=0;
    float genphistar=0;

    genphiacop=TMath::Pi()-fabs(genlep1->DeltaPhi(*genlep2));
    if(genq1<0) gencosthetastar=tanh(float((genlep1->Rapidity()-genlep2->Rapidity())/2));
    else gencosthetastar=tanh(float((genlep2->Rapidity()-genlep1->Rapidity())/2));
    genphistar=tan(genphiacop/2)*sqrt(1-pow(gencosthetastar,2));
    
    bool isReco=false;
    bool isGen=false;

    if(triggerDec&&goodPV&&matchTrigger&&nlep>=2&&lq1!=lq2&&dilep->M()>MASS_LOW&&dilep->M()<MASS_HIGH&&l1.Pt()>=PT_CUT&&l2.Pt()>=PT_CUT&&fabs(l1.Eta())<=ETA_CUT&&fabs(l2.Eta())<=ETA_CUT)
      {
	isReco=true;
      }
    if(ngenlep>=2&&genq1!=genq2&&gendilep->M()>MASS_LOW&&gendilep->M()<MASS_HIGH&&genlep1->Pt()>=PT_CUT&&genlep2->Pt()>=PT_CUT&&fabs(genlep1->Eta())<=ETA_CUT&&fabs(genlep2->Eta())<=ETA_CUT)
      {
	isGen=true;
      }

    Double_t genweight = 1;
    genweight *= scale1fbGen*lumi;
    Double_t weight = 1;
    weight *=scale1fb*lumi;

    sw_[1].Stop();

    if(isReco)
      {
	hMassMC ->Fill(dilep->M(),weight*corr);
	hZPtReco ->Fill(dilep->Pt(),weight*corr);
	hPhiStarReco ->Fill(phistar,weight*corr);
	hZRapReco ->Fill(fabs(dilep->Rapidity()),weight*corr);
	hLep1PtReco ->Fill(l1.Pt(),weight*corr);
	hLep2PtReco ->Fill(l2.Pt(),weight*corr);
	hLep1EtaReco ->Fill(fabs(l1.Eta()),weight*corr);
	hLep2EtaReco ->Fill(fabs(l2.Eta()),weight*corr);
	hZPtReco_EffBin ->Fill(dilep->Pt(),weight*corr2Bin);
	hPhiStarReco_EffBin ->Fill(phistar,weight*corr2Bin);
	hZRapReco_EffBin ->Fill(fabs(dilep->Rapidity()),weight*corr2Bin);
	hLep1PtReco_EffBin ->Fill(l1.Pt(),weight*corr2Bin);
	hLep2PtReco_EffBin ->Fill(l2.Pt(),weight*corr2Bin);
	hLep1EtaReco_EffBin ->Fill(fabs(l1.Eta()),weight*corr2Bin);
	hLep2EtaReco_EffBin ->Fill(fabs(l2.Eta()),weight*corr2Bin);
	hZPtReco_EffStatUp ->Fill(dilep->Pt(),weight*corrUp);
	hPhiStarReco_EffStatUp ->Fill(phistar,weight*corrUp);
	hZRapReco_EffStatUp ->Fill(fabs(dilep->Rapidity()),weight*corrUp);
	hLep1PtReco_EffStatUp ->Fill(l1.Pt(),weight*corrUp);
	hLep2PtReco_EffStatUp ->Fill(l2.Pt(),weight*corrUp);
	hLep1EtaReco_EffStatUp ->Fill(fabs(l1.Eta()),weight*corrUp);
	hLep2EtaReco_EffStatUp ->Fill(fabs(l2.Eta()),weight*corrUp);
	hZPtReco_EffStatDown ->Fill(dilep->Pt(),weight*corrDown);
	hPhiStarReco_EffStatDown ->Fill(phistar,weight*corrDown);
	hZRapReco_EffStatDown ->Fill(fabs(dilep->Rapidity()),weight*corrDown);
	hLep1PtReco_EffStatDown ->Fill(l1.Pt(),weight*corrDown);
	hLep2PtReco_EffStatDown ->Fill(l2.Pt(),weight*corrDown);
	hLep1EtaReco_EffStatDown ->Fill(fabs(l1.Eta()),weight*corrDown);
	hLep2EtaReco_EffStatDown ->Fill(fabs(l2.Eta()),weight*corrDown);
	hZPtReco_EffSigShape ->Fill(dilep->Pt(),weight*corrSigShape);
	hPhiStarReco_EffSigShape ->Fill(phistar,weight*corrSigShape);
	hZRapReco_EffSigShape ->Fill(fabs(dilep->Rapidity()),weight*corrSigShape);
	hLep1PtReco_EffSigShape ->Fill(l1.Pt(),weight*corrSigShape);
	hLep2PtReco_EffSigShape ->Fill(l2.Pt(),weight*corrSigShape);
	hLep1EtaReco_EffSigShape ->Fill(fabs(l1.Eta()),weight*corrSigShape);
	hLep2EtaReco_EffSigShape ->Fill(fabs(l2.Eta()),weight*corrSigShape);
	hZPtReco_EffBkgShape ->Fill(dilep->Pt(),weight*corrBkgShape);
	hPhiStarReco_EffBkgShape ->Fill(phistar,weight*corrBkgShape);
	hZRapReco_EffBkgShape ->Fill(fabs(dilep->Rapidity()),weight*corrBkgShape);
	hLep1PtReco_EffBkgShape ->Fill(l1.Pt(),weight*corrBkgShape);
	hLep2PtReco_EffBkgShape ->Fill(l2.Pt(),weight*corrBkgShape);
	hLep1EtaReco_EffBkgShape ->Fill(fabs(l1.Eta()),weight*corrBkgShape);
	hLep2EtaReco_EffBkgShape ->Fill(fabs(l2.Eta()),weight*corrBkgShape);
	
	if(lq1<0)
	  {
	    hLepNegPtReco ->Fill(l1.Pt(),weight*corr);
	    hLepPosPtReco ->Fill(l2.Pt(),weight*corr);
	    hLepNegPtReco_EffBin ->Fill(l1.Pt(),weight*corr2Bin);
	    hLepPosPtReco_EffBin ->Fill(l2.Pt(),weight*corr2Bin);
	    hLepNegPtReco_EffStatUp ->Fill(l1.Pt(),weight*corrUp);
	    hLepPosPtReco_EffStatUp ->Fill(l2.Pt(),weight*corrUp);
	    hLepNegPtReco_EffStatDown ->Fill(l1.Pt(),weight*corrDown);
	    hLepPosPtReco_EffStatDown ->Fill(l2.Pt(),weight*corrDown);
	    hLepNegPtReco_EffSigShape ->Fill(l1.Pt(),weight*corrSigShape);
	    hLepPosPtReco_EffSigShape ->Fill(l2.Pt(),weight*corrSigShape);
	    hLepNegPtReco_EffBkgShape ->Fill(l1.Pt(),weight*corrBkgShape);
	    hLepPosPtReco_EffBkgShape ->Fill(l2.Pt(),weight*corrBkgShape);
	  }
	else
	  {
	    hLepNegPtReco ->Fill(l2.Pt(),weight*corr);
	    hLepPosPtReco ->Fill(l1.Pt(),weight*corr);
	    hLepNegPtReco_EffBin ->Fill(l2.Pt(),weight*corr2Bin);
	    hLepPosPtReco_EffBin ->Fill(l1.Pt(),weight*corr2Bin);
	    hLepNegPtReco_EffStatUp ->Fill(l2.Pt(),weight*corrUp);
	    hLepPosPtReco_EffStatUp ->Fill(l1.Pt(),weight*corrUp);
	    hLepNegPtReco_EffStatDown ->Fill(l2.Pt(),weight*corrDown);
	    hLepPosPtReco_EffStatDown ->Fill(l1.Pt(),weight*corrDown);
	    hLepNegPtReco_EffSigShape ->Fill(l2.Pt(),weight*corrSigShape);
	    hLepPosPtReco_EffSigShape ->Fill(l1.Pt(),weight*corrSigShape);
	    hLepNegPtReco_EffBkgShape ->Fill(l2.Pt(),weight*corrBkgShape);
	    hLepPosPtReco_EffBkgShape ->Fill(l1.Pt(),weight*corrBkgShape);
	  }

      }
    if(isGen)
      {
	hZPtTruth ->Fill(gendilep->Pt(),genweight);
	hPhiStarTruth ->Fill(genphistar,genweight);
	hZRapTruth ->Fill(fabs(gendilep->Rapidity()),genweight);
	hLep1PtTruth ->Fill(genlep1->Pt(),genweight);
	hLep2PtTruth ->Fill(genlep2->Pt(),genweight);
	hLep1EtaTruth ->Fill(fabs(genlep1->Eta()),genweight);
	hLep2EtaTruth ->Fill(fabs(genlep2->Eta()),genweight);
	if(genq1<0)
	  {
	    hLepNegPtTruth ->Fill(genlep1->Pt(),genweight);
	    hLepPosPtTruth ->Fill(genlep2->Pt(),genweight);
	  }
	else
	  {
	    hLepNegPtTruth ->Fill(genlep2->Pt(),genweight);
	    hLepPosPtTruth ->Fill(genlep1->Pt(),genweight);
	  }
      }
    if(isReco&&isGen)
      {
	hZPtMatrix ->Fill(gendilep->Pt(),dilep->Pt(),weight*corr);
	hPhiStarMatrix ->Fill(genphistar,phistar,weight*corr);
	hZRapMatrix ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corr);
	hLep1PtMatrix ->Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	hLep2PtMatrix ->Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	hLep1EtaMatrix ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corr);
	hLep2EtaMatrix ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corr);
	hZPtMatrix_EffBin ->Fill(gendilep->Pt(),dilep->Pt(),weight*corr2Bin);
	hPhiStarMatrix_EffBin ->Fill(genphistar,phistar,weight*corr2Bin);
	hZRapMatrix_EffBin ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corr2Bin);
	hLep1PtMatrix_EffBin ->Fill(genlep1->Pt(),l1.Pt(),weight*corr2Bin);
	hLep2PtMatrix_EffBin ->Fill(genlep2->Pt(),l2.Pt(),weight*corr2Bin);
	hLep1EtaMatrix_EffBin ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corr2Bin);
	hLep2EtaMatrix_EffBin ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corr2Bin);
	hZPtMatrix_EffStatUp ->Fill(gendilep->Pt(),dilep->Pt(),weight*corrUp);
	hPhiStarMatrix_EffStatUp ->Fill(genphistar,phistar,weight*corrUp);

	hZRapMatrix_EffStatUp ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corrUp);
	hLep1PtMatrix_EffStatUp ->Fill(genlep1->Pt(),l1.Pt(),weight*corrUp);
	hLep2PtMatrix_EffStatUp ->Fill(genlep2->Pt(),l2.Pt(),weight*corrUp);
	hLep1EtaMatrix_EffStatUp ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corrUp);
	hLep2EtaMatrix_EffStatUp ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corrUp);
	hZPtMatrix_EffStatDown ->Fill(gendilep->Pt(),dilep->Pt(),weight*corrDown);
	hPhiStarMatrix_EffStatDown ->Fill(genphistar,phistar,weight*corrDown);
	hZRapMatrix_EffStatDown ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corrDown);
	hLep1PtMatrix_EffStatDown ->Fill(genlep1->Pt(),l1.Pt(),weight*corrDown);
	hLep2PtMatrix_EffStatDown ->Fill(genlep2->Pt(),l2.Pt(),weight*corrDown);
	hLep1EtaMatrix_EffStatDown ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corrDown);
	hLep2EtaMatrix_EffStatDown ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corrDown);
	hZPtMatrix_EffSigShape ->Fill(gendilep->Pt(),dilep->Pt(),weight*corrSigShape);
	hPhiStarMatrix_EffSigShape ->Fill(genphistar,phistar,weight*corrSigShape);
	hZRapMatrix_EffSigShape ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corrSigShape);
	hLep1PtMatrix_EffSigShape ->Fill(genlep1->Pt(),l1.Pt(),weight*corrSigShape);
	hLep2PtMatrix_EffSigShape ->Fill(genlep2->Pt(),l2.Pt(),weight*corrSigShape);
	hLep1EtaMatrix_EffSigShape ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corrSigShape);
	hLep2EtaMatrix_EffSigShape ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corrSigShape);
	hZPtMatrix_EffBkgShape ->Fill(gendilep->Pt(),dilep->Pt(),weight*corrBkgShape);
	hPhiStarMatrix_EffBkgShape ->Fill(genphistar,phistar,weight*corrBkgShape);
	hZRapMatrix_EffBkgShape ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corrBkgShape);
	hLep1PtMatrix_EffBkgShape ->Fill(genlep1->Pt(),l1.Pt(),weight*corrBkgShape);
	hLep2PtMatrix_EffBkgShape ->Fill(genlep2->Pt(),l2.Pt(),weight*corrBkgShape);
	hLep1EtaMatrix_EffBkgShape ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corrBkgShape);
	hLep2EtaMatrix_EffBkgShape ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corrBkgShape);

	if(lq1<0&&genq1<0)
	  {
	    hLepNegPtMatrix ->Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	    hLepNegPtMatrix_EffBin ->Fill(genlep1->Pt(),l1.Pt(),weight*corr2Bin);
	    hLepPosPtMatrix_EffBin ->Fill(genlep2->Pt(),l2.Pt(),weight*corr2Bin);
	    hLepNegPtMatrix_EffStatUp ->Fill(genlep1->Pt(),l1.Pt(),weight*corrUp);
	    hLepPosPtMatrix_EffStatUp ->Fill(genlep2->Pt(),l2.Pt(),weight*corrUp);
	    hLepNegPtMatrix_EffStatDown ->Fill(genlep1->Pt(),l1.Pt(),weight*corrDown);
	    hLepPosPtMatrix_EffStatDown ->Fill(genlep2->Pt(),l2.Pt(),weight*corrDown);
	    hLepNegPtMatrix_EffSigShape ->Fill(genlep1->Pt(),l1.Pt(),weight*corrSigShape);
	    hLepPosPtMatrix_EffSigShape ->Fill(genlep2->Pt(),l2.Pt(),weight*corrSigShape);
	    hLepNegPtMatrix_EffBkgShape ->Fill(genlep1->Pt(),l1.Pt(),weight*corrBkgShape);
	    hLepPosPtMatrix_EffBkgShape ->Fill(genlep2->Pt(),l2.Pt(),weight*corrBkgShape);
	  }
	else if(lq1<0&&genq1>0)
	  {
	    hLepNegPtMatrix ->Fill(genlep2->Pt(),l1.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep1->Pt(),l2.Pt(),weight*corr);
	    hLepNegPtMatrix_EffBin ->Fill(genlep2->Pt(),l1.Pt(),weight*corr2Bin);
	    hLepPosPtMatrix_EffBin ->Fill(genlep1->Pt(),l2.Pt(),weight*corr2Bin);
	    hLepNegPtMatrix_EffStatUp ->Fill(genlep2->Pt(),l1.Pt(),weight*corrUp);
	    hLepPosPtMatrix_EffStatUp ->Fill(genlep1->Pt(),l2.Pt(),weight*corrUp);
	    hLepNegPtMatrix_EffStatDown ->Fill(genlep2->Pt(),l1.Pt(),weight*corrDown);
	    hLepPosPtMatrix_EffStatDown ->Fill(genlep1->Pt(),l2.Pt(),weight*corrDown);
	    hLepNegPtMatrix_EffSigShape ->Fill(genlep2->Pt(),l1.Pt(),weight*corrSigShape);
	    hLepPosPtMatrix_EffSigShape ->Fill(genlep1->Pt(),l2.Pt(),weight*corrSigShape);
	    hLepNegPtMatrix_EffBkgShape ->Fill(genlep2->Pt(),l1.Pt(),weight*corrBkgShape);
	    hLepPosPtMatrix_EffBkgShape ->Fill(genlep1->Pt(),l2.Pt(),weight*corrBkgShape);
	  }
	else if(lq1>0&&genq1<0)
	  {
	    hLepNegPtMatrix ->Fill(genlep1->Pt(),l2.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep2->Pt(),l1.Pt(),weight*corr);
	    hLepNegPtMatrix_EffBin ->Fill(genlep1->Pt(),l2.Pt(),weight*corr2Bin);
	    hLepPosPtMatrix_EffBin ->Fill(genlep2->Pt(),l1.Pt(),weight*corr2Bin);
	    hLepNegPtMatrix_EffStatUp ->Fill(genlep1->Pt(),l2.Pt(),weight*corrUp);
	    hLepPosPtMatrix_EffStatUp ->Fill(genlep2->Pt(),l1.Pt(),weight*corrUp);
	    hLepNegPtMatrix_EffStatDown ->Fill(genlep1->Pt(),l2.Pt(),weight*corrDown);
	    hLepPosPtMatrix_EffStatDown ->Fill(genlep2->Pt(),l1.Pt(),weight*corrDown);
	    hLepNegPtMatrix_EffSigShape ->Fill(genlep1->Pt(),l2.Pt(),weight*corrSigShape);
	    hLepPosPtMatrix_EffSigShape ->Fill(genlep2->Pt(),l1.Pt(),weight*corrSigShape);
	    hLepNegPtMatrix_EffBkgShape ->Fill(genlep1->Pt(),l2.Pt(),weight*corrBkgShape);
	    hLepPosPtMatrix_EffBkgShape ->Fill(genlep2->Pt(),l1.Pt(),weight*corrBkgShape);
	  }
	else if(lq1>0&&genq1>0)
	  {
	    hLepNegPtMatrix ->Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	    hLepNegPtMatrix_EffBin ->Fill(genlep2->Pt(),l2.Pt(),weight*corr2Bin);
	    hLepPosPtMatrix_EffBin ->Fill(genlep1->Pt(),l1.Pt(),weight*corr2Bin);
	    hLepNegPtMatrix_EffStatUp ->Fill(genlep2->Pt(),l2.Pt(),weight*corrUp);
	    hLepPosPtMatrix_EffStatUp ->Fill(genlep1->Pt(),l1.Pt(),weight*corrUp);
	    hLepNegPtMatrix_EffStatDown ->Fill(genlep2->Pt(),l2.Pt(),weight*corrDown);
	    hLepPosPtMatrix_EffStatDown ->Fill(genlep1->Pt(),l1.Pt(),weight*corrDown);
	    hLepNegPtMatrix_EffSigShape ->Fill(genlep2->Pt(),l2.Pt(),weight*corrSigShape);
	    hLepPosPtMatrix_EffSigShape ->Fill(genlep1->Pt(),l1.Pt(),weight*corrSigShape);
	    hLepNegPtMatrix_EffBkgShape ->Fill(genlep2->Pt(),l2.Pt(),weight*corrBkgShape);
	    hLepPosPtMatrix_EffBkgShape ->Fill(genlep1->Pt(),l1.Pt(),weight*corrBkgShape);
	  }
      }
    delete gendilep;
    delete dilep;
  }
  delete infile;
  infile=0, intree=0; 

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  
  // string buffers
  char pname[100];       // plot name
  char xlabel[100];      // x-axis label
  char ylabel[100];      // y-axis label
  char plabel[100];      // plot label
  char lumitext[100];    // lumi label

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->cd();
  c->SetTopMargin(0.10);
  c->SetBottomMargin(0.15);
  c->SetLeftMargin(0.15);  
  c->SetRightMargin(0.10);  
  c->SetTickx(1);
  c->SetTicky(1); 
  
  if(snamev[isam].CompareTo("zmm",TString::kIgnoreCase)==0)
    {
      sprintf(plabel,"Migration matrix (aMC@NLO)");
    }
  else if(snamev[isam].CompareTo("zmmmg",TString::kIgnoreCase)==0)
    {
      sprintf(plabel,"Migration matrix (MADGRAPH)");
    }
  else sprintf(plabel,"Migration matrix");

  
  gStyle->SetTitleOffset(1.4,"Y");
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("4.2f");

  //
  // Z Pt
  //

  TH2D *hZPtMatrixNorm  = new TH2D("hZPtMatrixNorm","",34,ZPtBins,34,ZPtBins); hZPtMatrixNorm->Sumw2();
  TH1D *hZPtBinEfficiency  = new TH1D("hZPtBinEfficiency","",34,ZPtBins); hZPtBinEfficiency->Sumw2();
  TH1D *hZPtBinPurity  = new TH1D("hZPtBinPurity","",34,ZPtBins); hZPtBinPurity->Sumw2();
  double sumtot=0;
  for(int i=0;i<hZPtMatrix->GetNbinsX();++i)
    {
      double sum=0.;
      for(int j=0;j<hZPtMatrix->GetNbinsY();++j)
	{
	  if(hZPtMatrix->GetBinContent(i+1,j+1)/hZPtTruth->GetBinContent(i+1)>0)
	    {
	      hZPtMatrixNorm->SetBinContent(i+1,j+1,hZPtMatrix->GetBinContent(i+1,j+1)/hZPtTruth->GetBinContent(i+1));
	    }
	  if(hZPtMatrix->GetBinContent(i+1,j+1)<0) hZPtMatrix->SetBinContent(i+1,j+1,0);

	  sum+=hZPtMatrix->GetBinContent(i+1,j+1);

	}
      hZPtBinEfficiency->SetBinContent(i+1,hZPtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i<hZPtMatrix->GetNbinsY();++i)
    {
      double sum=0.;
      for(int j=0;j<hZPtMatrix->GetNbinsX();++j)
	{
	  sum+=hZPtMatrix->GetBinContent(j+1,i+1);
	}
      hZPtBinPurity->SetBinContent(i+1,hZPtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  CPlot plotZPtEff("ZPtBinEfficiency","",xlabel,ylabel);
  plotZPtEff.AddHist1D(hZPtBinEfficiency,"E");
  plotZPtEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZPtEff.SetLogx();
  plotZPtEff.SetLogy(0);
  plotZPtEff.SetYRange(0,1);
  plotZPtEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  CPlot plotZPtPurity("ZPtBinPurity","",xlabel,ylabel);
  plotZPtPurity.AddHist1D(hZPtBinPurity,"E");
  plotZPtPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZPtPurity.SetLogx();
  plotZPtPurity.SetLogy(0);
  plotZPtPurity.SetYRange(0,1);
  plotZPtPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"p_{T}^{#mu^{+}#mu^{-}} (detector level) [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{+}#mu^{-}} (hadron level) [GeV]");
  CPlot plotZPtMatrix("ZPtMatrixNorm","",xlabel,ylabel);
  plotZPtMatrix.AddHist2D(hZPtMatrixNorm,"COLZ");
  plotZPtMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZPtMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotZPtMatrix.SetLogx();
  plotZPtMatrix.SetLogy();
  plotZPtMatrix.Draw(c,kTRUE,format);


  //
  // PhiStar
  //

  TH2D *hPhiStarMatrixNorm  = new TH2D("hPhiStarMatrixNorm","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrixNorm->Sumw2();
  TH1D *hPhiStarBinEfficiency  = new TH1D("hPhiStarBinEfficiency","",27,PhiStarBins); hPhiStarBinEfficiency->Sumw2();
  TH1D *hPhiStarBinPurity  = new TH1D("hPhiStarBinPurity","",27,PhiStarBins); hPhiStarBinPurity->Sumw2();
  sumtot=0;
  for(int i=0;i<hPhiStarMatrix->GetNbinsX();++i)
    {
      double sum=0.;
      for(int j=0;j<hPhiStarMatrix->GetNbinsY();++j)
	{
	  if(hPhiStarMatrix->GetBinContent(i+1,j+1)/hPhiStarTruth->GetBinContent(i+1)>0)
	    {
	      hPhiStarMatrixNorm->SetBinContent(i+1,j+1,hPhiStarMatrix->GetBinContent(i+1,j+1)/hPhiStarTruth->GetBinContent(i+1));
	    }
	  if(hPhiStarMatrix->GetBinContent(i+1,j+1)<0) hPhiStarMatrix->SetBinContent(i+1,j+1,0);

	  sum+=hPhiStarMatrix->GetBinContent(i+1,j+1);

	}
      hPhiStarBinEfficiency->SetBinContent(i+1,hPhiStarMatrix->GetBinContent(i+1,i+1)/sum);
    }
  for(int i=0;i<hPhiStarMatrix->GetNbinsY();++i)
    {
      double sum=0.;
      for(int j=0;j<hPhiStarMatrix->GetNbinsX();++j)
	{
	  sum+=hPhiStarMatrix->GetBinContent(j+1,i+1);
	}
      hPhiStarBinPurity->SetBinContent(i+1,hPhiStarMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"#phi_{#eta}*");
  CPlot plotPhiStarEff("PhiStarBinEfficiency","",xlabel,ylabel);
  plotPhiStarEff.AddHist1D(hPhiStarBinEfficiency,"E");
  plotPhiStarEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotPhiStarEff.SetLogx();
  plotPhiStarEff.SetLogy(0);
  plotPhiStarEff.SetYRange(0,1);
  plotPhiStarEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"#phi_{#eta}*");
  CPlot plotPhiStarPurity("PhiStarBinPurity","",xlabel,ylabel);
  plotPhiStarPurity.AddHist1D(hPhiStarBinPurity,"E");
  plotPhiStarPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotPhiStarPurity.SetLogx();
  plotPhiStarPurity.SetLogy(0);
  plotPhiStarPurity.SetYRange(0,1);
  plotPhiStarPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"#phi_{#eta}* (detector level)");
  sprintf(xlabel,"#phi_{#eta}* (hadron level)");
  CPlot plotPhiStarMatrix("PhiStarMatrixNorm","",xlabel,ylabel);
  plotPhiStarMatrix.AddHist2D(hPhiStarMatrixNorm,"COLZ");
  plotPhiStarMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotPhiStarMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotPhiStarMatrix.SetLogx();
  plotPhiStarMatrix.SetLogy();
  plotPhiStarMatrix.Draw(c,kTRUE,format);

  //
  // Z Rap
  //
  
  TH2D *hZRapMatrixNorm  = new TH2D("hZRapMatrixNorm","",24,0,2.4,24,0,2.4); hZRapMatrixNorm->Sumw2();
  TH1D *hZRapBinEfficiency  = new TH1D("hZRapBinEfficiency","",24,0,2.4); hZRapBinEfficiency->Sumw2();
  TH1D *hZRapBinPurity  = new TH1D("hZRapBinPurity","",24,0,2.4); hZRapBinPurity->Sumw2();
  for(int i=0;i!=hZRapMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hZRapMatrix->GetNbinsY();++j)
	{
	  if(hZRapMatrix->GetBinContent(i+1,j+1)/hZRapTruth->GetBinContent(i+1)>0)
	    {
	      hZRapMatrixNorm->SetBinContent(i+1,j+1,hZRapMatrix->GetBinContent(i+1,j+1)/hZRapTruth->GetBinContent(i+1));
	    }
	  if(hZRapMatrix->GetBinContent(i+1,j+1)<0) hZRapMatrix->SetBinContent(i+1,j+1,0);
	  sum+=hZRapMatrix->GetBinContent(i+1,j+1);
	}
      hZRapBinEfficiency->SetBinContent(i+1,hZRapMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hZRapMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hZRapMatrix->GetNbinsX();++j)
	{
	  sum+=hZRapMatrix->GetBinContent(j+1,i+1);
	}
      hZRapBinPurity->SetBinContent(i+1,hZRapMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"|y^{#mu^{+}#mu^{-}}|");
  CPlot plotZRapEff("ZRapBinEfficiency","",xlabel,ylabel);
  plotZRapEff.AddHist1D(hZRapBinEfficiency,"E");
  plotZRapEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZRapEff.SetLogx(0);
  plotZRapEff.SetLogy(0);
  plotZRapEff.SetYRange(0,1);
  plotZRapEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"|y^{#mu^{+}#mu^{-}}|");
  CPlot plotZRapPurity("ZRapBinPurity","",xlabel,ylabel);
  plotZRapPurity.AddHist1D(hZRapBinPurity,"E");
  plotZRapPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZRapPurity.SetLogx(0);
  plotZRapPurity.SetLogy(0);
  plotZRapPurity.SetYRange(0,1);
  plotZRapPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"|y^{#mu^{+}#mu^{-}}| (detector level)");
  sprintf(xlabel,"|y^{#mu^{+}#mu^{-}}| (hadron level)");
  CPlot plotZRapMatrix("ZRapMatrixNorm","",xlabel,ylabel);
  plotZRapMatrix.AddHist2D(hZRapMatrixNorm,"COLZ");
  plotZRapMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZRapMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotZRapMatrix.SetLogx(0);
  plotZRapMatrix.SetLogy(0);
  plotZRapMatrix.Draw(c,kTRUE,format);
  
  //
  // Lep1 Pt
  //
  
  TH2D *hLep1PtMatrixNorm  = new TH2D("hLep1PtMatrixNorm","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrixNorm->Sumw2();
  TH1D *hLep1PtBinEfficiency  = new TH1D("hLep1PtBinEfficiency","",25,Lep1PtBins); hLep1PtBinEfficiency->Sumw2();
  TH1D *hLep1PtBinPurity  = new TH1D("hLep1PtBinPurity","",25,Lep1PtBins); hLep1PtBinPurity->Sumw2();
 
  for(int i=0;i!=hLep1PtMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep1PtMatrix->GetNbinsY();++j)
	{
	  if(hLep1PtMatrix->GetBinContent(i+1,j+1)/hLep1PtTruth->GetBinContent(i+1)>0)
	    {
	      hLep1PtMatrixNorm->SetBinContent(i+1,j+1,hLep1PtMatrix->GetBinContent(i+1,j+1)/hLep1PtTruth->GetBinContent(i+1));
	    }
	  if(hLep1PtMatrix->GetBinContent(i+1,j+1)<0) hLep1PtMatrix->SetBinContent(i+1,j+1,0);
	  sum+= hLep1PtMatrix->GetBinContent(i+1,j+1);
	}
      hLep1PtBinEfficiency->SetBinContent(i+1,hLep1PtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hLep1PtMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep1PtMatrix->GetNbinsX();++j)
	{
	  sum+= hLep1PtMatrix->GetBinContent(j+1,i+1);
	}
      hLep1PtBinPurity->SetBinContent(i+1,hLep1PtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"p_{T} (leading muon) [GeV]");
  CPlot plotLep1PtEff("Lep1PtBinEfficiency","",xlabel,ylabel);
  plotLep1PtEff.AddHist1D(hLep1PtBinEfficiency,"E");
  plotLep1PtEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep1PtEff.SetLogx();
  plotLep1PtEff.SetLogy(0);
  plotLep1PtEff.SetYRange(0,1);
  plotLep1PtEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"p_{T} (leading muon) [GeV]");
  CPlot plotLep1PtPurity("Lep1PtBinPurity","",xlabel,ylabel);
  plotLep1PtPurity.AddHist1D(hLep1PtBinPurity,"E");
  plotLep1PtPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep1PtPurity.SetLogx();
  plotLep1PtPurity.SetLogy(0);
  plotLep1PtPurity.SetYRange(0,1);
  plotLep1PtPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"p_{T} (leading muon, detector level) [GeV]");
  sprintf(xlabel,"p_{T} (leading muon, hadron level) [GeV]");
  CPlot plotLep1PtMatrix("Lep1PtMatrixNorm","",xlabel,ylabel);
  plotLep1PtMatrix.AddHist2D(hLep1PtMatrixNorm,"COLZ");
  plotLep1PtMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep1PtMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotLep1PtMatrix.SetLogx();
  plotLep1PtMatrix.SetLogy();
  plotLep1PtMatrix.Draw(c,kTRUE,format);
  
  
  //
  // Lep2 Pt
  //
  
  TH2D *hLep2PtMatrixNorm  = new TH2D("hLep2PtMatrixNorm","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrixNorm->Sumw2();
  TH1D *hLep2PtBinEfficiency  = new TH1D("hLep2PtBinEfficiency","",20,Lep2PtBins); hLep2PtBinEfficiency->Sumw2();
  TH1D *hLep2PtBinPurity  = new TH1D("hLep2PtBinPurity","",20,Lep2PtBins); hLep2PtBinPurity->Sumw2();
  for(int i=0;i!=hLep2PtMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep2PtMatrix->GetNbinsY();++j)
	{
	  if(hLep2PtMatrix->GetBinContent(i+1,j+1)/hLep2PtTruth->GetBinContent(i+1)>0)
	    {
	      hLep2PtMatrixNorm->SetBinContent(i+1,j+1,hLep2PtMatrix->GetBinContent(i+1,j+1)/hLep2PtTruth->GetBinContent(i+1));
	    }
	  if(hLep2PtMatrix->GetBinContent(i+1,j+1)<0) hLep2PtMatrix->SetBinContent(i+1,j+1,0);
	  sum+= hLep2PtMatrix->GetBinContent(i+1,j+1);
	}
      hLep2PtBinEfficiency->SetBinContent(i+1,hLep2PtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hLep2PtMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep2PtMatrix->GetNbinsX();++j)
	{
	  sum+= hLep2PtMatrix->GetBinContent(j+1,i+1);
	}
      hLep2PtBinPurity->SetBinContent(i+1,hLep2PtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"p_{T} (2nd leading muon) [GeV]");
  CPlot plotLep2PtEff("Lep2PtBinEfficiency","",xlabel,ylabel);
  plotLep2PtEff.AddHist1D(hLep2PtBinEfficiency,"E");
  plotLep2PtEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep2PtEff.SetLogx();
  plotLep2PtEff.SetLogy(0);
  plotLep2PtEff.SetYRange(0,1);
  plotLep2PtEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"p_{T} (2nd leading muon) [GeV]");
  CPlot plotLep2PtPurity("Lep2PtBinPurity","",xlabel,ylabel);
  plotLep2PtPurity.AddHist1D(hLep2PtBinPurity,"E");
  plotLep2PtPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep2PtPurity.SetLogx();
  plotLep2PtPurity.SetLogy(0);
  plotLep2PtPurity.SetYRange(0,1);
  plotLep2PtPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"p_{T} (2nd leading muon, detector level) [GeV]");
  sprintf(xlabel,"p_{T} (2nd leading muon, hadron level) [GeV]");
  CPlot plotLep2PtMatrix("Lep2PtMatrixNorm","",xlabel,ylabel);
  plotLep2PtMatrix.AddHist2D(hLep2PtMatrixNorm,"COLZ");
  plotLep2PtMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep2PtMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotLep2PtMatrix.SetLogx();
  plotLep2PtMatrix.SetLogy();
  plotLep2PtMatrix.Draw(c,kTRUE,format);

  //
  // LepNeg Pt
  //
  
  TH2D *hLepNegPtMatrixNorm  = new TH2D("hLepNegPtMatrixNorm","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrixNorm->Sumw2();
  TH1D *hLepNegPtBinEfficiency  = new TH1D("hLepNegPtBinEfficiency","",25,LepNegPtBins); hLepNegPtBinEfficiency->Sumw2();
  TH1D *hLepNegPtBinPurity  = new TH1D("hLepNegPtBinPurity","",25,LepNegPtBins); hLepNegPtBinPurity->Sumw2();
 
  for(int i=0;i!=hLepNegPtMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hLepNegPtMatrix->GetNbinsY();++j)
	{
	  if(hLepNegPtMatrix->GetBinContent(i+1,j+1)/hLepNegPtTruth->GetBinContent(i+1)>0)
	    {
	      hLepNegPtMatrixNorm->SetBinContent(i+1,j+1,hLepNegPtMatrix->GetBinContent(i+1,j+1)/hLepNegPtTruth->GetBinContent(i+1));
	    }
	  if(hLepNegPtMatrix->GetBinContent(i+1,j+1)<0) hLepNegPtMatrix->SetBinContent(i+1,j+1,0);
	  sum+= hLepNegPtMatrix->GetBinContent(i+1,j+1);
	}
      hLepNegPtBinEfficiency->SetBinContent(i+1,hLepNegPtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hLepNegPtMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hLepNegPtMatrix->GetNbinsX();++j)
	{
	  sum+= hLepNegPtMatrix->GetBinContent(j+1,i+1);
	}
      hLepNegPtBinPurity->SetBinContent(i+1,hLepNegPtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"p_{T}^{#mu^{-}} [GeV]");
  CPlot plotLepNegPtEff("LepNegPtBinEfficiency","",xlabel,ylabel);
  plotLepNegPtEff.AddHist1D(hLepNegPtBinEfficiency,"E");
  plotLepNegPtEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLepNegPtEff.SetLogx();
  plotLepNegPtEff.SetLogy(0);
  plotLepNegPtEff.SetYRange(0,1);
  plotLepNegPtEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"p_{T}^{#mu^{-}} [GeV]");
  CPlot plotLepNegPtPurity("LepNegPtBinPurity","",xlabel,ylabel);
  plotLepNegPtPurity.AddHist1D(hLepNegPtBinPurity,"E");
  plotLepNegPtPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLepNegPtPurity.SetLogx();
  plotLepNegPtPurity.SetLogy(0);
  plotLepNegPtPurity.SetYRange(0,1);
  plotLepNegPtPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"p_{T}^{#mu^{-}} (detector level) [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{-}} (hadron level) [GeV]");
  CPlot plotLepNegPtMatrix("LepNegPtMatrixNorm","",xlabel,ylabel);
  plotLepNegPtMatrix.AddHist2D(hLepNegPtMatrixNorm,"COLZ");
  plotLepNegPtMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLepNegPtMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotLepNegPtMatrix.SetLogx();
  plotLepNegPtMatrix.SetLogy();
  plotLepNegPtMatrix.Draw(c,kTRUE,format);
 
  //
  // LepPos Pt
  //
  
  TH2D *hLepPosPtMatrixNorm  = new TH2D("hLepPosPtMatrixNorm","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrixNorm->Sumw2();
  TH1D *hLepPosPtBinEfficiency  = new TH1D("hLepPosPtBinEfficiency","",25,LepPosPtBins); hLepPosPtBinEfficiency->Sumw2();
  TH1D *hLepPosPtBinPurity  = new TH1D("hLepPosPtBinPurity","",25,LepPosPtBins); hLepPosPtBinPurity->Sumw2();
 
  for(int i=0;i!=hLepPosPtMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hLepPosPtMatrix->GetNbinsY();++j)
	{
	  if(hLepPosPtMatrix->GetBinContent(i+1,j+1)/hLepPosPtTruth->GetBinContent(i+1)>0)
	    {
	      hLepPosPtMatrixNorm->SetBinContent(i+1,j+1,hLepPosPtMatrix->GetBinContent(i+1,j+1)/hLepPosPtTruth->GetBinContent(i+1));
	    }
	  if(hLepPosPtMatrix->GetBinContent(i+1,j+1)<0) hLepPosPtMatrix->SetBinContent(i+1,j+1,0);
	  sum+= hLepPosPtMatrix->GetBinContent(i+1,j+1);
	}
      hLepPosPtBinEfficiency->SetBinContent(i+1,hLepPosPtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hLepPosPtMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hLepPosPtMatrix->GetNbinsX();++j)
	{
	  sum+= hLepPosPtMatrix->GetBinContent(j+1,i+1);
	}
      hLepPosPtBinPurity->SetBinContent(i+1,hLepPosPtMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"p_{T}^{#mu^{+}} [GeV]");
  CPlot plotLepPosPtEff("LepPosPtBinEfficiency","",xlabel,ylabel);
  plotLepPosPtEff.AddHist1D(hLepPosPtBinEfficiency,"E");
  plotLepPosPtEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLepPosPtEff.SetLogx();
  plotLepPosPtEff.SetLogy(0);
  plotLepPosPtEff.SetYRange(0,1);
  plotLepPosPtEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"p_{T}^{#mu^{+}} [GeV]");
  CPlot plotLepPosPtPurity("LepPosPtBinPurity","",xlabel,ylabel);
  plotLepPosPtPurity.AddHist1D(hLepPosPtBinPurity,"E");
  plotLepPosPtPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLepPosPtPurity.SetLogx();
  plotLepPosPtPurity.SetLogy(0);
  plotLepPosPtPurity.SetYRange(0,1);
  plotLepPosPtPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"p_{T}^{#mu^{+}} (detector level) [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{+}} (hadron level) [GeV]");
  CPlot plotLepPosPtMatrix("LepPosPtMatrixNorm","",xlabel,ylabel);
  plotLepPosPtMatrix.AddHist2D(hLepPosPtMatrixNorm,"COLZ");
  plotLepPosPtMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLepPosPtMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotLepPosPtMatrix.SetLogx();
  plotLepPosPtMatrix.SetLogy();
  plotLepPosPtMatrix.Draw(c,kTRUE,format);

  
  //
  // Lep1 Eta
  //

  TH2D *hLep1EtaMatrixNorm  = new TH2D("hLep1EtaMatrixNorm","",24,0,2.4,24,0,2.4); hLep1EtaMatrixNorm->Sumw2();
  TH1D *hLep1EtaBinEfficiency  = new TH1D("hLep1EtaBinEfficiency","",24,0,2.4); hLep1EtaBinEfficiency->Sumw2();
  TH1D *hLep1EtaBinPurity  = new TH1D("hLep1EtaBinPurity","",24,0,2.4); hLep1EtaBinPurity->Sumw2();
  for(int i=0;i!=hLep1EtaMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep1EtaMatrix->GetNbinsY();++j)
	{
	  if(hLep1EtaMatrix->GetBinContent(i+1,j+1)/hLep1EtaTruth->GetBinContent(i+1)>0)
	    {
	      hLep1EtaMatrixNorm->SetBinContent(i+1,j+1,hLep1EtaMatrix->GetBinContent(i+1,j+1)/hLep1EtaTruth->GetBinContent(i+1));
	    }
	  if(hLep1EtaMatrix->GetBinContent(i+1,j+1)<0) hLep1EtaMatrix->SetBinContent(i+1,j+1,0);
	  sum+= hLep1EtaMatrix->GetBinContent(i+1,j+1);
	}
      hLep1EtaBinEfficiency->SetBinContent(i+1,hLep1EtaMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hLep1EtaMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep1EtaMatrix->GetNbinsX();++j)
	{
	  sum+= hLep1EtaMatrix->GetBinContent(j+1,i+1);
	}
      hLep1EtaBinPurity->SetBinContent(i+1,hLep1EtaMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"|#eta| (leading muon)");
  CPlot plotLep1EtaEff("Lep1EtaBinEfficiency","",xlabel,ylabel);
  plotLep1EtaEff.AddHist1D(hLep1EtaBinEfficiency,"E");
  plotLep1EtaEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep1EtaEff.SetLogx(0);
  plotLep1EtaEff.SetLogy(0);
  plotLep1EtaEff.SetYRange(0,1);
  plotLep1EtaEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"|#eta| (leading muon)");
  CPlot plotLep1EtaPurity("Lep1EtaBinPurity","",xlabel,ylabel);
  plotLep1EtaPurity.AddHist1D(hLep1EtaBinPurity,"E");
  plotLep1EtaPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep1EtaPurity.SetLogx(0);
  plotLep1EtaPurity.SetLogy(0);
  plotLep1EtaPurity.SetYRange(0,1);
  plotLep1EtaPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"|#eta| (leading muon, detector level)");
  sprintf(xlabel,"|#eta| (leading muon, hadron level)");
  CPlot plotLep1EtaMatrix("Lep1EtaMatrixNorm","",xlabel,ylabel);
  plotLep1EtaMatrix.AddHist2D(hLep1EtaMatrixNorm,"COLZ");
  plotLep1EtaMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep1EtaMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotLep1EtaMatrix.SetLogx(0);
  plotLep1EtaMatrix.SetLogy(0);
  plotLep1EtaMatrix.Draw(c,kTRUE,format);
  

  //
  // Lep2 Eta
  //

  TH2D *hLep2EtaMatrixNorm  = new TH2D("hLep2EtaMatrixNorm","",24,0,2.4,24,0,2.4); hLep2EtaMatrixNorm->Sumw2();
  TH1D *hLep2EtaBinEfficiency  = new TH1D("hLep2EtaBinEfficiency","",24,0,2.4); hLep2EtaBinEfficiency->Sumw2();
  TH1D *hLep2EtaBinPurity  = new TH1D("hLep2EtaBinPurity","",24,0,2.4); hLep2EtaBinPurity->Sumw2();
  for(int i=0;i!=hLep2EtaMatrix->GetNbinsX();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep2EtaMatrix->GetNbinsY();++j)
	{
	  if(hLep2EtaMatrix->GetBinContent(i+1,j+1)/hLep2EtaTruth->GetBinContent(i+1)>0)
	    {
	      hLep2EtaMatrixNorm->SetBinContent(i+1,j+1,hLep2EtaMatrix->GetBinContent(i+1,j+1)/hLep2EtaTruth->GetBinContent(i+1));
	    }
	  if(hLep2EtaMatrix->GetBinContent(i+1,j+1)<0) hLep2EtaMatrix->SetBinContent(i+1,j+1,0);
	  sum+= hLep2EtaMatrix->GetBinContent(i+1,j+1);
	}
      hLep2EtaBinEfficiency->SetBinContent(i+1,hLep2EtaMatrix->GetBinContent(i+1,i+1)/sum);
    }

  for(int i=0;i!=hLep2EtaMatrix->GetNbinsY();++i)
    {
      double sum=0;
      for(int j=0;j!=hLep2EtaMatrix->GetNbinsX();++j)
	{
	  sum+= hLep2EtaMatrix->GetBinContent(j+1,i+1);
	}
      hLep2EtaBinPurity->SetBinContent(i+1,hLep2EtaMatrix->GetBinContent(i+1,i+1)/sum);
    }

  sprintf(ylabel,"Efficiency [%%]");
  sprintf(xlabel,"|#eta| (2nd leading muon)");
  CPlot plotLep2EtaEff("Lep2EtaBinEfficiency","",xlabel,ylabel);
  plotLep2EtaEff.AddHist1D(hLep2EtaBinEfficiency,"E");
  plotLep2EtaEff.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep2EtaEff.SetLogx(0);
  plotLep2EtaEff.SetLogy(0);
  plotLep2EtaEff.SetYRange(0,1);
  plotLep2EtaEff.Draw(c,kTRUE,format);

  sprintf(ylabel,"Purity [%%]");
  sprintf(xlabel,"|#eta| (2nd leading muon)");
  CPlot plotLep2EtaPurity("Lep2EtaBinPurity","",xlabel,ylabel);
  plotLep2EtaPurity.AddHist1D(hLep2EtaBinPurity,"E");
  plotLep2EtaPurity.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep2EtaPurity.SetLogx(0);
  plotLep2EtaPurity.SetLogy(0);
  plotLep2EtaPurity.SetYRange(0,1);
  plotLep2EtaPurity.Draw(c,kTRUE,format);

  sprintf(ylabel,"|#eta| (2nd leading muon, detector level)");
  sprintf(xlabel,"|#eta| (2nd leading muon, hadron level)");
  CPlot plotLep2EtaMatrix("Lep2EtaMatrixNorm","",xlabel,ylabel);
  plotLep2EtaMatrix.AddHist2D(hLep2EtaMatrixNorm,"COLZ");
  plotLep2EtaMatrix.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotLep2EtaMatrix.AddTextBox(plabel,0.32,0.90,0.95,0.98,0);
  plotLep2EtaMatrix.SetLogx(0);
  plotLep2EtaMatrix.SetLogy(0);
  plotLep2EtaMatrix.Draw(c,kTRUE,format);

  c->Delete();
  
  outFile->cd();
  outFile->Write();
  outFile->Close();
  }

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

  gBenchmark->Show("plotZmmGen"); 
}
