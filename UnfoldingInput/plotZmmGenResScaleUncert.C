//================================================================================================
//
// Make plots of various distributions after Zee selection
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

#endif

//=== MAIN MACRO ================================================================================================= 

void plotZmmGenResScaleUncert(const TString  inputDir,        // input directory
			      const TString  outputDir,       // output directory
			      const Double_t lumi             // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmmGenResScaleUncert");

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
   
   //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
 
 
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
    
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

    
  // Read input file and get the TTrees
  TString infilename = inputDir + TString("/") + TString("zmm_select.raw.root");
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
  TString outfilename = outputDir + TString("/") + TString("zmm_UnfoldInputs_ResScale.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");

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

  TH1D *hPhiStarReco  = new TH1D("hPhiStarReco","",27,PhiStarBins); hPhiStarReco->Sumw2();
  TH1D *hPhiStarTruth  = new TH1D("hPhiStarTruth","",27,PhiStarBins); hPhiStarTruth->Sumw2();
  TH2D *hPhiStarMatrix  = new TH2D("hPhiStarMatrix","",27,PhiStarBins,27,PhiStarBins); hPhiStarMatrix->Sumw2();
  

  TH1D *hZRapReco  = new TH1D("hZRapReco","",24,0,2.4); hZRapReco->Sumw2();
  TH1D *hZRapTruth  = new TH1D("hZRapTruth","",24,0,2.4); hZRapTruth->Sumw2();
  TH2D *hZRapMatrix  = new TH2D("hZRapMatrix","",24,0,2.4,24,0,2.4); hZRapMatrix->Sumw2();
  
  
  TH1D *hLep1PtReco  = new TH1D("hLep1PtReco","",25,Lep1PtBins); hLep1PtReco->Sumw2();
  TH1D *hLep1PtTruth  = new TH1D("hLep1PtTruth","",25,Lep1PtBins); hLep1PtTruth->Sumw2();
  TH2D *hLep1PtMatrix  = new TH2D("hLep1PtMatrix","",25,Lep1PtBins,25,Lep1PtBins); hLep1PtMatrix->Sumw2();
  
  TH1D *hLep2PtReco  = new TH1D("hLep2PtReco","",20,Lep2PtBins); hLep2PtReco->Sumw2();
  TH1D *hLep2PtTruth  = new TH1D("hLep2PtTruth","",20,Lep2PtBins); hLep2PtTruth->Sumw2();
  TH2D *hLep2PtMatrix  = new TH2D("hLep2PtMatrix","",20,Lep2PtBins,20,Lep2PtBins); hLep2PtMatrix->Sumw2();

  TH1D *hLepNegPtReco  = new TH1D("hLepNegPtReco","",25,LepNegPtBins); hLepNegPtReco->Sumw2();
  TH1D *hLepNegPtTruth  = new TH1D("hLepNegPtTruth","",25,LepNegPtBins); hLepNegPtTruth->Sumw2();
  TH2D *hLepNegPtMatrix  = new TH2D("hLepNegPtMatrix","",25,LepNegPtBins,25,LepNegPtBins); hLepNegPtMatrix->Sumw2();

  TH1D *hLepPosPtReco  = new TH1D("hLepPosPtReco","",25,LepPosPtBins); hLepPosPtReco->Sumw2();
  TH1D *hLepPosPtTruth  = new TH1D("hLepPosPtTruth","",25,LepPosPtBins); hLepPosPtTruth->Sumw2();
  TH2D *hLepPosPtMatrix  = new TH2D("hLepPosPtMatrix","",25,LepPosPtBins,25,LepPosPtBins); hLepPosPtMatrix->Sumw2();
 
  TH1D *hLep1EtaReco  = new TH1D("hLep1EtaReco","",24,0,2.4); hLep1EtaReco->Sumw2();
  TH1D *hLep1EtaTruth  = new TH1D("hLep1EtaTruth","",24,0,2.4); hLep1EtaTruth->Sumw2();
  TH2D *hLep1EtaMatrix  = new TH2D("hLep1EtaMatrix","",24,0,2.4,24,0,2.4); hLep1EtaMatrix->Sumw2();
  

  TH1D *hLep2EtaReco  = new TH1D("hLep2EtaReco","",24,0,2.4); hLep2EtaReco->Sumw2();
  TH1D *hLep2EtaTruth  = new TH1D("hLep2EtaTruth","",24,0,2.4); hLep2EtaTruth->Sumw2();
  TH2D *hLep2EtaMatrix  = new TH2D("hLep2EtaMatrix","",24,0,2.4,24,0,2.4); hLep2EtaMatrix->Sumw2();
  
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

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

   
    if(triggerDec&&goodPV&&matchTrigger&&nlep>=2&&q1!=q2&&dilep->M()>MASS_LOW&&dilep->M()<MASS_HIGH&&l1.Pt()>=PT_CUT&&l2.Pt()>=PT_CUT&&fabs(l1.Eta())<=ETA_CUT&&fabs(l2.Eta())<=ETA_CUT)
      {
	isReco=true;
      }
    if(ngenlep>=2&&gendilep->M()>MASS_LOW&&gendilep->M()<MASS_HIGH&&genlep1->Pt()>=PT_CUT&&genlep2->Pt()>=PT_CUT&&fabs(genlep1->Eta())<=ETA_CUT&&fabs(genlep2->Eta())<=ETA_CUT)
      {
	isGen=true;
      }

    Double_t genweight = 1;
    genweight *= scale1fbGen*lumi;
    Double_t weight = 1;
    weight *=scale1fb*lumi;

    if(isReco)
      {
	hMassMC ->Fill(dilep->M(),weight*corr);
	hZPtReco ->Fill(dilep->Pt(),weight*corr);
	hPhiStarReco ->Fill(phistar,weight*corr);
	hZRapReco ->Fill(fabs(dilep->Rapidity()),weight*corr);
	hLep1PtReco ->Fill(l1.Pt(),weight*corr);
	hLep2PtReco ->Fill(l2.Pt(),weight*corr);
	if(lq1<0)
	  {
	    hLepNegPtReco ->Fill(l1.Pt(),weight*corr);
	    hLepPosPtReco ->Fill(l2.Pt(),weight*corr);
	  }
	else
	  {
	    hLepNegPtReco ->Fill(l2.Pt(),weight*corr);
	    hLepPosPtReco ->Fill(l1.Pt(),weight*corr);
	  }
	hLep1EtaReco ->Fill(fabs(l1.Eta()),weight*corr);
	hLep2EtaReco ->Fill(fabs(l2.Eta()),weight*corr);
      }
    if(isGen)
      {
	hZPtTruth ->Fill(gendilep->Pt(),genweight);
	hPhiStarTruth ->Fill(genphistar,genweight);
	hZRapTruth ->Fill(fabs(gendilep->Rapidity()),genweight);
	hLep1PtTruth ->Fill(genlep1->Pt(),genweight);
	hLep2PtTruth ->Fill(genlep2->Pt(),genweight);
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
	hLep1EtaTruth ->Fill(fabs(genlep1->Eta()),genweight);
	hLep2EtaTruth ->Fill(fabs(genlep2->Eta()),genweight);
      }
    if(isReco&&isGen)
      {
	hZPtMatrix ->Fill(gendilep->Pt(),dilep->Pt(),weight*corr);
	hPhiStarMatrix ->Fill(genphistar,phistar,weight*corr);
	hZRapMatrix ->Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corr);
	hLep1PtMatrix ->Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	hLep2PtMatrix ->Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	if(lq1<0&&genq1<0)
	  {
	    hLepNegPtMatrix ->Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	  }
	else if(lq1<0&&genq1>0)
	  {
	    hLepNegPtMatrix ->Fill(genlep2->Pt(),l1.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep1->Pt(),l2.Pt(),weight*corr);
	  }
	else if(lq1>0&&genq1<0)
	  {
	    hLepNegPtMatrix ->Fill(genlep1->Pt(),l2.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep2->Pt(),l1.Pt(),weight*corr);
	  }
	else if(lq1>0&&genq1>0)
	  {
	    hLepNegPtMatrix ->Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	    hLepPosPtMatrix ->Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	  }
	hLep1EtaMatrix ->Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corr);
	hLep2EtaMatrix ->Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corr);
      }
    delete gendilep;
    delete dilep;
  }
  delete infile;
  infile=0, intree=0; 
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  outFile->cd();
  outFile->Write();
  outFile->Close();
     
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl; 

  gBenchmark->Show("plotZmmGenResScaleUncert"); 
}
