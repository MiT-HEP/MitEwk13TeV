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

//=== FUNCTION DECLARATIONS ======================================================================================

void create( vector<TH1D> &v , const string name, int nbins, double xlow, double xhigh , const int ntoys) {for(int itoys=0;itoys!=ntoys;++itoys) {v.push_back(TH1D((name+Form("_%d",itoys)).c_str(),"",nbins,xlow,xhigh)); v[itoys].Sumw2();}}
void create( vector<TH1D> &v , const string name, int nbins, double* x , const int ntoys) {for(int itoys=0;itoys!=ntoys;++itoys) {v.push_back(TH1D((name+Form("_%d",itoys)).c_str(),"",nbins,x)); v[itoys].Sumw2();}}
void create( vector<TH2D> &v , const string name, int nbins, double xlow, double xhigh , const int ntoys) {for(int itoys=0;itoys!=ntoys;++itoys) {v.push_back(TH2D((name+Form("_%d",itoys)).c_str(),"",nbins,xlow,xhigh,nbins,xlow,xhigh)); v[itoys].Sumw2();}}
void create( vector<TH2D> &v , const string name, int nbins, double* x , const int ntoys) {for(int itoys=0;itoys!=ntoys;++itoys) {v.push_back(TH2D((name+Form("_%d",itoys)).c_str(),"",nbins,x,nbins,x)); v[itoys].Sumw2();}}

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

  //----- 
  const int NTOYS = 100;

  // efficiency files
  const TString baseDir = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/"; 
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

   
   //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
 
 
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
    
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  triggerDec;
  UInt_t  goodPV;
  UInt_t  matchTrigger;
  UInt_t  ngenlep;
  TLorentzVector *genlep1=0, *genlep2=0;
  Int_t   genq1, genq2;
  UInt_t nlep;
  TLorentzVector *lep1=0, *lep2=0;
  Int_t   q1, q2;
  Float_t scale1fbGen,scale1fb;

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

    
  // Read input file and get the TTrees
  TString infilename = inputDir + TString("/") + TString("zmm_select.raw.root");
  cout << "Processing " << infilename << "..." << endl;
  infile = TFile::Open(infilename);assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);

  intree -> SetBranchStatus("*",0);
  intree -> SetBranchStatus("runNum",1);
  intree -> SetBranchStatus("lumiSec",1);
  intree -> SetBranchStatus("evtNum",1);
  intree -> SetBranchStatus("triggerDec",1);
  intree -> SetBranchStatus("goodPV",1);
  intree -> SetBranchStatus("matchTrigger",1);
  intree -> SetBranchStatus("ngenlep",1);
  intree -> SetBranchStatus("genlep1",1);
  intree -> SetBranchStatus("genlep2",1);
  intree -> SetBranchStatus("genq1",1);
  intree -> SetBranchStatus("genq2",1);
  intree -> SetBranchStatus("nlep",1);
  intree -> SetBranchStatus("lep1",1);
  intree -> SetBranchStatus("lep2",1);
  intree -> SetBranchStatus("q1",1);
  intree -> SetBranchStatus("q2",1);
  intree -> SetBranchStatus("scale1fbGen",1);
  intree -> SetBranchStatus("scale1fb",1);
  
  
  intree->SetBranchAddress("runNum",   &runNum);     // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);     // event number
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
   
  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("zmm_UnfoldInputs_ResScale.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");

  //
  // Create histograms
  //
  double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
    double PhiStarBins[]={0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.01,0.012,0.014,0.016,0.018,0.021,0.024,0.027,0.030,0.034,0.038,0.044,0.050,0.058,0.066,0.076,0.088,0.10,0.12,0.14,0.16,0.18,0.20,0.24,0.28,0.34,0.42,0.52,0.64,0.8,1.0,1.5,2,3};
    double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
    double Lep2PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
    double LepNegPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
    double LepPosPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};

    const int nBinsZPt= sizeof(ZPtBins)/sizeof(double)-1;
  vector<TH1D> hZPtReco; create(hZPtReco,"hZPtReco",nBinsZPt,ZPtBins,NTOYS);
  vector<TH1D> hZPtTruth; create(hZPtTruth,"hZPtTruth",nBinsZPt,ZPtBins,NTOYS);
  vector<TH2D> hZPtMatrix; create(hZPtMatrix,"hZPtMatrix",nBinsZPt,ZPtBins,NTOYS);

  const int nBinsPhiStar= sizeof(PhiStarBins)/sizeof(double)-1;
  vector<TH1D> hPhiStarReco; create(hPhiStarReco,"hPhiStarReco",nBinsPhiStar,PhiStarBins,NTOYS);
  vector<TH1D> hPhiStarTruth; create(hPhiStarTruth,"hPhiStarTruth",nBinsPhiStar,PhiStarBins,NTOYS);
  vector<TH2D> hPhiStarMatrix; create(hPhiStarMatrix,"hPhiStarMatrix",nBinsPhiStar,PhiStarBins,NTOYS);
  
  vector<TH1D> hZRapReco; create(hZRapReco,"hZRapReco",24,0,2.4,NTOYS);
  vector<TH1D> hZRapTruth; create(hZRapTruth,"hZRapTruth",24,0,2.4,NTOYS); 
  vector<TH2D> hZRapMatrix; create(hZRapMatrix,"hZRapMatrix",24,0,2.4,NTOYS);
  
  const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
  vector<TH1D> hLep1PtReco; create(hLep1PtReco,"hLep1PtReco",nBinsLep1Pt,Lep1PtBins,NTOYS);
  vector<TH1D> hLep1PtTruth; create(hLep1PtTruth,"hLep1PtTruth",nBinsLep1Pt,Lep1PtBins,NTOYS); 
  vector<TH2D> hLep1PtMatrix; create(hLep1PtMatrix,"hLep1PtMatrix",nBinsLep1Pt,Lep1PtBins,NTOYS);
  
  const int nBinsLep2Pt= sizeof(Lep2PtBins)/sizeof(double)-1;
  vector<TH1D> hLep2PtReco; create(hLep2PtReco,"hLep2PtReco",nBinsLep2Pt,Lep2PtBins,NTOYS);
  vector<TH1D> hLep2PtTruth; create(hLep2PtTruth,"hLep2PtTruth",nBinsLep2Pt,Lep2PtBins,NTOYS);
  vector<TH2D> hLep2PtMatrix; create(hLep2PtMatrix,"hLep2PtMatrix",nBinsLep2Pt,Lep2PtBins,NTOYS);

  const int nBinsLepNegPt= sizeof(LepNegPtBins)/sizeof(double)-1;
  vector<TH1D> hLepNegPtReco; create(hLepNegPtReco,"hLepNegPtReco",nBinsLepNegPt,LepNegPtBins,NTOYS);
  vector<TH1D> hLepNegPtTruth; create(hLepNegPtTruth,"hLepNegPtTruth",nBinsLepNegPt,LepNegPtBins,NTOYS);
  vector<TH2D> hLepNegPtMatrix; create(hLepNegPtMatrix,"hLepNegPtMatrix",nBinsLepNegPt,LepNegPtBins,NTOYS);

  const int nBinsLepPosPt= sizeof(LepPosPtBins)/sizeof(double)-1;
  vector<TH1D> hLepPosPtReco; create(hLepPosPtReco,"hLepPosPtReco",nBinsLepPosPt,LepPosPtBins,NTOYS);
  vector<TH1D> hLepPosPtTruth; create(hLepPosPtTruth,"hLepPosPtTruth",nBinsLepPosPt,LepPosPtBins,NTOYS);
  vector<TH2D> hLepPosPtMatrix; create(hLepPosPtMatrix,"hLepPosPtMatrix",nBinsLepPosPt,LepPosPtBins,NTOYS);
 

  vector<TH1D> hLep1EtaReco; create(hLep1EtaReco,"hLep1EtaReco",24,0,2.4,NTOYS);
  vector<TH1D> hLep1EtaTruth; create(hLep1EtaTruth,"hLep1EtaTruth",24,0,2.4,NTOYS);
  vector<TH2D> hLep1EtaMatrix; create(hLep1EtaMatrix,"hLep1EtaMatrix",24,0,2.4,NTOYS);
  

  vector<TH1D> hLep2EtaReco; create(hLep2EtaReco,"hLep2EtaReco",24,0,2.4,NTOYS);
  vector<TH1D> hLep2EtaTruth; create(hLep2EtaTruth,"hLep2EtaTruth",24,0,2.4,NTOYS);
  vector<TH2D> hLep2EtaMatrix; create(hLep2EtaMatrix,"hLep2EtaMatrix",24,0,2.4,NTOYS);
  
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    Double_t genweight = 1;
    genweight *= scale1fbGen*lumi;
    Double_t weight = 1;
    weight *=scale1fb*lumi;	

    for(int itoys=0;itoys!=NTOYS;++itoys)
      {
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
	
	if(isReco)
	  {
	    hZPtReco[itoys].Fill(dilep->Pt(),weight*corr);
	    hPhiStarReco[itoys].Fill(phistar,weight*corr);
	    hZRapReco[itoys].Fill(fabs(dilep->Rapidity()),weight*corr);
	    hLep1PtReco[itoys].Fill(l1.Pt(),weight*corr);
	    hLep2PtReco[itoys].Fill(l2.Pt(),weight*corr);
	    if(lq1<0)
	      {
		hLepNegPtReco[itoys].Fill(l1.Pt(),weight*corr);
		hLepPosPtReco[itoys].Fill(l2.Pt(),weight*corr);
	      }
	    else
	      {
		hLepNegPtReco[itoys].Fill(l2.Pt(),weight*corr);
		hLepPosPtReco[itoys].Fill(l1.Pt(),weight*corr);
	      }
	    hLep1EtaReco[itoys].Fill(fabs(l1.Eta()),weight*corr);
	    hLep2EtaReco[itoys].Fill(fabs(l2.Eta()),weight*corr);
	  }
	if(isGen)
	  {
	    hZPtTruth[itoys].Fill(gendilep->Pt(),genweight);
	    hPhiStarTruth[itoys].Fill(genphistar,genweight);
	    hZRapTruth[itoys].Fill(fabs(gendilep->Rapidity()),genweight);
	    hLep1PtTruth[itoys].Fill(genlep1->Pt(),genweight);
	    hLep2PtTruth[itoys].Fill(genlep2->Pt(),genweight);
	    if(genq1<0)
	      {
		hLepNegPtTruth[itoys].Fill(genlep1->Pt(),genweight);
		hLepPosPtTruth[itoys].Fill(genlep2->Pt(),genweight);
	      }
	    else
	      {
		hLepNegPtTruth[itoys].Fill(genlep2->Pt(),genweight);
		hLepPosPtTruth[itoys].Fill(genlep1->Pt(),genweight);
	      }
	    hLep1EtaTruth[itoys].Fill(fabs(genlep1->Eta()),genweight);
	    hLep2EtaTruth[itoys].Fill(fabs(genlep2->Eta()),genweight);
	  }
	if(isReco&&isGen)
	  {
	    hZPtMatrix[itoys].Fill(gendilep->Pt(),dilep->Pt(),weight*corr);
	    hPhiStarMatrix[itoys].Fill(genphistar,phistar,weight*corr);
	    hZRapMatrix[itoys].Fill(fabs(gendilep->Rapidity()),fabs(dilep->Rapidity()),weight*corr);
	    hLep1PtMatrix[itoys].Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	    hLep2PtMatrix[itoys].Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	    if(lq1<0&&genq1<0)
	      {
		hLepNegPtMatrix[itoys].Fill(genlep1->Pt(),l1.Pt(),weight*corr);
		hLepPosPtMatrix[itoys].Fill(genlep2->Pt(),l2.Pt(),weight*corr);
	      }
	    else if(lq1<0&&genq1>0)
	      {
		hLepNegPtMatrix[itoys].Fill(genlep2->Pt(),l1.Pt(),weight*corr);
		hLepPosPtMatrix[itoys].Fill(genlep1->Pt(),l2.Pt(),weight*corr);
	      }
	    else if(lq1>0&&genq1<0)
	      {
		hLepNegPtMatrix[itoys].Fill(genlep1->Pt(),l2.Pt(),weight*corr);
		hLepPosPtMatrix[itoys].Fill(genlep2->Pt(),l1.Pt(),weight*corr);
	      }
	    else if(lq1>0&&genq1>0)
	      {
		hLepNegPtMatrix[itoys].Fill(genlep2->Pt(),l2.Pt(),weight*corr);
		hLepPosPtMatrix[itoys].Fill(genlep1->Pt(),l1.Pt(),weight*corr);
	      }
	    hLep1EtaMatrix[itoys].Fill(fabs(genlep1->Eta()),fabs(l1.Eta()),weight*corr);
	    hLep2EtaMatrix[itoys].Fill(fabs(genlep2->Eta()),fabs(l2.Eta()),weight*corr);
	  }
	delete gendilep;
	delete dilep;
      }
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
