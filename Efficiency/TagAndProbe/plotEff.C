//================================================================================================
//
// Signal Extraction
//-------------------
//  0: probe counting
//  1: Breit-Wigner convolved with Crystal Ball function
//  2: MC template convolved with Gaussian
//  3: Phil's Crystal Ball based "Voigtian" shape
//  4: Unbinned MC data convolved with Gaussian
//
// Background Model
//------------------
//  0: no background
//  1: exponential model
//  2: erfc*exp model
//  3: double exponential model
//  4: linear*exp model
//  5: quadratic*exp model
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <cmath>
#include <math.h>

#include "../../Utils/CPlot.hh"	            // helper class for plots
#include "../../Utils/MitStyleRemix.hh"         // style settings for drawing
#include "../../Utils/CEffUser1D.hh"            // class for handling efficiency graphs
#include "../../Utils/CEffUser2D.hh"            // class for handling efficiency tables

// structure for output ntuple
#include "EffData.hh"

#include "ZSignals.hh"
#include "ZBackgrounds.hh"
#endif

// RooFit headers
#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"

// bin size constants
#define BIN_SIZE_PASS 1
#define BIN_SIZE_FAIL 1

//=== FUNCTION DECLARATIONS ======================================================================================


// generate web page
void makeHTML(const TString outDir);
void makeHTML(const TString outDir, const TString name, const Int_t n);

// Make efficiency graph
TGraphAsymmErrors* makeEffGraph(const vector<Double_t> &edgesv, const vector<TTree*> &passv, const vector<TTree*> &failv, const Int_t method,
                                const TString name, const Double_t massLo, const Double_t massHi, 
				const TString format, const Bool_t doAbsEta,const double lumi);
TGraphAsymmErrors* makeEffGraph(const vector<Double_t> &edgesv, const vector<TTree*> &passv, const vector<TTree*> &failv,
                                const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 		                
				const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi, 
				const TString format, const Bool_t doAbsEta,const double lumi, const TString yaxislabel, const int charge);

// Make 2D efficiency map
void makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, const vector<TTree*> &passv, const vector<TTree*> &failv, const Int_t method,
                   const TString name, const Double_t massLo, const Double_t massHi,
		   const TString format, const Bool_t doAbsEta,const double lumi);
void makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, const vector<TTree*> &passv, const vector<TTree*> &failv,
                   const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		   const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi, 
		   const TString format, const Bool_t doAbsEta,const double lumi,const TString yaxislabel,const int charge);


// Generate MC-based signal templates
void generateHistTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const TH1D* puWeights, const int typeCode); 
void generateDataTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge); 
			   
// Perform count
void performCount(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                  const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		  TTree *passTree, TTree *failTree, const Int_t method, 
		  const TString name, const Double_t massLo, const Double_t massHi, const TString format, const Bool_t doAbsEta,
		  TCanvas *cpass, TCanvas *cfail,const double lumi);

void performFitBkgOnly(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		TTree *passTree, TTree *failTree,
		const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
		const TString format, const Bool_t doAbsEta, TCanvas *cpass, TCanvas *cfail,const double lumi,const TString yaxislabel,const int charge);

// Perform fit
void performFit(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		TTree *passTree, TTree *failTree,
		const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
		const TString format, const Bool_t doAbsEta, TCanvas *cpass, TCanvas *cfail,const double lumi,const TString yaxislabel,const int charge);

// Print correlations
void printCorrelations(ostream& os, RooFitResult *res);

// Parse fit results file
void parseFitResults(ifstream &ifs, double &eff, double &errl, double &errh);


//=== MAIN MACRO ================================================================================================= 

void plotEff(const TString conf,            // input binning file
             const Int_t   sigModPass,      // signal extraction method for PASS sample
	     const Int_t   bkgModPass,      // background model for PASS sample
	     const Int_t   sigModFail,      // signal extraction method for FAIL sample	     
	     const Int_t   bkgModFail,      // background model for FAIL sample
	     const TString infilename,      // ROOT file of probes
             const TString outputDir,       // output directory
             const TString format,          // plot format
	     const Bool_t  doAbsEta,        // bin in |eta| instead of eta?
	     const Int_t   doPU,            // PU re-weighting mode
	     const Int_t   charge,          // 0 (no charge requirement), -1, +1
             const TString xaxislabel="",   // 'Supercluster' or 'Muon'
             const TString yaxislabel="",   // use different labels for each efficiency
             const double ylow=0.5,         // lower limit of eff axis
             const double yhigh=1.02,       // upper limit of eff axis
             const double lumi=40.0,         // luminosity for plot label
	     const TString mcfilename="",   // ROOT file containing MC events to generate templates from
	     const UInt_t  runNumLo=0,      // lower bound of run range
	     const UInt_t  runNumHi=999999  // upper bound of run range
) {
  gBenchmark->Start("plotEff");
  // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // signal extraction mass region
  const Double_t massLo    = 60;
  const Double_t massHi    = 120;
  
  // fit mass region
  const Double_t fitMassLo = massLo;
  const Double_t fitMassHi = massHi;

  // Weights for PU-reweighting
  //const TString pufnameA   ("/data/blue/ksung/TagAndProbeExample/PileupReweighting.Summer11DYmm_To_Run2011A.root");
  //const TString pufnameB   ("/data/blue/ksung/TagAndProbeExample/PileupReweighting.Summer11DYmm_To_Run2011B.root");
  //const TString pufnameFull("/data/blue/ksung/TagAndProbeExample/PileupReweighting.Summer11DYmm_To_Full2011.root");
  
  // efficiency error calculation method
  // method: 0 -> Clopper-Pearson
  //         1 -> Feldman-Cousins
  const Int_t method=0;
  
  // bin edges for kinematic variables
  vector<Double_t> ptBinEdgesv;
  vector<Double_t> etaBinEdgesv;
  vector<Double_t> phiBinEdgesv;
  vector<Double_t> npvBinEdgesv;
  
  //
  // parse binning file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  Int_t opts[6];
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    Double_t edge;
    stringstream ss(line);
    if(state==0) {
      ss >> opts[0] >> opts[1] >> opts[2] >> opts[3] >> opts[4] >> opts[5];
    } else {
      ss >> edge;
      if(state==1)      { ptBinEdgesv.push_back(edge);  }
      else if(state==2) { etaBinEdgesv.push_back(edge); }
      else if(state==3) { phiBinEdgesv.push_back(edge); }
      else if(state==4) { npvBinEdgesv.push_back(edge); }
    }
  }
  ifs.close();
  
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir + TString("/plots");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  //  
  // Set up binning in kinematic variables
  //        
  const UInt_t ptNbins = ptBinEdgesv.size()-1;
  Double_t ptEdges[ptBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<ptBinEdgesv.size(); iedge++)
    ptEdges[iedge] = ptBinEdgesv[iedge];
    
  const UInt_t etaNbins = etaBinEdgesv.size()-1;
  Double_t etaEdges[etaBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<etaBinEdgesv.size(); iedge++)
    etaEdges[iedge] = etaBinEdgesv[iedge];
  
  const UInt_t phiNbins = phiBinEdgesv.size()-1;
  Double_t phiEdges[phiBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<phiBinEdgesv.size(); iedge++)
    phiEdges[iedge] = phiBinEdgesv[iedge];
  
  const UInt_t npvNbins = npvBinEdgesv.size()-1;
  Double_t npvEdges[npvBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<npvBinEdgesv.size(); iedge++)
    npvEdges[iedge] = npvBinEdgesv[iedge];
  
  char tname[50];
  Float_t mass;
  Double_t wgt;
  
  vector<TTree*> passTreePtv;
  vector<TTree*> failTreePtv;
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    sprintf(tname,"passPt_%i",ibin);
    passTreePtv.push_back(new TTree(tname,""));
    passTreePtv[ibin]->Branch("m",&mass,"m/F");
    passTreePtv[ibin]->Branch("w",&wgt,"w/D");
    passTreePtv[ibin]->SetDirectory(0);
    sprintf(tname,"failPt_%i",ibin);
    failTreePtv.push_back(new TTree(tname,""));
    failTreePtv[ibin]->Branch("m",&mass,"m/F");
    failTreePtv[ibin]->Branch("w",&wgt,"w/D");
    failTreePtv[ibin]->SetDirectory(0);
  }
  
  vector<TTree*> passTreeEtav;
  vector<TTree*> failTreeEtav;
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
    sprintf(tname,"passEta_%i",ibin);
    passTreeEtav.push_back(new TTree(tname,""));
    passTreeEtav[ibin]->Branch("m",&mass,"m/F");
    passTreeEtav[ibin]->Branch("w",&wgt,"w/D");
    passTreeEtav[ibin]->SetDirectory(0);
    sprintf(tname,"failEta_%i",ibin);
    failTreeEtav.push_back(new TTree(tname,""));
    failTreeEtav[ibin]->Branch("m",&mass,"m/F");
    failTreeEtav[ibin]->Branch("w",&wgt,"w/D");
    failTreeEtav[ibin]->SetDirectory(0);
  }
  
  vector<TTree*> passTreePhiv;  
  vector<TTree*> failTreePhiv;
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    sprintf(tname,"passPhi_%i",ibin);
    passTreePhiv.push_back(new TTree(tname,""));
    passTreePhiv[ibin]->Branch("m",&mass,"m/F");
    passTreePhiv[ibin]->Branch("w",&wgt,"w/D");
    passTreePhiv[ibin]->SetDirectory(0);
    sprintf(tname,"failPhi_%i",ibin);
    failTreePhiv.push_back(new TTree(tname,""));
    failTreePhiv[ibin]->Branch("m",&mass,"m/F");
    failTreePhiv[ibin]->Branch("w",&wgt,"w/D");
    failTreePhiv[ibin]->SetDirectory(0);
  }
  
  vector<TTree*> passTreeEtaPtv;  
  vector<TTree*> failTreeEtaPtv;
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    sprintf(tname,"passEtaPt_%i",ibin);
    passTreeEtaPtv.push_back(new TTree(tname,""));
    passTreeEtaPtv[ibin]->Branch("m",&mass,"m/F");
    passTreeEtaPtv[ibin]->Branch("w",&wgt,"w/D");
    passTreeEtaPtv[ibin]->SetDirectory(0);
    sprintf(tname,"failEtaPt_%i",ibin);
    failTreeEtaPtv.push_back(new TTree(tname,""));
    failTreeEtaPtv[ibin]->Branch("m",&mass,"m/F");
    failTreeEtaPtv[ibin]->Branch("w",&wgt,"w/D");
    failTreeEtaPtv[ibin]->SetDirectory(0);
  }
  
  vector<TTree*> passTreeEtaPhiv;
  vector<TTree*> failTreeEtaPhiv;
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    sprintf(tname,"passEtaPhi_%i",ibin); 
    passTreeEtaPhiv.push_back(new TTree(tname,""));
    passTreeEtaPhiv[ibin]->Branch("m",&mass,"m/F");
    passTreeEtaPhiv[ibin]->Branch("w",&wgt,"w/D");
    passTreeEtaPhiv[ibin]->SetDirectory(0);
    sprintf(tname,"failEtaPhi_%i",ibin);
    failTreeEtaPhiv.push_back(new TTree(tname,""));
    failTreeEtaPhiv[ibin]->Branch("m",&mass,"m/F");
    failTreeEtaPhiv[ibin]->Branch("w",&wgt,"w/D");
    failTreeEtaPhiv[ibin]->SetDirectory(0);
  }

  vector<TTree*> passTreeNPVv;  
  vector<TTree*> failTreeNPVv;
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    sprintf(tname,"passNPV_%i",ibin);
    passTreeNPVv.push_back(new TTree(tname,""));
    passTreeNPVv[ibin]->Branch("m",&mass,"m/F");
    passTreeNPVv[ibin]->Branch("w",&wgt,"w/D");
    passTreeNPVv[ibin]->SetDirectory(0);
    sprintf(tname,"failNPV_%i",ibin);
    failTreeNPVv.push_back(new TTree(tname,""));
    failTreeNPVv[ibin]->Branch("m",&mass,"m/F");
    failTreeNPVv[ibin]->Branch("w",&wgt,"w/D");
    failTreeNPVv[ibin]->SetDirectory(0);
  }  
  
  //
  // Pile-up reweighting functions 
  //
  //TFile *pufile    = 0;
  TH1D  *puWeights = 0;
  //if(abs(doPU)==1) pufile = new TFile(pufnameA);
  //if(abs(doPU)==2) pufile = new TFile(pufnameB);
  //if(abs(doPU)==3) pufile = new TFile(pufnameFull);
  //if(doPU!=0) {
  //  assert(pufile);
  //  puWeights = (TH1D*)pufile->Get("puWeights");
  //}
  
  //
  // Generate histogram templates from MC if necessary
  //
  int doType = 2; // set 5 to be the powheg+pythia and 6 to be powheg+photos
  if(sigModPass==2 || sigModFail==2 || sigModPass==5 || sigModFail==5 || sigModPass==6 || sigModFail==6) {
    doType = sigModPass;
    generateHistTemplates(mcfilename,ptBinEdgesv,etaBinEdgesv,phiBinEdgesv,npvBinEdgesv,fitMassLo,fitMassHi,doAbsEta,charge,puWeights,doType);
  }
  if(sigModPass==4 || sigModFail==4) {
    generateDataTemplates(mcfilename,ptBinEdgesv,etaBinEdgesv,phiBinEdgesv,npvBinEdgesv,fitMassLo,fitMassHi,doAbsEta,charge);
  }
    
  //
  // Read in probes data
  //
  TFile *infile    = new TFile(infilename);
  TTree *eventTree = (TTree*)infile->Get("Events");

  Float_t pt, eta, phi;
  Double_t weight;
  Int_t q;
  UInt_t npv, npu, pass, runNum, lumiSec, evtNum;
  Double_t weightPowPyth, weightPowPhot;
  eventTree->SetBranchAddress("mass",   &mass);
  eventTree->SetBranchAddress("pt",     &pt);
  eventTree->SetBranchAddress("eta",    &eta);
  eventTree->SetBranchAddress("phi",    &phi);
  eventTree->SetBranchAddress("weight", &weight);
  eventTree->SetBranchAddress("q",      &q);
  eventTree->SetBranchAddress("npv",    &npv);
  eventTree->SetBranchAddress("npu",    &npu);
  eventTree->SetBranchAddress("pass",   &pass);
  eventTree->SetBranchAddress("runNum", &runNum);
  eventTree->SetBranchAddress("lumiSec",&lumiSec);
  eventTree->SetBranchAddress("evtNum", &evtNum);
  eventTree->SetBranchAddress("weightPowPyth", &weightPowPyth);
  eventTree->SetBranchAddress("weightPowPhot", &weightPowPhot);
  std::cout << "data" << std::endl;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    eventTree->GetEntry(ientry);
    if((q)*charge < 0)    continue;
    if(mass < fitMassLo)  continue;
    if(mass > fitMassHi)  continue;
    if(runNum < runNumLo) continue;
    if(runNum > runNumHi) continue;
    wgt  = weight;
    if(doPU>0) wgt *= puWeights->GetBinContent(npu+1);
    // wgt*=weightPowPhot;
    
    Int_t ipt=-1;
    for(UInt_t ibin=0; ibin<ptNbins; ibin++)
      if((pt >= ptBinEdgesv[ibin]) && (pt < ptBinEdgesv[ibin+1]))
        ipt = ibin;
    if(ipt<0) continue;
    
    Int_t ieta=-1;
    for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
      if(doAbsEta) {
        assert(etaBinEdgesv[ibin]>=0);
        if((fabs(eta) >= etaBinEdgesv[ibin]) && (fabs(eta) < etaBinEdgesv[ibin+1]))
          ieta = ibin;
      } else {
        if((eta >= etaBinEdgesv[ibin]) && (eta < etaBinEdgesv[ibin+1]))
          ieta = ibin;
      }
    }
    if(ieta<0) continue;
	
    Int_t iphi=-1;
    for(UInt_t ibin=0; ibin<phiNbins; ibin++)
      if((phi >= phiBinEdgesv[ibin]) && (phi < phiBinEdgesv[ibin+1]))
        iphi = ibin;
    if(iphi<0) continue;

    Int_t inpv=-1;
    for(UInt_t ibin=0; ibin<npvNbins; ibin++)
      if((npv >= npvBinEdgesv[ibin]) && (npv < npvBinEdgesv[ibin+1]))
        inpv = ibin;
    if(inpv<0) continue;
        
       // fit only the background, skip signal
    
    
    if(pass) {
      passTreePtv[ipt]->Fill();
      passTreeEtav[ieta]->Fill();
      passTreePhiv[iphi]->Fill();
      passTreeEtaPtv[ipt*etaNbins + ieta]->Fill();
      passTreeEtaPhiv[iphi*etaNbins + ieta]->Fill();
      if(inpv>=0)      passTreeNPVv[inpv]->Fill();
    } else {
        // if (mass > 80 && mass < 100 ) continue;       
      failTreePtv[ipt]->Fill();
      failTreeEtav[ieta]->Fill();
      failTreePhiv[iphi]->Fill();
      failTreeEtaPtv[ipt*etaNbins + ieta]->Fill();
      failTreeEtaPhiv[iphi*etaNbins + ieta]->Fill();
      if(inpv>=0)      failTreeNPVv[inpv]->Fill();
    }    
  }  
  delete infile;
  infile=0, eventTree=0;
  std::cout << "done data" << std::endl;
  //
  // Compute efficiencies and make plots 
  // 
  TCanvas *c = MakeCanvas("c","c",800,600);
   
  TGraphAsymmErrors *grEffPt=0;
  TGraphAsymmErrors *grEffEta=0;
  TGraphAsymmErrors *grEffPhi=0;
  TGraphAsymmErrors *grEffNPV=0;
  TH2D *hEffEtaPt   = new TH2D("hEffEtaPt","",etaNbins,etaEdges,ptNbins,ptEdges);
  TH2D *hErrlEtaPt  = (TH2D*)hEffEtaPt->Clone("hErrlEtaPt");
  TH2D *hErrhEtaPt  = (TH2D*)hEffEtaPt->Clone("hErrhEtaPt");
  TH2D *hEffEtaPhi  = new TH2D("hEffEtaPhi","",etaNbins,etaEdges,phiNbins,phiEdges);
  TH2D *hErrlEtaPhi = (TH2D*)hEffEtaPhi->Clone("hErrlEtaPhi");
  TH2D *hErrhEtaPhi = (TH2D*)hEffEtaPhi->Clone("hErrhEtaPhi");
    
  char lumitext[100]; // lumi label
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);    
  std::cout << "do some fits" << std::endl;
  if(sigModPass==0 && sigModFail==0) {  // probe counting
    
    // efficiency in pT
    if(opts[0]) {
      grEffPt = makeEffGraph(ptBinEdgesv, passTreePtv, failTreePtv, method, "pt", massLo, massHi, format, doAbsEta,lumi);
      grEffPt->SetName("grEffPt");
      CPlot plotEffPt("effpt","",xaxislabel+" p_{T} [GeV]",yaxislabel+" efficiency");
      plotEffPt.AddGraph(grEffPt,"",kBlack);
      plotEffPt.SetYRange(ylow,yhigh);
      plotEffPt.SetXRange(0.9*(ptBinEdgesv[0]),1.1*(ptBinEdgesv[ptNbins-1]));
      //plotEffPt.AddTextBox(lumitext,0.55,0.78,0.90,0.84,0);
      plotEffPt.AddTextBox(lumitext,0.55,0.84,0.90,0.90,0);
      plotEffPt.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
      plotEffPt.Draw(c,kTRUE,format);    
    }

    // efficiency in eta
    if(opts[1]) {
      grEffEta = makeEffGraph(etaBinEdgesv, passTreeEtav, failTreeEtav, method, "eta", massLo, massHi, format, doAbsEta,lumi);
      grEffEta->SetName("grEffEta");
      CPlot plotEffEta("effeta","",xaxislabel+" #eta",yaxislabel+" efficiency");
      if(doAbsEta) plotEffEta.SetXTitle("probe |#eta|");
      plotEffEta.AddGraph(grEffEta,"",kBlack);
      plotEffEta.SetYRange(ylow,yhigh);
      //plotEffEta.AddTextBox(lumitext,0.55,0.70,0.90,0.76,0);
      plotEffEta.AddTextBox(lumitext,0.55,0.84,0.90,0.90,0);
      plotEffEta.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
      plotEffEta.Draw(c,kTRUE,format);
    
      CPlot plotEffEta2("effeta2","","probe #eta","#varepsilon");
      if(doAbsEta) plotEffEta2.SetXTitle("probe |#eta|");
      plotEffEta2.AddGraph(grEffEta,"",kBlack);
      plotEffEta2.SetYRange(ylow,yhigh);
      //plotEffEta2.AddTextBox(lumitext,0.55,0.70,0.90,0.76,0);
      plotEffEta2.AddTextBox(lumitext,0.55,0.84,0.90,0.90,0);
      plotEffEta2.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
      plotEffEta2.Draw(c,kTRUE,format);    
    }

    // efficiency in phi
    if(opts[2]) {
      grEffPhi = makeEffGraph(phiBinEdgesv, passTreePhiv, failTreePhiv, method, "phi", massLo, massHi, format, doAbsEta,lumi);
      grEffPhi->SetName("grEffPhi");
      CPlot plotEffPhi("effphi","","probe #phi","#varepsilon");
      plotEffPhi.AddGraph(grEffPhi,"",kBlack);
      plotEffPhi.SetYRange(ylow,yhigh);
      plotEffPhi.Draw(c,kTRUE,format);   
    }

    // efficiency in N_PV
    if(opts[3]) {
      grEffNPV = makeEffGraph(npvBinEdgesv, passTreeNPVv, failTreeNPVv, method, "npv", massLo, massHi, format, doAbsEta,lumi);
      grEffNPV->SetName("grEffNPV");
      CPlot plotEffNPV("effnpv","","N_{PV}","#varepsilon");
      plotEffNPV.AddGraph(grEffNPV,"",kBlack);
      plotEffNPV.SetYRange(ylow,yhigh);
      plotEffNPV.Draw(c,kTRUE,format);   
    }
        
    gStyle->SetPalette(1);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    
    //
    // eta-pT efficiency maps
    //
    if(opts[4]) {
      makeEffHist2D(hEffEtaPt, hErrlEtaPt, hErrhEtaPt, passTreeEtaPtv, failTreeEtaPtv, method, "etapt", massLo, massHi, format, doAbsEta,lumi);
      hEffEtaPt->SetTitleOffset(1.2,"Y");
      if(ptNbins>2)
        hEffEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-2]);
      CPlot plotEffEtaPt("effetapt","","probe #eta","probe p_{T} [GeV/c]");
      if(doAbsEta) plotEffEtaPt.SetXTitle("probe |#eta|");
      plotEffEtaPt.AddHist2D(hEffEtaPt,"COLZ");
      plotEffEtaPt.Draw(c,kTRUE,format);    

      hErrlEtaPt->SetTitleOffset(1.2,"Y");
      if(ptNbins>2)
        hErrlEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-2]);
      CPlot plotErrlEtaPt("errletapt","","probe #eta","probe p_{T} [GeV/c]");
      if(doAbsEta) plotErrlEtaPt.SetXTitle("probe |#eta|");
      plotErrlEtaPt.AddHist2D(hErrlEtaPt,"COLZ");
      plotErrlEtaPt.Draw(c,kTRUE,format);
  
      hErrhEtaPt->SetTitleOffset(1.2,"Y");
      if(ptNbins>2)
        hErrhEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-2]);
      CPlot plotErrhEtaPt("errhetapt","","probe #eta","probe p_{T} [GeV/c]");
      if(doAbsEta) plotErrhEtaPt.SetXTitle("probe |#eta|");
      plotErrhEtaPt.AddHist2D(hErrhEtaPt,"COLZ");
      plotErrhEtaPt.Draw(c,kTRUE,format);
    }
    
    //
    // eta-phi efficiency maps
    //
    if(opts[5]) {
      makeEffHist2D(hEffEtaPhi, hErrlEtaPhi, hErrhEtaPhi, passTreeEtaPhiv, failTreeEtaPhiv, method, "etaphi", massLo, massHi, format, doAbsEta,lumi);
      hEffEtaPhi->SetTitleOffset(1.2,"Y");
      CPlot plotEffEtaPhi("effetaphi","","probe #eta","probe #phi");
      if(doAbsEta) plotEffEtaPhi.SetXTitle("probe |#eta|");
      plotEffEtaPhi.AddHist2D(hEffEtaPhi,"COLZ");
      plotEffEtaPhi.Draw(c,kTRUE,format);   

      hErrlEtaPhi->SetTitleOffset(1.2,"Y");
      CPlot plotErrlEtaPhi("errletaphi","","probe #eta","probe #phi");
      plotErrlEtaPhi.AddHist2D(hErrlEtaPhi,"COLZ");
      if(doAbsEta) plotErrlEtaPhi.SetXTitle("probe |#eta|");
      plotErrlEtaPhi.Draw(c,kTRUE,format);

      hErrhEtaPhi->SetTitleOffset(1.2,"Y");
      CPlot plotErrhEtaPhi("errhetaphi","","probe #eta","probe #phi");
      if(doAbsEta) plotErrhEtaPhi.SetXTitle("probe |#eta|");
      plotErrhEtaPhi.AddHist2D(hErrhEtaPhi,"COLZ");
      plotErrhEtaPhi.Draw(c,kTRUE,format);
    }
       
  } else {
    // efficiency in pT
    if(opts[0]) {
      grEffPt = makeEffGraph(ptBinEdgesv, passTreePtv, failTreePtv, sigModPass, bkgModPass, sigModFail, bkgModFail, "pt", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);
      grEffPt->SetName("grEffPt");
      CPlot plotEffPt("effpt","","probe p_{T} [GeV/c]","#varepsilon");
      plotEffPt.AddGraph(grEffPt,"",kBlack);
      plotEffPt.SetYRange(ylow,yhigh);
      plotEffPt.SetXRange(0.9*(ptBinEdgesv[0]),1.1*(ptBinEdgesv[ptNbins-1]));
      plotEffPt.Draw(c,kTRUE,format);
    }
        
    // efficiency in eta
    if(opts[1]) {
      grEffEta = makeEffGraph(etaBinEdgesv, passTreeEtav, failTreeEtav, sigModPass, bkgModPass, sigModFail, bkgModFail, "eta", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);
      grEffEta->SetName("grEffEta");
      CPlot plotEffEta("effeta","","probe #eta","#varepsilon");
      if(doAbsEta) plotEffEta.SetXTitle("probe |#eta|");
      plotEffEta.AddGraph(grEffEta,"",kBlack);
      plotEffEta.SetYRange(0.2,1.04);
      plotEffEta.Draw(c,kTRUE,format);
    
      CPlot plotEffEta2("effeta2","","probe #eta","#varepsilon");
      if(doAbsEta) plotEffEta2.SetXTitle("probe |#eta|");
      plotEffEta2.AddGraph(grEffEta,"",kBlack);
      plotEffEta2.SetYRange(ylow,yhigh);
      plotEffEta2.Draw(c,kTRUE,format);
    }    
    
    // efficiency in phi
    if(opts[2]) {
      grEffPhi = makeEffGraph(phiBinEdgesv, passTreePhiv, failTreePhiv, sigModPass, bkgModPass, sigModFail, bkgModFail, "phi", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);
      grEffPhi->SetName("grEffPhi");
      CPlot plotEffPhi("effphi","","probe #phi","#varepsilon");
      plotEffPhi.AddGraph(grEffPhi,"",kBlack);
      plotEffPhi.SetYRange(ylow,yhigh);
      plotEffPhi.Draw(c,kTRUE,format);
    }
    
    // efficiency in N_PV
    if(opts[3]) {
      grEffNPV = makeEffGraph(npvBinEdgesv, passTreeNPVv, failTreeNPVv, sigModPass, bkgModPass, sigModFail, bkgModFail, "npv", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);
      grEffNPV->SetName("grEffNPV");
      CPlot plotEffNPV("effnpv","","N_{PV}","#varepsilon");
      plotEffNPV.AddGraph(grEffNPV,"",kBlack);
      plotEffNPV.SetYRange(ylow,yhigh);
      plotEffNPV.Draw(c,kTRUE,format);
    }
    
            
    gStyle->SetPalette(1);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    
    //
    // eta-pT efficiency maps
    //
    if(opts[4]) {
        std::cout << "making efficiency maps" << std::endl;
      makeEffHist2D(hEffEtaPt, hErrlEtaPt, hErrhEtaPt, passTreeEtaPtv, failTreeEtaPtv, sigModPass, bkgModPass, sigModFail, bkgModFail, "etapt", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);
      std::cout << "asf" << std::endl;
      hEffEtaPt->SetTitleOffset(1.2,"Y");
      if(ptBinEdgesv.size()<3) hEffEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-1]);
      else                     hEffEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-2]);
      CPlot plotEffEtaPt("effetapt","","probe #eta","probe p_{T} [GeV/c]");
      if(doAbsEta) plotEffEtaPt.SetXTitle("probe |#eta|");
      plotEffEtaPt.AddHist2D(hEffEtaPt,"COLZ");
      plotEffEtaPt.Draw(c,kTRUE,format);    

      hErrlEtaPt->SetTitleOffset(1.2,"Y");
      if(ptBinEdgesv.size()<3) hErrlEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-1]);
      else                     hErrlEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-2]);
      CPlot plotErrlEtaPt("errletapt","","probe #eta","probe p_{T} [GeV/c]");
      if(doAbsEta) plotErrlEtaPt.SetXTitle("probe |#eta|");
      plotErrlEtaPt.AddHist2D(hErrlEtaPt,"COLZ");
      plotErrlEtaPt.Draw(c,kTRUE,format);

      hErrhEtaPt->SetTitleOffset(1.2,"Y");
      if(ptBinEdgesv.size()<3) hErrhEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-1]);
      else                     hErrhEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins-2]);
      CPlot plotErrhEtaPt("errhetapt","","probe #eta","probe p_{T} [GeV/c]");
      if(doAbsEta) plotErrhEtaPt.SetXTitle("probe |#eta|");
      plotErrhEtaPt.AddHist2D(hErrhEtaPt,"COLZ");
      plotErrhEtaPt.Draw(c,kTRUE,format);
    }
    
    //
    // eta-phi efficiency maps
    //
    if(opts[5]) {
      makeEffHist2D(hEffEtaPhi, hErrlEtaPhi, hErrhEtaPhi, passTreeEtaPhiv, failTreeEtaPhiv, sigModPass, bkgModPass, sigModFail, bkgModFail, "etaphi", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);
      hEffEtaPhi->SetTitleOffset(1.2,"Y");
      CPlot plotEffEtaPhi("effetaphi","","probe #eta","probe #phi");
      if(doAbsEta) plotEffEtaPhi.SetXTitle("probe |#eta|");
      plotEffEtaPhi.AddHist2D(hEffEtaPhi,"COLZ");
      plotEffEtaPhi.Draw(c,kTRUE,format);    

      hErrlEtaPhi->SetTitleOffset(1.2,"Y");
      CPlot plotErrlEtaPhi("errletaphi","","probe #eta","probe #phi");
      plotErrlEtaPhi.AddHist2D(hErrlEtaPhi,"COLZ");
      if(doAbsEta) plotErrlEtaPhi.SetXTitle("probe |#eta|");
      plotErrlEtaPhi.Draw(c,kTRUE,format);
  
      hErrhEtaPhi->SetTitleOffset(1.2,"Y");
      CPlot plotErrhEtaPhi("errhetaphi","","probe #eta","probe #phi");
      if(doAbsEta) plotErrhEtaPhi.SetXTitle("probe |#eta|");
      plotErrhEtaPhi.AddHist2D(hErrhEtaPhi,"COLZ");
      plotErrhEtaPhi.Draw(c,kTRUE,format);
    }
  }

  // Undo scaling of axes before saving to file
  hEffEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins]);
  hErrlEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins]);
  hErrhEtaPt->GetYaxis()->SetRangeUser(ptBinEdgesv[0],ptBinEdgesv[ptNbins]);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  TFile *outfile = new TFile(outputDir + TString("/eff.root"), "RECREATE");
  if(grEffPt)  grEffPt->Write();
  if(grEffEta) grEffEta->Write();
  if(grEffPhi) grEffPhi->Write();
  if(grEffNPV) grEffNPV->Write();
  hEffEtaPt->Write();
  hErrlEtaPt->Write();
  hErrhEtaPt->Write();
  hEffEtaPhi->Write();
  hErrlEtaPhi->Write();
  hErrhEtaPhi->Write();
  outfile->Close();
  delete outfile;
  
  makeHTML(outputDir);
  makeHTML(outputDir, "pt", ptNbins);
  makeHTML(outputDir, "eta", etaNbins);
  makeHTML(outputDir, "phi", phiNbins);
  makeHTML(outputDir, "npv", npvNbins);  
  makeHTML(outputDir, "etapt", etaNbins*ptNbins);
  makeHTML(outputDir, "etaphi", etaNbins*phiNbins);

  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());

  // To make an eta-phi efficiency table using LaTex
  ofstream latexfile;
  char latexfname[100];    
  sprintf(latexfname,"%s/latex.txt",outputDir.Data());
  latexfile.open(latexfname);
  assert(latexfile.is_open());

  CEffUser2D effetapt;
  CEffUser2D effetaphi;
 
  CEffUser1D effpt;
  CEffUser1D effeta;
  CEffUser1D effphi;
  CEffUser1D effnpv;
  
  if(hEffEtaPt->GetEntries()>0) {
    effetapt.loadEff(hEffEtaPt,hErrlEtaPt,hErrhEtaPt);
    effetapt.printEff(txtfile);     txtfile << endl;
    effetapt.printErrLow(txtfile);  txtfile << endl;
    effetapt.printErrHigh(txtfile); txtfile << endl;
    txtfile << endl;

    effetapt.printEffLatex(latexfile); latexfile << endl;
  }
  
  if(hEffEtaPhi->GetEntries()>0) {
    effetaphi.loadEff(hEffEtaPhi,hErrlEtaPhi,hErrhEtaPhi);
    effetaphi.printEff(txtfile);     txtfile << endl;
    effetaphi.printErrLow(txtfile);  txtfile << endl;
    effetaphi.printErrHigh(txtfile); txtfile << endl;
    txtfile << endl;
  }
  
  if(grEffPt) {
    effpt.loadEff(grEffPt);
    effpt.printEff(txtfile);
    txtfile << endl;
    effpt.printErrLow(txtfile);
    txtfile << endl;
    effpt.printErrHigh(txtfile);
    txtfile << endl;
    txtfile << endl;
  }
  
  if(grEffEta) {
    effeta.loadEff(grEffEta);
    effeta.printEff(txtfile);
    txtfile << endl;
    effeta.printErrLow(txtfile);
    txtfile << endl;
    effeta.printErrHigh(txtfile);
    txtfile << endl;
    txtfile << endl;
  }
  
  if(grEffPhi) {
    effphi.loadEff(grEffPhi);
    effphi.printEff(txtfile);
    txtfile << endl;
    effphi.printErrLow(txtfile);
    txtfile << endl;
    effphi.printErrHigh(txtfile);
    txtfile << endl;
    txtfile << endl;
  }
  
  if(grEffNPV) {
    effnpv.loadEff(grEffNPV);
    effnpv.printEff(txtfile);
    txtfile << endl;
    effnpv.printErrLow(txtfile);
    txtfile << endl;
    effnpv.printErrHigh(txtfile);
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("plotEff"); 
  
}  


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;

  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effpt.png\"><img src=\"plots/effpt.png\" alt=\"plots/effpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effeta.png\"><img src=\"plots/effeta.png\" alt=\"plots/effeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effeta2.png\"><img src=\"plots/effeta2.png\" alt=\"plots/effeta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effphi.png\"><img src=\"plots/effphi.png\" alt=\"plots/effphi.png\" width=\"100%\"></a></td>" << endl;  
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"pt.html\">pT bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"eta.html\">&eta; bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"phi.html\">&phi; bins</a></td>" << endl;  
  htmlfile << "</tr>" << endl;    
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effetapt.png\"><img src=\"plots/effetapt.png\" alt=\"plots/effetapt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errletapt.png\"><img src=\"plots/errletapt.png\" alt=\"plots/errletapt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errhetapt.png\"><img src=\"plots/errhetapt.png\" alt=\"plots/errhetapt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;    
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"etapt.html\">&eta;-pT bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effetaphi.png\"><img src=\"plots/effetaphi.png\" alt=\"plots/effetaphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errletaphi.png\"><img src=\"plots/errletaphi.png\" alt=\"plots/errletaphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errhetaphi.png\"><img src=\"plots/errhetaphi.png\" alt=\"plots/errhetaphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;    
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"etaphi.html\">&eta;-&phi; bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effnpv.png\"><img src=\"plots/effnpv.png\" alt=\"plots/effnpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;  
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"npv.html\">N_PV bins</a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;  
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
    
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
}

void makeHTML(const TString outDir, const TString name, const Int_t n)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/%s.html",outDir.Data(),name.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;

  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;    
  Int_t i;
  for(i=0; i<n; i++) {
    if(i%2==0) htmlfile << "<tr>" << endl;    
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pass" << name << "_" << i << ".png\"><img src=\"plots/pass" << name << "_" << i << ".png\"alt=\"plots/pass" << name << "_" << i << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fail" << name << "_" << i << ".png\"><img src=\"plots/fail" << name << "_" << i << ".png\"alt=\"plots/fail" << name << "_" << i << ".png\" width=\"100%\"></a></td>" << endl;
    if(i%2) htmlfile << "</tr>" << endl;
  }
  if(i%2) {
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
    
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
}

//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* makeEffGraph(const vector<Double_t> &edgesv, const vector<TTree*> &passv, const vector<TTree*> &failv, const Int_t method, 
				const TString name, const Double_t massLo, const Double_t massHi, const TString format, const Bool_t doAbsEta,const double lumi)
{
  const UInt_t n = edgesv.size()-1;
  Double_t xval[n], xerr[n];
  Double_t yval[n], yerrl[n], yerrh[n];
  
  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);
    
  for(UInt_t ibin=0; ibin<n; ibin++) {
    xval[ibin] = 0.5*(edgesv[ibin+1] + edgesv[ibin]);
    xerr[ibin] = 0.5*(edgesv[ibin+1] - edgesv[ibin]);
    
    Double_t eff, errl, errh;
    performCount(eff, errl, errh, ibin, edgesv[ibin], edgesv[ibin+1], 0, 0,
	         passv[ibin], failv[ibin], method, 
	         name, massLo, massHi, format, doAbsEta,
	         cpass, cfail,lumi);
    
    yval[ibin]  = eff;
    yerrl[ibin] = errl;
    yerrh[ibin] = errh;
  }
  delete cpass;
  delete cfail;
  
  return new TGraphAsymmErrors(n,xval,yval,xerr,xerr,yerrl,yerrh);
}

TGraphAsymmErrors* makeEffGraph(const vector<Double_t> &edgesv, const vector<TTree*> &passv, const vector<TTree*> &failv,
                                const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		                const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
				const TString format, const Bool_t doAbsEta,const double lumi, const TString yaxislabel, const int charge)
{
  const UInt_t n = edgesv.size()-1;
  Double_t xval[n], xerr[n];
  Double_t yval[n], yerrl[n], yerrh[n];
  
  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);

  for(UInt_t ibin=0; ibin<n; ibin++) {
    xval[ibin] = 0.5*(edgesv[ibin+1] + edgesv[ibin]);
    xerr[ibin] = 0.5*(edgesv[ibin+1] - edgesv[ibin]);

    ifstream rfile;
    char rname[100];
    sprintf(rname,"%s/fitres%s_%i.txt",CPlot::sOutDir.Data(),name.Data(),ibin);
    rfile.open(rname);

    Double_t eff, errl, errh;
    if(rfile.is_open()) {	
      parseFitResults(rfile,eff,errl,errh);
      rfile.close();
	
    } else {
      performFit(eff, errl, errh, ibin, edgesv[ibin], edgesv[ibin+1], 0, 0,
	         passv[ibin], failv[ibin],
	         sigpass, bkgpass, sigfail, bkgfail, 
	         name, massLo, massHi, fitMassLo, fitMassHi, 
		 format, doAbsEta, cpass, cfail, lumi, yaxislabel, charge);
         
       // performFitBkgOnly(eff, errl, errh, ibin, edgesv[ibin], edgesv[ibin+1], 0, 0,
	         // passv[ibin], failv[ibin],
	         // sigpass, bkgpass, sigfail, bkgfail, 
	         // name, massLo, massHi, fitMassLo, fitMassHi, 
		 // format, doAbsEta, cpass, cfail, lumi, yaxislabel, charge);
    }
    
    yval[ibin]  = eff;
    yerrl[ibin] = errl;
    yerrh[ibin] = errh; 
  }
  delete cpass;
  delete cfail;    
  
  return new TGraphAsymmErrors(n,xval,yval,xerr,xerr,yerrl,yerrh);
}

//--------------------------------------------------------------------------------------------------
void makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, const vector<TTree*> &passv, const vector<TTree*> &failv, const Int_t method,
                   const TString name, const Double_t massLo, const Double_t massHi, const TString format, const Bool_t doAbsEta,const double lumi)
{
  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);
    
  for(Int_t iy=0; iy<hEff->GetNbinsY(); iy++) {
    for(Int_t ix=0; ix<hEff->GetNbinsX(); ix++) {
      Int_t ibin = iy*(hEff->GetNbinsX()) + ix;
      
      Double_t eff, errl, errh;
      performCount(eff, errl, errh, ibin, 
                   hEff->GetXaxis()->GetBinLowEdge(ix+1), hEff->GetXaxis()->GetBinLowEdge(ix+2),
		   hEff->GetYaxis()->GetBinLowEdge(iy+1), hEff->GetYaxis()->GetBinLowEdge(iy+2),
		   passv[ibin], failv[ibin], method, 
		   name, massLo, massHi, format, doAbsEta,
		   cpass, cfail,lumi);
      
      hEff ->SetBinContent(hEff ->GetBin(ix+1, iy+1), eff);
      hErrl->SetBinContent(hErrl->GetBin(ix+1, iy+1), errl);
      hErrh->SetBinContent(hErrh->GetBin(ix+1, iy+1), errh);
    }    
  }  
  delete cpass;
  delete cfail;  
}

void makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, const vector<TTree*> &passv, const vector<TTree*> &failv, 
                   const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		   const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
		   const TString format, const Bool_t doAbsEta, const double lumi, const TString yaxislabel, const int charge)
{ 

  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);

  for(Int_t iy=0; iy<hEff->GetNbinsY(); iy++) {
    for(Int_t ix=0; ix<hEff->GetNbinsX(); ix++) {
      Int_t ibin = iy*(hEff->GetNbinsX()) + ix;
      ifstream rfile;
      char rname[200];
      sprintf(rname,"%s/fitres%s_%i.txt",CPlot::sOutDir.Data(),name.Data(),ibin);
      rfile.open(rname);

      Double_t eff, errl, errh;
      if(rfile.is_open()) {	
	parseFitResults(rfile,eff,errl,errh);
	rfile.close();
	
      } else {
	performFit(eff, errl, errh, ibin, 
                   hEff->GetXaxis()->GetBinLowEdge(ix+1), hEff->GetXaxis()->GetBinLowEdge(ix+2),
		   hEff->GetYaxis()->GetBinLowEdge(iy+1), hEff->GetYaxis()->GetBinLowEdge(iy+2),
		   passv[ibin], failv[ibin],
		   sigpass, bkgpass, sigfail, bkgfail, 
		   name, massLo, massHi, fitMassLo, fitMassHi,
		   format, doAbsEta, cpass, cfail, lumi, yaxislabel, charge);
           
           	// performFitBkgOnly(eff, errl, errh, ibin, 
                   // hEff->GetXaxis()->GetBinLowEdge(ix+1), hEff->GetXaxis()->GetBinLowEdge(ix+2),
		   // hEff->GetYaxis()->GetBinLowEdge(iy+1), hEff->GetYaxis()->GetBinLowEdge(iy+2),
		   // passv[ibin], failv[ibin],
		   // sigpass, bkgpass, sigfail, bkgfail, 
		   // name, massLo, massHi, fitMassLo, fitMassHi,
		   // format, doAbsEta, cpass, cfail, lumi, yaxislabel, charge);
      }

      hEff ->SetBinContent(hEff ->GetBin(ix+1, iy+1), eff);
      hErrl->SetBinContent(hErrl->GetBin(ix+1, iy+1), errl);
      hErrh->SetBinContent(hErrh->GetBin(ix+1, iy+1), errh);
    }    
  }  
  delete cpass;
  delete cfail;
}

//--------------------------------------------------------------------------------------------------
void generateHistTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const TH1D* puWeights, int typeCode)
{
  cout << "Creating histogram templates... " << std::endl; cout.flush();

  char hname[50];
  
  const UInt_t ptNbins  = ptEdgesv.size()-1;
  const UInt_t etaNbins = etaEdgesv.size()-1;
  const UInt_t phiNbins = phiEdgesv.size()-1;
  const UInt_t npvNbins = npvEdgesv.size()-1;
  
  TH1D* passPt[ptNbins];
  TH1D* failPt[ptNbins];
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    sprintf(hname,"passpt_%i",ibin);
    passPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passPt[ibin]->SetDirectory(0);
    sprintf(hname,"failpt_%i",ibin);
    failPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failPt[ibin]->SetDirectory(0);
  }
  
  TH1D* passEta[etaNbins];
  TH1D* failEta[etaNbins];
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
    sprintf(hname,"passeta_%i",ibin);
    passEta[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passEta[ibin]->SetDirectory(0);
    sprintf(hname,"faileta_%i",ibin);
    failEta[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failEta[ibin]->SetDirectory(0);
  }
  
  TH1D* passPhi[phiNbins];  
  TH1D* failPhi[phiNbins];
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    sprintf(hname,"passphi_%i",ibin);
    passPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passPhi[ibin]->SetDirectory(0);
    sprintf(hname,"failphi_%i",ibin);
    failPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failPhi[ibin]->SetDirectory(0);
  }
  
  TH1D* passEtaPt[etaNbins*ptNbins];  
  TH1D* failEtaPt[etaNbins*ptNbins];
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    sprintf(hname,"passetapt_%i",ibin);
    passEtaPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passEtaPt[ibin]->SetDirectory(0);
    sprintf(hname,"failetapt_%i",ibin);
    failEtaPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failEtaPt[ibin]->SetDirectory(0);
  }
  
  TH1D* passEtaPhi[etaNbins*phiNbins];
  TH1D* failEtaPhi[etaNbins*phiNbins];
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    sprintf(hname,"passetaphi_%i",ibin); 
    passEtaPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passEtaPhi[ibin]->SetDirectory(0);
    sprintf(hname,"failetaphi_%i",ibin);
    failEtaPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failEtaPhi[ibin]->SetDirectory(0);
  }

  TH1D* passNPV[npvNbins];  
  TH1D* failNPV[npvNbins];
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    sprintf(hname,"passnpv_%i",ibin);
    passNPV[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passNPV[ibin]->SetDirectory(0);
    sprintf(hname,"failnpv_%i",ibin);
    failNPV[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failNPV[ibin]->SetDirectory(0);
  }
  
  std::cout << "blah" << std::endl;
    
  TFile infile(infilename);
  TTree *intree = (TTree*)infile.Get("Events");

  Float_t mass, pt, eta, phi;
  Double_t weight;
  Int_t q;
  UInt_t npv, npu, pass, runNum, lumiSec, evtNum;
  Double_t weightPowPyth, weightPowPhot;
  intree->SetBranchAddress("mass",   &mass);
  intree->SetBranchAddress("pt",     &pt);
  intree->SetBranchAddress("eta",    &eta);
  intree->SetBranchAddress("phi",    &phi);
  intree->SetBranchAddress("weight", &weight);
  intree->SetBranchAddress("q",      &q);
  intree->SetBranchAddress("npv",    &npv);
  intree->SetBranchAddress("npu",    &npu);
  intree->SetBranchAddress("pass",   &pass);
  intree->SetBranchAddress("runNum", &runNum);
  intree->SetBranchAddress("lumiSec",&lumiSec);
  intree->SetBranchAddress("evtNum", &evtNum);
  intree->SetBranchAddress("weightPowPyth", &weightPowPyth);
  intree->SetBranchAddress("weightPowPhot", &weightPowPhot);
  if(weightPowPyth>10||weightPowPyth<0) weightPowPyth=1;
  if(weightPowPhot>10||weightPowPhot<0) weightPowPhot=1; 
  // if(isnan(nan)) std::cout << "blah nanananana" << std::endl;
  std::cout << "hello" << std::endl;
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    Double_t puWgt=1;
    if(puWeights)
      puWgt = puWeights->GetBinContent(npu+1);
    if(typeCode==5) puWgt*=weightPowPyth;
    // std::cout << "weight " << puWgt << "  pyth weight " << weightPowPyth << std::endl;
    if(typeCode==6) puWgt*=weightPowPhot;
    // std::cout << "weight " << puWgt << "  pyth weight " << weightPowPhot << std::endl;
    
    if((q)*charge < 0) continue;
    
    Int_t ipt=-1;
    for(UInt_t ibin=0; ibin<ptNbins; ibin++)
      if((pt >= ptEdgesv[ibin]) && (pt < ptEdgesv[ibin+1]))
        ipt = ibin;
    if(ipt<0) continue;
    
    Int_t ieta=-1;
    for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
      if(doAbsEta) {
        assert(etaEdgesv[ibin]>=0);
        if((fabs(eta) >= etaEdgesv[ibin]) && (fabs(eta) < etaEdgesv[ibin+1]))
          ieta = ibin;
      } else {
        if((eta >= etaEdgesv[ibin]) && (eta < etaEdgesv[ibin+1]))
          ieta = ibin;
      }
    }
    if(ieta<0) continue;
	
    Int_t iphi=-1;
    for(UInt_t ibin=0; ibin<phiNbins; ibin++)
      if((phi >= phiEdgesv[ibin]) && (phi < phiEdgesv[ibin+1]))
        iphi = ibin;
    if(iphi<0) continue;

    Int_t inpv=-1;
    for(UInt_t ibin=0; ibin<npvNbins; ibin++)
      if((npv >= npvEdgesv[ibin]) && (npv < npvEdgesv[ibin+1]))
        inpv = ibin;
    if(inpv<0) continue;
        
        
    if(pass) {
      passPt[ipt]->Fill(mass,puWgt);
      passEta[ieta]->Fill(mass,puWgt);
      passPhi[iphi]->Fill(mass,puWgt);
      passEtaPt[ipt*etaNbins + ieta]->Fill(mass,puWgt);
      passEtaPhi[iphi*etaNbins + ieta]->Fill(mass,puWgt);
      passNPV[inpv]->Fill(mass,puWgt);
    } else {
      failPt[ipt]->Fill(mass,puWgt);
      failEta[ieta]->Fill(mass,puWgt);
      failPhi[iphi]->Fill(mass,puWgt);
      failEtaPt[ipt*etaNbins + ieta]->Fill(mass,puWgt);
      failEtaPhi[iphi*etaNbins + ieta]->Fill(mass,puWgt);
      failNPV[inpv]->Fill(mass,puWgt);
    }    
  }
  infile.Close();
  std::cout << "blah2" << std::endl;
 
  TFile outfile("histTemplates.root", "RECREATE");
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    passPt[ibin]->Write();
    failPt[ibin]->Write();
    delete passPt[ibin];
    delete failPt[ibin];
  }
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) { 
    passEta[ibin]->Write();
    failEta[ibin]->Write();
    delete passEta[ibin];
    delete failEta[ibin];
  }
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    passPhi[ibin]->Write();
    failPhi[ibin]->Write();
    delete passPhi[ibin];
    delete failPhi[ibin];
  }
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    passEtaPt[ibin]->Write();
    failEtaPt[ibin]->Write();
    delete passEtaPt[ibin];
    delete failEtaPt[ibin];
  }
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    passEtaPhi[ibin]->Write();
    failEtaPhi[ibin]->Write();
    delete passEtaPhi[ibin];
    delete failEtaPhi[ibin];
  }
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    passNPV[ibin]->Write();
    failNPV[ibin]->Write();
    delete passNPV[ibin];
    delete failNPV[ibin];
  }
  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
void generateDataTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge)
{
  cout << "Creating data templates... "; cout.flush();

  char tname[50];
  
  const UInt_t ptNbins  = ptEdgesv.size()-1;
  const UInt_t etaNbins = etaEdgesv.size()-1;
  const UInt_t phiNbins = phiEdgesv.size()-1;
  const UInt_t npvNbins = npvEdgesv.size()-1;
  
  Float_t mass;
  
  TTree* passPt[ptNbins];
  TTree* failPt[ptNbins];
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    sprintf(tname,"passpt_%i",ibin);
    passPt[ibin] = new TTree(tname,"");
    passPt[ibin]->Branch("m",&mass,"m/F");
    passPt[ibin]->SetDirectory(0);    
    sprintf(tname,"failpt_%i",ibin);
    failPt[ibin] = new TTree(tname,"");
    failPt[ibin]->Branch("m",&mass,"m/F");
    failPt[ibin]->SetDirectory(0);
  }
  
  TTree* passEta[etaNbins];
  TTree* failEta[etaNbins];
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
    sprintf(tname,"passeta_%i",ibin);
    passEta[ibin] = new TTree(tname,"");
    passEta[ibin]->Branch("m",&mass,"m/F");
    passEta[ibin]->SetDirectory(0); 
    sprintf(tname,"faileta_%i",ibin);
    failEta[ibin] = new TTree(tname,"");
    failEta[ibin]->Branch("m",&mass,"m/F");
    failEta[ibin]->SetDirectory(0); 
  }
  
  TTree* passPhi[phiNbins];  
  TTree* failPhi[phiNbins];
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    sprintf(tname,"passphi_%i",ibin);
    passPhi[ibin] = new TTree(tname,"");
    passPhi[ibin]->Branch("m",&mass,"m/F");
    passPhi[ibin]->SetDirectory(0); 
    sprintf(tname,"failphi_%i",ibin);
    failPhi[ibin] = new TTree(tname,"");
    failPhi[ibin]->Branch("m",&mass,"m/F");
    failPhi[ibin]->SetDirectory(0);
  }
  
  TTree* passEtaPt[etaNbins*ptNbins];  
  TTree* failEtaPt[etaNbins*ptNbins];
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    sprintf(tname,"passetapt_%i",ibin);
    passEtaPt[ibin] = new TTree(tname,"");
    passEtaPt[ibin]->Branch("m",&mass,"m/F");
    passEtaPt[ibin]->SetDirectory(0); 
    sprintf(tname,"failetapt_%i",ibin);
    failEtaPt[ibin] = new TTree(tname,"");
    failEtaPt[ibin]->Branch("m",&mass,"m/F");
    failEtaPt[ibin]->SetDirectory(0);
  }
  
  TTree* passEtaPhi[etaNbins*phiNbins];
  TTree* failEtaPhi[etaNbins*phiNbins];
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    sprintf(tname,"passetaphi_%i",ibin); 
    passEtaPhi[ibin] = new TTree(tname,"");
    passEtaPhi[ibin]->Branch("m",&mass,"m/F");
    passEtaPhi[ibin]->SetDirectory(0); 
    sprintf(tname,"failetaphi_%i",ibin);
    failEtaPhi[ibin] = new TTree(tname,"");
    failEtaPhi[ibin]->Branch("m",&mass,"m/F");
    failEtaPhi[ibin]->SetDirectory(0);
  }

  TTree* passNPV[npvNbins];  
  TTree* failNPV[npvNbins];
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    sprintf(tname,"passnpv_%i",ibin);
    passNPV[ibin] = new TTree(tname,"");
    passNPV[ibin]->Branch("m",&mass,"m/F");
    passNPV[ibin]->SetDirectory(0); 
    sprintf(tname,"failnpv_%i",ibin);
    failNPV[ibin] = new TTree(tname,"");
    failNPV[ibin]->Branch("m",&mass,"m/F");
    failNPV[ibin]->SetDirectory(0);
  }
  
  TFile infile(infilename);
  TTree *intree = (TTree*)infile.Get("Events");

  Float_t pt, eta, phi, weight;
  Int_t q;
  UInt_t npv, npu, pass, runNum, lumiSec, evtNum;

  intree->SetBranchAddress("mass",   &mass);
  intree->SetBranchAddress("pt",     &pt);
  intree->SetBranchAddress("eta",    &eta);
  intree->SetBranchAddress("phi",    &phi);
  intree->SetBranchAddress("weight", &weight);
  intree->SetBranchAddress("q",      &q);
  intree->SetBranchAddress("npv",    &npv);
  intree->SetBranchAddress("npu",    &npu);
  intree->SetBranchAddress("pass",   &pass);
  intree->SetBranchAddress("runNum", &runNum);
  intree->SetBranchAddress("lumiSec",&lumiSec);
  intree->SetBranchAddress("evtNum", &evtNum);
  
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if((q)*charge < 0) continue;
    if(mass < fitMassLo)  continue;
    if(mass > fitMassHi)  continue;
    
    mass = mass;
    
    Int_t ipt=-1;
    for(UInt_t ibin=0; ibin<ptNbins; ibin++)
      if((pt >= ptEdgesv[ibin]) && (pt < ptEdgesv[ibin+1]))
        ipt = ibin;
    if(ipt<0) continue;
    
    Int_t ieta=-1;
    for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
      if(doAbsEta) {
        assert(etaEdgesv[ibin]>=0);
        if((fabs(eta) >= etaEdgesv[ibin]) && (fabs(eta) < etaEdgesv[ibin+1]))
          ieta = ibin;
      } else {
        if((eta >= etaEdgesv[ibin]) && (eta < etaEdgesv[ibin+1]))
          ieta = ibin;
      }
    }
    if(ieta<0) continue;
	
    Int_t iphi=-1;
    for(UInt_t ibin=0; ibin<phiNbins; ibin++)
      if((phi >= phiEdgesv[ibin]) && (phi < phiEdgesv[ibin+1]))
        iphi = ibin;
    if(iphi<0) continue;
	
    Int_t inpv=-1;
    for(UInt_t ibin=0; ibin<npvNbins; ibin++)
      if((npv >= npvEdgesv[ibin]) && (npv < npvEdgesv[ibin+1]))
        inpv = ibin;
    if(inpv<0) continue;
        
    if(pass) {
      passPt[ipt]->Fill();
      passEta[ieta]->Fill();
      passPhi[iphi]->Fill();
      passEtaPt[ipt*etaNbins + ieta]->Fill();
      passEtaPhi[iphi*etaNbins + ieta]->Fill();
      passNPV[inpv]->Fill();
    } else {
      failPt[ipt]->Fill();
      failEta[ieta]->Fill();
      failPhi[iphi]->Fill();
      failEtaPt[ipt*etaNbins + ieta]->Fill();
      failEtaPhi[iphi*etaNbins + ieta]->Fill();
      failNPV[inpv]->Fill();
    }    
  }
  infile.Close();
 
  TFile outfile("dataTemplates.root", "RECREATE");
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    passPt[ibin]->Write();
    failPt[ibin]->Write();
    delete passPt[ibin];
    delete failPt[ibin];
  }
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) { 
    passEta[ibin]->Write();
    failEta[ibin]->Write();
    delete passEta[ibin];
    delete failEta[ibin];
  }
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    passPhi[ibin]->Write();
    failPhi[ibin]->Write();
    delete passPhi[ibin];
    delete failPhi[ibin];
  }
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    passEtaPt[ibin]->Write();
    failEtaPt[ibin]->Write();
    delete passEtaPt[ibin];
    delete failEtaPt[ibin];
  }
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    passEtaPhi[ibin]->Write();
    failEtaPhi[ibin]->Write();
    delete passEtaPhi[ibin];
    delete failEtaPhi[ibin];
  }
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    passNPV[ibin]->Write();
    failNPV[ibin]->Write();
    delete passNPV[ibin];
    delete failNPV[ibin];
  }
  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
void performCount(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                  const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		  TTree *passTree, TTree *failTree, const Int_t method, 
		  const TString name, const Double_t massLo, const Double_t massHi, const TString format, const Bool_t doAbsEta,
		  TCanvas *cpass, TCanvas *cfail,const double lumi)
{
  // skip ECAL gap region
  // if(xbinLo==1.4442 && xbinHi==1.566) return;
  // if(xbinLo==-1.566 && xbinHi==-1.4442) return;

  Float_t m;
  Double_t w;
  char pname[50];
  char binlabelx[100];
  char binlabely[100];
  char yield[50];
  char ylabel[50];
  char effstr[100];    
  
  Double_t npass=0, ntotal=0;
  passTree->SetBranchAddress("m",&m);
  passTree->SetBranchAddress("w",&w);
  for(UInt_t ientry=0; ientry<passTree->GetEntries(); ientry++) {
    passTree->GetEntry(ientry);
    if(m<massLo || m>massHi) continue;
    npass+=w;
    ntotal+=w;
  }
  failTree->SetBranchAddress("m",&m);
  failTree->SetBranchAddress("w",&w);
  for(UInt_t ientry=0; ientry<failTree->GetEntries(); ientry++) {
    failTree->GetEntry(ientry);
    if(m<massLo || m>massHi) continue;
    ntotal+=w;
  }
  if(ntotal < npass) ntotal = npass;
  resEff  = (ntotal>0) ? npass/ntotal : 0;
  if(method==0) {
    resErrl = resEff - TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kFALSE);
    resErrh = TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kTRUE) - resEff;
  }
   
  if(name.CompareTo("pt")==0) {
    sprintf(binlabelx,"%i GeV/c < p_{T} < %i GeV/c",Int_t(xbinLo),Int_t(xbinHi));
  
  } else if(name.CompareTo("eta")==0) { 
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else	 sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);
  
  } else if(name.CompareTo("phi")==0) { 
    sprintf(binlabelx,"%.1f < #phi < %.1f",xbinLo,xbinHi); 
  
  } else if(name.CompareTo("etapt")==0) {
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);    
    sprintf(binlabely,"%i GeV/c < p_{T} < %i GeV/c",Int_t(ybinLo),Int_t(ybinHi));
  
  } else if(name.CompareTo("etaphi")==0) {
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);					   
    sprintf(binlabely,"%.1f < #phi < %.1f",ybinLo,ybinHi);
  
  } else if(name.CompareTo("npv")==0) { 
    sprintf(binlabelx,"%i #leq N_{PV} < %i",(Int_t)xbinLo,(Int_t)xbinHi); 
  
  }
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);
  
  char lumitext[100]; // lumi label
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);    

  //
  // Plot passing probes
  //
  TH1D *hpass = new TH1D("hpass","",Int_t(massHi-massLo)/BIN_SIZE_PASS,massLo,massHi);
  passTree->SetBranchAddress("m",&m);
  passTree->SetBranchAddress("w",&w);
  hpass->Sumw2();
  for(UInt_t ientry=0; ientry<passTree->GetEntries(); ientry++) {
    passTree->GetEntry(ientry);
    hpass->Fill(m,w);
  }
  sprintf(pname,"pass%s_%i",name.Data(),ibin);
  sprintf(yield,"%i Events",(UInt_t)npass);
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_PASS);
  CPlot plotPass(pname,"Passing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotPass.AddHist1D(hpass,"E");
  plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  if((name.CompareTo("etapt")==0) || (name.CompareTo("etaphi")==0)) {
    plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);        
    plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);    
  } else {
    plotPass.AddTextBox(yield,0.21,0.81,0.51,0.85,0,kBlack,-1);
  }
  plotPass.AddTextBox(effstr,0.70,0.85,0.95,0.90,0,kBlack,-1);
  plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  plotPass.Draw(cpass,kTRUE,format);
  
  //
  // Plot failing probes
  //
  TH1D *hfail = new TH1D("hfail","",Int_t(massHi-massLo)/BIN_SIZE_FAIL,massLo,massHi);
  hfail->Sumw2();
  failTree->SetBranchAddress("m",&m);
  failTree->SetBranchAddress("w",&w);
  for(UInt_t ientry=0; ientry<failTree->GetEntries(); ientry++) {
    failTree->GetEntry(ientry);
    hfail->Fill(m,w);
  }
  sprintf(pname,"fail%s_%i",name.Data(),ibin);
  sprintf(yield,"%i Events",(UInt_t)(ntotal-npass));
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_FAIL);
  CPlot plotFail(pname,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotFail.AddHist1D(hfail,"E");
  plotFail.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  if((name.CompareTo("etapt")==0) || (name.CompareTo("etaphi")==0)) {
    plotFail.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);    
    plotFail.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);    
  } else {
    plotFail.AddTextBox(yield,0.21,0.81,0.51,0.85,0,kBlack,-1);
  }
  plotFail.AddTextBox(effstr,0.70,0.85,0.95,0.90,0,kBlack,-1);
  plotFail.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  plotFail.Draw(cfail,kTRUE,format);
  
  delete hpass;
  delete hfail;
}
//--------------------------------------------------------------------------------------------------
void performFitBkgOnly(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		TTree *passTree, TTree *failTree,
		const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
		const TString format, const Bool_t doAbsEta, TCanvas *cpass, TCanvas *cfail, const double lumi, const TString yaxislabel, const int charge)
{
  // skip ECAL gap region
  // if(xbinLo==1.4442 && xbinHi==1.566) return;
  // if(xbinLo==-1.566 && xbinHi==-1.4442) return;

  RooRealVar m("m","mass",fitMassLo,fitMassHi);
  m.setBins(60);

  char bkgpassname[50];
  char bkgfailname[50];
  char pname[50];
  char binlabelx[100];
  char binlabely[100];
  char yield[50];
  char ylabel[50];
  char effstr[100];
  char nsigstr[100];
  char nbkgstr[100];
  char chi2str[100];
 
  sprintf(bkgpassname, "backgroundPass_%d",ibin);
  sprintf(bkgfailname, "backgroundFail_%d",ibin); 

  Int_t nflpass=0, nflfail=0;
    
  TFile *histfile = 0;
  if(sigpass==2 || sigfail==2 || sigpass==5 || sigfail==5 || sigpass==6 || sigfail==6) {
    histfile = new TFile("histTemplates.root");
    assert(histfile);
  }
  TFile *datfile = 0;
  if(sigpass==4 || sigfail==4) {
    datfile = new TFile("dataTemplates.root");
    assert(datfile);
  }
  
  // Define categories
  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);
  
  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  TH1D histPass("histPass","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS,fitMassLo,fitMassHi); 
  TH1D histFail("histFail","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL,fitMassLo,fitMassHi);
  RooAbsData *dataCombined=0;
  
  const Bool_t doBinned = kTRUE;//(passTree->GetEntries()>1000 && failTree->GetEntries()>1000);
  
  if(doBinned) {
    passTree->Draw("m>>histPass","w");
    failTree->Draw("m>>histFail","w");
    dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),&histPass);
    dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),&histFail);
    //m.setBins(100);  
   
    dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                   RooFit::Index(sample),
				   RooFit::Import("Pass",*((RooDataHist*)dataPass)),
				   RooFit::Import("Fail",*((RooDataHist*)dataFail)));  
  
  } else {
    dataPass = new RooDataSet("dataPass","dataPass",passTree,RooArgSet(m));
    dataFail = new RooDataSet("dataFail","dataFail",failTree,RooArgSet(m));
    
    dataCombined = new RooDataSet("dataCombined","dataCombined",RooArgList(m),
      				  RooFit::Index(sample),
  	  			  RooFit::Import("Pass",*((RooDataSet*)dataPass)),
  				  RooFit::Import("Fail",*((RooDataSet*)dataFail))); 
  }
  
  // Define signal and background models
  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;
  
  

  
  if(sigpass==1) {
    sigPass = new CBreitWignerConvCrystalBall(m,kTRUE);
    nflpass += 4;
  
  } else if(sigpass==2||sigpass==5||sigpass==6) { 
    char hname[50];
    sprintf(hname,"pass%s_%i",name.Data(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    sigPass = new CMCTemplateConvGaussian(m,h,kTRUE,ibin);
    nflpass += 2;
  
  } else if(sigpass==3) {
    sigPass = new CVoigtianCBShape(m,kTRUE);
    nflpass += 4;
  
  } else if(sigpass==4) {
    char tname[50];
    sprintf(tname,"pass%s_%i",name.Data(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigPass = new CMCDatasetConvGaussian(m,t,kTRUE);
    nflpass += 2;
  }

  if(bkgpass==1) { 
    bkgPass = new CExponential(m,kTRUE,ibin);
    nflpass += 1;
  
  } else if(bkgpass==2) {
    bkgPass = new CErfExpo(m,kTRUE,ibin);
    nflpass += 3;
     
  } else if(bkgpass==3) {
    bkgPass = new CDoubleExp(m,kTRUE);
    nflpass += 3;
  
  } else if(bkgpass==4) {
    bkgPass = new CLinearExp(m,kTRUE);
    nflpass += 2;
  
  } else if(bkgpass==5) {
    bkgPass = new CQuadraticExp(m,kTRUE,ibin);
    nflpass += 3;  

  } else if(bkgpass==6) {
    bkgPass = new CQuadratic(m,kTRUE,ibin,0.,0.,0.,0.,0.,0.);
    nflpass += 3;

  } else if(bkgpass==7) {
    bkgPass = new CPower(m,kTRUE,ibin,0.,0.);
    nflfail += 1;
  }


 RooRealVar *sigmaFail = new RooRealVar("sigmaFail","sigmaFail",0.1,0,1.0);
  if(sigfail==1) {
    sigFail = new CBreitWignerConvCrystalBall(m,kFALSE);
    nflfail += 4;
  
  } else if(sigfail==2||sigfail==5||sigfail==6) {
    char hname[50];
    sprintf(hname,"fail%s_%i",name.Data(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    // if(ibin==7||ibin==11){
        // std::cout << "it's ibin " << ibin << std::endl;
        // for(int i=1; i < h->GetNbinsX(); i++){
            // std::cout << "h bin " << i << "  content "<< h->GetBinContent(i) << std::endl;
        // }
    // }
    sigFail = new CMCTemplateConvGaussian(m,h,kFALSE,ibin);//,((CMCTemplateConvGaussian*)sigPass)->sigma);
    nflfail += 2;
  } else if(sigfail==3) {
    sigFail = new CVoigtianCBShape(m,kFALSE);
    nflfail += 4;
  
  } else if(sigfail==4) {
    char tname[50];
    sprintf(tname,"fail%s_%i",name.Data(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigFail = new CMCDatasetConvGaussian(m,t,kFALSE);
    nflfail += 2;
  }

if(yaxislabel.CompareTo("stand-alone")==0       && bkgfail==6 && charge==0 && ibin==0){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-20.0,50.0,16.0,10.0,-0.1,10.0); // bin0 amcnlo 2,6
    // bkgFail = new CQuadratic(m,kFALSE,ibin,-70.0,20.0,13.0,5.0,-0.1,10.0); // bin0 minlo 2,6
    std::cout << "setting bins..." << std::endl;
    nflfail += 3;

  } else if(bkgfail==1) { 
    bkgFail = new CExponential(m,kFALSE,ibin);
    nflfail += 1;
  
  } else if(bkgfail==2) {
    bkgFail = new CErfExpo(m,kFALSE,ibin); 
    nflfail += 3;
  
  } else if(bkgfail==3) {
    bkgFail = new CDoubleExp(m,kFALSE);
    nflfail += 3;
    
  } else if(bkgfail==4) {
    bkgFail = new CLinearExp(m,kFALSE);
    nflfail += 2;
  
  } else if(bkgfail==5) {
    bkgFail = new CQuadraticExp(m,kFALSE,ibin);
    nflfail += 3;  

  } else if(bkgfail==6) {
    bkgFail = new CQuadratic(m,kFALSE,ibin,10.,30.,10.,30.,2.,15.);
    nflfail += 3;

  } else if(bkgpass==7) {
    bkgFail = new CPower(m,kFALSE,ibin,0.,0.);
    nflfail += 1;
  }


  // Define free parameters
  Double_t NsigMax     = doBinned ? histPass.Integral()+histFail.Integral() : passTree->GetEntries()+failTree->GetEntries();
  Double_t NbkgFailMax = doBinned ? histFail.Integral() : failTree->GetEntries();
  Double_t NbkgPassMax = doBinned ? histPass.Integral() : passTree->GetEntries();
  RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0.0,1.5*NsigMax);
  RooRealVar eff("eff","Efficiency",0.95,0.0,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.1*NbkgPassMax,0.01,NbkgPassMax);
  if(bkgpass==0) NbkgPass.setVal(0);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax,0.01,NbkgFailMax);
  if(yaxislabel.CompareTo("stand-alone")==0 && fabs(xbinLo) >= 0.9 ) {
      NbkgFail.setVal(0.95*NbkgFailMax);
      NbkgFail.setRange(0.85*NbkgFailMax,NbkgFailMax);
      eff.setVal(0.95);
  }
   else if(yaxislabel.CompareTo("stand-alone")==0       && bkgfail==6 && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==25 && ybinHi==40  ){
       NbkgFail.setVal(0.5*NbkgFailMax);
   }
  std::cout << NbkgFail.getVal() << std::endl;
  // if(bkgfail==6) NbkgFail.setVal(0.9*NbkgFailMax);
std::cout << NbkgFail.getVal() << std::endl;
std::cout << "-----------------" << std::endl;

  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
  RooAddPdf *modelPass=0, *modelFail=0;
  RooExtendPdf *esignalPass=0, *ebackgroundPass=0, *esignalFail=0, *ebackgroundFail=0;

  if(massLo!=fitMassLo || massHi!=fitMassHi) {
    m.setRange("signalRange",massLo,massHi);
    
    esignalPass     = new RooExtendPdf("esignalPass","esignalPass",*(sigPass->model),NsigPass,"signalRange");
    ebackgroundPass = new RooExtendPdf("ebackgroundPass","ebackgroundPass",(bkgpass>0) ? *(bkgPass->model) : *(sigPass->model),NbkgPass,"signalRange");
    modelPass       = new RooAddPdf("modelPass","Model for PASS sample",(bkgpass>0) ? RooArgList(*esignalPass,*ebackgroundPass) : RooArgList(*esignalPass));    

    esignalFail     = new RooExtendPdf("esignalFail","esignalFail",*(sigFail->model),NsigFail,"signalRange");
    ebackgroundFail = new RooExtendPdf("ebackgroundFail","ebackgroundFail",*(bkgFail->model),NbkgFail,"signalRange");
    // modelFail       = new RooAddPdf("modelFail","Model for FAIL sample", (bkgfail>0) ? RooArgList(*esignalFail,*ebackgroundFail) : RooArgList(*esignalFail));
    modelFail       = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*ebackgroundFail));
  
  } else {

    modelPass = new RooAddPdf("modelPass","Model for PASS sample",
                              (bkgpass>0) ? RooArgList(*(sigPass->model),*(bkgPass->model)) :  RooArgList(*(sigPass->model)),
		              (bkgpass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));

    // modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigFail->model),*(bkgFail->model)),RooArgList(NsigFail,NbkgFail));
    modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(bkgFail->model)),RooArgList(NbkgFail));
  }
  
  // eff.setVal(1.0);
  // eff.setConstant(kTRUE);
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");  
  totalPdf.addPdf(*modelFail,"Fail");

  int strategy = 0;

  // m.setR
  m.setRange("R1",60,75) ;
  m.setRange("R2",105,120) ;
  std::cout << "set ranges and fitting " << std::endl;
  RooFitResult *fitResult=0;
  RooMsgService::instance().setSilentMode(kTRUE);
  fitResult = modelFail->fitTo(*dataFail, 
             RooFit::Range("R1,R2"),
                             RooFit::PrintEvalErrors(-1),
			     RooFit::Extended(),
  			     RooFit::Strategy(strategy), // MINOS STRATEGY
  			     // RooFit::Minos(RooArgSet(eff)),
  			     RooFit::Save());
  std::cout << "fit status " << fitResult->status() << std::endl;
  
  if(fitResult->status()!=0){
      std::cout << "Backupt Fits!" << std::endl;
    // fitResult = modelFail->fitTo(*dataFail, RooFit::Range("R1,R2"), RooFit::PrintEvalErrors(-1), RooFit::Extended(),RooFit::Strategy(1), RooFit::Save());
    fitResult = modelFail->fitTo(*dataFail, RooFit::Range("R1,R2"), RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","scan"),RooFit::Strategy(2), RooFit::Save());
    // // fitResult = modelFail->fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","migrad"),RooFit::Range("R1,R2"), RooFit::Save());
    fitResult = modelFail->fitTo(*dataFail,RooFit::Range("R1,R2"), RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","improve"),RooFit::Strategy(2), RooFit::Save());
    fitResult = modelFail->fitTo(*dataFail,RooFit::Range("R1,R2"), RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","minimize"),RooFit::Strategy(2), RooFit::Save());
    // // // fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),RooFit::Strategy(2),RooFit::Minos(RooArgSet(eff)), RooFit::Save());
  }
  
    // do {
    // fitResult = modelFail->fitTo(*dataFail,
				 // NumCPU(4),
				 // //				 Minimizer("Minuit2","minimize"),
				 // Minimizer("Minuit2","scan"),
                 // RooFit::Range("R1,R2"),
				 // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // RooFit::Minos(),
				 // RooFit::Strategy(2),
				 // RooFit::Save());
                 
        // // nTries++;

    // fitResult = modelFail->fitTo(*dataFail,
				 // NumCPU(4),
				 // //				 Minimizer("Minuit2","minimize"),
				 // Minimizer("Minuit2","migrad"),
                 // RooFit::Range("R1,R2"),
				 // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // RooFit::Hesse(),
				 // RooFit::Strategy(2),
				 // RooFit::Save());
    // fitResult = modelFail->fitTo(*dataFail,
				 // NumCPU(4),
				 // //				 Minimizer("Minuit2","minimize"),
				 // Minimizer("Minuit2","improve"),
                 // RooFit::Range("R1,R2"),
				 // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // RooFit::Minos(),
				 // RooFit::Strategy(2),
				 // RooFit::Save());
				 
	// fitResult = modelFail->fitTo(*dataFail,
			       // NumCPU(4),
			       // Minimizer("Minuit2","minimize"),
                   // RooFit::Range("R1,R2"),
			       // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       // //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       // RooFit::Minos(),
			       // RooFit::Strategy(2),
	               // RooFit::Save());
        // nTries++;
    // }while((fitResult->status()>0)&&nTries < 10);
                 
                 
                 std::cout << "fitted " << std::endl;

  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();
    
  if(name.CompareTo("pt")==0) {
    sprintf(binlabelx,"%i GeV/c < p_{T} < %i GeV/c",Int_t(xbinLo),Int_t(xbinHi));
  
  } else if(name.CompareTo("eta")==0) { 
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else	 sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);
  
  } else if(name.CompareTo("phi")==0) { 
    sprintf(binlabelx,"%.1f < #phi < %.1f",xbinLo,xbinHi); 
  
  } else if(name.CompareTo("etapt")==0) {
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);    
    sprintf(binlabely,"%i GeV/c < p_{T} < %i GeV/c",Int_t(ybinLo),Int_t(ybinHi));
  
  } else if(name.CompareTo("etaphi")==0) {
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);					   
    sprintf(binlabely,"%.1f < #phi < %.1f",ybinLo,ybinHi);
  
  } else if(name.CompareTo("npv")==0) { 
    sprintf(binlabelx,"%i #leq N_{PV} < %i",(Int_t)xbinLo,(Int_t)xbinHi); 
  
  } 
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",eff.getVal(),fabs(eff.getErrorLo()),eff.getErrorHi());

std::cout << "blah1 " << std::endl;

  RooPlot *mframePass = m.frame(Bins(Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS));
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
  if(bkgpass>0)
    modelPass->plotOn(mframePass,Components(bkgpassname/*"backgroundPass"*/),LineStyle(kDashed),LineColor(kRed));
  modelPass->plotOn(mframePass);
std::cout << "blah2 " << std::endl;
  RooPlot *mframeFail = m.frame(Bins(Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL));
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail,Components(bkgfailname/*"backgroundFail"*/),LineStyle(kDashed),LineColor(kRed),Range("FULL"),NormRange("FULL"));
  modelFail->plotOn(mframeFail,Range("R1,R2"),NormRange("R1,R2"));
  std::cout << "blah3 " << std::endl;
  char lumitext[100]; // lumi label
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);    


  //
  // Plot failing probes
  //
  // double f = NsigFail.getVal(), d = NbkgFail.getVal();
  std::cout << "blah3a1" << std::endl;
  sprintf(pname,"fail%s_%i",name.Data(),ibin);
  sprintf(yield,"%u Events",(Int_t)failTree->GetEntries());
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_FAIL);
  // sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigFail.getVal(),NsigFail.getPropagatedError(*fitResult));
  std::cout << "blah3a2" << std::endl;
  // sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgFail.getVal(),NbkgFail.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframePass->chiSquare(nflfail));
  CPlot plotFail(pname,mframeFail,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  std::cout << "blah3a3" << std::endl;
  plotFail.AddTextBox(binlabelx,0.21,0.75,0.51,0.80,0,kBlack,-1);
  std::cout << "blah3b" << std::endl;
  if((name.CompareTo("etapt")==0) || (name.CompareTo("etaphi")==0)) {
    plotFail.AddTextBox(binlabely,0.21,0.70,0.51,0.75,0,kBlack,-1);    
    plotFail.AddTextBox(yield,0.21,0.66,0.51,0.70,0,kBlack,-1);    
  } else {
    plotFail.AddTextBox(yield,0.21,0.81,0.51,0.85,0,kBlack,-1);
  }
  std::cout << "blah3c" << std::endl;
  plotFail.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);  
  // plotFail.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
  plotFail.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,-1);
  plotFail.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  std::cout << "blah3d" << std::endl;
  if(yaxislabel.CompareTo("supercluster")==0 && xbinLo==0.0 && xbinHi==1.4442 && ybinLo==25 && ybinHi==100) plotFail.SetYRange(0,200);
  else if(yaxislabel.CompareTo("supercluster")==0 && xbinLo==1.566 && xbinHi==2.5 && ybinLo==25 && ybinHi==100) plotFail.SetYRange(0,100);
  plotFail.Draw(cfail,kTRUE,format);  
std::cout << "blah4" << std::endl;
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/fitres%s_%i.txt",CPlot::sOutDir.Data(),name.Data(),ibin);
  
std::cout << "blah4a" << std::endl;
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  
std::cout << "blah4b" << std::endl;
  txtfile << endl;
  // printCorrelations(txtfile, fitResult);
// std::cout << "blah4c" << std::endl;
  txtfile.close();

    char outFileName[100];
  sprintf(outFileName, "%s/%s_%i.root",CPlot::sOutDir.Data(),name.Data(),ibin);

  RooWorkspace *w = new RooWorkspace("w", "workspace");
  w->import(*modelFail);
  w->import(*fitResult);

  w->writeToFile(outFileName);

    TFile *fw = new TFile(outFileName, "UPDATE");
    TTree *t = new TTree("Bin", "Bin");

    UInt_t nEvents, nBkgFail, nBkgPass, iBin, absEta;
    Float_t ptLo, ptHi;
    Float_t etaLo, etaHi;
    Float_t phiLo, phiHi;
    Float_t npvLo, npvHi;


std::cout << "blah5" << std::endl;
    nEvents = NsigMax;
    nBkgFail = NbkgFailMax;
    nBkgPass = NbkgPassMax;
    iBin = ibin;
    if (doAbsEta) absEta=1;
    else absEta=0;

    ptLo = 999; ptHi = 999; etaLo = 999; etaHi = 999;
    phiLo = 999; phiHi = 999; npvLo = 999; npvHi = 999;

    if(name.CompareTo("pt")==0) {
      ptLo = xbinLo;
      ptHi = xbinHi;

    } else if(name.CompareTo("eta")==0) {
      etaLo = xbinLo;
      etaHi = xbinHi;

    } else if(name.CompareTo("phi")==0) {
      phiLo = xbinLo;
      phiHi = xbinHi;

    } else if(name.CompareTo("etapt")==0) {
      etaLo = xbinLo;
      etaHi = xbinHi;
      ptLo = ybinLo;
      ptHi = ybinHi;

    } else if(name.CompareTo("etaphi")==0) {
      etaLo = xbinLo;
      etaHi = xbinHi;
      phiLo = ybinLo;
      phiHi = ybinHi;

    } else if(name.CompareTo("npv")==0) {
      npvLo = xbinLo;
      npvHi = xbinHi;

    }

std::cout << "blah6" << std::endl;
    t->Fill();
    fw->Write();
    fw->Close();


  //
  // Clean up
  //
  delete esignalPass;
  delete ebackgroundPass;
  delete esignalFail;
  delete ebackgroundFail;
  delete modelPass;
  delete modelFail;  
  delete dataCombined;
  delete dataPass;
  delete dataFail;
  delete sigPass;
  delete bkgPass;
  delete sigFail;
  delete bkgFail;        
  delete histfile;
  delete datfile;   
  delete sigmaFail;

}

//--------------------------------------------------------------------------------------------------
void performFit(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		TTree *passTree, TTree *failTree,
		const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
		const TString format, const Bool_t doAbsEta, TCanvas *cpass, TCanvas *cfail, const double lumi, const TString yaxislabel, const int charge)
{
  // // skip ECAL gap region
  // if(xbinLo==1.4442 && xbinHi==1.566) return;
  // if(xbinLo==-1.566 && xbinHi==-1.4442) return;

  RooRealVar m("m","mass",fitMassLo,fitMassHi);
  m.setBins(60);

  char bkgpassname[50];
  char bkgfailname[50];
  char pname[50];
  char binlabelx[100];
  char binlabely[100];
  char yield[50];
  char ylabel[50];
  char effstr[100];
  char nsigstr[100];
  char nbkgstr[100];
  char chi2str[100];
 
  sprintf(bkgpassname, "backgroundPass_%d",ibin);
  sprintf(bkgfailname, "backgroundFail_%d",ibin); 

  Int_t nflpass=0, nflfail=0;
    
  TFile *histfile = 0;
  if(sigpass==2 || sigfail==2 || sigpass==5 || sigfail==5 || sigpass==6 || sigfail==6) {
    histfile = new TFile("histTemplates.root");
    assert(histfile);
  }
  TFile *datfile = 0;
  if(sigpass==4 || sigfail==4) {
    datfile = new TFile("dataTemplates.root");
    assert(datfile);
  }
  
  // Define categories
  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);
  
  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  TH1D histPass("histPass","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS,fitMassLo,fitMassHi); 
  TH1D histFail("histFail","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL,fitMassLo,fitMassHi);
  RooAbsData *dataCombined=0;
  
  const Bool_t doBinned = kTRUE;//(passTree->GetEntries()>1000 && failTree->GetEntries()>1000);
  
  if(doBinned) {
    passTree->Draw("m>>histPass","w");
    failTree->Draw("m>>histFail","w");
    dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),&histPass);
    dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),&histFail);
    //m.setBins(100);  
   
    dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                   RooFit::Index(sample),
				   RooFit::Import("Pass",*((RooDataHist*)dataPass)),
				   RooFit::Import("Fail",*((RooDataHist*)dataFail)));  
  
  } else {
    dataPass = new RooDataSet("dataPass","dataPass",passTree,RooArgSet(m));
    dataFail = new RooDataSet("dataFail","dataFail",failTree,RooArgSet(m));
    
    dataCombined = new RooDataSet("dataCombined","dataCombined",RooArgList(m),
      				  RooFit::Index(sample),
  	  			  RooFit::Import("Pass",*((RooDataSet*)dataPass)),
  				  RooFit::Import("Fail",*((RooDataSet*)dataFail))); 
  }
  
  // Define signal and background models
  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;
  
  

  
  if(sigpass==1) {
    sigPass = new CBreitWignerConvCrystalBall(m,kTRUE);
    nflpass += 4;
  
  } else if(sigpass==2||sigpass==5||sigpass==6) { 
    char hname[50];
    sprintf(hname,"pass%s_%i",name.Data(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    sigPass = new CMCTemplateConvGaussian(m,h,kTRUE,ibin);
    nflpass += 2;
  
  } else if(sigpass==3) {
    sigPass = new CVoigtianCBShape(m,kTRUE);
    nflpass += 4;
  
  } else if(sigpass==4) {
    char tname[50];
    sprintf(tname,"pass%s_%i",name.Data(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigPass = new CMCDatasetConvGaussian(m,t,kTRUE);
    nflpass += 2;
  }

  if(bkgpass==1) { 
    bkgPass = new CExponential(m,kTRUE,ibin);
    nflpass += 1;
  
  } else if(bkgpass==2) {
    bkgPass = new CErfExpo(m,kTRUE,ibin);
    nflpass += 3;
     
  } else if(bkgpass==3) {
    bkgPass = new CDoubleExp(m,kTRUE);
    nflpass += 3;
  
  } else if(bkgpass==4) {
    bkgPass = new CLinearExp(m,kTRUE);
    nflpass += 2;
  
  } else if(bkgpass==5) {
    bkgPass = new CQuadraticExp(m,kTRUE,ibin);
    nflpass += 3;  

  } else if(bkgpass==6) {
    bkgPass = new CQuadratic(m,kTRUE,ibin,0.,0.,0.,0.,0.,0.);
    nflpass += 3;

  } else if(bkgpass==7) {
    bkgPass = new CPower(m,kTRUE,ibin,0.,0.);
    nflfail += 1;
  }


 RooRealVar *sigmaFail = new RooRealVar("sigmaFail","sigmaFail",0.1,0,1.0);
  if(sigfail==1) {
    sigFail = new CBreitWignerConvCrystalBall(m,kFALSE);
    nflfail += 4;
  
  } else if(sigfail==2||sigfail==5||sigfail==6) {
    char hname[50];
    sprintf(hname,"fail%s_%i",name.Data(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
    // if(ibin==7||ibin==11){
        // std::cout << "it's ibin " << ibin << std::endl;
        // for(int i=1; i < h->GetNbinsX(); i++){
            // std::cout << "h bin " << i << "  content "<< h->GetBinContent(i) << std::endl;
        // }
    // }
    sigFail = new CMCTemplateConvGaussian(m,h,kFALSE,ibin);//,((CMCTemplateConvGaussian*)sigPass)->sigma);
    nflfail += 2;
  } else if(sigfail==3) {
    sigFail = new CVoigtianCBShape(m,kFALSE);
    nflfail += 4;
  
  } else if(sigfail==4) {
    char tname[50];
    sprintf(tname,"fail%s_%i",name.Data(),ibin);
    TTree *t = (TTree*)datfile->Get(tname);
    assert(t);
    sigFail = new CMCDatasetConvGaussian(m,t,kFALSE);
    nflfail += 2;
  }



// for alternative iso cone, 0.12
// {-2.4,-2.1,-1.2,-0.9,0,0.9,1.2,2.1,2.4}; // muon eta binning kms
/*  }else*/ 

// if(yaxislabel.CompareTo("stand-alone")==0       && bkgfail==6 && charge==0 && fabs(xbinLo)>2.0){
    // bkgFail = new CQuadratic(m,kFALSE,ibin,-20.0,20.0,16.0,5.0,-0.1,0.1); // bin0
    // std::cout << "setting bins..." << std::endl;
    // nflfail += 3;
  // }
if(yaxislabel.CompareTo("stand-alone")==0       && bkgfail==8 && charge==0 && ibin==0){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-20.0,50.0,16.0,10.0,-0.1,10.0); // bin0 amcnlo 2,6
    // bkgFail = new CQuadratic(m,kFALSE,ibin,-70.0,20.0,13.0,5.0,-0.1,10.0); // bin0 minlo 2,6
    std::cout << "setting bins..." << std::endl;
    nflfail += 3;
} else if (yaxislabel.CompareTo("stand-alone")==0       && bkgfail==6 && charge==0) {
    std::cout << "open" << std::endl;
    char templatename[200];
    sprintf(templatename,"/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_5TeV/results/Zmm/Data/MuStaEff_aMCxPythia_BKG_v2/Combined/plots/etapt_%d.root",ibin);
    // sprintf(templatename,"/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zmm/Data/MuStaEff_aMCxPythia_BKG/Combined/plots/etapt_%d.root",ibin);
    // sprintf(templatename,"/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zmm/Data/MuStaEff_aMCxPythia_tagPt_BKG/Combined/plots/etapt_%d.root",ibin);

    TFile *f = new TFile(templatename);
    std::cout << "opeend" << std::endl;
    RooWorkspace *w = (RooWorkspace*) f->Get("w");
    w->Print();
    std::cout << "wksp" << std::endl;
    sprintf(templatename,"backgroundFail_%d",ibin);
    std::cout << templatename << std::endl;
    // RooAbsPdf *bkgFail = w->pdf(templatename);
    // bkgFail->Print();
    // bkgFail->model = (RooGenericPdf*)w->pdf(templatename);
    std::cout << "bleh" << std::endl;
    // (bkgFail->model)->Print();
    std::cout << "blalblhab"<< std::endl;
    // set the starting values and upper and lower limits to the fit value +- error
    char varname[100];
    sprintf(varname,"a0Fail_%d",ibin);
    double pc0=w->var(varname)->getVal();
    double pcerrl0=fabs(w->var(varname)->getErrorLo());
    double pcerrh0=fabs(w->var(varname)->getErrorHi());
    
    std::cout << "a0"<< std::endl;
    sprintf(varname,"a1Fail_%d",ibin);
    double pc1=w->var(varname)->getVal();
    double pcerrl1=fabs(w->var(varname)->getErrorLo());
    double pcerrh1=fabs(w->var(varname)->getErrorHi());
    
    std::cout << "a1"<< std::endl;
    sprintf(varname,"a2Fail_%d",ibin);
    double pc2=w->var(varname)->getVal();
    double pcerrl2=fabs(w->var(varname)->getErrorLo());
    double pcerrh2=fabs(w->var(varname)->getErrorHi());
    std::cout << "a2"<< std::endl;
    
    bkgFail = new CQuadratic(m,kFALSE,ibin,pc0,pcerrh0,pc1,pcerrh1,pc2,pcerrh1);
 } else if (yaxislabel.CompareTo("stand-alone")==0       && bkgfail==7 && charge==0) {
    std::cout << "open" << std::endl;
    char templatename[200];
    sprintf(templatename,"/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zmm/Data/MuStaEff_POWBKG_BKG/Combined/plots/etapt_%d.root",ibin);
    // sprintf(templatename,"../LowPU2017ID_13TeV_v1/results/Zmm/Data/MuStaEff_POWBKG_v2_BkgFailOnly/Combined/plots/etapt_%d.root",ibin);
    TFile *f = new TFile(templatename);
    std::cout << "opeend" << std::endl;
    RooWorkspace *w = (RooWorkspace*) f->Get("w");
    w->Print();
    std::cout << "wksp" << std::endl;
    sprintf(templatename,"backgroundFail_%d",ibin);
    std::cout << templatename << std::endl;
    // RooAbsPdf *bkgFail = w->pdf(templatename);
    // bkgFail->Print();
    // bkgFail->model = (RooGenericPdf*)w->pdf(templatename);
    std::cout << "bleh" << std::endl;
    // (bkgFail->model)->Print();
    std::cout << "blalblhab"<< std::endl;
    // set the starting values and upper and lower limits to the fit value +- error
    char varname[100];
    sprintf(varname,"tFail_%d",ibin);
    double pc0=w->var(varname)->getVal();
    double pcerrl0=fabs(w->var(varname)->getErrorLo());
    double pcerrh0=fabs(w->var(varname)->getErrorHi());
    
    bkgFail = new CPower(m,kFALSE,ibin,pc0,pcerrh0);
 
  } else if(bkgfail==1) { 
    bkgFail = new CExponential(m,kFALSE,ibin);
    nflfail += 1;
  
  } else if(bkgfail==2) {
    bkgFail = new CErfExpo(m,kFALSE,ibin); 
    nflfail += 3;
  
  } else if(bkgfail==3) {
    bkgFail = new CDoubleExp(m,kFALSE);
    nflfail += 3;
    
  } else if(bkgfail==4) {
    bkgFail = new CLinearExp(m,kFALSE);
    nflfail += 2;
  
  } else if(bkgfail==5) {
    bkgFail = new CQuadraticExp(m,kFALSE,ibin);
    nflfail += 3;  

  } else if(bkgfail==6) {
    bkgFail = new CQuadratic(m,kFALSE,ibin,10.,30.,10.,30.,2.,15.);
    nflfail += 3;

  } else if(bkgpass==7) {
    bkgFail = new CPower(m,kFALSE,ibin,0.,0.);
    nflfail += 1;
  }

   std::cout << "hi " << std::endl;
  // Define free parameters
  Double_t NsigMax     = doBinned ? histPass.Integral()+histFail.Integral() : passTree->GetEntries()+failTree->GetEntries();
  Double_t NbkgFailMax = doBinned ? histFail.Integral() : failTree->GetEntries();
  Double_t NbkgPassMax = doBinned ? histPass.Integral() : passTree->GetEntries();
  RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0.0,1.5*NsigMax);
  RooRealVar eff("eff","Efficiency",0.95,0.0,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.1*NbkgPassMax,0.01,NbkgPassMax);
  if(bkgpass==0) NbkgPass.setVal(0);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax,0.01,NbkgFailMax);
  if(yaxislabel.CompareTo("stand-alone")==0 && fabs(xbinLo) >= 0.9 ) {
      NbkgFail.setVal(0.0);
      NbkgFail.setRange(0.0,NbkgFailMax*2.0);
      eff.setVal(0.98);
  }
   else if(yaxislabel.CompareTo("stand-alone")==0       && bkgfail==6 && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==25 && ybinHi==40  ){
       NbkgFail.setVal(0.5*NbkgFailMax);
   }
  std::cout << NbkgFail.getVal() << std::endl;
    std::cout << NbkgFail.getVal() << std::endl;
std::cout << "-----------------" << std::endl;
 

  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
  RooAddPdf *modelPass=0, *modelFail=0;
  RooExtendPdf *esignalPass=0, *ebackgroundPass=0, *esignalFail=0, *ebackgroundFail=0;

  if(massLo!=fitMassLo || massHi!=fitMassHi) {
    m.setRange("signalRange",massLo,massHi);
    
    esignalPass     = new RooExtendPdf("esignalPass","esignalPass",*(sigPass->model),NsigPass,"signalRange");
    ebackgroundPass = new RooExtendPdf("ebackgroundPass","ebackgroundPass",(bkgpass>0) ? *(bkgPass->model) : *(sigPass->model),NbkgPass,"signalRange");
    modelPass       = new RooAddPdf("modelPass","Model for PASS sample",(bkgpass>0) ? RooArgList(*esignalPass,*ebackgroundPass) : RooArgList(*esignalPass));    

    esignalFail     = new RooExtendPdf("esignalFail","esignalFail",*(sigFail->model),NsigFail,"signalRange");
    ebackgroundFail = new RooExtendPdf("ebackgroundFail","ebackgroundFail",*(bkgFail->model),NbkgFail,"signalRange");
    modelFail       = new RooAddPdf("modelFail","Model for FAIL sample", (bkgfail>0) ? RooArgList(*esignalFail,*ebackgroundFail) : RooArgList(*esignalFail));
  
  } else {

    modelPass = new RooAddPdf("modelPass","Model for PASS sample",
                              (bkgpass>0) ? RooArgList(*(sigPass->model),*(bkgPass->model)) :  RooArgList(*(sigPass->model)),
		              (bkgpass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));

    modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigFail->model),*(bkgFail->model)),RooArgList(NsigFail,NbkgFail));
  }
  std::cout << "here" << std::endl;
  
  // eff.setVal(1.0);
  // eff.setConstant(kTRUE);
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");  
  totalPdf.addPdf(*modelFail,"Fail");

  int strategy = 2;

  

// set some of the bins by hand
    strategy=1;

    int nTries = 0;
    
    if(yaxislabel.CompareTo("GSF")==0 && charge==-1&& (ibin ==18)){
      strategy=1;
      nTries=10;
  }

  RooFitResult *fitResult=0;
  RooMsgService::instance().setSilentMode(kTRUE);
  fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),
  			     RooFit::Strategy(strategy), RooFit::Minos(RooArgSet(eff)),RooFit::Save());
                 
   // if((fabs(eff.getErrorLo())<5e-4) || (eff.getErrorHi()<5e-4) || fitResult->status()!=0){
       // fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),RooFit::Strategy(1), RooFit::Minos(RooArgSet(eff)),RooFit::Save());
   // }
    do {
    fitResult = totalPdf.fitTo(*dataCombined,
				 NumCPU(4),
				 //				 Minimizer("Minuit2","minimize"),
				 Minimizer("Minuit2","scan"),
				 // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Minos(),
				 RooFit::Strategy(2),
				 RooFit::Save());
                 
        // nTries++;

    fitResult = totalPdf.fitTo(*dataCombined,
				 NumCPU(4),
				 //				 Minimizer("Minuit2","minimize"),
				 Minimizer("Minuit2","migrad"),
				 // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Hesse(),
				 RooFit::Strategy(2),
				 RooFit::Save());
    fitResult = totalPdf.fitTo(*dataCombined,
				 NumCPU(4),
				 //				 Minimizer("Minuit2","minimize"),
				 Minimizer("Minuit2","improve"),
				 // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 RooFit::Minos(),
				 RooFit::Strategy(2),
				 RooFit::Save());
				 
	fitResult = totalPdf.fitTo(*dataCombined,
			       NumCPU(4),
			       Minimizer("Minuit2","minimize"),
			       // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       RooFit::Minos(),
			       RooFit::Strategy(2),
	               RooFit::Save());
        nTries++;
    }while((fitResult->status()>0)&&nTries < 1);
                 
  // Refit w/o MINOS if MINOS errors are strange...
  if((fabs(eff.getErrorLo())<5e-4) || (eff.getErrorHi()<5e-4) || fitResult->status()!=0){
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","scan"),RooFit::Strategy(2), RooFit::Save());
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","migrad"),RooFit::Strategy(2), RooFit::Save());
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","improve"),RooFit::Strategy(2), RooFit::Save());
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),Minimizer("Minuit2","minimize"),RooFit::Strategy(2), RooFit::Save());
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(),RooFit::Strategy(2),RooFit::Minos(RooArgSet(eff)), RooFit::Save());
  }
  
  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();
    
  if(name.CompareTo("pt")==0) {
    sprintf(binlabelx,"%i GeV/c < p_{T} < %i GeV/c",Int_t(xbinLo),Int_t(xbinHi));
  
  } else if(name.CompareTo("eta")==0) { 
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else	 sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);
  
  } else if(name.CompareTo("phi")==0) { 
    sprintf(binlabelx,"%.1f < #phi < %.1f",xbinLo,xbinHi); 
  
  } else if(name.CompareTo("etapt")==0) {
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);    
    sprintf(binlabely,"%i GeV/c < p_{T} < %i GeV/c",Int_t(ybinLo),Int_t(ybinHi));
  
  } else if(name.CompareTo("etaphi")==0) {
    if(doAbsEta) sprintf(binlabelx,"%.1f < |#eta| < %.1f",xbinLo,xbinHi);
    else         sprintf(binlabelx,"%.1f < #eta < %.1f",xbinLo,xbinHi);					   
    sprintf(binlabely,"%.1f < #phi < %.1f",ybinLo,ybinHi);
  
  } else if(name.CompareTo("npv")==0) { 
    sprintf(binlabelx,"%i #leq N_{PV} < %i",(Int_t)xbinLo,(Int_t)xbinHi); 
  
  } 
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",eff.getVal(),fabs(eff.getErrorLo()),eff.getErrorHi());

  RooPlot *mframePass = m.frame(Bins(Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS));
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
  if(bkgpass>0)
    modelPass->plotOn(mframePass,Components(bkgpassname/*"backgroundPass"*/),LineStyle(kDashed),LineColor(kRed));
  modelPass->plotOn(mframePass);

  RooPlot *mframeFail = m.frame(Bins(Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL));
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail,Components(bkgfailname/*"backgroundFail"*/),LineStyle(kDashed),LineColor(kRed));
  modelFail->plotOn(mframeFail);
  
  char lumitext[100]; // lumi label
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);    

  //
  // Plot passing probes
  //
  double a = NsigPass.getVal(), b = NbkgPass.getVal();
  sprintf(pname,"pass%s_%i",name.Data(),ibin);
  sprintf(yield,"%u Events",(Int_t)passTree->GetEntries());
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_PASS);
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigPass.getVal(),NsigPass.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframePass->chiSquare(nflpass));
  if(bkgpass>0)
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgPass.getVal(),NbkgPass.getPropagatedError(*fitResult));
  CPlot plotPass(pname,mframePass,"Passing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  if((name.CompareTo("etapt")==0) || (name.CompareTo("etaphi")==0)) {
    plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);        
    plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);    
  } else {
    plotPass.AddTextBox(yield,0.21,0.81,0.51,0.85,0,kBlack,-1);
  }
  plotPass.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
  if(bkgpass>0) {
    plotPass.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
  } else {
    plotPass.AddTextBox(0.70,0.73,0.94,0.83,0,kBlack,-1,1,nsigstr);
    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
  }
  plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  if(yaxislabel.CompareTo("supercluster")==0 && xbinLo==0.0 && xbinHi==1.4442 && ybinLo==25 && ybinHi==100) plotPass.SetYRange(0,2000);
  else if(yaxislabel.CompareTo("supercluster")==0 && xbinLo==1.566 && xbinHi==2.5 && ybinLo==25 && ybinHi==100) plotPass.SetYRange(0,350);
  plotPass.Draw(cpass,kTRUE,format);

  //
  // Plot failing probes
  //
  double f = NsigFail.getVal(), d = NbkgFail.getVal();
  sprintf(pname,"fail%s_%i",name.Data(),ibin);
  sprintf(yield,"%u Events",(Int_t)failTree->GetEntries());
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_FAIL);
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigFail.getVal(),NsigFail.getPropagatedError(*fitResult));
  sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgFail.getVal(),NbkgFail.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframePass->chiSquare(nflfail));
  CPlot plotFail(pname,mframeFail,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotFail.AddTextBox(binlabelx,0.21,0.75,0.51,0.80,0,kBlack,-1);
  if((name.CompareTo("etapt")==0) || (name.CompareTo("etaphi")==0)) {
    plotFail.AddTextBox(binlabely,0.21,0.70,0.51,0.75,0,kBlack,-1);    
    plotFail.AddTextBox(yield,0.21,0.66,0.51,0.70,0,kBlack,-1);    
  } else {
    plotFail.AddTextBox(yield,0.21,0.81,0.51,0.85,0,kBlack,-1);
  }
  plotFail.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);  
  plotFail.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
  plotFail.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,-1);
  plotFail.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  if(yaxislabel.CompareTo("supercluster")==0 && xbinLo==0.0 && xbinHi==1.4442 && ybinLo==25 && ybinHi==100) plotFail.SetYRange(0,200);
  else if(yaxislabel.CompareTo("supercluster")==0 && xbinLo==1.566 && xbinHi==2.5 && ybinLo==25 && ybinHi==100) plotFail.SetYRange(0,100);
  plotFail.Draw(cfail,kTRUE,format);  

  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/fitres%s_%i.txt",CPlot::sOutDir.Data(),name.Data(),ibin);
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResult);
  txtfile.close();

  char outFileName[100];
  sprintf(outFileName, "%s/%s_%i.root",CPlot::sOutDir.Data(),name.Data(),ibin);

  RooWorkspace *w = new RooWorkspace("w", "workspace");
  w->import(totalPdf);
  w->import(*fitResult);

  w->writeToFile(outFileName);

    TFile *fw = new TFile(outFileName, "UPDATE");
    TTree *t = new TTree("Bin", "Bin");

    UInt_t nEvents, nBkgFail, nBkgPass, iBin, absEta;
    Float_t ptLo, ptHi;
    Float_t etaLo, etaHi;
    Float_t phiLo, phiHi;
    Float_t npvLo, npvHi;

    t->Branch("nEvents", &nEvents, "nEvents/i");
    t->Branch("nBkgFail", &nBkgFail, "nBkgFail/i");
    t->Branch("nBkgPass", &nBkgPass, "nBkgPass/i");
    t->Branch("iBin", &iBin, "iBin/i");
    t->Branch("absEta", &absEta, "absEta/i");
    t->Branch("ptLo",  &ptLo, "ptLo/F");
    t->Branch("ptHi",  &ptHi, "ptHi/F");
    t->Branch("etaLo",&etaLo, "etaLo/F");
    t->Branch("etaHi",&etaHi, "etaHi/F");
    t->Branch("phiLo",&phiLo, "phiLo/F");
    t->Branch("phiHi",&phiHi, "phiHi/F");
    t->Branch("npvLo",&npvLo, "npvLo/F");
    t->Branch("npvHi",&npvHi, "npvHi/F");

    nEvents = NsigMax;
    nBkgFail = NbkgFailMax;
    nBkgPass = NbkgPassMax;
    iBin = ibin;
    if (doAbsEta) absEta=1;
    else absEta=0;

    ptLo = 999; ptHi = 999; etaLo = 999; etaHi = 999;
    phiLo = 999; phiHi = 999; npvLo = 999; npvHi = 999;

    if(name.CompareTo("pt")==0) {
      ptLo = xbinLo;
      ptHi = xbinHi;

    } else if(name.CompareTo("eta")==0) {
      etaLo = xbinLo;
      etaHi = xbinHi;

    } else if(name.CompareTo("phi")==0) {
      phiLo = xbinLo;
      phiHi = xbinHi;

    } else if(name.CompareTo("etapt")==0) {
      etaLo = xbinLo;
      etaHi = xbinHi;
      ptLo = ybinLo;
      ptHi = ybinHi;

    } else if(name.CompareTo("etaphi")==0) {
      etaLo = xbinLo;
      etaHi = xbinHi;
      phiLo = ybinLo;
      phiHi = ybinHi;

    } else if(name.CompareTo("npv")==0) {
      npvLo = xbinLo;
      npvHi = xbinHi;

    }

    t->Fill();
    fw->Write();
    fw->Close();


  //
  // Clean up
  //
  delete esignalPass;
  delete ebackgroundPass;
  delete esignalFail;
  delete ebackgroundFail;
  delete modelPass;
  delete modelFail;  
  delete dataCombined;
  delete dataPass;
  delete dataFail;
  delete sigPass;
  delete bkgPass;
  delete sigFail;
  delete bkgFail;        
  delete histfile;
  delete datfile;   
  delete sigmaFail;

}

//--------------------------------------------------------------------------------------------------
void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList *parlist = res->correlation("eff");
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist->getSize(); i++) {
    for(Int_t j=0; j<parlist->getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void parseFitResults(ifstream &ifs, double &eff, double &errl, double &errh)
{
  string line;
  while(getline(ifs,line)) {
    size_t found = line.find("eff");
    if(found!=string::npos) {
      found = line.find("+/-");
      if(found!=string::npos) {
        string varname, initval, pmstr;
        stringstream ss(line);
        ss >> varname >> initval >> eff >> pmstr >> errl;
        errh = errl;
        
      } else {
        string varname, initval, errstr;
        stringstream ss(line);
        ss >> varname >> initval >> eff >> errstr;
        size_t ipos = errstr.find(","); 	
        string errlstr = errstr.substr(2,ipos-2);
        string errhstr = errstr.substr(ipos+2,errstr.length()-ipos-1);
        errl = atof(errlstr.c_str());
        errh = atof(errhstr.c_str());
      }
    }
  }
}
