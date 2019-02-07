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

#include "CPlot.hh"	            // helper class for plots
#include "MitStyleRemix.hh"         // style settings for drawing
#include "CEffUser1D.hh"            // class for handling efficiency graphs
#include "CEffUser2D.hh"            // class for handling efficiency tables

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
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const TH1D* puWeights); 
void generateDataTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge); 
			   
// Perform count
void performCount(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                  const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		  TTree *passTree, TTree *failTree, const Int_t method, 
		  const TString name, const Double_t massLo, const Double_t massHi, const TString format, const Bool_t doAbsEta,
		  TCanvas *cpass, TCanvas *cfail,const double lumi);

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
  if(sigModPass==2 || sigModFail==2) {
    generateHistTemplates(mcfilename,ptBinEdgesv,etaBinEdgesv,phiBinEdgesv,npvBinEdgesv,fitMassLo,fitMassHi,doAbsEta,charge,puWeights);
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

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    eventTree->GetEntry(ientry);
    if((q)*charge < 0)    continue;
    if(mass < fitMassLo)  continue;
    if(mass > fitMassHi)  continue;
    if(runNum < runNumLo) continue;
    if(runNum > runNumHi) continue;
    wgt  = weight;
    if(doPU>0) wgt *= puWeights->GetBinContent(npu+1);
    
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
        
    if(pass) {
      passTreePtv[ipt]->Fill();
      passTreeEtav[ieta]->Fill();
      passTreePhiv[iphi]->Fill();
      passTreeEtaPtv[ipt*etaNbins + ieta]->Fill();
      passTreeEtaPhiv[iphi*etaNbins + ieta]->Fill();
      if(inpv>=0)      passTreeNPVv[inpv]->Fill();
    } else {
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
      makeEffHist2D(hEffEtaPt, hErrlEtaPt, hErrhEtaPt, passTreeEtaPtv, failTreeEtaPtv, sigModPass, bkgModPass, sigModFail, bkgModFail, "etapt", massLo, massHi, fitMassLo, fitMassHi, format, doAbsEta, lumi, yaxislabel, charge);

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
		           const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const TH1D* puWeights)
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
    
    Double_t puWgt=1;
    if(puWeights)
      puWgt = puWeights->GetBinContent(npu+1);
    
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
  if(xbinLo==1.4442 && xbinHi==1.566) return;
  if(xbinLo==-1.566 && xbinHi==-1.4442) return;

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
void performFit(Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
                const Int_t ibin, const Double_t xbinLo, const Double_t xbinHi, const Double_t ybinLo, const Double_t ybinHi,
		TTree *passTree, TTree *failTree,
		const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail, 
		const TString name, const Double_t massLo, const Double_t massHi, const Double_t fitMassLo, const Double_t fitMassHi,
		const TString format, const Bool_t doAbsEta, TCanvas *cpass, TCanvas *cfail, const double lumi, const TString yaxislabel, const int charge)
{
  // skip ECAL gap region
  if(xbinLo==1.4442 && xbinHi==1.566) return;
  if(xbinLo==-1.566 && xbinHi==-1.4442) return;

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
  if(sigpass==2 || sigfail==2) {
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
  
  } else if(sigpass==2) { 
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


  if(sigfail==1) {
    sigFail = new CBreitWignerConvCrystalBall(m,kFALSE);
    nflfail += 4;
  
  } else if(sigfail==2) {
    char hname[50];
    sprintf(hname,"fail%s_%i",name.Data(),ibin);
    TH1D *h = (TH1D*)histfile->Get(hname);
    assert(h);
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

// only for finebin, must comment out afterwards
// 10-25GeV 3ptbin
// fine 
/*
  if(ibin==0){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-59.783,67.8937,8.3132,1.56861,-0.0596545,0.00861623);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==1){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-9.6466,51.4286,4.29752,1.18684,-0.0322013,0.0065136);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==2){
    bkgFail = new CQuadratic(m,kFALSE,ibin,690.78,139.598,13.2109,3.21826,-0.132308,0.0176632);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==3){
    bkgFail = new CQuadratic(m,kFALSE,ibin,851.445,118.326,1.31549,2.72538,-0.0508294,0.0149547);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==4){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-233.387,42.2047,7.02647,0.99123,-0.0385993,0.00551175);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==5){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-188.181,32.3922,5.34675,0.764585,-0.02936,0.004266);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==6){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-655.256,87.297,22.5921,2.04671,-0.127129,0.011375);
    bkgFail = new CQuadratic(m,kFALSE,ibin,);
    nflfail += 3;
  }else if(ibin==7){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-385.667,73.4868,14.4824,1.72042,-0.0827718,0.009553);
    nflfail += 3;
  }else if(ibin==8){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-52.104,27.0318,1.59234,0.635189,-0.00748997,0.00353835);
    nflfail += 3;
  }else if(ibin==9){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-148.242,19.2788,3.64503,0.461285,-0.0189668,0.00259245);
    nflfail += 3;
  }else if(ibin==10){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-506.523,54.5831,13.4538,1.29179,-0.0684939,0.00722625);
    nflfail += 3;
  }else if(ibin==11){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-293.452,45.2107,7.89477,1.07068,-0.0391149,0.00599699);
    nflfail += 3;
  }else if(ibin==12){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-27.268,17.4549,0.69916,0.414163,-0.00261404,0.00232734);
    nflfail += 3;
  }else if(ibin==13){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-39.7514,12.4144,0.953071,0.2954,-0.00460477,0.00165729);
    nflfail += 3;
  }else if(ibin==14){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-128.869,35.7934,3.37636,0.84775,-0.0142926,0.004754);
    nflfail += 3;
  }else if(ibin==15){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-81.6185,29.4768,2.15664,0.696177,-0.00905715,0.00389467);
    nflfail += 3;
  }else if(ibin==16){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-17.8746,11.5959,0.445289,0.273801,-0.0019419,0.00153015);
    nflfail += 3;
  }else if(ibin==17){
    bkgFail = new CQuadratic(m,kFALSE,ibin,8.31552,8.0714,-0.198305,0.190136,0.00145247,0.00106679);
    nflfail += 3;
  }else if(ibin==18){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-2.63915,25.3908,0.204934,0.597787,0.00122431,0.0033423);
    nflfail += 3;
  }else if(ibin==19){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-12.1681,20.2046,0.353465,0.479098,-0.000258999,0.0026943);
    nflfail += 3;
  }else if(ibin==20){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-12.285,9.98626,0.300958,0.236609,-0.00132876,0.00132831);
    nflfail += 3;
  }else if(ibin==21){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-6.29294,7.65293,0.167095,0.176793,-0.000819758,0.000970331);
    nflfail += 3;
  }else if(ibin==22){
    bkgFail = new CQuadratic(m,kFALSE,ibin,18.1732,19.5609,-0.396395,0.462074,0.00373596,0.00259343);
    nflfail += 3;
  }else if(ibin==23){
    bkgFail = new CQuadratic(m,kFALSE,ibin,14.6117,15.1915,-0.345972,0.35629,0.00302296,0.00198844);
    nflfail += 3;
  }else if(ibin==24){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-3.44895,8.70837,0.0929538,0.198519,-0.000349284,0.00108335);
    nflfail += 3;
  }else if(ibin==25){
    bkgFail = new CQuadratic(m,kFALSE,ibin,3.13921,9.4075,-0.0483007,0.210828,0.000303128,0.00113683);
    nflfail += 3;
  }else if(ibin==26){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-5.28589,13.9806,0.14074,0.331183,0.000134727,0.00186114);
    nflfail += 3;
  }else if(ibin==27){
    bkgFail = new CQuadratic(m,kFALSE,ibin,4.65854,11.4471,-0.11488,0.268156,0.00130256,0.00149423);
    nflfail += 3;
  }else if(ibin==28){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-3.24836,7.40525,0.0995607,0.17071,-0.000511867,0.000935773);
    nflfail += 3;
  }else if(ibin==29){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-10.4309,8.6913,0.271935,0.202471,-0.00147659,0.00110837);
    nflfail += 3;
  }else if(ibin==30){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-6.45202,10.8196,0.201684,0.251387,-0.000869546,0.00138669);
    nflfail += 3;
  }else if(ibin==31){
    bkgFail = new CQuadratic(m,kFALSE,ibin,3.68306,9.6837,-0.0776881,0.227405,0.000794782,0.00126951);
    nflfail += 3;
  }else if(ibin==32){
    bkgFail = new CQuadratic(m,kFALSE,ibin,7.51475,8.18432,-0.151074,0.193261,0.00103473,0.00108676);
    nflfail += 3;
  }else if(ibin==33){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-2.45849,7.04288,0.0874785,0.16379,-0.000438514,0.000905603);
    nflfail += 3;
  }else if(ibin==34){
    bkgFail = new CQuadratic(m,kFALSE,ibin,9.79961,17.8941,-0.17601,0.422011,0.00222819,0.00236496);
    nflfail += 3;
  }else if(ibin==35){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-8.02915,15.7827,0.26603,0.372445,-0.000736144,0.00208482);
    nflfail += 3;
*/




/*
// 10-25 GeV 1ptBin
  if(ibin==0){
    bkgFail = new CQuadratic(m,kFALSE,ibin,1709.76,84.0521,-26.0728,1.91147,0.108192,0.0104126);
    nflfail += 3;
  }else if(ibin==1){
    bkgFail = new CQuadratic(m,kFALSE,ibin,7642.42,164.126,-119.418,3.71331,0.498022,0.0201525);
    nflfail += 3;
  }else if(ibin==2){
    bkgFail = new CQuadratic(m,kFALSE,ibin,3429.28,101.034,-52.8942,2.26829,0.212252,0.0122385);
    nflfail += 3;
  }else if(ibin==3){
    bkgFail = new CQuadratic(m,kFALSE,ibin,9986.87,141.757,-161.386,3.13166,0.666215,0.0167113);
    nflfail += 3;
  }else if(ibin==4){
    bkgFail = new CQuadratic(m,kFALSE,ibin,1936.04,55.1018,-32.099,1.20266,0.134949,0.00636721);
    nflfail += 3;
  }else if(ibin==5){
    bkgFail = new CQuadratic(m,kFALSE,ibin,3750.44,76.0305,-61.9658,1.65552,0.259615,0.00874733);
    nflfail += 3;
  }else if(ibin==6){
    bkgFail = new CQuadratic(m,kFALSE,ibin,4037.26,78.1124,-67.7376,1.70417,0.288245,0.00902468);
    nflfail += 3;
  }else if(ibin==7){
    bkgFail = new CQuadratic(m,kFALSE,ibin,2157.42,57.4862,-36.5161,1.25784,0.157041,0.00667901);
    nflfail += 3;
  }else if(ibin==8){
    bkgFail = new CQuadratic(m,kFALSE,ibin,10113.2,144.163,-162.459,3.18749,0.666608,0.0170186);
    nflfail += 3;
  }else if(ibin==9){
    bkgFail = new CQuadratic(m,kFALSE,ibin,3389.89,101.456,-51.0329,2.27295,0.199101,0.0122369);
    nflfail += 3;
  }else if(ibin==10){
    bkgFail = new CQuadratic(m,kFALSE,ibin,7256.56,169.339,-109.357,3.83773,0.441685,0.0208428);
    nflfail += 3;
  }else if(ibin==11){
    bkgFail = new CQuadratic(m,kFALSE,ibin,1767.34,85.2565,-27.3228,1.9373,0.11478,0.0105444);
    nflfail += 3;
*/
/*
  if(ibin==1){
    bkgFail = new CQuadratic(m,kFALSE,ibin,371.124,75.402,0.20212,1.73545,-0.019815,0.00951763);
    nflfail += 3;
  }else if(ibin==6){
    bkgFail = new CQuadratic(m,kFALSE,ibin,6.4549,48.5684,3.43535,1.12065,-0.02622,0.00615057);
    nflfail += 3;
  }else if(ibin==10){
    bkgFail = new CQuadratic(m,kFALSE,ibin,427.356,76.9267,-1.16924,1.77106,-0.0114582,0.00971668);
    nflfail += 3;
*/

// for alternative iso cone, 0.12
// for nominal iso cone, 0.15
/*  }else*/ if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-2.4 && xbinHi==-2.1 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,14.688,38.4373,1.63443,0.890089,-0.0124674,0.00490208);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-2.1 && xbinHi==-1.2 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,350.161,84.9754,2.20028,1.96228,-0.0320498,0.01079);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-1.2 && xbinHi==-0.9 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,160.725,59.0307,1.1749,1.36233,-0.015835,0.00748612);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-0.9 && xbinHi==-0.3 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,79.7773,92.7792,10.1913,2.14619,-0.0782481,0.0118087);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-0.3 && xbinHi==-0.2 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-88.1536,40.6558,4.38496,0.940991,-0.0285021,0.00517186);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-0.2 && xbinHi== 0.0 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-148.153,53.0005,7.34116,1.22898,-0.0473752,0.00676667);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.0 && xbinHi== 0.2 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-53.0383,54.3729,5.33904,1.2581,-0.0369095,0.00691909);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.2 && xbinHi== 0.3 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-17.7026,41.248,2.85444,0.956444,-0.0203918,0.00527155);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.3 && xbinHi== 0.9 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,76.5751,93.4465,10.5841,2.16069,-0.0814634,0.0118835);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.9 && xbinHi== 1.2 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,69.0969,60.033,3.34439,1.38772,-0.027481,0.00763195);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,334.14,85.6047,2.54751,1.97506,-0.0334427,0.0108505);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 2.1 && xbinHi== 2.4 && ybinLo==25 && ybinHi==40  ){
    bkgFail = new CQuadratic(m,kFALSE,ibin,42.0962,38.015,1.1052,0.879851,-0.0101413,0.00484586);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-2.4 && xbinHi==-2.1 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-50.3944,21.5028,1.342,0.508744,-0.00594804,0.00284998);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-2.1 && xbinHi==-1.2 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-284.695,50.9969,8.17359,1.20505,-0.0406249,0.00674167);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-1.2 && xbinHi==-0.9 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-144.451,35.8615,4.11209,0.846398,-0.0205359,0.0047295);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-0.9 && xbinHi==-0.3 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-341.994,55.7471,9.55779,1.31615,-0.0461556,0.00735736);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-0.3 && xbinHi==-0.2 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-142.544,24.0544,3.53045,0.572755,-0.0175086,0.00321633);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo==-0.2 && xbinHi== 0.0 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-71.9995,32.1319,2.18304,0.760187,-0.00961905,0.00426188);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.0 && xbinHi== 0.2 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-159.72,33.2644,4.2551,0.784839,-0.0209857,0.0043811);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.2 && xbinHi== 0.3 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-111.643,25.1041,2.90704,0.592871,-0.0144992,0.00331016);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.3 && xbinHi== 0.9 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-300.41,57.8766,8.75623,1.36615,-0.0419159,0.00763852);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 0.9 && xbinHi== 1.2 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-169.412,35.74,4.56915,0.845013,-0.0222843,0.0047273);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-226.414,52.1268,6.74502,1.22804,-0.0320565,0.00685594);
    nflfail += 3;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==6 && charge==0 && xbinLo== 2.1 && xbinHi== 2.4 && ybinLo==40 && ybinHi==8000){
    bkgFail = new CQuadratic(m,kFALSE,ibin,-66.9827,21.6348,1.83467,0.509517,-0.00924994,0.00283993);
    nflfail += 3;

  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo==-2.4 && xbinHi==-2.1 && ybinLo==40 && ybinHi==8000) { //12
    bkgFail = new CPower(m,kFALSE,ibin,-0.658153,0.00805732);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo==-2.1 && xbinHi==-1.2 && ybinLo==40 && ybinHi==8000) { //13
    bkgFail = new CPower(m,kFALSE,ibin,-1.02923,0.00342985);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo==-1.2 && xbinHi==-0.9 && ybinLo==40 && ybinHi==8000) { //14
    bkgFail = new CPower(m,kFALSE,ibin,-0.869479,0.00494074);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo==-0.9 && xbinHi==-0.3 && ybinLo==40 && ybinHi==8000) { //15
    bkgFail = new CPower(m,kFALSE,ibin,-1.06908,0.00313672);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo==-0.2 && xbinHi== 0.0 && ybinLo==40 && ybinHi==8000) { //17
    bkgFail = new CPower(m,kFALSE,ibin,-0.833005,0.00538412);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo== 0.0 && xbinHi== 0.2 && ybinLo==40 && ybinHi==8000) { //18
    bkgFail = new CPower(m,kFALSE,ibin,-0.839799,0.00530299);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo== 0.3 && xbinHi== 0.9 && ybinLo==40 && ybinHi==8000) { //20
    bkgFail = new CPower(m,kFALSE,ibin,-1.07879,0.00306544);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo== 0.9 && xbinHi== 1.2 && ybinLo==40 && ybinHi==8000) { //21
    bkgFail = new CPower(m,kFALSE,ibin,-0.874991,0.00489383);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==40 && ybinHi==8000) { //22
    bkgFail = new CPower(m,kFALSE,ibin,-1.03631,0.0033772);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("Sta")==0       && bkgfail==7 && charge==0 && xbinLo== 2.1 && xbinHi== 2.4 && ybinLo==40 && ybinHi==8000) { //23
    bkgFail = new CPower(m,kFALSE,ibin,-0.655257,0.00805995);
    nflfail += 1;

  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo==-2.5    && xbinHi==-2.0    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.749614,0.00661016);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo==-2.0    && xbinHi==-1.566  && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.752473,0.00654215);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo==-1.4442 && xbinHi==-1.0    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.815651,0.00564858);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo==-1.0    && xbinHi==-0.5    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.896797,0.00470852);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo==-0.5    && xbinHi== 0.0    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.93297,0.00432771);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo== 0.0    && xbinHi== 0.5    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.922113,0.00443263);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo== 0.5    && xbinHi== 1.0    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.905108,0.00461378);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo== 1.0    && xbinHi== 1.4442 && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.819914,0.0056123);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo== 1.566  && xbinHi== 2.0    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.756879,0.00649284);
    nflfail += 1;
  } else if(yaxislabel.CompareTo("GsfSel")==0    && bkgfail==7 && charge==0 && xbinLo== 2.0    && xbinHi== 2.5    && ybinLo==55 && ybinHi==8000) {
    bkgFail = new CPower(m,kFALSE,ibin,-0.752264,0.00658988);
    nflfail += 1;

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
    bkgFail = new CQuadratic(m,kFALSE,ibin,0.,0.,0.,0.,0.,0.);
    nflfail += 3;

  } else if(bkgpass==7) {
    bkgFail = new CPower(m,kFALSE,ibin,0.,0.);
    nflfail += 1;
  }


  // Define free parameters
  Double_t NsigMax     = doBinned ? histPass.Integral()+histFail.Integral() : passTree->GetEntries()+failTree->GetEntries();
  Double_t NbkgFailMax = doBinned ? histFail.Integral() : failTree->GetEntries();
  Double_t NbkgPassMax = doBinned ? histPass.Integral() : passTree->GetEntries();
  RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0,1.5*NsigMax);
  RooRealVar eff("eff","Efficiency",0.9,0.0,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.1*NbkgPassMax,0.01,NbkgPassMax);
  if(bkgpass==0) NbkgPass.setVal(0);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax,0.01,NbkgFailMax);


//(2,1,2,2)
//  if(yaxislabel.CompareTo("GsfSel")==0) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}
 
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.5    && xbinHi==-2.0    && ybinLo==25 && ybinHi==  40) {}//0 
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.0    && xbinHi==-1.566  && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//1
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-1.4442 && xbinHi==-1.0    && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//3
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-1.0    && xbinHi==-0.5    && ybinLo==25 && ybinHi==  40) {}//4
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-0.5    && xbinHi== 0.0    && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//5
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.0    && xbinHi== 0.5    && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.75);Nsig.setRange(0,2.0*NsigMax);}//6
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.5    && xbinHi== 1.0    && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.75);Nsig.setRange(0,2.0*NsigMax);}//7
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.0    && xbinHi== 1.4442 && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//8
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.566  && xbinHi== 2.0    && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//10
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 2.0    && xbinHi== 2.5    && ybinLo==25 && ybinHi==  40) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//11
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.5    && xbinHi==-2.0    && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,5.0*NsigMax);}//12 
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.0    && xbinHi==-1.566  && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//13
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-1.4442 && xbinHi==-1.0    && ybinLo==40 && ybinHi==  55) {}//15
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-1.0    && xbinHi==-0.5    && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//16
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-0.5    && xbinHi== 0.0    && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,1.5*NsigMax);}//17
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.0    && xbinHi== 0.5    && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.1*NsigMax);}//18
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.5    && xbinHi== 1.0    && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,1.0*NsigMax);}//19
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.0    && xbinHi== 1.4442 && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.6*NsigMax);}//20
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.566  && xbinHi== 2.0    && ybinLo==40 && ybinHi==  55) {}//22
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 2.0    && xbinHi== 2.5    && ybinLo==40 && ybinHi==  55) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.72);Nsig.setRange(0,1.8*NsigMax);}//23

//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.5    && xbinHi==-2.0    && ybinLo==55 && ybinHi==8000) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//24
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.0    && xbinHi==-1.566  && ybinLo==55 && ybinHi==8000) {}//25
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-1.4442 && xbinHi==-1.0    && ybinLo==55 && ybinHi==8000) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.8);Nsig.setRange(0,2.0*NsigMax);}//27
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-1.0    && xbinHi==-0.5    && ybinLo==55 && ybinHi==8000) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.85);Nsig.setRange(0,2.0*NsigMax);}//28
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-0.5    && xbinHi== 0.0    && ybinLo==55 && ybinHi==8000) {}//29
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.0    && xbinHi== 0.5    && ybinLo==55 && ybinHi==8000) {}//30
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.5    && xbinHi== 1.0    && ybinLo==55 && ybinHi==8000) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.87);Nsig.setRange(0,2.0*NsigMax);}//31
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.0    && xbinHi== 1.4442 && ybinLo==55 && ybinHi==8000) {Nsig.setVal(0.8*NsigMax);eff.setVal(0.85);Nsig.setRange(0,2.0*NsigMax);}//32
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.566  && xbinHi== 2.0    && ybinLo==55 && ybinHi==8000) {}//34
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 2.0    && xbinHi== 2.5    && ybinLo==55 && ybinHi==8000) {}//35

// for nominal iso cone, 0.15
// for 76X
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo== 0.2 && xbinHi== 0.3 && ybinLo==25 && ybinHi==40  ) { Nsig.setRange(0,1.2*NsigMax);}
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo== 0.9 && xbinHi== 1.2 && ybinLo==25 && ybinHi==40  ) { Nsig.setRange(0,2.0*NsigMax);}
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==25 && ybinHi==40  ) { Nsig.setRange(0,1.0*NsigMax);}
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo==-0.9 && xbinHi==-0.3 && ybinLo==40 && ybinHi==8000) { Nsig.setRange(0,2.0*NsigMax);}
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo==-0.2 && xbinHi== 0.0 && ybinLo==40 && ybinHi==8000) { Nsig.setRange(0,1.2*NsigMax);}
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo== 0.9 && xbinHi== 1.2 && ybinLo==40 && ybinHi==8000) { Nsig.setRange(0,2.0*NsigMax);}
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==40 && ybinHi==8000) { Nsig.setRange(0,1.2*NsigMax);}

  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-2.1 && xbinHi==-1.2 && ybinLo==25 && ybinHi==40  ) { Nsig.setRange(0,1.0*NsigMax);}//1
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-0.2 && xbinHi== 0.0 && ybinLo==25 && ybinHi==40  ) { Nsig.setRange(0,1.0*NsigMax);}//5
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo== 0.0 && xbinHi== 0.2 && ybinLo==25 && ybinHi==40  ) { Nsig.setRange(0,1.0*NsigMax);}//6

// for 76X 
  if(yaxislabel.CompareTo("SIT")==0       && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==10 && ybinHi==  25) { eff.setVal(0.80);}
  if(yaxislabel.CompareTo("SIT")==0 	  && charge==0 && xbinLo==-2.4 && xbinHi==-2.1 && ybinLo==40 && ybinHi==8000) { eff.setVal(0.95);}
  if(yaxislabel.CompareTo("SIT")==0       && charge==0 && xbinLo== 0.0 && xbinHi== 0.2 && ybinLo==40 && ybinHi==8000) { Nsig.setRange(0,2.0*NsigMax);}
  if(yaxislabel.CompareTo("SIT")==0       && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==40 && ybinHi==8000) { Nsig.setRange(0,2.0*NsigMax);}

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

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");  
  totalPdf.addPdf(*modelFail,"Fail");

  int strategy = 2;
// for 76X
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo== 0.0 && xbinHi== 2.4 && ybinLo==25 && ybinHi==8000) strategy = 1;
  if(yaxislabel.CompareTo("Sta")==0       && charge==0 && xbinLo==-0.9 && xbinHi==-0.3 && ybinLo==40 && ybinHi==8000) strategy = 1; 
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-2.4 && xbinHi==-2.1 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//0
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-2.1 && xbinHi==-1.2 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//1
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-1.2 && xbinHi==-0.9 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//2
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-0.9 && xbinHi==-0.3 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//3
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo==-0.3 && xbinHi==-0.2 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//4
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo== 0.3 && xbinHi== 0.9 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//8
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo== 0.9 && xbinHi== 1.2 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//9
  if(yaxislabel.CompareTo("Sta")==0       && sigpass==1 && charge==0 && xbinLo== 1.2 && xbinHi== 2.1 && ybinLo==25 && ybinHi==40  ) { strategy = 1;}//10

  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.5    && xbinHi==-2.0    && ybinLo==25 && ybinHi==  40) strategy = 1;//1
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-0.5    && xbinHi== 0.0    && ybinLo==25 && ybinHi==  40) strategy = 1;//5
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.0    && xbinHi== 1.4442 && ybinLo==25 && ybinHi==  40) strategy = 1;//8
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.566  && xbinHi== 2.0    && ybinLo==25 && ybinHi==  40) strategy = 1;//10
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 2.0    && xbinHi== 2.5    && ybinLo==25 && ybinHi==  40) strategy = 1;//11
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 2.0    && xbinHi== 2.5    && ybinLo==40 && ybinHi==  55) strategy = 1;//23
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-2.0    && xbinHi==-1.566  && ybinLo==55 && ybinHi==8000) strategy = 1;//25
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo==-0.5    && xbinHi== 0.0    && ybinLo==55 && ybinHi==8000) strategy = 1;//29
  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 0.0    && xbinHi== 0.5    && ybinLo==55 && ybinHi==8000) strategy = 1;//30
//  if(yaxislabel.CompareTo("GsfSel")==0    && charge==0 && xbinLo== 1.566  && xbinHi== 2.0    && ybinLo==55 && ybinHi==8000) strategy = 1;//34

  RooFitResult *fitResult=0;
  RooMsgService::instance().setSilentMode(kTRUE);
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
			     RooFit::Extended(),
  			     RooFit::Strategy(strategy), // MINOS STRATEGY
  			     RooFit::Minos(RooArgSet(eff)),
  			     RooFit::Save());

    // int nTries = 0;
    // do {
      // fitResult = totalPdf.fitTo(*dataCombined,
				 // NumCPU(4),
				 // //				 Minimizer("Minuit2","minimize"),
				 // Minimizer("Minuit2","scan"),
				 // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // RooFit::Minos(),
				 // RooFit::Strategy(2),
				 // RooFit::Save());
                 
        // // nTries++;

      // fitResult = totalPdf.fitTo(*dataCombined,
				 // NumCPU(4),
				 // //				 Minimizer("Minuit2","minimize"),
				 // Minimizer("Minuit2","migrad"),
				 // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // RooFit::Hesse(),
				 // RooFit::Strategy(2),
				 // RooFit::Save());
       // fitResult = totalPdf.fitTo(*dataCombined,
				 // NumCPU(4),
				 // //				 Minimizer("Minuit2","minimize"),
				 // Minimizer("Minuit2","improve"),
				 // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
				 // RooFit::Minos(),
				 // RooFit::Strategy(2),
				 // RooFit::Save());
				 
	// fitResult = totalPdf.fitTo(*dataCombined,
			       // NumCPU(4),
			       // Minimizer("Minuit2","minimize"),
			       // // ExternalConstraints(constGauss1),ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       // //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
			       // RooFit::Minos(),
			       // RooFit::Strategy(2),
	               // RooFit::Save());
        // nTries++;
    // }while((fitResult->status()>0)&&nTries < 10);
                 
  // Refit w/o MINOS if MINOS errors are strange...
  if((fabs(eff.getErrorLo())<5e-4) || (eff.getErrorHi()<5e-4) || fitResult->status()>0)
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());
  
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
