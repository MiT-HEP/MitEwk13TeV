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

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

//=== MAIN MACRO ================================================================================================= 

void plotZmm(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmm");
  gStyle->SetTitleOffset(1.100,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  const Double_t mu_MASS  = 0.1057;
  //
  // input ntuple file names
  //
  enum { eData, eZmm, eEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;

  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_final_read/Zmumu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zmumu/ntuples/zmm_select.raw.root");   typev.push_back(eZmm);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zmumu/ntuples/ewk_select.raw.root");  typev.push_back(eEWK);
  fnamev.push_back("/data/blue/Bacon/Run2/wz_flat_app4/Zmumu/ntuples/top_select.raw.root");  typev.push_back(eEWK);

  //
  // Fit options
  //
  const Int_t    NBINS     = 60;
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;  
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;

  // efficiency files
  const TString dataHLTEffName_pos = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuHLTEff/eff.root";
  const TString dataHLTEffName_neg = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuHLTEff/eff.root";
  const TString zmmHLTEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuHLTEff/eff.root";
  const TString zmmHLTEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuHLTEff/eff.root";

  const TString dataSelEffName_pos = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuSelStaEff_iso/eff.root";
  const TString dataSelEffName_neg = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuSelStaEff_iso/eff.root";
  const TString zmmSelEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuSelStaEff_iso/eff.root";
  const TString zmmSelEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuSelStaEff_iso/eff.root";

  const TString dataTrkEffName_pos = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuTrkEff/eff.root";
  const TString dataTrkEffName_neg = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuTrkEff/eff.root";
  const TString zmmTrkEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuTrkEff/eff.root";
  const TString zmmTrkEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuTrkEff/eff.root";

  const TString dataStaEffName_pos = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuStaEff_iso/eff.root";
  const TString dataStaEffName_neg = "/data/blue/cmedlock/wz-efficiency-results/DataZmm_MuStaEff_iso/eff.root";
  const TString zmmStaEffName_pos  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuStaEff_iso/eff.root";
  const TString zmmStaEffName_neg  = "/data/blue/cmedlock/wz-efficiency-results/Zmm_MuStaEff_iso/eff.root";

  
  // plot output file format
  const TString format("png");

  TString pufname = "../Tools/pileup_weights_2015B.root";

  Int_t yield = 0;
  Double_t yield_zmm = 0, yield_zmm_unc=0;
  Double_t yield_ewk = 0, yield_ewk_unc=0;
   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // event category enumeration
  enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk }; // event category enum


  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_weights_2015B.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("npv_rw");

    
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  // histograms for full selection (MuMu2HLT + MuMu1HLT)
  TH1D *hData = new TH1D("hData","",NBINS,MASS_LOW,MASS_HIGH); hData->Sumw2();
  TH1D *hZmm  = new TH1D("hZmm", "",NBINS,MASS_LOW,MASS_HIGH); hZmm->Sumw2();
  TH1D *hEWK  = new TH1D("hEWK", "",NBINS,MASS_LOW,MASS_HIGH); hEWK->Sumw2();
  TH1D *hMC   = new TH1D("hMC",  "",NBINS,MASS_LOW,MASS_HIGH); hMC->Sumw2();
    
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, u1, u2;
  Float_t puWeight;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  Float_t pfCombIso1, pfCombIso2;

  
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

  TFile *infile=0;
  TTree *intree=0;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);	    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",     &runNum);      // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);      // event number
    intree->SetBranchAddress("matchGen",   &matchGen);    // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category",   &category);    // dilepton category
    intree->SetBranchAddress("npv",        &npv);	  // number of primary vertices
    intree->SetBranchAddress("npu",        &npu);	  // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",     &genVPt);      // GEN Z boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",    &genVPhi);     // GEN Z boson phi (signal MC)
    intree->SetBranchAddress("genVy",      &genVy);       // GEN Z boson rapidity (signal MC)
    intree->SetBranchAddress("genVMass",   &genVMass);    // GEN Z boson mass (signal MC)
    intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",        &met);	  // MET
    intree->SetBranchAddress("metPhi",     &metPhi);      // phi(MET)
    intree->SetBranchAddress("sumEt",      &sumEt);       // Sum ET
    intree->SetBranchAddress("u1",         &u1);	  // parallel component of recoil
    intree->SetBranchAddress("u2",         &u2);	  // perpendicular component of recoil
    intree->SetBranchAddress("q1",         &q1);	  // charge of tag lepton
    intree->SetBranchAddress("q2",         &q2);	  // charge of probe lepton
    intree->SetBranchAddress("dilep",      &dilep);       // dilepton 4-vector
    intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
    intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
    intree->SetBranchAddress("pfCombIso1", &pfCombIso1);  // combined PF isolation of tag lepton
    intree->SetBranchAddress("pfCombIso2", &pfCombIso2);  // combined PF isolation of probe lepton
    intree->SetBranchAddress("puWeight",      &puWeight);        // pu weight
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
   
      //if(dilep->M()        < MASS_LOW)  continue;
      //if(dilep->M()        > MASS_HIGH) continue;
      //if(lep1->Pt()        < PT_CUT)    continue;
      //if(lep2->Pt()        < PT_CUT)    continue;
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;
      //if(q1*q2<0) continue;
      
      Float_t mass = dilep->M();
      
      Double_t weight=1;
      if(typev[ifile]!=eData) {
	weight *= scale1fb*lumi;
	weight*= h_rw->GetBinContent(npv+1);
	//weight *= getMuScaleCorr(lep1->Eta(),0)*getMuScaleCorr(lep2->Eta(),0);
      }
      
      // fill Z events passing selection (MuMu2HLT + MuMu1HLT)
      if((category==eMuMu2HLT) || (category==eMuMu1HLT) || (category==eMuMu1HLT1L1)) {
	//if((category==eMuMu2HLT) || (category==eMuMu1HLT)) {
        if(typev[ifile]==eData) { 
	  if(dilep->M()        < MASS_LOW)  continue;
	  if(dilep->M()        > MASS_HIGH) continue;
	  if(lep1->Pt()        < PT_CUT)    continue;
	  if(lep2->Pt()        < PT_CUT)    continue;
	  hData->Fill(mass); 
	  
	  yield++;
	
	} else {
	  Double_t lp1 = gRandom->Gaus(lep1->Pt()*getMuScaleCorr(lep1->Eta(),0), getMuResCorr(lep1->Eta(),0));
	  Double_t lp2 = gRandom->Gaus(lep2->Pt()*getMuScaleCorr(lep2->Eta(),0), getMuResCorr(lep2->Eta(),0));
	  //Double_t lp1 = lep1->Pt();
	  //Double_t lp2 = lep2->Pt();
	  TLorentzVector l1, l2;
	  l1.SetPtEtaPhiM(lp1,lep1->Eta(),lep1->Phi(),mu_MASS);
	  l2.SetPtEtaPhiM(lp2,lep2->Eta(),lep2->Phi(),mu_MASS);
	  double mll=(l1+l2).M();
	  Double_t effdata, effmc;
	  Double_t corr=1;
	  if(mll       < MASS_LOW)  continue;
	  if(mll       > MASS_HIGH) continue;
	  if(lp1        < PT_CUT)    continue;
	  if(lp2        < PT_CUT)    continue;
	  effdata=1; effmc=1;    
          if(q1>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff((lep1->Eta()), lep1->Pt())); 
            effmc   *= (1.-zmmHLTEff_pos.getEff((lep1->Eta()), lep1->Pt())); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff((lep1->Eta()), lep1->Pt())); 
            effmc   *= (1.-zmmHLTEff_neg.getEff((lep1->Eta()), lep1->Pt())); 
          }
          if(q2>0) {
            effdata *= (1.-dataHLTEff_pos.getEff((lep2->Eta()), lep2->Pt())); 
            effmc   *= (1.-zmmHLTEff_pos.getEff((lep2->Eta()), lep2->Pt()));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff((lep2->Eta()), lep2->Pt())); 
            effmc   *= (1.-zmmHLTEff_neg.getEff((lep2->Eta()), lep2->Pt()));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corr *= effdata/effmc;
    
          effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= dataSelEff_pos.getEff((lep1->Eta()), lep1->Pt()); 
            effmc   *= zmmSelEff_pos.getEff((lep1->Eta()), lep1->Pt()); 
          } else {
            effdata *= dataSelEff_neg.getEff((lep1->Eta()), lep1->Pt()); 
            effmc   *= zmmSelEff_neg.getEff((lep1->Eta()), lep1->Pt()); 
          }
          if(q2>0) {
            effdata *= dataSelEff_pos.getEff((lep2->Eta()), lep2->Pt()); 
            effmc   *= zmmSelEff_pos.getEff((lep2->Eta()), lep2->Pt());
          } else {
            effdata *= dataSelEff_neg.getEff((lep2->Eta()), lep2->Pt()); 
            effmc   *= zmmSelEff_neg.getEff((lep2->Eta()), lep2->Pt());
          }
          corr *= effdata/effmc;
    
	  effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= dataStaEff_pos.getEff((lep1->Eta()), lep1->Pt()); 
            effmc   *= zmmStaEff_pos.getEff((lep1->Eta()), lep1->Pt()); 
          } else {
            effdata *= dataStaEff_neg.getEff((lep1->Eta()), lep1->Pt()); 
            effmc   *= zmmStaEff_neg.getEff((lep1->Eta()), lep1->Pt()); 
          }
          if(q2>0) {
            effdata *= dataStaEff_pos.getEff((lep2->Eta()), lep2->Pt()); 
            effmc   *= zmmStaEff_pos.getEff((lep2->Eta()), lep2->Pt());
          } else {
            effdata *= dataStaEff_neg.getEff((lep2->Eta()), lep2->Pt()); 
            effmc   *= zmmStaEff_neg.getEff((lep2->Eta()), lep2->Pt());
          }
	  //corr *= effdata/effmc; 
	  
          effdata=1; effmc=1;
          if(q1>0) { 
            effdata *= dataTrkEff_pos.getEff((lep1->Eta()), lep1->Pt()); 
            effmc   *= zmmTrkEff_pos.getEff((lep1->Eta()), lep1->Pt()); 
          } else {
            effdata *= dataTrkEff_neg.getEff((lep1->Eta()), lep1->Pt()); 
            effmc   *= zmmTrkEff_neg.getEff((lep1->Eta()), lep1->Pt()); 
          }
          if(q2>0) {
            effdata *= dataTrkEff_pos.getEff((lep2->Eta()), lep2->Pt()); 
            effmc   *= zmmTrkEff_pos.getEff((lep2->Eta()), lep2->Pt());
          } else {
            effdata *= dataTrkEff_neg.getEff((lep2->Eta()), lep2->Pt()); 
            effmc   *= zmmTrkEff_neg.getEff((lep2->Eta()), lep2->Pt());
          }
          corr *= effdata/effmc;
	  corr=1;
	  //TLorentzVector slep1 = (*lep1);
	  //slep1 *= gRandom->Gaus(slep1.Pt(), getMuResCorr(lep1->Eta(),0))/slep1.Pt();
	  
	  //TLorentzVector slep2 = (*lep2);
	  //slep2 *= gRandom->Gaus(slep2.Pt(), getMuResCorr(lep2->Eta(),0))/slep2.Pt();
	  
	  //mass = (slep1+slep2).M();
	  mass = (l1+l2).M();	
	  
	  if(typev[ifile]==eZmm) 
	    {
	      yield_zmm += weight*corr;
	      yield_zmm_unc += weight*weight*corr*corr;
	      //hZmm->Fill(mass,weight*(1.0402)*corr); 
	      //hMC->Fill(mass,weight*(1.0402)*corr);
	      hZmm->Fill(mass,weight*(1.0)*corr); 
	      hMC->Fill(mass,weight*(1.0)*corr);
	    }
	  if(typev[ifile]==eEWK) 
	    {
	      yield_ewk += weight*corr;
	      yield_ewk_unc += weight*weight*corr*corr;
	      hEWK->Fill(mass,weight*corr); 
	      hMC->Fill(mass,weight*corr);
	    }
	}
      }
    }
    
    delete infile;
    infile=0, intree=0;
  } 

  hZmm->Scale((hData->Integral()-hEWK->Integral())/hZmm->Integral());

  TH1D *hZmumuDiff = makeDiffHist(hData,hMC,"hZmumuDiff");
  hZmumuDiff->SetMarkerStyle(kFullCircle); 
  hZmumuDiff->SetMarkerSize(0.9);
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi*1000.);
  else         sprintf(lumitext,"%.0f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);  
  
  // plot colors
  Int_t linecolorZ   = kOrange-3;
  Int_t fillcolorZ   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t ratioColor   = kGray+2;

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
  TGaxis::SetMaxDigits(3);
  
  //
  // MuMu2HLT + MuMu1HLT categories
  //   
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hData->GetBinWidth(1));
  CPlot plotZmumu("zmm","","",ylabel);
  plotZmumu.AddHist1D(hData,"data","E");
  plotZmumu.AddToStack(hZmm,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  //plotZmumu.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  //plotZmumu.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotZmumu.AddTextBox("CMS Preliminary",0.55,0.80,0.90,0.86,0);
  plotZmumu.AddTextBox(lumitext,0.60,0.91,0.92,0.98,0);
  plotZmumu.SetYRange(0.01,1.2*(hData->GetMaximum() + sqrt(hData->GetMaximum())));
  plotZmumu.TransLegend(-0.35,-0.15);
  plotZmumu.Draw(c,kFALSE,format,1);

  CPlot plotZmumuDiff("zmm","","M(#mu^{+}#mu^{-}) [GeV/c^{2}]","#frac{Data-Pred}{Data}");
  plotZmumuDiff.AddHist1D(hZmumuDiff,"EX0",ratioColor);
  plotZmumuDiff.SetYRange(-0.2,0.2);
  plotZmumuDiff.AddLine(MASS_LOW, 0,MASS_HIGH, 0,kBlack,1);
  plotZmumuDiff.AddLine(MASS_LOW, 0.1,MASS_HIGH, 0.1,kBlack,3);
  plotZmumuDiff.AddLine(MASS_LOW,-0.1,MASS_HIGH,-0.1,kBlack,3);
  plotZmumuDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumu2("zmmlog","","",ylabel);
  plotZmumu2.AddHist1D(hData,"data","E");
  plotZmumu2.AddToStack(hEWK,"EWK",fillcolorEWK,linecolorEWK);
  plotZmumu2.AddToStack(hZmm,"Z#rightarrow#mu#mu",fillcolorZ,linecolorZ);
  plotZmumu2.AddTextBox(lumitext,0.60,0.91,0.92,0.98,0);plotZmumu2.SetName("zmmlog");
  plotZmumu2.AddTextBox("CMS Preliminary",0.55,0.80,0.90,0.86,0);
  plotZmumu2.SetLogy();
  plotZmumu2.SetYRange(1e-4*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plotZmumu2.TransLegend(-0.35,-0.15);
  plotZmumu2.Draw(c,kTRUE,format,1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;

  cout << " The Zmm event yield is " << yield << " +/-" << sqrt(yield) << "." << endl;
  cout << " The Zmm expected event yield is " << yield_zmm << " +/-" << sqrt(yield_zmm_unc) << "." << endl;
  cout << " The EWK event yield is " << yield_ewk << " +/-" << sqrt(yield_ewk_unc) << "." << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     

  gBenchmark->Show("plotZmm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
/*TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
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
  }*/
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
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

