//================================================================================================
//
// Compute Z->mumu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties no yet implemented
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
// #include <math.h>                  // mathematics
#include "TLorentzVector.h"         // 4-vector class

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions

// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"
#endif

//=== MAIN MACRO ================================================================================================= 

void computeAccSelZmmBinned_Sys(const TString conf,      // input file
			    const TString inputDir,
                const TString outputDir,  // output directory
			    const Int_t   doPU,
			    const TString sysFileSIT, // condense these into 1 file per type of eff (pos & neg into 1 file)
			    const TString sysFileSta,
                const bool is13TeV=1
) {
  gBenchmark->Start("computeAccSelZmmBinned");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t MUON_MASS  = 0.105658369;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 13;
  
  // efficiency files
  // const TString inputDir = "/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results/Zmm/";
  
  const TString dataHLTEffName_pos = inputDir + "Data/MuHLTEff_aMCxPythia/Positive/eff.root";
  const TString dataHLTEffName_neg = inputDir + "Data/MuHLTEff_aMCxPythia/Negative/eff.root";
  const TString zmmHLTEffName_pos  = inputDir + "MC/MuHLTEff_aMCxPythia/Positive/eff.root";
  const TString zmmHLTEffName_neg  = inputDir + "MC/MuHLTEff_aMCxPythia/Negative/eff.root";

  const TString dataSelEffName_pos = inputDir + "Data/MuSITEff_aMCxPythia/Positive/eff.root";
  const TString dataSelEffName_neg = inputDir + "Data/MuSITEff_aMCxPythia/Negative/eff.root";
  const TString zmmSelEffName_pos  = inputDir + "MC/MuSITEff_aMCxPythia/Positive/eff.root";
  const TString zmmSelEffName_neg  = inputDir + "MC/MuSITEff_aMCxPythia/Negative/eff.root";

  const TString dataStaEffName_pos = inputDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString dataStaEffName_neg = inputDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString zmmStaEffName_pos  = inputDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString zmmStaEffName_neg  = inputDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";

  // uncertainty files
  TFile *fSysSIT = TFile::Open(sysFileSIT);
  TFile *fSysSta = TFile::Open(sysFileSta);

  TH2D  *hSysSITSigFSRNeg = (TH2D*) fSysSIT->Get("hMuSITEffSigFSRNeg");
  TH2D  *hSysSITSigFSRPos = (TH2D*) fSysSIT->Get("hMuSITEffSigFSRPos"); 
  TH2D  *hSysSITSigMCNeg  = (TH2D*) fSysSIT->Get("hMuSITEffSigMCNeg"); 
  TH2D  *hSysSITSigMCPos  = (TH2D*) fSysSIT->Get("hMuSITEffSigMCPos"); 
  TH2D  *hSysSITBkgNeg    = (TH2D*) fSysSIT->Get("hMuSITEffBkgNeg"); 
  TH2D  *hSysSITBkgPos    = (TH2D*) fSysSIT->Get("hMuSITEffBkgPos"); 
  
  // the real one but i don't have the files done yet
  TH2D  *hSysStaSigFSRNeg = (TH2D*) fSysSta->Get("hMuStaEffSigFSRNeg");
  TH2D  *hSysStaSigFSRPos = (TH2D*) fSysSta->Get("hMuStaEffSigFSRPos"); 
  TH2D  *hSysStaSigMCNeg  = (TH2D*) fSysSta->Get("hMuStaEffSigMCNeg"); 
  TH2D  *hSysStaSigMCPos  = (TH2D*) fSysSta->Get("hMuStaEffSigMCPos"); 
  TH2D  *hSysStaBkgNeg    = (TH2D*) fSysSta->Get("hMuStaEffBkgNeg"); 
  TH2D  *hSysStaBkgPos    = (TH2D*) fSysSta->Get("hMuStaEffBkgPos"); 
  
  // // placeholder
  // TH2D  *hSysStaSigFSRNeg = (TH2D*) fSysSIT->Get("hMuSITEffSigFSRNeg");
  // TH2D  *hSysStaSigFSRPos = (TH2D*) fSysSIT->Get("hMuSITEffSigFSRPos"); 
  // TH2D  *hSysStaSigMCNeg  = (TH2D*) fSysSIT->Get("hMuSITEffSigMCNeg"); 
  // TH2D  *hSysStaSigMCPos  = (TH2D*) fSysSIT->Get("hMuSITEffSigMCPos"); 
  // TH2D  *hSysStaBkgNeg    = (TH2D*) fSysSIT->Get("hMuSITEffBkgNeg"); 
  // TH2D  *hSysStaBkgPos    = (TH2D*) fSysSIT->Get("hMuSITEffBkgPos"); 


  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
 
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   colorv;  // plot color per input file
  vector<Int_t>   linev;   // plot line style per input file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    ss >> fname >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
  }
  ifs.close();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);  

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

  h =(TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt");
  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
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

  h =(TH2D*)dataSelEffFile_pos->Get("hEffEtaPt");
  TH2D *hSelErr_pos = new TH2D("hSelErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErr_neg = new TH2D("hSelErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

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

  h =(TH2D*)dataStaEffFile_pos->Get("hEffEtaPt");
  TH2D *hStaErr_pos = new TH2D("hStaErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErr_neg = new TH2D("hStaErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			       h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr        = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
   
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nAllv;
  vector<Double_t> nSelCorrv, nSelCorrVarv;
  vector<Double_t> accv, accCorrv;
  vector<Double_t> accErrv, accErrCorrv;
  vector<Double_t> nSelCorrvFSR, nSelCorrvMC, nSelCorrvBkg;
  vector<Double_t> nSelCorrVarvFSR, nSelCorrVarvMC, nSelCorrVarvBkg;
  vector<Double_t> accCorrvFSR, accCorrvMC, accCorrvBkg;
  vector<Double_t> accErrCorrvFSR, accErrCorrvMC, accErrCorrvBkg;
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr  = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch* genPartBr = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Muon",          &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");   
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    nAllv.push_back(0);
    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelCorrv.push_back(0);
    nSelCorrVarv.push_back(0);
    nSelCorrvFSR.push_back(0);
    nSelCorrVarvFSR.push_back(0);
    nSelCorrvMC.push_back(0);
    nSelCorrVarvMC.push_back(0);
    nSelCorrvBkg.push_back(0);
    nSelCorrVarvBkg.push_back(0);

    //
    // loop over events
    //      
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; ientry<(uint)(0.1*eventTree->GetEntries()); ientry++) {
    // for(UInt_t ientry=895100; ientry<895102; ientry++) {
      if(ientry%1000000==0)   cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);
      // if(info->evtNum!=6768167)continue;
      Int_t glepq1=-99;
      Int_t glepq2=-99;
      bool alreadyDid=false;

      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      
      nAllv[ifile]+=gen->weight;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,&glepq1,&glepq2,1);
      if((vec->M()<MASS_LOW || vec->M()>MASS_HIGH)) continue;
      delete vec; delete lep1; delete lep2;

      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;
      
      // trigger requirement               
      // if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
      if (!isMuonTrigger(triggerMenu, info->triggerBits,kFALSE,is13TeV)) continue;

      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);

      for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++) {
  	const baconhep::TMuon *mu1 = (baconhep::TMuon*)((*muonArr)[i1]);

        if(mu1->pt	  < PT_CUT)  continue;  // lepton pT cut
        if(fabs(mu1->eta) > ETA_CUT) continue;  // lepton |eta| cut
        if(!passMuonID(mu1))	     continue;  // lepton selection
	
	TLorentzVector vMu1(0,0,0,0);
	vMu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, MUON_MASS);

        for(Int_t i2=i1+1; i2<muonArr->GetEntriesFast(); i2++) {
          const baconhep::TMuon *mu2 = (baconhep::TMuon*)((*muonArr)[i2]);
        
          if(mu1->q == mu2->q)	       continue;  // opposite charge requirement
          if(mu2->pt        < PT_CUT)  continue;  // lepton pT cut
          if(fabs(mu2->eta) > ETA_CUT) continue;  // lepton |eta| cut
	  if(!passMuonID(mu2))	       continue;  // lepton selection

          TLorentzVector vMu2(0,0,0,0);
	  vMu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, MUON_MASS);  

          // trigger match
	  // if(!isMuonTriggerObj(triggerMenu, mu1->hltMatchBits, kFALSE) && !isMuonTriggerObj(triggerMenu, mu2->hltMatchBits, kFALSE)) continue;
	  if(!isMuonTriggerObj(triggerMenu, mu1->hltMatchBits, kFALSE,is13TeV)&& !isMuonTriggerObj(triggerMenu, mu2->hltMatchBits, kFALSE,is13TeV)) continue;
	  // mass window
          TLorentzVector vDilep = vMu1 + vMu2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          /******** We have a Z candidate! HURRAY! ********/
          Double_t effdata, effmc;
          Double_t effdataFSR, effdataMC, effdataBkg;
	      Double_t corr=1;
	      Double_t corrFSR=1;
	      Double_t corrMC=1;
	      Double_t corrBkg=1;
	  
	      effdata=1; effmc=1;    
	      effdataFSR=1; effdataMC=1; effdataBkg=1;  
          if(mu1->q>0) { 
            effdata *= (1.-dataHLTEff_pos.getEff(mu1->eta, mu1->pt)); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(mu1->eta, mu1->pt)); 
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(mu1->eta, mu1->pt)); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(mu1->eta, mu1->pt)); 
          }
          if(mu2->q>0) {
            effdata *= (1.-dataHLTEff_pos.getEff(mu2->eta, mu2->pt)); 
            effmc   *= (1.-zmmHLTEff_pos.getEff(mu2->eta, mu2->pt));
          } else {
            effdata *= (1.-dataHLTEff_neg.getEff(mu2->eta, mu2->pt)); 
            effmc   *= (1.-zmmHLTEff_neg.getEff(mu2->eta, mu2->pt));
          }
          effdata = 1.-effdata;
          effmc   = 1.-effmc;
          corrFSR *= effdata/effmc; // alternate fsr model
          corrMC *= effdata/effmc; // alternate mc gen model
          corrBkg *= effdata/effmc; // alternate bkg model
          corr *= effdata/effmc; // orig
    
          effdata=1; effmc=1;
          if(mu1->q>0) {
            effdataFSR *= dataSelEff_pos.getEff(mu1->eta, mu1->pt) * hSysSITSigFSRPos->GetBinContent(hSysSITSigFSRPos->GetXaxis()->FindBin(mu1->eta), hSysSITSigFSRPos->GetYaxis()->FindBin(mu1->pt));
            effdataMC  *= dataSelEff_pos.getEff(mu1->eta, mu1->pt) * hSysSITSigMCPos->GetBinContent(hSysSITSigMCPos->GetXaxis()->FindBin(mu1->eta), hSysSITSigMCPos->GetYaxis()->FindBin(mu1->pt));
            effdataBkg *= dataSelEff_pos.getEff(mu1->eta, mu1->pt) * hSysSITBkgPos->GetBinContent(hSysSITBkgPos->GetXaxis()->FindBin(mu1->eta), hSysSITBkgPos->GetYaxis()->FindBin(mu1->pt));
            effdata *= dataSelEff_pos.getEff(mu1->eta, mu1->pt); // original
            // std::cout << "original " << effdata << "  unc " << effdataFSR << "  ratio " << effdata/effdataFSR << std::endl;
            effmc   *= zmmSelEff_pos.getEff(mu1->eta, mu1->pt);
          } else {
            effdataFSR *= dataSelEff_neg.getEff(mu1->eta, mu1->pt) * hSysSITSigFSRNeg->GetBinContent(hSysSITSigFSRNeg->GetXaxis()->FindBin(mu1->eta), hSysSITSigFSRNeg->GetYaxis()->FindBin(mu1->pt));
            effdataMC  *= dataSelEff_neg.getEff(mu1->eta, mu1->pt) * hSysSITSigMCNeg->GetBinContent(hSysSITSigMCNeg->GetXaxis()->FindBin(mu1->eta), hSysSITSigMCNeg->GetYaxis()->FindBin(mu1->pt));
            effdataBkg *= dataSelEff_neg.getEff(mu1->eta, mu1->pt) * hSysSITBkgNeg->GetBinContent(hSysSITBkgNeg->GetXaxis()->FindBin(mu1->eta), hSysSITBkgNeg->GetYaxis()->FindBin(mu1->pt));
            effdata *= dataSelEff_neg.getEff(mu1->eta, mu1->pt);//orig
            // std::cout << "original " << effdata << "  unc " << effdataFSR << std::endl;
            effmc   *= zmmSelEff_neg.getEff(mu1->eta, mu1->pt);
          }
          if(mu2->q>0) {
            effdataFSR *= dataSelEff_pos.getEff(mu2->eta, mu2->pt) * hSysSITSigFSRPos->GetBinContent(hSysSITSigFSRPos->GetXaxis()->FindBin(mu2->eta), hSysSITSigFSRPos->GetYaxis()->FindBin(mu2->pt));
            effdataMC  *= dataSelEff_pos.getEff(mu2->eta, mu2->pt) * hSysSITSigMCPos->GetBinContent(hSysSITSigMCPos->GetXaxis()->FindBin(mu2->eta), hSysSITSigMCPos->GetYaxis()->FindBin(mu2->pt));
            effdataBkg *= dataSelEff_pos.getEff(mu2->eta, mu2->pt) * hSysSITBkgPos->GetBinContent(hSysSITBkgPos->GetXaxis()->FindBin(mu2->eta), hSysSITBkgPos->GetYaxis()->FindBin(mu2->pt));
            effdata *= dataSelEff_pos.getEff(mu2->eta, mu2->pt);// orig
            effmc   *= zmmSelEff_pos.getEff(mu2->eta, mu2->pt);
          } else {
            effdataFSR *= dataSelEff_neg.getEff(mu2->eta, mu2->pt) * hSysSITSigFSRNeg->GetBinContent(hSysSITSigFSRNeg->GetXaxis()->FindBin(mu2->eta), hSysSITSigFSRNeg->GetYaxis()->FindBin(mu2->pt));
            effdataMC  *= dataSelEff_neg.getEff(mu2->eta, mu2->pt) * hSysSITSigMCNeg->GetBinContent(hSysSITSigMCNeg->GetXaxis()->FindBin(mu2->eta), hSysSITSigMCNeg->GetYaxis()->FindBin(mu2->pt));
            effdataBkg *= dataSelEff_neg.getEff(mu2->eta, mu2->pt) * hSysSITBkgNeg->GetBinContent(hSysSITBkgNeg->GetXaxis()->FindBin(mu2->eta), hSysSITBkgNeg->GetYaxis()->FindBin(mu2->pt));
            effdata *= dataSelEff_neg.getEff(mu2->eta, mu2->pt);// orig
            effmc   *= zmmSelEff_neg.getEff(mu2->eta, mu2->pt);
          }
          corrFSR *= effdataFSR/effmc; // alternate fsr model
          corrMC *= effdataMC/effmc; // alternate mc gen model
          corrBkg *= effdataBkg/effmc; // alternate bkg model
          corr *= effdata/effmc; // orig
    
    // std::cout << "corr/corrFSR1  " << (corr/corrFSR-1)*100 << std::endl;
    
          effdata=1; effmc=1;
	     effdataFSR=1; effdataMC=1; effdataBkg=1;  
          if(mu1->q>0) {
            effdataFSR *= dataStaEff_pos.getEff(mu1->eta, mu1->pt) * hSysStaSigFSRPos->GetBinContent(hSysStaSigFSRPos->GetXaxis()->FindBin(mu1->eta), hSysStaSigFSRPos->GetYaxis()->FindBin(mu1->pt));
            effdataMC  *= dataStaEff_pos.getEff(mu1->eta, mu1->pt) * hSysStaSigMCPos->GetBinContent(hSysStaSigMCPos->GetXaxis()->FindBin(mu1->eta), hSysStaSigMCPos->GetYaxis()->FindBin(mu1->pt));
            effdataBkg *= dataStaEff_pos.getEff(mu1->eta, mu1->pt) * hSysStaBkgPos->GetBinContent(hSysStaBkgPos->GetXaxis()->FindBin(mu1->eta), hSysStaBkgPos->GetYaxis()->FindBin(mu1->pt));
            effdata *= dataStaEff_pos.getEff(mu1->eta, mu1->pt);//orig
            effmc   *= zmmStaEff_pos.getEff(mu1->eta, mu1->pt);
          } else {
            effdataFSR *= dataStaEff_neg.getEff(mu1->eta, mu1->pt) * hSysStaSigFSRNeg->GetBinContent(hSysStaSigFSRNeg->GetXaxis()->FindBin(mu1->eta), hSysStaSigFSRNeg->GetYaxis()->FindBin(mu1->pt));
            effdataMC  *= dataStaEff_neg.getEff(mu1->eta, mu1->pt) * hSysStaSigMCNeg->GetBinContent(hSysStaSigMCNeg->GetXaxis()->FindBin(mu1->eta), hSysStaSigMCNeg->GetYaxis()->FindBin(mu1->pt));
            effdataBkg *= dataStaEff_neg.getEff(mu1->eta, mu1->pt) * hSysStaBkgNeg->GetBinContent(hSysStaBkgNeg->GetXaxis()->FindBin(mu1->eta), hSysStaBkgNeg->GetYaxis()->FindBin(mu1->pt));
            effdata *= dataStaEff_neg.getEff(mu1->eta, mu1->pt);//orig
            effmc   *= zmmStaEff_neg.getEff(mu1->eta, mu1->pt);
          }
          if(mu2->q>0) {
            effdataFSR *= dataStaEff_pos.getEff(mu2->eta, mu2->pt) * hSysStaSigFSRPos->GetBinContent(hSysStaSigFSRPos->GetXaxis()->FindBin(mu2->eta), hSysStaSigFSRPos->GetYaxis()->FindBin(mu2->pt));
            effdataMC *= dataStaEff_pos.getEff(mu2->eta, mu2->pt) * hSysStaSigMCPos->GetBinContent(hSysStaSigMCPos->GetXaxis()->FindBin(mu2->eta), hSysStaSigMCPos->GetYaxis()->FindBin(mu2->pt));
            effdataBkg *= dataStaEff_pos.getEff(mu2->eta, mu2->pt) * hSysStaBkgPos->GetBinContent(hSysStaBkgPos->GetXaxis()->FindBin(mu2->eta), hSysStaBkgPos->GetYaxis()->FindBin(mu2->pt));
            effdata *= dataStaEff_pos.getEff(mu2->eta, mu2->pt);//orig
            effmc   *= zmmStaEff_pos.getEff(mu2->eta, mu2->pt);
          } else {
            effdataFSR *= dataStaEff_neg.getEff(mu2->eta, mu2->pt) * hSysStaSigFSRNeg->GetBinContent(hSysStaSigFSRNeg->GetXaxis()->FindBin(mu2->eta), hSysStaSigFSRNeg->GetYaxis()->FindBin(mu2->pt));
            effdataMC *= dataStaEff_neg.getEff(mu2->eta, mu2->pt) * hSysStaSigMCNeg->GetBinContent(hSysStaSigMCNeg->GetXaxis()->FindBin(mu2->eta), hSysStaSigMCNeg->GetYaxis()->FindBin(mu2->pt));
            effdataBkg *= dataStaEff_neg.getEff(mu2->eta, mu2->pt) * hSysStaBkgNeg->GetBinContent(hSysStaBkgNeg->GetXaxis()->FindBin(mu2->eta), hSysStaBkgNeg->GetYaxis()->FindBin(mu2->pt));
            effdata *= dataStaEff_neg.getEff(mu2->eta, mu2->pt);//orig
            effmc   *= zmmStaEff_neg.getEff(mu2->eta, mu2->pt);
          }
          corrFSR *= effdataFSR/effmc;
          corrMC  *= effdataMC/effmc;
          corrBkg *= effdataBkg/effmc;
          corr *= effdata/effmc;
    
    // std::cout << "corr/corrFSR2  " << (corr/corrFSR-1)*100 << std::endl;
	  // scale factor uncertainties                                                                                                                                         
	  // STANDALONE
          if(mu1->q>0) {
            Double_t effdata = dataStaEff_pos.getEff(mu1->eta, mu1->pt);
            Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(mu1->eta, mu1->pt), dataStaEff_pos.getErrHigh(mu1->eta, mu1->pt));
            Double_t effmc   = zmmStaEff_pos.getEff(mu1->eta, mu1->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(mu1->eta, mu1->pt), zmmStaEff_pos.getErrHigh(mu1->eta, mu1->pt));
            isnan(errdata)?errdata=errmc:errdata=errdata; // stupid b/c fits for standalone are stupid
            Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            // std::cout << " pos effdata " << effdata << "  effmc " << effmc << "  errdata " << errdata << "  errmc " << errmc <<   std::endl;
            hStaErr_pos->Fill(mu1->eta, mu1->pt, errSta);
          } else {
            Double_t effdata = dataStaEff_neg.getEff(mu1->eta, mu1->pt);
            Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(mu1->eta, mu1->pt), dataStaEff_neg.getErrHigh(mu1->eta, mu1->pt));
            Double_t effmc   = zmmStaEff_neg.getEff(mu1->eta, mu1->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(mu1->eta, mu1->pt), zmmStaEff_neg.getErrHigh(mu1->eta, mu1->pt));
            isnan(errdata)?errdata=errmc:errdata=errdata; // stupid b/c fits for standalone are stupid
            Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            // std::cout << " neg effdata " << effdata << "  effmc " << effmc << std::endl;
            hStaErr_neg->Fill(mu1->eta, mu1->pt, errSta);
          }

          if(mu2->q>0) {
            Double_t effdata = dataStaEff_pos.getEff(mu2->eta, mu2->pt);
            Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(mu2->eta, mu2->pt), dataStaEff_pos.getErrHigh(mu2->eta, mu2->pt));
            Double_t effmc   = zmmStaEff_pos.getEff(mu2->eta, mu2->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(mu2->eta, mu2->pt), zmmStaEff_pos.getErrHigh(mu2->eta, mu2->pt));
            isnan(errdata)?errdata=errmc:errdata=errdata;
            Double_t errSta = ((effdata/effmc))*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            // std::cout << " pos effdata " << effdata << "  effmc " << effmc << "  errdata " << errdata << "  errmc " << errmc <<   std::endl;
            hStaErr_pos->Fill(mu2->eta, mu2->pt, errSta);
          } else {
            Double_t effdata = dataStaEff_neg.getEff(mu2->eta, mu2->pt);
            Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(mu2->eta, mu2->pt), dataStaEff_neg.getErrHigh(mu2->eta, mu2->pt));
            Double_t effmc   = zmmStaEff_neg.getEff(mu2->eta, mu2->pt);
            Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(mu2->eta, mu2->pt), zmmStaEff_neg.getErrHigh(mu2->eta, mu2->pt));
            isnan(errdata)?errdata=errmc:errdata=errdata; // stupid b/c fits for standalone are stupid
            Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            // std::cout << " neg effdata " << effdata << "  effmc " << effmc << std::endl;
            hStaErr_neg->Fill(mu2->eta, mu2->pt, errSta);
	  }

      
	  // SELECTION
          if(mu1->q>0) {
            Double_t effdata = dataSelEff_pos.getEff(mu1->eta, mu1->pt);
            Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(mu1->eta, mu1->pt), dataSelEff_pos.getErrHigh(mu1->eta, mu1->pt));
            Double_t effmc   = zmmSelEff_pos.getEff(mu1->eta, mu1->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(mu1->eta, mu1->pt), zmmSelEff_pos.getErrHigh(mu1->eta, mu1->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_pos->Fill(mu1->eta, mu1->pt, errSel);
          } else {
            Double_t effdata = dataSelEff_neg.getEff(mu1->eta, mu1->pt);
            Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(mu1->eta, mu1->pt), dataSelEff_neg.getErrHigh(mu1->eta, mu1->pt));
            Double_t effmc   = zmmSelEff_neg.getEff(mu1->eta, mu1->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(mu1->eta, mu1->pt), zmmSelEff_neg.getErrHigh(mu1->eta, mu1->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_neg->Fill(mu1->eta, mu1->pt, errSel);
          }

          if(mu2->q>0) {
            Double_t effdata = dataSelEff_pos.getEff(mu2->eta, mu2->pt);
            Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(mu2->eta, mu2->pt), dataSelEff_pos.getErrHigh(mu2->eta, mu2->pt));
            Double_t effmc   = zmmSelEff_pos.getEff(mu2->eta, mu2->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(mu2->eta, mu2->pt), zmmSelEff_pos.getErrHigh(mu2->eta, mu2->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_pos->Fill(mu2->eta, mu2->pt, errSel);
          } else {
            Double_t effdata = dataSelEff_neg.getEff(mu2->eta, mu2->pt);
            Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(mu2->eta, mu2->pt), dataSelEff_neg.getErrHigh(mu2->eta, mu2->pt));
            Double_t effmc   = zmmSelEff_neg.getEff(mu2->eta, mu2->pt);
            Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(mu2->eta, mu2->pt), zmmSelEff_neg.getErrHigh(mu2->eta, mu2->pt));
            Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hSelErr_neg->Fill(mu2->eta, mu2->pt, errSel);
	  }

      
	  //HLT
          if(mu1->q>0) {
            Double_t effdata = dataHLTEff_pos.getEff(mu1->eta, mu1->pt);
            Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(mu1->eta, mu1->pt), dataHLTEff_pos.getErrHigh(mu1->eta, mu1->pt));
            Double_t effmc   = zmmHLTEff_pos.getEff(mu1->eta, mu1->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(mu1->eta, mu1->pt), zmmHLTEff_pos.getErrHigh(mu1->eta, mu1->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_pos->Fill(mu1->eta, mu1->pt, errHLT);
          } else {
            Double_t effdata = dataHLTEff_neg.getEff(mu1->eta, mu1->pt);
            Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(mu1->eta, mu1->pt), dataHLTEff_neg.getErrHigh(mu1->eta, mu1->pt));
            Double_t effmc   = zmmHLTEff_neg.getEff(mu1->eta, mu1->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(mu1->eta, mu1->pt), zmmHLTEff_neg.getErrHigh(mu1->eta, mu1->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_neg->Fill(mu1->eta, mu1->pt, errHLT);
          }

          if(mu2->q>0) {
            Double_t effdata = dataHLTEff_pos.getEff(mu2->eta, mu2->pt);
            Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(mu2->eta, mu2->pt), dataHLTEff_pos.getErrHigh(mu2->eta, mu2->pt));
            Double_t effmc   = zmmHLTEff_pos.getEff(mu2->eta, mu2->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(mu2->eta, mu2->pt), zmmHLTEff_pos.getErrHigh(mu2->eta, mu2->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_pos->Fill(mu2->eta, mu2->pt, errHLT);
          } else {
            Double_t effdata = dataHLTEff_neg.getEff(mu2->eta, mu2->pt);
            Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(mu2->eta, mu2->pt), dataHLTEff_neg.getErrHigh(mu2->eta, mu2->pt));
            Double_t effmc   = zmmHLTEff_neg.getEff(mu2->eta, mu2->pt);
            Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(mu2->eta, mu2->pt), zmmHLTEff_neg.getErrHigh(mu2->eta, mu2->pt));
            Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
            hHLTErr_neg->Fill(mu2->eta, mu2->pt, errHLT);
          }
          if(alreadyDid) continue;
          alreadyDid=true;
	  // std::cout << info->evtNum << " " << corr << " " << std::endl;
	  nSelv[ifile]    +=weight;
	  nSelCorrvFSR[ifile]+=weight*corrFSR;
	  nSelCorrvMC[ifile]+=weight*corrMC;
	  nSelCorrvBkg[ifile]+=weight*corrBkg;
	  nSelCorrv[ifile]+=weight*corr;
      // std::cout << "corr " << corr << " corr FSR " << corrFSR << "  corr MC " << corrMC << "  corr Bkg " << corrBkg << std::endl;
	  nSelCorrVarvFSR[ifile]+=weight*weight*corrFSR*corrFSR;
	  nSelCorrVarvMC[ifile]+=weight*weight*corrMC*corrMC;
	  nSelCorrVarvBkg[ifile]+=weight*weight*corrBkg*corrBkg;
	  nSelCorrVarv[ifile]+=weight*weight*corr*corr;
        }
      }      
    }

    // std::cout << "nSelCorrVarv[ifile]  " <<  nSelCorrVarv[ifile] << std::endl;
    // std::cout << "var" << std::endl;
    Double_t var=0;
    for(Int_t iy=0; iy<=hHLTErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_pos->GetNbinsX(); ix++) {
        Double_t err=hHLTErr_pos->GetBinContent(ix,iy);
        var+=err*err;
        // std::cout << "hlt pos " << var << std::endl;
      }
    }
    for(Int_t iy=0; iy<=hHLTErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_neg->GetNbinsX(); ix++) {
        Double_t err=hHLTErr_neg->GetBinContent(ix,iy);
        var+=err*err;
        // std::cout << "hlt neg " << var << std::endl;
      }
    }
    for(Int_t iy=0; iy<=hSelErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_pos->GetNbinsX(); ix++) {
        Double_t err=hSelErr_pos->GetBinContent(ix,iy);
        var+=err*err;
        // std::cout << "sel pos " << var << std::endl;
      }
    }
    for(Int_t iy=0; iy<=hSelErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_neg->GetNbinsX(); ix++) {
        Double_t err=hSelErr_neg->GetBinContent(ix,iy);
        var+=err*err;
        // std::cout << "sel neg " << var << std::endl;
      }
    }
    for(Int_t iy=0; iy<=hStaErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_pos->GetNbinsX(); ix++) {
        Double_t err=hStaErr_pos->GetBinContent(ix,iy);
	    var+=err*err;
        // std::cout << "sta pos " << var << std::endl;
      }
    }
    for(Int_t iy=0; iy<=hStaErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_neg->GetNbinsX(); ix++) {
        Double_t err=hStaErr_neg->GetBinContent(ix,iy);
        // err=0;// for now
	    var+=err*err;
        // std::cout << "sta neg " << var << std::endl;
      }
    }
    // std::cout << "blah  " << std::endl;
    nSelCorrVarvFSR[ifile]+=var;
    nSelCorrVarvMC[ifile]+=var;
    nSelCorrVarvBkg[ifile]+=var;
    nSelCorrVarv[ifile]+=var;
    // std::cout << "comput acceptances" << std::endl;
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(accv[ifile]*sqrt((1.+accv[ifile])/nEvtsv[ifile]));
    accCorrvFSR.push_back(nSelCorrvFSR[ifile]/nEvtsv[ifile]); 
    accErrCorrvFSR.push_back(accCorrvFSR[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvFSR[ifile]*nSelCorrvFSR[ifile]) + 1./nEvtsv[ifile]));
    // std::cout << "acccorrvfsr " << accCorrvFSR[ifile] << "  nSelCorrVarv[ifile] " << nSelCorrVarv[ifile] << " nSelCorrvFSR[ifile] " << nSelCorrvFSR[ifile] << "  nEvtsv[ifile] " << nEvtsv[ifile] << std::endl;
    accCorrvMC.push_back(nSelCorrvMC[ifile]/nEvtsv[ifile]);   accErrCorrvMC.push_back(accCorrvMC[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvMC[ifile]*nSelCorrvMC[ifile]) + 1./nEvtsv[ifile]));
    accCorrvBkg.push_back(nSelCorrvBkg[ifile]/nEvtsv[ifile]); accErrCorrvBkg.push_back(accCorrvBkg[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvBkg[ifile]*nSelCorrvBkg[ifile]) + 1./nEvtsv[ifile]));
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); accErrCorrv.push_back(accCorrv[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  } 
  
    // std::cout << "clean up memeory sturff" << std::endl;
  delete info;
  delete gen;
  delete muonArr;

    
    std::cout << "print output" << std::endl;
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================    
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    cout << "          FSR unc: " << accCorrvFSR[ifile] << " +/- " << accErrCorrvFSR[ifile] << endl;
    cout << "           MC unc: " << accCorrvMC[ifile]  << " +/- " << accErrCorrvMC[ifile]  << endl;
    cout << "          Bkg unc: " << accCorrvBkg[ifile] << " +/- " << accErrCorrvBkg[ifile] << endl;
    cout << "  fraction passing gen cut: " << nEvtsv[ifile] << " / " << nAllv[ifile] << " = " << nEvtsv[ifile]/nAllv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/binned.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> mu mu" << endl;
  txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    txtfile << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    txtfile << "          FSR unc: " << accCorrvFSR[ifile] << " +/- " << accErrCorrvFSR[ifile] << endl;
    txtfile << "           MC unc: " << accCorrvMC[ifile]  << " +/- " << accErrCorrvMC[ifile]  << endl;
    txtfile << "          Bkg unc: " << accCorrvBkg[ifile] << " +/- " << accErrCorrvBkg[ifile] << endl;
    txtfile << "  fraction passing gen cut: " << nEvtsv[ifile] << " / " << nAllv[ifile] << " = " << nEvtsv[ifile]/nAllv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZmmBinned"); 
}
