//================================================================================================
//
// Compute W->munu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties need to be checked//
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
#include "TLorentzVector.h"

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
#include "CEffUser2D.hh"
#endif


//=== MAIN MACRO ================================================================================================= 

void computeAccSelWm(const TString conf,       // input file
                     const TString outputDir,  // output directory
		     const Int_t   charge,      // 0 = inclusive, +1 = W+, -1 = W-
		     const Int_t   doPU
) {
  gBenchmark->Start("computeAccSelWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;

  const Double_t VETO_PT   = 10;
  const Double_t VETO_ETA  = 2.4;

  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 13;
  
  // efficiency files
  TString dataHLTEffName("/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuHLTEff/MG/eff.root");
  TString zmmHLTEffName( "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuHLTEff/CT/eff.root");
  TString dataSelEffName("/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuSITEff/MG/eff.root");
  TString zmmSelEffName( "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuSITEff/CT/eff.root");
  TString dataTrkEffName("/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/MG/eff.root");
  TString zmmTrkEffName( "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/CT/eff.root");
  TString dataStaEffName("/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/MG/eff.root");
  TString zmmStaEffName( "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/CT/eff.root");
  if(charge==1) {
    dataHLTEffName ="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuHLTEff/MG/eff.root";
    zmmHLTEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuHLTEff/CT/eff.root";
    dataSelEffName ="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuSITEff/MG/eff.root";
    zmmSelEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuSITEff/CT/eff.root";
    dataTrkEffName ="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/MG/eff.root";
    zmmTrkEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/CT/eff.root";
    dataStaEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/MG/eff.root";
    zmmStaEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/CT/eff.root";
  }
  if(charge==-1) {
    dataHLTEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuHLTEff/MG/eff.root";
    zmmHLTEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuHLTEff/CT/eff.root";
    dataSelEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuSITEff/MG/eff.root";
    zmmSelEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuSITEff/CT/eff.root";
    dataTrkEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/MG/eff.root";
    zmmTrkEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/CT/eff.root";
    dataStaEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/MG/eff.root";
    zmmStaEffName = "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/MuStaEff/CT/eff.root";
  }

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_rw_76X.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("h_rw_golden");


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
  
  //
  // Get efficiency
  //
  TFile *dataHLTEffFile = new TFile(dataHLTEffName);
  CEffUser2D dataHLTEff;
  TH2D *hHLTErr=0, *hHLTErrB=0, *hHLTErrE=0;
  if(dataHLTEffFile) {    
    dataHLTEff.loadEff((TH2D*)dataHLTEffFile->Get("hEffEtaPt"), 
                       (TH2D*)dataHLTEffFile->Get("hErrlEtaPt"),
		       (TH2D*)dataHLTEffFile->Get("hErrhEtaPt"));
    
    TH2D* h =(TH2D*)dataHLTEffFile->Get("hEffEtaPt");
    hHLTErr  = new TH2D("hHLTErr", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hHLTErrB = new TH2D("hHLTErrB","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hHLTErrE = new TH2D("hHLTErrE","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  
  TFile *zmmHLTEffFile = new TFile(zmmHLTEffName);
  CEffUser2D zmmHLTEff;
  if(zmmHLTEffFile) {
    zmmHLTEff.loadEff((TH2D*)zmmHLTEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmHLTEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmHLTEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataSelEffFile = new TFile(dataSelEffName);
  CEffUser2D dataSelEff;
  TH2D *hSelErr=0, *hSelErrB=0, *hSelErrE=0;
  if(dataSelEffFile) {
    dataSelEff.loadEff((TH2D*)dataSelEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataSelEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataSelEffFile->Get("hErrhEtaPt"));
    
    TH2D* h =(TH2D*)dataSelEffFile->Get("hEffEtaPt");
    hSelErr  = new TH2D("hSelErr", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hSelErrB = new TH2D("hSelErrB","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hSelErrE = new TH2D("hSelErrE","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  
  TFile *zmmSelEffFile = new TFile(zmmSelEffName);
  CEffUser2D zmmSelEff;
  if(zmmSelEffFile) {
    zmmSelEff.loadEff((TH2D*)zmmSelEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmSelEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmSelEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataTrkEffFile = new TFile(dataTrkEffName);
  CEffUser2D dataTrkEff;
  TH2D *hTrkErr=0, *hTrkErrB=0, *hTrkErrE=0;
  if(dataTrkEffFile) {
    dataTrkEff.loadEff((TH2D*)dataTrkEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataTrkEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataTrkEffFile->Get("hErrhEtaPt"));
    
    TH2D* h =(TH2D*)dataTrkEffFile->Get("hEffEtaPt");
    hTrkErr  = new TH2D("hTrkErr", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hTrkErrB = new TH2D("hTrkErrB","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hTrkErrE = new TH2D("hTrkErrE","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  
  TFile *zmmTrkEffFile = new TFile(zmmTrkEffName);
  CEffUser2D zmmTrkEff;
  if(zmmTrkEffFile) {
    zmmTrkEff.loadEff((TH2D*)zmmTrkEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmTrkEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmTrkEffFile->Get("hErrhEtaPt"));
  }
  
  TFile *dataStaEffFile = new TFile(dataStaEffName);
  CEffUser2D dataStaEff;
  TH2D *hStaErr=0, *hStaErrB=0, *hStaErrE=0;
  if(dataStaEffFile) {
    dataStaEff.loadEff((TH2D*)dataStaEffFile->Get("hEffEtaPt"),
                       (TH2D*)dataStaEffFile->Get("hErrlEtaPt"),
                       (TH2D*)dataStaEffFile->Get("hErrhEtaPt"));
    
    TH2D* h =(TH2D*)dataStaEffFile->Get("hEffEtaPt");
    hStaErr  = new TH2D("hStaErr", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hStaErrB = new TH2D("hStaErrB","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
    hStaErrE = new TH2D("hStaErrE","",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                                      h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  
  TFile *zmmStaEffFile = new TFile(zmmStaEffName);
  CEffUser2D zmmStaEff;
  if(zmmStaEffFile) {
    zmmStaEff.loadEff((TH2D*)zmmStaEffFile->Get("hEffEtaPt"),
                      (TH2D*)zmmStaEffFile->Get("hErrlEtaPt"),
                      (TH2D*)zmmStaEffFile->Get("hErrhEtaPt"));
  }
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo   *info = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray     *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray        *muonArr = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;

  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBv, nSelEv;
  vector<Double_t> accv, accBv, accEv;
  vector<Double_t> accErrv, accErrBv, accErrEv;
  vector<Double_t> nSelCorrv, nSelBCorrv, nSelECorrv;
  vector<Double_t> nSelCorrVarv, nSelBCorrVarv, nSelECorrVarv;
  vector<Double_t> accCorrv, accBCorrv, accECorrv;
  vector<Double_t> accErrCorrv, accErrBCorrv, accErrECorrv;

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
    eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr    = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr     = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch *genPartBr = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Muon",          &muonArr); TBranch *muonBr    = eventTree->GetBranch("Muon"); 
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBv.push_back(0);
    nSelEv.push_back(0);
    nSelCorrv.push_back(0);
    nSelBCorrv.push_back(0);
    nSelECorrv.push_back(0);
    nSelCorrVarv.push_back(0);
    nSelBCorrVarv.push_back(0);
    nSelECorrVarv.push_back(0);
    
    for(Int_t iy=0; iy<=hHLTErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr->GetNbinsX(); ix++) {
        hHLTErr ->SetBinContent(ix,iy,0);
        hHLTErrB->SetBinContent(ix,iy,0);
        hHLTErrE->SetBinContent(ix,iy,0);
      }
    }
    for(Int_t iy=0; iy<=hSelErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr->GetNbinsX(); ix++) {
        hSelErr ->SetBinContent(ix,iy,0);
        hSelErrB->SetBinContent(ix,iy,0);
        hSelErrE->SetBinContent(ix,iy,0);
      }
    }
    for(Int_t iy=0; iy<=hTrkErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hTrkErr->GetNbinsX(); ix++) {
        hTrkErr ->SetBinContent(ix,iy,0);
        hTrkErrB->SetBinContent(ix,iy,0);
        hTrkErrE->SetBinContent(ix,iy,0);
      }
    }
    for(Int_t iy=0; iy<=hStaErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr->GetNbinsX(); ix++) {
        hStaErr ->SetBinContent(ix,iy,0);
        hStaErrB->SetBinContent(ix,iy,0);
        hStaErrE->SetBinContent(ix,iy,0);
      }
    }    
    
    //
    // loop over events
    //    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    //for(UInt_t ientry=0; ientry<1000000; ientry++) {
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);

      if (charge==-1 && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
      if (charge==1 && toolbox::flavor(genPartArr, BOSON_ID)!=-LEPTON_ID) continue;
      if (charge==0 && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      /*TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,1);*/
      
      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;
      
      // trigger requirement               
      if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
   
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);
      Int_t nLooseLep=0;
      const baconhep::TMuon *goodMuon=0;
      Bool_t passSel=kFALSE;
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
  	const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);

        if(fabs(mu->eta) > VETO_ETA) continue; // loose lepton |eta| cut
        if(mu->pt	 < VETO_PT)  continue; // loose lepton pT cut
        if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
        if(nLooseLep>1) {  // extra lepton veto
          passSel=kFALSE;
          break;
        }
        
        if(fabs(mu->eta) > ETA_CUT)         continue;  // lepton |eta| cut
        if(mu->pt < PT_CUT)		    continue;  // lepton pT cut	
        if(!passMuonID(mu))		    continue;  // lepton selection
	if(!isMuonTriggerObj(triggerMenu, mu->hltMatchBits, kFALSE)) continue;
	
	if(charge!=0 && mu->q!=charge) continue;  // check charge (if necessary)
	
	passSel=kTRUE;
	goodMuon=mu;
      }
      
      if(passSel) {
        
	/******** We have a W candidate! HURRAY! ********/
        
	Bool_t isBarrel = (fabs(goodMuon->eta)<ETA_BARREL) ? kTRUE : kFALSE;
        
	// data/MC scale factor corrections
	Double_t corr=1;
	if(dataHLTEffFile && zmmHLTEffFile) {
	  Double_t effdata = dataHLTEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmHLTEff.getEff(goodMuon->eta, goodMuon->pt);
	  corr *= effdata/effmc;
	}
	if(dataSelEffFile && zmmSelEffFile) {
	  Double_t effdata = dataSelEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmSelEff.getEff(goodMuon->eta, goodMuon->pt);
	  corr *= effdata/effmc;
	}
	if(dataTrkEffFile && zmmTrkEffFile) {
	  Double_t effdata = dataTrkEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmTrkEff.getEff(goodMuon->eta, goodMuon->pt);
	  //corr *= effdata/effmc;
	}
	if(dataStaEffFile && zmmStaEffFile) {
	  Double_t effdata = dataStaEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmStaEff.getEff(goodMuon->eta, goodMuon->pt);
	  corr *= effdata/effmc;
	}
	
	// scale factor uncertainties
	if(dataHLTEffFile && zmmHLTEffFile) {
	  Double_t effdata = dataHLTEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmHLTEff.getEff(goodMuon->eta, goodMuon->pt);	  
	  Double_t errdata = TMath::Max(dataHLTEff.getErrLow(goodMuon->eta, goodMuon->pt),dataHLTEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t errmc   = TMath::Max(zmmHLTEff.getErrLow(goodMuon->eta, goodMuon->pt), zmmHLTEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t err     = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);
	  hHLTErr->Fill(goodMuon->eta,goodMuon->pt,err);
	  if(isBarrel) hHLTErrB->Fill(goodMuon->eta,goodMuon->pt,err);
	  else         hHLTErrE->Fill(goodMuon->eta,goodMuon->pt,err);
	}
	if(dataSelEffFile && zmmSelEffFile) {
	  Double_t effdata = dataSelEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmSelEff.getEff(goodMuon->eta, goodMuon->pt);	  
	  Double_t errdata = TMath::Max(dataSelEff.getErrLow(goodMuon->eta, goodMuon->pt),dataSelEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t errmc   = TMath::Max(zmmSelEff.getErrLow(goodMuon->eta, goodMuon->pt), zmmSelEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t err     = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);
	  hSelErr->Fill(goodMuon->eta,goodMuon->pt,err);
	  if(isBarrel) hSelErrB->Fill(goodMuon->eta,goodMuon->pt,err);
	  else         hSelErrE->Fill(goodMuon->eta,goodMuon->pt,err);
	}
	if(dataTrkEffFile && zmmTrkEffFile) {
	  Double_t effdata = dataTrkEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmTrkEff.getEff(goodMuon->eta, goodMuon->pt);	  
	  Double_t errdata = TMath::Max(dataTrkEff.getErrLow(goodMuon->eta, goodMuon->pt),dataTrkEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t errmc   = TMath::Max(zmmTrkEff.getErrLow(goodMuon->eta, goodMuon->pt), zmmTrkEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t err     = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);
	  /* if(goodMuon->eta>1.2 && goodMuon->eta<2.1) 
	    {
	      err=0.0013;
	      }*/
	  hTrkErr->Fill(goodMuon->eta,goodMuon->pt,err);
	  if(isBarrel) hTrkErrB->Fill(goodMuon->eta,goodMuon->pt,err);
	  else         hTrkErrE->Fill(goodMuon->eta,goodMuon->pt,err);
	}
	if(dataStaEffFile && zmmStaEffFile) {
	  Double_t effdata = dataStaEff.getEff(goodMuon->eta, goodMuon->pt);
	  Double_t effmc   = zmmStaEff.getEff(goodMuon->eta, goodMuon->pt);	  
	  Double_t errdata = TMath::Max(dataStaEff.getErrLow(goodMuon->eta, goodMuon->pt),dataStaEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t errmc   = TMath::Max(zmmStaEff.getErrLow(goodMuon->eta, goodMuon->pt), zmmStaEff.getErrHigh(goodMuon->eta, goodMuon->pt));
	  Double_t err     = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata+errmc*errmc/effmc/effmc);
	  hStaErr->Fill(goodMuon->eta,goodMuon->pt,err);
	  if(isBarrel) hStaErrB->Fill(goodMuon->eta,goodMuon->pt,err);
	  else         hStaErrE->Fill(goodMuon->eta,goodMuon->pt,err);
	}
	
	nSelv[ifile]+=weight;
	nSelCorrv[ifile]+=weight*corr;
	nSelCorrVarv[ifile]+=weight*weight*corr*corr;
  	if(isBarrel) { 
	  nSelBv[ifile]+=weight;
	  nSelBCorrv[ifile]+=weight*corr;
	  nSelBCorrVarv[ifile]+=weight*weight*corr*corr;
	  	
	} else { 
	  nSelEv[ifile]+=weight;
	  nSelECorrv[ifile]+=weight*corr;
	  nSelECorrVarv[ifile]+=weight*weight*corr*corr;
	}
      }
    }
    
    Double_t var=0, varB=0, varE=0;
    for(Int_t iy=0; iy<=hHLTErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr->GetNbinsX(); ix++) {
        Double_t err;
	err=hHLTErr->GetBinContent(ix,iy);  var+=err*err;
        err=hHLTErrB->GetBinContent(ix,iy); varB+=err*err;
        err=hHLTErrE->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hSelErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr->GetNbinsX(); ix++) {
        Double_t err;
	err=hSelErr->GetBinContent(ix,iy);  var+=err*err;
	err=hSelErrB->GetBinContent(ix,iy); varB+=err*err;
	err=hSelErrE->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hTrkErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hTrkErr->GetNbinsX(); ix++) {
        Double_t err;
	err=hTrkErr->GetBinContent(ix,iy);  var+=0;//err*err;
        err=hTrkErrB->GetBinContent(ix,iy); varB+=0;//err*err;
        err=hTrkErrE->GetBinContent(ix,iy); varE+=0;//err*err;
      }
    }
    for(Int_t iy=0; iy<=hStaErr->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr->GetNbinsX(); ix++) {
        Double_t err;
	err=hStaErr->GetBinContent(ix,iy);  var+=err*err;
	err=hStaErrB->GetBinContent(ix,iy); varB+=err*err;
	err=hStaErrE->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    nSelCorrVarv[ifile]+=var;
    nSelBCorrVarv[ifile]+=varB;
    nSelECorrVarv[ifile]+=varE;
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.+accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.+accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.+accEv[ifile])/nEvtsv[ifile]));
    
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]);   accErrCorrv.push_back(accCorrv[ifile]*sqrt(nSelCorrVarv[ifile]/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));
    accBCorrv.push_back(nSelBCorrv[ifile]/nEvtsv[ifile]); accErrBCorrv.push_back(accBCorrv[ifile]*sqrt(nSelBCorrVarv[ifile]/nSelBCorrv[ifile]/nSelBCorrv[ifile] + 1./nEvtsv[ifile]));
    accECorrv.push_back(nSelECorrv[ifile]/nEvtsv[ifile]); accErrECorrv.push_back(accECorrv[ifile]*sqrt(nSelECorrVarv[ifile]/(nSelECorrv[ifile]*nSelECorrv[ifile]) + 1./nEvtsv[ifile]));
   
    delete infile;
    infile=0, eventTree=0;  
  }
  delete info;
  delete gen;
  delete muonArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  if(charge== 0) cout << " W -> mu nu"  << endl;
  if(charge==-1) cout << " W- -> mu nu" << endl;
  if(charge== 1) cout << " W+ -> mu nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    cout << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    cout << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    cout << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    cout << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    cout << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/sel.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  if(charge== 0) txtfile << " W -> mu nu"  << endl;
  if(charge==-1) txtfile << " W- -> mu nu" << endl;
  if(charge== 1) txtfile << " W+ -> mu nu" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    txtfile << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    txtfile << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    txtfile << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    txtfile << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    txtfile << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelWm"); 
}
