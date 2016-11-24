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

void computeAccSelWm_Inclusive(const TString conf,       // input file
		     const TString inputDir,
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
  const TString dataHLTEffName_pos = inputDir + "MuHLTEff/1MGpositive/eff.root";
  const TString dataHLTEffName_neg = inputDir + "MuHLTEff/1MGnegative/eff.root";
  const TString zmmHLTEffName_pos  = inputDir + "MuHLTEff/1CTpositive/eff.root";
  const TString zmmHLTEffName_neg  = inputDir + "MuHLTEff/1CTnegative/eff.root";

  const TString dataSelEffName_pos = inputDir + "MuSITEff/1MGpositive/eff.root";
  const TString dataSelEffName_neg = inputDir + "MuSITEff/1MGnegative/eff.root";
  const TString zmmSelEffName_pos  = inputDir + "MuSITEff/1CTpositive/eff.root";
  const TString zmmSelEffName_neg  = inputDir + "MuSITEff/1CTnegative/eff.root";

  const TString dataTrkEffName_pos = inputDir + "MuSITEff/1MGpositive/eff.root";
  const TString dataTrkEffName_neg = inputDir + "MuSITEff/1MGnegative/eff.root";
  const TString zmmTrkEffName_pos  = inputDir + "MuSITEff/1CTpositive/eff.root";
  const TString zmmTrkEffName_neg  = inputDir + "MuSITEff/1CTnegative/eff.root";

  const TString dataStaEffName_pos = inputDir + "MuStaEff/1MGpositive/eff.root";
  const TString dataStaEffName_neg = inputDir + "MuStaEff/1MGnegative/eff.root";
  const TString zmmStaEffName_pos  = inputDir + "MuStaEff/1CTpositive/eff.root";
  const TString zmmStaEffName_neg  = inputDir + "MuStaEff/1CTnegative/eff.root";

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
  
  //
  // Get efficiency
  //
  TH2D *h=0;

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
  TH2D *hHLTErrB_pos = new TH2D("hHLTErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrB_neg = new TH2D("hHLTErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrE_pos = new TH2D("hHLTErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hHLTErrE_neg = new TH2D("hHLTErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());


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
  TH2D *hSelErrB_pos = new TH2D("hSelErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrB_neg = new TH2D("hSelErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrE_pos = new TH2D("hSelErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hSelErrE_neg = new TH2D("hSelErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());


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
  TH2D *hStaErrB_pos = new TH2D("hStaErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrB_neg = new TH2D("hStaErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrE_pos = new TH2D("hStaErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hStaErrE_neg = new TH2D("hStaErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());


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

  h =(TH2D*)dataTrkEffFile_pos->Get("hEffEtaPt");
  TH2D *hTrkErr_pos = new TH2D("hTrkErr_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hTrkErr_neg = new TH2D("hTrkErr_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hTrkErrB_pos = new TH2D("hTrkErrB_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hTrkErrB_neg = new TH2D("hTrkErrB_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hTrkErrE_pos = new TH2D("hTrkErrE_pos", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  TH2D *hTrkErrE_neg = new TH2D("hTrkErrE_neg", "",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
                               h->GetNbinsY(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

 
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
         Double_t effdata, effmc;
         Double_t corr=1;

         effdata=1; effmc=1;
         if(goodMuon->q>0) {
           effdata *= dataHLTEff_pos.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmHLTEff_pos.getEff(goodMuon->eta, goodMuon->pt);
         } else {
           effdata *= dataHLTEff_neg.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmHLTEff_neg.getEff(goodMuon->eta, goodMuon->pt);
         }
	 corr *= effdata/effmc;

         effdata=1; effmc=1;
         if(goodMuon->q>0) {
           effdata *= dataSelEff_pos.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmSelEff_pos.getEff(goodMuon->eta, goodMuon->pt);
         } else {
           effdata *= dataSelEff_neg.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmSelEff_neg.getEff(goodMuon->eta, goodMuon->pt);
         }
         corr *= effdata/effmc;

         effdata=1; effmc=1;
         if(goodMuon->q>0) {
           effdata *= dataStaEff_pos.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmStaEff_pos.getEff(goodMuon->eta, goodMuon->pt);
         } else {
           effdata *= dataStaEff_neg.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmStaEff_neg.getEff(goodMuon->eta, goodMuon->pt);
         }
         corr *= effdata/effmc;

         effdata=1; effmc=1;
         if(goodMuon->q>0) {
           effdata *= dataTrkEff_pos.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmTrkEff_pos.getEff(goodMuon->eta, goodMuon->pt);
         } else {
           effdata *= dataTrkEff_neg.getEff(goodMuon->eta, goodMuon->pt);
           effmc   *= zmmTrkEff_neg.getEff(goodMuon->eta, goodMuon->pt);
         }
         //corr *= effdata/effmc;

	// scale factor uncertainties
	// TRACKER
        if(goodMuon->q>0) {
          Double_t effdata = dataTrkEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataTrkEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), dataTrkEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmTrkEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmTrkEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), zmmTrkEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hTrkErr_pos->Fill(goodMuon->eta, goodMuon->pt, errTrk);
          if(isBarrel) hTrkErrB_pos->Fill(goodMuon->eta,goodMuon->pt,errTrk);
          else         hTrkErrE_pos->Fill(goodMuon->eta,goodMuon->pt,errTrk);
        } else {
          Double_t effdata = dataTrkEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataTrkEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), dataTrkEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmTrkEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmTrkEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), zmmTrkEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errTrk = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hTrkErr_neg->Fill(goodMuon->eta, goodMuon->pt, errTrk);
          if(isBarrel) hTrkErrB_neg->Fill(goodMuon->eta,goodMuon->pt,errTrk);
          else         hTrkErrE_neg->Fill(goodMuon->eta,goodMuon->pt,errTrk);
        }

	// STANDALONE
        if(goodMuon->q>0) {
          Double_t effdata = dataStaEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), dataStaEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmStaEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), zmmStaEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hStaErr_pos->Fill(goodMuon->eta, goodMuon->pt, errSta);
	  if(isBarrel) hStaErrB_pos->Fill(goodMuon->eta,goodMuon->pt,errSta);
	  else	       hStaErrE_pos->Fill(goodMuon->eta,goodMuon->pt,errSta);
        } else {
          Double_t effdata = dataStaEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), dataStaEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmStaEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), zmmStaEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hStaErr_neg->Fill(goodMuon->eta, goodMuon->pt, errSta);
          if(isBarrel) hStaErrB_neg->Fill(goodMuon->eta,goodMuon->pt,errSta);
          else         hStaErrE_neg->Fill(goodMuon->eta,goodMuon->pt,errSta);
        }

	// SELECTION
        if(goodMuon->q>0) {
          Double_t effdata = dataSelEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), dataSelEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmSelEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), zmmSelEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hSelErr_pos->Fill(goodMuon->eta, goodMuon->pt, errSel);
          if(isBarrel) hSelErrB_pos->Fill(goodMuon->eta,goodMuon->pt,errSel);
          else         hSelErrE_pos->Fill(goodMuon->eta,goodMuon->pt,errSel);
        } else {
          Double_t effdata = dataSelEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), dataSelEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmSelEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), zmmSelEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hSelErr_neg->Fill(goodMuon->eta, goodMuon->pt, errSel);
          if(isBarrel) hSelErrB_neg->Fill(goodMuon->eta,goodMuon->pt,errSel);
          else         hSelErrE_neg->Fill(goodMuon->eta,goodMuon->pt,errSel);
        }

	//HLT
        if(goodMuon->q>0) {
          Double_t effdata = dataHLTEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), dataHLTEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmHLTEff_pos.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(goodMuon->eta, goodMuon->pt), zmmHLTEff_pos.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hHLTErr_pos->Fill(goodMuon->eta, goodMuon->pt, errHLT);
          if(isBarrel) hHLTErrB_pos->Fill(goodMuon->eta,goodMuon->pt,errHLT);
          else         hHLTErrE_pos->Fill(goodMuon->eta,goodMuon->pt,errHLT);
        } else {
          Double_t effdata = dataHLTEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), dataHLTEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t effmc   = zmmHLTEff_neg.getEff(goodMuon->eta, goodMuon->pt);
          Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(goodMuon->eta, goodMuon->pt), zmmHLTEff_neg.getErrHigh(goodMuon->eta, goodMuon->pt));
          Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
          hHLTErr_neg->Fill(goodMuon->eta, goodMuon->pt, errHLT);
          if(isBarrel) hHLTErrB_neg->Fill(goodMuon->eta,goodMuon->pt,errHLT);
          else         hHLTErrE_neg->Fill(goodMuon->eta,goodMuon->pt,errHLT);
        }

        nSelv[ifile]    +=weight;
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
    for(Int_t iy=0; iy<=hHLTErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_pos->GetNbinsX(); ix++) {
        Double_t err;
	err=hHLTErr_pos->GetBinContent(ix,iy);  var+=err*err;
        err=hHLTErrB_pos->GetBinContent(ix,iy); varB+=err*err;
        err=hHLTErrE_pos->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hHLTErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_neg->GetNbinsX(); ix++) {
        Double_t err;
        err=hHLTErr_neg->GetBinContent(ix,iy);  var+=err*err;
        err=hHLTErrB_neg->GetBinContent(ix,iy); varB+=err*err;
        err=hHLTErrE_neg->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hSelErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_pos->GetNbinsX(); ix++) {
        Double_t err;
        err=hSelErr_pos->GetBinContent(ix,iy);  var+=err*err;
        err=hSelErrB_pos->GetBinContent(ix,iy); varB+=err*err;
        err=hSelErrE_pos->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hSelErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_neg->GetNbinsX(); ix++) {
        Double_t err;
        err=hSelErr_neg->GetBinContent(ix,iy);  var+=err*err;
        err=hSelErrB_neg->GetBinContent(ix,iy); varB+=err*err;
        err=hSelErrE_neg->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hTrkErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hTrkErr_pos->GetNbinsX(); ix++) {
        Double_t err;
        err=hTrkErr_pos->GetBinContent(ix,iy);  var+=0.0;//err*err;
        err=hTrkErrB_pos->GetBinContent(ix,iy); varB+=0.0;//err*err;
        err=hTrkErrE_pos->GetBinContent(ix,iy); varE+=0.0;//err*err;
      }
    }
    for(Int_t iy=0; iy<=hTrkErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hTrkErr_neg->GetNbinsX(); ix++) {
        Double_t err;
        err=hTrkErr_neg->GetBinContent(ix,iy);  var+=0.0;//err*err;
        err=hTrkErrB_neg->GetBinContent(ix,iy); varB+=0.0;//err*err;
        err=hTrkErrE_neg->GetBinContent(ix,iy); varE+=0.0;//err*err;
      }
    }
    for(Int_t iy=0; iy<=hStaErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_pos->GetNbinsX(); ix++) {
        Double_t err;
        err=hStaErr_pos->GetBinContent(ix,iy);  var+=err*err;
        err=hStaErrB_pos->GetBinContent(ix,iy); varB+=err*err;
        err=hStaErrE_pos->GetBinContent(ix,iy); varE+=err*err;
      }
    }
    for(Int_t iy=0; iy<=hStaErr_neg->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_neg->GetNbinsX(); ix++) {
        Double_t err;
        err=hStaErr_neg->GetBinContent(ix,iy);  var+=err*err;
        err=hStaErrB_neg->GetBinContent(ix,iy); varB+=err*err;
        err=hStaErrE_neg->GetBinContent(ix,iy); varE+=err*err;
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
