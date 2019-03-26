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
#include "TLorentzVector.h"         // 4-vector class
#include "TCanvas.h"

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
// #include "CEffUser1D.hh"
// #include "CEffUser2D.hh"
#endif

const float fitMassLo = 60;
const float fitMassHi = 120;
const int massBin = 60;

// const int muEtaNB = 12;
// const float muEtaRange[muEtaNB+1] = {-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4};
// const int muPtNB = 8;
// const float muPtRange[muPtNB+1] = {25,30,35,40,45,50,60,80,8000};

  const int muEtaNB = 2;
  const float muEtaRange[muEtaNB+1] = {-2.4,0,2.4};
  const int muPtNB = 2;
  const float muPtRange[muPtNB+1] = {25,40,8000};

//=== MAIN MACRO ================================================================================================= 

void makeGenMllTemplate_forElectron_v2(
	const TString conf,      // input file
	const TString outputDir, // output directory
        const TString outputFile
) {
  gBenchmark->Start("makeGenMllTemplate_forMuon");

//  TCanvas *c = new TCanvas("c","",800,600);
  TH1D *h_mll_inclusive = new TH1D("h_mll_inclusive", "", massBin, fitMassLo, fitMassHi);
  TH1D *h_mll_perBin[muEtaNB][muPtNB];

  for(int iEtaBin = 0; iEtaBin < muEtaNB; iEtaBin ++){
    for(int iPtBin = 0; iPtBin < muPtNB; iPtBin ++){
      h_mll_perBin[iEtaBin][iPtBin] = new TH1D(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin), "", massBin, fitMassLo, fitMassHi);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t MUON_MASS  = 0.000511;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;
  
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

  // Data structures to store info from TTrees
  //baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  //TClonesArray *muonArr        = new TClonesArray("baconhep::TMuon");
  //TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
   
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv;
  vector<Double_t> nSelCorrv, nSelCorrVarv;
  vector<Double_t> accv, accCorrv;
  vector<Double_t> accErrv, accErrCorrv;

//  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    //eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr  = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch* genPartBr = eventTree->GetBranch("GenParticle");
    //eventTree->SetBranchAddress("Muon",          &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");   
    //eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    //
    // loop over events
    //      

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<10000; ientry++) {
      if(ientry%1000000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
      //infoBr->GetEntry(ientry);
      //muonArr->Clear(); muonBr->GetEntry(ientry);

      Int_t glepq1=-99;
      Int_t glepq2=-99;

      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,&glepq1,&glepq2,1);
      //if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;
      if(lep1->Pt() < PT_CUT || lep2->Pt() < PT_CUT) continue;
      if(abs(lep1->Eta()) > ETA_CUT || abs(lep2->Eta()) > ETA_CUT) continue;

//      delete vec; delete lep1; delete lep2;

      TLorentzVector vDilep = *lep1 + *lep2;
      Double_t mass = vDilep.M();
      if((mass<MASS_LOW) || (mass>MASS_HIGH)) continue;

      h_mll_inclusive->Fill(mass);

      Float_t eta = -99.;
      Float_t pt  = -999.;
      if(glepq1 > 0){
        eta = lep1->Eta();
        pt  = lep1->Pt();
      }
      else{
        eta = lep2->Eta();
        pt  = lep2->Pt();
      }

      bool fillFlag = false;
      for(int iEtaBin = 0; iEtaBin < muEtaNB; iEtaBin ++){
        for(int iPtBin = 0; iPtBin < muPtNB; iPtBin ++){
          if( eta > muEtaRange[iEtaBin] && eta < muEtaRange[iEtaBin+1] && pt > muPtRange[iPtBin] && pt < muPtRange[iPtBin+1]) {
            h_mll_perBin[iEtaBin][iPtBin]->Fill(mass);
            fillFlag = true;
          }
        }
      }
      if(!fillFlag) cout<<"[ERROR!] This should not happen, this didn't fill per-Bin mll distribution..."<<endl;
    }

    delete infile;
    infile=0, eventTree=0;  
  }  

  TFile *f = new TFile(TString(outputDir+"/"+outputFile), "RECREATE");
  h_mll_inclusive->Write();
  for(int iEtaBin = 0; iEtaBin < muEtaNB; iEtaBin ++){
    for(int iPtBin = 0; iPtBin < muPtNB; iPtBin ++){
      h_mll_perBin[iEtaBin][iPtBin]->Write();
    }
  }
  f->Close();
  
  //delete info;
  delete gen;
  //delete muonArr;

  gBenchmark->Show("makeGenMllTemplate_forMuon");
}
