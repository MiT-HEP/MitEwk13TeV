#if !defined(__CINT__) //|| defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include "TLorentzVector.h"

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace std;

#endif

void makeFlat(TString input="root://eoscms.cern.ch//store/user/jlawhorn/NNPDF30-GENSIM/wmm-pythia8-powheg1-bacon.root",
	      TString output="test.root") {
  
  TChain chain("Events");
  chain.Add(input);
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &part);        TBranch *partBr     = chain.GetBranch("GenParticle");

  //PDF info
  Double_t id_1,      id_2,       x_1,        x_2;
  Double_t xPDF_1,    xPDF_2,     scalePDF,   weight;

  //Generator level V+l info
  Double_t genV_id,   genL1_id,   genL2_id;
  Double_t genV_pt,   genV_eta,   genV_phi,   genV_m;
  Double_t genVf_pt,  genVf_eta,  genVf_phi,  genVf_m;
  Double_t genL1_pt,  genL1_eta,  genL1_phi,  genL1_m;
  Double_t genL2_pt,  genL2_eta,  genL2_phi,  genL2_m;
  Double_t genL1f_pt, genL1f_eta, genL1f_phi, genL1f_m;
  Double_t genL2f_pt, genL2f_eta, genL2f_phi, genL2f_m;
  
  TFile *ofile = new TFile(output, "recreate");
  TTree *otree = new TTree("Events", "Events");

  otree->Branch("id_1",       &id_1,       "id_1/d");
  otree->Branch("id_2",       &id_2,       "id_2/d");
  otree->Branch("x_1",        &x_1,        "x_1/d");
  otree->Branch("x_2",        &x_2,        "x_2/d");
  otree->Branch("xPDF_1",     &xPDF_1,     "xPDF_1/d");
  otree->Branch("xPDF_2",     &xPDF_2,     "xPDF_2/d");
  otree->Branch("scalePDF",   &scalePDF,   "scalePDF/d"); 
  otree->Branch("weight",     &weight,     "weight/d"); 
  otree->Branch("genV_pt",    &genV_pt,    "genV_pt/d");
  otree->Branch("genV_eta",   &genV_eta,   "genV_eta/d");
  otree->Branch("genV_phi",   &genV_phi,   "genV_phi/d");
  otree->Branch("genV_m",     &genV_m,     "genV_m/d");
  otree->Branch("genV_id",    &genV_id,    "genV_id/d");
  otree->Branch("genVf_pt",   &genVf_pt,   "genVf_pt/d");
  otree->Branch("genVf_eta",  &genVf_eta,  "genVf_eta/d");
  otree->Branch("genVf_phi",  &genVf_phi,  "genVf_phi/d");
  otree->Branch("genVf_m",    &genVf_m,    "genVf_m/d");
  otree->Branch("genL1_pt",   &genL1_pt,   "genL1_pt/d");
  otree->Branch("genL1_eta",  &genL1_eta,  "genL1_eta/d");
  otree->Branch("genL1_phi",  &genL1_phi,  "genL1_phi/d");
  otree->Branch("genL1_m",    &genL1_m,    "genL1_m/d");
  otree->Branch("genL1_id",   &genL1_id,   "genL1_id/d");
  otree->Branch("genL2_pt",   &genL2_pt,   "genL2_pt/d");
  otree->Branch("genL2_eta",  &genL2_eta,  "genL2_eta/d");
  otree->Branch("genL2_phi",  &genL2_phi,  "genL2_phi/d");
  otree->Branch("genL2_m",    &genL2_m,    "genL2_m/d");
  otree->Branch("genL2_id",   &genL2_id,   "genL2_id/d");
  otree->Branch("genL1f_pt",  &genL1f_pt,  "genL1f_pt/d");
  otree->Branch("genL1f_eta", &genL1f_eta, "genL1f_eta/d");
  otree->Branch("genL1f_phi", &genL1f_phi, "genL1f_phi/d");
  otree->Branch("genL1f_m",   &genL1f_m,   "genL1f_m/d");
  otree->Branch("genL2f_pt",  &genL2f_pt,  "genL2f_pt/d");
  otree->Branch("genL2f_eta", &genL2f_eta, "genL2f_eta/d");
  otree->Branch("genL2f_phi", &genL2f_phi, "genL2f_phi/d");
  otree->Branch("genL2f_m",   &genL2f_m,   "genL2f_m/d");

  for (Int_t i=0; i<chain.GetEntries(); i++) {
  //for (Int_t i=0; i<100; i++) {
    infoBr->GetEntry(i);
    part->Clear(); partBr->GetEntry(i);
    if (part->GetEntries()==0) continue;
    Int_t iv=-1, il1=-1, il2=-1;

    genV_id=0; genL1_id=0; genL2_id=0;
    genV_pt=0; genV_eta=0; genV_phi=0; genV_m=0;
    genVf_pt=0; genVf_eta=0; genVf_phi=0; genVf_m=0;
    genL1_pt=0; genL1_eta=0; genL1_phi=0; genL1_m=0;
    genL1f_pt=0; genL1f_eta=0; genL1f_phi=0; genL1f_m=0;
    genL2_pt=0; genL2_eta=0; genL2_phi=0; genL2_m=0;
    genL2f_pt=0; genL2f_eta=0; genL2f_phi=0; genL2f_m=0;
    
    for (Int_t j=0; j<part->GetEntries(); j++) { 
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
      
      if ((fabs(genloop->pdgId)==24||fabs(genloop->pdgId)==23) && (genloop->status==3||genloop->status==22)) {
	//cout << "found " << genloop->pdgId << " at " << j << endl;
	genV_pt  = genloop->pt;
	genV_eta = genloop->eta;
	genV_phi = genloop->phi;
	genV_m   = genloop->mass;
	genV_id  = genloop->pdgId;
	iv=j;
      }
      else if (iv!=-1 && genloop->parent==iv) {
	if (genloop->pdgId==genV_id) { 
	  genVf_pt  = genloop->pt;
	  genVf_eta = genloop->eta;
	  genVf_phi = genloop->phi;
	  genVf_m   = genloop->mass;
	  iv=j; //cout << "found daughter " << genloop->pdgId << " at " << j << " with parent " << genloop->parent << endl; 
	}
	else if (genloop->pdgId==-13 || genloop->pdgId==-11) {
	  //cout << "found -l at " << j << " with parent " << genloop->parent << endl;
	  genL1_pt  = genloop->pt;
	  genL1_eta = genloop->eta;
	  genL1_phi = genloop->phi;
	  genL1_m   = genloop->mass;
	  genL1_id  = genloop->pdgId;
	  il1=j;
	}
	else if (genloop->pdgId==13 || genloop->pdgId==11) {
	  //cout << "found +l at " << j << " with parent " << genloop->parent << endl;
	  genL2_pt  = genloop->pt;
	  genL2_eta = genloop->eta;
	  genL2_phi = genloop->phi;
	  genL2_m   = genloop->mass;
	  genL2_id  = genloop->pdgId;
	  il2=j;
	}
      }
      else if (il1!=-1 && genloop->parent==il1 && genloop->pdgId==genL1_id) {
	//cout << "found daughter -l at " << j << " with parent " << genloop->parent << endl;
	genL1f_pt  = genloop->pt;
	genL1f_eta = genloop->eta;
	genL1f_phi = genloop->phi;
	genL1f_m   = genloop->mass;
	il1=j;
      }
      else if (il2!=-1 && genloop->parent==il2 && genloop->pdgId==genL2_id) {
	//cout << "found daughter +l at " << j << " with parent " << genloop->parent << endl;
	genL2f_pt  = genloop->pt;
	genL2f_eta = genloop->eta;
	genL2f_phi = genloop->phi;
	genL2f_m   = genloop->mass;
	il2=j;
      }
    }

    if (genL1f_m==0 && genL1_m!=0) {
      genL1f_pt  = genL1_pt;
      genL1f_eta = genL1_eta;
      genL1f_phi = genL1_phi;
      genL1f_m   = genL1_m;
    }

    if (genL2f_m==0 && genL2_m!=0) {
      genL2f_pt  = genL2_pt;
      genL2f_eta = genL2_eta;
      genL2f_phi = genL2_phi;
      genL2f_m   = genL2_m;
    }

    id_1 = info->id_1;
    id_2 = info->id_2;
    x_1 = info->x_1;
    x_2 = info->x_2;
    xPDF_1 = info->xPDF_1;
    xPDF_2 = info->xPDF_2;
    scalePDF = info->scalePDF;
    weight = info->weight;

    otree->Fill();
    
  }
  
  ofile->Write();
  ofile->Close();
  
}
