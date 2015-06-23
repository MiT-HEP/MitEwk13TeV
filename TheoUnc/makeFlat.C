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

void makeFlat(TString input="root://eoscms.cern.ch//store/user/jlawhorn/Ewk13TeV/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
	      TString output="/afs/cern.ch/work/j/jlawhorn/theo-unc-06-10/WmJetsToLNu.root",
	      Int_t vid=-24) {
  
  TChain chain("Events");
  chain.Add(input);
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr      = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &genPartArr);  TBranch *partBr     = chain.GetBranch("GenParticle");

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

  for (Int_t ie=0; ie<chain.GetEntries(); ie++) {
  //for (Int_t ie=0; ie<10000000; ie++) {
  //for (Int_t ie=0; ie<100; ie++) {
    infoBr->GetEntry(ie);
    genPartArr->Clear(); partBr->GetEntry(ie);
    if (genPartArr->GetEntries()==0) continue;
    TLorentzVector *vec=0, *lepPos=0, *lepNeg=0;
    TLorentzVector *preVec=0, *preLepPos=0, *preLepNeg=0;
    Int_t flavor=0;
    Int_t iv=-1, iv1=-1, iv2=-1;
    Int_t vidLoop=vid;
    //cout << "-----" << endl;
    for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
      //cout << i << ", " << genloop->pdgId << ", " << genloop->parent << ", " << genloop->pt << ", " << genloop->status << endl;
      if (genloop->pdgId==-vidLoop) {
	vec=0; lepPos=0; lepNeg=0;
	vidLoop=-vid;
	break;
      }
      if (genloop->status==23 && (fabs(genloop->pdgId)==15 || fabs(genloop->pdgId)==13 || fabs(genloop->pdgId)==11)) {
	if (flavor==0) {
	  flavor=genloop->pdgId;
	}
	if (genloop->pdgId<0 && lepPos==0) {
	  lepPos=new TLorentzVector(0,0,0,0);
	  lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	  preLepPos=new TLorentzVector(0,0,0,0);
	  preLepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	  iv1=i;
	}
	else if (genloop->pdgId>0 && lepNeg==0) {
	  lepNeg=new TLorentzVector(0,0,0,0);
	  lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	  preLepNeg=new TLorentzVector(0,0,0,0);
	  preLepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	  iv2=i;
	}
      }
      else if (genloop->pdgId==vid && (genloop->status==3||genloop->status==22)) {
	preVec=new TLorentzVector(0,0,0,0);
	preVec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	vec=new TLorentzVector(0,0,0,0);
	vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	iv=i;
      }
      else if (iv!=-1 && genloop->parent==iv) {
	if (genloop->pdgId==vid) {
	  vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	  iv=i;
	}
	else if (fabs(genloop->pdgId)==15 || fabs(genloop->pdgId)==13 || fabs(genloop->pdgId)==11) {
	  if (flavor==0) {
	    flavor=genloop->pdgId;
	  }
	  if (genloop->pdgId<0 && lepPos==0) {
	    lepPos=new TLorentzVector(0,0,0,0);
	    lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    preLepPos=new TLorentzVector(0,0,0,0);
	    preLepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    iv1=i;
	  }
	  else if (genloop->pdgId>0 && lepNeg==0) {
	    lepNeg=new TLorentzVector(0,0,0,0);
	    lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    preLepNeg=new TLorentzVector(0,0,0,0);
	    preLepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	    iv2=i;
	  }
	}
      }
      else if (iv1!=-1 && genloop->parent==iv1) {
	lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	iv1=i;
      }
      else if (iv2!=-1 && genloop->parent==iv2) {
	lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	iv2=i;
      }
    }

    //if (vec==0 && lepPos==0 && lepNeg==0) continue;

    if (vec==0 && preLepNeg && preLepPos) {
      TLorentzVector temp = *preLepNeg + *preLepPos;
      vec=&temp;
    }
    genV_id=vidLoop;
    genL1_id=-fabs(flavor);
    genL2_id=fabs(flavor);
    
    if (preVec && preVec->M()>0) {
      genV_pt  = preVec->Pt();
      genV_eta = preVec->Eta();
      genV_phi = preVec->Phi();
      genV_m   = preVec->M();
    }
    else if (vec) {
      genV_pt  = vec->Pt();
      genV_eta = vec->Eta();
      genV_phi = vec->Phi();
      genV_m   = vec->M();
    }

    if (lepPos && lepNeg) {
      TLorentzVector temp = *lepPos + *lepNeg;
      genVf_pt  = temp.Pt();
      genVf_eta = temp.Eta();
      genVf_phi = temp.Phi();
      genVf_m   = temp.M();
    }

    if (lepPos) {
      genL1f_pt  = lepPos->Pt();
      genL1f_eta = lepPos->Eta();
      genL1f_phi = lepPos->Phi();
      genL1f_m   = lepPos->M();
    }
    else {
      genL1f_pt  = 0;
      genL1f_eta = 0;
      genL1f_phi = 0;
      genL1f_m   = 0;
    }

    if (lepNeg) {
      genL2f_pt  = lepNeg->Pt();
      genL2f_eta = lepNeg->Eta();
      genL2f_phi = lepNeg->Phi();
      genL2f_m   = lepNeg->M();
    }
    else {
      genL2f_pt  = 0;
      genL2f_eta = 0;
      genL2f_phi = 0;
      genL2f_m   = 0;
    }

    if (preLepPos) {
      genL1_pt  = preLepPos->Pt();
      genL1_eta = preLepPos->Eta();
      genL1_phi = preLepPos->Phi();
      genL1_m   = preLepPos->M();
    }
    else {
      genL1_pt  = 0;
      genL1_eta = 0;
      genL1_phi = 0;
      genL1_m   = 0;
    }

    if (preLepNeg) {
      genL2_pt  = preLepNeg->Pt();
      genL2_eta = preLepNeg->Eta();
      genL2_phi = preLepNeg->Phi();
      genL2_m   = preLepNeg->M();
    }
    else {
      genL2_pt  = 0;
      genL2_eta = 0;
      genL2_phi = 0;
      genL2_m   = 0;
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
