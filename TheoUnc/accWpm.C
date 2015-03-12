#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include "Math/LorentzVector.h"     // 4-vector class

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void accWpm(const TString input="/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN-SIM/Wpm_bacon.root") {
  
  TChain chain("Events");
  chain.Add(input);

  const Int_t    PDG_ID    = -13;
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &part);        TBranch *partBr     = chain.GetBranch("GenParticle");

  Float_t nPassPre=0;
  Float_t nPassPost=0;
  Float_t nTotal=0;
  
  for (Int_t i=0; i<chain.GetEntries(); i++) {
    infoBr->GetEntry(i);
    
    part->Clear(); partBr->GetEntry(i);
    
    Int_t iPre1=-1, iPost1=-1;

    if (part->GetEntries()==0) continue;
    
    for (Int_t j=0; j<part->GetEntries(); j++) { 
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);

      if (genloop->pdgId==PDG_ID) {
	if (genloop->status==3) iPre1=j;
	
	else if (genloop->status==1 && iPost1==-1) iPost1=j;
	else if (genloop->status==1 && genloop->parent==iPost1) iPost1=j;
      }       
    }

    nTotal+=1;

    if (iPre1==-1 || iPost1==-1) {
      cout << "?????" << endl;
      continue;
    }
    
    const baconhep::TGenParticle* pre1 = (baconhep::TGenParticle*) ((*part)[iPre1]);
    const baconhep::TGenParticle* post1 = (baconhep::TGenParticle*) ((*part)[iPost1]);

    if (pre1->pt>25 && fabs(pre1->eta)<2.1) {
      nPassPre+=1;
    }
    if (post1->pt>25 && fabs(post1->eta)<2.1) {
      nPassPost+=1;
    }
    
  }
  cout << "Pre-FSR:  " << nPassPre << ", " << nTotal << endl;
  cout << "Post-FSR: " << nPassPost << ", " << nTotal << endl;
  cout << "Pre-FSR:  " << nPassPre/nTotal << endl;
  cout << "Post-FSR: " << nPassPost/nTotal << endl;
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
