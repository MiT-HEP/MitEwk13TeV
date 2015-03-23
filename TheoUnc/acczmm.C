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
#include <TLorentzVector.h>
#include "LHAPDF/LHAPDF.h"

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace std;

#endif

void acczmm(TString input="/afs/cern.ch/work/j/jlawhorn/PYTHIA-GEN/CT10nlo_rw/zmm-bacon.root",
	    TString pdfName="CT10nlo_as_0117",
	    Int_t iPdfSet=0) {

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.1;
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;

  LHAPDF::setVerbosity(LHAPDF::SILENT);

  LHAPDF::PDF* nomPdf = LHAPDF::mkPDF("CT10nlo",0);
  LHAPDF::PDF* testPdf = LHAPDF::mkPDF(pdfName.Data(),iPdfSet);
  
  TChain chain("Events");
  chain.Add(input);

  const Int_t    PDG_ID    = 13;
  const Double_t MUON_MASS = 0.105658369;
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &part);        TBranch *partBr     = chain.GetBranch("GenParticle");

  Float_t nPostBB=0, nPostBE=0, nPostEE=0;
  Float_t nTotal=0;
  
  for (Int_t i=0; i<chain.GetEntries(); i++) {
    infoBr->GetEntry(i);
    
    part->Clear(); partBr->GetEntry(i);
    
    Int_t iPost1=-1;
    Int_t iPost2=-1;

    if (part->GetEntries()==0) continue;
    
    for (Int_t j=0; j<part->GetEntries(); j++) { 
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);

      if (genloop->pdgId==PDG_ID) {
	if (genloop->status==1 && iPost1==-1) iPost1=j;
	else if (genloop->status==1 && genloop->parent==iPost1) iPost1=j;
      }
      
      if (genloop->pdgId==-PDG_ID) {
        if (genloop->status==1 && iPost2==-1) iPost2=j;
        else if (genloop->status==1 && genloop->parent==iPost2) iPost2=j;
      }
    }

    Int_t id_1 = ( (info->id_1==0) ? 21 : info->id_1 );
    Int_t id_2 = ( (info->id_2==0) ? 21 : info->id_2 );
    
    Float_t weight = testPdf->xfxQ(id_1, info->x_1, info->scalePDF)*testPdf->xfxQ(id_2, info->x_2, info->scalePDF)/(nomPdf->xfxQ(id_1, info->x_1, info->scalePDF)*nomPdf->xfxQ(id_2, info->x_2, info->scalePDF));

    const baconhep::TGenParticle* post1 = (baconhep::TGenParticle*) ((*part)[iPost1]);
    const baconhep::TGenParticle* post2 = (baconhep::TGenParticle*) ((*part)[iPost2]);

    Bool_t isB1=(fabs(post1->eta)<ETA_BARREL) ? kTRUE : kFALSE;
    Bool_t isB2=(fabs(post2->eta)<ETA_BARREL) ? kTRUE : kFALSE;

    TLorentzVector mu1(0,0,0,0);
    mu1.SetPtEtaPhiM(post1->pt, post1->eta, post1->phi, MUON_MASS);
    TLorentzVector mu2(0,0,0,0);
    mu2.SetPtEtaPhiM(post2->pt, post2->eta, post2->phi, MUON_MASS);
    TLorentzVector zmm = mu1+mu2;

    nTotal+=weight;

    if (post1->pt < PT_CUT) continue;
    if (post2->pt < PT_CUT) continue;
    if (fabs(post1->eta) > ETA_CUT) continue;
    if (fabs(post2->eta) > ETA_CUT) continue;
    if (zmm.M()<MASS_LOW || zmm.M()>MASS_HIGH) continue;

    if (isB1 && isB2) {
      nPostBB+=weight;
    }
    else if (!isB1 && !isB2) {
      nPostEE+=weight;
    }
    else {
      nPostBE+=weight;
    }
  }
  
  cout << "Tot: " << (nPostBB+nPostBE+nPostEE)/nTotal << " " << nTotal/chain.GetEntries() << endl;
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
