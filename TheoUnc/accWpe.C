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
#include "LHAPDF/LHAPDF.h"

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace std;

#endif

void accWpe(TString input="root://eoscms.cern.ch//store/user/jlawhorn/PYTHIA-CT10nlo/Wpe_bacon.root",
	    TString pdfName="CT10nlo",
	    Int_t iPdfSet=0) {

  LHAPDF::setVerbosity(LHAPDF::SILENT);

  LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(pdfName.Data(),0);
  LHAPDF::PDF* testPdf = LHAPDF::mkPDF(pdfName.Data(),iPdfSet);
  
  TChain chain("Events");
  chain.Add(input);

  const Int_t    PDG_ID    = -11;
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &part);        TBranch *partBr     = chain.GetBranch("GenParticle");

  Float_t nPreBarr=0;
  Float_t nPreEnd=0;
  Float_t nPostBarr=0;
  Float_t nPostEnd=0;
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

    Int_t id_1 = ( (info->id_1==0) ? 21 : info->id_1 );
    Int_t id_2 = ( (info->id_2==0) ? 21 : info->id_2 );

    Float_t weight = testPdf->xfxQ(id_1, info->x_1, info->scalePDF)*testPdf->xfxQ(id_2, info->x_2, info->scalePDF)/(nomPdf->xfxQ(id_1, info->x_1, info->scalePDF)*nomPdf->xfxQ(id_2, info->x_2, info->scalePDF));

    nTotal+=weight;

    if (iPre1==-1 || iPost1==-1) {
      cout << "?????" << endl;
      continue;
    }
    
    const baconhep::TGenParticle* pre1 = (baconhep::TGenParticle*) ((*part)[iPre1]);
    const baconhep::TGenParticle* post1 = (baconhep::TGenParticle*) ((*part)[iPost1]);

    if (pre1->pt>25 && fabs(pre1->eta)<1.4442) nPreBarr+=weight;
    else if (pre1->pt>25 && fabs(pre1->eta)>1.566 && fabs(pre1->eta)<2.5) nPreEnd+=weight;
    
    if (post1->pt>25 && fabs(post1->eta)<1.4442) nPostBarr+=weight;
    else if (post1->pt>25 && fabs(post1->eta)>1.566 && fabs(post1->eta)<2.5) nPostEnd+=weight;

  }
  
  Float_t accB=nPostBarr/nTotal;
  Float_t accE=nPostEnd/nTotal;
  Float_t accT=(nPostBarr+nPostEnd)/nTotal;
  
  cout << "B: " << accB << " +/- " << sqrt((1-accB)*(accB)/nTotal) << " ( " << nPostBarr << "/" << nTotal << " )" << endl;
  cout << "E: " << accE << " +/- " << sqrt((1-accE)*(accE)/nTotal) << " ( " << nPostEnd << "/" << nTotal << " )" << endl;
  cout << "Tot:  "  << accT << " +/- " << sqrt((1-accT)*(accT)/nTotal) << " ( " << nPostBarr+nPostEnd << "/" << nTotal << " )" << endl;
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
