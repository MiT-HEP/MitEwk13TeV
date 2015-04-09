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

void accwmm_deta(TString input="root://eoscms.cern.ch//store/user/jlawhorn/PYTHIA-CT10nlo-RW/wmm-bacon.root",
		 TString outputDir="/afs/cern.ch/work/j/jlawhorn/public/wz-envelopes/",
		 TString pdfName="CT10nlo_as_0117",
		 Int_t iPdfSet=0) {

  TString output=outputDir+"wmm_"+pdfName+".root";

  LHAPDF::setVerbosity(LHAPDF::SILENT);

  LHAPDF::PDF* nomPdf = LHAPDF::mkPDF("CT10nlo",0);
  LHAPDF::PDF* testPdf = LHAPDF::mkPDF(pdfName.Data(),iPdfSet);
  
  TChain chain("Events");
  chain.Add(input);

  const Int_t    PDG_ID    = 13;
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenParticle", &part);        TBranch *partBr     = chain.GetBranch("GenParticle");

  TFile *outFile = new TFile(output, "update");
  
  char histname[100];
  sprintf(histname,"dEta_%s_%i",pdfName.Data(),iPdfSet);
  
  TH1D* dEta = new TH1D(histname, "", 20, -5, 5); dEta->Sumw2();

  sprintf(histname,"dPt_%s_%i",pdfName.Data(),iPdfSet);
  
  TH1D* dPt = new TH1D(histname, "", 25, 25, 100); dPt->Sumw2();

  Double_t nPostBarr=0;
  Double_t nPostEnd=0;
  Double_t nTotal=0;
  
  for (Int_t i=0; i<chain.GetEntries(); i++) {
    infoBr->GetEntry(i);
    
    part->Clear(); partBr->GetEntry(i);
    
    Int_t iPost1=-1;

    if (part->GetEntries()==0) continue;
    
    for (Int_t j=0; j<part->GetEntries(); j++) { 
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);

      if (genloop->pdgId==PDG_ID) {
	if (genloop->status==1 && iPost1==-1) iPost1=j;
	else if (genloop->status==1 && genloop->parent==iPost1) iPost1=j;
      }       
    }
    
    Int_t id_1 = ( (info->id_1==0) ? 21 : info->id_1 );
    Int_t id_2 = ( (info->id_2==0) ? 21 : info->id_2 );
    
    Double_t weight = testPdf->xfxQ(id_1, info->x_1, info->scalePDF)*testPdf->xfxQ(id_2, info->x_2, info->scalePDF)/(nomPdf->xfxQ(id_1, info->x_1, info->scalePDF)*nomPdf->xfxQ(id_2, info->x_2, info->scalePDF));

    nTotal+=weight;

    const baconhep::TGenParticle* post1 = (baconhep::TGenParticle*) ((*part)[iPost1]);
    
    if (post1->pt>25 && fabs(post1->eta)<1.2) nPostBarr+=weight;
    else if (post1->pt>25 && fabs(post1->eta)>1.2 && fabs(post1->eta)<2.1) nPostEnd+=weight;

    dEta->Fill(post1->eta, weight);
    dPt->Fill(post1->pt, weight);
    
  }

  cout << "Tot:  "  << (nPostBarr+nPostEnd)/nTotal << " " << nTotal/chain.GetEntries() << endl;

  outFile->Write();
  outFile->Close();

}
