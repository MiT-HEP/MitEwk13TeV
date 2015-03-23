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
#include "Eff/Eff.hh"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"

//#include "MitEwk13TeV/Ntupler/interface/TGenInfo.hh"
//#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
//#include "BaconAna/DataFormats/interface/TGenJet.hh"
//#include "BaconAna/DataFormats/interface/TGenParticle.hh"


using namespace std;

#endif

//#typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void accWme_test(TString input="/afs/cern.ch/work/j/jlawhorn/phil_wme_select_2.root",
	    TString pdfName="CT10nlo_as_0115",
	    Int_t iPdfSet=0) {

  LHAPDF::setVerbosity(LHAPDF::SILENT);

  //LHAPDF::PDF* nomPdf = LHAPDF::mkPDF("CT10nlo",0);
  LHAPDF::PDF* testPdf = LHAPDF::mkPDF(pdfName.Data(),iPdfSet);
  
  TChain chain("Events");
  chain.Add(input);

  const Int_t    PDG_ID    = 11;

  Int_t npv, q;
  Float_t x1,x2,w,weff;
  Int_t   id1,id2;
  Float_t scalePDF;
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >  *lep=0;
  
  chain.SetBranchAddress("npv"      ,&npv);
  chain.SetBranchAddress("q"        ,&q);
  chain.SetBranchAddress("lep"      ,&lep);
  chain.SetBranchAddress("w"        ,&w);
  chain.SetBranchAddress("weff"     ,&weff);
  chain.SetBranchAddress("x1"       ,&x1);
  chain.SetBranchAddress("x2"       ,&x2);
  chain.SetBranchAddress("id1"      ,&id1);
  chain.SetBranchAddress("id2"      ,&id2);
  chain.SetBranchAddress("scalePDF" ,&scalePDF);

  Float_t nPost=0;
  Float_t nTotal=0;

  Double_t nPost_d=0;
  Double_t nTotal_d=0;

  TH1D* tAll=new TH1D("tAll", "", 10, -10, 10000);
  TH1D* tPost=new TH1D("tPost", "", 10, -10, 10000);
  
  //for (Int_t i=0; i<chain.GetEntries(); i++) {
  for (Int_t i=0; i<10; i++) {
    chain.GetEntry(i);

    if (id1==0) id1=21;
    if (id2==0) id2=21;

    cout << w << " " << x1 << " " << x2 << " " << id1 << " " << id2 << " " << scalePDF << " " << testPdf->xfxQ(id1, x1, scalePDF)*testPdf->xfxQ(id2, x2, scalePDF)/(x1*x2) << endl;
    
    Float_t weight = testPdf->xfxQ(id1, x1, scalePDF)*testPdf->xfxQ(id2, x2, scalePDF)/(x1*x2)/w;
    Double_t weight_d = testPdf->xfxQ(id1, x1, scalePDF)*testPdf->xfxQ(id2, x2, scalePDF)/(x1*x2)/w;

    nTotal+=weight;
    nTotal_d+=weight_d;

    nPost+=weight*weff;
    nPost_d+=weight_d*weff;

    tAll->Fill(npv,weight);
    tPost->Fill(npv, weight*weff);

  }

  //cout << tPost->GetBinContent(0) << " " << tPost->GetBinContent(10000) << endl;
  //cout << nPost << " " << nPost_d << " " << nTotal << " " << nTotal_d << " " << chain.GetEntries() << endl;
  cout << "TotF:  "  << nPost/nTotal*chain.GetEntries() << " " << nTotal/chain.GetEntries() << endl;
  cout << "TotD:  "  << nPost_d/nTotal_d*chain.GetEntries() << " " << nTotal_d/chain.GetEntries() << endl;
  //cout << "Tot: " << tPost->Integral()/tAll->Integral()*chain.GetEntries() << " " << tAll->Integral()/chain.GetEntries() << endl;
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
