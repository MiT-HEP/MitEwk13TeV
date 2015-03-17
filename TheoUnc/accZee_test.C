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
//#include "LHAPDF/LHAPDF.h"

#include "MitEwk13TeV/Ntupler/interface/TGenInfo.hh"
//#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
//#include "BaconAna/DataFormats/interface/TGenJet.hh"
//#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace std;

#endif

void accZee_test(TString input="root://eoscms.cern.ch//store/user/jlawhorn/Ewk8TeV/s12-zee-v9-s8_scntuple.root") {//,
		 //TString input="root://eoscms.cern.ch//store/user/jlawhorn/PYTHIA-CT10nlo/Zee_bacon.root",
		 //TString pdfName="CT10nlo",
		 //Int_t iPdfSet=0) {

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.5;
  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;

  //LHAPDF::setVerbosity(LHAPDF::SILENT);
  
  //LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(pdfName.Data(),0);
  //LHAPDF::PDF* testPdf = LHAPDF::mkPDF(pdfName.Data(),iPdfSet);

  //const Int_t    PDG_ID    = 11;
  const Double_t ELE_MASS  = 0.000511;

  //Float_t nPreBB=0, nPreBE=0, nPreEE=0;
  Float_t nPostBB=0, nPostBE=0, nPostEE=0;
  Float_t nTotal=0;

  TChain *fChain = new TChain("Events");
  fChain->Add(input);

  mithep::TGenInfo *gen = new mithep::TGenInfo();
  fChain->SetBranchAddress("Gen", &gen);
  TBranch *genBr = fChain->GetBranch("Gen");

  cout << fChain->GetEntries() << endl;

  for (Int_t i=0; i<fChain->GetEntries(); i++) {
  //for (Int_t i=0; i<10; i++) {
    fChain->GetEntry(i);
    genBr->GetEntry(i);

    Float_t weight = 1;//testPdf->xfxQ(id_1, info->x_1, info->scalePDF)*testPdf->xfxQ(id_2, info->x_2, info->scalePDF)/(nomPdf->xfxQ(id_1, info->x_1, info->scalePDF)*nomPdf->xfxQ(id_2, info->x_2, info->scalePDF));

    Bool_t isB1=(fabs(gen->eta_1)<ETA_BARREL) ? kTRUE : kFALSE;
    Bool_t isB2=(fabs(gen->eta_2)<ETA_BARREL) ? kTRUE : kFALSE;
    //Bool_t isE1=(fabs(post1->eta)>ETA_ENDCAP && fabs(post1->eta)<ETA_CUT) ? kTRUE : kFALSE;
    //Bool_t isE2=(fabs(post2->eta)>ETA_ENDCAP && fabs(post2->eta)<ETA_CUT) ? kTRUE : kFALSE;

    TLorentzVector ele1(0,0,0,0);
    ele1.SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, ELE_MASS);
    TLorentzVector ele2(0,0,0,0);
    ele2.SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, ELE_MASS);
    TLorentzVector zee = ele1+ele2;

    if (gen->vmass<MASS_LOW || gen->vmass>MASS_HIGH) continue;

    if ( fabs(gen->eta_1)>ETA_BARREL && fabs(gen->eta_1)<ETA_ENDCAP ) continue;
    if ( fabs(gen->eta_2)>ETA_BARREL && fabs(gen->eta_2)<ETA_ENDCAP ) continue;

    nTotal+=weight;

    if (gen->pt_1 < PT_CUT) continue;
    if (gen->pt_2 < PT_CUT) continue;
    if (fabs(gen->eta_1) > ETA_CUT) continue;
    if (fabs(gen->eta_2) > ETA_CUT) continue;
    if (zee.M()<MASS_LOW || zee.M()>MASS_HIGH) continue;

    if (isB1 && isB2) {
      //cout << "BB " << post1->eta << " " << post2->eta << endl;
      nPostBB+=weight;
    }
    else if (!isB1 && !isB2) {
      //cout << "EE " << post1->eta << " " << post2->eta << endl;
      nPostEE+=weight;
    }
    else {
      //cout << "BE " << post1->eta << " " << post2->eta << endl;
      nPostBE+=weight;
    }
    /*    if (pre1->pt>25 && pre2->pt>25  && fabs(pre1->eta)<2.5 && fabs(pre2->eta)<2.5) {
      Bool_t isB1=(fabs(pre1->eta)<1.4442) ? kTRUE : kFALSE;
      Bool_t isB2=(fabs(pre2->eta)<1.4442) ? kTRUE : kFALSE;
      Bool_t isE1=(fabs(pre1->eta)>1.566 && fabs(pre1->eta)<2.5) ? kTRUE : kFALSE;
      Bool_t isE2=(fabs(pre2->eta)>1.566 && fabs(pre2->eta)<2.5) ? kTRUE : kFALSE;

      TLorentzVector ele1(pre1->pt, pre1->eta, pre1->phi, ELE_MASS);
      TLorentzVector ele2(pre2->pt, pre2->eta, pre2->phi, ELE_MASS);
      TLorentzVector zee = ele1+ele2;

      if (fabs(zee.M())>60 && fabs(zee.M())<120) {
	if (isB1 && isB2) nPreBB+=weight;
	else if (isE1 && isE2) nPreEE+=weight;
	else if ( (isB1&&isE2)||(isB2&&isE1) ) nPreBE+=weight;
      }
    }
    */
  }

  cout << "pre-FSR currently broken" << endl;

  Float_t accBB=nPostBB/nTotal; 
  Float_t accBE=nPostBE/nTotal;
  Float_t accEE=nPostEE/nTotal;
  Float_t accT = (nPostBB+nPostBE+nPostEE)/nTotal;
  
  cout << "BB: " << accBB << " +/- " << sqrt(accBB*(1-accBB)/nTotal) << endl;
  cout << "BE: " << accBE << " +/- " << sqrt(accBE*(1-accBE)/nTotal) << endl;
  cout << "EE: " << accEE << " +/- " << sqrt(accEE*(1-accEE)/nTotal) << endl;
  cout << "Tot: " << accT << " +/- " << sqrt(accT*(1-accT)/nTotal) << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {
  
  const Float_t pi = 3.14159265358979;
  
  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);
  
  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
