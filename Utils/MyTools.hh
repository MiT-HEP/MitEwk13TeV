#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <iostream>
#include "TLorentzVector.h"
namespace toolbox 
{
  Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2);
  
  Double_t deltaPhi(const Double_t phi1, const Double_t phi2);
  
  Int_t roundToInt(const Double_t x);
  
  void fillGen(TClonesArray *genPartArr, Int_t vid, Int_t lid, TLorentzVector* &vec, TLorentzVector* &fvec, TLorentzVector* &lep1, TLorentzVector* &lep2);

}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);
    
  Double_t deta = eta1-eta2;
  
  return sqrt(dphi*dphi + deta*deta);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

//------------------------------------------------------------------------------------------------------------------------
Int_t toolbox::roundToInt(Double_t x)
{
  if(x>0)
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
  else
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

//------------------------------------------------------------------------------------------------------------------------
void toolbox::fillGen(TClonesArray *genPartArr, Int_t vid, Int_t lid, TLorentzVector* &vec, TLorentzVector* &fvec, TLorentzVector* &lep1, TLorentzVector* &lep2) 
{
  Int_t iv=-1, iv1=-1, iv2=-1;
  TLorentzVector *lepPos=0, *lepNeg=0;

  for (Int_t i=0; i<genPartArr->GetEntries(); i++) {
    const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
    if (fabs(genloop->pdgId==vid) && (genloop->status==3||genloop->status==22)) {
      vec=new TLorentzVector(0,0,0,0);
      vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
      iv=i;
    }
    else if (iv!=-1 && genloop->parent==iv) {
      if (fabs(genloop->pdgId)==vid) {
        fvec=new TLorentzVector(0,0,0,0);
	fvec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
        iv=i;
      }
      if(fabs(genloop->pdgId)==lid) {
	if (genloop->pdgId<0 && lepPos==0) {
          lepPos=new TLorentzVector(0,0,0,0);
          lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
          iv1=i;
	}
        else if (lepNeg==0) {
          lepNeg=new TLorentzVector(0,0,0,0);
          lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
          iv2=i;
        }
      }
    }
    else if (fabs(genloop->pdgId)==lid) {
      if (genloop->parent==iv1) {
        lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
        iv1=i;
      }
      if (genloop->parent==iv2) {
        lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
        iv2=i;
      }
    }
  }
  if (lepPos && lepNeg) {
    if (lepPos->Pt()>lepNeg->Pt()) {
      lep1=lepPos;
      lep2=lepNeg;
    }
    else {
      lep1=lepNeg;
      lep2=lepPos;
    }
  }
  else if (lepPos) lep1=lepPos;
  else if (lepNeg) lep1=lepNeg;
  
}

#endif
