#ifndef EWKANA_UTILS_LEPTONIDCUTS_HH
#define EWKANA_UTILS_LEPTONIDCUTS_HH

#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include <TMath.h>
#include <cassert>
		   
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passAntiMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passEleID(const baconhep::TElectron *electron, const Double_t rho=0);
Bool_t passEleLooseID(const baconhep::TElectron *electron, const Double_t rho=0);
Bool_t passAntiEleID(const baconhep::TElectron *electron, const Double_t rho=0);
Double_t getEffArea(const Double_t eta);

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho)
{
  if(muon->nTkLayers  < 6)            return kFALSE;
  if(muon->nPixHits   < 1)            return kFALSE;
  if(fabs(muon->d0)   > 0.2)          return kFALSE;
  if(fabs(muon->dz)   > 0.5)          return kFALSE;
  if(muon->muNchi2    > 10)           return kFALSE;
  if(muon->nMatchStn  < 2)            return kFALSE;
  if(muon->nValidHits < 1)            return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kGlobal)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
  Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
  if(iso > 0.12*(muon->pt)) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passAntiMuonID(const baconhep::TMuon *muon, const Double_t rho)
{
  if(muon->nTkLayers  < 6)            return kFALSE;
  if(muon->nPixHits   < 1)            return kFALSE;
  if(fabs(muon->d0)   > 0.2)          return kFALSE;
  if(fabs(muon->dz)   > 0.5)          return kFALSE;
  if(muon->muNchi2    > 10)           return kFALSE;
  if(muon->nMatchStn  < 2)            return kFALSE;
  if(muon->nValidHits < 1)            return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kGlobal)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
  Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
  if(iso < 0.3*(muon->pt)) return kFALSE;

  return kTRUE;
}


//--------------------------------------------------------------------------------------------------
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho)
{
  if(!(muon->typeBits & baconhep::EMuType::kGlobal) && !(muon->typeBits & baconhep::EMuType::kTracker)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
//  Double_t iso = muon->pfChIso04 + TMath::Max(muon->pfNeuIso04 + muon->pfGamIso04 - 0.5*(muon->puIso04),Double_t(0));
//  if(iso > 0.20*(muon->pt)) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const baconhep::TElectron *electron, const Double_t rho)
{ // CSA14 Medium working point
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & baconhep::EEleType::kEcalDriven)) return kFALSE;
    
  // conversion rejection
  if(electron->isConv)            return kFALSE;
     
  Double_t ea = getEffArea(electron->scEta);
  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
    if(iso > 0.13*(electron->pt))  return kFALSE;
    if(electron->nMissingHits > 2) return kFALSE;     

    if(electron->sieie  	  > 0.011)                        return kFALSE;
    if(fabs(electron->dPhiIn)     > 0.033)                        return kFALSE;
    if(fabs(electron->dEtaIn)     > 0.0011)                       return kFALSE;
    if(electron->hovere 	  > 0.091)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.034*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) > 0.012)                                return kFALSE;
    if(fabs(electron->dz) > 0.39)                                 return kFALSE;

  } else {
    // endcap
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
    if(iso > 0.13*(electron->pt)) return kFALSE;
    if(electron->nMissingHits > 1) return kFALSE;
     
    if(electron->sieie  	  > 0.031)                        return kFALSE;
    if(fabs(electron->dPhiIn)     > 0.047)                        return kFALSE;
    if(fabs(electron->dEtaIn)     > 0.024)                        return kFALSE;
    if(electron->hovere 	  > 0.099)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.086*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) > 0.016)                                return kFALSE;
    if(fabs(electron->dz) > 0.78)                                 return kFALSE;

  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passAntiEleID(const baconhep::TElectron *electron, const Double_t rho)
{ // CSA14 Medium working point
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & baconhep::EEleType::kEcalDriven)) return kFALSE;
  
  // conversion rejection
  if(electron->isConv)            return kFALSE;
     
  Double_t ea = getEffArea(electron->scEta);

  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
    if(iso > 0.13*(electron->pt))  return kFALSE;
    if(electron->nMissingHits > 2) return kFALSE;
     
    if(electron->sieie  	  > 0.011)                        return kFALSE;
    if(fabs(electron->dPhiIn)     > 0.033)                        return kFALSE;
    if(fabs(electron->dEtaIn)     > 0.0011)                       return kFALSE;
    if(electron->hovere 	  > 0.091)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.034*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) > 0.012)                                return kFALSE;
    if(fabs(electron->dz) > 0.39)                                 return kFALSE;

  } else {
    // endcap
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
    if(iso > 0.13*(electron->pt))  return kFALSE;
    if(electron->nMissingHits > 1) return kFALSE;
     
    if(electron->sieie  	  > 0.031)                        return kFALSE;
    if(fabs(electron->dPhiIn)     > 0.047)                        return kFALSE;
    if(fabs(electron->dEtaIn)     > 0.024)                        return kFALSE;
    if(electron->hovere 	  > 0.099)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.086*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) > 0.016)                                return kFALSE;
    if(fabs(electron->dz) > 0.78)                                 return kFALSE;

  }
  return kTRUE;
}


//--------------------------------------------------------------------------------------------------
Bool_t passEleLooseID(const baconhep::TElectron *electron, const Double_t rho)
{ // CSA14 Veto working point
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & baconhep::EEleType::kEcalDriven)) return kFALSE;

  // conversion rejection
  if(electron->isConv)            return kFALSE;
       
  Double_t ea = getEffArea(electron->scEta);
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
    if(iso > 0.24*(electron->pt))  return kFALSE;
    if(electron->nMissingHits > 2) return kFALSE;

    if(electron->sieie        > 0.012)                           return kFALSE;
    if(fabs(electron->dPhiIn) > 0.17)                            return kFALSE;
    if(fabs(electron->dEtaIn) > 0.0015)                          return kFALSE;
    if(electron->hovere       > 0.29)                            return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.23*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0)     > 0.037)                           return kFALSE;
    if(fabs(electron->dz)     > 0.82)                            return kFALSE;
  
  } else {
    // endcap
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
    if(iso > 0.21*(electron->pt))  return kFALSE;
    if(electron->nMissingHits > 3) return kFALSE;

    if(electron->sieie        > 0.038)                           return kFALSE;
    if(fabs(electron->dPhiIn) > 0.25)                            return kFALSE;
    if(fabs(electron->dEtaIn) > 0.0026)                          return kFALSE;
    if(electron->hovere       > 0.21)                            return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.17*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0)     > 0.28)                           return kFALSE;
    if(fabs(electron->dz)     > 0.95)                            return kFALSE;
  
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Double_t getEffArea(const Double_t eta) {
  
  if     (fabs(eta) < 1.0)   return 0.100;
  else if(fabs(eta) < 1.479) return 0.120;
  else if(fabs(eta) < 2.0)   return 0.085;
  else if(fabs(eta) < 2.2)   return 0.110;
  else if(fabs(eta) < 2.3)   return 0.120;
  else if(fabs(eta) < 2.4)   return 0.120;
  else                       return 0.130;
}
#endif

