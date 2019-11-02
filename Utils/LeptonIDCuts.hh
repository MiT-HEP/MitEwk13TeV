#ifndef EWKANA_UTILS_LEPTONIDCUTS_HH
#define EWKANA_UTILS_LEPTONIDCUTS_HH

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include <TMath.h>
#include <cassert>
		   
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonID_2015(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passAntiMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passAntiMuonID_2015(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonLooseID_2015(const baconhep::TMuon *muon, const Double_t rho=0);

Bool_t passEleID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passEleTightID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passEleTightID_2015(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passEleLooseID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passEleLooseID_2015(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passAntiEleTightID(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho=0);
Bool_t passAntiEleTightID_2015(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho=0);
Bool_t passAntiEleID(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho=0);

Bool_t isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits,Bool_t isData);
Bool_t isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isData);

Bool_t isMuonTriggerNoIso(baconhep::TTrigger triggerMenu, TriggerBits hltBits);
Bool_t isMuonTriggerObjNoIso(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1);

Bool_t isEleTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Bool_t isData);
Bool_t isEleTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Bool_t isData);

Double_t getEffAreaEl(const Double_t eta);
Double_t getEffAreaEl_2015(const Double_t eta);
Double_t getEffAreaMu(const Double_t eta);
Double_t getEffAreaMu_2015(const Double_t eta);

//=== FUNCTION DEFINITIONS ======================================================================================
//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho)
{ // 2017 tight muon id is unchanged from 2015
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
  if(iso > 0.15*(muon->pt)) return kFALSE;

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
  if(iso < 0.15*(muon->pt)) return kFALSE;

  return kTRUE;
}


//--------------------------------------------------------------------------------------------------
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho)
{
  if(!(muon->typeBits & baconhep::EMuType::kGlobal) && !(muon->typeBits & baconhep::EMuType::kTracker)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
  return kTRUE;
}
// -------------------------------------------------------------------------------------------------
//            Muon ID from 2015
//--------------------------------------------------------------------------------------------------
Bool_t passMuonID_2015(const baconhep::TMuon *muon, const Double_t rho)
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
  if(iso > 0.15*(muon->pt)) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passAntiMuonID_2015(const baconhep::TMuon *muon, const Double_t rho)
{
  if(muon->nTkLayers  < 6)            return kFALSE;
  if(muon->nPixHits   < 1)            return kFALSE;
  if(fabs(muon->d0)   > 0.2)          return kFALSE;
  //if(fabs(muon->d0)   <= 0.2)          return kFALSE;
  if(fabs(muon->dz)   > 0.5)          return kFALSE;
  if(muon->muNchi2    > 10)           return kFALSE;
  if(muon->nMatchStn  < 2)            return kFALSE;
  if(muon->nValidHits < 1)            return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kGlobal)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
  Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
  if(iso < 0.15*(muon->pt)) return kFALSE;

  return kTRUE;
}


//--------------------------------------------------------------------------------------------------
Bool_t passMuonLooseID_2015(const baconhep::TMuon *muon, const Double_t rho)
{
  if(!(muon->typeBits & baconhep::EMuType::kGlobal) && !(muon->typeBits & baconhep::EMuType::kTracker)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
//                             Electron ID
//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Medium Electron ID for 2015

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    if(iso >= 0.0766*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.0103)                              return kFALSE;
    if(electron->hovere >= 0.0876)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0174*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0118)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.373)                                  return kFALSE;
  } else {
    if(iso >= 0.0678*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie  	  >= 0.0283)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.114)                           return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.00733)                         return kFALSE;
    if(electron->hovere 	  >= 0.0678)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0898*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0739)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.602)                                  return kFALSE;
  }

  return kTRUE;

}

Bool_t passAntiEleID(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho)
{// ele ID from 2015
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  // conversion rejection
  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());	
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  // barrel/endcap dependent requirements
  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    // barrel
    //if(iso >= 0.0766*(tag.Pt()))                                      return kFALSE;
    if(iso < 0.0766*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.0103)                              return kFALSE;
    if(electron->hovere >= 0.0876)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0174*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0118)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.373)                                  return kFALSE;
  } else {
    // endcap
    //if(iso >= 0.0678*(electron->pt))                                  return kFALSE;
    if(iso < 0.0678*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie  	  >= 0.0283)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.114)                           return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.00733)                         return kFALSE;
    if(electron->hovere 	  >= 0.0678)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0898*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0739)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.602)                                  return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passEleTightID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Tight Electron ID from 2017

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  // keeping these for now, can always remove at later step
  // if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    if(iso/(tag.Pt()) >= 0.0287+0.506/(tag.Pt()))                                      return kFALSE; // original
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie >= 0.0104)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.022)                              return kFALSE; // original
    if(fabs(electron->dEtaIn) >= 0.00255)                              return kFALSE; // original
    if(electron->hovere >= 0.026+1.15/electron->ecalEnergy+0.0324*rho/electron->ecalEnergy)          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.159*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.05)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.10)                                  return kFALSE;
  } else {
    if(iso/(tag.Pt()) >= 0.0445+0.963/(tag.Pt()))                                      return kFALSE; // original
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.0353)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.0236)                           return kFALSE; // original
    if(fabs(electron->dEtaIn)     >= 0.00501)                         return kFALSE; // original
    // std::cout << "----------------" << std::endl;
    if(electron->hovere           >= 0.0188+2.06/electron->ecalEnergy+0.183*rho/electron->ecalEnergy)         return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0197*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.10)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.20)                                 return kFALSE;
  }

  return kTRUE;

}

//--------------------------------------------------------------------------------------------------

Bool_t passAntiEleTightID(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho)
{// Tight Electron ID from 2017

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  // keeping these for now, can always remove later
  // if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
 // std::cout << "electron. dEta in " << fabs(electron->dEtaIn) << "  , dPhi in " << fabs(electron->dPhiIn)  << std::endl;
  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    // if(iso >= 0.0287+0.506/(tag.Pt()))                                      return kFALSE; // original
    if(iso/(tag.Pt()) < 0.0287+0.506/(tag.Pt()))                                      return kFALSE; // reverse
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie >= 0.0104)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.022)                              return kFALSE; // original
    if(fabs(electron->dEtaIn) >= 0.00255)                              return kFALSE; // original
    // if(fabs(electron->dPhiIn) < 0.022)                               return kFALSE; // reverse
    // if(fabs(electron->dEtaIn) < 0.00255)                              return kFALSE; // reverse
    if(electron->hovere >= 0.026+1.15/electron->ecalEnergy+0.0324*rho/electron->ecalEnergy)          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.159*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.05)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.10)                                  return kFALSE;
  } else {
    // if(iso >= 0.0445+0.963/(tag.Pt()))                                      return kFALSE; // original
    if(iso/(tag.Pt()) < 0.0445+0.963/(tag.Pt()))                                      return kFALSE;  // reverse
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.0353)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.0236)                           return kFALSE; // original
    if(fabs(electron->dEtaIn)     >= 0.00501)                         return kFALSE; // original
    // if(fabs(electron->dPhiIn)     < 0.0236)                           return kFALSE; // reverse
    // if(fabs(electron->dEtaIn)     < 0.00501)                          return kFALSE; // reverse 
    // std::cout << "----------------" << std::endl;
    if(electron->hovere           >= 0.0188+2.06/electron->ecalEnergy+0.183*rho/electron->ecalEnergy)         return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0197*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.10)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.20)                                 return kFALSE;
  }

  return kTRUE;
}

Bool_t passEleLooseID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Veto ID  with 2017 loose WP
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  // keeping these for now, can remove later?
  // if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;
  
  // conversion rejection
  if(electron->isConv) return kFALSE;
       
  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);     

  // barrel/endcap dependent requirements
  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    // barrel
    if(iso/(tag.Pt()) >= 0.198+0.506/(tag.Pt()))                                       return kFALSE;
    if(electron->nMissingHits >= 2)                                   return kFALSE;
    if(electron->sieie        >= 0.0126)                              return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.148)                               return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.00463 )                             return kFALSE;
    if(electron->hovere       >= 0.05+1.16/(electron->ecalEnergy)+0.0324*rho/(electron->ecalEnergy))  return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.209*(electron->ecalEnergy))    return kFALSE;
    if(fabs(electron->d0)     >= 0.05)                              return kFALSE;
    if(fabs(electron->dz)     >= 0.10)                                return kFALSE;
  } else {
    // endcap
    if(iso/(tag.Pt()) >= 0.203+0.963/(tag.Pt()))                                   return kFALSE;
    if(electron->nMissingHits >= 3)                                   return kFALSE;
    if(electron->sieie        >= 0.0457)                              return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.19)                               return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.00814)                             return kFALSE;
    if(electron->hovere       >= 0.05+2.54/(electron->ecalEnergy)+0.183*rho/(electron->ecalEnergy))  return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.132*(electron->ecalEnergy))    return kFALSE;
    if(fabs(electron->d0)     >= 0.10)                               return kFALSE;
    if(fabs(electron->dz)     >= 0.20)                               return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
//                       Ele 2015
//--------------------------------------------------------------------------------------------------
Bool_t passEleID_2015(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Medium Electron ID for 2015

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    if(iso >= 0.0766*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.0103)                              return kFALSE;
    if(electron->hovere >= 0.0876)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0174*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0118)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.373)                                  return kFALSE;
  } else {
    if(iso >= 0.0678*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie  	  >= 0.0283)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.114)                           return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.00733)                         return kFALSE;
    if(electron->hovere 	  >= 0.0678)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0898*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0739)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.602)                                  return kFALSE;
  }

  return kTRUE;

} 
Bool_t passEleTightID_2015(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Tight Electron ID from 2015

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    if(iso >= 0.0354*(tag.Pt()))                                      return kFALSE; // regular ISO
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.00926)                              return kFALSE;
    if(electron->hovere >= 0.0597)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.012*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0111)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.0466)                                  return kFALSE;
  } else {
    if(iso >= 0.0646*(tag.Pt()))                                      return kFALSE; // regular ISO
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.0279)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.0918)                           return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.00724)                         return kFALSE;
    if(electron->hovere           >= 0.0615)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.00999*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0351)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.417)                                  return kFALSE;
  }

  return kTRUE;

}
Bool_t passAntiEleTightID_2015(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho)
{// Tight Electron ID from 2015

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);
 // std::cout << "electron. dEta in " << fabs(electron->dEtaIn) << "  , dPhi in " << fabs(electron->dPhiIn)  << std::endl;
  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    if(iso >= 0.0354*(tag.Pt()))                                      return kFALSE; // original
    // if(iso < 0.0354*(tag.Pt()))                                      return kFALSE; // reverse
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    // if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE; // original
    // if(fabs(electron->dEtaIn) >= 0.00926)                              return kFALSE; // original
    if(fabs(electron->dPhiIn) < 0.0336)                               return kFALSE; // reverse
    if(fabs(electron->dEtaIn) < 0.00926)                              return kFALSE; // reverse
    // std::cout << "----------------" << std::endl;
    if(electron->hovere >= 0.0597)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.012*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0111)                                  std::cout << "FAIL d0 " << std::endl;//return kFALSE;
    if(fabs(electron->dz) >= 0.0466)                                  std::cout << "FAIL dz " << std::endl;//return kFALSE;
  } else {
    if(iso >= 0.0646*(tag.Pt()))                                      return kFALSE; // original
    // if(iso < 0.0646*(tag.Pt()))                                      return kFALSE;  // reverse
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.0279)                          return kFALSE;
    // if(fabs(electron->dPhiIn)     >= 0.0918)                           return kFALSE; // original
    // if(fabs(electron->dEtaIn)     >= 0.00724)                         return kFALSE; // original
    // if(fabs(electron->dPhiIn)     < 0.0918)                       std::cout << "FAIL PHI " << std::endl;//    return kFALSE; // reverse
    // if(fabs(electron->dEtaIn)     < 0.00724)                      std::cout << "FAIL ETA " << std::endl;//   return kFALSE; // reverse
    if(fabs(electron->dPhiIn)     < 0.0918)                           return kFALSE; // reverse
    if(fabs(electron->dEtaIn)     < 0.00724)                          return kFALSE; // reverse
    // std::cout << "----------------" << std::endl;
    if(electron->hovere           >= 0.0615)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.00999*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0351)                                  std::cout << "FAIL d0 " << std::endl;//return kFALSE;
    if(fabs(electron->dz) >= 0.417)                                  std::cout << "FAIL dz " << std::endl;//return kFALSE;
  }

  return kTRUE;
}

Bool_t passAntiEleID_2015(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho)
{
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  // if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

  // conversion rejection
  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());	
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  // barrel/endcap dependent requirements
  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    // barrel
    //if(iso >= 0.0766*(tag.Pt()))                                      return kFALSE;
    if(iso < 0.0766*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.0103)                              return kFALSE;
    if(electron->hovere >= 0.0876)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0174*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0118)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.373)                                  return kFALSE;
  } else {
    // endcap
    //if(iso >= 0.0678*(electron->pt))                                  return kFALSE;
    if(iso < 0.0678*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie  	  >= 0.0283)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.114)                           return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.00733)                         return kFALSE;
    if(electron->hovere 	  >= 0.0678)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0898*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0739)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.602)                                  return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------


Bool_t passEleLooseID_2015(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Veto ID  for 2015 Loose WP
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;
  
  // conversion rejection
  if(electron->isConv) return kFALSE;
       
  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);     

  // barrel/endcap dependent requirements
  if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    // barrel
    if(iso > 0.126*(tag.Pt()))                                       return kFALSE;
    if(electron->nMissingHits > 2)                                   return kFALSE;
    if(electron->sieie        > 0.0114)                              return kFALSE;
    if(fabs(electron->dPhiIn) > 0.216)                               return kFALSE;
    if(fabs(electron->dEtaIn) > 0.0152 )                             return kFALSE;
    if(electron->hovere       > 0.181)                               return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.207*(electron->ecalEnergy))    return kFALSE;
    if(fabs(electron->d0)     > 0.0564)                              return kFALSE;
    if(fabs(electron->dz)     > 0.472)                                return kFALSE;
  } else {
    // endcap
    if(iso > 0.144*(tag.Pt()))                                   return kFALSE;
    if(electron->nMissingHits > 3)                                   return kFALSE;
    if(electron->sieie        > 0.0352)                              return kFALSE;
    if(fabs(electron->dPhiIn) > 0.237)                               return kFALSE;
    if(fabs(electron->dEtaIn) > 0.0113)                             return kFALSE;
    if(electron->hovere       > 0.116)                              return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.174*(electron->ecalEnergy))    return kFALSE;
    if(fabs(electron->d0)     > 0.222)                               return kFALSE;
    if(fabs(electron->dz)     > 0.921)                               return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits,Bool_t isData, Bool_t is13TeV=true) {
  if(isData)
    return triggerMenu.pass("HLT_HIMu17_v*",hltBits);// both 5 & 13
  else{
    if(is13TeV) return triggerMenu.pass("HLT_Mu17_v*",hltBits); // 13 TeV
    else  return triggerMenu.pass("HLT_HIMu17_v*",hltBits); // 5 TeV
  }
}

Bool_t isMuonTriggerNoIso(baconhep::TTrigger triggerMenu, TriggerBits hltBits) {
  return triggerMenu.pass("HLT_HIMu17_v*",hltBits);
}

Bool_t isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isData, Bool_t is13TeV=true) {
  if (isData) return triggerMenu.passObj("HLT_HIMu17_v*","hltL3fL1sMu10lqL1f0L2f10L3Filtered17",hltMatchBits);// both 5 & 13
  else{ if(is13TeV) return triggerMenu.passObj("HLT_Mu17_v*","hltL3fL1sMu10lqL1f0L2f10L3Filtered17",hltMatchBits); // 13 TeV
    else return triggerMenu.passObj("HLT_HIMu17_v*","hltL3fL1sMu10lqL1f0L2f10L3Filtered17",hltMatchBits); // 5 TeV
  }
}

Bool_t isMuonTriggerObjNoIso(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1) {
  if (isL1) return triggerMenu.passObj("HLT_HIMu17_v*","hltL3fL1sMu10lqL1f0L2f10L3Filtered17",hltMatchBits);
  else return triggerMenu.passObj("HLT_Mu17_v*","hltL3fL1sMu10lqL1f0L2f10L3Filtered17",hltMatchBits);
}


Bool_t isEleTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Bool_t isData, Bool_t is13TeV=true) {
  if (isData) {
    return triggerMenu.pass("HLT_HIEle17_WPLoose_Gsf_v*",hltBits); // both 5 & 13
    // return triggerMenu.pass("HLT_HIEle20_WPLoose_Gsf_v*",hltBits); // both 5 & 13
    //return triggerMenu.pass("HLT_Ele27_WPLoose_Gsf_v*",hltBits);
  }
  else {
    if(is13TeV) return triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",hltBits); // 13 TeV
    // else return triggerMenu.pass("HLT_HIEle20_WPLoose_Gsf_v*",hltBits); // 5 TeV
    else return triggerMenu.pass("HLT_HIEle17_WPLoose_Gsf_v*",hltBits); // 5 TeV
  }
}
Bool_t isEleTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Bool_t isData, Bool_t is13TeV=true) {
  if (isData) {
    return triggerMenu.passObj("HLT_HIEle17_WPLoose_Gsf_v*","hltEle17WPLoose1GsfTrackIsoFilterForHI",hltMatchBits); // both 13 and 5
    // return triggerMenu.passObj("HLT_HIEle20_WPLoose_Gsf_v*","hltEle20WPLoose1GsfTrackIsoFilterForHI",hltMatchBits); // both 13 and 5
  }
  else {
    if(is13TeV){return triggerMenu.passObj("HLT_Ele27_WPTight_Gsf_v*","hltEle27WPTightGsfTrackIsoFilter",hltMatchBits);} // 13 TeV
    // else{ return triggerMenu.passObj("HLT_HIEle20_WPLoose_Gsf_v*","hltEle20WPLoose1GsfTrackIsoFilterForHI",hltMatchBits);} // 5 TeV
    else return triggerMenu.passObj("HLT_HIEle17_WPLoose_Gsf_v*","hltEle17WPLoose1GsfTrackIsoFilterForHI",hltMatchBits); // 5 TeV
    // return triggerMenu.passObj("HLT_Ele27_WPLoose_Gsf_v*","hltEle27noerWPLooseGsfTrackIsoFilter",hltMatchBits);
  }
}

//--------------------------------------------------------------------------------------------------
Double_t getEffAreaEl(const Double_t eta) { //2017
  if      (fabs(eta) < 1.0) return 0.1440;
  else if (fabs(eta) < 1.479)  return 0.1562;
  else if (fabs(eta) < 2.0)  return 0.1032;
  else if (fabs(eta) < 2.2)  return 0.0859;
  else if (fabs(eta) < 2.3)  return 0.1116;
  else if (fabs(eta) < 2.4)  return 0.1321;
  else if (fabs(eta) < 2.5)  return 0.1654;
  else return -1; // This should never happen, a cut on |eta|<2.5 is applied before this function is used.
}

Double_t getEffAreaEl_2015(const Double_t eta) { // 2015
  if      (fabs(eta) < 1.0) return 0.1752;
  else if (fabs(eta) < 1.479)  return 0.1862;
  else if (fabs(eta) < 2.0)  return 0.1411;
  else if (fabs(eta) < 2.2)  return 0.1534;
  else if (fabs(eta) < 2.3)  return 0.1903;
  else if (fabs(eta) < 2.4)  return 0.2243;
  else if (fabs(eta) < 2.5)  return 0.2687;
  else return -1; // This should never happen, a cut on |eta|<2.5 is applied before this function is used.
}

Double_t getEffAreaMu(const Double_t eta) {
  // not used, so i didn't update? 
  // well, ok this probably should be used....
  if      (fabs(eta) < 0.8) return 0.0913;
  else if (fabs(eta) < 1.3)  return 0.0765;
  else if (fabs(eta) < 2.0)  return 0.0546;
  else if (fabs(eta) < 2.2)  return 0.0728;
  else if (fabs(eta) < 2.5)  return 0.1177;
  else return -1; // This should never happen, a cut on |eta|<2.4 is applied in selection.
}

//---------------------------------------------------------------------------------------------------
Int_t getEtaBinLabel(const Double_t eta) {

  int ieta = -1;
  if(eta<-2.0)              return 0;
  else if(eta<-1.566)       return 1;
  else if(eta<-1.4442)      return 2;
  else if(eta<-1.0)         return 3;
  else if(eta<-0.5)         return 4;
  else if(eta< 0.0)         return 5;
  else if(eta< 0.5)         return 6;
  else if(eta< 1.0)         return 7;
  else if(eta< 1.4442)      return 8;
  else if(eta< 1.566)       return 9;
  else if(eta< 2.0)         return 10;
  else 			    return 11;
}
#endif

