#ifndef EWKANA_UTILS_LEPTONIDCUTS_HH
#define EWKANA_UTILS_LEPTONIDCUTS_HH

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include <TMath.h>
#include <cassert>
		   
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passAntiMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho=0);

Bool_t passEleID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passEleTightID(const baconhep::TElectron *electron, const Double_t rho);
Bool_t passEleLooseID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passAntiEleID(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho=0);

Bool_t isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits);
Bool_t isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1);

Bool_t isMuonTriggerNoIso(baconhep::TTrigger triggerMenu, TriggerBits hltBits);
Bool_t isMuonTriggerObjNoIso(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1);

Bool_t isEleTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Bool_t isData);
Bool_t isEleTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Bool_t isData);

Double_t getEffAreaEl(const Double_t eta);
Double_t getEffAreaMu(const Double_t eta);

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

//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Medium Electron ID for PU20 bx25

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;

  // conversion rejection
  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  // barrel/endcap dependent requirements
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
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
    // endcap
//    if(iso >= 0.0678*(electron->pt))                                  return kFALSE;
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

//--------------------------------------------------------------------------------------------------
Bool_t passAntiEleID(const baconhep::TElectron *electron,const TLorentzVector tag, const Double_t rho)
{
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;

  // conversion rejection
  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());	
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  // barrel/endcap dependent requirements
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
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
//    if(iso >= 0.0678*(electron->pt))                                  return kFALSE;
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
Bool_t passEleLooseID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Veto ID 
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  // conversion rejection
  if(electron->isConv) return kFALSE;
       
  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);     

  // barrel/endcap dependent requirements
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
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
Bool_t isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits) {
  return triggerMenu.pass("HLT_IsoMu20_v*",hltBits);
}

Bool_t isMuonTriggerNoIso(baconhep::TTrigger triggerMenu, TriggerBits hltBits) {
  return triggerMenu.pass("HLT_Mu20_v*",hltBits);
}

Bool_t isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1) {
  if (isL1) return triggerMenu.passObj("HLT_IsoMu20_v*","hltL1sL1SingleMu16",hltMatchBits);
  else return triggerMenu.passObj("HLT_IsoMu20_v*","hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09",hltMatchBits);
}

Bool_t isMuonTriggerObjNoIso(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1) {
  if (isL1) return triggerMenu.passObj("HLT_Mu20_v*","hltL1sL1SingleMu16",hltMatchBits);
  else return triggerMenu.passObj("HLT_Mu20_v*","hltL3fL1sMu16L1f0L2f10QL3Filtered20Q",hltMatchBits);
}


Bool_t isEleTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Bool_t isData) {
  if (isData) {
    return triggerMenu.pass("HLT_Ele23_WPLoose_Gsf_v*",hltBits);
  }
  else {
    return triggerMenu.pass("HLT_Ele23_WPLoose_Gsf_v*",hltBits);
  }
}

Bool_t isEleTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Bool_t isData) {
  if (isData) {
    if (isL1) {
      return triggerMenu.passObj("HLT_Ele23_WPLoose_Gsf_v*","hltEGL1SingleEG20ORL1SingleEG15Filter",hltMatchBits);
    }
    else {
      return triggerMenu.passObj("HLT_Ele23_WPLoose_Gsf_v*","hltEle23WPLooseGsfTrackIsoFilter",hltMatchBits);
    }
  }
  else if (isL1) {
    return triggerMenu.passObj("HLT_Ele23_WPLoose_Gsf_v*","hltL1sL1SingleEG20",hltMatchBits);
  }
  else {
    return triggerMenu.passObj("HLT_Ele23_WPLoose_Gsf_v*","hltEle23WPLooseGsfTrackIsoFilter",hltMatchBits);
  }
}

//--------------------------------------------------------------------------------------------------
Double_t getEffAreaEl(const Double_t eta) {
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
  else                      return 11;
}

#endif

