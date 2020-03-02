#ifndef Prefiring_Efficiency
#define Prefiring_Efficiency

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class

// The bacon particle interfaces
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"

// helper class to handle maps
#include "../Utils/CCorrUser2D.hh"
#include "../Utils/CCorrUser2D.hh"

// Make everything public for now while I test this
class PrefiringEfficiency {
  public:
    PrefiringEfficiency(){}
    ~PrefiringEfficiency(){}
    PrefiringEfficiency(TString fname, TString era){
      this->fname = fname;
      this->era = era;
      f = new TFile(fname);
      assert(f);
      loadJetMap();
      loadPhotonMap();
    }
    
    void loadJetMap(){
      cout << "Loading Prefiring Jet map for " << era.Data() << endl;
      mapJets.loadCorr((TH2D*)f->Get("L1prefiring_jetpt_"+era));
      return;
    }
    void loadPhotonMap(){
      cout << "Loading Prefiring Photon map for " << era.Data() << endl;
      mapPhot.loadCorr((TH2D*)f->Get("L1prefiring_photonpt_"+era));
      return;
    }
    
    void setObjects(TClonesArray *photons, TClonesArray *jets){
      scArr  = photons;
      jetArr = jets;
      return;
    }
    
    // could be done better, but i am lazy
    void computeJetsOnly(float &main, float &up, float &down){
      // loop through photons
      float uncTot = 0;
      for(Int_t ij=0; ij<jetArr->GetEntriesFast(); ij++) {
        const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[ij]);
        if(!etaCut(jet->eta)) continue;
        float eff = mapJets.getCorr(jet->eta, jet->pt);
        float unc = max(mapJets.getErr(jet->eta, jet->pt), (float)(eff*0.2));
        
        main   *= 1 - eff;
        uncTot += unc*unc;
      } 
      
      up = min(main+sqrt(uncTot),(float)1.0);
      down = max(main-sqrt(uncTot),(float)0.0);
      return;
    }
    
    void computePhotonsOnly(float &main, float &up, float &down){
      float uncTot = 0;
      for(Int_t ip=0; ip<scArr->GetEntriesFast(); ip++) {
        const baconhep::TPhoton *photon = (baconhep::TPhoton*)((*scArr)[ip]);
        if(!etaCut(photon->eta)) continue;
        
        float eff = mapPhot.getCorr(photon->eta, photon->pt);
        float unc = max(mapPhot.getErr(photon->eta, photon->pt), (float)(eff*0.2));
        
        main   *= 1 - eff;
        uncTot += unc*unc;
      } 
      
      up = min(main+sqrt(uncTot),(float)1.0);
      down = max(main-sqrt(uncTot),(float)0.0);
      return;
    }

    void computeFullPrefire(float &main, float &up, float &down){
      float unc = 0.;
      int njets = jetArr->GetEntriesFast();
      // in cases with no jets?
      if(njets == 0){
        cout << "No Jets - computing directly from photons " << endl;
        computePhotonsOnly(main,up,down);
        return;
      }
      // loop through jets!
      for(int ij=0; ij < njets; ij++) {
        const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[ij]);
        if(!etaCut(jet->eta)) continue;
        float effJet = mapJets.getCorr(jet->eta, jet->pt);
        float uncJet = max(mapJets.getErr( jet->eta, jet->pt), (float)(effJet*0.2) );
        float effObj = effJet;
        float uncObj = uncJet;
        if(effJet < 0 || effJet > 1) cout << "what?? " << endl;
        
        // check photons: 
        for(Int_t ip=0; ip<scArr->GetEntriesFast(); ip++) {
          const baconhep::TPhoton *photon = (baconhep::TPhoton*)((*scArr)[ip]);
          // only consider photons in correct eta region
          if(!etaCut(photon->eta)) continue;
          // check if the photons overlap with the jet
          if(toolbox::deltaR(jet->eta, jet->phi, photon->eta, photon->phi)>0.4) continue;
          // get the efficiency & uncertainty for this photon
          float effPho = mapPhot.getCorr(photon->eta, photon->pt);
          float uncPho = max(mapPhot.getErr( photon->eta, photon->pt), (float)(effPho*0.2) );
          
          // if we have an overlapping photon & jet, take the maximum of the jetcorr & photoncorr
          if(effPho > effJet){
            effObj = effPho;
            uncObj = uncPho;
          }
          
        } // end photon loop
        main *= 1 - effObj;
        unc += uncObj*uncObj; 
      } // end jet loop
      up   = min( main + sqrt(unc), (float)1.0 );
      down = max( main - sqrt(unc), (float)0.0 );
      // cout << "prefire eff: " << main << "  +1sig: " << up << "  -1sig: " << down << endl;
      return;
    }
    
  private:
    TString fname, era;
    CCorrUser2D mapPhot, mapJets;
    TClonesArray *scArr;
    TClonesArray *jetArr;
    TFile *f;
    
    bool etaCut(double eta){
      if(fabs(eta) < 2 || fabs(eta) > 3) return false;
      return true;
    }
  
};

#endif