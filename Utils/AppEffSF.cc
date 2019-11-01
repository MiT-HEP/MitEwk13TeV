#ifndef App_Eff_SF
#define App_Eff_SF

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
// #include "AppEffSF.hh"
// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"


struct basicEff {
  // data pos
  CEffUser2D dataPos;
  // data neg
  CEffUser2D dataNeg;
  // mc pos
  CEffUser2D mcPos;
  // mc neg
  CEffUser2D mcNeg;
  
};


struct effUnc {
  // fsr pos
  CEffUser2D fsrPos;
  // fsr neg
  CEffUser2D fsrNeg;
  // mc pos
  CEffUser2D mcPos;
  // mc neg
  CEffUser2D mcNeg;
  // bkg pos
  CEffUser2D bkgPos;
  // bkg neg
  CEffUser2D bkgNeg;
  // tag pt pos
  CEffUser2D tagPos;
  // tag pt neg
  CEffUser2D tagNeg;
  
};

// Make everything public for now while I test this
class AppEffSF {
  public:
    AppEffSF(){}
    ~AppEffSF(){}
    AppEffSF(TString dirName){
      this->dirName = dirName;
      isMuon=false;
    }
    // some functions to load the individual correction locations
    void loadHLT(TString dir, TString q1, TString q2){
      loadEfficiencyFile( hlt, dir,  q1,  q2);
      return;
    }
    
    void loadSta( TString dir, TString q1, TString q2){
      loadEfficiencyFile( sta, dir,  q1,  q2);
      isMuon = true;
      return;
    }
    
    void loadSel( TString dir, TString q1, TString q2){
      loadEfficiencyFile( sel, dir,  q1,  q2);
      return;
    }
    
    void loadUncSel(TString dir){
      cout << "Loading Selection Uncertainties " << endl;
      loadUncertaintyFile(u_sel, dir);
      return;
    }
    
    void loadUncSta(TString dir){
      loadUncertaintyFile(u_sta, dir);
      return;
    }
    
    // do i also need for hlt?
    vector<double> getUncSta(TLorentzVector *l1, int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      return computeAllUncertainties(sta, u_sta, l1, q1, l2, q2);
    }
    
    vector<double> getUncSel(TLorentzVector *l1, int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      return computeAllUncertainties(sel, u_sel, l1, q1, l2, q2);
    }
    
    // return a vector with each of the uncertainties propagated to final acc values
    vector<double> computeAllUncertainties(basicEff &effs, effUnc &uncs, TLorentzVector *l1, int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      vector<double> res;
      // cout << "FSR " << endl;
      res.push_back(computeWithUnc(effs, uncs.fsrPos, uncs.fsrNeg, l1, q1, l2, q2));
      
      // cout << "mc " << endl;
      res.push_back(computeWithUnc(effs, uncs.mcPos,  uncs.mcNeg,  l1, q1, l2, q2));
      // cout << "bkg " << endl;
      res.push_back(computeWithUnc(effs, uncs.bkgPos, uncs.bkgNeg, l1, q1, l2, q2));
      // cout << "tag " << endl;
      res.push_back(computeWithUnc(effs, uncs.tagPos, uncs.tagNeg, l1, q1, l2, q2));
      
      return res;
    }
    
    double computeWithUnc(basicEff &effs, CEffUser2D &varPos, CEffUser2D &varNeg, TLorentzVector *l1, int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      double effdata = 1, effmc = 1;
      // cout << "hello" << endl;
      if(q1>0) {
          effmc   *= effs.mcPos.getEff(l1->Eta(), l1->Pt());
          effdata *= effs.dataPos.getEff(l1->Eta(), l1->Pt()) * varPos.getEff(l1->Eta(), l1->Pt());
        } else {
          effmc   *= effs.mcNeg.getEff(l1->Eta(), l1->Pt());
          effdata *= effs.dataNeg.getEff(l1->Eta(), l1->Pt()) * varNeg.getEff(l1->Eta(), l1->Pt());
        }
        if(!l2) {
          return  effdata/effmc;
        }
        if(q2>0) {
          effmc   *= effs.mcPos.getEff(l2->Eta(), l2->Pt());
          effdata *= effs.dataPos.getEff(l2->Eta(), l2->Pt()) * varPos.getEff(l2->Eta(), l2->Pt());
        } else {
          effmc   *= effs.mcNeg.getEff(l2->Eta(), l2->Pt());
          effdata *= effs.dataNeg.getEff(l2->Eta(), l2->Pt()) * varNeg.getEff(l2->Eta(), l2->Pt());
        }
      return effdata/effmc;
    }
    
    // some functions to load the systmatic uncertainty files
    // syst unc only applies to the GSF (electrons), and Sta / SIT (muons)
    void addSystematic(){
      return;
    }
    
    // compute efficiencies
    double fullEfficiencies(TLorentzVector *l1, int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      double effmc = 1, effdata = 1;
      // vector<double>
      // if there is a second lepton we need to calculate twice
      double corr = 1;
      corr *= computeHLTSF(l1, q1, l2, q2);
      corr *= computeSelSF(l1, q1, l2, q2);
      if(isMuon){
        // cout << "Muon" << endl;
        corr *= computeStaSF(l1, q1, l2, q2);
      }
     
      
      return corr;
    }
  

    
    
    double statUncHLT(TLorentzVector *l1,  int q1, TH2D *pos, TH2D *neg, double wgt){
      return statUncTrigger(hlt, l1, q1, pos, neg, wgt);
    }
    double statUncSel(TLorentzVector *l1,  int q1, TH2D *pos, TH2D *neg, double wgt){
      return statUnc(sel, l1, q1, pos, neg, wgt);
    }
    double statUncSta(TLorentzVector *l1,  int q1, TH2D *pos, TH2D *neg, double wgt){
      return statUnc(sta, l1, q1, pos, neg, wgt);
    }

    
  //private:
  
    // structs i'm storying the corrections in
    basicEff hlt;
    basicEff sta;
    basicEff sel;
    
    effUnc u_hlt;
    effUnc u_sta;
    effUnc u_sel;
    
    TString dirName; // name of primary directory
    bool isMuon;
    
    void loadEfficiencyFile(basicEff &effs, TString dir, TString q1, TString q2){
      
      TString dataPosName = dirName + "Data/" + dir +"/"+ q1 + "/eff.root"; 
      TString dataNegName = dirName + "Data/" + dir +"/"+ q2 + "/eff.root"; 
      TString mcPosName   = dirName + "MC/"   + dir +"/"+ q1 + "/eff.root"; 
      TString mcNegName   = dirName + "MC/"   + dir +"/"+ q2 + "/eff.root"; 
      cout << "attempting to load " << dataPosName << endl;
      
      TFile *dp = new TFile(dataPosName);
      effs.dataPos.loadEff((TH2D*)dp->Get("hEffEtaPt"), (TH2D*)dp->Get("hErrlEtaPt"), (TH2D*)dp->Get("hErrhEtaPt"));
      
      TFile *dn = new TFile(dataNegName);
      effs.dataNeg.loadEff((TH2D*)dn->Get("hEffEtaPt"), (TH2D*)dn->Get("hErrlEtaPt"), (TH2D*)dn->Get("hErrhEtaPt"));
      
      TFile *mp = new TFile(mcPosName);
      effs.mcPos.loadEff((TH2D*)mp->Get("hEffEtaPt"), (TH2D*)mp->Get("hErrlEtaPt"), (TH2D*)mp->Get("hErrhEtaPt"));
      
      TFile *mn = new TFile(mcNegName);
      effs.mcNeg.loadEff((TH2D*)mn->Get("hEffEtaPt"), (TH2D*)mn->Get("hErrlEtaPt"), (TH2D*)mn->Get("hErrhEtaPt"));
      
      return;
    }
    
    void loadUncertaintyFile(effUnc &uncs, TString filename){
      TFile *fSys = new TFile(filename);
      
      cout << "Loading ... " << filename.Data() << endl;
      
      uncs.fsrNeg.loadEff((TH2D*)fSys->Get("hFSRNeg"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      uncs.fsrPos.loadEff((TH2D*)fSys->Get("hFSRPos"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      
      uncs.mcNeg.loadEff((TH2D*)fSys->Get("hMCNeg"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      uncs.mcPos.loadEff((TH2D*)fSys->Get("hMCPos"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      
      uncs.bkgNeg.loadEff((TH2D*)fSys->Get("hBkgNeg"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      uncs.bkgPos.loadEff((TH2D*)fSys->Get("hBkgPos"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      
      uncs.tagNeg.loadEff((TH2D*)fSys->Get("hTagNeg"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      uncs.tagPos.loadEff((TH2D*)fSys->Get("hTagPos"), (TH2D*)fSys->Get(""), (TH2D*)fSys->Get(""));
      
      // cout << uncs.fsrNeg.getEff(0, 30) << endl; 
      
      return;
    }
    
       
    double statUnc(basicEff eff, TLorentzVector *l1,  int q1, TH2D *pos, TH2D *neg, double wgt){
      double var = 0.0;
      if(q1>0) {
        Double_t effdata = eff.dataPos.getEff(l1->Eta(), l1->Pt());
        Double_t errdata = TMath::Max(eff.dataPos.getErrLow(l1->Eta(), l1->Pt()), eff.dataPos.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t effmc   = eff.mcPos.getEff(l1->Eta(), l1->Pt());
        Double_t errmc   = TMath::Max(eff.mcPos.getErrLow(l1->Eta(), l1->Pt()), eff.mcPos.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        pos->Fill(l1->Eta(), l1->Pt(), errSta*wgt);
        var+=errSta*errSta;
      } else {
        Double_t effdata = eff.dataNeg.getEff(l1->Eta(), l1->Pt());
        Double_t errdata = TMath::Max(eff.dataNeg.getErrLow(l1->Eta(), l1->Pt()), eff.dataNeg.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t effmc   = eff.mcNeg.getEff(l1->Eta(), l1->Pt());
        Double_t errmc   = TMath::Max(eff.mcNeg.getErrLow(l1->Eta(), l1->Pt()), eff.mcNeg.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        neg->Fill(l1->Eta(), l1->Pt(), errSta*wgt);
        var+=errSta*errSta;
      }
      return var;
    }
    
    double statUncTrigger(basicEff eff, TLorentzVector *l1,  int q1, TH2D *pos, TH2D *neg, double wgt){
      double var = 0.0;
      if(q1>0) {
        Double_t effdata = 1-eff.dataPos.getEff(l1->Eta(), l1->Pt());
        Double_t errdata = TMath::Max(eff.dataPos.getErrLow(l1->Eta(), l1->Pt()), eff.dataPos.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t effmc   = 1-eff.mcPos.getEff(l1->Eta(), l1->Pt());
        Double_t errmc   = TMath::Max(eff.mcPos.getErrLow(l1->Eta(), l1->Pt()), eff.mcPos.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t errSta = ((1 - effdata)/(1 - effmc))*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        pos->Fill(l1->Eta(), l1->Pt(), errSta*wgt);
        var+=errSta*errSta;
      } else {
        Double_t effdata = 1-eff.dataNeg.getEff(l1->Eta(), l1->Pt());
        Double_t errdata = TMath::Max(eff.dataNeg.getErrLow(l1->Eta(), l1->Pt()), eff.dataNeg.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t effmc   = 1-eff.mcNeg.getEff(l1->Eta(), l1->Pt());
        Double_t errmc   = TMath::Max(eff.mcNeg.getErrLow(l1->Eta(), l1->Pt()), eff.mcNeg.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t errSta = ((1 - effdata)/(1 - effmc))*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        neg->Fill(l1->Eta(), l1->Pt(), errSta*wgt);
        var+=errSta*errSta;
      }
      return var;
    }
    
    double computeHLTSF(TLorentzVector *l1,  int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      double effdata = 1.0; 
      double effmc = 1.0;
      double corr = 1.0;
      if(q1>0) { 
        effdata *= (1.-hlt.dataPos.getEff((l1->Eta()), l1->Pt())); 
        effmc   *= (1.-hlt.mcPos.getEff((l1->Eta()), l1->Pt())); 
      } else {
        effdata *= (1.-hlt.dataNeg.getEff((l1->Eta()), l1->Pt())); 
        effmc   *= (1.-hlt.mcNeg.getEff((l1->Eta()), l1->Pt())); 
      }
      if(!l2) {
        effdata = 1.-effdata;
        effmc   = 1.-effmc;
        corr *= effdata/effmc;
        return corr;
      }
      if(q2>0) {
        effdata *= (1.-hlt.dataPos.getEff((l2->Eta()), l2->Pt())); 
        effmc   *= (1.-hlt.mcPos.getEff((l2->Eta()), l2->Pt()));
      } else {
        effdata *= (1.-hlt.dataNeg.getEff((l2->Eta()), l2->Pt())); 
        effmc   *= (1.-hlt.mcNeg.getEff((l2->Eta()), l2->Pt()));
      }
      effdata = 1.-effdata;
      effmc   = 1.-effmc;
      corr *= effdata/effmc;
      return corr;
    }
    
    double computeStaSF(TLorentzVector *l1,  int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      double effdata = 1.0; 
      double effmc = 1.0;
      double corr = 1.0;
      if(q1>0) { 
        effdata *= sta.dataPos.getEff((l1->Eta()), l1->Pt()); 
        effmc   *= sta.mcPos.getEff((l1->Eta()), l1->Pt()); 
      } else {
        effdata *= sta.dataNeg.getEff((l1->Eta()), l1->Pt()); 
        effmc   *= sta.mcNeg.getEff((l1->Eta()), l1->Pt()); 
      }
      if(!l2) {
        effdata = effdata;
        effmc   = effmc;
        corr *= effdata/effmc;
        return corr;
      }
      if(q2>0) {
        effdata *= sta.dataPos.getEff((l2->Eta()), l2->Pt()); 
        effmc   *= sta.mcPos.getEff((l2->Eta()), l2->Pt());
      } else {
        effdata *= sta.dataNeg.getEff((l2->Eta()), l2->Pt()); 
        effmc   *= sta.mcNeg.getEff((l2->Eta()), l2->Pt());
      }
      effdata = effdata;
      effmc   = effmc;
      corr *= effdata/effmc;
      return corr;
    }
    
    double computeSelSF(TLorentzVector *l1,  int q1, TLorentzVector* l2 = nullptr, int q2 = 0){
      double effdata = 1.0; 
      double effmc = 1.0;
      double corr = 1.0;
      if(q1>0) { 
        effdata *= sel.dataPos.getEff((l1->Eta()), l1->Pt()); 
        effmc   *= sel.mcPos.getEff((l1->Eta()), l1->Pt()); 
      } else {
        effdata *= sel.dataNeg.getEff((l1->Eta()), l1->Pt()); 
        effmc   *= sel.mcNeg.getEff((l1->Eta()), l1->Pt()); 
      }
      if(!l2) {
        effdata = effdata;
        effmc   = effmc;
        corr *= effdata/effmc;
        return corr;
      }
      if(q2>0) {
        effdata *= sel.dataPos.getEff((l2->Eta()), l2->Pt()); 
        effmc   *= sel.mcPos.getEff((l2->Eta()), l2->Pt());
      } else {
        effdata *= sel.dataNeg.getEff((l2->Eta()), l2->Pt()); 
        effmc   *= sel.mcNeg.getEff((l2->Eta()), l2->Pt());
      }
      effdata = effdata;
      effmc   = effmc;
      corr *= effdata/effmc;
      return corr;
    }
    
  
};

#endif