//================================================================================================
//
// Compute W->munu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties need to be checked//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TLorentzVector.h"

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions

// helper class to handle efficiency tables
#include "../Utils/CEffUser2D.hh"
#include "../Utils/AppEffSF.cc"
#endif


//=== MAIN MACRO ================================================================================================= 

void computeAccSelWm_Sys(const TString conf,       // input file
          const TString inputDir,
          const TString outputDir,  // output directory
		     const Int_t   charge,      // 0 = inclusive, +1 = W+, -1 = W-
		     const Int_t   doPU,
			    const TString sysFileSIT, // condense these into 1 file per type of eff (pos & neg into 1 file)
			    const TString sysFileSta,
          const bool is13TeV=1
) {
  gBenchmark->Start("computeAccSelWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  const Double_t mu_MASS  = 0.1057;

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  // const Double_t ETA_BARREL = 1.2;
  // const Double_t ETA_ENDCAP = 1.2;

  const Double_t ETA_BARREL = 1.442;
  const Double_t ETA_ENDCAP = 1.442;
  const Double_t VETO_PT   = 10;
  const Double_t VETO_ETA  = 2.4;

  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 13;
  
  const int NBptSta  = 3;
  const float ptrangeSta[NBptSta +1]   = {25., 35,  50., 10000.};

  const int NBeta = 12;
  const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};
  const int NBptHLT = 12;
  const float ptrangeHLT[NBptHLT+1] = {25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 10000};
  
  AppEffSF effs(inputDir);
  effs.loadHLT("MuHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("MuSITEff_aMCxPythia","Combined","Combined");
  effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  effs.loadUncSel(sysFileSIT);
  effs.loadUncSta(sysFileSta);
  
  TH2D *hSelErr_pos = new TH2D("hSelErr_pos", "",NBeta,etarange,NBptSta,ptrangeSta);
  TH2D *hSelErr_neg = new TH2D("hSelErr_neg", "",NBeta,etarange,NBptSta,ptrangeSta);
  
  TH2D *hStaErr_pos = new TH2D("hStaErr_pos", "",NBeta,etarange,NBptSta,ptrangeSta);
  TH2D *hStaErr_neg = new TH2D("hStaErr_neg", "",NBeta,etarange,NBptSta,ptrangeSta);

  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",NBeta,etarange,NBptHLT,ptrangeHLT);
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",NBeta,etarange,NBptHLT,ptrangeHLT);

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/pileup_rw_76X.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("h_rw_golden");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   colorv;  // plot color per input file
  vector<Int_t>   linev;   // plot line style per input file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    ss >> fname >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
  }
  ifs.close();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  // Data structures to store info from TTrees
  baconhep::TEventInfo   *info = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray     *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray        *muonArr = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;

  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBv, nSelEv;
  vector<Double_t> accv, accBv, accEv;
  vector<Double_t> accErrv, accErrBv, accErrEv;
  vector<Double_t> nSelCorrv, nSelBCorrv, nSelECorrv;
  vector<Double_t> nSelCorrVarv, nSelBCorrVarv, nSelECorrVarv;
  vector<Double_t> accCorrv, accBCorrv, accECorrv;
  vector<Double_t> accErrCorrv, accErrBCorrv, accErrECorrv;
  
  vector<Double_t> nSelCorrvFSR, nSelCorrvMC, nSelCorrvBkg, nSelCorrvTag;//, nSelCorrvStat;
  vector<Double_t> nSelCorrvFSR_I, nSelCorrvMC_I, nSelCorrvBkg_I, nSelCorrvTag_I;//, nSelCorrvStat_I;
  vector<Double_t> nSelCorrvFSR_S, nSelCorrvMC_S, nSelCorrvBkg_S, nSelCorrvTag_S;//, nSelCorrvStat_S;
  vector<Double_t> nSelCorrVarvFSR, nSelCorrVarvMC, nSelCorrVarvBkg, nSelCorrVarvTag;//, nSelCorrVarvStat;
  vector<Double_t> pctDiffvFSR, pctDiffvMC, pctDiffvBkg, pctDiffvTag;//, accCorrvStat;
  vector<Double_t> accCorrvFSR, accCorrvMC, accCorrvBkg, accCorrvTag;//, accCorrvStat;
  vector<Double_t> accCorrvFSR_I, accCorrvMC_I, accCorrvBkg_I, accCorrvTag_I;//, accCorrvStat_I;
  vector<Double_t> accCorrvFSR_S, accCorrvMC_S, accCorrvBkg_S, accCorrvTag_S;//, accCorrvStat_S;
  vector<Double_t> accErrCorrvFSR, accErrCorrvMC, accErrCorrvBkg, accErrCorrvTag;//, accErrCorrvStat;

  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr    = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr     = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch *genPartBr = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Muon",          &muonArr); TBranch *muonBr    = eventTree->GetBranch("Muon"); 
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBv.push_back(0);
    nSelEv.push_back(0);
    nSelCorrv.push_back(0);
    nSelBCorrv.push_back(0);
    nSelECorrv.push_back(0);
    nSelCorrVarv.push_back(0);
    nSelBCorrVarv.push_back(0);
    nSelECorrVarv.push_back(0);
    nSelCorrvFSR.push_back(0);  nSelCorrVarvFSR.push_back(0);
    nSelCorrvMC.push_back(0);   nSelCorrVarvMC.push_back(0);
    nSelCorrvBkg.push_back(0);  nSelCorrVarvBkg.push_back(0);
    nSelCorrvTag.push_back(0);  nSelCorrVarvTag.push_back(0);
    
    nSelCorrvFSR_I.push_back(0);  nSelCorrvFSR_S.push_back(0);
    nSelCorrvMC_I.push_back(0);   nSelCorrvMC_S.push_back(0);
    nSelCorrvBkg_I.push_back(0);  nSelCorrvBkg_S.push_back(0);
    nSelCorrvTag_I.push_back(0);  nSelCorrvTag_S.push_back(0);
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<(0.25)*((uint)eventTree->GetEntries()); ientry++) {
      if(ientry%100000==0)   cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
    //for(UInt_t ientry=0; ientry<1000000; ientry++) {
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);

      if (charge==-1 && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
      if (charge==1 && toolbox::flavor(genPartArr, BOSON_ID)!=-LEPTON_ID) continue;
      if (charge==0 && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      
      
      
      /*TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,1);*/
      Int_t glepq1=-99;
      Int_t glepq2=-99;
      TLorentzVector *gvec=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *glep2=new TLorentzVector(0,0,0,0);
	    toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2,&glepq1,&glepq2,1);
   
      // TLorentzVector tvec=*glep1+*glep2;
      // TLorentzVector* genV=new TLorentzVector(0,0,0,0);
      // genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
      // genVPt   = tvec.Pt();
      // genVPhi  = tvec.Phi();
      // genVy    = tvec.Rapidity();
      // double genVMass = tvec.M();
      double mtgen  = sqrt( 2.0 * (glep1->Pt()) * (glep2->Pt()) * (1.0-cos(toolbox::deltaPhi(glep1->Phi(),glep2->Phi()))) );
      // if(mtgen > 40) continue;
      // if(mtgen < 40 || mtgen > 140) continue;
     // cout << "mass " << genVMass <<  "   mt " << mt << endl;
      
      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;
      
      TLorentzVector vMu(0,0,0,0);
      // trigger requirement               
      // if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
      if (!isMuonTrigger(triggerMenu, info->triggerBits,kFALSE,is13TeV)) continue;
   
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);
      Int_t nLooseLep=0;
      const baconhep::TMuon *goodMuon=0;
      Bool_t passSel=kFALSE;
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
        const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);

        if(fabs(mu->eta) > VETO_ETA) continue; // loose lepton |eta| cut
        if(mu->pt	 < VETO_PT)  continue; // loose lepton pT cut
        if(passMuonLooseID(mu)) nLooseLep++;   // loose lepton selection
        if(nLooseLep>1) {  // extra lepton veto
          passSel=kFALSE;
          break;
        }
        
        if(fabs(mu->eta) > ETA_CUT)         continue;  // lepton |eta| cut
        if(mu->pt < PT_CUT)		    continue;  // lepton pT cut	
        if(!passMuonID(mu))		    continue;  // lepton selection
        if(!isMuonTriggerObj(triggerMenu, mu->hltMatchBits, kFALSE,is13TeV)) continue;
	
        if(charge!=0 && mu->q!=charge) continue;  // check charge (if necessary)
	
        double mtreco  = sqrt( 2.0 * (mu->pt) * (info->pfMETC) * (1.0-cos(toolbox::deltaPhi(mu->phi,info->pfMETCphi))) );
        
        if(mtreco < 40) continue;
        // if(mtreco > 40 && mtreco < 140) continue;
        passSel=kTRUE;
        goodMuon=mu;
        vMu.SetPtEtaPhiM(mu->pt, mu->eta, mu->phi, mu_MASS);
      }
      
      if(passSel) {
        
        /******** We have a W candidate! HURRAY! ********/
              
        Bool_t isBarrel = (fabs(goodMuon->eta)<ETA_BARREL) ? kTRUE : kFALSE;
              
        // data/MC scale factor corrections
        Double_t effdata, effmc, emTag;
        Double_t edFSR, edMC, edBkg, edTag;//, edStat;
	      Double_t corr=1;
	      Double_t corrFSR=1, corrMC=1, corrBkg=1, corrTag=1;//, corrStat=1;
	      Double_t corrFSR_I=1, corrMC_I=1, corrBkg_I=1, corrTag_I=1;//, corrStat_I=1;
	      Double_t corrFSR_S=1, corrMC_S=1, corrBkg_S=1, corrTag_S=1;//, corrStat_S=1;
        
        effdata=1; effmc=1;   emTag=1; 
	      edFSR=1; edMC=1; edBkg=1; edTag=1;//edStat=1;  
	
        int q = goodMuon->q;
        
        corr = effs.fullEfficiencies(&vMu,q);
        vector<double> uncs_sta = effs.getUncSta(&vMu,q);
        vector<double> uncs_sit = effs.getUncSel(&vMu,q);
        
        corrFSR *= uncs_sta[0]*uncs_sit[0]*effs.computeHLTSF(&vMu,q); // alternate fsr model
        corrMC  *= uncs_sta[1]*uncs_sit[1]*effs.computeHLTSF(&vMu,q); // alternate mc gen model
        corrBkg *= uncs_sta[2]*uncs_sit[2]*effs.computeHLTSF(&vMu,q); // alternate bkg model
        corrTag *= uncs_sta[3]*uncs_sit[3]*effs.computeHLTSF(&vMu,q); // alternate bkg model
        // corr *= effdata/effmc; // orig
        
        corrFSR_I *= uncs_sit[0]*effs.computeHLTSF(&vMu,q)*effs.computeStaSF(&vMu,q); // alternate fsr model
        corrMC_I  *= uncs_sit[1]*effs.computeHLTSF(&vMu,q)*effs.computeStaSF(&vMu,q); // alternate mc gen model
        corrBkg_I *= uncs_sit[2]*effs.computeHLTSF(&vMu,q)*effs.computeStaSF(&vMu,q); // alternate bkg model
        corrTag_I *= uncs_sit[3]*effs.computeHLTSF(&vMu,q)*effs.computeStaSF(&vMu,q); // alternate bkg model
         
        corrFSR_S *= uncs_sta[0]*effs.computeHLTSF(&vMu,q)*effs.computeSelSF(&vMu,q); // alternate fsr model
        corrMC_S  *= uncs_sta[1]*effs.computeHLTSF(&vMu,q)*effs.computeSelSF(&vMu,q); // alternate mc gen model
        corrBkg_S *= uncs_sta[2]*effs.computeHLTSF(&vMu,q)*effs.computeSelSF(&vMu,q); // alternate bkg model
        corrTag_S *= uncs_sta[3]*effs.computeHLTSF(&vMu,q)*effs.computeSelSF(&vMu,q); // alternate bkg model
         
        double var=0.;        
        // var += effs.statUncSta(&l1, q) + effs.statUncSta(&l2, q2);
        var += effs.statUncSta(&vMu, q, hStaErr_pos, hStaErr_neg, fabs(weight)*corr);
        var += effs.statUncSel(&vMu, q, hSelErr_pos, hSelErr_neg, fabs(weight)*corr);
        var += effs.statUncHLT(&vMu, q, hHLTErr_pos, hHLTErr_neg, fabs(weight)*corr);
  
        nSelv[ifile]    +=weight;
        nSelCorrvFSR[ifile] +=weight*corrFSR;  nSelCorrvFSR_I[ifile] +=weight*corrFSR_I;  nSelCorrvFSR_S[ifile] +=weight*corrFSR_S;
        nSelCorrvMC[ifile]  +=weight*corrMC;   nSelCorrvMC_I[ifile]  +=weight*corrMC_I;   nSelCorrvMC_S[ifile]  +=weight*corrMC_S;
        nSelCorrvBkg[ifile] +=weight*corrBkg;  nSelCorrvBkg_I[ifile] +=weight*corrBkg_I;  nSelCorrvBkg_S[ifile] +=weight*corrBkg_S;
        nSelCorrvTag[ifile] +=weight*corrTag;  nSelCorrvTag_I[ifile] +=weight*corrTag_I;  nSelCorrvTag_S[ifile] +=weight*corrTag_S;
        // /*nSelCorrvStat[ifile]+=weight*corrStat; */nSelCorrvStat_I[ifile]+=weight*corrStat_I; nSelCorrvStat_S[ifile]+=weight*corrStat_S;
        
        
        nSelCorrv[ifile]+=weight*corr;
        nSelCorrVarvFSR[ifile]+=weight*weight*corrFSR*corrFSR;
        nSelCorrVarvMC[ifile]+=weight*weight*corrMC*corrMC;
        nSelCorrVarvBkg[ifile]+=weight*weight*corrBkg*corrBkg;
        nSelCorrVarvTag[ifile]+=weight*weight*corrTag*corrTag;
        nSelCorrVarv[ifile]+=weight*weight*corr*corr;
        
        if(isBarrel) { 
        nSelBv[ifile]+=weight;
        nSelBCorrv[ifile]+=weight*corr;
        nSelBCorrVarv[ifile]+=weight*weight*corr*corr;
          
        } else { 
          nSelEv[ifile]+=weight;
          nSelECorrv[ifile]+=weight*corr;
          nSelECorrVarv[ifile]+=weight*weight*corr*corr;
        }
      }
    }
    
    Double_t var=0, varB=0, varE=0;
    for(Int_t iy=0; iy<=hHLTErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hHLTErr_pos->GetNbinsX(); ix++) {
        Double_t err=hHLTErr_pos->GetBinContent(ix,iy);
        var+=err*err;
        err=hHLTErr_neg->GetBinContent(ix,iy);
        var+=err*err;
        // std::cout << "hlt pos " << var << std::endl;
      }
    }

    for(Int_t iy=0; iy<=hSelErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hSelErr_pos->GetNbinsX(); ix++) {
        Double_t err=hSelErr_pos->GetBinContent(ix,iy);
        var+=err*err;
        err=hSelErr_neg->GetBinContent(ix,iy);
        var+=err*err;
        // std::cout << "sel pos " << var << std::endl;
      }
    }

    for(Int_t iy=0; iy<=hStaErr_pos->GetNbinsY(); iy++) {
      for(Int_t ix=0; ix<=hStaErr_pos->GetNbinsX(); ix++) {
        Double_t err=hStaErr_pos->GetBinContent(ix,iy);
	      var+=err*err;
        err=hStaErr_neg->GetBinContent(ix,iy);
	      var+=err*err;
      }
    }
    nSelCorrVarv[ifile]+=var;
    nSelBCorrVarv[ifile]+=varB;
    nSelECorrVarv[ifile]+=varE;
    nSelCorrVarvFSR[ifile]+=var;
    nSelCorrVarvMC[ifile]+=var;
    nSelCorrVarvBkg[ifile]+=var;
    nSelCorrVarvTag[ifile]+=var;
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.+accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.+accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.+accEv[ifile])/nEvtsv[ifile]));
    
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]);   accErrCorrv.push_back(accCorrv[ifile]*sqrt(nSelCorrVarv[ifile]/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));
    accBCorrv.push_back(nSelBCorrv[ifile]/nEvtsv[ifile]); accErrBCorrv.push_back(accBCorrv[ifile]*sqrt(nSelBCorrVarv[ifile]/nSelBCorrv[ifile]/nSelBCorrv[ifile] + 1./nEvtsv[ifile]));
    accECorrv.push_back(nSelECorrv[ifile]/nEvtsv[ifile]); accErrECorrv.push_back(accECorrv[ifile]*sqrt(nSelECorrVarv[ifile]/(nSelECorrv[ifile]*nSelECorrv[ifile]) + 1./nEvtsv[ifile]));
    
    accCorrvFSR.push_back(nSelCorrvFSR[ifile]/nEvtsv[ifile]);
    accCorrvMC.push_back(nSelCorrvMC[ifile]/nEvtsv[ifile]);
    accCorrvBkg.push_back(nSelCorrvBkg[ifile]/nEvtsv[ifile]);
    accCorrvTag.push_back(nSelCorrvTag[ifile]/nEvtsv[ifile]);
    
    accCorrvFSR_I.push_back(nSelCorrvFSR_I[ifile]/nEvtsv[ifile]);
    accCorrvMC_I.push_back(nSelCorrvMC_I[ifile]/nEvtsv[ifile]);
    accCorrvBkg_I.push_back(nSelCorrvBkg_I[ifile]/nEvtsv[ifile]);
    accCorrvTag_I.push_back(nSelCorrvTag_I[ifile]/nEvtsv[ifile]);
    
    accCorrvFSR_S.push_back(nSelCorrvFSR_S[ifile]/nEvtsv[ifile]);
    accCorrvMC_S.push_back(nSelCorrvMC_S[ifile]/nEvtsv[ifile]);
    accCorrvBkg_S.push_back(nSelCorrvBkg_S[ifile]/nEvtsv[ifile]);
    accCorrvTag_S.push_back(nSelCorrvTag_S[ifile]/nEvtsv[ifile]);

    accErrCorrvFSR.push_back(accCorrvFSR[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvFSR[ifile]*nSelCorrvFSR[ifile]) + 1./nEvtsv[ifile]));
    accErrCorrvMC.push_back(accCorrvMC[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvMC[ifile]*nSelCorrvMC[ifile]) + 1./nEvtsv[ifile]));
    accErrCorrvBkg.push_back(accCorrvBkg[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvBkg[ifile]*nSelCorrvBkg[ifile]) + 1./nEvtsv[ifile]));
    accErrCorrvTag.push_back(accCorrvTag[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrvTag[ifile]*nSelCorrvTag[ifile]) + 1./nEvtsv[ifile]));
    
    accCorrv.push_back(nSelCorrv[ifile]/nEvtsv[ifile]); 
    accErrCorrv.push_back(accCorrv[ifile]*sqrt((nSelCorrVarv[ifile])/(nSelCorrv[ifile]*nSelCorrv[ifile]) + 1./nEvtsv[ifile]));
    
   
    delete infile;
    infile=0, eventTree=0;  
  }
  delete info;
  delete gen;
  delete muonArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  if(charge== 0) cout << " W -> mu nu"  << endl;
  if(charge==-1) cout << " W- -> mu nu" << endl;
  if(charge== 1) cout << " W+ -> mu nu" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    cout << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    cout << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    cout << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    cout << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    cout << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    cout << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    cout << "          pct: " << 100*accErrCorrv[ifile] /accCorrv[ifile] << endl;
    cout << "          FSR unc: " << accCorrvFSR[ifile]  << " / Sel: " << accCorrvFSR_I[ifile] << " / Sta: " << accCorrvFSR_S[ifile] << endl;
    cout << "           MC unc: " << accCorrvMC[ifile]   << " / Sel: " << accCorrvMC_I[ifile]  << " / Sta: " << accCorrvMC_S[ifile]  << endl;
    cout << "          Bkg unc: " << accCorrvBkg[ifile]  << " / Sel: " << accCorrvBkg_I[ifile] << " / Sta: " << accCorrvBkg_S[ifile] << endl;
    cout << "          Tag unc: " << accCorrvTag[ifile]  << " / Sel: " << accCorrvTag_I[ifile] << " / Sta: " << accCorrvTag_S[ifile] << endl;
    cout << "Acc (FSR/MC/Bkg/Tag): " << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[500];
  sprintf(txtfname,"%s/sel.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  if(charge== 0) txtfile << " W -> mu nu"  << endl;
  if(charge==-1) txtfile << " W- -> mu nu" << endl;
  if(charge== 1) txtfile << " W+ -> mu nu" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "     barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile];
    txtfile << "  ==eff corr==> " << accBCorrv[ifile] << " +/- " << accErrBCorrv[ifile] << endl;
    txtfile << "     endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile];
    txtfile << "  ==eff corr==> " << accECorrv[ifile] << " +/- " << accErrECorrv[ifile] << endl;
    txtfile << "      total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile];
    txtfile << "  ==eff corr==> " << accCorrv[ifile]  << " +/- " << accErrCorrv[ifile] << endl;
    txtfile << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    txtfile << "          FSR unc: " << accCorrvFSR[ifile]  << " / Sel: " << accCorrvFSR_I[ifile] << " / Sta: " << accCorrvFSR_S[ifile] << endl;
    txtfile << "           MC unc: " << accCorrvMC[ifile]   << " / Sel: " << accCorrvMC_I[ifile]  << " / Sta: " << accCorrvMC_S[ifile]  << endl;
    txtfile << "          Bkg unc: " << accCorrvBkg[ifile]  << " / Sel: " << accCorrvBkg_I[ifile] << " / Sta: " << accCorrvBkg_S[ifile] << endl;
    txtfile << "          Tag unc: " << accCorrvTag[ifile]  << " / Sel: " << accCorrvTag_I[ifile] << " / Sta: " << accCorrvTag_S[ifile] << endl;
    txtfile << "Acc (FSR/MC/Bkg/Tag): " << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

    txtfile << endl;
  }
  txtfile.close();  
  
  
  // char txtfname[100];
  sprintf(txtfname,"%s/sel_nums_only.txt",outputDir.Data());
  ofstream txtfile2;
  txtfile2.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {//accv[ifile]
    txtfile2 << "uncorrected: " << accv[ifile]    << endl;
    txtfile2 << accCorrv[ifile]     << " " << accErrCorrv[ifile]    << endl;
    txtfile2 << accCorrvFSR[ifile]  << " " << accCorrvFSR_I[ifile] << " " << accCorrvFSR_S[ifile] << endl;
    txtfile2 << accCorrvMC[ifile]   << " " << accCorrvMC_I[ifile]  << " " << accCorrvMC_S[ifile]  << endl;
    txtfile2 << accCorrvBkg[ifile]  << " " << accCorrvBkg_I[ifile] << " " << accCorrvBkg_S[ifile] << endl;
    txtfile2 << accCorrvTag[ifile]  << " " << accCorrvTag_I[ifile] << " " << accCorrvTag_S[ifile] << endl;
    // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

    txtfile2 << endl;
  }
  txtfile2.close();  
  
    // char txtfname[100];
  sprintf(txtfname,"%s/sit_unc.txt",outputDir.Data());
  ofstream txtfile3;
  txtfile3.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile3 << accCorrv[ifile]      << endl;
    txtfile3 << accCorrvFSR_I[ifile] << endl;
    txtfile3 << accCorrvMC_I[ifile]  << endl;
    txtfile3 << accCorrvBkg_I[ifile] << endl;
    txtfile3 << accCorrvTag_I[ifile] << endl;
    // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

    // txtfile3 << endl;
  }
  txtfile3.close();
  
      // char txtfname[100];
  sprintf(txtfname,"%s/sta_unc.txt",outputDir.Data());
  ofstream txtfile4;
  txtfile4.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile4 << accCorrv[ifile]   << endl;
    txtfile4 << accCorrvFSR_S[ifile] << endl;
    txtfile4 << accCorrvMC_S[ifile]  << endl;
    txtfile4 << accCorrvBkg_S[ifile] << endl;
    txtfile4 << accCorrvTag_S[ifile] << endl;
    // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

    // txtfile4 << endl;
  }
  txtfile4.close();
  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelWm"); 
}
