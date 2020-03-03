//================================================================================================
//
// Compute Z->mumu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties no yet implemented
//
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
// #include <math.h>                  // mathematics
#include "TLorentzVector.h"         // 4-vector class

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
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"
#include "../Utils/AppEffSF.cc"
#endif

//=== MAIN MACRO ================================================================================================= 

void computeAccSelZmm(const TString conf,      // input file
			    const TString inputDir,
                const TString outputDir,  // output directory
			    const Int_t   doPU,
			    const TString sysFileSIT, // condense these into 1 file per type of eff (pos & neg into 1 file)
			    const TString sysFileSta,
                const bool is13TeV=1
) {
  gBenchmark->Start("computeAccSelZmmBinned");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t MUON_MASS  = 0.105658369;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 13;
  
  const int NBptSta  = 3;
  const float ptrangeSta[NBptSta +1]   = {25., 35,  50., 10000.};

  const int NBeta = 12;
  const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};
  const int NBptSIT = 3;
  const float ptrangeSIT[NBptSIT+1] = {25,35,50,10000};
  
      const int NBptHLT = 12;
  const float ptrangeHLT[NBptHLT+1] = {25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 10000};
  
  AppEffSF effs(inputDir);
  effs.loadHLT("MuHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("MuSITEff_aMCxPythia","Combined","Combined");
  effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  effs.loadUncSel(sysFileSIT);
  effs.loadUncSta(sysFileSta);
  

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
 
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

  TH2D *h=0;

  TH2D *hSelErr_pos = new TH2D("hSelErr_pos", "",NBeta,etarange,NBptSIT,ptrangeSIT);
  TH2D *hSelErr_neg = new TH2D("hSelErr_neg", "",NBeta,etarange,NBptSIT,ptrangeSIT);
  
  TH2D *hStaErr_pos = new TH2D("hStaErr_pos", "",NBeta,etarange,NBptSta,ptrangeSta);
  TH2D *hStaErr_neg = new TH2D("hStaErr_neg", "",NBeta,etarange,NBptSta,ptrangeSta);

  TH2D *hHLTErr_pos = new TH2D("hHLTErr_pos", "",NBeta,etarange,NBptHLT,ptrangeHLT);
  TH2D *hHLTErr_neg = new TH2D("hHLTErr_neg", "",NBeta,etarange,NBptHLT,ptrangeHLT);

  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr        = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;
   
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nAllv;
  vector<Double_t> nSelCorrv, nSelCorrVarv;
  vector<Double_t> accv, accCorrv;
  vector<Double_t> accErrv, accErrCorrv;
  vector<Double_t> nSelCorrvFSR, nSelCorrvMC, nSelCorrvBkg, nSelCorrvTag;//, nSelCorrvStat;
  vector<Double_t> nSelCorrvFSR_I, nSelCorrvMC_I, nSelCorrvBkg_I, nSelCorrvTag_I;//, nSelCorrvStat_I;
  vector<Double_t> nSelCorrvFSR_S, nSelCorrvMC_S, nSelCorrvBkg_S, nSelCorrvTag_S;//, nSelCorrvStat_S;
  vector<Double_t> nSelCorrVarvFSR, nSelCorrVarvMC, nSelCorrVarvBkg, nSelCorrVarvTag;//, nSelCorrVarvStat;
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
    eventTree->SetBranchAddress("Info",             &info); TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("GenEvtInfo",        &gen); TBranch *genBr  = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle",&genPartArr); TBranch* genPartBr = eventTree->GetBranch("GenParticle");
    eventTree->SetBranchAddress("Muon",          &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");   
    eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");

    nAllv.push_back(0);
    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelCorrv.push_back(0);
    nSelCorrVarv.push_back(0);
    nSelCorrvFSR.push_back(0);  nSelCorrVarvFSR.push_back(0);
    nSelCorrvMC.push_back(0);   nSelCorrVarvMC.push_back(0);
    nSelCorrvBkg.push_back(0);  nSelCorrVarvBkg.push_back(0);
    nSelCorrvTag.push_back(0);  nSelCorrVarvTag.push_back(0);
    nSelCorrvFSR_I.push_back(0);  nSelCorrvFSR_S.push_back(0);
    nSelCorrvMC_I.push_back(0);   nSelCorrvMC_S.push_back(0);
    nSelCorrvBkg_I.push_back(0);  nSelCorrvBkg_S.push_back(0);
    nSelCorrvTag_I.push_back(0);  nSelCorrvTag_S.push_back(0);
    //
    // loop over events
    //      
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; ientry<(uint)(0.10*eventTree->GetEntries()); ientry++) {
      if(ientry%100000==0)   cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); genPartBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      Int_t glepq1=-99;
      Int_t glepq2=-99;
      bool alreadyDid=false;

      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      
      nAllv[ifile]+=gen->weight;
      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,&glepq1,&glepq2,1);
      // if((vec->M()<MASS_LOW || vec->M()>MASS_HIGH)) continue;
      delete vec; delete lep1; delete lep2;

      vertexArr->Clear();
      vertexBr->GetEntry(ientry);
      double npv  = vertexArr->GetEntries();
      Double_t weight=gen->weight;
      if(doPU>0) weight*=h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

      nEvtsv[ifile]+=weight;
      
      // trigger requirement               
      if (!isMuonTrigger(triggerMenu, info->triggerBits,kFALSE,is13TeV)) continue;

      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      muonArr->Clear();
      muonBr->GetEntry(ientry);

      for(Int_t i1=0; i1<muonArr->GetEntriesFast(); i1++) {
        const baconhep::TMuon *mu1 = (baconhep::TMuon*)((*muonArr)[i1]);

        if(mu1->pt	  < PT_CUT)  continue;  // lepton pT cut
        if(fabs(mu1->eta) > ETA_CUT) continue;  // lepton |eta| cut
        if(!passMuonID(mu1))	     continue;  // lepton selection
	
        TLorentzVector vMu1(0,0,0,0);
        vMu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, MUON_MASS);
       

        for(Int_t i2=i1+1; i2<muonArr->GetEntriesFast(); i2++) {
          const baconhep::TMuon *mu2 = (baconhep::TMuon*)((*muonArr)[i2]);
        
          if(mu1->q == mu2->q)	       continue;  // opposite charge requirement
          if(mu2->pt        < PT_CUT)  continue;  // lepton pT cut
          if(fabs(mu2->eta) > ETA_CUT) continue;  // lepton |eta| cut
          if(!passMuonID(mu2))	       continue;  // lepton selection

          TLorentzVector vMu2(0,0,0,0);
          vMu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, MUON_MASS);  

          // trigger match
           if(!isMuonTriggerObj(triggerMenu, mu1->hltMatchBits, kFALSE,is13TeV) && !isMuonTriggerObj(triggerMenu, mu2->hltMatchBits, kFALSE,is13TeV)) continue;
          // mass window
          TLorentzVector vDilep = vMu1 + vMu2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          /******** We have a Z candidate! HURRAY! ********/
          Double_t effdata, effmc, emTag;
          Double_t edFSR, edMC, edBkg, edTag;//, edStat;
          Double_t corr=1;
          Double_t corrFSR=1, corrMC=1, corrBkg=1, corrTag=1;//, corrStat=1;
          Double_t corrFSR_I=1, corrMC_I=1, corrBkg_I=1, corrTag_I=1;//, corrStat_I=1;
          Double_t corrFSR_S=1, corrMC_S=1, corrBkg_S=1, corrTag_S=1;//, corrStat_S=1;
      
          effdata=1; effmc=1;   emTag=1; 
          edFSR=1; edMC=1; edBkg=1; edTag=1;//edStat=1;  
          // }
          if(alreadyDid) continue;
          alreadyDid=true;
          
          
          
          int q1 = mu1->q;
          int q2 = mu2->q;
          
          corr = effs.fullEfficiencies(&vMu1,q1,&vMu2,q2);
          // corr = effs.dataOnly(&vMu1,q1,&vMu2,q2);
          vector<double> uncs_sta = effs.getUncSta(&vMu1,q1,&vMu2,q2);
          vector<double> uncs_sit = effs.getUncSel(&vMu1,q1,&vMu2,q2);
          
          corrFSR *= uncs_sta[0]*uncs_sit[0]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2); // alternate fsr model
          corrMC  *= uncs_sta[1]*uncs_sit[1]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2); // alternate mc gen model
          corrBkg *= uncs_sta[2]*uncs_sit[2]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2); // alternate bkg model
          corrTag *= uncs_sta[3]*uncs_sit[3]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2); // alternate bkg model
          // corr *= effdata/effmc; // orig
          
          corrFSR_I *= uncs_sit[0]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeStaSF(&vMu1,q1,&vMu2,q2); 
          corrMC_I  *= uncs_sit[1]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeStaSF(&vMu1,q1,&vMu2,q2); 
          corrBkg_I *= uncs_sit[2]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeStaSF(&vMu1,q1,&vMu2,q2); 
          corrTag_I *= uncs_sit[3]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeStaSF(&vMu1,q1,&vMu2,q2); 
           
          corrFSR_S *= uncs_sta[0]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeSelSF(&vMu1,q1,&vMu2,q2); 
          corrMC_S  *= uncs_sta[1]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeSelSF(&vMu1,q1,&vMu2,q2); 
          corrBkg_S *= uncs_sta[2]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeSelSF(&vMu1,q1,&vMu2,q2); 
          corrTag_S *= uncs_sta[3]*effs.computeHLTSF(&vMu1,q1,&vMu2,q2)*effs.computeSelSF(&vMu1,q1,&vMu2,q2); 
           
          double var=0.;        
          // var += effs.statUncSta(&l1, q1) + effs.statUncSta(&l2, q2);
          var += effs.statUncSta(&vMu1, q1, hStaErr_pos, hStaErr_neg, fabs(weight)*corr);
          var += effs.statUncSta(&vMu2, q2, hStaErr_pos, hStaErr_neg, fabs(weight)*corr);
          var += effs.statUncSel(&vMu1, q1, hSelErr_pos, hSelErr_neg, fabs(weight)*corr);
          var += effs.statUncSel(&vMu2, q2, hSelErr_pos, hSelErr_neg, fabs(weight)*corr);
          var += effs.statUncHLT(&vMu1, q1, hHLTErr_pos, hHLTErr_neg, fabs(weight)*corr);
          var += effs.statUncHLT(&vMu2, q2, hHLTErr_pos, hHLTErr_neg, fabs(weight)*corr);
          
          // std::cout << info->evtNum << " " << corr << " " << std::endl;
          nSelv[ifile]    +=weight;
          nSelCorrvFSR[ifile] +=weight*corrFSR;  nSelCorrvFSR_I[ifile] +=weight*corrFSR_I;  nSelCorrvFSR_S[ifile] +=weight*corrFSR_S;
          nSelCorrvMC[ifile]  +=weight*corrMC;   nSelCorrvMC_I[ifile]  +=weight*corrMC_I;   nSelCorrvMC_S[ifile]  +=weight*corrMC_S;
          nSelCorrvBkg[ifile] +=weight*corrBkg;  nSelCorrvBkg_I[ifile] +=weight*corrBkg_I;  nSelCorrvBkg_S[ifile] +=weight*corrBkg_S;
          nSelCorrvTag[ifile] +=weight*corrTag;  nSelCorrvTag_I[ifile] +=weight*corrTag_I;  nSelCorrvTag_S[ifile] +=weight*corrTag_S;
          
          nSelCorrv[ifile]+=weight*corr;
          nSelCorrVarvFSR[ifile]+=weight*weight*corrFSR*corrFSR;
          nSelCorrVarvMC[ifile]+=weight*weight*corrMC*corrMC;
          nSelCorrVarvBkg[ifile]+=weight*weight*corrBkg*corrBkg;
          nSelCorrVarvTag[ifile]+=weight*weight*corrTag*corrTag;
          // nSelCorrVarv[ifile]+=weight*weight*corr*corr;
        }
      }      
    }

    // std::cout << "nSelCorrVarv[ifile]  " <<  nSelCorrVarv[ifile] << std::endl;
    // std::cout << "var" << std::endl;
    Double_t var=0;
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

    nSelCorrVarvFSR[ifile]+=var;
    nSelCorrVarvMC[ifile]+=var;
    nSelCorrVarvBkg[ifile]+=var;
    nSelCorrVarvTag[ifile]+=var;
    nSelCorrVarv[ifile]+=var;
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(accv[ifile]*sqrt((1.+accv[ifile])/nEvtsv[ifile]));
    
        
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

    
  std::cout << "print output" << std::endl;
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================    
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << "     SF corrected: " << accCorrv[ifile]     << " +/- " << accErrCorrv[ifile]    << endl;
    cout << "  ==total efficiency==> " <<  setw(4) << accCorrv[ifile]/accv[ifile] << endl;
    cout << "          pct: " << 100*accErrCorrv[ifile] /accCorrv[ifile] << endl;
    cout << "          FSR unc: " << accCorrvFSR[ifile]  << " / Sel: " << accCorrvFSR_I[ifile] << " / Sta: " << accCorrvFSR_S[ifile] << endl;
    cout << "           MC unc: " << accCorrvMC[ifile]   << " / Sel: " << accCorrvMC_I[ifile]  << " / Sta: " << accCorrvMC_S[ifile]  << endl;
    cout << "          Bkg unc: " << accCorrvBkg[ifile]  << " / Sel: " << accCorrvBkg_I[ifile] << " / Sta: " << accCorrvBkg_S[ifile] << endl;
    cout << "          Tag unc: " << accCorrvTag[ifile]  << " / Sel: " << accCorrvTag_I[ifile] << " / Sta: " << accCorrvTag_S[ifile] << endl;
    // cout << "         Stat unc: " << accCorrvStat[ifile] << " / Sel: " << accCorrvStat_I[ifile]<< " / Sta: " << accCorrvStat_S[ifile] << endl;
    cout << "  fraction passing gen cut: " << nEvtsv[ifile] << " / " << nAllv[ifile] << " = " << nEvtsv[ifile]/nAllv[ifile] << endl;
    cout << endl;
  }
  
  char txtfname[500];
  sprintf(txtfname,"%s/binned.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> mu mu" << endl;
  txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "          nominal: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    txtfile << "  ==total efficiency==> " <<  setw(4) << accCorrv[ifile]/accv[ifile] << endl;
    txtfile << "     SF corrected: " << accCorrv[ifile]    << " +/- " << accErrCorrv[ifile]    << endl;
    txtfile << "          FSR unc: " << accCorrvFSR[ifile]  << " / Sel: " << accCorrvFSR_I[ifile] << " / Sta: " << accCorrvFSR_S[ifile] << endl;
    txtfile << "           MC unc: " << accCorrvMC[ifile]   << " / Sel: " << accCorrvMC_I[ifile]  << " / Sta: " << accCorrvMC_S[ifile]  << endl;
    txtfile << "          Bkg unc: " << accCorrvBkg[ifile]  << " / Sel: " << accCorrvBkg_I[ifile] << " / Sta: " << accCorrvBkg_S[ifile] << endl;
    txtfile << "          Tag unc: " << accCorrvTag[ifile]  << " / Sel: " << accCorrvTag_I[ifile] << " / Sta: " << accCorrvTag_S[ifile] << endl;
    // txtfile << "         Stat unc: " << accCorrvStat[ifile] << " / Sel: " << accCorrvStat_I[ifile]<< " / Sta: " << accCorrvStat_S[ifile] << endl;
    txtfile << "  fraction passing gen cut: " << nEvtsv[ifile] << " / " << nAllv[ifile] << " = " << nEvtsv[ifile]/nAllv[ifile] << endl;
    txtfile << endl;
  }
  txtfile.close();
  
    // char txtfname[100];
  sprintf(txtfname,"%s/sel_nums_only.txt",outputDir.Data());
  ofstream txtfile2;
  txtfile2.open(txtfname);
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
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

    txtfile3 << endl;
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

    txtfile4 << endl;
  }
  txtfile4.close();
  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZmmBinned"); 
}
