//================================================================================================
//
// Select probes for electron efficiencies with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include "TLorentzVector.h"               // 4-vector class
#include "TRandom.h"
#include "TH1.h"

#include "../Utils/LeptonCorr.hh"

// structure for output ntuple
#include "EffData.hh" 
#endif

//=== MAIN MACRO ================================================================================================= 

void selectProbesEleEff(const TString infilename,           // input ntuple
                        const TString outputDir,            // output directory
			const Int_t   effType,              // type of efficiency to compute
		        const Bool_t  doGenMatch = kFALSE,  // match to generator leptons
			const Bool_t  doWeighted = kFALSE,  // store events with weights
                        const UInt_t  desiredrunNum = 0     // select a specific run (0 for all runs)
) {
  gBenchmark->Start("selectProbesEleEff");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
 
  const Double_t TAG_PT_CUT = 25;
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  enum { eHLTEff, eL1Eff, eSelEff, eGsfEff, eGsfSelEff, eSCEff, eIDEff, eIsoEff };  // efficiency type enum
  if(effType > eIsoEff) {
    cout << "Invalid effType option! Exiting..." << endl;
    return; 
  }
  enum { eEleEle2HLT=1, eEleEle1HLT1L1, eEleEle1HLT, eEleEleNoSel, eEleSC, eTrkSC, eTrkNoSC, eEleIso, eEleNoIso };  // event category enum
  
  Double_t nProbes = 0;

  const float fitMassLo = 60;
  const float fitMassHi = 120;
  const int massBin = 60;

  const int muEtaNB = 12;
  const float muEtaRange[muEtaNB+1] = {-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4};
  const int muPtNB = 3;
  const float muPtRange[muPtNB+1] = {25,35,50,10000};
  
  TH1D *h_mll_amcnlo[muEtaNB][muPtNB];
  TH1D *h_mll_pythia[muEtaNB][muPtNB];
  TH1D *h_mll_photos[muEtaNB][muPtNB];
  
  TString rwDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/probes/";
  // Set up to load both the photos and pythia extra weighting
  TFile *inAMCNLO = new TFile(rwDir+TString("GenReweightMassHist/EleGSF/zee-amc-pyth-mll.root"),"OPEN");
  TFile *inPhotos = new TFile(rwDir+TString("GenReweightMassHist/EleGSF/zee-pow-phot-mll.root"),"OPEN");
  // TH1D * = (TH1F*)f.Get(“h1”);
  TFile *inPythia = new TFile(rwDir+TString("GenReweightMassHist/EleGSF/zee-pow-pyth-mll.root"),"OPEN"); 
    
  for(int iEtaBin = 0; iEtaBin < muEtaNB; iEtaBin ++){
    for(int iPtBin = 0; iPtBin < muPtNB; iPtBin ++){
      h_mll_amcnlo[iEtaBin][iPtBin] = (TH1D*)inAMCNLO->Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
      double norm = h_mll_amcnlo[iEtaBin][iPtBin]->Integral(); h_mll_amcnlo[iEtaBin][iPtBin]->Scale(1/norm); norm=0;
      h_mll_photos[iEtaBin][iPtBin] = (TH1D*)inPhotos->Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
      norm = h_mll_photos[iEtaBin][iPtBin]->Integral(); h_mll_photos[iEtaBin][iPtBin]->Scale(1/norm); norm=0;
      h_mll_pythia[iEtaBin][iPtBin] = (TH1D*)inPythia->Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
      norm = h_mll_pythia[iEtaBin][iPtBin]->Integral(); h_mll_pythia[iEtaBin][iPtBin]->Scale(1/norm); norm=0;
    }
  }

  
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  //EffData data;
  //outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum");
  Float_t mass, pt, eta, phi;
  Double_t weight, weightPowPhot, weightPowPyth;
  Int_t q;
  UInt_t npv, npu, passes, runNum, lumiSec, evtNum;
  outTree->Branch("mass",   &mass,   "mass/F");
  outTree->Branch("pt",     &pt,     "pt/F");
  outTree->Branch("eta",    &eta,    "eta/F");
  outTree->Branch("phi",    &phi,    "phi/F");
  outTree->Branch("weight", &weight, "weight/D");
  outTree->Branch("q",      &q,      "q/I");
  outTree->Branch("npv",    &npv,    "npv/i");
  outTree->Branch("npu",    &npu,    "npu/i");
  outTree->Branch("pass",   &passes, "pass/i");
  outTree->Branch("runNum", &runNum, "runNum/i");
  outTree->Branch("lumiSec",&lumiSec,"lumiSec/i");
  outTree->Branch("evtNum", &evtNum, "evtNum/i");
  outTree->Branch("weightPowPhot", &weightPowPhot, "weightPowPhot/D"); // powheg+photos/amcnlo+pythia weight
  outTree->Branch("weightPowPyth", &weightPowPyth, "weightPowPyth/D"); // powheg+pythia/amcnlo+pythia weight
 
  //
  // Declare input ntuple variables
  //
  //UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen, matchGenSCEff;
  UInt_t  category;
  //UInt_t  npv, npu;
  Float_t scale1fb;
  Float_t genWeight, PUWeight, prefireWeight=1;
  Float_t met, metPhi, sumEt, u1, u2;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Int_t   q1, q2;
  Int_t   glepq1, glepq2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0, *genlep1=0, *genlep2=0;
  TLorentzVector *sc1=0, *sc2=0;
  
  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  TFile *infile = new TFile(infilename);	 assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",   &runNum);     // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);    // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);     // event number
  intree->SetBranchAddress("matchGen", &matchGen);   // event has both leptons matched to MC Z->ll
  if(effType==eSCEff) intree->SetBranchAddress("matchGenSCEff", &matchGenSCEff);
  intree->SetBranchAddress("category", &category);   // dilepton category
  intree->SetBranchAddress("npv",      &npv);	     // number of primary vertices
  intree->SetBranchAddress("npu",      &npu);	     // number of in-time PU events (MC)
  intree->SetBranchAddress("genWeight",  &genWeight);
  intree->SetBranchAddress("PUWeight",   &PUWeight);
  // intree->SetBranchAddress("prefireWeight",   &prefireWeight);
  intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",      &met);	     // MET
  intree->SetBranchAddress("metPhi",   &metPhi);     // phi(MET)
  intree->SetBranchAddress("sumEt",    &sumEt);	     // Sum ET
  intree->SetBranchAddress("u1",       &u1);	     // parallel component of recoil
  intree->SetBranchAddress("u2",       &u2);	     // perpendicular component of recoil
  intree->SetBranchAddress("q1",       &q1);	     // charge of tag lepton
  intree->SetBranchAddress("q2",       &q2);	     // charge of probe lepton
  intree->SetBranchAddress("glepq1",       &glepq1);	     // charge of probe lepton
  intree->SetBranchAddress("glepq2",       &glepq2);	     // charge of probe lepton
  intree->SetBranchAddress("dilep",    &dilep);      // dilepton 4-vector
  intree->SetBranchAddress("lep1",     &lep1);       // tag lepton 4-vector
  intree->SetBranchAddress("lep2",     &lep2);       // probe lepton 4-vector
  intree->SetBranchAddress("genlep1",  &genlep1);       // tag lepton 4-vector
  intree->SetBranchAddress("genlep2",  &genlep2);       // probe lepton 4-vector
  intree->SetBranchAddress("sc1",      &sc1);        // tag Supercluster 4-vector
  intree->SetBranchAddress("sc2",      &sc2);        // probe Supercluster 4-vector 
  intree->SetBranchAddress("genVPt",   &genVPt);     // generator dilepton pT
  intree->SetBranchAddress("genVPhi",  &genVPhi);    // generator dilepton pT
  intree->SetBranchAddress("genVy",    &genVy);      // generator dilepton pT
  intree->SetBranchAddress("genVMass", &genVMass);   // generator dilepton pT
  
  /*
  	  genVPt   = tvec.Pt();
	  genVPhi  = tvec.Phi();
	  genVy    = tvec.Rapidity();
	  genVMass = tvec.M();
  */
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
  //for(UInt_t ientry=0; ientry<15; ientry++) {
    intree->GetEntry(ientry);

    if(desiredrunNum!=0 && runNum!=desiredrunNum) continue;

    if(lep1->Pt() < TAG_PT_CUT) continue;

    // check GEN match if necessary
    if(effType==eSCEff) {
      if(doGenMatch && !matchGenSCEff) continue;
    } else {
      if(doGenMatch && !matchGen) continue;
    }

    Bool_t  pass=kFALSE;
    
    if(effType==eHLTEff) {
      //
      // probe = electron passing selection
      // pass  = matched to HLT
      // * EleEle2HLT event means a passing probe, EleEle1HLT(1L1) event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eEleEle2HLT)    { pass=kTRUE; }
      else if(category==eEleEle1HLT1L1) { pass=kFALSE; }
      else if(category==eEleEle1HLT)    { pass=kFALSE; }
      else if(category==eEleEleNoSel)   { continue; }
      else                              { continue; }

    } else if(effType==eL1Eff) {
      //
      // probe = electron passing selection
      // pass  = matched to L1
      // * EleEle2HLT, EleEle1HLT1L1  event means a passing probe, EleEle1HLT event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eEleEle2HLT)    { pass=kTRUE; }
      else if(category==eEleEle1HLT1L1) { pass=kTRUE; }
      else if(category==eEleEle1HLT)    { pass=kFALSE; }
      else if(category==eEleEleNoSel)   { continue; }
      else                              { continue; }
    
    } else if(effType==eSelEff) {
      //
      // probe = electron
      // pass  = passing selection
      // * EleEle2HLT, EleEle1HLT(1L1) event means a passing probe, EleEleNoSel event means a failing probe,
      //   EleSC event does not satisfy probe requirements
      //    
      if     (category==eEleEle2HLT)    { pass=kTRUE; }
      else if(category==eEleEle1HLT1L1) { pass=kTRUE; }
      else if(category==eEleEle1HLT)    { pass=kTRUE; }
      else if(category==eEleEleNoSel)   { pass=kFALSE; }
      else                              { continue; }
    } else if(effType==eIDEff) {
      //
      // probe = electron
      // pass  = passing selection without isolation requirement
      // * EleIso, EleNoIso event means a passing probe, EleEleNoSel event means a failing probe,
      //   EleSC event does not satisfy probe requirements
      //    
      if     (category==eEleIso)      { pass=kTRUE; }
      else if(category==eEleNoIso)    { pass=kTRUE; }
      else if(category==eEleEleNoSel) { pass=kFALSE; }
      else                            { continue; }
    
    } else if(effType==eIsoEff) {
      //
      // probe = electron passing selection without the isolation requirement
      // pass  = passing the isolation requirement
      // * EleIso event means a passing probe, EleNoIso event means a failing probe,
      //   EleEleNoSel, EleSC event does not satisfy probe requirements
      //    
      if     (category==eEleIso)    { pass=kTRUE; }
      else if(category==eEleNoIso)  { pass=kFALSE; }
      else                          { continue; }
        
    } else if(effType==eGsfEff) {
      //
      // probe = supercluster
      // pass  = passing selection
      // * EleEle2HLT, EleEle1HLT(1L1), EleEleNoSel event means a passing probe,
      //   EleSC event means a failing probe
      //    
      if     (category==eEleEle2HLT)    { pass=kTRUE; }
      else if(category==eEleEle1HLT1L1) { pass=kTRUE; }
      else if(category==eEleEle1HLT)    { pass=kTRUE; }
      else if(category==eEleEleNoSel)   { pass=kTRUE; }
      else                              { pass=kFALSE; }  
    
    } else if(effType==eGsfSelEff) {
      //
      // probe = supercluster
      // pass  = passing selection
      // * EleEle2HLT, EleEle1HLT(1L1), EleEleNoSel event means a passing probe,
      //   EleSC event means a failing probe
      //    
      if     (category==eEleEle2HLT)    { pass=kTRUE; }
      else if(category==eEleEle1HLT1L1) { pass=kTRUE; }
      else if(category==eEleEle1HLT)    { pass=kTRUE; }
      else if(category==eEleEleNoSel)   { pass=kFALSE; }
      else                              { pass=kFALSE; }  
    } else if(effType==eSCEff) {
      //
      // probe = supercluster
      // pass  = passing selection
      // * EleEle2HLT, EleEle1HLT(1L1), EleEleNoSel event means a passing probe,
      //   EleSC event means a failing probe
      //    
      if     (category==eTrkSC)    { pass=kTRUE; }
      else if(category==eTrkNoSC)  { pass=kFALSE; }
      else                         { continue; }
    }
    weightPowPhot=1;
    weightPowPyth=1;
    nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
    if(doWeighted){
     Float_t geneta = -99.;
     Float_t genpt  = -999.;
     if(glepq1 > 0){
       geneta = genlep1->Eta();
       genpt  = genlep1->Pt();
     }
     else{
       geneta = genlep2->Eta();
       genpt  = genlep2->Pt();
     }
    
        // here we get the weights for the different gen shit
        // std::cout << "eta " << eta << "  pt " << pt << "  mass " << genVMass << std::endl; 
      for(int iEtaBin = 0; iEtaBin < muEtaNB; iEtaBin ++){
        if(geneta < muEtaRange[iEtaBin] || geneta > muEtaRange[iEtaBin+1]) continue;
        // std::cout << "eta bin = " << iEtaBin << std::endl;
        for(int iPtBin = 0; iPtBin < muPtNB; iPtBin ++){
          if(genpt < muPtRange[iPtBin] || genpt > muPtRange[iPtBin+1]) continue;
          // std::cout << "pt bin = " << iPtBin << std::endl;
          for(int iMassBin = 1; iMassBin < h_mll_photos[iEtaBin][iPtBin]->GetNbinsX()+1; iMassBin ++){
              // std::cout << "stupid Mass iterator " << h_mll_photos[iEtaBin][iPtBin]->GetBinLowEdge(iMassBin) << std::endl;
              if(genVMass < h_mll_photos[iEtaBin][iPtBin]->GetBinLowEdge(iMassBin)||genVMass > h_mll_photos[iEtaBin][iPtBin]->GetBinLowEdge(iMassBin+1)) continue;
              // std::cout << "mass bin = " << iMassBin << std::endl;
              if(h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin) !=0 ){
                weightPowPhot = h_mll_photos[iEtaBin][iPtBin]->GetBinContent(iMassBin)/h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin);
                // std::cout << h_mll_photos[iEtaBin][iPtBin]->GetBinContent(iMassBin) << "  " << h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin) << std::endl;
                // std::cout << "blah = " << weightPowPhot << std::endl;
                weightPowPyth = h_mll_pythia[iEtaBin][iPtBin]->GetBinContent(iMassBin)/h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin);
              }
                break;
            // h_mll_powheg[iEtaBin][iPtBin] = (TH1D*)inPowheg.Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
            // h_mll_pythia[iEtaBin][iPtBin] = (TH1D*)inPythia.Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
          }
        }
      }
    }
    
    // Fill tree
    mass    = dilep->M();
    pt      = (effType==eGsfSelEff && !pass) ? sc2->Pt()  : lep2->Pt();
    eta     = (effType==eGsfSelEff && !pass) ? sc2->Eta() : lep2->Eta();
    phi     = (effType==eGsfSelEff && !pass) ? sc2->Phi() : lep2->Phi();
    weight  = doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
    q	    = q2;
    npv	    = npv;
    npu	    = npu;
    passes  = (pass) ? 1 : 0;
    runNum  = runNum;
    lumiSec = lumiSec;
    evtNum  = evtNum;
    outTree->Fill();

    if(category==eEleEle2HLT) {
      if(lep2->Pt() < TAG_PT_CUT) continue;

      nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
      
      mass    = dilep->M();
      pt      = lep1->Pt();
      eta     = lep1->Eta();
      phi     = (effType==eGsfEff) ? sc1->Phi() : lep1->Phi();
      weight  = doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
      q	      = q1;
      npv     = npv;
      npu     = npu;
      passes  = 1;
      runNum  = runNum;
      lumiSec = lumiSec;
      evtNum  = evtNum;
      outTree->Fill();
    }
  }
  delete infile;
  infile=0, intree=0;	   


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectProbesEleEff"); 
}
