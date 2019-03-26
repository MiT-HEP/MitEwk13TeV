//================================================================================================
//
// Select probes for muon efficiencies with Tag&Probe method
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
#include "TH1.h"
// structure for output ntuple
#include "EffData.hh" 
#endif

//=== MAIN MACRO ================================================================================================= 

void selectProbesMuEff(const TString infilename,           // input ntuple
                       const TString outputDir, 	   // output directory
		       const Int_t   effType, 	           // type of efficiency to compute
		       const Bool_t  doGenMatch = kFALSE,  // match to generator leptons
		       const Bool_t  doWeighted = kFALSE   // store events with weights
) {
  gBenchmark->Start("selectProbesMuEff");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
 
  const Double_t TAG_PT_CUT = 25;
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  enum { eHLTEff, eL1Eff, eSelEff, eTrkEff, eStaEff, eStaEff_iso, ePOGIDEff, ePOGIsoEff, eSelIsoTrkEff, eHLTSelStaEff, eHLTSelStaEff_iso };
  if(effType > eHLTSelStaEff_iso) {
    cout << "Invalid effType option! Exiting..." << endl;
    return;
  }
  
  enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum
  
  
    const float fitMassLo = 60;
  const float fitMassHi = 120;
  const int massBin = 60;

const int muEtaNB = 8;
const float muEtaRange[muEtaNB+1] = {-2.4,-2.1,-1.2,-0.9,0,0.9,1.2,2.1,2.4};
const int muPtNB = 8;
const float muPtRange[muPtNB+1] = {25,30,35,40,45,50,60,80,8000};
  
  
      TH1D *h_mll_amcnlo[muEtaNB][muPtNB];
    TH1D *h_mll_pythia[muEtaNB][muPtNB];
    TH1D *h_mll_photos[muEtaNB][muPtNB];

      // Set up to load both the photos and pythia extra weighting
    TFile *inAMCNLO = new TFile(TString("GenReweightMassHist/zmm-amc-pyth-mll.root"),"OPEN");
    TFile *inPhotos = new TFile(TString("GenReweightMassHist/zmm-pow-phot-mll.root"),"OPEN");
    // TH1D * = (TH1F*)f.Get(“h1”);
    TFile *inPythia = new TFile(TString("GenReweightMassHist/zmm-pow-pyth-mll.root"),"OPEN"); 
    
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
  
  Double_t nProbes = 0;
  
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
  // Declare output ntuple variables
  //
  //UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  //UInt_t  npv, npu;
  Float_t scale1fb;
  Float_t genWeight, PUWeight;
    Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  Int_t   glepq1, glepq2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0, *genlep1=0, *genlep2=0;
  TLorentzVector *sta1=0, *sta2=0;
  Float_t pfCombIso1, pfCombIso2;
  Float_t d01, dz1, d02, dz2;
  Float_t muNchi21,  muNchi22;
  UInt_t nPixHits1, nTkLayers1, nPixHits2, nTkLayers2;
  UInt_t nValidHits1, nMatch1, nValidHits2, nMatch2;
  UInt_t typeBits1, typeBits2;

  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  TFile *infile = new TFile(infilename);	 assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",     &runNum);      // event run number
  intree->SetBranchAddress("lumiSec",    &lumiSec);     // event lumi section
  intree->SetBranchAddress("evtNum",     &evtNum);      // event number
  intree->SetBranchAddress("matchGen",   &matchGen);    // event has both leptons matched to MC Z->ll
  intree->SetBranchAddress("category",   &category);    // dilepton category
  intree->SetBranchAddress("npv",        &npv);	        // number of primary vertices
  intree->SetBranchAddress("npu",        &npu);	        // number of in-time PU events (MC)
  intree->SetBranchAddress("genWeight",  &genWeight);
  intree->SetBranchAddress("PUWeight",   &PUWeight);
  intree->SetBranchAddress("scale1fb",   &scale1fb);    // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",        &met);	        // MET
  intree->SetBranchAddress("metPhi",     &metPhi);      // phi(MET)
  intree->SetBranchAddress("sumEt",      &sumEt);       // Sum ET
  intree->SetBranchAddress("u1",         &u1);	        // parallel component of recoil
  intree->SetBranchAddress("u2",         &u2);	        // perpendicular component of recoil
  intree->SetBranchAddress("q1",         &q1);	        // charge of tag lepton
  intree->SetBranchAddress("q2",         &q2);	        // charge of probe lepton
  intree->SetBranchAddress("glepq1",       &glepq1);	     // charge of probe lepton
  intree->SetBranchAddress("glepq2",       &glepq2);	     // charge of probe lepton
  intree->SetBranchAddress("dilep",      &dilep);       // dilepton 4-vector
  intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
  intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
  intree->SetBranchAddress("genlep1",  &genlep1);       // tag lepton 4-vector
  intree->SetBranchAddress("genlep2",  &genlep2);       // probe lepton 4-vector
  intree->SetBranchAddress("sta1",       &sta1);        // tag STA muon 4-vector
  intree->SetBranchAddress("sta2",       &sta2);        // probe STA muon 4-vector 
  intree->SetBranchAddress("pfCombIso1", &pfCombIso1);  // PF combined isolation of tag lepton
  intree->SetBranchAddress("pfCombIso2", &pfCombIso2);  // PF combined isolation of probe lepton    
  intree->SetBranchAddress("d01",        &d01); 	// transverse impact parameter of tag lepton
  intree->SetBranchAddress("d02",	 &d02); 	// transverse impact parameter of probe lepton      
  intree->SetBranchAddress("dz1",	 &dz1); 	// longitudinal impact parameter of tag lepton
  intree->SetBranchAddress("dz2",	 &dz2); 	// longitudinal impact parameter of probe lepton    
  intree->SetBranchAddress("muNchi21",	 &muNchi21);	// muon fit normalized chi^2 of tag lepton
  intree->SetBranchAddress("muNchi22",	 &muNchi22);	// muon fit normalized chi^2 of probe lepton	 
  intree->SetBranchAddress("nPixHits1",  &nPixHits1);   // number of pixel hits of tag muon
  intree->SetBranchAddress("nPixHits2",  &nPixHits2);   // number of pixel hits of probe muon
  intree->SetBranchAddress("nTkLayers1", &nTkLayers1);  // number of tracker layers of tag muon
  intree->SetBranchAddress("nTkLayers2", &nTkLayers2);  // number of tracker layers of probe muon
  intree->SetBranchAddress("nMatch1",	 &nMatch1);	// number of matched segments of tag muon
  intree->SetBranchAddress("nMatch2",	 &nMatch2);	// number of matched segments of probe muon    
  intree->SetBranchAddress("nValidHits1",&nValidHits1); // number of valid muon hits of tag muon
  intree->SetBranchAddress("nValidHits2",&nValidHits2); // number of valid muon hits of probe muon
  intree->SetBranchAddress("typeBits1",  &typeBits1);   // muon type of tag muon
  intree->SetBranchAddress("typeBits2",  &typeBits2);   // muon type of probe muon
  intree->SetBranchAddress("genVPt",   &genVPt);     // generator dilepton pT
  intree->SetBranchAddress("genVPhi",  &genVPhi);    // generator dilepton pT
  intree->SetBranchAddress("genVy",    &genVy);      // generator dilepton pT
  intree->SetBranchAddress("genVMass", &genVMass);   // generator dilepton pT

  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    if(lep1->Pt() < TAG_PT_CUT) continue;
    
    // check GEN match if necessary
    if(doGenMatch && !matchGen) continue;
    
    Bool_t  pass=kFALSE;
    Float_t m=0;
    
    if(effType==eHLTEff) {
      //
      // probe = muon passing selection
      // pass  = matched to HLT
      // * MuMu2HLT event means a passing probe, MuMu1HLT(1L1) event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kFALSE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { continue; }
      else if(category==eMuSta)       { continue; }
      else                            { continue; }
      
      m = dilep->M();

    } else if(effType==eL1Eff) {
      //
      // probe = muon passing selection
      // pass  = matched to HLT
      // * MuMu2HLT, MuMu1HLT1L1 event means a passing probe, MuMu1HLT event means a failing probe
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { continue; }
      else if(category==eMuSta)       { continue; }
      else                            { continue; }
      
      m = dilep->M();
    
    } else if(effType==eSelEff) {
      //
      // probe = GLB muon
      // pass  = passing selection
      // * MuMu2HLT, MuMu1HLT(1L1) event means a passing probe, MuMuNoSel event means a failing probe,
      //   all other categories do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kFALSE; }
      else if(category==eMuSta)       { continue; }
      else                            { continue; }
      
      m = dilep->M();
    
    } else if(effType==eTrkEff) {
      //
      // probe = STA muon
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), MuMuNoSel event means a passing probe, MuSta event means a failing probe, 
      //   MuTrk event does not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { pass=kFALSE; }
      else                            { continue; }
      
      // compute mass using probe STA muon pT
      TLorentzVector tp = *lep1 + *sta2;
      m = tp.M();    
    
    } else if(effType==eStaEff) {
      //
      // probe = tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), MuMuNoSel event means a passing probe, MuTrk event means a failing probe, 
      //   MuSta event does not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { continue; }
      else                            { pass=kFALSE; }
      
      m = dilep->M();
    
    } else if(effType==eStaEff_iso) {
      //
      // probe = isolated tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), isolated MuMuNoSel event means a passing probe, isolated MuTrk event means a failing probe, 
      //   MuSta, non-isolated MuMuNoSel, and non-isolated MuTrk events do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { continue; }
      else                            { if(pfCombIso2>0.12*(lep2->Pt())) continue; else pass=kFALSE; }
      
      m = dilep->M();
    } else if(effType==eSelIsoTrkEff) {

      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { if(pfCombIso2>0.15*(lep2->Pt())) pass=kFALSE; else pass=kTRUE; }
      else if(category==eMuMuNoSel)   { pass=kFALSE; }
      else if(category==eMuSta)       { pass=kFALSE; }
      else                            { continue; }

      m = dilep->M();

    }
 else if(effType==eHLTSelStaEff_iso) {
      //
      // probe = isolated tracker track
      // pass  = is also a GLB muon
      // * MuMu2HLT, MuMu1HLT(1L1), isolated MuMuNoSel event means a passing probe, isolated MuTrk event means a failing probe, 
      //   MuSta, non-isolated MuMuNoSel, and non-isolated MuTrk events do not satisfy probe requirements
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kFALSE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { if(pfCombIso2>0.12*(lep2->Pt())) continue; else pass=kFALSE; }
      else if(category==eMuSta)       { continue; }
      else                            { if(pfCombIso2>0.12*(lep2->Pt())) continue; else pass=kFALSE; }
      
      m = dilep->M();    
    } else if(effType==eHLTSelStaEff) {
      //
      // probe = tracker track
      // pass = is also a GLB muon
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kFALSE; }
      else if(category==eMuMu1HLT)    { pass=kFALSE; }
      else if(category==eMuMuNoSel)   { pass=kTRUE; }
      else if(category==eMuSta)       { continue; }
      else                            { pass=kFALSE; }
      
      m = dilep->M();
    } else if(effType==ePOGIDEff) {
      //
      // probe = tracker track
      // pass  = passes "tight" muon ID
      // * 
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   { 
        
	pass=kTRUE; 
        if(nTkLayers2  < 6)   pass=kFALSE;
        if(nPixHits2   < 1)   pass=kFALSE;
        if(fabs(d02)   > 0.2) pass=kFALSE;
        if(fabs(dz2)   > 0.5) pass=kFALSE;
        if(muNchi22    > 10)  pass=kFALSE;
        if(nMatch2     < 2)   pass=kFALSE;
        if(nValidHits2 < 1)   pass=kFALSE;
        if(!(typeBits2 & 1))  pass=kFALSE;
        if(!(typeBits2 & 8))  pass=kFALSE;
      }
      else if(category==eMuSta)     { continue; }
      else                          { pass=kFALSE; }
      
      m = dilep->M();    
    
    } else if(effType==ePOGIsoEff) {
      //
      // probe = "tight" muon
      // pass  = passes isolation
      // * 
      //    
      if     (category==eMuMu2HLT)    { pass=kTRUE; }
      else if(category==eMuMu1HLT1L1) { pass=kTRUE; }
      else if(category==eMuMu1HLT)    { pass=kTRUE; }
      else if(category==eMuMuNoSel)   {
        
	if(nTkLayers2  < 6)   continue;
        if(nPixHits2   < 1)   continue;
        if(fabs(d02)   > 0.2) continue;
        if(fabs(dz2)   > 0.5) continue;
        if(muNchi22    > 10)  continue;
        if(nMatch2     < 2)   continue;
        if(nValidHits2 < 1)   continue;
        if(!(typeBits2 & 1))  continue;
        if(!(typeBits2 & 8))  continue;

	pass = (pfCombIso2 < 0.12*(lep2->Pt()));      
      }
      else if(category==eMuSta)     { continue; }
      else                          { continue; }
      
      m = dilep->M();    
    }
    
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
                weightPowPhot = h_mll_photos[iEtaBin][iPtBin]->GetBinContent(iMassBin)/h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin);
                // std::cout << h_mll_photos[iEtaBin][iPtBin]->GetBinContent(iMassBin) << "  " << h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin) << std::endl;
                // std::cout << "blah = " << weightPowPhot << std::endl;
                weightPowPyth = h_mll_pythia[iEtaBin][iPtBin]->GetBinContent(iMassBin)/h_mll_amcnlo[iEtaBin][iPtBin]->GetBinContent(iMassBin);
                break;
            // h_mll_powheg[iEtaBin][iPtBin] = (TH1D*)inPowheg.Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
            // h_mll_pythia[iEtaBin][iPtBin] = (TH1D*)inPythia.Get(Form("h_mll_etaBin%d_ptBin%d", iEtaBin, iPtBin));
          }
        }
      }
    }

    // Fill tree
    mass = m;
    pt	 = (effType==eTrkEff) ? sta2->Pt()  : lep2->Pt();
    eta	 = (effType==eTrkEff) ? sta2->Eta() : lep2->Eta();
    phi	 = (effType==eTrkEff) ? sta2->Phi() : lep2->Phi();
    weight  = doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
    q	    = q2;
    npv	    = npv;
    npu	    = npu;
    passes  = (pass) ? 1 : 0;
    runNum  = runNum;
    lumiSec = lumiSec;
    evtNum  = evtNum;
    outTree->Fill();
    
    if(category==eMuMu2HLT) {
      if(lep2->Pt() < TAG_PT_CUT) continue;

      nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;

      mass	   = m;
      pt	   = (effType==eTrkEff) ? sta1->Pt()  : lep1->Pt();
      eta	   = (effType==eTrkEff) ? sta1->Eta() : lep1->Eta();
      phi	   = (effType==eTrkEff) ? sta1->Phi() : lep1->Phi();
      weight  = doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
      q	   = q1;
      npv	   = npv;
      npu	   = npu;
      passes	   = 1;
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
      
  gBenchmark->Show("selectProbesMuEff"); 
}
