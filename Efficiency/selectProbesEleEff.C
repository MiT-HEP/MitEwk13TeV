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

#include "../Utils/LeptonCorr.hh"
#include "../EleScale/EnergyScaleCorrection_class.hh"

// structure for output ntuple
#include "EffData.hh" 
#endif

//=== MAIN MACRO ================================================================================================= 

void selectProbesEleEff(const TString infilename,           // input ntuple
                        const TString outputDir,            // output directory
			const Int_t   effType,              // type of efficiency to compute
		        const Bool_t  doGenMatch = kFALSE,  // match to generator leptons
			const Bool_t  doWeighted = kFALSE,  // store events with weights
			const Bool_t  doScaleAndSmear = kTRUE,
			const TString dataType = "",
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

  const TString corrFiles = "../EleScale/76X_16DecRereco_2015";

  //data
  EnergyScaleCorrection_class eleCorr( corrFiles.Data()); eleCorr.doScale= true; eleCorr.doSmearings =true;  

  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  //EffData data;
  //outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum");
  Float_t mass, pt, eta, phi;
  Double_t weight;
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
  //
  // Declare input ntuple variables
  //
  //UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen, matchGenSCEff;
  UInt_t  category;
  //UInt_t  npv, npu;
  Float_t scale1fb;
  Float_t genWeight, PUWeight;
  Float_t met, metPhi, sumEt, u1, u2;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  TLorentzVector *sc1=0, *sc2=0;
  Float_t r91,r92;
  
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
  intree->SetBranchAddress("scale1fb", &scale1fb);   // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",      &met);	     // MET
  intree->SetBranchAddress("metPhi",   &metPhi);     // phi(MET)
  intree->SetBranchAddress("sumEt",    &sumEt);	     // Sum ET
  intree->SetBranchAddress("u1",       &u1);	     // parallel component of recoil
  intree->SetBranchAddress("u2",       &u2);	     // perpendicular component of recoil
  intree->SetBranchAddress("q1",       &q1);	     // charge of tag lepton
  intree->SetBranchAddress("q2",       &q2);	     // charge of probe lepton
  intree->SetBranchAddress("dilep",    &dilep);      // dilepton 4-vector
  intree->SetBranchAddress("lep1",     &lep1);       // tag lepton 4-vector
  intree->SetBranchAddress("lep2",     &lep2);       // probe lepton 4-vector
  intree->SetBranchAddress("sc1",      &sc1);        // tag Supercluster 4-vector
  intree->SetBranchAddress("sc2",      &sc2);        // probe Supercluster 4-vector 
  intree->SetBranchAddress("r91",      &r91);               // r9
  intree->SetBranchAddress("r92",      &r92);               // r9

  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
  //for(UInt_t ientry=0; ientry<15; ientry++) {
    intree->GetEntry(ientry);

    if(desiredrunNum!=0 && runNum!=desiredrunNum) continue;

    if(sc1->Pt() < TAG_PT_CUT) continue;

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
    
    nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;

      if(doScaleAndSmear && !(effType==eGsfSelEff && !pass)){

        float smear1 = 0.0, scale1 = 1.0;
        float aeta1= fabs(lep1->Eta());
        float et1 = lep1->E() / cosh(aeta1);
        bool eb1 = aeta1 < 1.4442; 
        float smear2 = 0.0, scale2 = 1.0;
        float aeta2= fabs(lep2->Eta());
        float et2 = lep2->E() / cosh(aeta2);
        bool eb2 = aeta2 < 1.4442;
        float error1;
        float error2;

        if(dataType.CompareTo("data")==0){// DATA
                scale1 = eleCorr.ScaleCorrection( runNum, eb1, r91, aeta1, et1);
                scale2 = eleCorr.ScaleCorrection( runNum, eb2, r92, aeta2, et2);

                error1 = eleCorr.ScaleCorrectionUncertainty(runNum, eb1, r91, aeta1, et1);
                error2 = eleCorr.ScaleCorrectionUncertainty(runNum, eb2, r92, aeta2, et2);
         
                (*lep1) *= scale1;
                (*lep2) *= scale2;
        }else{ // MC
                smear1 = eleCorr.getSmearingSigma( runNum, eb1, r91, aeta1, et1 ,0.,0.);
                smear2 = eleCorr.getSmearingSigma( runNum, eb2, r92, aeta2, et2 ,0.,0.);

                float smearE1 = smear1 + std::hypot( eleCorr.getSmearingSigma( runNum, eb1, r91, aeta1,et1,1.,0.) -smear1,  eleCorr.getSmearingSigma( runNum,  eb1,r91, aeta1,et1,0.,1. ) -smear1 ) ;
                float smearE2 = smear2 + std::hypot( eleCorr.getSmearingSigma( runNum, eb2, r92, aeta2,et2,1.,0.), eleCorr.getSmearingSigma( runNum, eb2, r92, aeta2,et2,0.,1.) );
                double r1= gRandom->Gaus(0,1);
                double r2= gRandom->Gaus(0,1);
                double corr1 = 1.0 + smear1 * r1;
                double corr2 = 1.0 + smear2 * r2;
                error1 = lep1->E() * (1.0 + smearE1 * r1);
                error2 = lep2->E() * (1.0 + smearE2 * r2);
                (*lep1) *= corr1;
                (*lep2) *= corr2;
        } 
      } // do scale and smear*/

    // Fill tree
    mass    = dilep->M();
    pt      = (effType==eGsfSelEff && !pass) ? sc2->Pt()  : lep2->Pt();
    eta     = (effType==eGsfSelEff && !pass) ? sc2->Eta() : lep2->Eta();
    phi     = (effType==eGsfSelEff && !pass) ? sc2->Phi() : lep2->Phi();
//    pt	    = (effType==eSCEff) ? lep2->Pt()  : sc2->Pt();
//    eta	    = (effType==eSCEff) ? lep2->Eta() : sc2->Eta();
//    phi	    = (effType==eGsfEff || effType==eGsfSelEff) ? sc2->Phi() : lep2->Phi();
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
      if(sc2->Pt() < TAG_PT_CUT) continue;

      nProbes += doWeighted ? genWeight*PUWeight/std::abs(genWeight) : 1;
      
      mass    = dilep->M();
      pt      = lep1->Pt();//sc1->Pt();
      eta     = lep1->Eta();//sc1->Eta();
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
