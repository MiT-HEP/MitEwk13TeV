#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <TRandom3.h>
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include "Math/LorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void makeTemplatesWm() 
{
  gBenchmark->Start("makeTemplatesWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.1;

  // input W signal file
  TString infilename("/data/blue/ksung/EWKAna/8TeV/Selection/Wmunu/ntuples/wm_select.root");
  
  // file name with recoil correction
  TString recoilfname("ZmmDataMay23_2/fits.root");
  
  // output file name
  TString outfilename("/data/blue/ksung/EWKAna/8TeV/RecoilSyst/WmTemplates_May23_2.root");
  

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
    
  // Access recoil corrections
  RecoilCorrector recoilCorr(recoilfname);

  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  LorentzVector *lep=0;
    
  //
  // Set up output TTrees
  //
  Float_t out_met;
  TFile *outFile  = new TFile(outfilename,"RECREATE"); 
  
  TTree *rawWmTree  = new TTree("RawWm","RawWm");
  rawWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  rawWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWmpTree  = new TTree("RawWmp","RawWmp");
  rawWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  rawWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWmmTree  = new TTree("RawWmm","RawWmm");
  rawWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  rawWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrWmTree = new TTree("CorrWm","corrWm");
  corrWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWmpTree = new TTree("CorrWmp","corrWmp");
  corrWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWmmTree = new TTree("CorrWmm","corrWmm");
  corrWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrUpWmTree = new TTree("CorrUpWm","corrUpWm");
  corrUpWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrUpWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWmpTree = new TTree("CorrUpWmp","corrUpWmp");
  corrUpWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrUpWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWmmTree = new TTree("CorrUpWmm","corrUpWmm");
  corrUpWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrUpWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrDownWmTree = new TTree("CorrDownWm","corrDownWm");
  corrDownWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrDownWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWmpTree = new TTree("CorrDownWmp","corrDownWmp");
  corrDownWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrDownWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWmmTree = new TTree("CorrDownWmm","corrDownWmm");
  corrDownWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrDownWmmTree->Branch("out_met",  &out_met,  "out_met/F");
  
  TTree *lepScaleUpWmTree = new TTree("LepScaleUpWm","lepScaleUpWm");
  lepScaleUpWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleUpWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWmpTree = new TTree("LepScaleUpWmp","lepScaleUpWmp");
  lepScaleUpWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleUpWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWmmTree = new TTree("LepScaleUpWmm","lepScaleUpWmm");
  lepScaleUpWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleUpWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepScaleDownWmTree = new TTree("LepScaleDownWm","lepScaleDownWm");
  lepScaleDownWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleDownWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWmpTree = new TTree("LepScaleDownWmp","lepScaleDownWmp");
  lepScaleDownWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleDownWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWmmTree = new TTree("LepScaleDownWmm","lepScaleDownWmm");
  lepScaleDownWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleDownWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResUpWmTree = new TTree("LepResUpWm","lepResUpWm");
  lepResUpWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResUpWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWmpTree = new TTree("LepResUpWmp","lepResUpWmp");
  lepResUpWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResUpWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWmmTree = new TTree("LepResUpWmm","lepResUpWmm");
  lepResUpWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResUpWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResDownWmTree = new TTree("LepResDownWm","lepResDownWm");
  lepResDownWmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResDownWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWmpTree = new TTree("LepResDownWmp","lepResDownWmp");
  lepResDownWmpTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResDownWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWmmTree = new TTree("LepResDownWmm","lepResDownWmm");
  lepResDownWmmTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResDownWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  infile = new TFile(infilename);	  assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",   &runNum);    // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);    // event number
  intree->SetBranchAddress("npv",      &npv);	    // number of primary vertices
  intree->SetBranchAddress("npu",      &npu);	    // number of in-time PU events (MC)
  intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
  intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)   
  intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("met",      &met);	    // MET
  intree->SetBranchAddress("metPhi",   &metPhi);    // phi(MET)
  intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
  intree->SetBranchAddress("mt",       &mt);	    // transverse mass
  intree->SetBranchAddress("u1",       &u1);	    // parallel component of recoil
  intree->SetBranchAddress("u2",       &u2);	    // perpendicular component of recoil
  intree->SetBranchAddress("q",        &q);	    // lepton charge
  intree->SetBranchAddress("lep",      &lep);	    // lepton 4-vector
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(lep->Pt()	< PT_CUT)  continue;  
    if(fabs(lep->Eta()) > ETA_CUT) continue;

    out_met = met;
    rawWmTree->Fill();
    if(q>0) rawWmpTree->Fill();
    else    rawWmmTree->Fill();
  	  
    Double_t corrMet=met, corrMetPhi=metPhi;
    Double_t lepPt=lep->Pt(), lepPhi=lep->Phi();
    
    Double_t lepScale=1;
    Double_t lepScaleErr=0;
    Double_t lepRes=0;
    Double_t lepResErr=0;
    
    lepScale=1; lepScaleErr=0.005;
    lepRes=0.5; lepResErr=0.5;
    
    // apply recoil corrections with nominal lepton scale and resolution corrections
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    corrWmTree->Fill();
    if(q>0) corrWmpTree->Fill();
    else    corrWmmTree->Fill();    
    
    // recoil corrections "up"
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,1);
    out_met = corrMet;
    corrUpWmTree->Fill();
    if(q>0) corrUpWmpTree->Fill();
    else    corrUpWmmTree->Fill();
  	  
    // recoil corrections "down"
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,-1);
    out_met = corrMet;
    corrDownWmTree->Fill();
    if(q>0) corrDownWmpTree->Fill();
    else    corrDownWmmTree->Fill();    
    
    // lepton scale "up"
    lepPt  = gRandom->Gaus(lep->Pt()*(lepScale+lepScaleErr), lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepScaleUpWmTree->Fill();
    if(q>0) lepScaleUpWmpTree->Fill();
    else    lepScaleUpWmmTree->Fill();
  	  
    // lepton scale "down"
    lepPt  = gRandom->Gaus(lep->Pt()*(lepScale-lepScaleErr), lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepScaleDownWmTree->Fill();
    if(q>0) lepScaleDownWmpTree->Fill();
    else    lepScaleDownWmmTree->Fill();
     
    
    // lepton resolution "up"
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes+lepResErr);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepResUpWmTree->Fill();
    if(q>0) lepResUpWmpTree->Fill();
    else    lepResUpWmmTree->Fill();
  	  
    // lepton resolution "down"
    lepPt  = (lepRes-lepResErr>0) ? gRandom->Gaus(lep->Pt()*lepScale, lepRes-lepResErr) : lep->Pt()*lepScale;
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepResDownWmTree->Fill();
    if(q>0) lepResDownWmpTree->Fill();
    else    lepResDownWmmTree->Fill();     
  }   
  delete infile;
  infile=0, intree=0;   
  

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output: " << outfilename << endl;    
  cout << endl;     
  
  gBenchmark->Show("makeTemplatesWm");
}
