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

void makeTemplatesWe() 
{
  gBenchmark->Start("makeTemplatesWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.5;

  // input W signal file
  TString infilename("/data/blue/ksung/EWKAna/8TeV/Selection/Wenu/ntuples/we_select.root");
  
  // file name with recoil correction
  TString recoilfname("ZeeDataMay23_2/fits.root");
  
  // output file name
  TString outfilename("/data/blue/ksung/EWKAna/8TeV/RecoilSyst/WeTemplates_May23_2.root");


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
  LorentzVector *sc=0;

  //
  // Set up output TTrees
  //
  Float_t out_met;
  TFile *outFile  = new TFile(outfilename,"RECREATE"); 
  
  TTree *rawWeTree  = new TTree("RawWe","RawWe");
  rawWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  rawWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWepTree  = new TTree("RawWep","RawWep");
  rawWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  rawWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWemTree  = new TTree("RawWem","RawWem");
  rawWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  rawWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrWeTree = new TTree("CorrWe","corrWe");
  corrWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWepTree = new TTree("CorrWep","corrWep");
  corrWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWemTree = new TTree("CorrWem","corrWem");
  corrWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrUpWeTree = new TTree("CorrUpWe","corrUpWe");
  corrUpWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrUpWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWepTree = new TTree("CorrUpWep","corrUpWep");
  corrUpWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrUpWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWemTree = new TTree("CorrUpWem","corrUpWem");
  corrUpWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrUpWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrDownWeTree = new TTree("CorrDownWe","corrDownWe");
  corrDownWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrDownWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWepTree = new TTree("CorrDownWep","corrDownWep");
  corrDownWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrDownWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWemTree = new TTree("CorrDownWem","corrDownWem");
  corrDownWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  corrDownWemTree->Branch("out_met",  &out_met,  "out_met/F");
  
  TTree *lepScaleUpWeTree = new TTree("LepScaleUpWe","lepScaleUpWe");
  lepScaleUpWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleUpWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWepTree = new TTree("LepScaleUpWep","lepScaleUpWep");
  lepScaleUpWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleUpWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWemTree = new TTree("LepScaleUpWem","lepScaleUpWem");
  lepScaleUpWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleUpWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepScaleDownWeTree = new TTree("LepScaleDownWe","lepScaleDownWe");
  lepScaleDownWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleDownWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWepTree = new TTree("LepScaleDownWep","lepScaleDownWep");
  lepScaleDownWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleDownWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWemTree = new TTree("LepScaleDownWem","lepScaleDownWem");
  lepScaleDownWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepScaleDownWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResUpWeTree = new TTree("LepResUpWe","lepResUpWe");
  lepResUpWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResUpWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWepTree = new TTree("LepResUpWep","lepResUpWep");
  lepResUpWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResUpWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWemTree = new TTree("LepResUpWem","lepResUpWem");
  lepResUpWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResUpWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResDownWeTree = new TTree("LepResDownWe","lepResDownWe");
  lepResDownWeTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResDownWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWepTree = new TTree("LepResDownWep","lepResDownWep");
  lepResDownWepTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResDownWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWemTree = new TTree("LepResDownWem","lepResDownWem");
  lepResDownWemTree->Branch("scale1fb", &scale1fb, "scale1fb/F");
  lepResDownWemTree->Branch("out_met",  &out_met,  "out_met/F");
    
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
  intree->SetBranchAddress("sc",       &sc);	    // electron Supercluster 4-vector
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(sc->Pt()        < PT_CUT)  continue;   
    if(fabs(sc->Eta()) > ETA_CUT) continue;
  
    out_met = met;
    rawWeTree->Fill();
    if(q>0) rawWepTree->Fill();
    else    rawWemTree->Fill();
  	  
    Double_t corrMet=met, corrMetPhi=metPhi;
    Double_t lepPt=lep->Pt(), lepPhi=lep->Phi();
    
    Double_t lepScale=0;
    Double_t lepScaleErr=0;
    Double_t lepRes=0;
    Double_t lepResErr=0;
/* *** PR ***
    if     (fabs(sc->Eta()) < 0.4)    { lepScale = 1.0/0.99968;  lepScaleErr = 0.00141821*lepScale*lepScale; lepRes = 0.511222; lepResErr = 0.188401; }
    else if(fabs(sc->Eta()) < 0.8)    { lepScale = 1.0/1.00188;  lepScaleErr = 0.00142247*lepScale*lepScale; lepRes = 0.379167; lepResErr = 0.197064; }
    else if(fabs(sc->Eta()) < 1.2)    { lepScale = 1.0/1.00239;  lepScaleErr = 0.00185596*lepScale*lepScale; lepRes = 0.732184; lepResErr = 0.289305; }
    else if(fabs(sc->Eta()) < 1.4442) { lepScale = 1.0/1.00401;  lepScaleErr = 0.00272305*lepScale*lepScale; lepRes = 0.601056; lepResErr = 0.340967; }
    else if(fabs(sc->Eta()) < 1.566)  { lepScale = 1.0/0.997551; lepScaleErr = 0.00276288*lepScale*lepScale; lepRes = 1.385;    lepResErr = 0.260686; }
    else                              { lepScale = 1.0/0.984124; lepScaleErr = 0.00247995*lepScale*lepScale; lepRes = 1.19933;  lepResErr = 0.267988; }
//*/ 
//* *** May23 ***
    if     (fabs(sc->Eta()) < 0.4)    { lepScale = 1.0/1.00243;  lepScaleErr = 0.00138066*lepScale*lepScale; lepRes = 0.464061;   lepResErr = 0.184892; }
    else if(fabs(sc->Eta()) < 0.8)    { lepScale = 1.0/1.00549;  lepScaleErr = 0.00117978*lepScale*lepScale; lepRes = 0.00985329; lepResErr = 0.820835; }
    else if(fabs(sc->Eta()) < 1.2)    { lepScale = 1.0/1.00634;  lepScaleErr = 0.00192368*lepScale*lepScale; lepRes = 0.822958;   lepResErr = 0.2599;   }
    else if(fabs(sc->Eta()) < 1.4442) { lepScale = 1.0/1.00669;  lepScaleErr = 0.00279439*lepScale*lepScale; lepRes = 0.71369;    lepResErr = 0.282484; }
    else if(fabs(sc->Eta()) < 1.566)  { lepScale = 1.0/0.998547; lepScaleErr = 0.00275456*lepScale*lepScale; lepRes = 1.35987;    lepResErr = 0.257267; }
    else                              { lepScale = 1.0/0.983544; lepScaleErr = 0.00248148*lepScale*lepScale; lepRes = 1.27686;    lepResErr = 0.248619; }
//*/    

    // apply recoil corrections with nominal lepton scale and resolution corrections
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    corrWeTree->Fill();
    if(q>0) corrWepTree->Fill();
    else    corrWemTree->Fill();
    
    // recoil corrections "up"
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,1);
    out_met = corrMet;
    corrUpWeTree->Fill();
    if(q>0) corrUpWepTree->Fill();
    else    corrUpWemTree->Fill();
  	  
    // recoil corrections "down"
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,-1);
    out_met = corrMet;
    corrDownWeTree->Fill();
    if(q>0) corrDownWepTree->Fill();
    else    corrDownWemTree->Fill();    
    
    // lepton scale "up"
    lepPt  = gRandom->Gaus(lep->Pt()*(lepScale+lepScaleErr), lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepScaleUpWeTree->Fill();
    if(q>0) lepScaleUpWepTree->Fill();
    else    lepScaleUpWemTree->Fill();
  	  
    // lepton scale "down"
    lepPt  = gRandom->Gaus(lep->Pt()*(lepScale-lepScaleErr), lepRes);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepScaleDownWeTree->Fill();
    if(q>0) lepScaleDownWepTree->Fill();
    else    lepScaleDownWemTree->Fill();
     
    
    // lepton resolution "up"
    lepPt  = gRandom->Gaus(lep->Pt()*lepScale, lepRes+lepResErr);
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepResUpWeTree->Fill();
    if(q>0) lepResUpWepTree->Fill();
    else    lepResUpWemTree->Fill();
  	  
    // lepton resolution "down"
    lepPt  = (lepRes-lepResErr>0) ? gRandom->Gaus(lep->Pt()*lepScale, lepRes-lepResErr) : lep->Pt()*lepScale;
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi);
    out_met = corrMet;
    lepResDownWeTree->Fill();
    if(q>0) lepResDownWepTree->Fill();
    else    lepResDownWemTree->Fill();         
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
  
  gBenchmark->Show("makeTemplatesWe");
}
