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
#include <TLorentzVector.h>               // 4-vector class

#include "../Utils/LeptonCorr.hh"	  // lepton corrections
#include "../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET
#endif

//=== MAIN MACRO ================================================================================================= 
void makeTemplatesWm() 
{
  gBenchmark->Start("makeTemplatesWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;

  // input W signal file
  TString infilename("/data/blue/Bacon/Run2/wz_flat_07_23/Wmunu/ntuples/wm_select.root");
  
  // file name with Zll data
  TString datafname("ZmmData/fits_mva.root");
  // file name with Zl MC
  TString zllMCfname("ZeeMC/fits.root");
  // file name with Wp MC
  TString wpMCfname("WmpMC/fits.root");
    // file name with Wm MC
  TString wmMCfname("WmmMC/fits.root");
  

  
  // output file name
  TString outfilename("./testWm.root");


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // Access recoil corrections
  //RecoilCorrector recoilCorr(datafname,zllMCfname,wpMCfname,wmMCfname);
  RecoilCorrector recoilCorr(datafname);

  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t weight, scale1fb, puWeight;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0;
    
  //
  // Set up output TTrees
  //
  Float_t out_met;
  TFile *outFile  = new TFile(outfilename,"RECREATE"); 
  
  TTree *rawWmTree  = new TTree("RawWm","RawWm");
  rawWmTree->Branch("weight", &weight, "weight/F");
  rawWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWmpTree  = new TTree("RawWmp","RawWmp");
  rawWmpTree->Branch("weight", &weight, "weight/F");
  rawWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWmmTree  = new TTree("RawWmm","RawWmm");
  rawWmmTree->Branch("weight", &weight, "weight/F");
  rawWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrWmTree = new TTree("CorrWm","corrWm");
  corrWmTree->Branch("weight", &weight, "weight/F");
  corrWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWmpTree = new TTree("CorrWmp","corrWmp");
  corrWmpTree->Branch("weight", &weight, "weight/F");
  corrWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWmmTree = new TTree("CorrWmm","corrWmm");
  corrWmmTree->Branch("weight", &weight, "weight/F");
  corrWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrUpWmTree = new TTree("CorrUpWm","corrUpWm");
  corrUpWmTree->Branch("weight", &weight, "weight/F");
  corrUpWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWmpTree = new TTree("CorrUpWmp","corrUpWmp");
  corrUpWmpTree->Branch("weight", &weight, "weight/F");
  corrUpWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWmmTree = new TTree("CorrUpWmm","corrUpWmm");
  corrUpWmmTree->Branch("weight", &weight, "weight/F");
  corrUpWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrDownWmTree = new TTree("CorrDownWm","corrDownWm");
  corrDownWmTree->Branch("weight", &weight, "weight/F");
  corrDownWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWmpTree = new TTree("CorrDownWmp","corrDownWmp");
  corrDownWmpTree->Branch("weight", &weight, "weight/F");
  corrDownWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWmmTree = new TTree("CorrDownWmm","corrDownWmm");
  corrDownWmmTree->Branch("weight", &weight, "weight/F");
  corrDownWmmTree->Branch("out_met",  &out_met,  "out_met/F");
  
  TTree *lepScaleUpWmTree = new TTree("LepScaleUpWm","lepScaleUpWm");
  lepScaleUpWmTree->Branch("weight", &weight, "weight/F");
  lepScaleUpWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWmpTree = new TTree("LepScaleUpWmp","lepScaleUpWmp");
  lepScaleUpWmpTree->Branch("weight", &weight, "weight/F");
  lepScaleUpWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWmmTree = new TTree("LepScaleUpWmm","lepScaleUpWmm");
  lepScaleUpWmmTree->Branch("weight", &weight, "weight/F");
  lepScaleUpWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepScaleDownWmTree = new TTree("LepScaleDownWm","lepScaleDownWm");
  lepScaleDownWmTree->Branch("weight", &weight, "weight/F");
  lepScaleDownWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWmpTree = new TTree("LepScaleDownWmp","lepScaleDownWmp");
  lepScaleDownWmpTree->Branch("weight", &weight, "weight/F");
  lepScaleDownWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWmmTree = new TTree("LepScaleDownWmm","lepScaleDownWmm");
  lepScaleDownWmmTree->Branch("weight", &weight, "weight/F");
  lepScaleDownWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResUpWmTree = new TTree("LepResUpWm","lepResUpWm");
  lepResUpWmTree->Branch("weight", &weight, "weight/F");
  lepResUpWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWmpTree = new TTree("LepResUpWmp","lepResUpWmp");
  lepResUpWmpTree->Branch("weight", &weight, "weight/F");
  lepResUpWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWmmTree = new TTree("LepResUpWmm","lepResUpWmm");
  lepResUpWmmTree->Branch("weight", &weight, "weight/F");
  lepResUpWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResDownWmTree = new TTree("LepResDownWm","lepResDownWm");
  lepResDownWmTree->Branch("weight", &weight, "weight/F");
  lepResDownWmTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWmpTree = new TTree("LepResDownWmp","lepResDownWmp");
  lepResDownWmpTree->Branch("weight", &weight, "weight/F");
  lepResDownWmpTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWmmTree = new TTree("LepResDownWmm","lepResDownWmm");
  lepResDownWmmTree->Branch("weight", &weight, "weight/F");
  lepResDownWmmTree->Branch("out_met",  &out_met,  "out_met/F");

  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  infile = TFile::Open(infilename);	  assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("runNum",   &runNum);    // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);    // event number
  intree->SetBranchAddress("npv",      &npv);	    // number of primary vertices
  intree->SetBranchAddress("npu",      &npu);	    // number of in-time PU events (MC)
  intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
  intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)   
  intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("puWeight", &puWeight);  // pileup reweighting
  intree->SetBranchAddress("mvaMet",   &met);	    // MET
  intree->SetBranchAddress("mvaMetPhi",&metPhi);    // phi(MET)
  intree->SetBranchAddress("mvaSumEt", &sumEt);     // Sum ET
  intree->SetBranchAddress("mvaMt",    &mt);	    // transverse mass
  intree->SetBranchAddress("mvaU1",    &u1);	    // parallel component of recoil
  intree->SetBranchAddress("mvaU2",    &u2);	    // perpendicular component of recoil
  intree->SetBranchAddress("q",        &q);	    // lepton charge
  intree->SetBranchAddress("lep",      &lep);	    // lepton 4-vector
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    weight=scale1fb*puWeight;
    
    if(lep->Pt()	< PT_CUT)  continue;  
    if(fabs(lep->Eta()) > ETA_CUT) continue;

    // uncorrected MET
    out_met = met;

    rawWmTree->Fill();
    if(q>0) rawWmpTree->Fill();
    else    rawWmmTree->Fill();
  	  
    Double_t corrMet=met, corrMetPhi=metPhi;
    Double_t lepPt=lep->Pt(), lepPhi=lep->Phi();
    
    // apply recoil corrections with nominal lepton scale and resolution corrections
    lepPt  = gRandom->Gaus(lep->Pt()*getMuScaleCorr(lep->Eta(),0), getMuResCorr(lep->Eta(),0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    corrWmTree->Fill();
    if(q>0) corrWmpTree->Fill();
    else    corrWmmTree->Fill();

    // recoil corrections "up"
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,1,0);
    out_met = corrMet;

    corrUpWmTree->Fill();    
    if(q>0) corrUpWmpTree->Fill();
    else    corrUpWmmTree->Fill();

    // recoil corrections "down"
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,-1,0);
    out_met = corrMet;

    corrDownWmTree->Fill();
    if(q>0) corrDownWmpTree->Fill();
    else    corrDownWmmTree->Fill();   
    
    // lepton scale "up"
    lepPt  = gRandom->Gaus(lep->Pt()*getMuScaleCorr(lep->Eta(),1), getMuResCorr(lep->Eta(),0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    lepScaleUpWmTree->Fill();
    if(q>0) lepScaleUpWmpTree->Fill();
    else    lepScaleUpWmmTree->Fill();

    // lepton scale "down"
    lepPt  = gRandom->Gaus(lep->Pt()*getMuScaleCorr(lep->Eta(),-1), getMuResCorr(lep->Eta(),0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    lepScaleDownWmTree->Fill();
    if(q>0) lepScaleDownWmpTree->Fill();
    else    lepScaleDownWmmTree->Fill();
    
    // lepton resolution "up"
    lepPt  = gRandom->Gaus(lep->Pt()*getMuScaleCorr(lep->Eta(),0), getMuResCorr(lep->Eta(),1));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    lepResUpWmTree->Fill();
    if(q>0) lepResUpWmpTree->Fill();
    else    lepResUpWmmTree->Fill();

    // lepton resolution "down"
    lepPt  = gRandom->Gaus(lep->Pt()*getMuScaleCorr(lep->Eta(),0), TMath::Max(getMuResCorr(lep->Eta(),-1),0.0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
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
