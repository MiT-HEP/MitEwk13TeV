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
#include <TLorentzVector.h>

#include "../Utils/LeptonCorr.hh"         // lepton corrections
#include "../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET
#endif


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
  TString infilename("/data/blue/Bacon/Run2/wz_flat_07_23/Wenu/ntuples/we_select.root");
  
  // file name with Zll data
  TString datafname("ZeeData/fits_mva.root");
  // file name with Zl MC
  TString zllMCfname("ZeeMC/fits.root");
  // file name with Wp MC
  TString wpMCfname("WepMC/fits.root");
    // file name with Wm MC
  TString wmMCfname("WemMC/fits.root");
  
  // output file name
  TString outfilename("./testWe.root");


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
  TLorentzVector *sc=0;

  //
  // Set up output TTrees
  //
  Float_t out_met;
  TFile *outFile  = new TFile(outfilename,"RECREATE"); 
  
  TTree *rawWeTree  = new TTree("RawWe","RawWe");
  rawWeTree->Branch("weight", &weight, "weight/F");
  rawWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWepTree  = new TTree("RawWep","RawWep");
  rawWepTree->Branch("weight", &weight, "weight/F");
  rawWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *rawWemTree  = new TTree("RawWem","RawWem");
  rawWemTree->Branch("weight", &weight, "weight/F");
  rawWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrWeTree = new TTree("CorrWe","corrWe");
  corrWeTree->Branch("weight", &weight, "weight/F");
  corrWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWepTree = new TTree("CorrWep","corrWep");
  corrWepTree->Branch("weight", &weight, "weight/F");
  corrWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrWemTree = new TTree("CorrWem","corrWem");
  corrWemTree->Branch("weight", &weight, "weight/F");
  corrWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrUpWeTree = new TTree("CorrUpWe","corrUpWe");
  corrUpWeTree->Branch("weight", &weight, "weight/F");
  corrUpWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWepTree = new TTree("CorrUpWep","corrUpWep");
  corrUpWepTree->Branch("weight", &weight, "weight/F");
  corrUpWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrUpWemTree = new TTree("CorrUpWem","corrUpWem");
  corrUpWemTree->Branch("weight", &weight, "weight/F");
  corrUpWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *corrDownWeTree = new TTree("CorrDownWe","corrDownWe");
  corrDownWeTree->Branch("weight", &weight, "weight/F");
  corrDownWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWepTree = new TTree("CorrDownWep","corrDownWep");
  corrDownWepTree->Branch("weight", &weight, "weight/F");
  corrDownWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *corrDownWemTree = new TTree("CorrDownWem","corrDownWem");
  corrDownWemTree->Branch("weight", &weight, "weight/F");
  corrDownWemTree->Branch("out_met",  &out_met,  "out_met/F");
  
  TTree *lepScaleUpWeTree = new TTree("LepScaleUpWe","lepScaleUpWe");
  lepScaleUpWeTree->Branch("weight", &weight, "weight/F");
  lepScaleUpWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWepTree = new TTree("LepScaleUpWep","lepScaleUpWep");
  lepScaleUpWepTree->Branch("weight", &weight, "weight/F");
  lepScaleUpWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleUpWemTree = new TTree("LepScaleUpWem","lepScaleUpWem");
  lepScaleUpWemTree->Branch("weight", &weight, "weight/F");
  lepScaleUpWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepScaleDownWeTree = new TTree("LepScaleDownWe","lepScaleDownWe");
  lepScaleDownWeTree->Branch("weight", &weight, "weight/F");
  lepScaleDownWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWepTree = new TTree("LepScaleDownWep","lepScaleDownWep");
  lepScaleDownWepTree->Branch("weight", &weight, "weight/F");
  lepScaleDownWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepScaleDownWemTree = new TTree("LepScaleDownWem","lepScaleDownWem");
  lepScaleDownWemTree->Branch("weight", &weight, "weight/F");
  lepScaleDownWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResUpWeTree = new TTree("LepResUpWe","lepResUpWe");
  lepResUpWeTree->Branch("weight", &weight, "weight/F");
  lepResUpWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWepTree = new TTree("LepResUpWep","lepResUpWep");
  lepResUpWepTree->Branch("weight", &weight, "weight/F");
  lepResUpWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResUpWemTree = new TTree("LepResUpWem","lepResUpWem");
  lepResUpWemTree->Branch("weight", &weight, "weight/F");
  lepResUpWemTree->Branch("out_met",  &out_met,  "out_met/F");

  TTree *lepResDownWeTree = new TTree("LepResDownWe","lepResDownWe");
  lepResDownWeTree->Branch("weight", &weight, "weight/F");
  lepResDownWeTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWepTree = new TTree("LepResDownWep","lepResDownWep");
  lepResDownWepTree->Branch("weight", &weight, "weight/F");
  lepResDownWepTree->Branch("out_met",  &out_met,  "out_met/F");
  TTree *lepResDownWemTree = new TTree("LepResDownWem","lepResDownWem");
  lepResDownWemTree->Branch("weight", &weight, "weight/F");
  lepResDownWemTree->Branch("out_met",  &out_met,  "out_met/F");
    
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
  intree->SetBranchAddress("sc",       &sc);	    // electron Supercluster 4-vector
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    weight=scale1fb*puWeight;
    
    if(sc->Pt()        < PT_CUT)  continue;   
    if(fabs(sc->Eta()) > ETA_CUT) continue;

    // uncorrected MET  
    out_met = met;

    rawWeTree->Fill();
    if(q>0) rawWepTree->Fill();
    else    rawWemTree->Fill();
  	  
    Double_t corrMet=met, corrMetPhi=metPhi;
    Double_t lepPt=lep->Pt(), lepPhi=lep->Phi();
    
    // apply recoil corrections with nominal lepton scale and resolution corrections
    lepPt  = gRandom->Gaus(lep->Pt()*getEleScaleCorr(lep->Eta(),0), getEleResCorr(lep->Eta(),0)); 
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    corrWeTree->Fill();
    if(q>0) corrWepTree->Fill();
    else    corrWemTree->Fill();
    
    // recoil corrections "up"
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,1,0);
    out_met = corrMet;

    corrUpWeTree->Fill();
    if(q>0) corrUpWepTree->Fill();
    else    corrUpWemTree->Fill();

    // recoil corrections "down"
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,-1,0);
    out_met = corrMet;

    corrDownWeTree->Fill();
    if(q>0) corrDownWepTree->Fill();
    else    corrDownWemTree->Fill();    

    // lepton scale "up"
    lepPt  = gRandom->Gaus(lep->Pt()*getEleScaleCorr(lep->Eta(),1), getEleResCorr(lep->Eta(),0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    lepScaleUpWeTree->Fill();
    if(q>0) lepScaleUpWepTree->Fill();
    else    lepScaleUpWemTree->Fill();

    // lepton scale "down"
    lepPt  = gRandom->Gaus(lep->Pt()*getEleScaleCorr(lep->Eta(),-1), getEleResCorr(lep->Eta(),0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    lepScaleDownWeTree->Fill();
    if(q>0) lepScaleDownWepTree->Fill();
    else    lepScaleDownWemTree->Fill();

    // lepton resolution "up"
    lepPt  = gRandom->Gaus(lep->Pt()*getEleScaleCorr(lep->Eta(),0), getEleResCorr(lep->Eta(),1));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
    out_met = corrMet;

    lepResUpWeTree->Fill();
    if(q>0) lepResUpWepTree->Fill();
    else    lepResUpWemTree->Fill();

    // lepton resolution "down"
    lepPt  = gRandom->Gaus(lep->Pt()*getEleScaleCorr(lep->Eta(),0), TMath::Max(getEleResCorr(lep->Eta(),-1),0.0));
    recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lepPhi,0,0);
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
