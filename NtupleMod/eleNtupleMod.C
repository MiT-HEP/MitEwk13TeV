//================================================================================================
//
//  Perform recoil corrections, rochester corrections, make new branches for the info
//
//  * outputs another ntuple but with new branches
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
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

#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_asym2.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"

//helper class to handle rochester corrections
#include <../RochesterCorr/RoccoR.cc>



#endif

//=== MAIN MACRO ================================================================================================= 

void eleNtupleMod(const TString  outputDir,   // output directory 
                   const TString  inputDir,    // input directory
                   const TString  fileName    // both the input and output final file name i.e. data_select.root
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  // flage to control applying the recoil corrections
  bool doInclusive = true; // This should be the standard recoil correction: 3-gaussian inclusive eta
  bool doKeys = true; // RooKeysPDF instead of 3-Gaus
  bool doEta = true; // eta-binned 3-Gaus fit
  bool doStat = true; //  Statistical Uncertainty
  
  // which MET type we use
  bool doPF = true;
  
  std::string u1_name; std::string u2_name;
  std::string met_name; std::string metPhi_name;
//   std::string recoilType;


  if(doPF){
    u1_name = "u1";
    u2_name = "u2";
    met_name = "met";
    metPhi_name = "metPhi";
//     recoilType = "PF";
  } else {
    u1_name = "puppiU1";
    u2_name = "puppiU2";
    met_name = "puppiMet";
    metPhi_name = "puppiMetPhi";
//     recoilType = "Puppi";
  }
  
  // don't think these are really necessary but leaving them for now
  
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t mu_MASS = 0.1057;
 
 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Point to the Efficiency SF
 // ------------------------------------------------------------------------------------------------------------------------------------------
 // These are the official recoil corrections path, I will update corrections by putting new ones in this folder
  const TString baseDir = "/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV/results/Zee/";
  const TString dataHLTEffName_pos = baseDir + "Data/EleHLTEff_aMCxPythia/Positive/eff.root";
  const TString dataHLTEffName_neg = baseDir + "Data/EleHLTEff_aMCxPythia/Negative/eff.root";
  const TString zeeHLTEffName_pos  = baseDir + "MC/EleHLTEff_aMCxPythia/Positive/eff.root";
  const TString zeeHLTEffName_neg  = baseDir + "MC/EleHLTEff_aMCxPythia/Negative/eff.root";

  const TString dataGsfSelEffName_pos = baseDir + "Data/EleGSFSelEff_aMCxPythia/Positive/eff.root";
  const TString dataGsfSelEffName_neg = baseDir + "Data/EleGSFSelEff_aMCxPythia/Negative/eff.root";
  const TString zeeGsfSelEffName_pos  = baseDir + "MC/EleGSFSelEff_aMCxPythia/Positive/eff.root";
  const TString zeeGsfSelEffName_neg  = baseDir + "MC/EleGSFSelEff_aMCxPythia/Negative/eff.root";

  
   const TString dataHLTEff2BinName_pos = baseDir + "Data/EleHLTEff_aMCxPythia/Positive/eff.root";
  const TString dataHLTEff2BinName_neg = baseDir + "Data/EleHLTEff_aMCxPythia/Negative/eff.root";
  const TString zeeHLTEff2BinName_pos  = baseDir + "MC/EleHLTEff_aMCxPythia/Positive/eff.root";
  const TString zeeHLTEff2BinName_neg  = baseDir + "MC/EleHLTEff_aMCxPythia/Negative/eff.root";

  const TString dataGsfSelEff2BinName_pos = baseDir + "Data/EleGSFSelEff_aMCxPythia/Positive/eff.root";
  const TString dataGsfSelEff2BinName_neg = baseDir + "Data/EleGSFSelEff_aMCxPythia/Negative/eff.root";
  const TString zeeGsfSelEff2BinName_pos  = baseDir + "MC/EleGSFSelEff_aMCxPythia/Positive/eff.root";
  const TString zeeGsfSelEff2BinName_neg  = baseDir + "MC/EleGSFSelEff_aMCxPythia/Negative/eff.root";
  
  // Efficiency Uncertainties
  TFile *fSysGSFSel = TFile::Open(SysFileGSFSel);
  TH2D  *hSysGSFSelSigFSRNeg = (TH2D*) fSysGSFSel->Get("hEleGSFSelEffSigFSRNeg");
  TH2D  *hSysGSFSelSigFSRPos = (TH2D*) fSysGSFSel->Get("hEleGSFSelEffSigFSRPos"); 
  TH2D  *hSysGSFSelSigMCNeg  = (TH2D*) fSysGSFSel->Get("hEleGSFSelEffSigMCNeg"); 
  TH2D  *hSysGSFSelSigMCPos  = (TH2D*) fSysGSFSel->Get("hEleGSFSelEffSigMCPos"); 
  TH2D  *hSysGSFSelBkgNeg    = (TH2D*) fSysGSFSel->Get("hEleGSFSelEffBkgNeg"); 
  TH2D  *hSysGSFSelBkgPos    = (TH2D*) fSysGSFSel->Get("hEleGSFSelEffBkgPos"); 


  Bool_t isData = (fileName.CompareTo("data_select.root")==0);
  std::cout << fileName.CompareTo("data_select.root",TString::kIgnoreCase) << std::endl;
  Bool_t isRecoil = (fileName.CompareTo("we_select.root")==0||fileName.CompareTo("we0_select.root")==0||fileName.CompareTo("we1_select.root")==0||fileName.CompareTo("we2_select.root")==0||fileName.CompareTo("wx_select.root")==0||fileName.CompareTo("zxx_select.root")==0);
  
  std::cout << "isData" << isData << std::endl;
  std::cout << "isRecoil " << isRecoil << std::endl;

  if(isData) {doInclusive = false; doKeys = false; doEta = false; doStat = false;}

 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Load the Recoil Correction Files
 // ------------------------------------------------------------------------------------------------------------------------------------------
  // ===================== Recoil correction files ============================
  const TString directory("/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Recoil");
 
  // New Recoil Correctors for everything
  RecoilCorrector *rcMainWp    = new  RecoilCorrector("",""); RecoilCorrector *rcMainWm    = new  RecoilCorrector("","");
  RecoilCorrector *rcStatW     = new  RecoilCorrector("",""); 
  RecoilCorrector *rcKeysWp    = new  RecoilCorrector("",""); RecoilCorrector *rcKeysWm    = new  RecoilCorrector("","");
  RecoilCorrector *rcEta05Wp   = new  RecoilCorrector("",""); RecoilCorrector *rcEta05Wm   = new  RecoilCorrector("","");
  RecoilCorrector *rcEta051Wp  = new  RecoilCorrector("",""); RecoilCorrector *rcEta051Wm  = new  RecoilCorrector("","");
  RecoilCorrector *rcEta1Wp    = new  RecoilCorrector("",""); RecoilCorrector *rcEta1Wm    = new  RecoilCorrector("","");
  // also make sure to add the Wm stuff
  if(doInclusive){
    rcMainWp->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_WmpMCPF/",directory.Data()));
    rcMainWp->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()));
    rcMainWp->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
    
    rcMainWm->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_WmmMCPF/",directory.Data()));
    rcMainWm->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()));
    rcMainWm->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
  } 
  if (doStat){
    int rec_sig = 1;
    rcStatW->loadRooWorkspacesDiagMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_2G/",directory.Data()), rec_sig);
    rcStatW->loadRooWorkspacesDiagData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_2G/",directory.Data()), rec_sig);
    rcStatW->loadRooWorkspacesDiagMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_2G/",directory.Data()), rec_sig);
  } 
  if (doEta){
    rcEta05Wp->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_0-05_v1/",directory.Data()));
    rcEta05Wp->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_0-05_v1/",directory.Data()));
    rcEta05Wp->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_0-05_v1/",directory.Data()));
    
    rcEta051Wp->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_05-10_v1/",directory.Data()));
    rcEta051Wp->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_05-10_v1/",directory.Data()));
    rcEta051Wp->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_05-10_v1/",directory.Data()));
      
    rcEta1Wp->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_10-24_v1/",directory.Data()));
    rcEta1Wp->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_10-24_v1/",directory.Data()));
    rcEta1Wp->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_10-24_v1/",directory.Data()));

    rcEta05Wm->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_0-05_v1/",directory.Data()));
    rcEta05Wm->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_0-05_v1/",directory.Data()));
    rcEta05Wm->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_0-05_v1/",directory.Data()));
    
    rcEta051Wm->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_05-10_v1/",directory.Data()));
    rcEta051Wm->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_05-10_v1/",directory.Data()));
    rcEta051Wm->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_05-10_v1/",directory.Data()));
    
    rcEta1Wm->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_10-24_v1/",directory.Data()));
    rcEta1Wm->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_10-24_v1/",directory.Data()));
      rcEta1Wm->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_10-24_v1/",directory.Data()));
  }
  if (doKeys){
    rcKeysWp->loadRooWorkspacesMCtoCorrectKeys(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_keys_v1/",directory.Data()));
    rcKeysWp->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()));
    rcKeysWp->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
    
    rcKeysWm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_keys_v1/",directory.Data()));
    rcKeysWm->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()));
    rcKeysWm->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
  }
  


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Load the Efficiency SF
 // ------------------------------------------------------------------------------------------------------------------------------------------
 //
  // HLT efficiency
  //
  cout << "Loading trigger efficiencies..." << endl;

  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_pos = new TFile(zeeHLTEffName_pos);
  CEffUser2D zeeHLTEff_pos;
  zeeHLTEff_pos.loadEff((TH2D*)zeeHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEffFile_neg = new TFile(zeeHLTEffName_neg);
  CEffUser2D zeeHLTEff_neg;
  zeeHLTEff_neg.loadEff((TH2D*)zeeHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeHLTEffFile_neg->Get("hErrhEtaPt"));
   
  TFile *dataHLTEff2BinFile_pos = new TFile(dataHLTEff2BinName_pos);
  CEffUser2D dataHLTEff2Bin_pos;
  dataHLTEff2Bin_pos.loadEff((TH2D*)dataHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEff2BinFile_neg = new TFile(dataHLTEff2BinName_neg);
  CEffUser2D dataHLTEff2Bin_neg;
  dataHLTEff2Bin_neg.loadEff((TH2D*)dataHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrhEtaPt"));
    
  TFile *zeeHLTEff2BinFile_pos = new TFile(zeeHLTEff2BinName_pos);
  CEffUser2D zeeHLTEff2Bin_pos;
  zeeHLTEff2Bin_pos.loadEff((TH2D*)zeeHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zeeHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zeeHLTEff2BinFile_neg = new TFile(zeeHLTEff2BinName_neg);
  CEffUser2D zeeHLTEff2Bin_neg;
  zeeHLTEff2Bin_neg.loadEff((TH2D*)zeeHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zeeHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeHLTEff2BinFile_neg->Get("hErrhEtaPt"));

  //
  // Selection efficiency
  //
  cout << "Loading GSF+selection efficiencies..." << endl;
  
  TFile *dataGsfSelEffFile_pos = new TFile(dataGsfSelEffName_pos);
  CEffUser2D dataGsfSelEff_pos;
  dataGsfSelEff_pos.loadEff((TH2D*)dataGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataGsfSelEffFile_neg = new TFile(dataGsfSelEffName_neg);
  CEffUser2D dataGsfSelEff_neg;
  dataGsfSelEff_neg.loadEff((TH2D*)dataGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_pos = new TFile(zeeGsfSelEffName_pos);
  CEffUser2D zeeGsfSelEff_pos;
  zeeGsfSelEff_pos.loadEff((TH2D*)zeeGsfSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEffFile_neg = new TFile(zeeGsfSelEffName_neg);
  CEffUser2D zeeGsfSelEff_neg;
  zeeGsfSelEff_neg.loadEff((TH2D*)zeeGsfSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataGsfSelEff2BinFile_pos = new TFile(dataGsfSelEff2BinName_pos);
  CEffUser2D dataGsfSelEff2Bin_pos;
  dataGsfSelEff2Bin_pos.loadEff((TH2D*)dataGsfSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataGsfSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataGsfSelEff2BinFile_neg = new TFile(dataGsfSelEff2BinName_neg);
  CEffUser2D dataGsfSelEff2Bin_neg;
  dataGsfSelEff2Bin_neg.loadEff((TH2D*)dataGsfSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataGsfSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataGsfSelEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zeeGsfSelEff2BinFile_pos = new TFile(zeeGsfSelEff2BinName_pos);
  CEffUser2D zeeGsfSelEff2Bin_pos;
  zeeGsfSelEff2Bin_pos.loadEff((TH2D*)zeeGsfSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_pos->Get("hErrhEtaPt"));

  TFile *zeeGsfSelEff2BinFile_neg = new TFile(zeeGsfSelEff2BinName_neg);
  CEffUser2D zeeGsfSelEff2Bin_neg;
  zeeGsfSelEff2Bin_neg.loadEff((TH2D*)zeeGsfSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zeeGsfSelEff2BinFile_neg->Get("hErrhEtaPt"));

  // RoccoR  rc("../RochesterCorr/RoccoR2017.txt");
  
  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  cout << "Processing " << fileName.Data() << "..." << endl;
  infile = new TFile((inputDir+TString("/")+fileName).Data());    assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
    
  // Variables to get some of the branches out of the tree
  Float_t genVPt, genVPhi, genVy;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0, *genV=0, *genLep=0;
  TLorentzVector *lep_raw=0;
  Float_t pfCombIso;
  intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
  intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
  intree->SetBranchAddress("genVy",    &genVy);   // GEN W boson phi (signal MC)
  intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
  intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
  intree->SetBranchAddress("prefireWeight", &prefireWeight);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("scale1fb",     &scale1fb);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
  intree->SetBranchAddress(met_name.c_str(),      &met);       // MET
  intree->SetBranchAddress(metPhi_name.c_str(),   &metPhi);    // phi(MET)
  intree->SetBranchAddress(u1_name.c_str(),       &u1);        // parallel component of recoil
  intree->SetBranchAddress(u2_name.c_str(),       &u2);        // perpendicular component of recoil
  intree->SetBranchAddress("q",           &q);         // lepton charge
  intree->SetBranchAddress("lep",         &lep);       // lepton 4-vector
  intree->SetBranchAddress("lep_raw",         &lep_raw);       // lepton 4-vector (raw)
  intree->SetBranchAddress("genLep",      &genLep);       // lepton 4-vector
  intree->SetBranchAddress("genV",        &genV);       // lepton 4-vector
  intree->SetBranchAddress("pfCombIso",   &pfCombIso);       // lepton 4-vector

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + fileName;
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);
  // Actually just need to clone the tree: 
  TTree *outTree = intree->CloneTree(0);
  // Variables for the new branches:
  Double_t metCorrLep=0, metCorrMain=0, metCorrEta=0, metCorrStat=0, metCorrKeys=0, mtCorr=0;
  Double_t metCorrLepPhi=0, metCorrMainPhi=0, metCorrEtaPhi=0, metCorrStatPhi=0, metCorrKeysPhi=0;
  Double_t totalEvtWeight=1, effSFweight=1, relIso=0;
  
  outFile->cd();
  outTree->Branch("metCorrLep",      &metCorrLep,       "metCorrLep/d");      // corrected MET with only lepton corrections
  outTree->Branch("metCorrLepPhi",   &metCorrLepPhi,    "metCorrLepPhi/d");   // corrected METphi with only lepton corrections
  outTree->Branch("metCorrMain",     &metCorrMain,      "metCorrMain/d");     // corrected MET with main corrections
  outTree->Branch("metCorrMainPhi",  &metCorrMainPhi,   "metCorrMainPhi/d");  // corrected METphi with main corrections
  outTree->Branch("metCorrEta",      &metCorrEta,       "metCorrEta/d");      // corrected MET with eta corrections
  outTree->Branch("metCorrEtaPhi",   &metCorrEtaPhi,    "metCorrEtaPhi/d");   // corrected METphi with eta corrections
  outTree->Branch("metCorrStat",     &metCorrStat,      "metCorrStat/d");     // corrected MET with stat corrections
  outTree->Branch("metCorrStatPhi",  &metCorrStatPhi,   "metCorrStatPhi/d");  // corrected METphi with stat corrections
  outTree->Branch("metCorrKeys",     &metCorrKeys,      "metCorrKeys/d");     // corrected MET with keys corrections
  outTree->Branch("metCorrKeysPhi",  &metCorrKeysPhi,   "metCorrKeysPhi/d");  // corrected METphi with keys corrections
  outTree->Branch("relIso",          &relIso,           "relIso/d");          // scaled isolation variable that needs calculation
  outTree->Branch("mtCorr",          &mtCorr,           "mtCorr/d");          // corrected MET with keys corrections
  outTree->Branch("totalEvtWeight",  &totalEvtWeight,   "totalEvtWeight/d");  // total event weight
  outTree->Branch("effSFweight",     &effSFweight,      "effSFweight/d");     // scale factors weight
  // outTree->Branch("lep_raw",         "TLorentzVector",  &lep_raw);            // uncorrected lepton vector
  


  // Double_t mt=-999;
  
  
  UInt_t iterator=15;
  // UInt_t iterator=1;
  if(isData)iterator=1;
  //
  // loop over events
  //
  std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
  // for(UInt_t ientry=0; ientry<((int)intree->GetEntries())*0.01; ientry+=iterator) {
    intree->GetEntry(ientry);
    if(ientry%100000==0) cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;

    // vector containing raw lepton info for correcting MET
    TVector2 vLepRaw((lep_raw->Pt())*cos(lep_raw->Phi()),(lep_raw->Pt())*sin(lep_raw->Phi()));
    // (*lep_raw)=(*lep); // is this legit

    double pU1         = 0;  //--
    double pU2         = 0;  //--

    Double_t effdata =1, effmc=1;
    Double_t corr=1;

    if(fabs(lep->Eta()) > ETA_CUT) continue;
      
  // ========================================================================================
  //   Calculate the Efficiency SF Correction weight
  // ----------------------------------------------------------------------------------------
    if(q>0) { 
      effdata *= (1.-dataHLTEff_pos.getEff(lep->Eta(), lep->Pt())); 
      effmc   *= (1.-zeeHLTEff_pos.getEff(lep->Eta(), lep->Pt())); 
    } else {
      effdata *= (1.-dataHLTEff_neg.getEff(lep->Eta(), lep->Pt())); 
      effmc   *= (1.-zeeHLTEff_neg.getEff(lep->Eta(), lep->Pt())); 
    }

    effdata = 1.-effdata;
    effmc   = 1.-effmc;
    corr *= effdata/effmc;
    if(q>0) { 
      effdata *= dataGsfSelEff_pos.getEff(lep->Eta(), lep->Pt()); 
      effmc   *= zeeGsfSelEff_pos.getEff(lep->Eta(), lep->Pt());
    } else {
      effdata *= dataGsfSelEff_neg.getEff(lep->Eta(), lep->Pt()); 
      effmc   *= zeeGsfSelEff_neg.getEff(lep->Eta(), lep->Pt()); 
    }
    corr *= effdata/effmc;
    // std::cout << "got efficiencies" << std::endl;
    if(isData){
      // corrected (smear/scale) lepton for MET correction
      TVector2 vLepCor((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
      // calculate the corrected MET
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));  //move the declaration elsewhere
      Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod(); // calculate the MET corrected for lepton scale
      Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();// calculate the MET corrected for lepton scale
      mt     = sqrt( 2.0 * (lep->Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),corrMetPhiLepton))) );
      
      metCorrLep=corrMetWithLepton;
      mtCorr=mt;
    } else {
      totalEvtWeight=corr*scale1fb*prefireWeight;
      
      TVector2 vLepCor((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
      Double_t lepPt = vLepCor.Mod();
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
      Double_t corrMet = (vMetCorr + vLepRaw - vLepCor).Mod();
      Double_t corrMetPhi = (vMetCorr + vLepRaw - vLepCor).Phi();
      
      metCorrLep=corrMet; metCorrLepPhi=corrMetPhi;
        
        
    // std::cout << "start recoil" << std::endl;
      if(isRecoil) {
        if(q>0) {
          if(doKeys) {
            metCorrKeys=corrMet; metCorrKeysPhi=corrMetPhi;
            rcKeysWp->CorrectInvCdf(metCorrKeys,metCorrKeysPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kTRUE,kFALSE);
          }
          if(doEta) {
            metCorrEta=corrMet; metCorrEtaPhi=corrMetPhi;
            if(fabs(genVy)<0.5)
              rcEta05Wp->CorrectInvCdf(metCorrEta,metCorrEtaPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              rcEta051Wp->CorrectInvCdf(metCorrEta,metCorrEtaPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else
              rcEta1Wp->CorrectInvCdf(metCorrEta,metCorrEtaPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE); 
          }
          if(doInclusive){
            metCorrMain=corrMet; metCorrMainPhi=corrMetPhi;
            rcMainWp->CorrectInvCdf(metCorrMain,metCorrMainPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
          }
          if(doStat){
            metCorrStat=corrMet; metCorrStatPhi=corrMetPhi;
            rcStatW->CorrectInvCdf(metCorrStat,metCorrStatPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kTRUE);
          }

        } else {
          if(doKeys) {
            metCorrKeys=corrMet; metCorrKeysPhi=corrMetPhi;
            rcKeysWm->CorrectInvCdf(metCorrKeys,metCorrKeysPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kTRUE,kFALSE);
          }
          if(doEta) {
            metCorrEta=corrMet; metCorrEtaPhi=corrMetPhi;
            if(fabs(genVy)<0.5)
              rcEta05Wm->CorrectInvCdf(metCorrEta,metCorrEtaPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              rcEta051Wm->CorrectInvCdf(metCorrEta,metCorrEtaPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else
              rcEta1Wm->CorrectInvCdf(metCorrEta,metCorrEtaPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE); 
          }
          if(doInclusive){
            metCorrMain=corrMet; metCorrMainPhi=corrMetPhi;
            rcMainWm->CorrectInvCdf(metCorrMain,metCorrMainPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
          }
          if(doStat){
            metCorrStat=corrMet; metCorrStatPhi=corrMetPhi;
            rcStatW->CorrectInvCdf(metCorrStat,metCorrStatPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kTRUE);
          }
        }
      }
    }
    
    // std::cout << "done w recoil" << std::endl;
    // set electron relIso to be these quantities but scaled to a constant 0.15 cutoff
    // iso= 0.0287+0.506/lep->Pt()// barrel
    // iso/(tag.Pt()) >= 0.0445+0.963/(tag.Pt())
    if(fabs(lep->Eta())<ECAL_GAP_LOW)
      relIso=(pfCombIso-0.506)/lep->Pt()+(0.15-0.0287);
    if(fabs(lep->Eta())>ECAL_GAP_HIGH)
      relIso=(pfCombIso-0.963)/lep->Pt()+(0.15-0.0445);
    // relIso = pfCombIso/lep_raw->Pt(); // change this shit too
    effSFweight=corr;
    outTree->Fill(); // add new info per event to the new tree
  }//end of loop over events
  // std::cout << "end loop over events"<< std::endl;
  delete rcMainWp;
  delete rcMainWm;
  delete rcKeysWp;
  delete rcKeysWm;
  delete rcStatW;
  delete rcEta05Wp;
  delete rcEta051Wp;
  delete rcEta1Wp;
  delete rcEta05Wm;
  delete rcEta051Wm;
  delete rcEta1Wm;
  // std::cout << "clean up memory" << std::endl;
    
  outFile->cd();
  outFile->Write();
  std::cout << "wrote outfile" << std::endl;
  
  delete intree;
  delete infile;
  return;
} // end of function