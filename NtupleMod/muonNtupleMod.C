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

void muonNtupleMod(const TString  outputDir,   // output directory 
                   const TString  inputDir,    // input directory
                   const TString  sqrts,      // 13 or 5 TeV string specifier
                   const TString  fileName,    // both the input and output final file name i.e. data_select.root
                   const TString  sysFileSIT, // constains the uncertainty info for selection/id/trk efficiency
                   const TString  sysFileSta,  // contains the alternate shape info for standalone efficiencies
                   const TString  sysFileSta_Alt // alternate file to get the direct comparison for the one shape
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  // flage to control applying the recoil corrections
  bool doInclusive = true; // This should be the standard recoil correction: 3-gaussian inclusive eta
  bool doKeys = true; // RooKeysPDF instead of 3-Gaus
  bool doEta = true; // eta-binned 3-Gaus fit
  bool doStat = false; //  Statistical Uncertainty
  int nNV = 10;
  // which MET type we use
  bool doPF = true;
  
  std::string u1_name; std::string u2_name;
  std::string met_name; std::string metPhi_name;
//   std::string recoilType;


  // Control the types of uncertainties
  enum{no,cent,eta,keys,ru,rd,stat0,stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8,stat9};
  const string vMET[]={"no","cent","eta","keys","ru","rd","stat0","stat1","stat2","stat3","stat4","stat5","stat6","stat7","stat8","stat9"};
  int nMET = sizeof(vMET)/sizeof(vMET[0]);
  int ns=nMET-nNV;
  // front half should be nMET-nNV
  
  enum{main,mc,fsr,bkg,tagpt,statu,statd,pfireu,pfired};
  const string vWeight[]={"eff","mc","fsr","bkg","tagpt","statu","statd","pfireu","pfired"};
  int nWeight = sizeof(vWeight)/sizeof(vWeight[0]);

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
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t mu_MASS = 0.1057;
 
 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Point to the Efficiency SF
 // ------------------------------------------------------------------------------------------------------------------------------------------
 // These are the official recoil corrections path, I will update corrections by putting new ones in this folder
  const TString baseDir = "/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zmm/";
  // "/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV/results/Zmm/";
  // path needs to be updated ^
  const TString dataHLTEffName_pos = baseDir + "Data/MuHLTEff_aMCxPythia/Positive/eff.root";
  const TString dataHLTEffName_neg = baseDir + "Data/MuHLTEff_aMCxPythia/Negative/eff.root";
  const TString zmmHLTEffName_pos  = baseDir + "MC/MuHLTEff_aMCxPythia/Positive/eff.root";
  const TString zmmHLTEffName_neg  = baseDir + "MC/MuHLTEff_aMCxPythia/Negative/eff.root";

  const TString dataSelEffName_pos = baseDir + "Data/MuSITEff_aMCxPythia/Positive/eff.root";
  const TString dataSelEffName_neg = baseDir + "Data/MuSITEff_aMCxPythia/Negative/eff.root";
  const TString zmmSelEffName_pos  = baseDir + "MC/MuSITEff_aMCxPythia/Positive/eff.root";
  const TString zmmSelEffName_neg  = baseDir + "MC/MuSITEff_aMCxPythia/Negative/eff.root";

  const TString dataStaEffName_pos = baseDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString dataStaEffName_neg = baseDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString zmmStaEffName_pos  = baseDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString zmmStaEffName_neg  = baseDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";

  // Include the uncertainty maps for the efficiencies
  // uncertainty files
  // TFile *fSysSIT = TFile::Open(sysFileSIT);
  // TFile *fSysSta = TFile::Open(sysFileSta);
  // TFile *fSysSta2 = TFile::Open(sysFileSta_Alt);

  // TH2D  *hSysSITSigFSRNeg = (TH2D*) fSysSIT->Get("hMuSITEffSigFSRNeg");
  // TH2D  *hSysSITSigFSRPos = (TH2D*) fSysSIT->Get("hMuSITEffSigFSRPos"); 
  // TH2D  *hSysSITSigMCNeg  = (TH2D*) fSysSIT->Get("hMuSITEffSigMCNeg"); 
  // TH2D  *hSysSITSigMCPos  = (TH2D*) fSysSIT->Get("hMuSITEffSigMCPos"); 
  // TH2D  *hSysSITBkgNeg    = (TH2D*) fSysSIT->Get("hMuSITEffBkgNeg"); 
  // TH2D  *hSysSITBkgPos    = (TH2D*) fSysSIT->Get("hMuSITEffBkgPos"); 
  
  // // fsr needs to come from the 
  // TH2D  *hSysStaSigFSRNeg = (TH2D*) fSysSta2->Get("hMuStaEffSigFSRNeg");
  // TH2D  *hSysStaSigFSRPos = (TH2D*) fSysSta2->Get("hMuStaEffSigFSRPos"); 
  // TH2D  *hSysStaSigMCNeg  = (TH2D*) fSysSta->Get("hMuStaEffSigMCNeg"); 
  // TH2D  *hSysStaSigMCPos  = (TH2D*) fSysSta->Get("hMuStaEffSigMCPos"); 
  // TH2D  *hSysStaBkgNeg    = (TH2D*) fSysSta->Get("hMuStaEffBkgNeg"); 
  // TH2D  *hSysStaBkgPos    = (TH2D*) fSysSta->Get("hMuStaEffBkgPos"); 
  
  TFile *fSysSIT = new TFile(sysFileSIT);
  TFile *fSysSta = new TFile(sysFileSta);
  CEffUser2D sysSIT_FSR_Neg, sysSIT_MC_Neg, sysSIT_Bkg_Neg;
  CEffUser2D sysSIT_FSR_Pos, sysSIT_MC_Pos, sysSIT_Bkg_Pos;
  CEffUser2D sysSta_FSR_Neg, sysSta_MC_Neg, sysSta_Bkg_Neg;
  CEffUser2D sysSta_FSR_Pos, sysSta_MC_Pos, sysSta_Bkg_Pos;
  
  sysSIT_FSR_Neg.loadEff((TH2D*)fSysSIT->Get("hMuSITEffSigFSRNeg"), (TH2D*)fSysSIT->Get(""), (TH2D*)fSysSIT->Get(""));
  sysSIT_FSR_Pos.loadEff((TH2D*)fSysSIT->Get("hMuSITEffSigFSRPos"), (TH2D*)fSysSIT->Get(""), (TH2D*)fSysSIT->Get(""));
  
  sysSIT_MC_Neg.loadEff((TH2D*)fSysSIT->Get("hMuSITEffSigMCNeg"), (TH2D*)fSysSIT->Get(""), (TH2D*)fSysSIT->Get(""));
  sysSIT_MC_Pos.loadEff((TH2D*)fSysSIT->Get("hMuSITEffSigMCPos"), (TH2D*)fSysSIT->Get(""), (TH2D*)fSysSIT->Get(""));
  
  sysSIT_Bkg_Neg.loadEff((TH2D*)fSysSIT->Get("hMuSITEffBkgNeg"), (TH2D*)fSysSIT->Get(""), (TH2D*)fSysSIT->Get(""));
  sysSIT_Bkg_Pos.loadEff((TH2D*)fSysSIT->Get("hMuSITEffBkgPos"), (TH2D*)fSysSIT->Get(""), (TH2D*)fSysSIT->Get(""));
  
  sysSta_FSR_Neg.loadEff((TH2D*)fSysSta->Get("hMuStaEffSigFSRNeg"), (TH2D*)fSysSta->Get(""), (TH2D*)fSysSta->Get(""));
  sysSta_FSR_Pos.loadEff((TH2D*)fSysSta->Get("hMuStaEffSigFSRPos"), (TH2D*)fSysSta->Get(""), (TH2D*)fSysSta->Get(""));
  
  sysSta_MC_Neg.loadEff((TH2D*)fSysSta->Get("hMuStaEffSigMCNeg"), (TH2D*)fSysSta->Get(""), (TH2D*)fSysSta->Get(""));
  sysSta_MC_Pos.loadEff((TH2D*)fSysSta->Get("hMuStaEffSigMCPos"), (TH2D*)fSysSta->Get(""), (TH2D*)fSysSta->Get(""));
  
  sysSta_Bkg_Neg.loadEff((TH2D*)fSysSta->Get("hMuStaEffBkgNeg"), (TH2D*)fSysSta->Get(""), (TH2D*)fSysSta->Get(""));
  sysSta_Bkg_Pos.loadEff((TH2D*)fSysSta->Get("hMuStaEffBkgPos"), (TH2D*)fSysSta->Get(""), (TH2D*)fSysSta->Get(""));


  Bool_t isData = (fileName.CompareTo("data_select.root")==0);
  std::cout << fileName.CompareTo("data_select.root",TString::kIgnoreCase) << std::endl;
  Bool_t isRecoil = (fileName.CompareTo("wm_select.raw.root")==0||fileName.CompareTo("wm0_select.raw.root")==0||fileName.CompareTo("wm1_select.raw.root")==0||fileName.CompareTo("wm2_select.raw.root")==0||fileName.CompareTo("wx_select.raw.root")==0||fileName.CompareTo("zxx_select.raw.root")==0||fileName.CompareTo("wm_select.root")==0||fileName.CompareTo("wm0_select.root")==0||fileName.CompareTo("wm1_select.root")==0||fileName.CompareTo("wm2_select.root")==0||fileName.CompareTo("wx_select.root")==0||fileName.CompareTo("zxx_select.root")==0);
  std::cout << "isData " << isData << std::endl;
  std::cout << "isRecoil " << isRecoil << std::endl;

  if(isData) {doInclusive = false; doKeys = false; doEta = false; doStat = false;}

 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Load the Recoil Correction Files
 // ------------------------------------------------------------------------------------------------------------------------------------------
  // ===================== Recoil correction files ============================
  const TString directory("/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Recoil");
 
  // New Recoil Correctors for everything
  RecoilCorrector *rcMainWp    = new  RecoilCorrector("",""); RecoilCorrector *rcMainWm    = new  RecoilCorrector("","");
  vector<RecoilCorrector*> rcStatW;
  for(int i=0; i < nNV; ++i) {
    RecoilCorrector *tempStatW     = new  RecoilCorrector("","");
    rcStatW.push_back(tempStatW);
    
  }
  // RecoilCorrector *rcStatW     = new  RecoilCorrector("","");
  RecoilCorrector *rcKeysWp    = new  RecoilCorrector("",""); RecoilCorrector *rcKeysWm    = new  RecoilCorrector("","");
  RecoilCorrector *rcEta05Wp   = new  RecoilCorrector("",""); RecoilCorrector *rcEta05Wm   = new  RecoilCorrector("","");
  RecoilCorrector *rcEta051Wp  = new  RecoilCorrector("",""); RecoilCorrector *rcEta051Wm  = new  RecoilCorrector("","");
  RecoilCorrector *rcEta1Wp    = new  RecoilCorrector("",""); RecoilCorrector *rcEta1Wm    = new  RecoilCorrector("","");
  // also make sure to add the Wm stuff
  if(doInclusive){
    rcMainWp->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMC_PF_%s_2G/",directory.Data(),sqrts.Data()));
    rcMainWp->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G/",directory.Data(),sqrts.Data()));
    rcMainWp->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()));
    
    rcMainWm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()));
    rcMainWm->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G/",directory.Data(),sqrts.Data()));
    rcMainWm->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()));
  } 
  if (doStat){
    int rec_sig = 1;
    for(int i = 0; i < nNV; i++){
      rcStatW[i]->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()), i, rec_sig);
      rcStatW[i]->loadRooWorkspacesDiagData(Form("%s/ZmmData_PF_%s_2G/",directory.Data(),sqrts.Data()), i, rec_sig);
      rcStatW[i]->loadRooWorkspacesDiagMC(Form("%s/ZmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()), i, rec_sig);
    }
    
  } 
  if (doEta){
    rcEta05Wp->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMC_PF_%s_Eta1/",directory.Data(),sqrts.Data()));
    rcEta05Wp->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta1/",directory.Data(),sqrts.Data()));
    rcEta05Wp->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta1/",directory.Data(),sqrts.Data()));
    
    rcEta051Wp->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMC_PF_%s_Eta2/",directory.Data(),sqrts.Data()));
    rcEta051Wp->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta2/",directory.Data(),sqrts.Data()));
    rcEta051Wp->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta2/",directory.Data(),sqrts.Data()));
      
    rcEta1Wp->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMC_PF_%s_Eta3/",directory.Data(),sqrts.Data()));
    rcEta1Wp->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta3/",directory.Data(),sqrts.Data()));
    rcEta1Wp->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta3/",directory.Data(),sqrts.Data()));

    rcEta05Wm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_Eta1/",directory.Data(),sqrts.Data()));
    rcEta05Wm->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta1/",directory.Data(),sqrts.Data()));
    rcEta05Wm->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta1/",directory.Data(),sqrts.Data()));
    
    rcEta051Wm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_Eta2/",directory.Data(),sqrts.Data()));
    rcEta051Wm->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta2/",directory.Data(),sqrts.Data()));
    rcEta051Wm->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta2/",directory.Data(),sqrts.Data()));
    
    rcEta1Wm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_Eta3/",directory.Data(),sqrts.Data()));
    rcEta1Wm->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta3/",directory.Data(),sqrts.Data()));
      rcEta1Wm->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta3/",directory.Data(),sqrts.Data()));
  }
  if (doKeys){
    rcKeysWp->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMC_PF_%s_Keys/",directory.Data(),sqrts.Data()));
    rcKeysWp->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Keys/",directory.Data(),sqrts.Data()));
    rcKeysWp->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Keys/",directory.Data(),sqrts.Data()));
    
    rcKeysWm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMC_PF_%s_Keys/",directory.Data(),sqrts.Data()));
    rcKeysWm->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Keys/",directory.Data(),sqrts.Data()));
    rcKeysWm->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Keys/",directory.Data(),sqrts.Data()));
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
  cout << "Loading trigger efficiencies..." << endl;
  
  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
    
  TFile *zmmHLTEffFile_pos = new TFile(zmmHLTEffName_pos);  CEffUser2D zmmHLTEff_pos;
  zmmHLTEff_pos.loadEff((TH2D*)zmmHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEffFile_neg = new TFile(zmmHLTEffName_neg);  CEffUser2D zmmHLTEff_neg;
  zmmHLTEff_neg.loadEff((TH2D*)zmmHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrhEtaPt"));

  //
  // Selection efficiency
  //
  cout << "Loading selection efficiencies..." << endl;
  
  TFile *dataSelEffFile_pos = new TFile(dataSelEffName_pos);  CEffUser2D dataSelEff_pos;
  dataSelEff_pos.loadEff((TH2D*)dataSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEffFile_neg = new TFile(dataSelEffName_neg);  CEffUser2D dataSelEff_neg;
  dataSelEff_neg.loadEff((TH2D*)dataSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEffFile_pos = new TFile(zmmSelEffName_pos);  CEffUser2D zmmSelEff_pos;
  zmmSelEff_pos.loadEff((TH2D*)zmmSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEffFile_neg = new TFile(zmmSelEffName_neg);  CEffUser2D zmmSelEff_neg;
  zmmSelEff_neg.loadEff((TH2D*)zmmSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrhEtaPt"));

  
  //
  // Standalone efficiency
  //
  cout << "Loading standalone efficiencies..." << endl;
  
  TFile *dataStaEffFile_pos = new TFile(dataStaEffName_pos);  CEffUser2D dataStaEff_pos;
  dataStaEff_pos.loadEff((TH2D*)dataStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEffFile_neg = new TFile(dataStaEffName_neg);  CEffUser2D dataStaEff_neg;
  dataStaEff_neg.loadEff((TH2D*)dataStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_pos = new TFile(zmmStaEffName_pos);  CEffUser2D zmmStaEff_pos;
  zmmStaEff_pos.loadEff((TH2D*)zmmStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_neg = new TFile(zmmStaEffName_neg);  CEffUser2D zmmStaEff_neg;
  zmmStaEff_neg.loadEff((TH2D*)zmmStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrhEtaPt")); 

  RoccoR  rc("../RochesterCorr/RoccoR2017.txt");
  
  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  cout << "Processing " << fileName.Data() << "..." << endl;
  infile = new TFile((inputDir+TString("/")+fileName).Data());    assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
    
  // Variables to get some of the branches out of the tree
  Float_t genVPt, genVPhi, genVy;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight, prefireUp, prefireDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0, *genV=0, *genLep=0;
  Float_t pfCombIso;
  intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
  intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
  intree->SetBranchAddress("genVy",    &genVy);   // GEN W boson phi (signal MC)
  intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
  intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
  intree->SetBranchAddress("prefireWeight", &prefireWeight);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("prefireUp",   &prefireUp);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("prefireDown", &prefireDown);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("scale1fb",     &scale1fb);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("scale1fbUp",   &scale1fbUp);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
  intree->SetBranchAddress(met_name.c_str(),      &met);       // MET
  intree->SetBranchAddress(metPhi_name.c_str(),   &metPhi);    // phi(MET)
  intree->SetBranchAddress(u1_name.c_str(),       &u1);        // parallel component of recoil
  intree->SetBranchAddress(u2_name.c_str(),       &u2);        // perpendicular component of recoil
  intree->SetBranchAddress("q",           &q);         // lepton charge
  intree->SetBranchAddress("lep",         &lep);       // lepton 4-vector
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
  // Double_t metVars[no]=0, metVars[cent]=0, metVars[eta]=0,metVars[keys]=0, mtCorr=0;
  // Double_t metVarsPhi[no]=0, metVarsPhi[cent]=0, metVarsPhi[eta]=0, metVarsPhi[keys]=0;
  Double_t effSFweight=1, relIso=0;
  Double_t evtWeightSysFSR=1,evtWeightSysMC=1,evtWeightSysBkg=1;
  Double_t mtCorr=0;
  vector<Double_t>  metVars, metVarsPhi;
  vector<Double_t>  evtWeight;
  
  for(int i=0; i < nMET; i++) {metVars.push_back(0); metVarsPhi.push_back(0);}
  for(int i=0; i < nWeight; i++) evtWeight.push_back(0);
  
  TLorentzVector *lep_raw=0;
  outFile->cd();
  // outTree->Branch("metVars[no]",      &metVars[no],       "metVars[no]/d");      // corrected MET with only lepton corrections
  // outTree->Branch("metVarsPhi[no]",   &metVarsPhi[no],    "metVarsPhi[no]/d");   // corrected METphi with only lepton corrections
  // outTree->Branch("metVars[main]",     &metVars[main],      "metVars[main]/d");     // corrected MET with main corrections
  // outTree->Branch("metVarsPhi[main]",  &metVarsPhi[main],   "metVarsPhi[main]/d");  // corrected METphi with main corrections
  // outTree->Branch("metVars[eta]",      &metVars[eta],       "metVars[eta]/d");      // corrected MET with eta corrections
  // outTree->Branch("metVarsPhi[eta]",   &metVarsPhi[eta],    "metVarsPhi[eta]/d");   // corrected METphi with eta corrections
  // outTree->Branch("metCorrStat",     &metCorrStat,      "metCorrStat/d");     // corrected MET with stat corrections
  // outTree->Branch("metCorrStatPhi",  &metCorrStatPhi,   "metCorrStatPhi/d");  // corrected METphi with stat corrections
  // outTree->Branch("metVars[keys]",     &metVars[keys],      "metVars[keys]/d");     // corrected MET with keys corrections
  // outTree->Branch("metVarsPhi[keys]",  &metVarsPhi[keys],   "metVarsPhi[keys]/d");  // corrected METphi with keys corrections
  outTree->Branch("relIso",          &relIso,           "relIso/d");          // scaled isolation variable that needs calculation
  outTree->Branch("mtCorr",          &mtCorr,           "mtCorr/d");          // corrected MET with keys corrections
  // outTree->Branch("totalEvtWeight",  &totalEvtWeight,   "totalEvtWeight/d");  // total event weight
  // outTree->Branch("evtWeightSysFSR",  &evtWeightSysFSR,  "evtWeightSysFSR/d");  // total event weight
  // outTree->Branch("evtWeightSysMC",   &evtWeightSysMC,   "evtWeightSysMC/d");  // total event weight
  // outTree->Branch("evtWeightSysBkg",  &evtWeightSysBkg,  "evtWeightSysBkg/d");  // total event weight
  outTree->Branch("evtWeight",       "vector<Double_t>", &evtWeight); // event weight vector
  outTree->Branch("effSFweight",     &effSFweight,      "effSFweight/d");     // scale factors weight
  outTree->Branch("lep_raw",         "TLorentzVector",   &lep_raw);            // uncorrected lepton vector
  outTree->Branch("metVars",     "vector<Double_t>",  &metVars);            // uncorrected lepton vector
  outTree->Branch("metVarsPhi",  "vector<Double_t>",  &metVarsPhi);         // uncorrected lepton vector
  

  // Double_t mt=-999;
  
  
  UInt_t iterator=15;
  // UInt_t iterator=1;
  if(isData)iterator=1;
  //
  // loop over events
  //
  std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
   for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
  // for(UInt_t ientry=0; ientry<((int)intree->GetEntries())*0.1; ientry+=iterator) {
    intree->GetEntry(ientry);
    if(ientry%1000==0)  cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;

    // vector containing raw lepton info for correcting MET
    TVector2 vLepRaw((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
    (*lep_raw)=(*lep); // is this legit

    double pU1         = 0;  //--
    double pU2         = 0;  //--

    // data/MC scale factor corrections
    Double_t effdata, effmc;
    Double_t effdataFSR, effdataMC, effdataBkg;
    Double_t corr=1;
    Double_t corrFSR=1;
    Double_t corrMC=1;
    Double_t corrBkg=1;

    if(fabs(lep->Eta()) > ETA_CUT) continue;
      
  // ========================================================================================
  //   Calculate the Efficiency SF Correction weight
  // ----------------------------------------------------------------------------------------

    // HLT efficiency. 
    // HLT doesn't have associated systematic uncertainties
    effdata=1;effmc=1;
    if(q>0) {
      effdata *= (1.-dataHLTEff_pos.getEff((lep->Eta()), lep->Pt())); 
      effmc   *= (1.-zmmHLTEff_pos.getEff((lep->Eta()), lep->Pt())); 
    } else {
      effdata *= (1.-dataHLTEff_neg.getEff((lep->Eta()), lep->Pt())); 
      effmc   *= (1.-zmmHLTEff_neg.getEff((lep->Eta()), lep->Pt())); 
    }
    effdata = 1.-effdata;
    effmc   = 1.-effmc;
    corr *= effdata/effmc;
    corrFSR *= effdata/effmc;
    corrMC  *= effdata/effmc;
    corrBkg *= effdata/effmc;
    
    // if(lep->Pt() < 25) continue;
   
    // std::cout << q << std::endl;
    effdata=1; effmc=1;
    effdataMC=1; effdataBkg=1; effdataFSR=1; 
    if(q>0) {
      effdataFSR *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()) * sysSIT_FSR_Pos.getEff(lep->Eta(), lep->Pt());
      effdataMC *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt())  * sysSIT_MC_Pos.getEff(lep->Eta(), lep->Pt());
      effdataBkg *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()) * sysSIT_Bkg_Pos.getEff(lep->Eta(), lep->Pt());
      effdata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
    } else {
      effdataFSR *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()) * sysSIT_FSR_Neg.getEff(lep->Eta(), lep->Pt()); 
      effdataMC *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt())  * sysSIT_MC_Neg.getEff(lep->Eta(), lep->Pt()); 
      effdataBkg *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()) * sysSIT_Bkg_Neg.getEff(lep->Eta(), lep->Pt()); 
      effdata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
    }
    // std::cout << "fsr " << effdataFSR << std::endl;
    // std::cout << "-----------------------" << std::endl;
    // std::cout << "fsr " << effdataFSR << std::endl;
    // std::cout << "mc " << effdataMC << std::endl;
    // std::cout << "bkg " << effdataBkg << std::endl;
    // std::cout << lep->Pt() << "  " << lep->Eta() << std::endl;
    corr *= effdata/effmc;
    corrFSR *= effdataFSR/effmc;
    corrMC  *= effdataMC/effmc;
    corrBkg *= effdataBkg/effmc;
    
    effdata=1; effmc=1;
    effdataFSR=1; effdataMC=1; effdataBkg=1; 
    // effSigShapedata=1;
    // effBkgShapedata=1;
    if(q>0) {
      effdataFSR *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()) * sysSta_FSR_Pos.getEff(lep->Eta(), lep->Pt()); 
      effdataMC *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()) * sysSta_MC_Pos.getEff(lep->Eta(), lep->Pt()); 
      effdataBkg *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()) * sysSta_Bkg_Pos.getEff(lep->Eta(), lep->Pt()); 
      effdata *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
    } else {
      effdataFSR *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()) * sysSta_FSR_Neg.getEff(lep->Eta(), lep->Pt()); 
      effdataMC *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()) * sysSta_MC_Neg.getEff(lep->Eta(), lep->Pt()); 
      effdataBkg *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()) * sysSta_Bkg_Neg.getEff(lep->Eta(), lep->Pt()); 
      effdata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
    }
    // std::cout << "fsr " << effdataFSR << std::endl;
    // std::cout << "mc " << effdataMC << std::endl;
    // std::cout << "bkg " << effdataBkg << std::endl;
    // std::cout << lep->Pt() << "  " << lep->Eta() << std::endl;
    corr *= effdata/effmc; 
    corrFSR *= effdataFSR/effmc;
    corrMC  *= effdataMC/effmc;
    corrBkg *= effdataBkg/effmc;
    
    // std::cout << "did the efficiencies" << std::endl;

    
    if(isData){
      // Apply the Rochester Corrections to data
      TLorentzVector mu1;
      mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      double dtSF1 = rc.kScaleDT(q, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
      mu1*=dtSF1;
      (*lep)*=dtSF1; // is this legit lol
      
      if(mu1.Pt()        < PT_CUT)  continue;
      
      // std::cout << "lepPt " << mu1.Pt() << "  lep->Pt() " << lep->Pt() << "  lep_raw->Pt() " << lep_raw->Pt() << std::endl;
      // corrected (smear/scale) lepton for MET correction
      TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      // calculate the corrected MET
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));  //move the declaration elsewhere
      Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod(); // calculate the MET corrected for lepton scale
      Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();// calculate the MET corrected for lepton scale
      mt     = sqrt( 2.0 * (mu1.Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(mu1.Phi(),corrMetPhiLepton))) );
      
      metVars[no]=corrMetWithLepton;
      mtCorr=mt;
    } else {
      
      evtWeight[main]=corr*scale1fb*prefireWeight;
      evtWeight[fsr]=corrFSR*scale1fb*prefireWeight;
      evtWeight[mc] =corrMC*scale1fb*prefireWeight;
      evtWeight[bkg]=corrBkg*scale1fb*prefireWeight;
      // evtWeight[tagpt]
      // evtWeight[statu]
      // evtWeight[statd]
      evtWeight[pfireu]=corr*scale1fb*prefireUp;
      evtWeight[pfired]=corr*scale1fb*prefireDown;
      // std::cout << "storight the wweights" << std::endl;
      // std::cout << "weights and alternates" << totalEvtWeight << "  " << evtWeightSysFSR << "  " << evtWeightSysMC << "  " << evtWeightSysBkg << std::endl;
       
      // if(q<0) continue;

      // Do some Rochester corrections for MC
      TLorentzVector mu1;
      TLorentzVector mu1u;
      TLorentzVector mu1d;
      mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      // for the rochest up/down
      mu1u.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      mu1d.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      
      double mcSF1 = rc.kSpreadMC(q, mu1.Pt(), mu1.Eta(), mu1.Phi(), genLep->Pt());
      mu1*=mcSF1;
      (*lep)*=mcSF1; // is this legit lol
      
      // double deltaDtSF = rc.kScaleDTerror(Q, pt, eta, phi);
      double deltaMcSF = rc.kSpreadMCerror(q, mu1u.Pt(), mu1u.Eta(), mu1u.Phi(), genLep->Pt());
      // std::cout << deltaMcSF << std::endl; *=(1+delta)
      // double deltaMcSF = rc.kSmearMCerror(Q, pt, eta, phi, nl, u);
      
      TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      TVector2 vLepCorU((mu1u.Pt())*cos(mu1u.Phi()),(mu1u.Pt())*sin(mu1u.Phi()));
      TVector2 vLepCorD((mu1d.Pt())*cos(mu1d.Phi()),(mu1d.Pt())*sin(mu1d.Phi()));
      Double_t lepPt = vLepCor.Mod();
      // std::cout << "lepPt " << mu1.Pt() << "  lep->Pt() " << lep->Pt() << "  lep_raw->Pt() " << lep_raw->Pt() << std::endl;
      // change to have rochester corrected muon and raw lepton with MET corrected same way as electron channel 
    
            // Double_t corrMet=met, corrMetPhi=metPhi;
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
      Double_t corrMet = (vMetCorr + vLepRaw - vLepCor).Mod();
      Double_t corrMetPhi = (vMetCorr + vLepRaw - vLepCor).Phi();
      
      Double_t corrMetU = (vMetCorr + vLepRaw - vLepCorU).Mod();
      Double_t corrMetPhiU = (vMetCorr + vLepRaw - vLepCorU).Phi();
      
      Double_t corrMetD = (vMetCorr + vLepRaw - vLepCorD).Mod();
      Double_t corrMetPhiD = (vMetCorr + vLepRaw - vLepCorD).Phi();
      // Double_t corrMet
      metVars[no]=corrMet; metVarsPhi[no]=corrMetPhi;
      // std::cout << "blah " << std::endl;
        
      if(isRecoil) {
        
        if(q>0) {
          if(doKeys) {
            // std::cout << "keys" << std::endl;
            metVars[keys]=corrMet; metVarsPhi[keys]=corrMetPhi;
            rcKeysWp->CorrectInvCdf(metVars[keys],metVarsPhi[keys],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kTRUE,kFALSE);
          
            // std::cout << "keys done" << std::endl;
          }
          if(doEta) {
            // std::cout << "eta" << std::endl;
            metVars[eta]=corrMet; metVarsPhi[eta]=corrMetPhi;
            if(fabs(genVy)<0.5)
              rcEta05Wp->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              rcEta051Wp->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else
              rcEta1Wp->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE); 
            // std::cout << "eta done" << std::endl;
          }
          if(doInclusive){
            // std::cout << "inclusive" << std::endl;
            metVars[cent]=corrMet; metVarsPhi[cent]=corrMetPhi;
            rcMainWp->CorrectInvCdf(metVars[cent],metVarsPhi[cent],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            metVars[ru]=corrMetU; metVarsPhi[ru]=corrMetPhiU;
            rcMainWp->CorrectInvCdf(metVars[ru],metVarsPhi[ru],genVPt,genVPhi,mu1u.Pt(),mu1u.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            metVars[rd]=corrMetD; metVarsPhi[rd]=corrMetPhiD;
            rcMainWp->CorrectInvCdf(metVars[rd],metVarsPhi[rd],genVPt,genVPhi,mu1d.Pt(),mu1d.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
          
            // std::cout << "inclusive done" << std::endl;
          }
          if(doStat){
            for(int i = 0; i < nNV; i++){
              int ofs=i+ns;
              // std::cout << "filling stat " << i << std::endl;
              metVars[ofs]=corrMet; metVarsPhi[ofs]=corrMetPhi;
              rcStatW[i]->CorrectInvCdf(metVars[ofs],metVarsPhi[ofs],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kTRUE);
            }
            // std::cout << "stat done " << std::endl; 
          }

        } else {
          if(doKeys) {
            metVars[keys]=corrMet; metVarsPhi[keys]=corrMetPhi;
            rcKeysWm->CorrectInvCdf(metVars[keys],metVarsPhi[keys],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kTRUE,kFALSE);
          }
          if(doEta) {
            metVars[eta]=corrMet; metVarsPhi[eta]=corrMetPhi;
            if(fabs(genVy)<0.5)
              rcEta05Wm->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              rcEta051Wm->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else
              rcEta1Wm->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE); 
          }
          if(doInclusive){
            metVars[cent]=corrMet; metVarsPhi[cent]=corrMetPhi;
            rcMainWm->CorrectInvCdf(metVars[cent],metVarsPhi[cent],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            metVars[ru]=corrMetU; metVarsPhi[ru]=corrMetPhiU;
            rcMainWm->CorrectInvCdf(metVars[ru],metVarsPhi[ru],genVPt,genVPhi,mu1u.Pt(),mu1u.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            metVars[rd]=corrMetD; metVarsPhi[rd]=corrMetPhiD;
            rcMainWm->CorrectInvCdf(metVars[rd],metVarsPhi[rd],genVPt,genVPhi,mu1d.Pt(),mu1d.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);            
          }
          if(doStat){
            for(int i =0; i < nNV; i++){
              int ofs=i+ns;
              metVars[ofs]=corrMet; metVarsPhi[ofs]=corrMetPhi;
              rcStatW[i]->CorrectInvCdf(metVars[ofs],metVarsPhi[ofs],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kTRUE);
            }
          }
        }
        mtCorr  = sqrt( 2.0 * (lep->Pt()) * (metVars[cent]) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metVarsPhi[cent]))) );
      }
    }
    relIso = pfCombIso/lep_raw->Pt(); 
    effSFweight=corr;
    outTree->Fill(); // add new info per event to the new tree
  }//end of loop over events
  // std::cout << "end loop over events"<< std::endl;
  delete rcMainWp;
  delete rcMainWm;
  delete rcKeysWp;
  delete rcKeysWm;
  for(int i = 0; i < nNV; i ++)delete rcStatW[i];
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
