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
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t mu_MASS = 0.1057;
 
 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Point to the Efficiency SF
 // ------------------------------------------------------------------------------------------------------------------------------------------
 // These are the official recoil corrections path, I will update corrections by putting new ones in this folder
  const TString baseDir = "/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV/results/Zmm/";
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



  Bool_t isData = (fileName.CompareTo("data_select.root")==0);
  std::cout << fileName.CompareTo("data_select.root",TString::kIgnoreCase) << std::endl;
  Bool_t isRecoil = (fileName.CompareTo("wm_select.raw.root")==0||fileName.CompareTo("wm0_select.raw.root")==0||fileName.CompareTo("wm1_select.raw.root")==0||fileName.CompareTo("wm2_select.raw.root")==0||fileName.CompareTo("wx_select.raw.root")==0||fileName.CompareTo("zxx_select.raw.root")==0||fileName.CompareTo("wm_select.root")==0||fileName.CompareTo("wm0_select.root")==0||fileName.CompareTo("wm1_select.root")==0||fileName.CompareTo("wm2_select.root")==0||fileName.CompareTo("wx_select.root")==0||fileName.CompareTo("zxx_select.root")==0);
  std::cout << "isData" << isData << std::endl;

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
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight;
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
  Double_t metCorrLep=0, metCorrMain=0, metCorrEta=0, metCorrStat=0, metCorrKeys=0, mtCorr=0;
  Double_t metCorrLepPhi=0, metCorrMainPhi=0, metCorrEtaPhi=0, metCorrStatPhi=0, metCorrKeysPhi=0;
  Double_t totalEvtWeight=1, effSFweight=1, relIso=0;
  TLorentzVector *lep_raw=0;
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
  outTree->Branch("lep_raw",         "TLorentzVector",  &lep_raw);            // uncorrected lepton vector
  


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
    TVector2 vLepRaw((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
    (*lep_raw)=(*lep); // is this legit

    double pU1         = 0;  //--
    double pU2         = 0;  //--

    Double_t effdata =1, effmc=1;
    Double_t corr=1;

    if(fabs(lep->Eta()) > ETA_CUT) continue;
      
  // ========================================================================================
  //   Calculate the Efficiency SF Correction weight
  // ----------------------------------------------------------------------------------------

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
  
    effdata=1; effmc=1;
    // effSigShapedata=1;
    // effBkgShapedata=1;
    if(q>0) {
      effdata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
    } else {
      effdata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
    }
    corr *= effdata/effmc;
    
    effdata=1; effmc=1;
    // effSigShapedata=1;
    // effBkgShapedata=1;
    if(q>0) { 
      effdata *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
    } else {
      effdata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
      effmc   *= zmmStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
    }
    corr *= effdata/effmc; 
    
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
      
      metCorrLep=corrMetWithLepton;
      mtCorr=mt;
    } else {
      totalEvtWeight=corr*scale1fb*prefireWeight;
        

      // Do some Rochester corrections for MC
      TLorentzVector mu1;
      mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      double mcSF1 = rc.kSpreadMC(q, mu1.Pt(), mu1.Eta(), mu1.Phi(), genLep->Pt());
      mu1*=mcSF1;
      // (*lep)*=mcSF1; // is this legit lol
      
      TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      Double_t lepPt = vLepCor.Mod();
      // std::cout << "lepPt " << mu1.Pt() << "  lep->Pt() " << lep->Pt() << "  lep_raw->Pt() " << lep_raw->Pt() << std::endl;
      // change to have rochester corrected muon and raw lepton with MET corrected same way as electron channel 
    
            // Double_t corrMet=met, corrMetPhi=metPhi;
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
      Double_t corrMet = (vMetCorr + vLepRaw - vLepCor).Mod();
      Double_t corrMetPhi = (vMetCorr + vLepRaw - vLepCor).Phi();
      // Double_t corrMet
      metCorrLep=corrMet; metCorrLepPhi=corrMetPhi;
        
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
        mtCorr  = sqrt( 2.0 * (lep->Pt()) * (metCorrMain) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metCorrMainPhi))) );
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
