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


//helper class to handle rochester corrections
#include <../RochesterCorr/RoccoR.cc>

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

#endif

//=== MAIN MACRO ================================================================================================= 

void fitWm_lumis_2d(const TString  outputDir,   // output directory 
           const Double_t lumi,        // integrated luminosity (/fb)'
           const Double_t lumi2, // lumi for the anti-isolation trigger
       const Double_t nsigma=0,     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
       const TString input_section = "1"
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  bool doInclusive = true; // This should be the standard recoil correction: 3-gaussian
  bool doKeys = false; // RooKeysPDF instead of 3-Gaus
  bool doEta = false; // eta-binned 3-Gaus fit
  bool doStat = false; // 
  
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

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);

 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Load the Recoil Correction Files
 // ------------------------------------------------------------------------------------------------------------------------------------------
  // ===================== Recoil correction files ============================
  const TString directory("/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Recoil");

  // for PF, 13 TeV Low PU, inclusive
  RecoilCorrector *recoilCorr = new  RecoilCorrector("","");
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("","");
  
     // also make sure to add the Wm stuff
  if(doInclusive && !doStat){
  
    // // PF MET
  recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
  recoilCorr->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()));
  recoilCorr->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
  
    // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_bigBin_v1/",directory2.Data()));
    // recoilCorr->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_bigBin_v1/",directory2.Data()));
    // recoilCorr->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_bigBin_v1/",directory2.Data()));
  
  } else if (doInclusive && doStat){
      recoilCorr->loadRooWorkspacesDiagMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()), rec_sig);
      recoilCorr->loadRooWorkspacesDiagData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()), rec_sig);
     recoilCorr->loadRooWorkspacesDiagMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()), rec_sig);
     
           // recoilCorr->loadRooWorkspacesDiagMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_bigBin_v1/",directory2.Data()), rec_sig);
      // recoilCorr->loadRooWorkspacesDiagData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_bigBin_v1/",directory2.Data()), rec_sig);
     // recoilCorr->loadRooWorkspacesDiagMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_bigBin_v1/",directory2.Data()), rec_sig);
      
  } else if (doEta){
      recoilCorr05->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_0-05_v1/",directory.Data()));
      recoilCorr05->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_0-05_v1/",directory.Data()));
      recoilCorr05->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_0-05_v1/",directory.Data()));
      
      recoilCorr051->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_05-10_v1/",directory.Data()));
      recoilCorr051->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_05-10_v1/",directory.Data()));
      recoilCorr051->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_05-10_v1/",directory.Data()));
      
      recoilCorr1->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_10-24_v1/",directory.Data()));
      recoilCorr1->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF_10-24_v1/",directory.Data()));
      recoilCorr1->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_10-24_v1/",directory.Data()));
  
  } else if (doKeys){
       
  recoilCorrKeys->loadRooWorkspacesMCtoCorrectKeys(Form("%s/LowPU2017ID_13TeV_ZmmMCPF_keys_v1/",directory.Data()));
  recoilCorrKeys->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory.Data()));
  recoilCorrKeys->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory.Data()));
  }
  
  
   // *** THis part also needs to be fixed
  // input ntuple file names
  //
  enum { eData, eWmunu, eEWK, eBKG, eZxx, eWx, eTtb, eDib, eQCD, eAntiData, eAntiWmunu, eAntiEWK, eAntiQCD, eAntiTtb, eAntiDib, eAntiWx, eAntiZxx };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/data_select.root"); typev.push_back(eData);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/wm_select.raw.root");   typev.push_back(eWmunu);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/wx_select.raw.root");  typev.push_back(eWx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/zxx_select.raw.root");  typev.push_back(eZxx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/zz_select.raw.root");  typev.push_back(eDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/ww_select.raw.root");  typev.push_back(eDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/wz_select.raw.root");  typev.push_back(eDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/Wmunu/ntuples/top_select.raw.root");  typev.push_back(eTtb);

  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/data_select.root"); typev.push_back(eAntiData);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/wx_select.root"); typev.push_back(eAntiWx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/zxx_select.root"); typev.push_back(eAntiZxx);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/ww_select.root"); typev.push_back(eAntiDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/wz_select.root"); typev.push_back(eAntiDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/zz_select.root"); typev.push_back(eAntiDib);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/wm_select.root"); typev.push_back(eAntiWmunu);
  fnamev.push_back("/afs/cern.ch/work/s/sabrandt/public/LowPU2017ID_13TeV/AntiWmunu/ntuples/top_select.root");  typev.push_back(eAntiTtb);


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

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    // Actually just need to clone the tree: 
    TTree *newtree = intree->CloneTree(0);
    // Variables for the new branches:
    Double_t metCorrMain, metCorrEta, metCorrStat, metCorrKeys, mtCorr;
    Double_t totalEvtWeight, effSFweight;
    TLorentzVector *lep_raw=0;
    
    outTree->Branch("metCorrMain",     &metCorrMain,      "metCorrMain/d");      // corrected MET with main corrections
    outTree->Branch("metCorrEta",      &metCorrEta,       "metCorrEta/d");      // corrected MET with eta corrections
    outTree->Branch("metCorrStat",     &metCorrStat,      "metCorrStat/d");      // corrected MET with stat corrections
    outTree->Branch("metCorrKeys",     &metCorrKeys,      "metCorrKeys/d");      // corrected MET with keys corrections
    outTree->Branch("mtCorr",          &mtCorr,           "mtCorr/d");      // corrected MET with keys corrections
    outTree->Branch("totalEvtWeight",  &totalEvtWeight,   "totalEvtWeight/d");      // total event weight
    outTree->Branch("effSFweight",     &effSFweight,      "effSFweight/d");      // scale factors weight
    outTree->Branch("lep_raw",         "TLorentzVector",  &lep_raw);      // corrected MET with main corrections
    
    
    // Variables to get some of the branches out of the tree
    // UInt_t  runNum, lumiSec, evtNum;
    // UInt_t  npv, npu;
    Float_t genVPt, genVPhi, genVy;
    Float_t genLepPt, genLepPhi;
    Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight;
    Float_t met, metPhi, sumEt, mt, u1, u2;
    Int_t   q;
    TLorentzVector *lep=0;
    Float_t pfChIso, pfGamIso, pfNeuIso;
    Float_t pfCombIso;
    // intree->SetBranchAddress("runNum",   &runNum);    // event run number
    // intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    // intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    // intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    // intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",  &genVy);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireWeight", &prefireWeight);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress(met_name.c_str(),      &met);       // MET
    intree->SetBranchAddress(metPhi_name.c_str(),   &metPhi);    // phi(MET)
    // intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    // intree->SetBranchAddress("mt",       &mt);        // transverse mass
    intree->SetBranchAddress(u1_name.c_str(),       &u1);        // parallel component of recoil
    intree->SetBranchAddress(u2_name.c_str(),       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep_raw_in);       // lepton 4-vector
    intree->SetBranchAddress("genLep",      &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV",     &genV);       // lepton 4-vector
    // intree->SetBranchAddress("pfChIso",  &pfChIso);
    // intree->SetBranchAddress("pfGamIso", &pfGamIso);
    // intree->SetBranchAddress("pfNeuIso", &pfNeuIso);
    intree->SetBranchAddress("pfCombIso",      &pfCombIso);       // lepton 4-vector
  
    // Double_t mt=-999;
    
    
    UInt_t iterator=15;
    // UInt_t iterator=1;
    if(typev[ifile]==eData||typev[ifile]==eAntiData)iterator=1;
    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    // for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; ientry<((int)intree->GetEntries()); ientry+=iterator) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " percent done with this file." << endl;

      // vector containing raw lepton info for correcting MET
      TVector2 vLepRaw((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));

      double pU1         = 0;  //--
      double pU2         = 0;  //--

      Double_t effdata, effmc;
      Double_t corr=1;
      Double_t eff2Bindata, eff2Binmc;
      Double_t corr2Bin=1;
      Double_t corrUp=1;
      Double_t corrDown=1;
      Double_t effSigShapedata;
      Double_t corrSigShape=1;
      Double_t effBkgShapedata;
      Double_t corrBkgShape=1;

      if(fabs(lep->Eta()) > ETA_CUT) continue;
      
        
    // ========================================================================================
    //   Calculate the Efficiency SF Correction weight
    // ----------------------------------------------------------------------------------------
  
      effdata=1; effmc=1;    
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
      effSigShapedata=1;
      effBkgShapedata=1;
      if(q>0) {
        effdata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
      } else {
        effdata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
      }
      corr *= effdata/effmc;
      
      effdata=1; effmc=1;
      effSigShapedata=1;
      effBkgShapedata=1;
      if(q>0) { 
        effdata *= dataStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmStaEff_pos.getEff((lep->Eta()), lep->Pt()); 
      } else {
        effdata *= dataStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmStaEff_neg.getEff((lep->Eta()), lep->Pt()); 
      }
      corr *= effdata/effmc; 
      
      

          
      if(typev[ifile]==eData || typev[ifile]==eAntiData){
        // Apply the Rochester Corrections to data
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        double dtSF1 = rc.kScaleDT(q, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
        mu1*=dtSF1;
        
        if(mu1.Pt()        < PT_CUT)  continue;
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
        // calculate the corrected MET
        TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));  //move the declaration elsewhere
        Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod(); // calculate the MET corrected for lepton scale
        Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();// calculate the MET corrected for lepton scale
        mt     = sqrt( 2.0 * (mu1.Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(mu1.Phi(),corrMetPhiLepton))) );

      } else {
        Double_t weight = 1;Double_t weightUp = 1;Double_t weightDown = 1;
        Double_t weight2 =1;
        //corr = 1.0;

        weight *= scale1fb*lumi*corr*prefireWeight*iterator;

        // Do some Rochester corrections for MC
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        double mcSF1 = rc.kSpreadMC(q, mu1.Pt(), mu1.Eta(), mu1.Phi(), genLep->Pt());
        mu1*=mcSF1;
         // std::cout << "pt corr " << mu1.Pt() << "  no corr "  << lep->Pt() << std::endl;
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
		// vLepCor+=vLepCor*eleCorr;
        Double_t lepPt = vLepCor.Mod();
        // change to have rochester corrected muon and raw lepton with MET corrected same way as electron channel 
        if(typev[ifile]==eWmunu || typev[ifile]==eWx || typev[ifile]==eZxx) {
          Double_t corrMet=met, corrMetPhi=metPhi;
            
            corrMet=met, corrMetPhi=metPhi;
                          // Compute the corrected MET value
              TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
              Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
              Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
            if(q>0) {
	      if(doKeys) {
                  recoilCorrKeys->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doStat);
          } else if(doEta) {
                if(fabs(dilep->Eta())<0.5)
                  recoilCorr05->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doStat);
                else if (fabs(dilep->Eta())>=0.5 && fabs(dilep->Eta())<1.0)
                  recoilCorr051->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doStat);
                else
                  recoilCorr1->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doStat); 
          }else if(doInclusive){

	      recoilCorr->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,vDilepCor.Mod(),vDilepCor.Phi(),vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doStat);
	      // // recoilCorr->CorrectInvCdf(corrMetWithLepton,corrMetPhiLepton,genVPt,genVPhi,vDilepCor.Mod(),vDilepCor.Phi(),pU1,pU2,0,0,0,doKeys,doStat);
	      }


              corrMet=met, corrMetPhi=metPhi;
            } else {
              hMuonEtaMCm->Fill(fabs(mu1.Eta()),weight);
              if(doInclusive) recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doStat);
              // else if(doPF) recoilCorrPFm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doStat);
//               else if(pileupUp)    recoilCorrPuUpm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
//               else if(pileupDown)  recoilCorrPuDownm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
              // else if(doKeys){
                  // recoilCorrKeysm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
// //                 if(fabs(genVy)<0.5)
// //                   recoilCorrKeysm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //                 else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
// //                   recoilCorrKeysm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //                 else
// //                   recoilCorrKeysm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys); 
              // } else if(doEta) {
                // if(fabs(genVy)<0.5)
                  // recoilCorrm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                // else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                  // recoilCorrm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                // else
                  // recoilCorrm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
              // }
              // Compute the corrected MET value
              TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
              Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
              Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();

              corrMet=met, corrMetPhi=metPhi;
            }
            corrMet=met, corrMetPhi=metPhi;
          }
      }
      newTree->Fill(); // add new info per event to the new tree
    }//end of loop over events
    
  }// End of loop over files
  delete infile;
  infile=0, intree=0;   
} // end of function