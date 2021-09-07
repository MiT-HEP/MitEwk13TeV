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

#include "TRandom.h"
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/RecoilCorrector.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
#include "../Utils/CEffUser2D.hh"

//helper class to handle rochester corrections
#include "../RochesterCorr/RoccoR.cc"
#include "../Utils/AppEffSF.cc"



#endif

//=== MAIN MACRO ================================================================================================= 

void muonNtupleMod(const TString  outputDir,   // output directory 
                   const TString  inputDir,    // input directory
                   const TString  sqrts,      // 13 or 5 TeV string specifier
                   const TString  fileName,    // both the input and output final file name i.e. data_select.root
                   const TString  sysFileSIT, // constains the uncertainty info for selection/id/trk efficiency
                   const TString  sysFileSta  // contains the alternate shape info for standalone efficiencies
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
  
  enum{main,mc,fsr,bkg,tagpt,effstat,pfireu,pfired};
  const string vWeight[]={"eff","mc","fsr","bkg","tagpt","effstat","pfireu","pfired"};
  int nWeight = sizeof(vWeight)/sizeof(vWeight[0]);

  if(doPF){
    u1_name = "u1";
    u2_name = "u2";
    met_name = "met";
    metPhi_name = "metPhi";
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
 
 // -----------------------------------------------------------
 //   Point to the Efficiency SF
 // -----------------------------------------------------------
  TString effDir = "/afs/cern.ch/work/a/arapyan/Run2/test/CMSSW_10_2_13/src/MitEwk13TeV/data/Efficiency/lowpu_"+sqrts+"/results/Zmm/";
  std::cout << effDir << std::endl;
  AppEffSF effs(effDir);
  effs.loadHLT("MuHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("MuSITEff_aMCxPythia","Combined","Combined");
  effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  
  effs.loadUncSel(sysFileSIT);
  effs.loadUncSta(sysFileSta);
  TH2D *hErr  = new TH2D("hErr", "",10,0,10,20,0,20);


  Bool_t isData = (fileName.CompareTo("data_select.root")==0);
  std::cout << fileName.CompareTo("data_select.root",TString::kIgnoreCase) << std::endl;
  
  Bool_t isRecoil = (fileName.CompareTo("wm_select.raw.root" )==0||
                     fileName.CompareTo("wm0_select.raw.root")==0||
                     fileName.CompareTo("wm1_select.raw.root")==0||
                     fileName.CompareTo("wm2_select.raw.root")==0||
                     fileName.CompareTo("wx_select.raw.root" )==0||
                     fileName.CompareTo("wx0_select.raw.root")==0||
                     fileName.CompareTo("wx1_select.raw.root")==0||
                     fileName.CompareTo("wx2_select.raw.root")==0||
                     fileName.CompareTo("zxx_select.raw.root")==0||
                     fileName.CompareTo("wm_select.root"     )==0||
                     fileName.CompareTo("wm0_select.root"    )==0||
                     fileName.CompareTo("wm1_select.root"    )==0||
                     fileName.CompareTo("wm2_select.root"    )==0||
                     fileName.CompareTo("wx_select.root"     )==0||
                     fileName.CompareTo("wx0_select.root"    )==0||
                     fileName.CompareTo("wx1_select.root"    )==0||
                     fileName.CompareTo("wx2_select.root"    )==0||
                     fileName.CompareTo("zxx_select.root"    )==0);
                     
  if(inputDir.Contains("Anti") && isRecoil) {
    doInclusive = true; 
    doKeys      = false; 
    doEta       = false; 
    doStat      = false;
  }
  
  std::cout << "isData " << isData << std::endl;
  std::cout << "isRecoil " << isRecoil << std::endl;

  if(isData||(!isRecoil)) {
    doInclusive = false; 
    doKeys      = false; 
    doEta       = false; 
    doStat      = false;
  }

 // ------------------------------------------------------------------------------------------------------------------------------------------
 //   Load the Recoil Correction Files
 // ------------------------------------------------------------------------------------------------------------------------------------------
  // ===================== Recoil correction files ============================
  const TString directory("/afs/cern.ch/work/a/arapyan/Run2/test/CMSSW_10_2_13/src/MitEwk13TeV/data/Recoil");
 
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
    rcMainWp->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G_bkg_fixRoch/",directory.Data(),sqrts.Data()));
    rcMainWp->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()));
    
    rcMainWm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_2G/",directory.Data(),sqrts.Data()));
    rcMainWm->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G_bkg_fixRoch/",directory.Data(),sqrts.Data()));
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
  
  RoccoR  rc("/afs/cern.ch/work/a/arapyan/Run2/test/CMSSW_10_2_13/src/MitEwk13TeV/RochesterCorr/RoccoR2017.txt");
  
  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  cout << "Processing " << fileName.Data() << "..." << endl;
  infile = new TFile((inputDir+TString("/")+fileName).Data());    assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
    
  TH1D* hGenWeights = (TH1D*)infile->Get("hGenWeights");
  // Variables to get some of the branches out of the tree
  Float_t genVPt, genVPhi, genVy;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight, prefireUp, prefireDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  UInt_t nTkLayers;
  TLorentzVector *lep=0, *genV=0, *genLep=0;
  Float_t genMuonPt=0;
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
  intree->SetBranchAddress("nTkLayers",   &nTkLayers);       // lepton 4-vector
  intree->SetBranchAddress("genMuonPt",   &genMuonPt);       // lepton 4-vector

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + fileName;
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);
  // Actually just need to clone the tree: 
  TTree *outTree = intree->CloneTree(0);
  Double_t effSFweight=1, relIso=0;
  Double_t evtWeightSysFSR=1,evtWeightSysMC=1,evtWeightSysBkg=1;
  Double_t mtCorr=0;
  vector<Double_t>  metVars, metVarsPhi;
  vector<Double_t>  evtWeight;
  
  for(int i=0; i < nMET; i++) {metVars.push_back(0); metVarsPhi.push_back(0);}
  for(int i=0; i < nWeight; i++) evtWeight.push_back(0);
  
  TLorentzVector *lep_raw=0;
  outFile->cd();
  outTree->Branch("relIso",          &relIso,           "relIso/d");          // scaled isolation variable that needs calculation
  outTree->Branch("mtCorr",          &mtCorr,           "mtCorr/d");          // corrected MET with keys corrections
  outTree->Branch("evtWeight",       "vector<Double_t>", &evtWeight); // event weight vector
  outTree->Branch("effSFweight",     &effSFweight,      "effSFweight/d");     // scale factors weight
  outTree->Branch("lep_raw",         "TLorentzVector",   &lep_raw);            // uncorrected lepton vector
  outTree->Branch("metVars",     "vector<Double_t>",  &metVars);            // uncorrected lepton vector
  outTree->Branch("metVarsPhi",  "vector<Double_t>",  &metVarsPhi);         // uncorrected lepton vector
  
  //
  // loop over events
  //
  std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
   for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) { 
    intree->GetEntry(ientry);
    if(ientry%10000==0)  cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;

    // vector containing raw lepton info for correcting MET
    TVector2 vLepRaw((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
    (*lep_raw)=(*lep); // is this legit

    double pU1         = 0;  //--
    double pU2         = 0;  //--

    // data/MC scale factor corrections
    Double_t effdata, effmc;
    Double_t effdatah, effmch, effdatal, effmcl;
    Double_t effdataFSR, effdataMC, effdataBkg;
    Double_t edTag, emTag;    
    Double_t corr=1, corrdu=1, corrdd=1, corrmu=1, corrmd=1;
    Double_t corrFSR=1;
    Double_t corrMC=1;
    Double_t corrBkg=1;
    Double_t corrTag=1;


    if(fabs(lep->Eta()) > ETA_CUT) continue;
      
  // ========================================================================================
  //   Calculate the Efficiency SF Correction weight
  // ----------------------------------------------------------------------------------------

    // HLT efficiency. 
    // HLT doesn't have associated systematic uncertainties
    effdata=1;effmc=1;
    effdatah=1;effdatal=1; effmch=1;effmcl=1;
    
    if(isData){
      // Apply the Rochester Corrections to data
      TLorentzVector mu1;
      TLorentzVector mu1u;
      TLorentzVector mu1d;
      mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      mu1u.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      mu1d.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      double dtSF1 = rc.kScaleDT(q, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
      mu1*=dtSF1;
      (*lep)*=dtSF1; // is this legit lol
      
      if(mu1.Pt()        < PT_CUT)  continue;
      
      // corrected (smear/scale) lepton for MET correction
      TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      TVector2 vLepCoru((mu1u.Pt())*cos(mu1u.Phi()),(mu1u.Pt())*sin(mu1u.Phi()));
      TVector2 vLepCord((mu1d.Pt())*cos(mu1d.Phi()),(mu1d.Pt())*sin(mu1d.Phi()));
      // calculate the corrected MET
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));  //move the declaration elsewhere
      Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod(); // calculate the MET corrected for lepton scale
      Double_t corrMetPhiLepton  = (vMetCorr + vLepRaw - vLepCor).Phi();// calculate the MET corrected for lepton scale
      mt     = sqrt( 2.0 * (mu1.Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(mu1.Phi(),corrMetPhiLepton))) );
      
      Double_t corrMetWithLeptonu = (vMetCorr + vLepRaw - vLepCoru).Mod(); // calculate the MET corrected for lepton scale
      Double_t corrMetPhiLeptonu  = (vMetCorr + vLepRaw - vLepCoru).Phi();// calculate the MET corrected for lepton scale
      Double_t corrMetWithLeptond = (vMetCorr + vLepRaw - vLepCord).Mod(); // calculate the MET corrected for lepton scale
      Double_t corrMetPhiLeptond  = (vMetCorr + vLepRaw - vLepCord).Phi();// calculate the MET corrected for lepton scale
      
      metVars[no]=corrMetWithLepton;
      metVars[ru]=corrMetWithLeptonu;
      metVars[rd]=corrMetWithLeptond;
      metVarsPhi[no]=corrMetPhiLepton;
      metVarsPhi[ru]=corrMetPhiLeptonu;
      metVarsPhi[rd]=corrMetPhiLeptond;
      
      mtCorr=mt;
    } else {
      
      // Do some Rochester corrections for MC
      TLorentzVector mu1;
      TLorentzVector mu1u;
      TLorentzVector mu1d;
      mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      // for the rochest up/down
      mu1u.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      mu1d.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
      
    
      corr = effs.fullEfficiencies(&mu1,q);
      vector<double> uncs_sta = effs.getUncSta(&mu1,q);
      vector<double> uncs_sit = effs.getUncSel(&mu1,q);
      
      corrFSR *= uncs_sta[0]*uncs_sit[0]*effs.computeHLTSF(&mu1,q); // alternate fsr model
      corrMC  *= uncs_sta[1]*uncs_sit[1]*effs.computeHLTSF(&mu1,q); // alternate mc gen model
      corrBkg *= uncs_sta[2]*uncs_sit[2]*effs.computeHLTSF(&mu1,q); // alternate bkg model
      corrTag *= uncs_sta[3]*uncs_sit[3]*effs.computeHLTSF(&mu1,q); // alternate bkg model

      double var=0.;        
      // var += effs.statUncSta(&l1, q) + effs.statUncSta(&l2, q2);
      var += effs.statUncSta(&mu1, q, hErr, hErr, 1.0);
      var += effs.statUncSel(&mu1, q, hErr, hErr, 1.0);
      var += effs.statUncHLT(&mu1, q, hErr, hErr, 1.0);
    
      evtWeight[main]    =corr    *scale1fb*prefireWeight;
      evtWeight[fsr]     =corrFSR *scale1fb*prefireWeight;
      evtWeight[mc]      =corrMC  *scale1fb*prefireWeight;
      evtWeight[bkg]     =corrBkg *scale1fb*prefireWeight;
      evtWeight[tagpt]   =corrTag *scale1fb*prefireWeight;
      evtWeight[effstat] =var     *scale1fb*prefireWeight*scale1fb*prefireWeight;
      evtWeight[pfireu]  =corr    *scale1fb*prefireUp;
      evtWeight[pfired]  =corr    *scale1fb*prefireDown;

      double rand = gRandom->Uniform(1);
      double mcSF1 = 1;
      if(genMuonPt > 0){
        mcSF1 = rc.kSpreadMC(q, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt);
      } else {
        mcSF1 = rc.kSmearMC(q, mu1.Pt(),  mu1.Eta(), mu1.Phi(), nTkLayers, rand);
      }
      mu1*=mcSF1;
      (*lep)*=mcSF1; 
      
      double deltaMcSF = 1;
      if(genMuonPt>0){
        deltaMcSF = rc.kSpreadMCerror(q, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt);
      } else {
        deltaMcSF = rc.kSmearMCerror(q, mu1.Pt(),  mu1.Eta(), mu1.Phi(), nTkLayers, rand);
      }
      mu1u*=(1+deltaMcSF);
      mu1d*=(1-deltaMcSF);
      
      TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      TVector2 vLepCorU((mu1u.Pt())*cos(mu1u.Phi()),(mu1u.Pt())*sin(mu1u.Phi()));
      TVector2 vLepCorD((mu1d.Pt())*cos(mu1d.Phi()),(mu1d.Pt())*sin(mu1d.Phi()));
      Double_t lepPt = vLepCor.Mod();
    
      TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
      Double_t corrMet    = (vMetCorr + vLepRaw - vLepCor).Mod();
      Double_t corrMetPhi = (vMetCorr + vLepRaw - vLepCor).Phi();
      
      Double_t corrMetU    = (vMetCorr + vLepRaw - vLepCorU).Mod();
      Double_t corrMetPhiU = (vMetCorr + vLepRaw - vLepCorU).Phi();
      
      Double_t corrMetD    = (vMetCorr + vLepRaw - vLepCorD).Mod();
      Double_t corrMetPhiD = (vMetCorr + vLepRaw - vLepCorD).Phi();
      
      metVars[no]   =corrMet;  metVarsPhi[no]   =corrMetPhi;
      metVars[keys] =corrMet;  metVarsPhi[keys] =corrMetPhi;
      metVars[eta]  =corrMet;  metVarsPhi[eta]  =corrMetPhi;
      metVars[cent] =corrMet;  metVarsPhi[cent] =corrMetPhi;
      metVars[ru]   =corrMetU; metVarsPhi[ru]   =corrMetPhiU;
      metVars[rd]   =corrMetD; metVarsPhi[rd]   =corrMetPhiD;
      for(int i = 0; i < nNV; i++){
        int ofs=i+ns;
        metVars[ofs]=corrMet; metVarsPhi[ofs]=corrMetPhi;
      }
  
      if(isRecoil) {
        
        if(q>0) {
          if(doKeys) {
            rcKeysWp->CorrectInvCdf(metVars[keys],metVarsPhi[keys],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kTRUE,kFALSE);
          }
          if(doEta) {
            if(fabs(genVy)<0.5)
              rcEta05Wp->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              rcEta051Wp->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else
              rcEta1Wp->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE); 
          }
          if(doInclusive){
            rcMainWp->CorrectInvCdf(metVars[cent],metVarsPhi[cent],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            rcMainWp->CorrectInvCdf(metVars[ru],  metVarsPhi[ru],  genVPt,genVPhi,mu1u.Pt(),mu1u.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            rcMainWp->CorrectInvCdf(metVars[rd],  metVarsPhi[rd],  genVPt,genVPhi,mu1d.Pt(),mu1d.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
          }
          if(doStat){
            for(int i = 0; i < nNV; i++){
              int ofs=i+ns;
              rcStatW[i]->CorrectInvCdf(metVars[ofs],metVarsPhi[ofs],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kTRUE);
            }
          }

        } else {
          if(doKeys) {
            rcKeysWm->CorrectInvCdf(metVars[keys],metVarsPhi[keys],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kTRUE,kFALSE);
          }
          if(doEta) {
            if(fabs(genVy)<0.5)
              rcEta05Wm->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
              rcEta051Wm->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            else
              rcEta1Wm->CorrectInvCdf(metVars[eta],metVarsPhi[eta],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE); 
          }
          if(doInclusive){
            rcMainWm->CorrectInvCdf(metVars[cent],metVarsPhi[cent],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            rcMainWm->CorrectInvCdf(metVars[ru],  metVarsPhi[ru],  genVPt,genVPhi,mu1u.Pt(),mu1u.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);
            rcMainWm->CorrectInvCdf(metVars[rd],  metVarsPhi[rd],  genVPt,genVPhi,mu1d.Pt(),mu1d.Phi(),pU1,pU2,0,0,0,kFALSE,kFALSE);            
          }
          if(doStat){
            for(int i =0; i < nNV; i++){
              int ofs=i+ns;
              rcStatW[i]->CorrectInvCdf(metVars[ofs],metVarsPhi[ofs],genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,kFALSE,kTRUE);
            }
          }
        }
      }
      mtCorr  = sqrt( 2.0 * (lep->Pt()) * (metVars[cent]) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metVarsPhi[cent]))) );
    }
    relIso = pfCombIso/lep_raw->Pt(); 
    effSFweight=corr;
    outTree->Fill(); // add new info per event to the new tree
  }//end of loop over events
  // std::cout << "end loop over events"<< std::endl;
  delete rcMainWp;  delete rcMainWm;
  delete rcKeysWp;  delete rcKeysWm;
  for(int i = 0; i < nNV; i ++)delete rcStatW[i];
  delete rcEta05Wp;  delete rcEta051Wp;  delete rcEta1Wp;
  delete rcEta05Wm;  delete rcEta051Wm;  delete rcEta1Wm;
  // std::cout << "clean up memory" << std::endl;
    
  outFile->cd();
  hGenWeights->Write();
  outFile->Write();
  std::cout << "wrote outfile" << std::endl;

  delete intree;
  delete infile;
} // end of function
