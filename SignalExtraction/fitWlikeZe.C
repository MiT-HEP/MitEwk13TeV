//================================================================================================
//
// Perform fit to extract W->munu signal
//
//  * outputs plots and fit results summary
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
#include "../Utils/AppEffSF.cc"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_asym2.hh"
// #include "../Utils/RecoilCorrector_addJets.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// #include "ZBackgrounds.hh"

// helper class to handle efficiency tables
#include "../Utils/CEffUser1D.hh"
#include "../Utils/CEffUser2D.hh"

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
               const Double_t ksprob, const Double_t ksprobpe);

// make webpage
// void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void fitWlikeZe(const TString  outputDir,   // output directory
           const TString ntupleDir,
           const TString sqrts,
           const Double_t lumi,        // integrated luminosity (/fb)'
       const Double_t nsigma=0,     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
       const TString input_section = "1"
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // some flags to handle Recoil corrections
  bool doKeys = false;
  bool doInclusive = true;
  bool doEta = false;
  bool doDiago = false;
  bool doFootprint = false;
  bool doPF = false;
  bool doShittyRecoil=false;
  // some flags to handle the pileup Up/Down systematics
  bool pileupUp = false;
  bool pileupDown = false;
  
  
  // MET histogram binning and range
  const Int_t    NBINS   = 50;
  const Double_t METMAX  = 100;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  
  const Double_t ele_MASS  = 0.000511;
  const int NTOYS = 100;
  
  
   TString baseDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zee/";
  // TString baseDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results/Zee/";
  AppEffSF effs(baseDir);
  // effs.loadHLT("EleHLTEff_aMCxPythia","Combined","Combined");
  effs.loadHLT("EleHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("EleGSFSelEff_aMCxPythia","Combined","Combined");
  // effs.loadSel("EleGSFSelEff_aMCxPythia","Positive","Negative");
  // effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  string sysDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics/";
  string SysFileGSFSel = sysDir + "SysUnc_EleGSFSelEff.root";
  effs.loadUncSel(SysFileGSFSel);
  //
  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);

  // TFile *_rat2 = new TFile("TEST_Zmm_13TeV_incl_ZpTrw_v0/zPt_Normal13TeV.root");
  // TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  // hh_diff = (TH1D*)_rat2->Get("hZptRatio");
  

  // file format for output plots
  const TString format("png"); 


// const TString directory1("/home/sabrandt/SM/InclusiveMaster/CMSSW_7_6_3_patch2/src/MitEwk13TeV/Recoil/");
  // const TString directory2("../Recoil/");
  const TString directory2("/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Recoil");
  // const TString directory("/afs/cern.ch/user/d/dalfonso/public/WZ/JULY5");
  // // for Puppi, inclusive
  int rec_sig = 1;
  RecoilCorrector *recoilCorr = new  RecoilCorrector("","");
  // RecoilCorrector *recoilCorrm = new  RecoilCorrector("","");
  RecoilCorrector *recoilCorr05 = new  RecoilCorrector("","");
  RecoilCorrector *recoilCorr051 = new  RecoilCorrector("","");
  RecoilCorrector *recoilCorr1 = new  RecoilCorrector("","");
   RecoilCorrector *recoilCorrKeys = new  RecoilCorrector("","");
  
  
  if(doInclusive && !doDiago){
  
    // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s_Orig/WemMC_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    // recoilCorr->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    
    // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMC_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    // recoilCorr->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    recoilCorr->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G_bkg_fixRoch/",directory2.Data(),sqrts.Data()));
    recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_2G/",directory2.Data(),sqrts.Data()));
  
  } else if (doInclusive && doDiago){
    
    recoilCorr->loadRooWorkspacesDiagMCtoCorrect(Form("%s/WmmMC_PF_%s_2G/",directory2.Data(),sqrts.Data()), rec_sig,1);
    recoilCorr->loadRooWorkspacesDiagData(Form("%s/ZmmData_PF_%s_2G/",directory2.Data(),sqrts.Data()), rec_sig,1);
    recoilCorr->loadRooWorkspacesDiagMC(Form("%s/ZmmMC_PF_%s_2G/",directory2.Data(),sqrts.Data()), rec_sig,1);
  
      
  } else if (doEta){
    recoilCorr05->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMC_PF_%s_Eta1/",directory2.Data(),sqrts.Data()));
    recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta1/",directory2.Data(),sqrts.Data()));
    recoilCorr05->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta1/",directory2.Data(),sqrts.Data()));
    
    recoilCorr051->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMC_PF_%s_Eta2/",directory2.Data(),sqrts.Data()));
    recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta2/",directory2.Data(),sqrts.Data()));
    recoilCorr051->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta2/",directory2.Data(),sqrts.Data()));
    
    recoilCorr1->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMC_PF_%s_Eta3/",directory2.Data(),sqrts.Data()));
    recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Eta3/",directory2.Data(),sqrts.Data()));
    recoilCorr1->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Eta3/",directory2.Data(),sqrts.Data()));
  
  } else if (doKeys){
       
  recoilCorrKeys->loadRooWorkspacesMCtoCorrectKeys(Form("%s/ZmmMC_PF_%s_Keys/",directory2.Data(),sqrts.Data()));
  recoilCorrKeys->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_Keys/",directory2.Data(),sqrts.Data()));
  recoilCorrKeys->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_Keys/",directory2.Data(),sqrts.Data()));
  }

  //
  // input ntuple file names
  //
  enum { eData, eWmunu, eEWK, eBKG, eAntiData, eAntiWmunu, eAntiEWK };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  TString flav = "Zee";
   
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"));    typev.push_back(eData);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zee_select.root"));  typev.push_back(eWmunu);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx0_select.root"));  typev.push_back(eBKG);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx1_select.root"));  typev.push_back(eBKG);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx2_select.root"));  typev.push_back(eBKG);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eBKG);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.root"));  typev.push_back(eEWK);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.root"));  typev.push_back(eEWK);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.root"));  typev.push_back(eEWK);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top1_select.root")); typev.push_back(eEWK);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top2_select.root")); typev.push_back(eEWK);
  // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top_select.root")); typev.push_back(eEWK);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet   = new TH1D("hDataMet","",  NBINS,0,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm  = new TH1D("hDataMetm","", NBINS,0,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp  = new TH1D("hDataMetp","", NBINS,0,METMAX); hDataMetp->Sumw2();
  TH1D *hWmunuMet  = new TH1D("hWmunuMet","", NBINS,0,METMAX); hWmunuMet->Sumw2();
  TH1D *hWmunuMetp = new TH1D("hWmunuMetp","",NBINS,0,METMAX); hWmunuMetp->Sumw2();
  TH1D *hWmunuMetm = new TH1D("hWmunuMetm","",NBINS,0,METMAX); hWmunuMetm->Sumw2();
  
  
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,0,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,0,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,0,METMAX); hEWKMetm->Sumw2();
  

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet","",  NBINS,0,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,0,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,0,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWmunuMet  = new TH1D("hAntiWmunuMet","", NBINS,0,METMAX); hAntiWmunuMet->Sumw2();
  TH1D *hAntiWmunuMetp = new TH1D("hAntiWmunuMetp","",NBINS,0,METMAX); hAntiWmunuMetp->Sumw2();
  TH1D *hAntiWmunuMetm = new TH1D("hAntiWmunuMetm","",NBINS,0,METMAX); hAntiWmunuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet", "",  NBINS,0,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,0,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,0,METMAX); hAntiEWKMetm->Sumw2();
 
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy;
  // Float_t genlep->Pt(), genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q, q1, q2;
  TLorentzVector *lep=0, *lep_raw = 0, *genV=0, *metLep =0, *metLepRaw =0;
  TLorentzVector *lep2=0, *lep1 = 0, *lep1_raw=0, *lep2_raw=0;
  TLorentzVector *dilep=0;
  Float_t prefireWeight;
  Float_t lep1error, lep2error;
  
  TLorentzVector *genlep1=0, *genlep2=0;
//   TLorentzVector *lep2=0;
  // Float_t pfChIso, pfGamIso, pfNeuIso;
  UInt_t category;
    
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

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",  &genVy);   // GEN W boson phi (signal MC)
    // intree->SetBranchAddress("genlep->Pt()",   &genlep->Pt());    // GEN lepton pT (signal MC)
    // intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
    // intree->SetBranchAddress("puppiMet",      &met);       // MET
    // intree->SetBranchAddress("puppiMetPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("met",      &met);       // MET
    intree->SetBranchAddress("metPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("u1",       &u1);        // parallel component of recoil
    intree->SetBranchAddress("u2",       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q1",        &q1);         // lepton charge
    intree->SetBranchAddress("q2",        &q2);         // lepton charge
    intree->SetBranchAddress("lep1",      &lep1);
    intree->SetBranchAddress("lep2",      &lep2);
    intree->SetBranchAddress("lep1_raw",      &lep1_raw);
    intree->SetBranchAddress("lep2_raw",      &lep2_raw);
    intree->SetBranchAddress("lep1error",      &lep1error);
    intree->SetBranchAddress("lep2error",      &lep2error);
    intree->SetBranchAddress("genlep1",      &genlep1);
    intree->SetBranchAddress("genlep2",      &genlep2);
    intree->SetBranchAddress("dilep",      &dilep);
    intree->SetBranchAddress("genV",     &genV);       // lepton 4-vector
    intree->SetBranchAddress("category", &category);
    intree->SetBranchAddress("prefireWeight", &prefireWeight);

  
    Double_t mt=-999;
    TH1D* hGenWeights;
    double totalNorm = 1.0;
    cout << "Hello " << endl;
    if(typev[ifile] != eData ){
      cout << "get gen weights" << endl;
      hGenWeights = (TH1D*)infile->Get("hGenWeights");
      totalNorm = hGenWeights->Integral();
      cout << totalNorm << endl;
    }
    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<(int)(intree->GetEntries()*0.1); ientry++) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) std::cout << "On Entry.... " << ientry << std::endl;

      if(category != 1 && category != 2 && category != 3) continue;
      if(lep1->Pt() < PT_CUT || lep2->Pt() < PT_CUT) continue;
      if(fabs(lep1->Eta()) > ETA_CUT || fabs(lep2->Eta()) > ETA_CUT) continue;
      if(dilep->M() < MASS_LOW || dilep->M() > MASS_HIGH) continue;
      if(q1==q2) continue;

      // std::cout << "blah blah" << std::endl;
      double error;
      TVector2 vOtherLep;
      if(q1 > 0){
        q = q1;
        lep = lep1;
        lep_raw = lep1_raw;
        error = lep1error;
        metLep = lep2;
        metLepRaw = lep2_raw;
      } else {
        q = q2;
        lep = lep2;
        lep_raw = lep2_raw;
        error = lep2error;
        metLep = lep1;
        metLepRaw = lep1_raw;
      }

      
      double pU1         = 0;  //--
      double pU2         = 0;  //--
// cout << "lldl" << endl;
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

      q=1;
      if(fabs(lep->Eta()) > ETA_CUT) continue;
      // if(dilep->Pt() > 15 || dilep->Pt() < 5)continue;
  
      // mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );
      // cout << "lllaa" << endl;
      corr = effs.fullEfficiencies(lep1,q1,lep2,q2);

      TLorentzVector mu1;
      mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),ele_MASS);
      TLorentzVector mu2;
      mu2.SetPtEtaPhiM(metLep->Pt(),metLep->Eta(),metLep->Phi(),ele_MASS);
      if(mu1.Pt()        < PT_CUT)  continue;
      if(mu2.Pt()        < PT_CUT)  continue;
      // corrected (smear/scale) lepton for MET correction
            // cout << "lll" << endl;
      TVector2 vLepRaw((lep_raw->Pt())*cos(lep_raw->Phi()),(lep_raw->Pt())*sin(lep_raw->Phi()));
      TVector2 vOtherLepRaw((metLepRaw->Pt())*cos(metLepRaw->Phi()),(metLepRaw->Pt())*sin(metLepRaw->Phi()));
      TVector2 vMet((met)*cos(metPhi), (met)*sin(metPhi));
      TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
      vOtherLep.Set((mu2.Pt())*cos(mu2.Phi()),(mu2.Pt())*sin(mu2.Phi()));
      TVector2 vMetCorr = vMet + vOtherLepRaw;
      Double_t corrMetWithLepton = (vMetCorr).Mod();// + vLepRaw - vLepCor).Mod();
      
      if(typev[ifile]==eData || typev[ifile]==eAntiData){
        // corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor*(1-error)).Mod();
        if(typev[ifile]==eData) {
          hDataMet->Fill(corrMetWithLepton);
          if(q>0) { hDataMetp->Fill(corrMetWithLepton); }
          else    { hDataMetm->Fill(corrMetWithLepton); }
        } else if(typev[ifile]==eAntiData) {
          hAntiDataMet->Fill(corrMetWithLepton);
          if(q>0) { hAntiDataMetp->Fill(corrMetWithLepton); } 
          else    { hAntiDataMetm->Fill(corrMetWithLepton); }   
        }
      } else {
        Double_t weight = 1;Double_t weightUp = 1;Double_t weightDown = 1;
        Double_t weight2 =1;
        weight2*=scale1fb*lumi*corr*prefireWeight/totalNorm;
        weight *= scale1fb*lumi*corr*prefireWeight/totalNorm;

        // change to have rochester corrected muon and raw lepton with MET corrected same way as electron channel 
        if(typev[ifile]==eWmunu || typev[ifile]==eBKG) {

          if(lep->Pt()        > PT_CUT) {
            double bin = 0;
            // double w2 = 1;//hh_diff->GetBinContent(bin);
            // double recoilWeight = 1;
            // corrMetWithLepton=vMetCorr.Mod(), 
            double corrMetPhi=vMetCorr.Phi();
            hWmunuMet->Fill(corrMetWithLepton,weight);

            if(q>0) {
              double bin = 0;
              if(doKeys) {
                recoilCorrKeys->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              } else if(doEta) {
                if(fabs(dilep->Eta())<0.5)
                  recoilCorr05->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                else if (fabs(dilep->Eta())>=0.5 && fabs(dilep->Eta())<1.0)
                    recoilCorr051->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                else 
                  recoilCorr1->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago); 
              } else if(doInclusive){
                recoilCorr->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              }
              if(typev[ifile]==eWmunu) {
                hWmunuMetp->Fill(corrMetWithLepton,weight);
              } else {
                hEWKMetp->Fill(corrMetWithLepton,weight);
              }
            }
          }
        }
        // if(typev[ifile]==eAntiWmunu) {
          // if(lep->Pt()        < PT_CUT)  continue;
          // Double_t corrMetWithLepton=met, corrMetPhi=metPhi;
// //        this histogram isn't actually used in any results / fits
          // hAntiWmunuMet->Fill(corrMetWithLepton,weight2);
          // if(q>0) {              
            // pU1 = 0; pU2 = 0; 
            // double bin = 0;
            // // for(int i = 0; i <= hh_diff->GetNbinsX();++i){
              // // if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            // // }
              // double w2 = 1;// hh_diff->GetBinContent(bin);
              // double recoilWeight = 1;
              // if(doKeys) {
                // recoilCorrKeys->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              // } else if(doEta) {
                // if(fabs(dilep->Eta())<0.5)
                  // recoilCorr05->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                // else if (fabs(dilep->Eta())>=0.5 && fabs(dilep->Eta())<1.0)
                    // recoilCorr051->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                // else 
                  // recoilCorr1->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago); 
              // } else if(doInclusive){
                // recoilCorr->CorrectInvCdf(corrMetWithLepton,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              // }
            // TVector2 vMetCorr = vMet + vOtherLep;
            // Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            // hAntiWmunuMetp->Fill(corrMetWithLepton,weight2); 
            // corrMetWithLepton = met; corrMetPhi = metPhi;
          // } 
        // }
        if(typev[ifile]==eEWK) {
          if(lep->Pt()        < PT_CUT)  continue;
//           TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          TVector2 vMetCorr = vMet + vOtherLep;
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          hEWKMet->Fill(corrMetWithLepton,weight);
//           hEWKMet_PileupUp->Fill(corrMetWithLepton,weightUp);
//           hEWKMet_PileupDown->Fill(corrMetWithLepton,weightDown);
          if(q>0) {
            hEWKMetp->Fill(corrMetWithLepton,weight); 
            // RooFit doesn't see these histograms, only for combine, I will remove completely later
//             hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp); 
//             hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown); 
          }
          else { 
            hEWKMetm->Fill(corrMetWithLepton,weight); 
//             hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); 
//             hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); 
          }
        }
        // if(typev[ifile]==eAntiEWK) {
          // if(lep->Pt()        < PT_CUT)  continue;
// //           TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          // TVector2 vMetCorr = vMet + vOtherLep;
          // Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          // hAntiEWKMet->Fill(corrMetWithLepton,weight2);
          // if(q>0) { hAntiEWKMetp->Fill(corrMetWithLepton,weight2); }
          // else    { hAntiEWKMetm->Fill(corrMetWithLepton,weight2); }
        // }
      }
    }
  }
  delete infile;
  infile=0, intree=0;   
 
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",1.0*(hWmunuMetp->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  // cewk.setVal(hEWKMet->Integral()/hWmunuMet->Integral());
  cewk.setVal(0);
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWmunuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
//   dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",1.0*(hWmunuMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral(),0,hAntiDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWmunuMetp->Integral());
  // cewkp.setVal(0);
  cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  //RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral()*1.0,0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
//   dewkp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",1.0*(hWmunuMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWmunuMetm->Integral());
  cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWmunuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
//   dewkm.setConstant(kTRUE);
  RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));
  
  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,wmunuMetm,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  

    RooGaussian constm("constm","constm",nEWKm,RooConst(hEWKMetm->Integral()),RooConst(0.15*hEWKMetm->Integral()));
    RooGaussian constp("constp","constp",nEWKp,RooConst(hEWKMetp->Integral()),RooConst(0.15*hEWKMetp->Integral()));
 
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK),   RooArgList(nSig,nEWK));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp),RooArgList(nSigp,nEWKp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm),RooArgList(nSigm,nEWKm));
  

  
  // Anti-Signal PDFs
  RooDataHist awmunuMet ("awmunuMET", "awmunuMET", RooArgSet(pfmet),hAntiWmunuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,awmunuMet, 1);
  RooDataHist awmunuMetp("awmunuMETp","awmunuMETp",RooArgSet(pfmet),hAntiWmunuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,awmunuMetp,1);
  RooDataHist awmunuMetm("awmunuMETm","awmunuMETm",RooArgSet(pfmet),hAntiWmunuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,awmunuMetm,1); 
  
  // Anti-EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  

  
  // PDF for simultaneous fit
  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("Select");
  rooCat.defineType("Anti");
  
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(pdfMet, "Select");
  //pdfTotal.addPdf(apdfMet,"Anti");
  
  RooSimultaneous pdfTotalp("pdfTotalp","pdfTotalp",rooCat);
  pdfTotalp.addPdf(pdfMetp, "Select");
  // pdfTotalp.addPdf(apdfMetp,"Anti");
  
  RooSimultaneous pdfTotalm("pdfTotalm","pdfTotalm",rooCat);
  pdfTotalm.addPdf(pdfMetm, "Select");
  // pdfTotalm.addPdf(apdfMetm,"Anti");

  //
  // Perform fits
  //
  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);
 
  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
  RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
  RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);

  cout << "Starting values for Wmunu yields: " << endl;
  cout << "Selected: " << hDataMet->Integral() << endl;
  cout << "   sig: " << hWmunuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;

  cout << "Starting values for Wmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;

  cout << "Starting values for Wmunu_m yields: " << endl;
  cout << "   sig: " << hWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  
  
  cout << "Starting values for AntiWmunu yields: " << endl;
  cout << "Selected: " << hAntiDataMet->Integral() << endl;
  cout << "   sig: " << hAntiWmunuMet->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMet->Integral() << endl;

  cout << "Starting values for AntiWmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;

  cout << "Starting values for AntiWmunu_m yields: " << endl;
  cout << "   sig: " << hAntiWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);

  combine_workspace.import(pdfWm);
  combine_workspace.import(pdfWmp);
  combine_workspace.import(pdfWmm);
  // combine_workspace.import(pdfWm_RecoilUp);
  // combine_workspace.import(pdfWmp_RecoilUp);
  // combine_workspace.import(pdfWmm_RecoilUp);
  // combine_workspace.import(pdfWm_RecoilDown);
  // combine_workspace.import(pdfWmp_RecoilDown);
  // combine_workspace.import(pdfWmm_RecoilDown);
  // combine_workspace.import(pdfWm_ScaleUp);
  // combine_workspace.import(pdfWmp_ScaleUp);
  // combine_workspace.import(pdfWmm_ScaleUp);
  // combine_workspace.import(pdfWm_ScaleDown);
  // combine_workspace.import(pdfWmp_ScaleDown);
  // combine_workspace.import(pdfWmm_ScaleDown);
  // combine_workspace.import(pdfWm_PileupUp);
  // combine_workspace.import(pdfWmp_PileupUp);
  // combine_workspace.import(pdfWmm_PileupUp);
  // combine_workspace.import(pdfWm_PileupDown);
  // combine_workspace.import(pdfWmp_PileupDown);
  // combine_workspace.import(pdfWmm_PileupDown);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKm);
  // combine_workspace.import(pdfEWK_PileupUp);
  // combine_workspace.import(pdfEWKp_PileupUp);
  // combine_workspace.import(pdfEWKm_PileupUp);
  // combine_workspace.import(pdfEWK_PileupDown);
  // combine_workspace.import(pdfEWKp_PileupDown);
  // combine_workspace.import(pdfEWKm_PileupDown);
  
  combine_workspace.import(antiMet);
  combine_workspace.import(antiMetp);
  combine_workspace.import(antiMetm);
  
  combine_workspace.import(apdfWm);
  combine_workspace.import(apdfWmp);
  combine_workspace.import(apdfWmm);
  combine_workspace.import(apdfEWK);
  combine_workspace.import(apdfEWKp);
  combine_workspace.import(apdfEWKm);

  char nname[100];
  sprintf(nname,"%s/Wmunu_pdfTemplates.root",CPlot::sOutDir.Data());
  combine_workspace.writeToFile(nname);
  
  RooDataHist dataTotal("dataTotal","dataTotal", RooArgList(pfmet), Index(rooCat),
            Import("Select", dataMet),
            Import("Anti",   antiMet));
  RooFitResult *fitRes = 0;//pdfMet.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));//dataTotal.fitTo(dataMet,Extended(),Minos(kTRUE),Save(kTRUE));
  
 
  RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
             Import("Select", dataMetp),
             Import("Anti",   antiMetp));
             
  RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
  
  RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
             Import("Select", dataMetm),
             Import("Anti", antiMetm));
  RooFitResult *fitResm = 0;//pdfMetm.fitTo(dataMetm,Extended(),Minos(kTRUE),Save(kTRUE));
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal())/hPdfMetp->Integral());
  std::cout << "Scale = " << (nSigp.getVal()+nEWKp.getVal())/hPdfMetp->Integral() << std::endl;
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
   
   
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");
  TGaxis::SetMaxDigits(3);
  
  char ylabel[100];  // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 8 TeV",lumi*1000.);
  else         sprintf(lumitext,"%.1f fb^{-1}  at  #sqrt{s} = 13 TeV",lumi/1000.);
  
  // plot colors
  Int_t linecolorW   = kOrange-3;
  Int_t fillcolorW   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t ratioColor   = kGray+2;
  
  //
  // Dummy histograms for TLegend
  // (I can't figure out how to properly pass RooFit objects...)
  //
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  
  TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  hDummyW->SetLineColor(linecolorW);
  hDummyW->SetFillColor(fillcolorW);
  hDummyW->SetFillStyle(1001);
  
  TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  hDummyEWK->SetLineColor(linecolorEWK);
  hDummyEWK->SetFillColor(fillcolorEWK);
  hDummyEWK->SetFillStyle(1001);
  

   double yrange = 0.05;
   
   //
  // W MET plot
  //
  RooPlot *weframe = pfmet.frame(Bins(NBINS));
  weframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(weframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(weframe,LineColor(linecolorW));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK)),LineColor(linecolorEWK));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("wmunu_fitmet",weframe,"","mT [GeV]",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.AddTextBox("CMS",0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox(lumitext,0.63,0.92,0.95,0.99,0);
  plotMet.Draw(c,kTRUE,format,1);

  CPlot plotMetDiff("wmunu_fitmet","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-yrange,yrange);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, yrange*0.5,METMAX, yrange*0.5,kBlack,3);
  plotMetDiff.AddLine(0,-yrange*0.5,METMAX,-yrange*0.5,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  plotMetDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMet.SetName("wmunu_fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
    
  //
  // W+ MET plot
  //
  RooPlot *wepframe = pfmet.frame(Bins(NBINS));    
  wepframe->GetYaxis()->SetNdivisions(505);
  wepframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wepframe,LineColor(linecolorW));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp)),LineColor(linecolorEWK));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("wmunu_fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("wmunu_fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-yrange,yrange);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, yrange*0.5,METMAX, yrange*0.5,kBlack,3);
  plotMetpDiff.AddLine(0,-yrange*0.5,METMAX,-yrange*0.5,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("wmunu_fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  plotMetp.Draw(c,kTRUE,"pdf",1);

  //
  // W- MET plot
  //
  RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
  wemframe->GetYaxis()->SetNdivisions(505);
  wemframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wemframe,LineColor(linecolorW));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm)),LineColor(linecolorEWK));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("wmunu_fitmetm",wemframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("wmunu_fitmetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetmDiff.SetYRange(-yrange,yrange);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, yrange*0.5,METMAX, yrange*0.5,kBlack,3);
  plotMetmDiff.AddLine(0,-yrange*0.5,METMAX,-yrange*0.5,kBlack,3);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("wmunu_fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
  plotMetm.Draw(c,kTRUE,"pdf",1);
  
 
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  //

  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitWm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    if(hData->GetBinContent(ibin) == 0) diff = 0;
    // std::cout << "data " << hData->GetBinContent(ibin) << std::endl;
    // std::cout << "fits " << hFit->GetBinContent(ibin) << std::endl;
//     Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    if(hData->GetBinContent(ibin) == 0) err = 0;
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    // std::cout << "err = " << err << std::endl;
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);  
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  
  return hDiff;
}

//--------------------------------------------------------------------------------------------------
void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList parlist = res->floatParsFinal();
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist.getSize(); i++) {
    for(Int_t j=0; j<parlist.getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
               const Double_t ksprob, const Double_t ksprobpe)
{
  ios_base::fmtflags flags = os.flags();
  
  os << "  Chi2 Test" << endl;
  os << " -----------" << endl;
  os << "       prob = " << chi2prob << endl;
  os << "   chi2/ndf = " << chi2ndf << endl;
  os << endl;
  os << "  KS Test" << endl;
  os << " ---------" << endl;
  os << "   prob = " << ksprob << endl;
  os << "   prob = " << ksprobpe << " with 1000 pseudo-experiments" << endl;
  os << endl;
 
  os.flags(flags);
}
