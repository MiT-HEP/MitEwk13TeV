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

#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector_asym2.hh"
// #include "../Utils/RecoilCorrector_addJets.hh"
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// #include "ZBackgrounds.hh"

//helper class to handle rochester corrections
// #include <rochcor2015r.h>
// #include <muresolution_run2r.h>
#include <../RochesterCorr/RoccoR.cc>
#include "../Utils/AppEffSF.cc"

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
void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void fitWlikeZm(const TString  outputDir,   // output directory
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
  
  const Double_t mu_MASS = 0.1057;
  const int NTOYS = 100;
  
  
  
  // Double_t vEtaBinsFP[] = {2.5,2.1,1.5,0,-1.5,-2.1,-2.5};
  // Double_t vEtaBinsFP[] = {2.5,2.3,2.1,1.9,1.7,1.5,1.2,0.8,0.4,0,-0.4,-0.8,-1.2,-1.5,-1.7,-1.9,-2.1,-2.3,-2.5};
  Double_t vEtaBinsFP[] = {2.5,2.4,2.3,2.2,2.1,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.0,-1.1,-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5};
  Int_t nBinsFP = sizeof(vEtaBinsFP)/sizeof(vEtaBinsFP[0])-1;
  vector<Float_t> upper(nBinsFP);
  vector<Float_t> lower(nBinsFP);
  std::cout << "nEtaBins = " << nBinsFP << std::endl;
  
  // vector<Float_t> upper;
  for(int i = 0; i < nBinsFP; ++i){
	  std::cout << "bin val = " << vEtaBinsFP[i] << std::endl;
	  upper[i] = vEtaBinsFP[i];
	  lower[i] = vEtaBinsFP[i+1];
	  // lower.push_back(vEtaBinsFP[i+1]);
  }
   
  // vector<Float_t> upper; 
  // // upper.push_back(2.5); upper.push_back(2.3); upper.push_back(2.1); 
  // // upper.push_back(1.9); upper.push_back(1.7); upper.push_back(1.5);
  // // upper.push_back(1.2); upper.push_back(0.8); upper.push_back(0.4);
  // upper.push_back(2.5); upper.push_back(2.1); upper.push_back(1.5);upper.push_back(0);upper.push_back(-1.5); upper.push_back(-2.1);
  // vector<Float_t> lower; 
  // // lower.push_back(2.3);  lower.push_back(2.1);  lower.push_back(1.9);
  // // lower.push_back(1.7);  lower.push_back(1.5);  lower.push_back(1.2);
  // // lower.push_back(0.8); lower.push_back(0.4);   lower.push_back(0);
  // lower.push_back(2.1); lower.push_back(1.5);lower.push_back(0);lower.push_back(-1.5); lower.push_back(-2.1); lower.push_back(-2.5); 
  // // for electrons // Double_t vListOfFPs[] = {-0.00785963,-0.00993125,-0.0148858,-0.0130871,-0.0250669,-0.0173248,-0.0146016,-0.0136755,-0.0108911,-0.0116117,-0.0157723,-0.0194445,-0.0163345,-0.0230202,-0.0145553,-0.0160693,-0.00196283,-0.00298555};
  // // Double_t vListOfFPs[] = {-0.00109938, -0.00570463, -0.0101305,-0.00100094,-0.0110923, -0.00800751,-0.00763989,-0.00888777, -0.0113909,-0.0118738,-0.00842714,-0.00805976,-0.00790106,-0.00644149,-0.0106873,-0.00519317,-0.00817862, -0.0278066};
  // Double_t vListOfFPs[] = {0, -0.00109938, -0.00267664, -0.00935168, -0.0116393, -0.00877814, -0.00073066, -0.00168777, -0.0114005, -0.0102761, -0.010612, -0.00345421, -0.0092375, -0.00850389, -0.00988786, -0.00506354, -0.0064118, -0.00408221, -0.00398614, -0.0172608, -0.0104789, -0.007598, -0.0177603, -0.0129501, -0.00846278, -0.010801, -0.0179897, -0.014265, -0.004627, -0.00198961, -0.00940792, -0.0137091, -0.00842374, -0.00451072, -0.0134227, -0.013647, -0.000902762, -0.00728768, -0.00516646, -0.0112141, -0.000557311, -0.0124703, -0.0169047, -0.00473276, -0.00623816, -0.00367185, -0.0133091, -0.00340216, -0.0278066, 0};
  // vector<Float_t> fpcorr;
  // // std::cout << "size of list " << 
  // for(int i = 0; i < nBinsFP; ++i){
	  // std::cout << "loop bin " << i << std::endl;
	  // std::cout << vListOfFPs[i] << std::endl;
	  // fpcorr.push_back(vListOfFPs[i]);
  // }
  // std::cout << "done? " << std::endl;
  // fpcorr.push_back(-0.00923615); fpcorr.push_back(-0.0165217); fpcorr.push_back(-0.0136799); fpcorr.push_back(-0.0157187); fpcorr.push_back(-0.0172452); fpcorr.push_back(-0.00253613);
  // fpcorr.push_back(-0.000826797); fpcorr.push_back(0.00520995); fpcorr.push_back(0.0120186); fpcorr.push_back(0.0121537); fpcorr.push_back(0.0215729); fpcorr.push_back(0.0192083); fpcorr.push_back(0.014545); fpcorr.push_back(0.0155789);  fpcorr.push_back(0.0115706);
   
  TString baseDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zmm/";
  AppEffSF effs(baseDir);
  effs.loadHLT("MuHLTEff_aMCxPythia","Positive","Negative");
  effs.loadSel("MuSITEff_aMCxPythia","Combined","Combined");
  effs.loadSta("MuStaEff_aMCxPythia","Combined","Combined");
  
  string sysDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/Systematics/";
  string sysFileSIT = sysDir + "SysUnc_MuSITEff.root";
  string sysFileSta = sysDir + "SysUnc_MuStaEff.root";
  effs.loadUncSel(sysFileSIT);
  effs.loadUncSta(sysFileSta);
   
  // const TString baseDir = "/afs/cern.ch/work/s/sabrandt/public/FilesSM2017GH/Efficiency/LowPU2017ID_"+sqrts+"/results/Zmm/";
  // // "/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV/results/Zmm/";
  // // path needs to be updated ^
  // const TString dataHLTEffName_pos = baseDir + "Data/MuHLTEff_aMCxPythia/Positive/eff.root";
  // const TString dataHLTEffName_neg = baseDir + "Data/MuHLTEff_aMCxPythia/Negative/eff.root";
  // const TString zmmHLTEffName_pos  = baseDir + "MC/MuHLTEff_aMCxPythia/Positive/eff.root";
  // const TString zmmHLTEffName_neg  = baseDir + "MC/MuHLTEff_aMCxPythia/Negative/eff.root";

  // const TString dataSelEffName_pos = baseDir + "Data/MuSITEff_aMCxPythia/Positive/eff.root";
  // const TString dataSelEffName_neg = baseDir + "Data/MuSITEff_aMCxPythia/Negative/eff.root";
  // const TString zmmSelEffName_pos  = baseDir + "MC/MuSITEff_aMCxPythia/Positive/eff.root";
  // const TString zmmSelEffName_neg  = baseDir + "MC/MuSITEff_aMCxPythia/Negative/eff.root";

  // const TString dataStaEffName_pos = baseDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  // const TString dataStaEffName_neg = baseDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  // const TString zmmStaEffName_pos  = baseDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";
  // const TString zmmStaEffName_neg  = baseDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";
  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);

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
  
  
    recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMC_PF_%s_2G/",directory2.Data(),sqrts.Data()));
    recoilCorr->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G_bkg/",directory2.Data(),sqrts.Data()));
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
  
  TString flav = "Zmumu";

  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"));    typev.push_back(eData);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zmm_select.raw.root"));  typev.push_back(eWmunu);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx_select.raw.root"));  typev.push_back(eBKG);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.raw.root")); typev.push_back(eBKG);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.raw.root"));  typev.push_back(eEWK);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.raw.root"));  typev.push_back(eEWK);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.raw.root"));  typev.push_back(eEWK);
  fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top_select.raw.root")); typev.push_back(eEWK);


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
  
  // TH1D *hEWKMet_PileupUp   = new TH1D("hEWKMet_PileupUp",  "",NBINS,0,METMAX); hEWKMet_PileupUp->Sumw2();
  // TH1D *hEWKMetp_PileupUp  = new TH1D("hEWKMetp_PileupUp", "",NBINS,0,METMAX); hEWKMetp_PileupUp->Sumw2();
  // TH1D *hEWKMetm_PileupUp  = new TH1D("hEWKMetm_PileupUp", "",NBINS,0,METMAX); hEWKMetm_PileupUp->Sumw2();
  
  // TH1D *hEWKMet_PileupDown   = new TH1D("hEWKMet_PileupDown",  "",NBINS,0,METMAX); hEWKMet_PileupDown->Sumw2();
  // TH1D *hEWKMetp_PileupDown  = new TH1D("hEWKMetp_PileupDown", "",NBINS,0,METMAX); hEWKMetp_PileupDown->Sumw2();
  // TH1D *hEWKMetm_PileupDown  = new TH1D("hEWKMetm_PileupDown", "",NBINS,0,METMAX); hEWKMetm_PileupDown->Sumw2();
  // TH1D *hWmunuMet_RecoilUp  = new TH1D("hWmunuMet_RecoilUp", "",NBINS,0,METMAX); hWmunuMet_RecoilUp->Sumw2();
  // TH1D *hWmunuMetp_RecoilUp = new TH1D("hWmunuMetp_RecoilUp","",NBINS,0,METMAX); hWmunuMetp_RecoilUp->Sumw2();
  // TH1D *hWmunuMetm_RecoilUp = new TH1D("hWmunuMetm_RecoilUp","",NBINS,0,METMAX); hWmunuMetm_RecoilUp->Sumw2();
  // TH1D *hWmunuMet_RecoilDown  = new TH1D("hWmunuMet_RecoilDown", "",NBINS,0,METMAX); hWmunuMet_RecoilDown->Sumw2();
  // TH1D *hWmunuMetp_RecoilDown = new TH1D("hWmunuMetp_RecoilDown","",NBINS,0,METMAX); hWmunuMetp_RecoilDown->Sumw2();
  // TH1D *hWmunuMetm_RecoilDown = new TH1D("hWmunuMetm_RecoilDown","",NBINS,0,METMAX); hWmunuMetm_RecoilDown->Sumw2();

  // TH1D *hWmunuMet_RecoilCUp  = new TH1D("hWmunuMet_RecoilCUp", "",NBINS,0,METMAX); hWmunuMet_RecoilCUp->Sumw2();
  // TH1D *hWmunuMetp_RecoilCUp = new TH1D("hWmunuMetp_RecoilCUp","",NBINS,0,METMAX); hWmunuMetp_RecoilCUp->Sumw2();
  // TH1D *hWmunuMetm_RecoilCUp = new TH1D("hWmunuMetm_RecoilCUp","",NBINS,0,METMAX); hWmunuMetm_RecoilCUp->Sumw2();
  // TH1D *hWmunuMet_RecoilCDown  = new TH1D("hWmunuMet_RecoilCDown", "",NBINS,0,METMAX); hWmunuMet_RecoilCDown->Sumw2();
  // TH1D *hWmunuMetp_RecoilCDown = new TH1D("hWmunuMetp_RecoilCDown","",NBINS,0,METMAX); hWmunuMetp_RecoilCDown->Sumw2();
  // TH1D *hWmunuMetm_RecoilCDown = new TH1D("hWmunuMetm_RecoilCDown","",NBINS,0,METMAX); hWmunuMetm_RecoilCDown->Sumw2();

  // TH1D *hWmunuMet_ScaleUp  = new TH1D("hWmunuMet_ScaleUp", "",NBINS,0,METMAX); hWmunuMet_ScaleUp->Sumw2();
  // TH1D *hWmunuMetp_ScaleUp = new TH1D("hWmunuMetp_ScaleUp","",NBINS,0,METMAX); hWmunuMetp_ScaleUp->Sumw2();
  // TH1D *hWmunuMetm_ScaleUp = new TH1D("hWmunuMetm_ScaleUp","",NBINS,0,METMAX); hWmunuMetm_ScaleUp->Sumw2();
  // TH1D *hWmunuMet_ScaleDown  = new TH1D("hWmunuMet_ScaleDown", "",NBINS,0,METMAX); hWmunuMet_ScaleDown->Sumw2();
  // TH1D *hWmunuMetp_ScaleDown = new TH1D("hWmunuMetp_ScaleDown","",NBINS,0,METMAX); hWmunuMetp_ScaleDown->Sumw2();
  // TH1D *hWmunuMetm_ScaleDown = new TH1D("hWmunuMetm_ScaleDown","",NBINS,0,METMAX); hWmunuMetm_ScaleDown->Sumw2();
  
  // TH1D *hWmunuMet_PileupUp  = new TH1D("hWmunuMet_PileupUp", "",NBINS,0,METMAX); hWmunuMet_PileupUp->Sumw2();
  // TH1D *hWmunuMetp_PileupUp = new TH1D("hWmunuMetp_PileupUp","",NBINS,0,METMAX); hWmunuMetp_PileupUp->Sumw2();
  // TH1D *hWmunuMetm_PileupUp = new TH1D("hWmunuMetm_PileupUp","",NBINS,0,METMAX); hWmunuMetm_PileupUp->Sumw2();
  // TH1D *hWmunuMet_PileupDown  = new TH1D("hWmunuMet_PileupDown", "",NBINS,0,METMAX); hWmunuMet_PileupDown->Sumw2();
  // TH1D *hWmunuMetp_PileupDown = new TH1D("hWmunuMetp_PileupDown","",NBINS,0,METMAX); hWmunuMetp_PileupDown->Sumw2();
  // TH1D *hWmunuMetm_PileupDown = new TH1D("hWmunuMetm_PileupDown","",NBINS,0,METMAX); hWmunuMetm_PileupDown->Sumw2();

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet","",  NBINS,0,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,0,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,0,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWmunuMet  = new TH1D("hAntiWmunuMet","", NBINS,0,METMAX); hAntiWmunuMet->Sumw2();
  TH1D *hAntiWmunuMetp = new TH1D("hAntiWmunuMetp","",NBINS,0,METMAX); hAntiWmunuMetp->Sumw2();
  TH1D *hAntiWmunuMetm = new TH1D("hAntiWmunuMetm","",NBINS,0,METMAX); hAntiWmunuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet", "",  NBINS,0,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,0,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,0,METMAX); hAntiEWKMetm->Sumw2();
  
//    //
//   // HLT efficiency
//   //
  cout << "Loading trigger efficiencies..." << endl;
 
  TH2D *hErr  = new TH2D("hErr", "",10,0,10,20,0,20);
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy;
  // Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q, q1, q2;
  TLorentzVector *lep=0, *genV=0, *metLep =0;
  TLorentzVector *lep2=0, *lep1 = 0;
  TLorentzVector *dilep=0;
  Float_t prefireWeight;
  
  TLorentzVector *genlep1=0, *genlep2=0;
//   TLorentzVector *lep2=0;
  // Float_t pfChIso, pfGamIso, pfNeuIso;
  UInt_t category;
    
  // rochcor2015 *rmcor = new rochcor2015();
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

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",  &genVy);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",      &met);       // MET
    intree->SetBranchAddress("metPhi",   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("u1",       &u1);        // parallel component of recoil
    intree->SetBranchAddress("u2",       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q1",        &q1);         // lepton charge
    intree->SetBranchAddress("q2",        &q2);         // lepton charge
    intree->SetBranchAddress("lep1",      &lep1);
    intree->SetBranchAddress("lep2",      &lep2);
    intree->SetBranchAddress("genlep1",      &genlep1);
    intree->SetBranchAddress("genlep2",      &genlep2);
    intree->SetBranchAddress("dilep",      &dilep);
//     intree->SetBranchAddress(otherlep.c_str(),      &lep2);       // lepton 4-vector
//     intree->SetBranchAddress("lep1",      &lep1);       // lepton 4-vector
//     intree->SetBranchAddress("lep2",      &lep2);       // lepton 4-vector
    intree->SetBranchAddress("genV",     &genV);       // lepton 4-vector
    // intree->SetBranchAddress("pfChIso",  &pfChIso);
    // intree->SetBranchAddress("pfGamIso", &pfGamIso);
    // intree->SetBranchAddress("pfNeuIso", &pfNeuIso);
    intree->SetBranchAddress("category", &category);
    intree->SetBranchAddress("prefireWeight", &prefireWeight);

  
    Double_t mt=-999;

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

//       std::cout << "blah blah" << std::endl;
      // vector containing raw lepton info for correcting MET
      
      TVector2 vOtherLep;
      // TVector2 vOtherLepRaw;
      // use the + lepton as our main lepton
      if(q1 < 0){
        q = q1;
        lep = lep1;
        metLep = lep2;
        
      } else {
        q = q2;
        lep = lep2;
        metLep = lep1;
        // vOtherLepRaw.Set((lep1->Pt())*cos(lep1->Phi()),(lep1->Pt())*sin(lep1->Phi()));
      }
      TVector2 vLepRaw((lep->Pt())*cos(lep->Phi()),(lep->Pt())*sin(lep->Phi()));
      TVector2 vOtherLepRaw((metLep->Pt())*cos(metLep->Phi()),(metLep->Pt())*sin(metLep->Phi()));
      TVector2 vMet((met)*cos(metPhi), (met)*sin(metPhi));
      
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
      
      
      q=1;
      if(fabs(lep->Eta()) > ETA_CUT) continue;
  
      mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );

      // corr = effs.computeHLTSF(lep,q)*effs.computeSelSF(lep,q)*effs.computeStaSF(lep,q);
      corr = effs.fullEfficiencies(lep,q);
      // cout << corr1 << " " << corr << endl;
      
      vector<double> uncs_sta = effs.getUncSta(lep,q);
      vector<double> uncs_sit = effs.getUncSel(lep,q);
      
      // corrFSR *= uncs_sta[0]*uncs_sit[0]*effs.computeHLTSF(lep,q); // alternate fsr model
      // corrMC  *= uncs_sta[1]*uncs_sit[1]*effs.computeHLTSF(lep,q); // alternate mc gen model
      // corrBkg *= uncs_sta[2]*uncs_sit[2]*effs.computeHLTSF(lep,q); // alternate bkg model
      // corrTag *= uncs_sta[3]*uncs_sit[3]*effs.computeHLTSF(lep,q); // alternate bkg model
      // corr *= effdata/effmc; // orig
      
      
      // double var=0.;        
      // var += effs.statUncSta(&mu1, q1, hErr, hErr, fabs(weight)*corr);
      // var += effs.statUncSta(&mu2, q2, hErr, hErr, fabs(weight)*corr);
      // var += effs.statUncSel(&mu1, q1, hErr, hErr, fabs(weight)*corr);
      // var += effs.statUncSel(&mu2, q2, hErr, hErr, fabs(weight)*corr);
      // var += effs.statUncHLTDilep(&mu1, q1, &mu2, q2);
	  

      
      
      
      if(typev[ifile]==eData || typev[ifile]==eAntiData){
        TLorentzVector mu1;        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        TLorentzVector mu2;        mu2.SetPtEtaPhiM(metLep->Pt(),metLep->Eta(),metLep->Phi(),mu_MASS);
        float qter1=1.0;
        double dtSF1 = rc.kScaleDT(q, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
        mu1*=dtSF1;
        double dtSF2 = rc.kScaleDT(-q, mu2.Pt(), mu2.Eta(), mu2.Phi());//, s=0, m=0);
        mu1*=dtSF2;
        // rmcor->momcor_data(mu1,q,0,qter1);
        // rmcor->momcor_data(mu2,-q,0,qter1);
        if(mu1.Pt()        < PT_CUT)  continue;
        if(mu2.Pt()        < PT_CUT)  continue;
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
        vOtherLep.Set((mu2.Pt())*cos(mu2.Phi()),(mu2.Pt())*sin(mu2.Phi()));
        TVector2 vMetCorr = vMet + vOtherLepRaw;
        Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
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

        // if(pileupUp) {
          // weight *= scale1fbUp*lumi*corr;
          // weight2 *= scale1fbUp*lumi*corr;
        // } else if(pileupDown){
          // weight *= scale1fbDown*lumi*corr;
          // weight2 *= scale1fbDown*lumi*corr;
        // } else {
          weight2*=scale1fb*lumi*corr*prefireWeight;
          weight *= scale1fb*lumi*corr*prefireWeight;
        // }
        // Double_t eleCorr=0.0;
        // Double_t eleCorrOther=0.0;
        // Do some Rochester corrections for MC
        TLorentzVector mu1;
        TLorentzVector mu2;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        mu2.SetPtEtaPhiM(metLep->Pt(),metLep->Eta(),metLep->Phi(),mu_MASS);
        float qter1=1.0;
        double mcSF1 = 1;
        genlep1->Pt()>0? mcSF1 = rc.kSpreadMC(q, mu1.Pt(), mu1.Eta(), mu1.Phi(), genlep1->Pt()):mcSF1=1;
        mu1*=mcSF1;
        double mcSF2 = 1;
        genlep2->Pt()>0?mcSF2 = rc.kSpreadMC(q, mu2.Pt(), mu2.Eta(), mu2.Phi(), genlep2->Pt()) : mcSF2=1;
        mu1*=mcSF2;
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
	    	// vLepCor=vLepCor*(1-eleCorr);
        Double_t lepPt = vLepCor.Mod();
        
        // change to have rochester corrected muon and raw lepton with MET corrected same way as electron channel 
        if(typev[ifile]==eWmunu || typev[ifile]==eBKG) {
		       vOtherLep.Set((mu2.Pt())*cos(mu2.Phi()),(mu2.Pt())*sin(mu2.Phi()));
          TVector2 vMetCorr = vMet + vOtherLepRaw;
          Double_t corrMet=met, corrMetPhi=metPhi;
          if(lep->Pt()        > PT_CUT) {
            double bin = 0;
            double w2 = 1;//hh_diff->GetBinContent(bin);
            double recoilWeight = 1;
            corrMet=vMetCorr.Mod(), corrMetPhi=vMetCorr.Phi();
            hWmunuMet->Fill(corrMet,weight);
            if(q>0) {
              double bin = 0;
              double w2 =  1;//hh_diff->GetBinContent(bin);
              double recoilWeight = 1;
              if(doKeys) {
                recoilCorrKeys->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              } else if(doEta) {
                if(fabs(dilep->Eta())<0.5)
                  recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                else if (fabs(dilep->Eta())>=0.5 && fabs(dilep->Eta())<1.0)
                    recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                else 
                  recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago); 
              } else if(doInclusive){
                recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              }
              TVector2 vMetCorrAfter((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
              
              Double_t corrMetWithLepton = (vMetCorrAfter + vLepRaw - vLepCor).Mod();
              
        // std::cout << "hello5 "<< std::endl;
              if(typev[ifile]==eWmunu)
              {
                hWmunuMetp->Fill(corrMetWithLepton,weight*w2*recoilWeight);
//                 hWmunuMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
//                 hWmunuMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
              }
              else
              {
                hEWKMetp->Fill(corrMetWithLepton,weight*recoilWeight);
              }
              corrMet=vMetCorr.Mod(), corrMetPhi=vMetCorr.Phi();
//               corrMet=met, corrMetPhi=metPhi;
            } 
            corrMet=met, corrMetPhi=metPhi;
          }
          //             recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),1,q);
          // unused?
          // hWmunuMet_RecoilUp->Fill(corrMet,weight);
          // if(q>0) {
            // pU1 = 0; pU2 = 0; 
            // hWmunuMetp_RecoilUp->Fill(corrMet,weight); 
            // corrMet=met, corrMetPhi=metPhi;
          // } else { 
            // pU1 = 0; pU2 = 0; 
            // hWmunuMetm_RecoilUp->Fill(corrMet,weight);
            // corrMet=met, corrMetPhi=metPhi;
          // }
          // hWmunuMet_RecoilDown->Fill(corrMet,weight);
          // if(q>0) {
            // pU1 = 0; pU2 = 0; 
            // hWmunuMetp_RecoilDown->Fill(corrMet,weight);
            // corrMet=met, corrMetPhi=metPhi;
          // } else {
            // pU1 = 0; pU2 = 0; 
            // hWmunuMetm_RecoilDown->Fill(corrMet,weight);
            // corrMet=met, corrMetPhi=metPhi;
          // }
        }
        // if(typev[ifile]==eAntiWmunu) {
          // if(lep->Pt()        < PT_CUT)  continue;
          // Double_t corrMet=met, corrMetPhi=metPhi;
// //        this histogram isn't actually used in any results / fits
          // hAntiWmunuMet->Fill(corrMet,weight2);
          // if(q>0) {              
            // pU1 = 0; pU2 = 0; 
            // double bin = 0;
            // // for(int i = 0; i <= hh_diff->GetNbinsX();++i){
              // // if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            // // }
              // double w2 = 1;// hh_diff->GetBinContent(bin);
              // double recoilWeight = 1;
              // if(doKeys) {
                // recoilCorrKeys->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              // } else if(doEta) {
                // if(fabs(dilep->Eta())<0.5)
                  // recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                // else if (fabs(dilep->Eta())>=0.5 && fabs(dilep->Eta())<1.0)
                    // recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
                // else 
                  // recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago); 
              // } else if(doInclusive){
                // recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              // }
            // TVector2 vMetCorr = vMet + vOtherLep;
            // Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            // hAntiWmunuMetp->Fill(corrMetWithLepton,weight2); 
            // corrMet = met; corrMetPhi = metPhi;
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
 
  // // Here rescale the up/down whatever
// //   Calculate the shapes for W+
  // TH1D *corrP = (TH1D*) hWmunuMetp->Clone("up");
  // *corrP = (*hh_diffp)*(*hWmunuMetp);
  // hWmunuMetp_RecoilUp->Add(corrP,hWmunuMetp,-1);
  // hWmunuMetp_RecoilDown->Add(corrP,hWmunuMetp,1);
// //   hWmunuMetp_ScaleUp->Add(corrP,hWmunuMetp,-1);
// //   hWmunuMetp_ScaleDown->Add(corrP,hWmunuMetp,1);
  // // Calculate the shapes for W-
  // TH1D *corrM = (TH1D*) hWmunuMetm->Clone("up");
  // *corrM = (*hh_diffm)*(*hWmunuMetm);
  // hWmunuMetm_RecoilUp->Add(corrM,hWmunuMetm,-1);
  // hWmunuMetm_RecoilDown->Add(corrM,hWmunuMetm,1);
// //   hWmunuMetm_ScaleUp->Add(corrM,hWmunuMetm,-1);
// //   hWmunuMetm_ScaleDown->Add(corrM,hWmunuMetm,1);
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",1.0*(hWmunuMetp->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0*(hDataMet->Integral()),0,hDataMet->Integral());
  nQCD.setConstant(kTRUE);
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  // cewk.setVal(hEWKMet->Integral()/hWmunuMet->Integral());
  cewk.setVal(0);
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  // RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWmunuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  // RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  // RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  // dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
//   dewk.setConstant(kTRUE);
  // RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",1.0*(hWmunuMetp->Integral()),0,hDataMetp->Integral());
  // RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral(),0,hAntiDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",0,0,hDataMetp->Integral());
  nQCDp.setConstant(kTRUE);
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWmunuMetp->Integral());
  // cewkp.setVal(0);
  cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  //RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral()*1.0,0,hAntiDataMetp->Integral());
  // RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  // RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  // dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
//   dewkp.setConstant(kTRUE);
  // RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",1.0*(hWmunuMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",0,0,hDataMetm->Integral());
  nQCDm.setConstant(kTRUE);
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWmunuMetm->Integral());
  cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  // RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWmunuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  // RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  // RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  // dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
//   dewkm.setConstant(kTRUE);
  // RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));
  
  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,wmunuMetm,1); 
  // RooDataHist wmunuMet_RecoilUp("wmunuMET_RecoilUp", "wmunuMET_RecoilUp", RooArgSet(pfmet),hWmunuMet_RecoilUp);  RooHistPdf pdfWm_RecoilUp("wm_RecoilUp", "wm_RecoilUp", pfmet,wmunuMet_RecoilUp, 1);
  // RooDataHist wmunuMetp_RecoilUp("wmunuMETp_RecoilUp","wmunuMETp_RecoilUp",RooArgSet(pfmet),hWmunuMetp_RecoilUp); RooHistPdf pdfWmp_RecoilUp("wmp_RecoilUp","wmp_RecoilUp",pfmet,wmunuMetp_RecoilUp,1);
  // RooDataHist wmunuMetm_RecoilUp("wmunuMETm_RecoilUp","wmunuMETm_RecoilUp",RooArgSet(pfmet),hWmunuMetm_RecoilUp); RooHistPdf pdfWmm_RecoilUp("wmm_RecoilUp","wmm_RecoilUp",pfmet,wmunuMetm_RecoilUp,1); 
  // RooDataHist wmunuMet_RecoilDown("wmunuMET_RecoilDown", "wmunuMET_RecoilDown", RooArgSet(pfmet),hWmunuMet_RecoilDown);  RooHistPdf pdfWm_RecoilDown("wm_RecoilDown", "wm_RecoilDown", pfmet,wmunuMet_RecoilDown, 1);
  // RooDataHist wmunuMetp_RecoilDown("wmunuMETp_RecoilDown","wmunuMETp_RecoilDown",RooArgSet(pfmet),hWmunuMetp_RecoilDown); RooHistPdf pdfWmp_RecoilDown("wmp_RecoilDown","wmp_RecoilDown",pfmet,wmunuMetp_RecoilDown,1);
  // RooDataHist wmunuMetm_RecoilDown("wmunuMETm_RecoilDown","wmunuMETm_RecoilDown",RooArgSet(pfmet),hWmunuMetm_RecoilDown); RooHistPdf pdfWmm_RecoilDown("wmm_RecoilDown","wmm_RecoilDown",pfmet,wmunuMetm_RecoilDown,1); 
   // RooDataHist wmunuMet_ScaleUp("wmunuMET_ScaleUp", "wmunuMET_ScaleUp", RooArgSet(pfmet),hWmunuMet_ScaleUp);  RooHistPdf pdfWm_ScaleUp("wm_ScaleUp", "wm_ScaleUp", pfmet,wmunuMet_ScaleUp, 1);
  // RooDataHist wmunuMetp_ScaleUp("wmunuMETp_ScaleUp","wmunuMETp_ScaleUp",RooArgSet(pfmet),hWmunuMetp_ScaleUp); RooHistPdf pdfWmp_ScaleUp("wmp_ScaleUp","wmp_ScaleUp",pfmet,wmunuMetp_ScaleUp,1);
  // RooDataHist wmunuMetm_ScaleUp("wmunuMETm_ScaleUp","wmunuMETm_ScaleUp",RooArgSet(pfmet),hWmunuMetm_ScaleUp); RooHistPdf pdfWmm_ScaleUp("wmm_ScaleUp","wmm_ScaleUp",pfmet,wmunuMetm_ScaleUp,1); 
  // RooDataHist wmunuMet_ScaleDown("wmunuMET_ScaleDown", "wmunuMET_ScaleDown", RooArgSet(pfmet),hWmunuMet_ScaleDown);  RooHistPdf pdfWm_ScaleDown("wm_ScaleDown", "wm_ScaleDown", pfmet,wmunuMet_ScaleDown, 1);
  // RooDataHist wmunuMetp_ScaleDown("wmunuMETp_ScaleDown","wmunuMETp_ScaleDown",RooArgSet(pfmet),hWmunuMetp_ScaleDown); RooHistPdf pdfWmp_ScaleDown("wmp_ScaleDown","wmp_ScaleDown",pfmet,wmunuMetp_ScaleDown,1);
  // RooDataHist wmunuMetm_ScaleDown("wmunuMETm_ScaleDown","wmunuMETm_ScaleDown",RooArgSet(pfmet),hWmunuMetm_ScaleDown); RooHistPdf pdfWmm_ScaleDown("wmm_ScaleDown","wmm_ScaleDown",pfmet,wmunuMetm_ScaleDown,1); 
  
    // RooDataHist wmunuMet_PileupUp("wmunuMET_PileupUp", "wmunuMET_PileupUp", RooArgSet(pfmet),hWmunuMet_PileupUp);  RooHistPdf pdfWm_PileupUp("wm_PileupUp", "wm_PileupUp", pfmet,wmunuMet_PileupUp, 1);
  // RooDataHist wmunuMetp_PileupUp("wmunuMETp_PileupUp","wmunuMETp_PileupUp",RooArgSet(pfmet),hWmunuMetp_PileupUp); RooHistPdf pdfWmp_PileupUp("wmp_PileupUp","wmp_PileupUp",pfmet,wmunuMetp_PileupUp,1);
  // RooDataHist wmunuMetm_PileupUp("wmunuMETm_PileupUp","wmunuMETm_PileupUp",RooArgSet(pfmet),hWmunuMetm_PileupUp); RooHistPdf pdfWmm_PileupUp("wmm_PileupUp","wmm_PileupUp",pfmet,wmunuMetm_PileupUp,1); 
  // RooDataHist wmunuMet_PileupDown("wmunuMET_PileupDown", "wmunuMET_PileupDown", RooArgSet(pfmet),hWmunuMet_PileupDown);  RooHistPdf pdfWm_PileupDown("wm_PileupDown", "wm_PileupDown", pfmet,wmunuMet_PileupDown, 1);
  // RooDataHist wmunuMetp_PileupDown("wmunuMETp_PileupDown","wmunuMETp_PileupDown",RooArgSet(pfmet),hWmunuMetp_PileupDown); RooHistPdf pdfWmp_PileupDown("wmp_PileupDown","wmp_PileupDown",pfmet,wmunuMetp_PileupDown,1);
  // RooDataHist wmunuMetm_PileupDown("wmunuMETm_PileupDown","wmunuMETm_PileupDown",RooArgSet(pfmet),hWmunuMetm_PileupDown); RooHistPdf pdfWmm_PileupDown("wmm_PileupDown","wmm_PileupDown",pfmet,wmunuMetm_PileupDown,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
    // RooDataHist ewkMet_PileupUp("ewkMET_PileupUp", "ewkMET_PileupUp", RooArgSet(pfmet),hEWKMet_PileupUp);  RooHistPdf pdfEWK_PileupUp("ewk_PileupUp", "ewk_PileupUp", pfmet,ewkMet_PileupUp, 1);
  // RooDataHist ewkMetp_PileupUp("ewkMETp_PileupUp","ewkMETp_PileupUp",RooArgSet(pfmet),hEWKMetp_PileupUp); RooHistPdf pdfEWKp_PileupUp("ewkp_PileupUp","ewkp_PileupUp",pfmet,ewkMetp_PileupUp,1);
  // RooDataHist ewkMetm_PileupUp("ewkMETm_PileupUp","ewkMETm_PileupUp",RooArgSet(pfmet),hEWKMetm_PileupUp); RooHistPdf pdfEWKm_PileupUp("ewkm_PileupUp","ewkm_PileupUp",pfmet,ewkMetm_PileupUp,1); 
  // RooDataHist ewkMet_PileupDown("ewkMET_PileupDown", "ewkMET_PileupDown", RooArgSet(pfmet),hEWKMet_PileupDown);  RooHistPdf pdfEWK_PileupDown("ewk_PileupDown", "ewk_PileupDown", pfmet,ewkMet_PileupDown, 1);
  // RooDataHist ewkMetp_PileupDown("ewkMETp_PileupDown","ewkMETp_PileupDown",RooArgSet(pfmet),hEWKMetp_PileupDown); RooHistPdf pdfEWKp_PileupDown("ewkp_PileupDown","ewkp_PileupDown",pfmet,ewkMetp_PileupDown,1);
  // RooDataHist ewkMetm_PileupDown("ewkMETm_PileupDown","ewkMETm_PileupDown",RooArgSet(pfmet),hEWKMetm_PileupDown); RooHistPdf pdfEWKm_PileupDown("ewkm_PileupDown","ewkm_PileupDown",pfmet,ewkMetm_PileupDown,1); 
  
//   // test using the reversed dEta, dPhi cuts as background
//   RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
//   RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
//   RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  // QCD Pdfs
  
  //CExponential qcd(pfmet,kTRUE);
  //CExponential qcdp(pfmet,kTRUE);
  //CExponential qcdm(pfmet,kTRUE);
  // comment back in for qcd functional form
  // CPepeModel2 qcd("qcd",pfmet);
  // CPepeModel2 qcdp("qcdp",pfmet);
  // CPepeModel2 qcdm("qcdm",pfmet);
  //CPepeModel2 qcdm("qcdm",pfmet);
  //CPepeModel1 qcdm("qcdm",pfmet);
//  qcdp.a1->setConstant(kTRUE);
//   qcdp.a1->setVal(0.24);
//    qcdm.a1->setVal(0.24);

//   RooRealVar a1ConstMeanP("a1ConstMeanP","a1ConstMeanP",0.18);
//   RooRealVar a1ConstSigmaP("a1ConstSigmaP","a1ConstSigmaP",0.003);
// //   RooRealVar a1ConstMeanM("a1ConstMeanM","a1ConstMeanM",0.215);
// //   RooRealVar a1ConstSigmaM("a1ConstSigmaM","a1ConstSigmaM",0.007);
// 
//   RooRealVar sigConstMeanP("sigConstMeanP","sigConstMeanP",17.6);
//   RooRealVar sigConstSigmaP("sigConstSigmaP","sigConstSigmaP",0.2);
//   RooRealVar sigConstMeanM("sigConstMeanM","sigConstMeanM",16.9);
//   RooRealVar sigConstSigmaM("sigConstSigmaM","sigConstSigmaM",0.5);
//   
// // RooRealVar f("f","f",0.5,0.,1.) ;
//   RooGaussian fconsta1p("fconsta1p","fconsta1p",*(qcdp.a1),a1ConstMeanP,a1ConstSigmaP);
// //   RooGaussian fconsta1m("fconsta1m","fconsta1m",*(qcdm.a1),a1ConstMeanM,a1ConstSigmaM);
//   
//   RooGaussian fconstsigp("fconstsigp","fconstsigp",*(qcdp.sigma),sigConstMeanP,sigConstSigmaP);
//   RooGaussian fconstsigm("fconstsigm","fconstsigm",*(qcdm.sigma),sigConstMeanM,sigConstSigmaM);
//  
//   RooProdPdf qcdpc("qcdpc","qcdpc",RooArgSet(*(qcdp.model),fconsta1p,fconstsigp));
//   RooProdPdf qcdmc("qcdmc","qcdmc",RooArgSet(*(qcdm.model),fconstsigm));
//   RooProdPdf qcdmc("qcdmc","qcdmc",RooArgSet(qcdp,fconstraint));

    RooGaussian constm("constm","constm",nEWKm,RooConst(hEWKMetm->Integral()),RooConst(0.15*hEWKMetm->Integral()));
    RooGaussian constp("constp","constp",nEWKp,RooConst(hEWKMetp->Integral()),RooConst(0.15*hEWKMetp->Integral()));
 
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK),   RooArgList(nSig,nEWK));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp),RooArgList(nSigp,nEWKp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm),RooArgList(nSigm,nEWKm));
  
//   // constrained background?
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,qcdpc),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,qcdmc),RooArgList(nSigm,nEWKm,nQCDm));
  
  
    
  
  // // Anti-Signal PDFs
  // RooDataHist awmunuMet ("awmunuMET", "awmunuMET", RooArgSet(pfmet),hAntiWmunuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,awmunuMet, 1);
  // RooDataHist awmunuMetp("awmunuMETp","awmunuMETp",RooArgSet(pfmet),hAntiWmunuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,awmunuMetp,1);
  // RooDataHist awmunuMetm("awmunuMETm","awmunuMETm",RooArgSet(pfmet),hAntiWmunuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,awmunuMetm,1); 
  
  // // Anti-EWK+top PDFs
  // RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  // RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  // RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  
// // //   // Anti-QCD Pdfs
  // CPepeModel2 aqcd("aqcd",pfmet,qcd.a1);
  // CPepeModel2 aqcdp("aqcdp",pfmet,qcdp.a1);
  // CPepeModel2 aqcdm("aqcdm",pfmet,qcdm.a1);

//   CPepeModel2 aqcd("aqcd",pfmet);
//   CPepeModel2 aqcdp("aqcdp",pfmet);
//   CPepeModel2 aqcdm("aqcdm",pfmet);
  
  // // Anti-selection PDFs
  // RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  // RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  // RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
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
  
  cout << "Starting values for Wmunu yields: " << endl;
  cout << "Selected: " << hDataMet->Integral() << endl;
  cout << "   sig: " << hWmunuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;
  cout << "   qcd: " << hDataMet->Integral()-hWmunuMet->Integral()-hEWKMet->Integral() << endl;

  cout << "Starting values for Wmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hDataMetp->Integral()-hWmunuMetp->Integral()-hEWKMetp->Integral() << endl;

  cout << "Starting values for Wmunu_m yields: " << endl;
  cout << "   sig: " << hWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  cout << "   qcd: " << hDataMetm->Integral()-hWmunuMetm->Integral()-hEWKMetm->Integral() << endl;
  
  
  // cout << "Starting values for AntiWmunu yields: " << endl;
  // cout << "Selected: " << hAntiDataMet->Integral() << endl;
  // cout << "   sig: " << hAntiWmunuMet->Integral() << endl;
  // cout << "   EWK: " << hAntiEWKMet->Integral() << endl;
  // cout << "   qcd: " << hAntiDataMet->Integral()-hAntiWmunuMet->Integral()-hAntiEWKMet->Integral() << endl;

  // cout << "Starting values for AntiWmunu_p yields: " << endl;
  // cout << "   sig: " << hWmunuMetp->Integral() << endl;
  // cout << "   EWK: " << hEWKMetp->Integral() << endl;
  // cout << "   qcd: " << hAntiDataMetp->Integral()-hAntiWmunuMetp->Integral()-hAntiEWKMetp->Integral() << endl;

  // cout << "Starting values for AntiWmunu_m yields: " << endl;
  // cout << "   sig: " << hAntiWmunuMetm->Integral() << endl;
  // cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;
  // cout << "   qcd: " << hAntiDataMetm->Integral()-hAntiWmunuMetm->Integral()-hAntiEWKMetm->Integral() << endl;

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());

  // RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  // RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  
  // RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  // RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  
  
//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",40000,0,100000);
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",40000,0,100000);

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);
  combine_workspace.import(pdfWm);
  combine_workspace.import(pdfWmp);
  combine_workspace.import(pdfWmm);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKm);

  char nname[100];
  sprintf(nname,"%s/Wmunu_pdfTemplates.root",CPlot::sOutDir.Data());
  combine_workspace.writeToFile(nname);
  RooFitResult *fitResp = pdfMetp.fitTo(dataMetp,Extended(),Minos(kTRUE),Save(kTRUE));
  RooFitResult *fitResm = 0;
  
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  std::cout << "Scale = " << (nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral() << std::endl;
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
   
  // TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  // hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  // TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  // hAntiMetDiff->SetMarkerStyle(kFullCircle);
  // hAntiMetDiff->SetMarkerSize(0.9);
   
  // TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
  // hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  // TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  // hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  // hAntiMetpDiff->SetMarkerSize(0.9);
    
  // TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
  // hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  // TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  // hAntiMetmDiff->SetMarkerStyle(kFullCircle); 
  // hAntiMetmDiff->SetMarkerSize(0.9);
   
  
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
  Int_t linecolorQCD = kViolet+2;
  Int_t fillcolorQCD = kViolet-5;
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
  
  TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  hDummyQCD->SetLineColor(linecolorQCD);
  hDummyQCD->SetFillColor(fillcolorQCD);
  hDummyQCD->SetFillStyle(1001);
   
   double yrange = 0.05;
   
   //
  // W MET plot
  //
  RooPlot *weframe = pfmet.frame(Bins(NBINS));
  weframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(weframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(weframe,LineColor(linecolorW));
  // pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  // pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
  // pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  // pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("wmunu_fitmet",weframe,"","mT [GeV]",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
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
  // pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  // pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
  // pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  // pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("wmunu_fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
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
  
  // // Anti-selection W+ background fits
  // RooPlot *awepframe = pfmet.frame(Bins(NBINS));    
  // awepframe->GetYaxis()->SetNdivisions(505);
  // awepframe->GetXaxis()->SetLabelOffset(2.0);
  // antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  // apdfMetp.plotOn(awepframe,FillColor(fillcolorW),DrawOption("F"));
  // apdfMetp.plotOn(awepframe,LineColor(linecolorW));
  // apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  // apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),LineColor(linecolorEWK));
  // apdfMetp.plotOn(awepframe,Components(RooArgSet(*(aqcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  // apdfMetp.plotOn(awepframe,Components(RooArgSet(*(aqcdp.model))),LineColor(linecolorQCD));
  // apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfWmp)),LineColor(linecolorW),LineStyle(2));
  // antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  // sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
  // CPlot plotAntiMetp("wmunu_fitantimetp",awepframe,"","",ylabel);
  // plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
  // plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  // plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  // plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  // plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  // plotAntiMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  // plotAntiMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  // plotAntiMetp.Draw(c,kFALSE,format,1);

  // CPlot plotAntiMetpDiff("wmunu_fitantimetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  // hAntiMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  // hAntiMetpDiff->GetYaxis()->SetLabelSize(0.11);
  // plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
  // plotAntiMetpDiff.SetYRange(-0.20,0.20);
  // plotAntiMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  // plotAntiMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  // plotAntiMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  // plotAntiMetpDiff.Draw(c,kTRUE,format,2);
  // plotAntiMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  // plotAntiMetp.SetName("wmunu_fitantimetplog");
  // plotAntiMetp.SetLogy();
  // plotAntiMetp.SetYRange(1e-3*(hAntiDataMetp->GetMaximum()),10*(hAntiDataMetp->GetMaximum()));
  // plotAntiMetp.Draw(c,kTRUE,format,1);
  // plotAntiMetp.Draw(c,kTRUE,"pdf",1);
  
  //
  // W- MET plot
  //
  RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
  wemframe->GetYaxis()->SetNdivisions(505);
  wemframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wemframe,LineColor(linecolorW));
  // pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  // pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
  // pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  // pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("wmunu_fitmetm",wemframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
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
  
  // // Anti-selection W- background fits
  // // fix these
  // RooPlot *awemframe = pfmet.frame(Bins(NBINS)); 
  // awemframe->GetYaxis()->SetNdivisions(505);
  // awemframe->GetXaxis()->SetLabelOffset(2.0);
  // antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  // apdfMetm.plotOn(awemframe,FillColor(fillcolorW),DrawOption("F"));
  // apdfMetm.plotOn(awemframe,LineColor(linecolorW));
  // apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  // apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),LineColor(linecolorEWK));
  // apdfMetm.plotOn(awemframe,Components(RooArgSet(*(aqcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  // apdfMetm.plotOn(awemframe,Components(RooArgSet(*(aqcdm.model))),LineColor(linecolorQCD));
  // apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfWmm)),LineColor(linecolorW),LineStyle(2));
  // antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  // sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetm->GetBinWidth(1));
  // CPlot plotAntiMetm("wmunu_fitantimetm",awemframe,"","",ylabel);
  // plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
  // plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  // plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  // plotAntiMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  // plotAntiMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  // plotAntiMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  // plotAntiMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  // plotAntiMetm.Draw(c,kFALSE,format,1);

  // CPlot plotAntiMetmDiff("wmunu_fitantimetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  // hAntiMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  // hAntiMetmDiff->GetYaxis()->SetLabelSize(0.11);
  // plotAntiMetmDiff.SetYRange(-0.2,0.2);
  // plotAntiMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  // plotAntiMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  // plotAntiMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  // plotAntiMetmDiff.AddHist1D(hAntiMetmDiff,"EX0",ratioColor);
  // plotAntiMetmDiff.Draw(c,kTRUE,format,2);
  // plotAntiMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  // plotAntiMetm.SetName("wmunu_fitantimetmlog");
  // plotAntiMetm.SetLogy();
  // plotAntiMetm.SetYRange(1e-3*(hAntiDataMetm->GetMaximum()),10*(hAntiDataMetm->GetMaximum()));
  // plotAntiMetm.Draw(c,kTRUE,format,1);
  // plotAntiMetm.Draw(c,kTRUE,"pdf",1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  
  ios_base::fmtflags flags;
  
  // Double_t chi2prob, chi2ndf;
  // Double_t ksprob, ksprobpe;
  
  // chi2prob = hDataMet->Chi2Test(hPdfMet,"PUW");
  // chi2ndf  = hDataMet->Chi2Test(hPdfMet,"CHI2/NDFUW");
  // ksprob   = hDataMet->KolmogorovTest(hPdfMet);
  // ksprobpe = hDataMet->KolmogorovTest(hPdfMet,"DX");
  // sprintf(txtfname,"%s/fitresWm.txt",CPlot::sOutDir.Data());
  // txtfile.open(txtfname);
  // assert(txtfile.is_open());
  
  // chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  // chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  // ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
  // ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
  // // sprintf(txtfname,"%s/fitresWmp.txt",CPlot::sOutDir.Data());
  // // txtfile.open(txtfname);
  // // assert(txtfile.is_open());
  
  // flags = txtfile.flags();
  // txtfile << setprecision(10);
  // txtfile << " *** Yields *** " << endl;
  // txtfile << "Selected: " << hDataMetp->Integral() << endl;
  // txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  // // txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  // txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  // // txtfile << "AntiSelected: " << hAntiDataMetp->Integral() << endl;
  // // txtfile << "  AntiSignal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResp) << endl;
  // // txtfile << "     AntiQCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResp) << endl;
  // // txtfile << "   AntiOther: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResp) << endl;
  // txtfile << endl; 
  // txtfile.flags(flags);
  
  // fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  // txtfile << endl;
  // printCorrelations(txtfile, fitResp);
  // txtfile << endl;
  // printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  
  // // fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  // // txtfile << endl;
  // // printCorrelations(txtfile, fitRes);
  // // txtfile << endl;
  // // printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  // txtfile.close();
  
  // chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  // chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  // ksprob   = hDataMetp->KolmogorovTest(hPdfMetp);
  // ksprobpe = hDataMetp->KolmogorovTest(hPdfMetp,"DX");  
  sprintf(txtfname,"%s/fitresWmp.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << hDataMetp->Integral() << endl;
  txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  // txtfile << "AntiSelected: " << hAntiDataMetp->Integral() << endl;
  // txtfile << "  AntiSignal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResp) << endl;
  // txtfile << "     AntiQCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResp) << endl;
  // txtfile << "   AntiOther: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResp) << endl;
  txtfile << endl; 
  txtfile.flags(flags);
  
  // fitResp->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  // txtfile << endl;
  // printCorrelations(txtfile, fitResp);
  // txtfile << endl;
  // printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();

  // chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  // chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  // ksprob   = hDataMetm->KolmogorovTest(hPdfMetm);
  // ksprobpe = hDataMetm->KolmogorovTest(hPdfMetm,"DX");  
  // sprintf(txtfname,"%s/fitresWmm.txt",CPlot::sOutDir.Data());
  // txtfile.open(txtfname);
  // assert(txtfile.is_open());
  
  // flags = txtfile.flags();
  // txtfile << setprecision(10);
  // txtfile << " *** Yields *** " << endl;
  // txtfile << "Selected: " << hDataMetm->Integral() << endl;
  // txtfile << "  Signal: " << nSigm.getVal() << " +/- " << nSigm.getPropagatedError(*fitResm) << endl;
  // txtfile << "     QCD: " << nQCDm.getVal() << " +/- " << nQCDm.getPropagatedError(*fitResm) << endl;
  // txtfile << "   Other: " << nEWKm.getVal() << " +/- " << nEWKm.getPropagatedError(*fitResm) << endl;
  // // txtfile << "AntiSelected: " << hAntiDataMetm->Integral() << endl;
  // // txtfile << "  Signal: " << nAntiSigm.getVal() << " +/- " << nAntiSigm.getPropagatedError(*fitResm) << endl;
  // // txtfile << "     QCD: " << nAntiQCDm.getVal() << " +/- " << nAntiQCDm.getPropagatedError(*fitResm) << endl;
  // // txtfile << "   Other: " << nAntiEWKm.getVal() << " +/- " << nAntiEWKm.getPropagatedError(*fitResm) << endl;
  // txtfile << endl;
  // txtfile.flags(flags);

  // fitResm->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  // txtfile << endl;
  // printCorrelations(txtfile, fitResm);
  // txtfile << endl;
  // printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  // txtfile.close();

  makeHTML(outputDir);
  
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
               
//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/WmunuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wmunu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmet.png\"><img src=\"fitmet.png\" alt=\"fitmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetp.png\"><img src=\"fitmetp.png\" alt=\"fitmetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetm.png\"><img src=\"fitmetm.png\" alt=\"fitmetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetlog.png\"><img src=\"fitmetlog.png\" alt=\"fitmetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetplog.png\"><img src=\"fitmetplog.png\" alt=\"fitmetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetmlog.png\"><img src=\"fitmetmlog.png\" alt=\"fitmetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimet.png\"><img src=\"fitantimet.png\" alt=\"fitantimet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetp.png\"><img src=\"fitantimetp.png\" alt=\"fitantimetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetm.png\"><img src=\"fitantimetm.png\" alt=\"fitantimetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetlog.png\"><img src=\"fitantimetlog.png\" alt=\"fitantimetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetplog.png\"><img src=\"fitantimetplog.png\" alt=\"fitantimetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitantimetmlog.png\"><img src=\"fitantimetmlog.png\" alt=\"fitantimetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}
