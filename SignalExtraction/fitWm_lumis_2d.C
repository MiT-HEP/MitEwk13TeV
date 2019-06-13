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
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"

// #include "ZBackgrounds.hh"

//helper class to handle rochester corrections
// #include <rochcor2015r.h>
// #include <muresolution_run2r.h>
#include <../RochesterCorr/RoccoR.cc>

// helper class to handle efficiency tables
#include "CEffUser1D.hh"
#include "CEffUser2D.hh"

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
#include "RooMinuit.h"
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
  bool doMTCut = false;
  bool doMET = true;
  bool doTemplate = false;
  // some flags to handle Recoil corrections
  bool doInclusive = true; // This should be the standard recoil correction: 3-gaussian
  // some flags for the recoil alternate models (not operational right now)
  bool doKeys = false; // RooKeysPDF instead of 3-Gaus
  bool doEta = false; // eta-binned 3-Gaus fit
  bool doDiago = false; // 
  // some flags to handle the pileup Up/Down systematics
  bool pileupUp = false;
  bool pileupDown = false;
  double yscale=0.5;
  
  // which MET type we use
  bool doPF = true;
  
  // Double_t vIsoBins[] = {0.0,0.15,0.45,0.55,0.65,0.75};
  Double_t vIsoBins[] = {0.0,0.2,0.25,0.35,0.45,0.55,0.65,0.75};
  // Double_t vIsoBins[] = {0.0,0.15,0.65,1.5};
  // Double_t vIsoBins[] = {0.0,0.25,0.35,0.45,0.55,0.65,0.75};
  // Double_t vIsoBins[] = {0.0,0.25,0.35,0.45,0.55};
  int nIsoBins = sizeof(vIsoBins)/sizeof(vIsoBins[0])-1;
  std::cout << "size of isobin array is " << nIsoBins << std::endl;
  
  
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
  
  
  // MET histogram binning and range
  // const Int_t    NBINS   = 75*4;
  // const Int_t    NBINS   = 125;
  // const Int_t    NBINS   = 50;
  // const Int_t    NBINS   = 100;
  const Int_t    NBINS   = 75;
  const Double_t METMIN  = 0;
  // const Double_t METMIN  = 25;
  // const Double_t METMIN  = 75;
  const Double_t METMAX  = 150;
  
  const Double_t PT_CUT  = 25;
  // const Double_t PT_CUT  = 35;
  const Double_t ETA_CUT = 2.4;
  // const Double_t ETA_CUT = 1.4;
  
  const Double_t mu_MASS = 0.1057;
  const int NTOYS = 100;
  const Double_t MT_CUT = 50.0;
    // efficiency files
 
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

  // // efficiency files 2Bins
  
  const TString dataHLTEff2BinName_pos = baseDir + "Data/MuHLTEff_aMCxPythia/Positive/eff.root";
  const TString dataHLTEff2BinName_neg = baseDir + "Data/MuHLTEff_aMCxPythia/Negative/eff.root";
  const TString zmmHLTEff2BinName_pos  = baseDir + "MC/MuHLTEff_aMCxPythia/Positive/eff.root";
  const TString zmmHLTEff2BinName_neg  = baseDir + "MC/MuHLTEff_aMCxPythia/Negative/eff.root";

  const TString dataSelEff2BinName_pos = baseDir + "Data/MuSITEff_aMCxPythia/Positive/eff.root";
  const TString dataSelEff2BinName_neg = baseDir + "Data/MuSITEff_aMCxPythia/Negative/eff.root";
  const TString zmmSelEff2BinName_pos  = baseDir + "MC/MuSITEff_aMCxPythia/Positive/eff.root";
  const TString zmmSelEff2BinName_neg  = baseDir + "MC/MuSITEff_aMCxPythia/Negative/eff.root";


  const TString dataStaEff2BinName_pos = baseDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString dataStaEff2BinName_neg = baseDir + "Data/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString zmmStaEff2BinName_pos  = baseDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";
  const TString zmmStaEff2BinName_neg  = baseDir + "MC/MuStaEff_aMCxPythia/Combined/eff.root";

  TString StaEffSignalShapeSys     = baseDir + "Results/MuStaSigSys.root";
  TString StaEffBackgroundShapeSys = baseDir + "Results/MuStaBkgSys.root";
  TString SelEffSignalShapeSys     = baseDir + "Results/MuSITSigSys.root";
  TString SelEffBackgroundShapeSys = baseDir + "Results/MuSITBkgSys.root";

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);


  

  // file format for output plots
  const TString format("png"); 


  // ===================== Recoil correction files ============================
  const TString directory1("/afs/cern.ch/user/d/dalfonso/public/WZ/JUNE25");
  // const TString directory("/afs/cern.ch/user/d/dalfonso/public/WZ/JULY5");
  const TString directory2("/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Recoil");
  //  const TString directory1("/afs/cern.ch/user/d/dalfonso/public/WZ/dec5");
  const TString directory("/afs/cern.ch/user/d/dalfonso/public/WZ/JULY5");

  // for PF, 13 TeV Low PU, inclusive
  RecoilCorrector *recoilCorr = new  RecoilCorrector("","");

  
  RecoilCorrector *recoilCorrm = new  RecoilCorrector("","");

  // for Puppi, inclusive
  // RecoilCorrector *recoilCorr = new  RecoilCorrector("","");
  
  if(doDiago){
    recoilCorr->loadRooWorkspacesDiagMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory2.Data()),1);
    recoilCorr->loadRooWorkspacesDiagData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory2.Data()),1);
    recoilCorr->loadRooWorkspacesDiagMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory2.Data()),1);
  } else {
    recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_WmpMCPF/",directory2.Data()));
    recoilCorr->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory2.Data()));
    recoilCorr->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory2.Data()));
  }
  // // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/JUNE25_ZmmMCPuppi/",directory1.Data()));
  // // // recoilCorr->loadRooWorkspacesData(Form("%s/JUNE25_ZmmDataPuppi_bkg/",directory1.Data()));
  // // // recoilCorr->loadRooWorkspacesMC(Form("%s/JUNE25_ZmmMCPuppi/",directory1.Data()));
  // // // // // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/APR4_ZmmMCPuppi/",directory1.Data()));

  // // make the diagonalized and central consistent
  // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/JUNE25_ZmmMCPuppi/",directory1.Data()));
  // //recoilCorr->loadRooWorkspacesData(Form("%s/JUNE25_ZmmDataPuppi_bkg/",directory1.Data()));
  // //recoilCorr->loadRooWorkspacesMC(Form("%s/JUNE25_ZmmMCPuppi/",directory1.Data()));
  
    // // // // // // recoilCorr->loadRooWorkspacesData(Form("%s/ZeeDataPuppi_bkgTopEWK/",directory.Data()));
  // // // the ones i'm actually using
    // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_genFix/",directory.Data()));
   // recoilCorr->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
    // recoilCorr->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  // // // // // // // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  // // just test these
      // // recoilCorr->loadRooWorkspacesMCtoCorrect(Form("%s/ZeeMCPuppi/",directory.Data()));
  // // recoilCorr->loadRooWorkspacesData(Form("%s/ZeeDataPuppi_bkgTopEWK//",directory.Data()));
    // // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  // // // // // recoilCorr->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  
  // }
  
  if(doDiago){
    recoilCorrm->loadRooWorkspacesDiagMCtoCorrect(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory2.Data()),1);
    recoilCorrm->loadRooWorkspacesDiagData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory2.Data()),1);
    recoilCorrm->loadRooWorkspacesDiagMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory2.Data()),1);
  }else {  
    // recoilCorrm->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_WmmMCPF_mT50_v1/",directory2.Data()));
    recoilCorrm->loadRooWorkspacesMCtoCorrect(Form("%s/LowPU2017ID_13TeV_WmmMCPF/",directory2.Data()));
    recoilCorrm->loadRooWorkspacesData(Form("%s/LowPU2017ID_13TeV_ZmmDataPF/",directory2.Data()));
    recoilCorrm->loadRooWorkspacesMC(Form("%s/LowPU2017ID_13TeV_ZmmMCPF/",directory2.Data()));
  }
  
    // // // // The ones I'm actually using:
  // recoilCorrm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_genFix/",directory.Data()));
   // recoilCorrm->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
   // recoilCorrm->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
    // // // just test
  // // recoilCorrm->loadRooWorkspacesMCtoCorrect(Form("%s/ZeeMCPuppi/",directory.Data()));
    // // recoilCorrm->loadRooWorkspacesData(Form("%s/ZeeDataPuppi_bkgTopEWK/",directory.Data()));
  // // recoilCorrm->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  
    // // // // recoilCorrm->loadRooWorkspacesMC(Form("%s/ZeeMCPuppi/",directory.Data()));
  // // // // recoilCorrm->loadRooWorkspacesData(Form("%s/ZeeDataPuppi_bkgTopEWK/",directory.Data()));
  
  // }
  
// //   // --------------------- Eta-binned recoil corrections -----------------------
  // RecoilCorrector *recoilCorr05 = new  RecoilCorrector("","");
  // recoilCorr05->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_rap05_genFix/",directory.Data()));
  // // recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  // recoilCorr05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05/",directory.Data()));
  // recoilCorr05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));
  
  // RecoilCorrector *recoilCorrm05 = new  RecoilCorrector("","");
  // recoilCorrm05->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_rap05_genFix/",directory.Data()));
  // // recoilCorrm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
  // recoilCorrm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05/",directory.Data()));
  // recoilCorrm05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));

  // RecoilCorrector *recoilCorr051 = new  RecoilCorrector("","");
  // recoilCorr051->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_rap05-1_genFix/",directory.Data()));
  // // recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  // recoilCorr051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05-1/",directory.Data()));
  // recoilCorr051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));

  // RecoilCorrector *recoilCorrm051 = new  RecoilCorrector("","");
  // recoilCorrm051->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_rap05-1_genFix/",directory.Data()));
  // // recoilCorrm051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
  // recoilCorrm051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap05-1/",directory.Data()));
  // recoilCorrm051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));

  // RecoilCorrector *recoilCorr1 = new  RecoilCorrector("","");
  // recoilCorr1->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPuppi_rap1_genFix/",directory.Data()));
  // // recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  // recoilCorr1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap1/",directory.Data()));
  // recoilCorr1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data()));

  // RecoilCorrector *recoilCorrm1 = new  RecoilCorrector("","");
  // recoilCorrm1->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPuppi_rap1_genFix/",directory.Data()));
  // // recoilCorrm1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
  // recoilCorrm1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK_rap1/",directory.Data()));
  // recoilCorrm1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data())); 
// //   
// // //   // ---------------------- KEYS - inclusive-------------
  // std::cout << "load inclusive RooKeys files" << std::endl;
  // RecoilCorrector *recoilCorrKeys = new  RecoilCorrector("","");
  // // recoilCorrKeys->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys/",directory.Data()));
  // recoilCorrKeys->loadRooWorkspacesMCtoCorrectKeys(Form("%s/ZmmMCPuppi_keys/",directory.Data()));
  // recoilCorrKeys->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
  // recoilCorrKeys->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));
  
  // RecoilCorrector *recoilCorrKeysm = new  RecoilCorrector("","");
  // // recoilCorrKeysm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys/",directory.Data()));
  // recoilCorrKeysm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/ZmmMCPuppi_keys/",directory.Data()));
  // recoilCorrKeysm->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkgTopEWK/",directory.Data()));
  // recoilCorrKeysm->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi/",directory.Data()));


// //   // ---------------------- KEYS -------------
// //   RecoilCorrector *recoilCorrKeys05 = new  RecoilCorrector("","");
// //   recoilCorrKeys05->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys_rap05/",directory.Data()));
// //   recoilCorrKeys05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
// //   recoilCorrKeys05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));
// //   
// //   RecoilCorrector *recoilCorrKeysm05 = new  RecoilCorrector("","");
// //   recoilCorrKeysm05->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys_rap05/",directory.Data()));
// //   recoilCorrKeysm05->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05/",directory.Data()));
// //   recoilCorrKeysm05->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05/",directory.Data()));
// // 
// //   RecoilCorrector *recoilCorrKeys051 = new  RecoilCorrector("","");
// //   recoilCorrKeys051->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys_rap05-1/",directory.Data()));
// //   recoilCorrKeys051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
// //   recoilCorrKeys051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));
// // 
// //   RecoilCorrector *recoilCorrKeysm051 = new  RecoilCorrector("","");
// //   recoilCorrKeysm051->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys_rap05-1/",directory.Data()));
// //   recoilCorrKeysm051->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap05-1/",directory.Data()));
// //   recoilCorrKeysm051->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap05-1/",directory.Data()));
// // 
// //   RecoilCorrector *recoilCorrKeys1 = new  RecoilCorrector("","");
// //   recoilCorrKeys1->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmpMCPuppi_keys_rap1/",directory.Data()));
// //   recoilCorrKeys1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
// //   recoilCorrKeys1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data()));
// // 
// //   RecoilCorrector *recoilCorrKeysm1 = new  RecoilCorrector("","");
// //   recoilCorrKeysm1->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WmmMCPuppi_keys_rap1/",directory.Data()));
// //   recoilCorrKeysm1->loadRooWorkspacesData(Form("%s/ZmmDataPuppi_bkg_rap1/",directory.Data()));
// //   recoilCorrKeysm1->loadRooWorkspacesMC(Form("%s/ZmmMCPuppi_rap1/",directory.Data())); 
  
  // // ==========================================================================
  // // ---------------- Recoil corrections for PF MET ---------------------
  // const TString directoryPF("/afs/cern.ch/user/d/dalfonso/public/WZ/jan23");
  // RecoilCorrector *recoilCorrPF = new  RecoilCorrector("","");
  // recoilCorrPF->loadRooWorkspacesMCtoCorrect(Form("%s/WmpMCPF/",directoryPF.Data()));
  // recoilCorrPF->loadRooWorkspacesData(Form("%s/ZmmDataPF_bkg/",directoryPF.Data()));
  // recoilCorrPF->loadRooWorkspacesMC(Form("%s/ZmmMCPF/",directoryPF.Data()));
  
  // RecoilCorrector *recoilCorrPFm = new  RecoilCorrector("","");
  // recoilCorrPFm->loadRooWorkspacesMCtoCorrect(Form("%s/WmmMCPF/",directoryPF.Data()));
  // recoilCorrPFm->loadRooWorkspacesData(Form("%s/ZmmDataPF_bkg/",directoryPF.Data()));
  // recoilCorrPFm->loadRooWorkspacesMC(Form("%s/ZmmMCPF/",directoryPF.Data()));
  // // ==========================================================================
  
  // ----------------------------------------------------
  // Load the plots of relative difference, to create the up/down shapes
  TFile *_rdWmp = new TFile("shapeDiff/Wmp_relDiff.root");
  TFile *_rdWmm = new TFile("shapeDiff/Wmm_relDiff.root");
  
  TH1D *hh_diffm = new TH1D("hh_diffm","hh_diffm",75,0,150);
  TH1D *hh_diffp = new TH1D("hh_diffp","hh_diffp",75,0,150);
  
  hh_diffm = (TH1D*)_rdWmm->Get("hh_diff");
  hh_diffp = (TH1D*)_rdWmp->Get("hh_diff");
  // -----------------------------------------------------
  
  // Load the Z data and Z MC Pt spectra
  TFile *_rat1 = new TFile("shapeDiff/zmm_PDFUnc.root");
  TH1D *hh_mc;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_mc = (TH1D*)_rat1->Get("hZPtTruthNominal"); 
  hh_mc->Scale(1/hh_mc->Integral()); // normalize
  
  TFile *_rat2 = new TFile("shapeDiff/UnfoldingOutputZPt.root");
  TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_diff = (TH1D*)_rat2->Get("hUnfold");
  hh_diff->Scale(1/hh_diff->Integral()); // normalize
  hh_diff->Divide(hh_mc);
   
  //
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
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet   = new TH1D("hDataMet","",  NBINS,METMIN,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm  = new TH1D("hDataMetm","", NBINS,METMIN,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp  = new TH1D("hDataMetp","", NBINS,METMIN,METMAX); hDataMetp->Sumw2();
  TH1D *hWmunuMet  = new TH1D("hWmunuMet","", NBINS,METMIN,METMAX); hWmunuMet->Sumw2();
  TH1D *hWmunuMetp = new TH1D("hWmunuMetp","",NBINS,METMIN,METMAX); hWmunuMetp->Sumw2();
  TH1D *hWmunuMetm = new TH1D("hWmunuMetm","",NBINS,METMIN,METMAX); hWmunuMetm->Sumw2();
  
  TH1D *hQCDMet    = new TH1D("hQCDMet", "",  NBINS,METMIN,METMAX); hQCDMet->Sumw2();
  TH1D *hQCDMetp   = new TH1D("hQCDMetp", "", NBINS,METMIN,METMAX); hQCDMetp->Sumw2();
  TH1D *hQCDMetm   = new TH1D("hQCDMetm", "", NBINS,METMIN,METMAX); hQCDMetm->Sumw2();
  
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,METMIN,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,METMIN,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,METMIN,METMAX); hEWKMetm->Sumw2();
  
  TH1D *hEWKMet_PileupUp   = new TH1D("hEWKMet_PileupUp",  "",NBINS,METMIN,METMAX); hEWKMet_PileupUp->Sumw2();
  TH1D *hEWKMetp_PileupUp  = new TH1D("hEWKMetp_PileupUp", "",NBINS,METMIN,METMAX); hEWKMetp_PileupUp->Sumw2();
  TH1D *hEWKMetm_PileupUp  = new TH1D("hEWKMetm_PileupUp", "",NBINS,METMIN,METMAX); hEWKMetm_PileupUp->Sumw2();
  
  TH1D *hEWKMet_PileupDown   = new TH1D("hEWKMet_PileupDown",  "",NBINS,METMIN,METMAX); hEWKMet_PileupDown->Sumw2();
  TH1D *hEWKMetp_PileupDown  = new TH1D("hEWKMetp_PileupDown", "",NBINS,METMIN,METMAX); hEWKMetp_PileupDown->Sumw2();
  TH1D *hEWKMetm_PileupDown  = new TH1D("hEWKMetm_PileupDown", "",NBINS,METMIN,METMAX); hEWKMetm_PileupDown->Sumw2();
  TH1D *hWmunuMet_RecoilUp  = new TH1D("hWmunuMet_RecoilUp", "",NBINS,METMIN,METMAX); hWmunuMet_RecoilUp->Sumw2();
  TH1D *hWmunuMetp_RecoilUp = new TH1D("hWmunuMetp_RecoilUp","",NBINS,METMIN,METMAX); hWmunuMetp_RecoilUp->Sumw2();
  TH1D *hWmunuMetm_RecoilUp = new TH1D("hWmunuMetm_RecoilUp","",NBINS,METMIN,METMAX); hWmunuMetm_RecoilUp->Sumw2();
  TH1D *hWmunuMet_RecoilDown  = new TH1D("hWmunuMet_RecoilDown", "",NBINS,METMIN,METMAX); hWmunuMet_RecoilDown->Sumw2();
  TH1D *hWmunuMetp_RecoilDown = new TH1D("hWmunuMetp_RecoilDown","",NBINS,METMIN,METMAX); hWmunuMetp_RecoilDown->Sumw2();
  TH1D *hWmunuMetm_RecoilDown = new TH1D("hWmunuMetm_RecoilDown","",NBINS,METMIN,METMAX); hWmunuMetm_RecoilDown->Sumw2();

  TH1D *hWmunuMet_RecoilCUp  = new TH1D("hWmunuMet_RecoilCUp", "",NBINS,METMIN,METMAX); hWmunuMet_RecoilCUp->Sumw2();
  TH1D *hWmunuMetp_RecoilCUp = new TH1D("hWmunuMetp_RecoilCUp","",NBINS,METMIN,METMAX); hWmunuMetp_RecoilCUp->Sumw2();
  TH1D *hWmunuMetm_RecoilCUp = new TH1D("hWmunuMetm_RecoilCUp","",NBINS,METMIN,METMAX); hWmunuMetm_RecoilCUp->Sumw2();
  TH1D *hWmunuMet_RecoilCDown  = new TH1D("hWmunuMet_RecoilCDown", "",NBINS,METMIN,METMAX); hWmunuMet_RecoilCDown->Sumw2();
  TH1D *hWmunuMetp_RecoilCDown = new TH1D("hWmunuMetp_RecoilCDown","",NBINS,METMIN,METMAX); hWmunuMetp_RecoilCDown->Sumw2();
  TH1D *hWmunuMetm_RecoilCDown = new TH1D("hWmunuMetm_RecoilCDown","",NBINS,METMIN,METMAX); hWmunuMetm_RecoilCDown->Sumw2();

  TH1D *hWmunuMet_ScaleUp  = new TH1D("hWmunuMet_ScaleUp", "",NBINS,METMIN,METMAX); hWmunuMet_ScaleUp->Sumw2();
  TH1D *hWmunuMetp_ScaleUp = new TH1D("hWmunuMetp_ScaleUp","",NBINS,METMIN,METMAX); hWmunuMetp_ScaleUp->Sumw2();
  TH1D *hWmunuMetm_ScaleUp = new TH1D("hWmunuMetm_ScaleUp","",NBINS,METMIN,METMAX); hWmunuMetm_ScaleUp->Sumw2();
  TH1D *hWmunuMet_ScaleDown  = new TH1D("hWmunuMet_ScaleDown", "",NBINS,METMIN,METMAX); hWmunuMet_ScaleDown->Sumw2();
  TH1D *hWmunuMetp_ScaleDown = new TH1D("hWmunuMetp_ScaleDown","",NBINS,METMIN,METMAX); hWmunuMetp_ScaleDown->Sumw2();
  TH1D *hWmunuMetm_ScaleDown = new TH1D("hWmunuMetm_ScaleDown","",NBINS,METMIN,METMAX); hWmunuMetm_ScaleDown->Sumw2();
  
  TH1D *hWmunuMet_PileupUp  = new TH1D("hWmunuMet_PileupUp", "",NBINS,METMIN,METMAX); hWmunuMet_PileupUp->Sumw2();
  TH1D *hWmunuMetp_PileupUp = new TH1D("hWmunuMetp_PileupUp","",NBINS,METMIN,METMAX); hWmunuMetp_PileupUp->Sumw2();
  TH1D *hWmunuMetm_PileupUp = new TH1D("hWmunuMetm_PileupUp","",NBINS,METMIN,METMAX); hWmunuMetm_PileupUp->Sumw2();
  TH1D *hWmunuMet_PileupDown  = new TH1D("hWmunuMet_PileupDown", "",NBINS,METMIN,METMAX); hWmunuMet_PileupDown->Sumw2();
  TH1D *hWmunuMetp_PileupDown = new TH1D("hWmunuMetp_PileupDown","",NBINS,METMIN,METMAX); hWmunuMetp_PileupDown->Sumw2();
  TH1D *hWmunuMetm_PileupDown = new TH1D("hWmunuMetm_PileupDown","",NBINS,METMIN,METMAX); hWmunuMetm_PileupDown->Sumw2();

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet","",  NBINS,METMIN,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,METMIN,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,METMIN,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWmunuMet  = new TH1D("hAntiWmunuMet","", NBINS,METMIN,METMAX); hAntiWmunuMet->Sumw2();
  TH1D *hAntiWmunuMetp = new TH1D("hAntiWmunuMetp","",NBINS,METMIN,METMAX); hAntiWmunuMetp->Sumw2();
  TH1D *hAntiWmunuMetm = new TH1D("hAntiWmunuMetm","",NBINS,METMIN,METMAX); hAntiWmunuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet", "",  NBINS,METMIN,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,METMIN,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,METMIN,METMAX); hAntiEWKMetm->Sumw2();
  
  TH1D *hAntiQCDMet    = new TH1D("hAntiQCDMet", "",  NBINS,METMIN,METMAX); hAntiQCDMet->Sumw2();
  TH1D *hAntiQCDMetp   = new TH1D("hAntiQCDMetp", "", NBINS,METMIN,METMAX); hAntiQCDMetp->Sumw2();
  TH1D *hAntiQCDMetm   = new TH1D("hAntiQCDMetm", "", NBINS,METMIN,METMAX); hAntiQCDMetm->Sumw2();
  
  TH1D *hMuonEtaDatap = new TH1D("hMuonEtaDatap","",25,0,2.5); hMuonEtaDatap->Sumw2();
  TH1D *hMuonEtaDatam = new TH1D("hMuonEtaDatam","",25,0,2.5); hMuonEtaDatam->Sumw2();
  TH1D *hMuonEtaMCp = new TH1D("hMuonEtaMCp","",25,0,2.5); hMuonEtaMCp->Sumw2();
  TH1D *hMuonEtaMCm = new TH1D("hMuonEtaMCm","",25,0,2.5); hMuonEtaMCm->Sumw2();
  TH1D *hMuonEtaAntiDatap = new TH1D("hMuonEtaAntiDatap","",25,0,2.5); hMuonEtaAntiDatap->Sumw2();
  TH1D *hMuonEtaAntiDatam = new TH1D("hMuonEtaAntiDatam","",25,0,2.5); hMuonEtaAntiDatam->Sumw2();
  
  TH1D *hDataMetpPhi   = new TH1D("hDataMetpPhi","",  100,-3.15, 6.30); hDataMetpPhi->Sumw2();
  TH1D *hDataMetmPhi   = new TH1D("hDataMetmPhi","",  100,-3.15, 6.30); hDataMetmPhi->Sumw2();
  TH1D *hWmunuMetpPhi   = new TH1D("hWmunuMetpPhi","",  100,-3.15, 6.30); hWmunuMetpPhi->Sumw2();
  TH1D *hWmunuMetmPhi   = new TH1D("hWmunuMetmPhi","",  100,-3.15, 6.30); hWmunuMetmPhi->Sumw2();
  
  
  
  // For the background isolation bins
  // iso 0.25-0.35, 0.35-0.45, 0.45-0.55, 0.55-0.65. 0.65-0.75 (5 bins total)
  // TH1F **hM2M0_HB_5_10 = new TH1F*[5]; // example
  // TH1D **hAntiDataMetIsoBins   = new TH1D*[5];// hAntiDataMet->Sumw2();
  TH1D **hDataMetmIsoBins  = new TH1D*[nIsoBins];// hAntiDataMetm->Sumw2();  
  TH1D **hDataMetpIsoBins  = new TH1D*[nIsoBins];// hAntiDataMetp->Sumw2();
  // TH1D **hAntiWmunuMetIsoBins   = new TH1D*[5];// hAntiWmunuMet->Sumw2();
  TH1D **hWmunuMetpIsoBins  = new TH1D*[nIsoBins];// hAntiWmunuMetp->Sumw2();
  TH1D **hWmunuMetmIsoBins  = new TH1D*[nIsoBins];// hAntiWmunuMetm->Sumw2();
  // TH1D **hAntiEWKMetIsoBins    = new TH1D*[5];// hAntiEWKMet->Sumw2();
  TH1D **hEWKMetpIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hEWKMetmIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hDibMetpIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hDibMetmIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hTtbMetpIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hTtbMetmIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hWxMetpIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hWxMetmIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hZxxMetpIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hZxxMetmIsoBins   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hMetpIsoValues = new TH1D*[nIsoBins];
  TH1D **hMetmIsoValues = new TH1D*[nIsoBins];
  // Create a histogram pointer in each space in the array
  for(int i = 0; i < nIsoBins; i++){
    hDataMetmIsoBins[i]  = new TH1D(("hDataMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDataMetpIsoBins[i]  = new TH1D(("hDataMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWmunuMetpIsoBins[i]  = new TH1D(("hWmunuMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWmunuMetmIsoBins[i]  = new TH1D(("hWmunuMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hEWKMetpIsoBins[i]  = new TH1D(("hEwkMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hEWKMetmIsoBins[i]  = new TH1D(("hEwkMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	hDibMetpIsoBins[i]  = new TH1D(("hDibMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDibMetmIsoBins[i]  = new TH1D(("hDibMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	hTtbMetpIsoBins[i]  = new TH1D(("hTtbMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hTtbMetmIsoBins[i]  = new TH1D(("hTtbMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	hWxMetpIsoBins[i]  = new TH1D(("hWxMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWxMetmIsoBins[i]  = new TH1D(("hWxMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
	hZxxMetpIsoBins[i]  = new TH1D(("hZxxMetpIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hZxxMetmIsoBins[i]  = new TH1D(("hZxxMetmIsoBin_"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    //FileService->make<TH1F>(("hM2M0_HB_5_10_p"+int2string(i+1)).c_str(),  ("M2/M0, HB part"+int2string(i+1)+" 5-10GeV").c_str(),  50, 0.0, 2.0);
    if(i==0){
     hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  1000,vIsoBins[i],0.15);
     hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  1000,vIsoBins[i],0.15);
    }else{
	 hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
     hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
	}
  }
  
  
//    //
//   // HLT efficiency
//   //
  cout << "Loading trigger efficiencies..." << endl;
  
  TFile *dataHLTEffFile_pos = new TFile(dataHLTEffName_pos);
  CEffUser2D dataHLTEff_pos;
  dataHLTEff_pos.loadEff((TH2D*)dataHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEffFile_neg = new TFile(dataHLTEffName_neg);
  CEffUser2D dataHLTEff_neg;
  dataHLTEff_neg.loadEff((TH2D*)dataHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEffFile_neg->Get("hErrhEtaPt"));
    
  TFile *zmmHLTEffFile_pos = new TFile(zmmHLTEffName_pos);
  CEffUser2D zmmHLTEff_pos;
  zmmHLTEff_pos.loadEff((TH2D*)zmmHLTEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEffFile_neg = new TFile(zmmHLTEffName_neg);
  CEffUser2D zmmHLTEff_neg;
  zmmHLTEff_neg.loadEff((TH2D*)zmmHLTEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataHLTEff2BinFile_pos = new TFile(dataHLTEff2BinName_pos);
  CEffUser2D dataHLTEff2Bin_pos;
  dataHLTEff2Bin_pos.loadEff((TH2D*)dataHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataHLTEff2BinFile_neg = new TFile(dataHLTEff2BinName_neg);
  CEffUser2D dataHLTEff2Bin_neg;
  dataHLTEff2Bin_neg.loadEff((TH2D*)dataHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataHLTEff2BinFile_neg->Get("hErrhEtaPt"));
    
  TFile *zmmHLTEff2BinFile_pos = new TFile(zmmHLTEff2BinName_pos);
  CEffUser2D zmmHLTEff2Bin_pos;
  zmmHLTEff2Bin_pos.loadEff((TH2D*)zmmHLTEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmHLTEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmHLTEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmHLTEff2BinFile_neg = new TFile(zmmHLTEff2BinName_neg);
  CEffUser2D zmmHLTEff2Bin_neg;
  zmmHLTEff2Bin_neg.loadEff((TH2D*)zmmHLTEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmHLTEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmHLTEff2BinFile_neg->Get("hErrhEtaPt"));

 
  
  //
  // Selection efficiency
  //
  cout << "Loading selection efficiencies..." << endl;
  
  TFile *dataSelEffFile_pos = new TFile(dataSelEffName_pos);
  CEffUser2D dataSelEff_pos;
  dataSelEff_pos.loadEff((TH2D*)dataSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEffFile_neg = new TFile(dataSelEffName_neg);
  CEffUser2D dataSelEff_neg;
  dataSelEff_neg.loadEff((TH2D*)dataSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEffFile_pos = new TFile(zmmSelEffName_pos);
  CEffUser2D zmmSelEff_pos;
  zmmSelEff_pos.loadEff((TH2D*)zmmSelEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEffFile_neg = new TFile(zmmSelEffName_neg);
  CEffUser2D zmmSelEff_neg;
  zmmSelEff_neg.loadEff((TH2D*)zmmSelEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataSelEff2BinFile_pos = new TFile(dataSelEff2BinName_pos);
  CEffUser2D dataSelEff2Bin_pos;
  dataSelEff2Bin_pos.loadEff((TH2D*)dataSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataSelEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataSelEff2BinFile_neg = new TFile(dataSelEff2BinName_neg);
  CEffUser2D dataSelEff2Bin_neg;
  dataSelEff2Bin_neg.loadEff((TH2D*)dataSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataSelEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmSelEff2BinFile_pos = new TFile(zmmSelEff2BinName_pos);
  CEffUser2D zmmSelEff2Bin_pos;
  zmmSelEff2Bin_pos.loadEff((TH2D*)zmmSelEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmSelEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmSelEff2BinFile_pos->Get("hErrhEtaPt"));

  TFile *zmmSelEff2BinFile_neg = new TFile(zmmSelEff2BinName_neg);
  CEffUser2D zmmSelEff2Bin_neg;
  zmmSelEff2Bin_neg.loadEff((TH2D*)zmmSelEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmSelEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmSelEff2BinFile_neg->Get("hErrhEtaPt"));

  
    //
  // Standalone efficiency
  //
  cout << "Loading standalone efficiencies..." << endl;
  
  TFile *dataStaEffFile_pos = new TFile(dataStaEffName_pos);
  CEffUser2D dataStaEff_pos;
  dataStaEff_pos.loadEff((TH2D*)dataStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEffFile_neg = new TFile(dataStaEffName_neg);
  CEffUser2D dataStaEff_neg;
  dataStaEff_neg.loadEff((TH2D*)dataStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEffFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_pos = new TFile(zmmStaEffName_pos);
  CEffUser2D zmmStaEff_pos;
  zmmStaEff_pos.loadEff((TH2D*)zmmStaEffFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEffFile_neg = new TFile(zmmStaEffName_neg);
  CEffUser2D zmmStaEff_neg;
  zmmStaEff_neg.loadEff((TH2D*)zmmStaEffFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEffFile_neg->Get("hErrhEtaPt"));

  TFile *dataStaEff2BinFile_pos = new TFile(dataStaEff2BinName_pos);
  CEffUser2D dataStaEff2Bin_pos;
  dataStaEff2Bin_pos.loadEff((TH2D*)dataStaEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)dataStaEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)dataStaEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *dataStaEff2BinFile_neg = new TFile(dataStaEff2BinName_neg);
  CEffUser2D dataStaEff2Bin_neg;
  dataStaEff2Bin_neg.loadEff((TH2D*)dataStaEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)dataStaEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)dataStaEff2BinFile_neg->Get("hErrhEtaPt"));
  
  TFile *zmmStaEff2BinFile_pos = new TFile(zmmStaEff2BinName_pos);
  CEffUser2D zmmStaEff2Bin_pos;
  zmmStaEff2Bin_pos.loadEff((TH2D*)zmmStaEff2BinFile_pos->Get("hEffEtaPt"), (TH2D*)zmmStaEff2BinFile_pos->Get("hErrlEtaPt"), (TH2D*)zmmStaEff2BinFile_pos->Get("hErrhEtaPt"));
  
  TFile *zmmStaEff2BinFile_neg = new TFile(zmmStaEff2BinName_neg);
  CEffUser2D zmmStaEff2Bin_neg;
  zmmStaEff2Bin_neg.loadEff((TH2D*)zmmStaEff2BinFile_neg->Get("hEffEtaPt"), (TH2D*)zmmStaEff2BinFile_neg->Get("hErrlEtaPt"), (TH2D*)zmmStaEff2BinFile_neg->Get("hErrhEtaPt"));
 


  double tolerance = ROOT::Math::MinimizerOptions::DefaultTolerance();
  string algo = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
  string type = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
  // ROOT::Math::MinimizerOptions::SetStrategy(2);
  int strategy= ROOT::Math::MinimizerOptions::DefaultStrategy();
  // int strategy= ROOT::Math::MinimizerOptions::Strategy();

  int precision= ROOT::Math::MinimizerOptions::DefaultPrecision();
  int MaxFunctionCalls= ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls();
  int MaxIterations= ROOT::Math::MinimizerOptions::DefaultMaxIterations();

  cout << "DEFAULTS: algo " << algo.c_str() << " type " << type.c_str() << " tolerance " << tolerance << " strategy " << strategy << " precision " << precision << " MaxIterations " << MaxIterations << " MaxFunctionCalls " << MaxFunctionCalls << endl;
  

//   
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0, *genV=0, *genLep=0;
  Float_t pfChIso, pfGamIso, pfNeuIso;
    Float_t pfCombIso;
    
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
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireWeight", &prefireWeight);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress(met_name.c_str(),      &met);       // MET
    intree->SetBranchAddress(metPhi_name.c_str(),   &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",       &mt);        // transverse mass
    intree->SetBranchAddress(u1_name.c_str(),       &u1);        // parallel component of recoil
    intree->SetBranchAddress(u2_name.c_str(),       &u2);        // perpendicular component of recoil
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("genLep",      &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV",     &genV);       // lepton 4-vector
    intree->SetBranchAddress("pfChIso",  &pfChIso);
    intree->SetBranchAddress("pfGamIso", &pfGamIso);
    intree->SetBranchAddress("pfNeuIso", &pfNeuIso);
    intree->SetBranchAddress("pfCombIso",      &pfCombIso);       // lepton 4-vector
  
    Double_t mt=-999;
    
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
      // if(fabs(lep->Eta()) < 1.6 || fabs(lep->Eta()) > 2.4) continue;
      // if(fabs(lep->Eta()) > 1.4) continue; 
      // if(fabs(lep->Pt()) <30) continue; 
  
      mt     = sqrt( 2.0 * (lep->Pt()) * (met) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metPhi))) );
 
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
//       corrSigShape *= effdata/effmc;
//       corrBkgShape *= effdata/effmc;
    
      effdata=1; effmc=1;
      effSigShapedata=1;
      effBkgShapedata=1;
      if(q>0) {
        effdata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmSelEff_pos.getEff((lep->Eta()), lep->Pt()); 
//         effSigShapedata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(lep->Eta()), hSelSigSys->GetYaxis()->FindBin(lep->Pt())); 
//         effBkgShapedata *= dataSelEff_pos.getEff((lep->Eta()), lep->Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(lep->Eta()), hSelBkgSys->GetYaxis()->FindBin(lep->Pt())); 
      } else {
        effdata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
        effmc   *= zmmSelEff_neg.getEff((lep->Eta()), lep->Pt()); 
//         effSigShapedata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt())*hSelSigSys->GetBinContent(hSelSigSys->GetXaxis()->FindBin(lep->Eta()), hSelSigSys->GetYaxis()->FindBin(lep->Pt()));
//         effBkgShapedata *= dataSelEff_neg.getEff((lep->Eta()), lep->Pt())*hSelBkgSys->GetBinContent(hSelBkgSys->GetXaxis()->FindBin(lep->Eta()), hSelBkgSys->GetYaxis()->FindBin(lep->Pt())); 
      }
      corr *= effdata/effmc;
//       corrSigShape *= effSigShapedata/effmc;
//       corrBkgShape *= effBkgShapedata/effmc;
      
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
      
	  // if(corr < 0.5) std::cout << "corr  " << corr <<  "  lep eta " << lep->Eta() <<  " let pt" << lep->Pt() << std::endl;
	  // if(corr > 1.0) std::cout << "corr  " << corr <<  "  lep eta " << lep->Eta() <<  " let pt" << lep->Pt() << std::endl;
	  
      effdata=1; effmc=1;
      double var=0.;      
      // STANDALONE
      if(q>0) {
        Double_t effdata = dataStaEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataStaEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataStaEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmStaEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmStaEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmStaEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSta*errSta;
      } else {
        Double_t effdata = dataStaEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataStaEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataStaEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmStaEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmStaEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmStaEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSta = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSta*errSta;
      }
      // SELECTION
      if(q>0) {
        Double_t effdata = dataSelEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataSelEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataSelEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmSelEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmSelEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmSelEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSel*errSel;
      } else {
        Double_t effdata = dataSelEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataSelEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataSelEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmSelEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmSelEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmSelEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errSel = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errSel*errSel;
      }
      //HLT
      if(q>0) {
        Double_t effdata = dataHLTEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataHLTEff_pos.getErrLow(lep->Eta(), lep->Pt()), dataHLTEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmHLTEff_pos.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmHLTEff_pos.getErrLow(lep->Eta(), lep->Pt()), zmmHLTEff_pos.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errHLT*errHLT;
      } else {
        Double_t effdata = dataHLTEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errdata = TMath::Max(dataHLTEff_neg.getErrLow(lep->Eta(), lep->Pt()), dataHLTEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t effmc   = zmmHLTEff_neg.getEff(lep->Eta(), lep->Pt());
        Double_t errmc   = TMath::Max(zmmHLTEff_neg.getErrLow(lep->Eta(), lep->Pt()), zmmHLTEff_neg.getErrHigh(lep->Eta(), lep->Pt()));
        Double_t errHLT = (effdata/effmc)*sqrt(errdata*errdata/effdata/effdata + errmc*errmc/effmc/effmc);
        var+=errHLT*errHLT;
      }

      
      corrUp=corr+sqrt(var);
      corrDown=corr-sqrt(var);  
// std::cout << "made it 0!" << std::endl;
      eff2Bindata=1; eff2Binmc=1;    
          if(q>0) { 
            eff2Bindata *= (1.-dataHLTEff2Bin_pos.getEff((lep->Eta()), lep->Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_pos.getEff((lep->Eta()), lep->Pt())); 
          } else {
            eff2Bindata *= (1.-dataHLTEff2Bin_neg.getEff((lep->Eta()), lep->Pt())); 
            eff2Binmc   *= (1.-zmmHLTEff2Bin_neg.getEff((lep->Eta()), lep->Pt())); 
          }
          eff2Bindata = 1.-eff2Bindata;
          eff2Binmc   = 1.-eff2Binmc;
          corr2Bin *= eff2Bindata/eff2Binmc;
//     std::cout << "made it 1!" << std::endl;
          eff2Bindata=1; eff2Binmc=1;
          if(q>0) { 
            eff2Bindata *= dataSelEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
          } else {
            eff2Bindata *= dataSelEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmSelEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
          }
          corr2Bin *= eff2Bindata/eff2Binmc;
//     std::cout << "made it 2!" << std::endl;
      eff2Bindata=1; eff2Binmc=1;
          if(q>0) { 
            eff2Bindata *= dataStaEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_pos.getEff((lep->Eta()), lep->Pt()); 
          } else {
            eff2Bindata *= dataStaEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
            eff2Binmc   *= zmmStaEff2Bin_neg.getEff((lep->Eta()), lep->Pt()); 
          }
      corr2Bin *= eff2Bindata/eff2Binmc; 
//       std::cout << "made it 3!" << std::endl;
          eff2Bindata=1; eff2Binmc=1;
          
      if(typev[ifile]==eData || typev[ifile]==eAntiData){
        // Apply the Rochester Corrections to data (anti-isolation)
        TLorentzVector mu1;
        mu1.SetPtEtaPhiM(lep->Pt(),lep->Eta(),lep->Phi(),mu_MASS);
        double dtSF1 = rc.kScaleDT(q, mu1.Pt(), mu1.Eta(), mu1.Phi());//, s=0, m=0);
        mu1*=dtSF1;
        
        // std::cout << "pt corr " << mu1.Pt() << "  no corr "  << lep->Pt() << std::endl;
        if(mu1.Pt()        < PT_CUT)  continue;
        // corrected (smear/scale) lepton for MET correction
        TVector2 vLepCor((mu1.Pt())*cos(mu1.Phi()),(mu1.Pt())*sin(mu1.Phi()));
        // calculate the corrected MET
        TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
        Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
        Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
        mt     = sqrt( 2.0 * (mu1.Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(mu1.Phi(),corrMetPhiLepton))) );
        if(doMTCut&&(mt<MT_CUT)) continue;
        if(typev[ifile]==eData) {
          hDataMet->Fill(corrMetWithLepton);
          if(q>0) { 
              if(doMET){hDataMetp->Fill(corrMetWithLepton); }else {hDataMetp->Fill(mt);}
              hDataMetpPhi->Fill(corrMetPhiLepton);
              hMuonEtaDatap->Fill(fabs(mu1.Eta()));
              for(int it=0; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                if(pfCombIso/(lep->Pt()) >= vIsoBins[it] && pfCombIso/(lep->Pt()) < vIsoBins[it+1]) {
                    // if(pfCombIso/(lep->Pt()) <= vIsoBins[0]) std::cout << "pfcomb = 0?" << std::endl;
                    if(doMET){hDataMetpIsoBins[it]->Fill(corrMetWithLepton);}else{hDataMetpIsoBins[it]->Fill(mt);}
					hMetpIsoValues[it]->Fill(pfCombIso/(lep->Pt()));
                    break;
                }
              }
          }
          else    { 
              if(doMET){hDataMetm->Fill(corrMetWithLepton); }else{hDataMetm->Fill(mt); }
              hMuonEtaDatam->Fill(fabs(mu1.Eta()));
              hDataMetmPhi->Fill(corrMetPhiLepton);
              for(int it=0; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                if(pfCombIso/(lep->Pt()) >= vIsoBins[it] && pfCombIso/(lep->Pt()) < vIsoBins[it+1]) {
                    if(doMET){hDataMetmIsoBins[it]->Fill(corrMetWithLepton);}else{hDataMetmIsoBins[it]->Fill(mt);}
					hMetmIsoValues[it]->Fill(pfCombIso/(lep->Pt()));
                    break;
                }
              }
          }
        } else if(typev[ifile]==eAntiData) {
          hAntiDataMet->Fill(corrMetWithLepton);
          if(q>0) { 
              if(doMET){hAntiDataMetp->Fill(corrMetWithLepton);  }else{hAntiDataMetp->Fill(mt); }
              hMuonEtaAntiDatap->Fill(fabs(mu1.Eta()));
              for(int it=1; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                if(pfCombIso/(lep->Pt()) >= vIsoBins[it] && pfCombIso/(lep->Pt()) < vIsoBins[it+1]) {
                    if(doMET){hDataMetpIsoBins[it]->Fill(corrMetWithLepton); }else{hDataMetpIsoBins[it]->Fill(mt);}
					hMetpIsoValues[it]->Fill(pfCombIso/(lep->Pt()));
                    break;
                }
            }
          } 
          else    { 
              if(doMET){hAntiDataMetm->Fill(corrMetWithLepton);  }else{hAntiDataMetm->Fill(mt); }
              hMuonEtaAntiDatap->Fill(fabs(mu1.Eta()));
              for(int it=1; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                if(pfCombIso/(lep->Pt()) >= vIsoBins[it] && pfCombIso/(lep->Pt()) < vIsoBins[it+1]) {
                    if(doMET){hDataMetmIsoBins[it]->Fill(corrMetWithLepton); }else{hDataMetmIsoBins[it]->Fill(mt);}
					hMetmIsoValues[it]->Fill(pfCombIso/(lep->Pt()));
                    break;
                }
            }
          }   
        }
      } else {
        Double_t weight = 1;Double_t weightUp = 1;Double_t weightDown = 1;
        Double_t weight2 =1;
        //corr = 1.0;

        if(pileupUp) {
          weight *= scale1fbUp*lumi*corr;
          weight2 *= scale1fbUp*lumi2*corr;
        } else if(pileupDown){
          weight *= scale1fbDown*lumi*corr;
          weight2 *= scale1fbDown*lumi2*corr;
        } else {
          weight2*=scale1fb*lumi2*corr*prefireWeight*iterator;
          weight *= scale1fb*lumi*corr*prefireWeight*iterator;
        }

        // apply recoil corrections to W MC
        // Double_t lepPt = lep->Pt();
        // lepPt = lepPt + eleCorr*lepPt;
//         Double_t lepPtup = lep->Pt();
//         Double_t lepPtdown = lep->Pt();
		// std::cout << "after fill " << std::endl; 
        
        // Prepare 2-d Vector for raw Muon

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
          if(mu1.Pt()        > PT_CUT) {
            double bin = 0;
            for(int i = 1; i <= hh_diff->GetNbinsX();++i){
              if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            }
            double w2 = 1.0;//hh_diff->GetBinContent(bin);
            
            corrMet=met, corrMetPhi=metPhi;
            hWmunuMet->Fill(corrMet,weight);
//             hWmunuMet_PileupUp->Fill(corrMet,weightUp);
//             hWmunuMet_PileupDown->Fill(corrMet,weightDown);
            if(q>0) {
              hMuonEtaMCp->Fill(fabs(mu1.Eta()),weight);
              // if(doInclusive) recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              if(doInclusive) recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              // else if(doPF) recoilCorrPF->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
// //               else if(pileupUp)    recoilCorrPuUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
// //               else if(pileupDown)  recoilCorrPuDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
              // else if(doKeys){
                  // recoilCorrKeys->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
// //                 if(fabs(genVy)<0.5)
// //                   recoilCorrKeys05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //                 else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
// //                   recoilCorrKeys051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //                 else
// //                   recoilCorrKeys1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys); 
              // } else if(doEta) {
                // if(fabs(genVy)<0.5)
                  // recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                // else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                  // recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
                // else
                  // recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
              // } 
              // Compute the corrected MET value
              TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
              Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
              Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
			   // std::cout << "recCorrMet = " << corrMet <<  "  corrMetLep = " << corrMetWithLepton << "  u1 = "  << pU1 << "  u2 = " << pU2 << "  lep raw pt = " << lep->Pt() << "  lep corr pt = "<< mu1.Pt() << std::endl;
              mt     = sqrt( 2.0 * (mu1.Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(mu1.Phi(),corrMetPhiLepton))) );
              if(doMTCut&&(mt<MT_CUT)) continue;
              if(typev[ifile]==eWmunu)  {
                 if(doMET){hWmunuMetp->Fill(corrMetWithLepton,weight);}else{hWmunuMetp->Fill(mt,weight);}
                 hWmunuMetpPhi->Fill(corrMetPhiLepton);
//                 hWmunuMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
//                 hWmunuMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
                for(int it=0; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                    // if(pfCombIso/(mu1.Pt()) <= vIsoBins[0]) std::cout << "pfcomb = 0?" << std::endl;
                     if(doMET){hWmunuMetpIsoBins[it]->Fill(corrMetWithLepton, weight);}else{hWmunuMetpIsoBins[it]->Fill(mt, weight);}
                     break;
                  }
                }
              }
              else if (typev[ifile]==eWx)
              {
                 hEWKMet->Fill(corrMetWithLepton,weight);
                 if(doMET){hEWKMetp->Fill(corrMetWithLepton,weight);}else{hEWKMetp->Fill(mt,weight);}
                 hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
                 hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
                 for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                     if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                     if(doMET){
                         hWxMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                         hEWKMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                         }else{
                             hWxMetpIsoBins[it]->Fill(mt,weight);
                             hEWKMetpIsoBins[it]->Fill(mt,weight);
                             }
                     // hWxMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                     break;
                  }
               }
             } 
		     else 
		     {
		    	hEWKMet->Fill(corrMetWithLepton,weight);
 				if(doMET){hEWKMetp->Fill(corrMetWithLepton,weight);}else{hEWKMetp->Fill(mt,weight);}
                hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp);
                hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown);
               for(int it=0; it < nIsoBins; ++it){
                   // check the isolation value for this event and fill appropriate histogram
                   if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                       if(doMET){
                           hZxxMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                           hEWKMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                           }else{
                               hZxxMetpIsoBins[it]->Fill(mt,weight);
                               hEWKMetpIsoBins[it]->Fill(mt,weight);
                               }
                       // hZxxMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                  }
                }
		     }

              corrMet=met, corrMetPhi=metPhi;
            } else {
              hMuonEtaMCm->Fill(fabs(mu1.Eta()),weight);
              if(doInclusive) recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
              // else if(doPF) recoilCorrPFm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
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
              mt     = sqrt( 2.0 * (mu1.Pt()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(mu1.Phi(),corrMetPhiLepton))) );
              if(doMTCut&&(mt<MT_CUT)) continue;
              if(typev[ifile]==eWmunu)
              {
                if(doMET){hWmunuMetm->Fill(corrMetWithLepton,weight); }else{hWmunuMetm->Fill(mt,weight);}
                hWmunuMetmPhi->Fill(corrMetPhiLepton);
//                 hWmunuMetm_PileupUp->Fill(corrMetWithLepton,weightUp);
//                 hWmunuMetm_PileupDown->Fill(corrMetWithLepton,weightDown);
                for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                      if(doMET){hWmunuMetmIsoBins[it]->Fill(corrMetWithLepton,weight);}else{hWmunuMetmIsoBins[it]->Fill(mt,weight);}
                      break;
                  }
                }
              }
            else if (typev[ifile]==eWx)
            {
              hEWKMet->Fill(corrMetWithLepton,weight);
              if(doMET){hEWKMetm->Fill(corrMetWithLepton,weight); }else{hEWKMetm->Fill(mt,weight);} corrMet=met, corrMetPhi=metPhi; 
              hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); corrMet=met, corrMetPhi=metPhi; 
              hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); corrMet=met, corrMetPhi=metPhi; 
              for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                      if(doMET){
                          hWxMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                          hEWKMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                          }else{
                              hWxMetmIsoBins[it]->Fill(mt,weight);
                              hEWKMetmIsoBins[it]->Fill(mt,weight);
                              }
                      // hWxMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                 }
              }
            }
		    else
            {
              hEWKMet->Fill(corrMetWithLepton,weight);
              if(doMET){hEWKMetm->Fill(corrMetWithLepton,weight);}else{ hEWKMetm->Fill(mt,weight);} corrMet=met, corrMetPhi=metPhi; 
              hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); corrMet=met, corrMetPhi=metPhi; 
              hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); corrMet=met, corrMetPhi=metPhi; 
              for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                  if(doMET){
                      hZxxMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                      hEWKMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                      }else{
                          hZxxMetmIsoBins[it]->Fill(mt,weight);
                          hEWKMetmIsoBins[it]->Fill(mt,weight);
                          }
                      // hZxxMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                  }
               }
            }
              corrMet=met, corrMetPhi=metPhi;
            }
            corrMet=met, corrMetPhi=metPhi;
          }
          //             recoilCorr.Correct(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,mu1.Phi(),1,q);
          // unused?
          hWmunuMet_RecoilUp->Fill(corrMet,weight);
          if(q>0) {
            pU1 = 0; pU2 = 0; 
            hWmunuMetp_RecoilUp->Fill(corrMet,weight); 
            corrMet=met, corrMetPhi=metPhi;
          } else { 
            pU1 = 0; pU2 = 0; 
            hWmunuMetm_RecoilUp->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          }
          hWmunuMet_RecoilDown->Fill(corrMet,weight);
          if(q>0) {
            pU1 = 0; pU2 = 0; 
            hWmunuMetp_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          } else {
            pU1 = 0; pU2 = 0; 
            hWmunuMetm_RecoilDown->Fill(corrMet,weight);
            corrMet=met, corrMetPhi=metPhi;
          }
        }
        if(typev[ifile]==eAntiWmunu) {
          if(mu1.Pt()        < PT_CUT)  continue;
          Double_t corrMet=met, corrMetPhi=metPhi;
//           TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
//           Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
//        this histogram isn't actually used in any results / fits
          hAntiWmunuMet->Fill(corrMet,weight2);
          if(q>0) {              
            pU1 = 0; pU2 = 0; 
            if(doInclusive) recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
            // else if(doPF) recoilCorrPF->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
// //             else if(pileupUp)    recoilCorrPuUp->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
// //             else if(pileupDown)  recoilCorrPuDown->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
            // else if(doKeys){
                // recoilCorrKeys->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
// //               if(fabs(genVy)<0.5)
// //                 recoilCorrKeys05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //               else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
// //                 recoilCorrKeys051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //               else
// //                 recoilCorrKeys1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys); 
            // } else if(doEta){
              // if(fabs(genVy)<0.5)
                // recoilCorr05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              // else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                // recoilCorr051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              // else
                // recoilCorr1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            // }
            //recoilCorr->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,mu1.Phi(),pU1,pU2,0);
            TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
            Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
            mt     = sqrt( 2.0 * (vLepCor.Mod()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(vLepCor.Phi(),corrMetPhiLepton))) );
            if(doMTCut&&(mt<MT_CUT)) continue;
            if(doMET){hAntiWmunuMetp->Fill(corrMetWithLepton,weight2); }else{hAntiWmunuMetp->Fill(mt,weight2); }
            
            for(int it=1; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                    if(doMET){hWmunuMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);}else{hWmunuMetpIsoBins[it]->Fill(mt, weight2);}
                    break;
                }
            }
            corrMet = met; corrMetPhi = metPhi;
          } 
          else {
            pU1 = 0; pU2 = 0; 
            if(doInclusive) recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
            // else if(doPF) recoilCorrPFm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys,doDiago);
// //             else if(pileupUp)    recoilCorrPuUpm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
// //             else if(pileupDown)  recoilCorrPuDownm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0);
            // else if(doKeys){
                // recoilCorrKeysm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
// //               if(fabs(genVy)<0.5)
// //                 recoilCorrKeysm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //               else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
// //                 recoilCorrKeysm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys);
// //               else
// //                 recoilCorrKeysm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,mu1.Pt(),mu1.Phi(),pU1,pU2,0,0,0,doKeys); 
            // } else if(doEta){
              // if(fabs(genVy)<0.5)
                // recoilCorrm05->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              // else if (fabs(genVy)>=0.5 && fabs(genVy)<1.0)
                // recoilCorrm051->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys);
              // else
                // recoilCorrm1->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lep->Pt(),lep->Phi(),pU1,pU2,0,0,0,doKeys); 
            // }
            //recoilCorrm->CorrectInvCdf(corrMet,corrMetPhi,genVPt,genVPhi,lepPt,mu1.Phi(),pU1,pU2,0);
            TVector2 vMetCorr((corrMet)*cos(corrMetPhi),(corrMet)*sin(corrMetPhi));
            Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
            Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
            mt     = sqrt( 2.0 * (vLepCor.Mod()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(vLepCor.Phi(),corrMetPhiLepton))) );
            if(doMTCut&&(mt<MT_CUT)) continue;
            if(doMET){hAntiWmunuMetm->Fill(corrMetWithLepton,weight2);}else{hAntiWmunuMetm->Fill(mt,weight2);}
            for(int it=1; it < nIsoBins; ++it){
                // check the isolation value for this event and fill appropriate histogram
                if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                    if(doMET){hWmunuMetmIsoBins[it]->Fill(corrMetWithLepton, weight2); }else{hWmunuMetmIsoBins[it]->Fill(mt, weight2);}
                    break;
                }
            }
            corrMet = met; corrMetPhi = metPhi; 
          }
        }
        if(typev[ifile]==eDib || typev[ifile]==eTtb ) {
          if(mu1.Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
          mt     = sqrt( 2.0 * (vLepCor.Mod()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(vLepCor.Phi(),corrMetPhiLepton))) );
          if(doMTCut&&(mt<MT_CUT)) continue;
          if(doMET){hEWKMet->Fill(corrMetWithLepton,weight);}else{hEWKMet->Fill(mt,weight);}
          hEWKMet_PileupUp->Fill(corrMetWithLepton,weightUp);
          hEWKMet_PileupDown->Fill(corrMetWithLepton,weightDown);
          if(q>0) {
            if(doMET){hEWKMetp->Fill(corrMetWithLepton,weight); }else{hEWKMetp->Fill(mt,weight); }
            hEWKMetp_PileupUp->Fill(corrMetWithLepton,weightUp); 
            hEWKMetp_PileupDown->Fill(corrMetWithLepton,weightDown); 
			if(typev[ifile]==eDib){
              for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                      if(doMET){
                          hDibMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                          hEWKMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                          }else{
                              hDibMetpIsoBins[it]->Fill(mt,weight);
                              hEWKMetpIsoBins[it]->Fill(mt,weight);
                              }
                      // hDibMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                  }
              }
			}
			else
			{
              for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                      if(doMET){
                          hTtbMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                          hEWKMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                          }else{
                              hTtbMetpIsoBins[it]->Fill(mt,weight);
                              hEWKMetpIsoBins[it]->Fill(mt,weight);
                              }
                      // hTtbMetpIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                  }
              }
			}
          }
          else { 
            hEWKMetm->Fill(corrMetWithLepton,weight); 
            hEWKMetm_PileupUp->Fill(corrMetWithLepton,weightUp); 
            hEWKMetm_PileupDown->Fill(corrMetWithLepton,weightDown); 
			if(typev[ifile]==eDib){
              for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                      if(doMET){
                          hDibMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                          hEWKMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                          }else{
                              hDibMetmIsoBins[it]->Fill(mt,weight);
                              hEWKMetmIsoBins[it]->Fill(mt,weight);
                              }
                      // hDibMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                  }
              }
			}
			else 
			{
	          for(int it=0; it < nIsoBins; ++it){
                  // check the isolation value for this event and fill appropriate histogram
                  if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                      if(doMET){
                          hTtbMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                          hEWKMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                          }else{
                              hTtbMetmIsoBins[it]->Fill(mt,weight);
                              hEWKMetmIsoBins[it]->Fill(mt,weight);
                              }
                      // hTtbMetmIsoBins[it]->Fill(corrMetWithLepton,weight);
                      break;
                  }
              }
			}
          }
        }
        if(typev[ifile]==eAntiZxx || typev[ifile]==eAntiWx) {
          if(mu1.Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
          mt     = sqrt( 2.0 * (vLepCor.Mod()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(vLepCor.Phi(),corrMetPhiLepton))) );
          if(doMTCut&&(mt<MT_CUT)) continue;
          if(doMET){hAntiEWKMet->Fill(corrMetWithLepton,weight2);}else{hAntiEWKMet->Fill(mt,weight2);}
          if(q>0) { 
              hAntiEWKMetp->Fill(corrMetWithLepton,weight2); 
              if(typev[ifile]==eAntiWx){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hWxMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hWxMetpIsoBins[it]->Fill(mt, weight2);
                                hEWKMetpIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }  
              } if(typev[ifile]==eAntiZxx){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hZxxMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hZxxMetpIsoBins[it]->Fill(mt, weight2);
                                hEWKMetpIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }  
              }
          }
          else { 
              hAntiEWKMetm->Fill(corrMetWithLepton,weight2); 
              if(typev[ifile]==eAntiWx){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hWxMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hWxMetmIsoBins[it]->Fill(mt, weight2);
                                hEWKMetmIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }
              } if(typev[ifile]==eAntiZxx){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hZxxMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hZxxMetmIsoBins[it]->Fill(mt, weight2);
                                hEWKMetmIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }
              }
          }
        }if(typev[ifile]==eAntiDib || typev[ifile]==eAntiTtb) {
          if(mu1.Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          Double_t corrMetPhiLepton = (vMetCorr + vLepRaw - vLepCor).Phi();
          mt     = sqrt( 2.0 * (vLepCor.Mod()) * (corrMetWithLepton) * (1.0-cos(toolbox::deltaPhi(vLepCor.Phi(),corrMetPhiLepton))) );
          if(doMTCut&&(mt<MT_CUT)) continue;
          if(doMET){hAntiEWKMet->Fill(corrMetWithLepton,weight2);}else{hAntiEWKMet->Fill(mt,weight2);}
          if(q>0) { 
              hAntiEWKMetp->Fill(corrMetWithLepton,weight2); 
              if(typev[ifile]==eAntiDib){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hDibMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hDibMetpIsoBins[it]->Fill(mt, weight2);
                                hEWKMetpIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }  
              } if(typev[ifile]==eAntiTtb){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hTtbMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hTtbMetpIsoBins[it]->Fill(mt, weight2);
                                hEWKMetpIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetpIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }  
              }
          }
          else { 
              hAntiEWKMetm->Fill(corrMetWithLepton,weight2); 
              if(typev[ifile]==eAntiDib){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hDibMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hDibMetmIsoBins[it]->Fill(mt, weight2);
                                hEWKMetmIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }
              } if(typev[ifile]==eAntiTtb){
                  for(int it=1; it < nIsoBins; ++it){
                    // check the isolation value for this event and fill appropriate histogram
                    if(pfCombIso/(mu1.Pt()) >= vIsoBins[it] && pfCombIso/(mu1.Pt()) < vIsoBins[it+1]) {
                        if(doMET){
                            hTtbMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            hEWKMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                            }else{
                                hTtbMetmIsoBins[it]->Fill(mt, weight2);
                                hEWKMetmIsoBins[it]->Fill(mt, weight2);
                                }
                        // hTtbMetmIsoBins[it]->Fill(corrMetWithLepton, weight2);
                        break;
                    }
                }
              }
          }
        }
        if(typev[ifile]==eQCD) {
          if(mu1.Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          hQCDMet->Fill(corrMetWithLepton,weight);
          if(q>0) {
            hQCDMetp->Fill(corrMetWithLepton,weight); 
          }
          else { 
            hQCDMetm->Fill(corrMetWithLepton,weight); 
          }
        }
        if(typev[ifile]==eAntiQCD) {
          if(mu1.Pt()        < PT_CUT)  continue;
          TVector2 vMetCorr((met)*cos(metPhi),(met)*sin(metPhi));
          Double_t corrMetWithLepton = (vMetCorr + vLepRaw - vLepCor).Mod();
          hAntiQCDMet->Fill(corrMetWithLepton,weight2);
          if(q>0) { hAntiQCDMetp->Fill(corrMetWithLepton,weight2); }
          else    { hAntiQCDMetm->Fill(corrMetWithLepton,weight2); }
        }
      }
    }
  }
  delete infile;
  infile=0, intree=0;   
 
  // Here rescale the up/down whatever
//   Calculate the shapes for W+
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
   
      	ofstream txtfile2;
    char txtfname2[100];
    std::cout << "Printing We+. " << std::endl;
    sprintf(txtfname2,"%s/isoAvg.txt",CPlot::sOutDir.Data());
    txtfile2.open(txtfname2);
    assert(txtfile2.is_open());
 
     // flags = txtfile.flags();
    // txtfile << setprecision(10);
    // txtfile << " *** Yields *** " << endl;
   char plotname[100];
  TCanvas *c2 = MakeCanvas("c2","c2",800,600);
  for(int j = 0; j < nIsoBins; ++j){

     int binmaxWp = hMetpIsoValues[j]->GetMaximumBin(); 
	 double medianWp = hMetpIsoValues[j]->GetXaxis()->GetBinCenter(binmaxWp);
	 
	 int binmaxWm = hMetmIsoValues[j]->GetMaximumBin(); 
	 double medianWm = hMetmIsoValues[j]->GetXaxis()->GetBinCenter(binmaxWm);
	 
     txtfile2 << "bin+ : " << j << " ; average : " << hMetpIsoValues[j]->GetMean() << " ; median: " << medianWp <<  "  entries: " << hMetpIsoValues[j]->GetEntries() << endl;//<< " low edge " << hMetpIsoValues[j]->GetBinLowEdge(0) << " low edge " << hMetpIsoValues[j]->GetBinLowEdge(1000) << endl;
     txtfile2 << "bin- : " << j << " ; average : " << hMetmIsoValues[j]->GetMean() << " ; median: " << medianWm << "  entries: " << hMetmIsoValues[j]->GetEntries() << endl;
	 sprintf(plotname,"%s/hMetpIso_%d.png",CPlot::sOutDir.Data(),j);
     c2->Clear();//hMetpIsoValues[j]->Draw();
     if(j==0){
		 c2->SetLogy(1);
		 }else{
	     c2->SetLogy(0);
		 }
	 hMetpIsoValues[j]->Draw();
	 c2->Update();c2->SaveAs(plotname);
	 sprintf(plotname,"%s/hMetmIso_%d.png",CPlot::sOutDir.Data(),j);
     c2->Clear();
	  // if(j==0){
		 // c2->SetLogy(kTrue);
		 // }else{
	     // hMetmIsoValues[j]->GetYaxis()->SetLogy(kFalse);
		 // }
	 hMetmIsoValues[j]->Draw();c2->Update();c2->SaveAs(plotname);
  }
  
    txtfile2.close();
    sprintf(plotname,"%s/hMetpPhi_SR.png",CPlot::sOutDir.Data());
    c2->Clear();
    double norm = hDataMetpPhi->Integral()/hWmunuMetpPhi->Integral();
    hWmunuMetpPhi->Scale(norm);
    hWmunuMetpPhi->Draw();
    hDataMetpPhi->Draw("same");
	c2->Update();c2->SaveAs(plotname);
    sprintf(plotname,"%s/hMetmPhi_SR.png",CPlot::sOutDir.Data());
    c2->Clear();
    norm = hDataMetmPhi->Integral()/hWmunuMetmPhi->Integral();
    hWmunuMetmPhi->Scale(norm);
    hWmunuMetmPhi->Draw();
    hDataMetmPhi->Draw("same");
	c2->Update();c2->SaveAs(plotname);
    
  
    ofstream txtfile3;
    char txtfname3[100];
    std::cout << "Printing We+. " << std::endl;
    sprintf(txtfname3,"%s/isoBins_Yields.txt",CPlot::sOutDir.Data());
    txtfile3.open(txtfname3);
    assert(txtfile3.is_open());
	txtfile3 << "W+ " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  nEntries " << hDataMetpIsoBins[i]->Integral() << std::endl;
	}
	txtfile3 << "signal MC " << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  nEntries " << hWmunuMetpIsoBins[i]->Integral() << std::endl;
	}
	txtfile3 << "EWK MC " << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  nEntries " << hEWKMetpIsoBins[i]->Integral() << std::endl;
	}
	
	txtfile3 << "W- " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  nEntries " << hDataMetmIsoBins[i]->Integral() << std::endl;
	}
	txtfile3 << "signal MC " << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  nEntries " << hWmunuMetmIsoBins[i]->Integral() << std::endl;
	}
	txtfile3 << "EWK MC " << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  nEntries " << hEWKMetmIsoBins[i]->Integral() << std::endl;
	}

	txtfile3.close();
  
  
  std::cout << " ==== nEvents === " << std::endl;
  std::cout << "sig single  = " << hWmunuMetm->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWmunuMetmIsoBins[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetm->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetmIsoBins[0]->Integral() << std::endl;
  
  std::cout << "sig single  = " << hWmunuMetp->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWmunuMetpIsoBins[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetp->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetpIsoBins[0]->Integral() << std::endl;
  
  
 
 std::cout << "blah" << std::endl;
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWmunuMet->Integral());
//   cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWmunuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWmunuMet->Integral());
  dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",0.7*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral(),0,hAntiDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nQCDp("nQCDp","nQCDp",hDataMetp->Integral()*0.3,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWmunuMetp->Integral());
//   cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  //RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWmunuMetp->Integral()*1.0,0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWmunuMetp->Integral());
  dewkp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",0.7*(hDataMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm->Integral()*0.3,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWmunuMetm->Integral());
//   cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWmunuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWmunuMetm->Integral());
  dewkm.setConstant(kTRUE);
  RooFormulaVar nAntiEWKm("nAntiEWKm","nAntiEWKm","dewkm*nAntiSigm",RooArgList(nAntiSigm,dewkm));
  
    vector<RooRealVar*> nSigp_(nIsoBins), nQCDp_(nIsoBins), dewkp_(nIsoBins), cewkp_(nIsoBins);
   vector<RooRealVar*> nSigm_(nIsoBins), nQCDm_(nIsoBins), dewkm_(nIsoBins), cewkm_(nIsoBins);
   vector<RooRealVar*> nEWKp_(nIsoBins), nEWKm_(nIsoBins);
   vector<RooRealVar*> nTtbp_(nIsoBins), nTtbm_(nIsoBins);
   vector<RooRealVar*> nDibp_(nIsoBins), nDibm_(nIsoBins);
   vector<RooRealVar*> dwxp_(nIsoBins),  dwxm_(nIsoBins);
   vector<RooFormulaVar*> nWxp_(nIsoBins),  nWxm_(nIsoBins);
   vector<RooRealVar*> nZxxp_(nIsoBins), nZxxm_(nIsoBins);
   std::cout << "blah 2" << std::endl;
   // vector<RooFormulaVar*> nEWKp_(nIsoBins), nEWKm_(nIsoBins);
   // RooRealVar* nEWKpSR = new RooRealVar("nEWKpSR","nEWKpSR",hEWKMetpIsoBins[0]->Integral());
   // RooRealVar* nEWKmSR = new RooRealVar("nEWKmSR","nEWKmSR",hEWKMetmIsoBins[0]->Integral());
  char nname[50];
  char formula[50];
  for (int j = 0; j < nIsoBins; ++j){
	  // W+
	  double qcdFac = 0.9;
      if(j==0) qcdFac = 0.15;
      sprintf(nname, "nAntiSigp%d",j);
      nSigp_[j] = new RooRealVar(nname,nname,hWmunuMetpIsoBins[j]->Integral(),0,hDataMetpIsoBins[j]->Integral());      
      // sprintf(nname, "nAntiQCDp%d",j);
      // nQCDp_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetpIsoBins[j]->Integral(),0,hDataMetpIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiTtbp%d",j);
      nTtbp_[j] = new RooRealVar(nname,nname,hTtbMetpIsoBins[j]->Integral(),0.8*hTtbMetpIsoBins[j]->Integral(),hDataMetpIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiDibp%d",j);
      nDibp_[j] = new RooRealVar(nname,nname,hDibMetpIsoBins[j]->Integral(),0.8*hDibMetpIsoBins[j]->Integral(),hDataMetpIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiZxxp%d",j);
      nZxxp_[j] = new RooRealVar(nname,nname,hZxxMetpIsoBins[j]->Integral(),0.8*hZxxMetpIsoBins[j]->Integral(),hDataMetpIsoBins[j]->Integral());
	  
	  double qcd_remainder = hDataMetpIsoBins[j]->Integral() - hWmunuMetpIsoBins[j]->Integral() - hEWKMetpIsoBins[j]->Integral() ;
	  if(j==0){
         sprintf(nname, "nAntiQCDp%d",j);
		 // nQCDp_[j] = new RooRealVar(nname,nname,0.4*hDataMetpIsoBins[j]->Integral(),0.0*hDataMetpIsoBins[j]->Integral(),0.6*hDataMetpIsoBins[j]->Integral());
         nQCDp_[j] = new RooRealVar(nname,nname,qcd_remainder,0,2.0*qcd_remainder);
         // nQCDp_[j]->setConstant(kTRUE);
	  } else {
         sprintf(nname, "nAntiQCDp%d",j);
         // nQCDp_[j] = new RooRealVar(nname,nname,hDataMetpIsoBins[j]->Integral(),0.95*hDataMetpIsoBins[j]->Integral(),1.05*hDataMetpIsoBins[j]->Integral());
         nQCDp_[j] = new RooRealVar(nname,nname,hDataMetpIsoBins[j]->Integral(),0.95*hDataMetpIsoBins[j]->Integral(),1.05*hDataMetpIsoBins[j]->Integral());
	     // nQCDp_[j]->setConstant(kTRUE);
      }
	  // nQCDp_[j]->setConstant(kTRUE);
	  
	  	  // the ratio between W->enu and W->tau nu
	  sprintf(nname, "dwxp%d",j);
      dwxp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dwxp_[j]->setVal(hWxMetpIsoBins[j]->Integral()/hWmunuMetpIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiWxp%d",j); sprintf(formula,"dwxp%d*nAntiSigp%d",j,j);
	  nWxp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*dwxp_[j],*nSigp_[j]));
	  
      // sprintf(nname, "dewkp%d",j);
      // dewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkp_[j]->setVal(hEWKMetpIsoBins[j]->Integral()/hWmunuMetpIsoBins[j]->Integral());
	  
      sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"dewkp%d*nAntiSigp%d",j,j);
      // // nEWKp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigp_[j],*dewkp_[j]));
	  nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetpIsoBins[j]->Integral(),0.0,(2.0)*hEWKMetpIsoBins[j]->Integral());
	
	  // sprintf(nname, "nAntiDibp%d",j); 
	  // nDibp_[j] = new RooRealVar(nname,nname,hDibMetpIsoBins[j]->Integral(),0.0,(2.0)*hDibMetpIsoBins[j]->Integral());	  
	  
	  // sprintf(nname, "nAntiTtbp%d",j); 
	  // nTtbp_[j] = new RooRealVar(nname,nname,hTtbMetpIsoBins[j]->Integral(),0.0,(2.0)*hTtbMetpIsoBins[j]->Integral());
	
	 
      // W-	 
	  sprintf(nname, "nAntiSigm%d",j);
      nSigm_[j] = new RooRealVar(nname,nname,hWmunuMetmIsoBins[j]->Integral(),0,hDataMetmIsoBins[j]->Integral());
	  
	  	  
	  sprintf(nname, "nAntiTtbm%d",j);
      nTtbm_[j] = new RooRealVar(nname,nname,hTtbMetmIsoBins[j]->Integral(),0.8*hTtbMetmIsoBins[j]->Integral(),hDataMetmIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiDibm%d",j);
      nDibm_[j] = new RooRealVar(nname,nname,hDibMetmIsoBins[j]->Integral(),0.8*hDibMetmIsoBins[j]->Integral(),hDataMetmIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiZxxm%d",j);
      nZxxm_[j] = new RooRealVar(nname,nname,hZxxMetmIsoBins[j]->Integral(),0.8*hZxxMetmIsoBins[j]->Integral(),hDataMetmIsoBins[j]->Integral());

      // sprintf(nname, "nAntiQCDm%d",j);
      // nQCDm_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetpIsoBins[j]->Integral(),0,hDataMetmIsoBins[j]->Integral());
	  
	  qcd_remainder = hDataMetmIsoBins[j]->Integral() - hWmunuMetmIsoBins[j]->Integral() - hEWKMetmIsoBins[j]->Integral() ;
	  if(j==0){
         sprintf(nname, "nAntiQCDm%d",j);
         // nQCDm_[j] = new RooRealVar(nname,nname,0.4*hDataMetmIsoBins[j]->Integral(),0.0*hDataMetmIsoBins[j]->Integral(),0.6*hDataMetmIsoBins[j]->Integral());
         nQCDm_[j] = new RooRealVar(nname,nname,qcd_remainder,0,qcd_remainder*2.0);
         // nQCDm_[j]->setConstant(kTRUE);
	  } else {
		 sprintf(nname, "nAntiQCDm%d",j);
         // nQCDm_[j] = new RooRealVar(nname,nname,hDataMetmIsoBins[j]->Integral(),0.95*hDataMetmIsoBins[j]->Integral(),1.05*hDataMetmIsoBins[j]->Integral());
         nQCDm_[j] = new RooRealVar(nname,nname,hDataMetmIsoBins[j]->Integral(),0.95*hDataMetmIsoBins[j]->Integral(),1.05*hDataMetmIsoBins[j]->Integral());
	     // nQCDm_[j]->setConstant(kTRUE);
      }
	  
	  // nQCDm_[j]->setConstant(kTRUE);
	  
	  	  // the W->tau nu / W->e nu ratio
	  sprintf(nname, "dwxm%d",j);
      dwxm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dwxm_[j]->setVal(hWxMetmIsoBins[j]->Integral()/hWmunuMetmIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiWxm%d",j); sprintf(formula,"dwxm%d*nAntiSigm%d",j,j);
	  nWxm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*dwxm_[j],*nSigm_[j]));
	  
      // sprintf(nname, "dewkm%d",j);
      // dewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkm_[j]->setVal(hEWKMetmIsoBins[j]->Integral()/hWmunuMetmIsoBins[j]->Integral());
	  
      sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      // // nEWKm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigm_[j],*dewkm_[j]));
	  nEWKm_[j] = new RooRealVar(nname,nname,hEWKMetmIsoBins[j]->Integral(),0.0,(2.0)*hEWKMetmIsoBins[j]->Integral());
	  
	  // sprintf(nname, "nAntiDibm%d",j); 
	  // nDibm_[j] = new RooRealVar(nname,nname,hDibMetmIsoBins[j]->Integral(),0.0,(2.0)*hDibMetmIsoBins[j]->Integral());
	  
	  // sprintf(nname, "nAntiTtbm%d",j); 
	  // nTtbm_[j] = new RooRealVar(nname,nname,hTtbMetmIsoBins[j]->Integral(),0.0,(2.0)*hTtbMetmIsoBins[j]->Integral());
	  
	  /*
      double qcdFac = 0.9;
      if(j==0) qcdFac = 0.15;
	  double QCDlims = 0.5;
	  if(j==0) QCDlims = 0.0;
      sprintf(nname, "nAntiSigp%d",j);
      nSigp_[j] = new RooRealVar(nname,nname,hWmunuMetpIsoBins[j]->Integral(),0.5*hDataMetpIsoBins[j]->Integral(),hDataMetpIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiQCDp%d",j);
      nQCDp_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetpIsoBins[j]->Integral(),(0+QCDlims)*hDataMetpIsoBins[j]->Integral(),(1+QCDlims)*hDataMetpIsoBins[j]->Integral());
	  
      sprintf(nname, "dewkp%d",j);
      dewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dewkp_[j]->setVal(hEWKMetpIsoBins[j]->Integral()/hWmunuMetpIsoBins[j]->Integral());
	  
      // sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"dewkp%d*nAntiSigp%d",j,j);
      // nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetpIsoBins[j]->Integral(),(0.5)*hEWKMetpIsoBins[j]->Integral(),(1.5)*hEWKMetpIsoBins[j]->Integral());
	  sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"dewkp%d*nAntiSigp%d",j,j);
      nEWKp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigp_[j],*dewkp_[j]));
	
	 
      // W-	 
	  sprintf(nname, "nAntiSigm%d",j);
      nSigm_[j] = new RooRealVar(nname,nname,hWmunuMetmIsoBins[j]->Integral(),0,hDataMetmIsoBins[j]->Integral());

      sprintf(nname, "nAntiQCDm%d",j);
      nQCDm_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetmIsoBins[j]->Integral(),(0+QCDlims)*hDataMetmIsoBins[j]->Integral(),(1+QCDlims)*hDataMetmIsoBins[j]->Integral());
	  
      sprintf(nname, "dewkm%d",j);
      dewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dewkm_[j]->setVal(hEWKMetmIsoBins[j]->Integral()/hWmunuMetmIsoBins[j]->Integral());
	  
      // sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      // nEWKm_[j] = new RooRealVar(nname,nname,hEWKMetmIsoBins[j]->Integral(),(0.5)*hEWKMetmIsoBins[j]->Integral(),(1.5)*hEWKMetmIsoBins[j]->Integral());
	  
	  sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      nEWKm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigm_[j],*dewkm_[j]));
	  */
	 
	  // if(j == 0){
         // sprintf(nname, "nAntiQCDp%d",j);
         // nQCDp_[j] = new RooRealVar(nname,nname,0.15*hDataMetpIsoBins[j]->Integral(),0.01*hDataMetpIsoBins[j]->Integral(),0.4*hDataMetpIsoBins[j]->Integral());
		 // sprintf(nname, "nAntiSigp%d",j);
	     // nSigp_[j] = new RooRealVar(nname,nname,hWmunuMetpIsoBins[j]->Integral(),0.5*hDataMetpIsoBins[j]->Integral(),hDataMetpIsoBins[j]->Integral());
		 // // sprintf(nname, "cewkp%d",j);
         // // cewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
         // // cewkp_[j]->setVal(1);
		 // // sprintf(nname, "nAntiEWKp%d",j); //sprintf(formula,"cewkp%d*nAntiSigp%d",j);
         // // nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetpIsoBins[j]->Integral(),0.01*hDataMetpIsoBins[j]->Integral(),0.5*hDataMetpIsoBins[j]->Integral());
		 
	     // // nEWKp_[j] = new RooFormulaVar(nname, nname, "nEWKpSR", RooArgList(*nEWKpSR));
	  // } else {
         // sprintf(nname, "nAntiQCDp%d",j);
         // nQCDp_[j] = new RooRealVar(nname,nname,hDataMetpIsoBins[j]->Integral(),0.5*hDataMetpIsoBins[j]->Integral(),1.1*hDataMetpIsoBins[j]->Integral());
		 // sprintf(nname, "nAntiSigp%d",j);
	     // nSigp_[j] = new RooRealVar(nname,nname,hWmunuMetpIsoBins[j]->Integral(),0,0.5*hDataMetpIsoBins[j]->Integral());
		 // sprintf(nname, "cewkp%d",j);
         // // cewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
         // // cewkp_[j]->setVal(hEWKMetpIsoBins[j]->Integral()/hEWKMetpIsoBins[0]->Integral());
		 // // cewkp_[j]->setConstant(kTRUE);
		 // // sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"cewkp%d*nEWKpSR",j);
		 // // nEWKp_[j] = new RooFormulaVar(nname, nname, formula, RooArgList(*cewkp_[j],*nEWKpSR));
	  // }
	  
      // sprintf(nname, "dewkp%d",j);
      // dewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkp_[j]->setVal(hEWKMetpIsoBins[j]->Integral()/hWmunuMetpIsoBins[j]->Integral());
	  
	  // // sprintf
	  
      // // sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"cewkp%d*nAntiSigp%d",j,j);
      // // nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetpIsoBins[j]->Integral(),0.01*hDataMetpIsoBins[j]->Integral(),0.5*hDataMetpIsoBins[j]->Integral());
	  // // nEWKp_[j] = new RooFormulaVar(nname, nname, formula, RooArgList(cewkp_[j],));
	  
	  // sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"dewkp%d*nAntiSigp%d",j,j);
      // nEWKp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigp_[j],*dewkp_[j]));
	  
	  
	
	 
      // W-	 
	  // sprintf(nname, "nAntiSigm%d",j);
      // nSigm_[j] = new RooRealVar(nname,nname,hWmunuMetmIsoBins[j]->Integral(),0.5*hDataMetmIsoBins[j]->Integral(),hDataMetmIsoBins[j]->Integral());

	  // if(j==0){
         // sprintf(nname, "nAntiQCDm%d",j);
         // nQCDm_[j] = new RooRealVar(nname,nname,0.15*hDataMetmIsoBins[j]->Integral(),0.01*hDataMetmIsoBins[j]->Integral(),0.4*hDataMetmIsoBins[j]->Integral());
		 // sprintf(nname, "nAntiSigm%d",j);
	     // nSigm_[j] = new RooRealVar(nname,nname,hWmunuMetmIsoBins[j]->Integral(),0.5*hDataMetmIsoBins[j]->Integral(),hDataMetmIsoBins[j]->Integral());
		 // // sprintf(nname, "cewkm%d",j);
         // // cewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
         // // cewkm_[j]->setVal(1);
		 // // sprintf(nname, "nAntiEWKm%d",j); //sprintf(formula,"cewkp%d*nAntiSigp%d",j);
         // // nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetpIsoBins[j]->Integral(),0.01*hDataMetpIsoBins[j]->Integral(),0.5*hDataMetpIsoBins[j]->Integral());
		 
	     // // nEWKm_[j] = new RooFormulaVar(nname, nname, "nEWKmSR", RooArgList(*nEWKmSR));
	  // } else {
		 // sprintf(nname, "nAntiQCDm%d",j);
         // nQCDm_[j] = new RooRealVar(nname,nname,hDataMetmIsoBins[j]->Integral(),0.5*hDataMetmIsoBins[j]->Integral(),1.1*hDataMetmIsoBins[j]->Integral());
		 // sprintf(nname, "nAntiSigm%d",j);
	     // nSigm_[j] = new RooRealVar(nname,nname,hWmunuMetmIsoBins[j]->Integral(),0,0.5*hDataMetmIsoBins[j]->Integral());
		 // sprintf(nname, "cewkm%d",j);
         // // cewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
         // // cewkm_[j]->setVal(hEWKMetmIsoBins[j]->Integral()/hEWKMetmIsoBins[0]->Integral());
		 // // cewkm_[j]->setConstant(kTRUE);
		 // // sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"cewkm%d*nEWKmSR",j);
		 // // nEWKm_[j] = new RooFormulaVar(nname, nname, formula, RooArgList(*cewkm_[j],*nEWKmSR));
	  // }
	  
      // sprintf(nname, "dewkm%d",j);
      // dewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkm_[j]->setVal(hEWKMetmIsoBins[j]->Integral()/hWmunuMetmIsoBins[j]->Integral());
	  
	  // sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      // nEWKm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigm_[j],*dewkm_[j]));
	  
      // sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      // nEWKm_[j] = new RooRealVar(nname,nname,hEWKMetmIsoBins[j]->Integral(),0.01*hDataMetmIsoBins[j]->Integral(),0.5*hDataMetmIsoBins[j]->Integral());
	  
	  // for now we are setting the control region EWK and Signal yields to be zero
	  // if(j >= 1) {
		  // nSigp_[j]->setVal(0); nSigp_[j]->setConstant(kTRUE);
		  // nSigm_[j]->setVal(0); nSigm_[j]->setConstant(kTRUE);
		  // nEWKp_[j]->setVal(0.0); nEWKp_[j]->setConstant(kTRUE);
		  // nEWKm_[j]->setVal(0.0); nEWKm_[j]->setConstant(kTRUE);
		  // nDibp_[j]->setVal(0.0); nDibp_[j]->setConstant(kTRUE);
		  // nDibm_[j]->setVal(0.0); nDibm_[j]->setConstant(kTRUE);
		  // nTtbp_[j]->setVal(0.0); nTtbp_[j]->setConstant(kTRUE);
		  // nTtbm_[j]->setVal(0.0); nTtbm_[j]->setConstant(kTRUE);
	  // }
	  
	  std::cout << "blah 3." << j << std::endl;
	  
  }
  
   std::cout << "finished normalizations" << std::endl;
  
  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,wmunuMetm,1); 
  RooDataHist wmunuMet_RecoilUp("wmunuMET_RecoilUp", "wmunuMET_RecoilUp", RooArgSet(pfmet),hWmunuMet_RecoilUp);  RooHistPdf pdfWm_RecoilUp("wm_RecoilUp", "wm_RecoilUp", pfmet,wmunuMet_RecoilUp, 1);
  RooDataHist wmunuMetp_RecoilUp("wmunuMETp_RecoilUp","wmunuMETp_RecoilUp",RooArgSet(pfmet),hWmunuMetp_RecoilUp); RooHistPdf pdfWmp_RecoilUp("wmp_RecoilUp","wmp_RecoilUp",pfmet,wmunuMetp_RecoilUp,1);
  RooDataHist wmunuMetm_RecoilUp("wmunuMETm_RecoilUp","wmunuMETm_RecoilUp",RooArgSet(pfmet),hWmunuMetm_RecoilUp); RooHistPdf pdfWmm_RecoilUp("wmm_RecoilUp","wmm_RecoilUp",pfmet,wmunuMetm_RecoilUp,1); 
  RooDataHist wmunuMet_RecoilDown("wmunuMET_RecoilDown", "wmunuMET_RecoilDown", RooArgSet(pfmet),hWmunuMet_RecoilDown);  RooHistPdf pdfWm_RecoilDown("wm_RecoilDown", "wm_RecoilDown", pfmet,wmunuMet_RecoilDown, 1);
  RooDataHist wmunuMetp_RecoilDown("wmunuMETp_RecoilDown","wmunuMETp_RecoilDown",RooArgSet(pfmet),hWmunuMetp_RecoilDown); RooHistPdf pdfWmp_RecoilDown("wmp_RecoilDown","wmp_RecoilDown",pfmet,wmunuMetp_RecoilDown,1);
  RooDataHist wmunuMetm_RecoilDown("wmunuMETm_RecoilDown","wmunuMETm_RecoilDown",RooArgSet(pfmet),hWmunuMetm_RecoilDown); RooHistPdf pdfWmm_RecoilDown("wmm_RecoilDown","wmm_RecoilDown",pfmet,wmunuMetm_RecoilDown,1); 
   RooDataHist wmunuMet_ScaleUp("wmunuMET_ScaleUp", "wmunuMET_ScaleUp", RooArgSet(pfmet),hWmunuMet_ScaleUp);  RooHistPdf pdfWm_ScaleUp("wm_ScaleUp", "wm_ScaleUp", pfmet,wmunuMet_ScaleUp, 1);
  RooDataHist wmunuMetp_ScaleUp("wmunuMETp_ScaleUp","wmunuMETp_ScaleUp",RooArgSet(pfmet),hWmunuMetp_ScaleUp); RooHistPdf pdfWmp_ScaleUp("wmp_ScaleUp","wmp_ScaleUp",pfmet,wmunuMetp_ScaleUp,1);
  RooDataHist wmunuMetm_ScaleUp("wmunuMETm_ScaleUp","wmunuMETm_ScaleUp",RooArgSet(pfmet),hWmunuMetm_ScaleUp); RooHistPdf pdfWmm_ScaleUp("wmm_ScaleUp","wmm_ScaleUp",pfmet,wmunuMetm_ScaleUp,1); 
  RooDataHist wmunuMet_ScaleDown("wmunuMET_ScaleDown", "wmunuMET_ScaleDown", RooArgSet(pfmet),hWmunuMet_ScaleDown);  RooHistPdf pdfWm_ScaleDown("wm_ScaleDown", "wm_ScaleDown", pfmet,wmunuMet_ScaleDown, 1);
  RooDataHist wmunuMetp_ScaleDown("wmunuMETp_ScaleDown","wmunuMETp_ScaleDown",RooArgSet(pfmet),hWmunuMetp_ScaleDown); RooHistPdf pdfWmp_ScaleDown("wmp_ScaleDown","wmp_ScaleDown",pfmet,wmunuMetp_ScaleDown,1);
  RooDataHist wmunuMetm_ScaleDown("wmunuMETm_ScaleDown","wmunuMETm_ScaleDown",RooArgSet(pfmet),hWmunuMetm_ScaleDown); RooHistPdf pdfWmm_ScaleDown("wmm_ScaleDown","wmm_ScaleDown",pfmet,wmunuMetm_ScaleDown,1); 
  
    RooDataHist wmunuMet_PileupUp("wmunuMET_PileupUp", "wmunuMET_PileupUp", RooArgSet(pfmet),hWmunuMet_PileupUp);  RooHistPdf pdfWm_PileupUp("wm_PileupUp", "wm_PileupUp", pfmet,wmunuMet_PileupUp, 1);
  RooDataHist wmunuMetp_PileupUp("wmunuMETp_PileupUp","wmunuMETp_PileupUp",RooArgSet(pfmet),hWmunuMetp_PileupUp); RooHistPdf pdfWmp_PileupUp("wmp_PileupUp","wmp_PileupUp",pfmet,wmunuMetp_PileupUp,1);
  RooDataHist wmunuMetm_PileupUp("wmunuMETm_PileupUp","wmunuMETm_PileupUp",RooArgSet(pfmet),hWmunuMetm_PileupUp); RooHistPdf pdfWmm_PileupUp("wmm_PileupUp","wmm_PileupUp",pfmet,wmunuMetm_PileupUp,1); 
  RooDataHist wmunuMet_PileupDown("wmunuMET_PileupDown", "wmunuMET_PileupDown", RooArgSet(pfmet),hWmunuMet_PileupDown);  RooHistPdf pdfWm_PileupDown("wm_PileupDown", "wm_PileupDown", pfmet,wmunuMet_PileupDown, 1);
  RooDataHist wmunuMetp_PileupDown("wmunuMETp_PileupDown","wmunuMETp_PileupDown",RooArgSet(pfmet),hWmunuMetp_PileupDown); RooHistPdf pdfWmp_PileupDown("wmp_PileupDown","wmp_PileupDown",pfmet,wmunuMetp_PileupDown,1);
  RooDataHist wmunuMetm_PileupDown("wmunuMETm_PileupDown","wmunuMETm_PileupDown",RooArgSet(pfmet),hWmunuMetm_PileupDown); RooHistPdf pdfWmm_PileupDown("wmm_PileupDown","wmm_PileupDown",pfmet,wmunuMetm_PileupDown,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
  RooDataHist ewkMet_PileupUp("ewkMET_PileupUp", "ewkMET_PileupUp", RooArgSet(pfmet),hEWKMet_PileupUp);  RooHistPdf pdfEWK_PileupUp("ewk_PileupUp", "ewk_PileupUp", pfmet,ewkMet_PileupUp, 1);
  RooDataHist ewkMetp_PileupUp("ewkMETp_PileupUp","ewkMETp_PileupUp",RooArgSet(pfmet),hEWKMetp_PileupUp); RooHistPdf pdfEWKp_PileupUp("ewkp_PileupUp","ewkp_PileupUp",pfmet,ewkMetp_PileupUp,1);
  RooDataHist ewkMetm_PileupUp("ewkMETm_PileupUp","ewkMETm_PileupUp",RooArgSet(pfmet),hEWKMetm_PileupUp); RooHistPdf pdfEWKm_PileupUp("ewkm_PileupUp","ewkm_PileupUp",pfmet,ewkMetm_PileupUp,1); 
  RooDataHist ewkMet_PileupDown("ewkMET_PileupDown", "ewkMET_PileupDown", RooArgSet(pfmet),hEWKMet_PileupDown);  RooHistPdf pdfEWK_PileupDown("ewk_PileupDown", "ewk_PileupDown", pfmet,ewkMet_PileupDown, 1);
  RooDataHist ewkMetp_PileupDown("ewkMETp_PileupDown","ewkMETp_PileupDown",RooArgSet(pfmet),hEWKMetp_PileupDown); RooHistPdf pdfEWKp_PileupDown("ewkp_PileupDown","ewkp_PileupDown",pfmet,ewkMetp_PileupDown,1);
  RooDataHist ewkMetm_PileupDown("ewkMETm_PileupDown","ewkMETm_PileupDown",RooArgSet(pfmet),hEWKMetm_PileupDown); RooHistPdf pdfEWKm_PileupDown("ewkm_PileupDown","ewkm_PileupDown",pfmet,ewkMetm_PileupDown,1); 
  
//   // test using the reversed dEta, dPhi cuts as background
  RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
  RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
  RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  
  vector<RooDataHist*> wmunuMetp_(nIsoBins), wmunuMetm_(nIsoBins);
  vector<RooDataHist*> ewkMetp_(nIsoBins), ewkMetm_(nIsoBins);
  vector<RooDataHist*> qcdMetp_(nIsoBins), qcdMetm_(nIsoBins);
  vector<RooDataHist*> ttbMetp_(nIsoBins), dibMetp_(nIsoBins), wxMetp_(nIsoBins), zxxMetp_(nIsoBins);
  vector<RooDataHist*> ttbMetm_(nIsoBins), dibMetm_(nIsoBins), wxMetm_(nIsoBins), zxxMetm_(nIsoBins);
  
  vector<RooHistPdf*> pdfWep_(nIsoBins), pdfEWKp_(nIsoBins);
  vector<RooHistPdf*> pdfQCDp_(nIsoBins), pdfQCDm_(nIsoBins);
  vector<RooHistPdf*> pdfWem_(nIsoBins), pdfEWKm_(nIsoBins);
  vector<RooHistPdf*> pdfTtbp_(nIsoBins), pdfDibp_(nIsoBins), pdfWxp_(nIsoBins), pdfZxxp_(nIsoBins);
  vector<RooHistPdf*> pdfTtbm_(nIsoBins), pdfDibm_(nIsoBins), pdfWxm_(nIsoBins), pdfZxxm_(nIsoBins);
  
  
  TH1D *hDatapFirstBinSubtracted = (TH1D*) hDataMetpIsoBins[1]->Clone("hDatapFirstBinSubtracted");
  TH1D *hDatamFirstBinSubtracted = (TH1D*) hDataMetmIsoBins[1]->Clone("hDatamFirstBinSubtracted");
  hDatapFirstBinSubtracted->Add(hWmunuMetpIsoBins[1],-1);
  hDatapFirstBinSubtracted->Add(hEWKMetpIsoBins[1],-1);
  hDatamFirstBinSubtracted->Add(hWmunuMetmIsoBins[1],-1);
  hDatamFirstBinSubtracted->Add(hEWKMetmIsoBins[1],-1);
  
   std::cout << "making pdfs from histograms" << std::endl;
  // char hname[50];
  // char pname[50];
  for (int j = 0; j < nIsoBins; ++j){
      // signal pdfs
      sprintf(nname, "wmunuMETp%d",j); wmunuMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWmunuMetpIsoBins[j]);
      sprintf(nname, "wep%d",j); pdfWep_[j] = new RooHistPdf(nname,nname,pfmet,*wmunuMetp_[j],1);      
      // ewk pdfs
      sprintf(nname, "ewkMETp%d",j); ewkMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hEWKMetpIsoBins[j]);
      sprintf(nname, "ewkp%d",j); pdfEWKp_[j] = new RooHistPdf(nname,nname,pfmet,*ewkMetp_[j],1);
      
      if(j==0){
          sprintf(nname, "qcdMetp%d",j); qcdMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDatapFirstBinSubtracted);
      } else{
          sprintf(nname, "qcdMetp%d",j); qcdMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDataMetpIsoBins[j]);
      }
      sprintf(nname, "qcdp%d",j); pdfQCDp_[j] = new RooHistPdf(nname,nname,pfmet,*qcdMetp_[j],1);
      
	  // Split EWK into W, Z, ttbar, and di-boson
	  sprintf(nname, "ttbMETp%d",j); ttbMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hTtbMetpIsoBins[j]);
      sprintf(nname, "ttbp%d",j); pdfTtbp_[j] = new RooHistPdf(nname,nname,pfmet,*ttbMetp_[j],1);
	  
	  sprintf(nname, "dibMETp%d",j); dibMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDibMetpIsoBins[j]);
      sprintf(nname, "dibp%d",j); pdfDibp_[j] = new RooHistPdf(nname,nname,pfmet,*dibMetp_[j],1);
	  
	  sprintf(nname, "wxMETp%d",j); wxMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWxMetpIsoBins[j]);
      sprintf(nname, "wxp%d",j); pdfWxp_[j] = new RooHistPdf(nname,nname,pfmet,*wxMetp_[j],1);
	  
	  sprintf(nname, "zxxMETp%d",j); zxxMetp_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hZxxMetpIsoBins[j]);
      sprintf(nname, "zxxp%d",j); pdfZxxp_[j] = new RooHistPdf(nname,nname,pfmet,*zxxMetp_[j],1);
	  
	  // ----------------------------------------- W- ---------------------------------
	  
	  sprintf(nname, "wmunuMETm%d",j); wmunuMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWmunuMetmIsoBins[j]);
      sprintf(nname, "wem%d",j); pdfWem_[j] = new RooHistPdf(nname,nname,pfmet,*wmunuMetm_[j],1);      
      // ewk pdfs
      sprintf(nname, "ewkMETm%d",j); ewkMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hEWKMetmIsoBins[j]);
      sprintf(nname, "ewkm%d",j); pdfEWKm_[j] = new RooHistPdf(nname,nname,pfmet,*ewkMetm_[j],1);
	  
	  if(j==0){
          sprintf(nname, "qcdMetm%d",j); qcdMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDatamFirstBinSubtracted);
      } else{
          sprintf(nname, "qcdMetm%d",j); qcdMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDataMetmIsoBins[j]);
      }
      sprintf(nname, "qcdm%d",j); pdfQCDm_[j] = new RooHistPdf(nname,nname,pfmet,*qcdMetm_[j],1);
      
	  // Split the EWK into W, Z, ttbar, and di-boson
	  sprintf(nname, "ttbMETm%d",j); ttbMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hTtbMetmIsoBins[j]);
      sprintf(nname, "ttbm%d",j); pdfTtbm_[j] = new RooHistPdf(nname,nname,pfmet,*ttbMetm_[j],1);
	  
	  sprintf(nname, "dibMETm%d",j); dibMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hDibMetmIsoBins[j]);
      sprintf(nname, "dibm%d",j); pdfDibm_[j] = new RooHistPdf(nname,nname,pfmet,*dibMetm_[j],1);
	  
	  sprintf(nname, "wxMETm%d",j); wxMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hWxMetmIsoBins[j]);
      sprintf(nname, "wxm%d",j); pdfWxm_[j] = new RooHistPdf(nname,nname,pfmet,*wxMetm_[j],1);
	  
	  sprintf(nname, "zxxMETm%d",j); zxxMetm_[j] = new RooDataHist(nname, nname, RooArgSet(pfmet),hZxxMetmIsoBins[j]);
      sprintf(nname, "zxxm%d",j); pdfZxxm_[j] = new RooHistPdf(nname,nname,pfmet,*zxxMetm_[j],1);
  }
  
  
  // QCD Pdfs
  
  //CExponential qcd(pfmet,kTRUE);
  //CExponential qcdp(pfmet,kTRUE);
  //CExponential qcdm(pfmet,kTRUE);
  // comment back in for qcd functional form
  CPepeModel2 qcd("qcd",pfmet);
  CPepeModel2 qcdp("qcdp",pfmet);
  CPepeModel2 qcdm("qcdm",pfmet);
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
    
    RooGaussian constantim("constantim","constantim",nAntiSigm,RooConst(hAntiWmunuMetm->Integral()),RooConst(0.15*hAntiWmunuMetm->Integral()));
    RooGaussian constantip("constantip","constantip",nAntiSigp,RooConst(hAntiWmunuMetp->Integral()),RooConst(0.15*hAntiWmunuMetp->Integral()));
	

 
//   // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));

//   // Signal + Background PDFs
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,pdfQCD),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,pdfQCDp),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,pdfQCDm),RooArgList(nSigm,nEWKm,nQCDm));
  
//   // constrained background?
//   RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
//   RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,qcdpc),RooArgList(nSigp,nEWKp,nQCDp));
//   RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,qcdmc),RooArgList(nSigm,nEWKm,nQCDm));
  
  
    
  
  // Anti-Signal PDFs
  RooDataHist awmunuMet ("awmunuMET", "awmunuMET", RooArgSet(pfmet),hAntiWmunuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,awmunuMet, 1);
  RooDataHist awmunuMetp("awmunuMETp","awmunuMETp",RooArgSet(pfmet),hAntiWmunuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,awmunuMetp,1);
  RooDataHist awmunuMetm("awmunuMETm","awmunuMETm",RooArgSet(pfmet),hAntiWmunuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,awmunuMetm,1); 
  
  // Anti-EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  
//   // Anti-QCD Pdfs
  CPepeModel2 aqcd("aqcd",pfmet, qcd.a1);
  CPepeModel2 aqcdp("aqcdp",pfmet, qcdp.a1);
  CPepeModel2 aqcdm("aqcdm",pfmet, qcdm.a1);

//   CPepeModel2 aqcd("aqcd",pfmet);
//   CPepeModel2 aqcdp("aqcdp",pfmet);
//   CPepeModel2 aqcdm("aqcdm",pfmet);
  
  // Anti-selection PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
  
  // regular linear dependence
  // vector <CPepeModel2isobinsMuons*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  // qcdp_[0] = new CPepeModel2isobinsMuons("qcdp2d0",pfmet, 0.075);
  // qcdm_[0] = new CPepeModel2isobinsMuons("qcdm2d0",pfmet, 0.075);
  
  // quadratic dependence
  // vector <CPepeModel2isobinsQuad*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  // qcdp_[0] = new CPepeModel2isobinsQuad("qcdp2d0",pfmet, 0.075);
  // qcdm_[0] = new CPepeModel2isobinsQuad("qcdm2d0",pfmet, 0.075);
  
  // totally uncorrelated
  vector <CPepeModel2*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  qcdp_[0] = new CPepeModel2("qcdp2d0",pfmet);
  qcdm_[0] = new CPepeModel2("qcdm2d0",pfmet);
  
  for (int j = 1; j < nIsoBins; ++j){
	  
      // quadratic isolation dependence
	  // sprintf(nname, "qcdp2d%d",j);
	  // qcdp_[j] = new CPepeModel2isobinsQuad(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->b1, qcdp_[0]->b2, qcdp_[0]->b3, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  // printf(nname, "qcdm2d%d",j);
	  // qcdm_[j] = new CPepeModel2isobinsQuad(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2,  qcdm_[0]->b1, qcdm_[0]->b2, qcdm_[0]->b3, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
	  
	  // // Original isolation linear dependence
      // sprintf(nname, "qcdp2d%d",j);
      // qcdp_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  // sprintf(nname, "qcdm2d%d",j);
      // qcdm_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
  
      sprintf(nname, "qcdp2d%d",j);
      qcdp_[j] = new CPepeModel2(nname,pfmet);
	  
	  sprintf(nname, "qcdm2d%d",j);
      qcdm_[j] = new CPepeModel2(nname,pfmet);
  
  }
   std::cout << "just finished making the pepes" << std::endl;
   
   // create the appropriate gaussian constraint for the signal region electroweak+ttbar 
    RooGaussian constm_sr("constm_sr","constm_sr",*nEWKm_[0],RooConst(nEWKm_[0]->getVal()),RooConst(0.10*nEWKm_[0]->getVal()));
    RooGaussian constp_sr("constp_sr","constp_sr",*nEWKp_[0],RooConst(nEWKp_[0]->getVal()),RooConst(0.10*nEWKp_[0]->getVal()));
	
	RooGaussian constm_cr("constm_cr","constm_cr",*nEWKm_[1],RooConst(nEWKm_[1]->getVal()),RooConst(0.10*nEWKm_[1]->getVal()));
    RooGaussian constp_cr("constp_cr","constp_cr",*nEWKp_[1],RooConst(nEWKp_[1]->getVal()),RooConst(0.10*nEWKp_[1]->getVal()));
	
	// RooGaussian constm_sr("constm_sr","constm_sr",*dewkm_[0],RooConst(dewkm_[0]->getVal()),RooConst(0.15*dewkm_[0]->getVal()));
    // RooGaussian constp_sr("constp_sr","constp_sr",*dewkp_[0],RooConst(dewkp_[0]->getVal()),RooConst(0.15*dewkp_[0]->getVal()));
	
	
	     std::cout << "nsigp " << nSigp_[0]->getVal() << std::endl;
		 std::cout << "nttbp " << nTtbp_[0]->getVal() << std::endl;
		 std::cout << "nqcdp " << nQCDp_[0]->getVal() << std::endl;
		 std::cout << "nwxp " << nWxp_[0]->getVal() << std::endl;
		 std::cout << "nzxxp " << nZxxp_[0]->getVal() << std::endl;
	
		// Fill the actual values later
	// Wx sample - constrain dewk values w/ Gaussian
	RooGaussian const_wxm("const_wxm","const_wxm",*dwxm_[0],RooConst(dwxm_[0]->getVal()),RooConst(0.009*dwxm_[0]->getVal()));
    RooGaussian const_wxp("const_wxp","const_wxp",*dwxp_[0],RooConst(dwxp_[0]->getVal()),RooConst(0.009*dwxp_[0]->getVal()));
	
	RooGaussian const_ttbm("const_ttbm","const_ttbm",*nTtbm_[0],RooConst(nTtbm_[0]->getVal()*815/830),RooConst(0.0686*nTtbm_[0]->getVal()));
    RooGaussian const_ttbp("const_ttbp","const_ttbp",*nTtbp_[0],RooConst(nTtbp_[0]->getVal()*815/830),RooConst(0.0686*nTtbp_[0]->getVal()));
	
	RooGaussian const_dibm("const_dibm","const_dibm",*nDibm_[0],RooConst(nDibm_[0]->getVal()),RooConst(0.009*nDibm_[0]->getVal()));
    RooGaussian const_dibp("const_dibp","const_dibp",*nDibp_[0],RooConst(nDibp_[0]->getVal()),RooConst(0.009*nDibp_[0]->getVal()));
	
	RooGaussian const_zxxm("const_zxxm","const_zxxm",*nZxxm_[0],RooConst(nZxxm_[0]->getVal()*5688/5833),RooConst(0.0449*nZxxm_[0]->getVal()));
    RooGaussian const_zxxp("const_zxxp","const_zxxp",*nZxxp_[0],RooConst(nZxxp_[0]->getVal()*5688/5833),RooConst(0.0449*nZxxp_[0]->getVal()));
	
	
	
	std::cout << "make constraints" << std::endl;

	// RooGaussian constm_cr1("constm_cr1","constm_cr1",*nEWKm_[1],RooConst(nEWKm_[1]->getVal()),RooConst(0.05*nEWKm_[1]->getVal()));
    // RooGaussian constp_cr1("constp_cr1","constp_cr1",*nEWKp_[1],RooConst(nEWKp_[1]->getVal()),RooConst(0.05*nEWKp_[1]->getVal()));

	// RooGaussian constm_cr2("constm_cr2","constm_cr2",*nEWKm_[2],RooConst(nEWKm_[2]->getVal()),RooConst(0.05*nEWKm_[2]->getVal()));
    // RooGaussian constp_cr2("constp_cr2","constp_cr2",*nEWKp_[2],RooConst(nEWKp_[2]->getVal()),RooConst(0.05*nEWKp_[2]->getVal()));

	// RooGaussian constm_cr3("constm_cr3","constm_cr3",*nEWKm_[3],RooConst(nEWKm_[3]->getVal()),RooConst(0.05*nEWKm_[3]->getVal()));
    // RooGaussian constp_cr3("constp_cr3","constp_cr3",*nEWKp_[3],RooConst(nEWKp_[3]->getVal()),RooConst(0.05*nEWKp_[3]->getVal()));

	// RooGaussian constm_cr4("constm_cr4","constm_cr4",*nEWKm_[4],RooConst(nEWKm_[4]->getVal()),RooConst(0.05*nEWKm_[4]->getVal()));
    // RooGaussian constp_cr4("constp_cr4","constp_cr4",*nEWKp_[4],RooConst(nEWKp_[4]->getVal()),RooConst(0.05*nEWKp_[4]->getVal()));

	// RooGaussian constm_cr5("constm_cr5","constm_cr5",*nEWKm_[5],RooConst(nEWKm_[5]->getVal()),RooConst(0.05*nEWKm_[5]->getVal()));
    // RooGaussian constp_cr5("constp_cr5","constp_cr5",*nEWKp_[5],RooConst(nEWKp_[5]->getVal()),RooConst(0.05*nEWKp_[5]->getVal()));

	
     // contstuct the PDFs for the 2 test isolation bins
  // Test the vector map 
  // // map <string, RooAddPdf> pdfMetp_;
   // map <int, RooAddPdf*> pdfMetp_;
   vector < RooAddPdf*> pdfMetp_(nIsoBins), pdfMetm_(nIsoBins);
  // make the pdfs?
  for(int j = 0; j < nIsoBins; ++j){
	 // if(j==0 || j==1){
	 if(j==0||j==1||j==2||j==3||j==4){
		 sprintf(nname,"pdfWep%d",j);

		 if(doMET&&(!doTemplate)){
             pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j]));
         } else {
		     pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*pdfQCDp_[j]),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j]));
         }
		 // pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfWxp_[j],*pdfZxxp_[j],*pdfDibp_[j],*pdfTtbp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nWxp_[j],*nZxxp_[j],*nDibp_[j],*nTtbp_[j],*nQCDp_[j]));
		 std::cout << "blah" << std::endl;
		 
		 sprintf(nname,"pdfWem%d",j);
		 if(doMET&&(!doTemplate)){
             pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j]));
         }else {
             pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*pdfQCDm_[j]),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j]));
         }		
        // pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfWxm_[j],*pdfZxxm_[j],*pdfDibm_[j],*pdfTtbm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nWxm_[j],*nZxxm_[j],*nDibm_[j],*nTtbm_[j],*nQCDm_[j]));
		 // sprintf(nname,"pdfWep%d",j);
		 // pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfDibp_[j],*pdfTtbp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nDibp_[j],*nTtbp_[j],*nQCDp_[j]));
		 
		 // sprintf(nname,"pdfWem%d",j);
		 // pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfDibm_[j],*pdfTtbm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nDibm_[j],*nTtbm_[j]*nQCDm_[j]));
	 } else {
         sprintf(nname,"pdfWep%d",j);
		 if(doMET&&(!doTemplate)){
             pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*(qcdp_[j]->model)),RooArgList(*nQCDp_[j]));
		 } else {
             pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfQCDp_[j]),RooArgList(*nQCDp_[j]));
		 }
         
		 sprintf(nname,"pdfWem%d",j);
		 if (doMET&&(!doTemplate) ) {
             pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*(qcdm_[j]->model)),RooArgList(*nQCDm_[j]));
         } else {
		     pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfQCDm_[j]),RooArgList(*nQCDm_[j]));
         }
		 
	 }
  }
  
  std::cout << "made pdfs" << std::endl;
  
  
  // PDF for simultaneous fit  
  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("Select");
  rooCat.defineType("Anti");
  
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(pdfMet, "Select");
  //pdfTotal.addPdf(apdfMet,"Anti");
  
  RooSimultaneous pdfTotalp("pdfTotalp","pdfTotalp",rooCat);
  pdfTotalp.addPdf(pdfMetp, "Select");
  pdfTotalp.addPdf(apdfMetp,"Anti");
  
  RooSimultaneous pdfTotalm("pdfTotalm","pdfTotalm",rooCat);
  pdfTotalm.addPdf(pdfMetm, "Select");
  pdfTotalm.addPdf(apdfMetm,"Anti");

  //
  // Perform fits
  //
  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);
 
  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
  RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
  RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);

//   RooDataHist dataTotalp("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat),
// 			 Import("Select", dataMetp),
// 			 Import("Anti",   antiMetp));
// 
//   RooDataHist dataTotalm("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat),
// 			 Import("Select", dataMetm),
// 			 Import("Anti", antiMetm));


  
    // ------------------------------------------------------------------------
  // Make the data histograms for the RooCategory Fit
  
  // put them into a vector map, make sure the names match the RooCategory names
  map<string, RooDataHist*> dataMetp_, dataMetm_;
  std::vector<RooDataHist*> dataMetpHist_(nIsoBins),dataMetmHist_(nIsoBins);
  // map<string, RooDataHist*> dataMetm_;
  for(int i = 0; i < nIsoBins; ++i){
	  sprintf(nname, "isop%d",i);
	  dataMetp_[nname]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetpIsoBins[i]);
      sprintf(nname, "isom%d",i);
	  dataMetm_[nname]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetmIsoBins[i]);
      sprintf(nname, "isop%d",i);
      dataMetpHist_[i]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetpIsoBins[i]);
	  sprintf(nname, "isom%d",i);
      dataMetmHist_[i]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetmIsoBins[i]);
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  // set up RooCategory for data
  // test only 2 bins at first
  
  RooCategory rooCat2dTest("rooCat2dTest","rooCat2dTest");
  for (int i = 0; i < nIsoBins; ++ i){
	  sprintf(nname,"iso%d",i); rooCat2dTest.defineType(nname);
  }

  RooDataHist combDatap("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetp_));
  RooDataHist combDatam("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetm_));
  
  RooSimultaneous simPdfp("simPdfp","simultaneous pdf W+",rooCat2dTest) ;
  RooSimultaneous simPdfm("simPdfm","simultaneous pdf W-",rooCat2dTest) ;
  for(int i = 0; i < nIsoBins; ++i ){
	  sprintf(nname, "isop%d", i); simPdfp.addPdf(*pdfMetp_[i],nname);
	  sprintf(nname, "isom%d", i); simPdfm.addPdf(*pdfMetm_[i],nname);
  }

  
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
  
  
  cout << "Starting values for AntiWmunu yields: " << endl;
  cout << "Selected: " << hAntiDataMet->Integral() << endl;
  cout << "   sig: " << hAntiWmunuMet->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMet->Integral() << endl;
  cout << "   qcd: " << hAntiDataMet->Integral()-hAntiWmunuMet->Integral()-hAntiEWKMet->Integral() << endl;

  cout << "Starting values for AntiWmunu_p yields: " << endl;
  cout << "   sig: " << hWmunuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetp->Integral()-hAntiWmunuMetp->Integral()-hAntiEWKMetp->Integral() << endl;

  cout << "Starting values for AntiWmunu_m yields: " << endl;
  cout << "   sig: " << hAntiWmunuMetm->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetm->Integral()-hAntiWmunuMetm->Integral()-hAntiEWKMetm->Integral() << endl;

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());

  RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  
  RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  
  
  
  std::vector<RooRealVar*> pepe2Pdf_qcdp_norms_(nIsoBins), pepe2Pdf_qcdm_norms_(nIsoBins);
  for(int j = 0; j < nIsoBins; ++j){
     // sprintf(nname,"pepe2Pdf_qcdp_norm");
     if(j==0){
       sprintf(nname,"pepe2Pdf_qcdp2d%i_norm",j);
       pepe2Pdf_qcdp_norms_[j] = new RooRealVar(nname,nname,0.4*hDataMetpIsoBins[j]->Integral(),0,0.6*hDataMetpIsoBins[j]->Integral());
       sprintf(nname,"pepe2Pdf_qcdm2d%i_norm",j);
       pepe2Pdf_qcdm_norms_[j] = new RooRealVar(nname,nname,0.4*hDataMetmIsoBins[j]->Integral(),0,0.6*hDataMetmIsoBins[j]->Integral());
     } else {
       sprintf(nname,"pepe2Pdf_qcdp2d%i_norm",j);
       pepe2Pdf_qcdp_norms_[j] = new RooRealVar(nname,nname,hDataMetpIsoBins[j]->Integral(),0.95*hDataMetpIsoBins[j]->Integral(),1.05*hDataMetpIsoBins[j]->Integral());
       sprintf(nname,"pepe2Pdf_qcdm2d%i_norm",j);
       pepe2Pdf_qcdm_norms_[j] = new RooRealVar(nname,nname,hDataMetmIsoBins[j]->Integral(),0.95*hDataMetmIsoBins[j]->Integral(),1.05*hDataMetmIsoBins[j]->Integral());
     }
  }
  std::cout << "mark 3 " << std::endl;
  
//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",40000,0,100000);
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",40000,0,100000);

  RooWorkspace combine_workspace("combine_workspace");
  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);
  combine_workspace.import(pepe2Pdf_qcdp_norm);
  combine_workspace.import(pepe2Pdf_qcdm_norm);

  combine_workspace.import(pdfWm);
  combine_workspace.import(pdfWmp);
  combine_workspace.import(pdfWmm);
  combine_workspace.import(pdfWm_RecoilUp);
  combine_workspace.import(pdfWmp_RecoilUp);
  combine_workspace.import(pdfWmm_RecoilUp);
  combine_workspace.import(pdfWm_RecoilDown);
  combine_workspace.import(pdfWmp_RecoilDown);
  combine_workspace.import(pdfWmm_RecoilDown);
  combine_workspace.import(pdfWm_ScaleUp);
  combine_workspace.import(pdfWmp_ScaleUp);
  combine_workspace.import(pdfWmm_ScaleUp);
  combine_workspace.import(pdfWm_ScaleDown);
  combine_workspace.import(pdfWmp_ScaleDown);
  combine_workspace.import(pdfWmm_ScaleDown);
  combine_workspace.import(pdfWm_PileupUp);
  combine_workspace.import(pdfWmp_PileupUp);
  combine_workspace.import(pdfWmm_PileupUp);
  combine_workspace.import(pdfWm_PileupDown);
  combine_workspace.import(pdfWmp_PileupDown);
  combine_workspace.import(pdfWmm_PileupDown);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKm);
  combine_workspace.import(pdfEWK_PileupUp);
  combine_workspace.import(pdfEWKp_PileupUp);
  combine_workspace.import(pdfEWKm_PileupUp);
  combine_workspace.import(pdfEWK_PileupDown);
  combine_workspace.import(pdfEWKp_PileupDown);
  combine_workspace.import(pdfEWKm_PileupDown);
//   combine_workspace.import(*(qcd.model));
//   combine_workspace.import(*(qcdp.model));
//   combine_workspace.import(*(qcdm.model));

  combine_workspace.import(pdfQCD);
  combine_workspace.import(pdfQCDp);
  combine_workspace.import(pdfQCDm);
  
  combine_workspace.import(antiMet);
  combine_workspace.import(antiMetp);
  combine_workspace.import(antiMetm);
  combine_workspace.import(pepe2Pdf_aqcdp_norm);
  combine_workspace.import(pepe2Pdf_aqcdm_norm);
  
  combine_workspace.import(apdfWm);
  combine_workspace.import(apdfWmp);
  combine_workspace.import(apdfWmm);
  combine_workspace.import(apdfEWK);
  combine_workspace.import(apdfEWKp);
  combine_workspace.import(apdfEWKm);
  combine_workspace.import(*(aqcd.model));
  //combine_workspace.import(qcdpn);
  combine_workspace.import(*(aqcdp.model));
  combine_workspace.import(*(aqcdm.model));

  sprintf(nname,"%s/Wmunu_pdfTemplates.root",CPlot::sOutDir.Data());
  combine_workspace.writeToFile(nname);
  
  
   // separate file for the binned workspace just so i don't fuck up the original one
  RooWorkspace binned_workspace("binned_workspace");
  // loop through the number of bins, import the appropriate pdf for each one
  for(int j = 0; j < nIsoBins; ++j){
      // binned_workspace.import(dataMet);
      binned_workspace.import(*dataMetpHist_[j]);
      binned_workspace.import(*dataMetmHist_[j]);
	  // QCD normalization RooRealVars
      binned_workspace.import(*pepe2Pdf_qcdp_norms_[j]);
      binned_workspace.import(*pepe2Pdf_qcdm_norms_[j]);
	  // QCD shapes
      // binned_workspace.import(*(qcd.model));
      binned_workspace.import(*(qcdp_[j]->model));
      binned_workspace.import(*(qcdm_[j]->model));
      std::cout << "importing   " << qcdm_[j]->model->GetName() << std::endl;
	  //The DataHist
	  // combine_workspace.import(*wenuMet_[j]);
	  binned_workspace.import(*wmunuMetp_[j]);
	  binned_workspace.import(*wmunuMetm_[j]);
      // binned_workspace.import(pdfWe);
      // binned_workspace.import(*pdfWep_[j]);
      // binned_workspace.import(*pdfWem_[j]);
	  // The DataHist
	  // combine_workspace.import(*ewkMet_[j]);
	  binned_workspace.import(*ewkMetp_[j]);
	  binned_workspace.import(*ewkMetm_[j]);
      
      binned_workspace.import(*wxMetp_[j]);
	  binned_workspace.import(*wxMetm_[j]);
      
      binned_workspace.import(*zxxMetp_[j]);
	  binned_workspace.import(*zxxMetm_[j]);
      
      binned_workspace.import(*ttbMetp_[j]);
	  binned_workspace.import(*ttbMetm_[j]);
      
      binned_workspace.import(*dibMetp_[j]);
	  binned_workspace.import(*dibMetm_[j]);
	  //The Pdfs
      // binned_workspace.import(pdfEWK);
      // binned_workspace.import(*pdfEWKp_[j]);
      // binned_workspace.import(*pdfEWKm_[j]);
	  
	  binned_workspace.import(*hDataMetpIsoBins[j]);
	  binned_workspace.import(*hDataMetmIsoBins[j]);
	  binned_workspace.import(*hWmunuMetpIsoBins[j]);
	  binned_workspace.import(*hWmunuMetmIsoBins[j]);
	  binned_workspace.import(*hEWKMetpIsoBins[j]);
	  binned_workspace.import(*hEWKMetmIsoBins[j]);
      binned_workspace.import(*hWxMetpIsoBins[j]);
	  binned_workspace.import(*hWxMetmIsoBins[j]);
      binned_workspace.import(*hZxxMetpIsoBins[j]);
	  binned_workspace.import(*hZxxMetmIsoBins[j]);
      binned_workspace.import(*hDibMetpIsoBins[j]);
	  binned_workspace.import(*hDibMetmIsoBins[j]);
      binned_workspace.import(*hTtbMetpIsoBins[j]);
	  binned_workspace.import(*hTtbMetmIsoBins[j]);
	  // binned_workspace.import(*hDataMetpIsoBins[j]);
	  

  }// end of loop, now save the workspace
  sprintf(nname, "%s/Wmunu_pdfTemplates_binned.root",CPlot::sOutDir.Data());
  binned_workspace.writeToFile(nname);
  
  

  
  
     // had commented out the Min2/Strat2 when running some of the other options
  // RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),RooFit::Strategy(2)/*,Minimizer("Minuit2","minimize")*/,ExternalConstraints(RooArgList(const_wxp,const_zxxp,const_dibp,const_ttbp)),PrintEvalErrors(-1));
  // RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),RooFit::Strategy(2)/*,Minimizer("Minuit2","minimize")*/,ExternalConstraints(RooArgList(const_wxm,const_zxxm,const_dibm,const_ttbm)),PrintEvalErrors(-1));
  
  RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),ExternalConstraints(RooArgSet(constp_sr/*,constp_cr*/)),RooFit::Strategy(2),Minos(kTRUE),/*Minimizer("Minuit2","minimize"),*/PrintEvalErrors(-1));
  RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),ExternalConstraints(RooArgSet(constm_sr/*,constm_cr*/)),RooFit::Strategy(2),Minos(kTRUE),/*Minimizer("Minuit2","minimize"),*/PrintEvalErrors(-1));
  
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle);
  hMetmDiff->SetMarkerSize(0.9);
  
  // the diff hists for the anti-selection
  
  TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  
  hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetp->GetNbinsX(); ++ibin){hPdfAntiMetp->SetBinError(ibin, hAntiWmunuMetp->GetBinError(ibin));}
  hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetm->GetNbinsX(); ++ibin){hPdfAntiMetm->SetBinError(ibin, hAntiWmunuMetm->GetBinError(ibin));}
  hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  hAntiMetmDiff->SetMarkerStyle(kFullCircle);
  hAntiMetmDiff->SetMarkerSize(0.9);
    
    char plotname2[100];
    //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
    // label for lumi
  // char lumitext[100];
  // if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  (8 TeV)",lumi*1000.);
  // else         sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000.);
  
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
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.02);
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
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  (8 TeV)",lumi*1000.);
  else         sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000.);

  
    
  Double_t chi2probm, chi2ndfm, chi2probp, chi2ndfp;
  Double_t ksprobm, ksprobpem, ksprobp, ksprobpep;
  
    ofstream txtfile4;
    char txtfname4[100];
    std::cout << "Printing We+. " << std::endl;
    sprintf(txtfname4,"%s/chi2.txt",CPlot::sOutDir.Data());
    txtfile4.open(txtfname4);
    assert(txtfile4.is_open());
  
  for(int i = 0; i < nIsoBins; ++i){
      
    // set up the diff plot
    // turn this into its own function later?
    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetp = (TH1D*)(pdfMetp_[i]->createHistogram("hPdfMetp", pfmet));
    for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetpIsoBins[i]->GetBinError(ibin));}
    hPdfMetp->Scale((nSigp_[i]->getVal()+nEWKp_[i]->getVal()+nQCDp_[i]->getVal())/hPdfMetp->Integral());
    TH1D *hMetpDiff = makeDiffHist(hDataMetpIsoBins[i],hPdfMetp,"hMetpDiff");
    hMetpDiff->SetMarkerStyle(kFullCircle);
    hMetpDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    // turn this part also into its own function later
    RooPlot *wepframe = pfmet.frame(Bins(NBINS)); 
    wepframe->GetYaxis()->SetNdivisions(505);
    wepframe->GetXaxis()->SetLabelOffset(2.0);
    sprintf(nname,"isop%d",i);
    sprintf(plotname2, "wep_fitmetp_bin%i",i);
    dataMetp_[nname]->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMetp_[i]->plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMetp_[i]->plotOn(wepframe,LineColor(linecolorW));
    if(!doMET||doTemplate){
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*pdfQCDp_[i])),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*pdfQCDp_[i])),LineColor(linecolorEWK));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfQCDp_[i])),LineColor(linecolorQCD));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfQCDp_[i])),FillColor(fillcolorQCD),DrawOption("F"));
    } else {
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*(qcdp_[i]->model))),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfEWKp_[i],*(qcdp_[i]->model))),LineColor(linecolorEWK));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*(qcdp_[i]->model))),FillColor(fillcolorQCD),DrawOption("F"));
      pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*(qcdp_[i]->model))),LineColor(linecolorQCD));
    } 
    pdfMetp_[i]->plotOn(wepframe,Components(RooArgSet(*pdfWep_[i])),LineColor(linecolorW),LineStyle(2));
    dataMetp_[nname]->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
   
    sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
    CPlot plotMetp(plotname2,wepframe,"","",ylabel);
    plotMetp.SetLegend(0.68,0.57,0.93,0.77);
    plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMetp.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
    plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMetp.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
    plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
    plotMetp.Draw(c,kFALSE,format,1);
    
    std::cout << "Draw the W plot diff" << std::endl;
    CPlot plotMetpDiff(plotname2,"","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
    hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
    hMetpDiff->GetYaxis()->SetLabelSize(0.11);
    plotMetpDiff.SetYRange(-yscale,yscale);
    plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
    plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
    plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
    plotMetpDiff.Draw(c,kTRUE,format,2);
    plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
    std::cout << "Draw the W plot log" << std::endl;
    sprintf(plotname2,"wep_fitmetp_bin%i_log",i);
    plotMetp.SetName(plotname2);
    plotMetp.SetLogy();
    plotMetp.SetYRange(1e-5*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
    plotMetp.Draw(c,kTRUE,format,1);
    plotMetp.Draw(c,kTRUE,"pdf",1);
   
    chi2probp = hDataMetpIsoBins[i]->Chi2Test(hPdfMetp,"PUW");
    chi2ndfp  = hDataMetpIsoBins[i]->Chi2Test(hPdfMetp,"CHI2/NDFUW");
    ksprobp   = hDataMetpIsoBins[i]->KolmogorovTest(hPdfMetp);
    ksprobpep = hDataMetpIsoBins[i]->KolmogorovTest(hPdfMetp,"DX"); 
   
    // Do the W-  plots here
    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetm = (TH1D*)(pdfMetm_[i]->createHistogram("hPdfMetm", pfmet));
    for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetmIsoBins[i]->GetBinError(ibin));}
    hPdfMetm->Scale((nSigm_[i]->getVal()+nEWKm_[i]->getVal()+nQCDm_[i]->getVal())/hPdfMetm->Integral());
    TH1D *hMetmDiff = makeDiffHist(hDataMetmIsoBins[i],hPdfMetm,"hMetmDiff");
    hMetmDiff->SetMarkerStyle(kFullCircle);
    hMetmDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    // turn this part also into its own function later
    RooPlot *wemframe = pfmet.frame(Bins(NBINS)); 
    wemframe->GetYaxis()->SetNdivisions(505);
    wemframe->GetXaxis()->SetLabelOffset(2.0);
    sprintf(nname,"isom%d",i);
    sprintf(plotname2, "wem_fitmetm_bin%i",i);
    dataMetm_[nname]->plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    pdfMetm_[i]->plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
    pdfMetm_[i]->plotOn(wemframe,LineColor(linecolorW));
    if(!doMET||doTemplate){
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*pdfQCDm_[i])),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*pdfQCDm_[i])),LineColor(linecolorEWK));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfQCDm_[i])),FillColor(fillcolorQCD),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfQCDm_[i])),LineColor(linecolorQCD));
    } else {
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*(qcdm_[i]->model))),FillColor(fillcolorEWK),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfEWKm_[i],*(qcdm_[i]->model))),LineColor(linecolorEWK));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*(qcdm_[i]->model))),FillColor(fillcolorQCD),DrawOption("F"));
      pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*(qcdm_[i]->model))),LineColor(linecolorQCD));
    }
    pdfMetm_[i]->plotOn(wemframe,Components(RooArgSet(*pdfWem_[i])),LineColor(linecolorW),LineStyle(2));
    dataMetm_[nname]->plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
   
    sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
    CPlot plotMetm(plotname2,wemframe,"","",ylabel);
    plotMetm.SetLegend(0.68,0.57,0.93,0.77);
    plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
    plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotMetm.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
    plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
    plotMetm.Draw(c,kFALSE,format,1);
    
    std::cout << "Draw the W plot diff" << std::endl;
    CPlot plotMetmDiff(plotname2,"","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
    hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
    hMetmDiff->GetYaxis()->SetLabelSize(0.11);
    plotMetmDiff.SetYRange(-yscale,yscale);
    plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
    plotMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
    plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
    plotMetmDiff.Draw(c,kTRUE,format,2);
    plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
    std::cout << "Draw the W plot log" << std::endl;
    sprintf(plotname2,"wem_fitmetm_bin%i_log",i);
    plotMetm.SetName(plotname2);
    plotMetm.SetLogy();
    plotMetm.SetYRange(1e-5*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
    plotMetm.Draw(c,kTRUE,format,1);
    plotMetm.Draw(c,kTRUE,"pdf",1);
	
	chi2probm = hDataMetmIsoBins[i]->Chi2Test(hPdfMetm,"PUW");
    chi2ndfm  = hDataMetmIsoBins[i]->Chi2Test(hPdfMetm,"CHI2/NDFUW");
    ksprobm   = hDataMetmIsoBins[i]->KolmogorovTest(hPdfMetm);
    ksprobpem = hDataMetmIsoBins[i]->KolmogorovTest(hPdfMetm,"DX"); 
	
    txtfile4 << "bin+ : " << i << " ; chi2 : " << chi2ndfp << endl;
	txtfile4 << "bin- : " << i << " ; chi2 : " << chi2ndfm << endl;
	
  
   std::cout << "blah " << i << std::endl;
  
  }
  
  txtfile4.close();
  
  std::cout << "hi" << std::endl;
  ofstream txtfile;
  char txtfname1[100];
  std::cout << "Printing We+. " << std::endl;
  sprintf(txtfname1,"%s/fitresWe2dCatTest.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname1);
  assert(txtfile.is_open());
  
  // flags = txtfile.flags();
  // txtfile << setprecision(10);
  // txtfile << " *** Yields *** " << endl;
  txtfile << "Data+: " << hDataMetpIsoBins[0]->Integral() << endl;
  txtfile << "Data-: " << hDataMetmIsoBins[0]->Integral() << endl;
  // txtfile << "  Signal: " << nSigp.getVal() << " +/- " << nSigp.getPropagatedError(*fitResp) << endl;
  // txtfile << "     QCD: " << nQCDp.getVal() << " +/- " << nQCDp.getPropagatedError(*fitResp) << endl;
  // txtfile << "   Other: " << nEWKp.getVal() << " +/- " << nEWKp.getPropagatedError(*fitResp) << endl;
  // txtfile << "AntiSelected: " << hAntiDataMetp->Integral() << endl;
  // txtfile << "  AntiSignal: " << nAntiSigp.getVal() << " +/- " << nAntiSigp.getPropagatedError(*fitResp) << endl;
  // txtfile << "     AntiQCD: " << nAntiQCDp.getVal() << " +/- " << nAntiQCDp.getPropagatedError(*fitResp) << endl;
  // txtfile << "   AntiOther: " << nAntiEWKp.getVal() << " +/- " << nAntiEWKp.getPropagatedError(*fitResp) << endl;
  // txtfile << endl;  
  // txtfile.flags(flags);
  
   
  
  fitResp2dCatTest->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose); 
  txtfile << endl;
  printCorrelations(txtfile, fitResp2dCatTest);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2probp, chi2ndfp, ksprobp, ksprobpep);
  txtfile<< endl;
  
  fitResm2dCatTest->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitResm2dCatTest);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2probm, chi2ndfm, ksprobm, ksprobpem);
  txtfile <<endl;
  
  txtfile.close();
  
  
  // TCanvas *c2 = MakeCanvas("c2","c2",800,600);
  // xframe->GetYaxis()->SetMaximum(0)
  // // // test likelihood plotting --------------
    // // // // Construct Likelihood
    // RooAbsReal* nllm = simPdfm.createNLL(combDatam,NumCPU(2)) ;
    // RooAbsReal* nllp = simPdfp.createNLL(combDatap,NumCPU(2)) ;

   // // // // Minimize likelihood w.r.t all parameters before making plots
    // // // // RooMinuit(*nllm).SetPrintLevel(-1) ;
    // // // // RooMinuit(*nllp).SetPrintLevel(-1) ;
    // RooMinuit(*nllm).migrad() ;
    // RooMinuit(*nllp).migrad() ;

   // // // // // Plot likelihood scan frac 
   // RooPlot* framesigmll = nSigm_[0]->frame(Bins(10),Range(0.9*nSigm_[0]->getVal(),1.1*nSigm_[0]->getVal()),Title("")) ; 
   // RooPlot* frameewkmll = nEWKm_[0]->frame(Bins(10),Range(0.9*nEWKm_[0]->getVal(),1.1*nEWKm_[0]->getVal()),Title("")) ; //frameewkmll->SetMinimum(0); frameewkmll->SetMaximum(100);
   // RooPlot* frameqcdmll = nQCDm_[0]->frame(Bins(10),Range(0.9*nQCDm_[0]->getVal(),1.1*nQCDm_[0]->getVal()),Title("")) ; //frameqcdmll->SetMinimum(0); frameqcdmll->SetMaximum(100);
   // RooPlot* frameqcdmc1ll = (qcdm_[0]->c1)->frame(Bins(10),Title("")) ; //frameqcdmc1ll->SetMinimum(0); frameqcdmc1ll->SetMaximum(100);
   // RooPlot* frameqcdmc2ll = (qcdm_[0]->c2)->frame(Bins(10),Title("")) ; //frameqcdmc2ll->SetMinimum(0); frameqcdmc2ll->SetMaximum(100);
   // RooPlot* frameqcdmc3ll = (qcdm_[0]->c3)->frame(Bins(10),Title("")) ; //frameqcdmc3ll->SetMinimum(0); frameqcdmc3ll->SetMaximum(100);
   // RooPlot* frameqcdmd1ll = (qcdm_[0]->d1)->frame(Bins(10),Title("")) ; //frameqcdmd1ll->SetMinimum(0); frameqcdmd1ll->SetMaximum(100);
   // RooPlot* frameqcdmd2ll = (qcdm_[0]->d2)->frame(Bins(10),Title("")) ; //frameqcdmd2ll->SetMinimum(0); frameqcdmd2ll->SetMaximum(100);
   // RooPlot* frameqcdmd3ll = (qcdm_[0]->d3)->frame(Bins(10),Title("")) ; //frameqcdmd3ll->SetMinimum(0); frameqcdmd3ll->SetMaximum(100);
   
   // RooPlot* frameqcdmll1 = nQCDm_[1]->frame(Bins(10),Range(0.9*nQCDm_[1]->getVal(),1.1*nQCDm_[1]->getVal()),Title("")) ; //frameqcdmll1->SetMinimum(0); frameqcdmll1->SetMaximum(100);
   // RooPlot* frameqcdmll2 = nQCDm_[2]->frame(Bins(10),Range(0.9*nQCDm_[2]->getVal(),1.1*nQCDm_[2]->getVal()),Title("")) ;// frameqcdmll2->SetMinimum(0); frameqcdmll2->SetMaximum(100);
   // RooPlot* frameqcdmll3 = nQCDm_[3]->frame(Bins(10),Range(0.9*nQCDm_[3]->getVal(),1.1*nQCDm_[3]->getVal()),Title("")) ; //frameqcdmll3->SetMinimum(0); frameqcdmll3->SetMaximum(100);
   // RooPlot* frameqcdmll4 = nQCDm_[4]->frame(Bins(10),Range(0.9*nQCDm_[4]->getVal(),1.1*nQCDm_[4]->getVal()),Title("")) ; //frameqcdmll4->SetMinimum(0); frameqcdmll4->SetMaximum(100);
   // RooPlot* frameqcdmll5 = nQCDm_[5]->frame(Bins(10),Range(0.9*nQCDm_[5]->getVal(),1.1*nQCDm_[5]->getVal()),Title("")) ; //frameqcdmll5->SetMinimum(0); frameqcdmll5->SetMaximum(100);
   // // //nll->plotOn(framell,ShiftToZero()) ;
   

   // // // //nll->plotOn(framell,ShiftToZero()) ;
   
    // RooAbsReal* pll_nsigm = nllm->createProfile(*nSigm_[0]) ;
    // RooAbsReal* pll_newkm = nllm->createProfile(*nEWKm_[0]) ;
    // RooAbsReal* pll_nqcdm = nllm->createProfile(*nQCDm_[0]) ;
    // RooAbsReal* pll_qcdmc1 = nllm->createProfile(*(qcdm_[0]->c1)) ;
    // RooAbsReal* pll_qcdmc2 = nllm->createProfile(*(qcdm_[0]->c2)) ;
    // RooAbsReal* pll_qcdmc3 = nllm->createProfile(*(qcdm_[0]->c3)) ;
    // RooAbsReal* pll_qcdmd1 = nllm->createProfile(*(qcdm_[0]->d1)) ;
    // RooAbsReal* pll_qcdmd2 = nllm->createProfile(*(qcdm_[0]->d2)) ;
    // RooAbsReal* pll_qcdmd3 = nllm->createProfile(*(qcdm_[0]->d3)) ;
	
	// RooAbsReal* pll_nqcdm1 = nllm->createProfile(*nQCDm_[1]) ;
	// RooAbsReal* pll_nqcdm2 = nllm->createProfile(*nQCDm_[2]) ;
	// RooAbsReal* pll_nqcdm3 = nllm->createProfile(*nQCDm_[3]) ;
	// RooAbsReal* pll_nqcdm4 = nllm->createProfile(*nQCDm_[4]) ;
	// RooAbsReal* pll_nqcdm5 = nllm->createProfile(*nQCDm_[5]) ;
	
	// pll_nsigm->plotOn(framesigmll,LineColor(kRed)) ; framesigmll->SetMinimum(0); framesigmll->SetMaximum(10);
    // pll_newkm->plotOn(frameewkmll,LineColor(kRed)) ; frameewkmll->SetMinimum(0); frameewkmll->SetMaximum(10);
    // pll_nqcdm->plotOn(frameqcdmll,LineColor(kRed)) ; frameqcdmll->SetMinimum(0); frameqcdmll->SetMaximum(10);
    // pll_qcdmc1->plotOn(frameqcdmc1ll,LineColor(kRed)) ; frameqcdmc1ll->SetMinimum(0); frameqcdmc1ll->SetMaximum(10);
    // pll_qcdmc2->plotOn(frameqcdmc2ll,LineColor(kRed)) ; frameqcdmc2ll->SetMinimum(0); frameqcdmc2ll->SetMaximum(10);
    // pll_qcdmc3->plotOn(frameqcdmc3ll,LineColor(kRed)) ; frameqcdmc3ll->SetMinimum(0); frameqcdmc3ll->SetMaximum(10);
    // pll_qcdmd1->plotOn(frameqcdmd1ll,LineColor(kRed)) ; frameqcdmd1ll->SetMinimum(0); frameqcdmd1ll->SetMaximum(10);
    // pll_qcdmd2->plotOn(frameqcdmd2ll,LineColor(kRed)) ; frameqcdmd2ll->SetMinimum(0); frameqcdmd2ll->SetMaximum(10);
    // pll_qcdmd3->plotOn(frameqcdmd3ll,LineColor(kRed)) ; frameqcdmd3ll->SetMinimum(0); frameqcdmd3ll->SetMaximum(10);
	
	// pll_nqcdm1->plotOn(frameqcdmll1,LineColor(kRed)) ; frameqcdmll1->SetMinimum(0); frameqcdmll1->SetMaximum(10);
	// pll_nqcdm2->plotOn(frameqcdmll2,LineColor(kRed)) ; frameqcdmll2->SetMinimum(0); frameqcdmll2->SetMaximum(10);
	// pll_nqcdm3->plotOn(frameqcdmll3,LineColor(kRed)) ; frameqcdmll3->SetMinimum(0); frameqcdmll3->SetMaximum(10);
	// pll_nqcdm4->plotOn(frameqcdmll3,LineColor(kRed)) ; frameqcdmll3->SetMinimum(0); frameqcdmll3->SetMaximum(10);
	// pll_nqcdm5->plotOn(frameqcdmll5,LineColor(kRed)) ; frameqcdmll5->SetMinimum(0); frameqcdmll5->SetMaximum(10);
	
	   // char plotname[100];
    // // // Adjust frame maximum for visual clarity
   // // // // framell->SetMinimum(0) ;
   // // // // framell->SetMaximum(3) ;
   // sprintf(plotname,"%s/LogLikelihood_nSigm.png",CPlot::sOutDir.Data());
   // c2->Clear();framesigmll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nEWKm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameewkmll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_c1QCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmc1ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_c2QCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmc2ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_c3QCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmc3ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_d1QCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmd1ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_d2QCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmd2ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_d3QCDm.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmd3ll->Draw();c2->Update();c2->SaveAs(plotname);
   
   // sprintf(plotname,"%s/LogLikelihood_nQCDm1.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmll1->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDm2.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmll2->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDm3.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmll3->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDm4.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmll4->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDm5.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdmll5->Draw();c2->Update();c2->SaveAs(plotname);
   
       // // // // Plot likelihood scan frac 
   // RooPlot* framesigpll = nSigp_[0]->frame(Bins(10),Range(0.9*nSigp_[0]->getVal(),1.1*nSigp_[0]->getVal()),Title("")) ; 
   // RooPlot* frameewkpll = nEWKp_[0]->frame(Bins(10),Range(0.9*nEWKp_[0]->getVal(),1.1*nEWKp_[0]->getVal()),Title("")) ; 
   // RooPlot* frameqcdpll = nQCDp_[0]->frame(Bins(10),Range(0.9*nQCDp_[0]->getVal(),1.1*nQCDp_[0]->getVal()),Title("")) ; 
   // RooPlot* frameqcdpc1ll = (qcdp_[0]->c1)->frame(Bins(10),Title("")) ; 
   // RooPlot* frameqcdpc2ll = (qcdp_[0]->c2)->frame(Bins(10),Title("")) ; 
   // RooPlot* frameqcdpc3ll = (qcdp_[0]->c3)->frame(Bins(10),Title("")) ; 
   // RooPlot* frameqcdpd1ll = (qcdp_[0]->d1)->frame(Bins(10),Title("")) ; 
   // RooPlot* frameqcdpd2ll = (qcdp_[0]->d2)->frame(Bins(10),Title("")) ; 
   // RooPlot* frameqcdpd3ll = (qcdp_[0]->d3)->frame(Bins(10),Title("")) ; 

   // RooPlot* frameqcdpll1 = nQCDp_[1]->frame(Bins(10),Range(0.9*nQCDp_[1]->getVal(),1.1*nQCDp_[1]->getVal()),Title("")) ; 
   // RooPlot* frameqcdpll2 = nQCDp_[2]->frame(Bins(10),Range(0.9*nQCDp_[2]->getVal(),1.1*nQCDp_[2]->getVal()),Title("")) ; 
   // RooPlot* frameqcdpll3 = nQCDp_[3]->frame(Bins(10),Range(0.9*nQCDp_[3]->getVal(),1.1*nQCDp_[3]->getVal()),Title("")) ; 
   // RooPlot* frameqcdpll4 = nQCDp_[4]->frame(Bins(10),Range(0.9*nQCDp_[4]->getVal(),1.1*nQCDp_[4]->getVal()),Title("")) ; 
   // RooPlot* frameqcdpll5 = nQCDp_[5]->frame(Bins(10),Range(0.9*nQCDp_[5]->getVal(),1.1*nQCDp_[5]->getVal()),Title("")) ; 
	
    // RooAbsReal* pll_nsigp = nllp->createProfile(*nSigp_[0]) ;
    // RooAbsReal* pll_newkp = nllp->createProfile(*nEWKp_[0]) ;
    // RooAbsReal* pll_nqcdp = nllp->createProfile(*nQCDp_[0]) ;
    // RooAbsReal* pll_qcdpc1 = nllp->createProfile(*(qcdp_[0]->c1)) ;
    // RooAbsReal* pll_qcdpc2 = nllp->createProfile(*(qcdp_[0]->c2)) ;
    // RooAbsReal* pll_qcdpc3 = nllp->createProfile(*(qcdp_[0]->c3)) ;
    // RooAbsReal* pll_qcdpd1 = nllp->createProfile(*(qcdp_[0]->d1)) ;
    // RooAbsReal* pll_qcdpd2 = nllp->createProfile(*(qcdp_[0]->d2)) ;
    // RooAbsReal* pll_qcdpd3 = nllp->createProfile(*(qcdp_[0]->d3)) ;

	// RooAbsReal* pll_nqcdp1 = nllp->createProfile(*nQCDp_[1]) ;
	// RooAbsReal* pll_nqcdp2 = nllp->createProfile(*nQCDp_[2]) ;
	// RooAbsReal* pll_nqcdp3 = nllp->createProfile(*nQCDp_[3]) ;
	// RooAbsReal* pll_nqcdp4 = nllp->createProfile(*nQCDp_[4]) ;
	// RooAbsReal* pll_nqcdp5 = nllp->createProfile(*nQCDp_[5]) ;
	
    // // // // // Plot the profile likelihood in frac


    // // // // Plot the profile likelihood in frac
    // pll_nsigp->plotOn(framesigpll,LineColor(kRed)) ; framesigpll->SetMinimum(0); framesigpll->SetMaximum(10);
    // pll_newkp->plotOn(frameewkpll,LineColor(kRed)) ; frameewkpll->SetMinimum(0); frameewkpll->SetMaximum(10);
    // pll_nqcdp->plotOn(frameqcdpll,LineColor(kRed)) ; frameqcdpll->SetMinimum(0); frameqcdpll->SetMaximum(10);
    // pll_qcdpc1->plotOn(frameqcdpc1ll,LineColor(kRed)) ; frameqcdpc1ll->SetMinimum(0); frameqcdpc1ll->SetMaximum(10);
    // pll_qcdpc2->plotOn(frameqcdpc2ll,LineColor(kRed)) ; frameqcdpc2ll->SetMinimum(0); frameqcdpc2ll->SetMaximum(10);
    // pll_qcdpc3->plotOn(frameqcdpc3ll,LineColor(kRed)) ; frameqcdpc3ll->SetMinimum(0); frameqcdpc3ll->SetMaximum(10);
    // pll_qcdpd1->plotOn(frameqcdpd1ll,LineColor(kRed)) ; frameqcdpd1ll->SetMinimum(0); frameqcdpd1ll->SetMaximum(10);
    // pll_qcdpd2->plotOn(frameqcdpd2ll,LineColor(kRed)) ; frameqcdpd2ll->SetMinimum(0); frameqcdpd2ll->SetMaximum(10);
    // pll_qcdpd3->plotOn(frameqcdpd3ll,LineColor(kRed)) ; frameqcdpd3ll->SetMinimum(0); frameqcdpd3ll->SetMaximum(10);
	
	// pll_nqcdp1->plotOn(frameqcdpll1,LineColor(kRed)) ; frameqcdpll1->SetMinimum(0); frameqcdpll1->SetMaximum(10);
	// pll_nqcdp2->plotOn(frameqcdpll2,LineColor(kRed)) ; frameqcdpll2->SetMinimum(0); frameqcdpll2->SetMaximum(10);
	// pll_nqcdp3->plotOn(frameqcdpll3,LineColor(kRed)) ; frameqcdpll3->SetMinimum(0); frameqcdpll3->SetMaximum(10);
	// pll_nqcdp4->plotOn(frameqcdpll4,LineColor(kRed)) ; frameqcdpll4->SetMinimum(0); frameqcdpll4->SetMaximum(10);
	// pll_nqcdp5->plotOn(frameqcdpll5,LineColor(kRed)) ; frameqcdpll5->SetMinimum(0); frameqcdpll5->SetMaximum(10);
	

   
   
    // sprintf(plotname,"%s/LogLikelihood_nSigp.png",CPlot::sOutDir.Data());
   // c2->Clear();framesigpll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nEWKp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameewkpll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_c1QCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpc1ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_c2QCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpc2ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_c3QCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpc3ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_d1QCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpd1ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_d2QCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpd2ll->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_d3QCDp.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpd3ll->Draw();c2->Update();c2->SaveAs(plotname);
   
   // sprintf(plotname,"%s/LogLikelihood_nQCDp1.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpll1->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDp2.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpll2->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDp3.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpll3->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDp4.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpll4->Draw();c2->Update();c2->SaveAs(plotname);
   // sprintf(plotname,"%s/LogLikelihood_nQCDp5.png",CPlot::sOutDir.Data());
   // c2->Clear();frameqcdpll5->Draw();c2->Update();c2->SaveAs(plotname);
  
  
  
  
  
  

  
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
    std::cout << "bin " << ibin << std::endl;
    std::cout << "data " << hData->GetBinContent(ibin) << std::endl;
    std::cout << "fits " << hFit->GetBinContent(ibin) << std::endl;
//     Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    if(hData->GetBinContent(ibin) == 0) err = 0;
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    std::cout << "err = " << err << std::endl;
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
  }
  std::cout << "bin " << hData->GetNbinsX() + 1 << std::endl;
  std::cout << "data " << hData->GetBinContent(hData->GetNbinsX() + 1) << std::endl;
  std::cout << "fits " << hFit->GetBinContent(hData->GetNbinsX() + 1) << std::endl;
  
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
