//================================================================================================
// 
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>   
#include <TGraphAsymmErrors.h>   
#include <TColor.h>             // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing


#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1);
TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);
TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

//=== MAIN MACRO ================================================================================================= 

void plotZmmCorrelations(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmmCorrelations");
  gStyle->SetTitleOffset(0.75,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  //
  // input ntuple file names
  //
vector<TFile*> file;
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStar.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRap.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPt.root", "OPEN"));

  vector<TFile*> fileLumiUp;
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtLumiUp.root", "OPEN"));

  vector<TFile*> fileLumiDown;
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtLumiDown.root", "OPEN"));

  vector<TFile*> fileEWKBkgUp;
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEWKBkgUp.root", "OPEN"));

  vector<TFile*> fileEWKBkgDown;
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEWKBkgDown.root", "OPEN"));

  vector<TFile*> fileTopBkgUp;
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtTopBkgUp.root", "OPEN"));

  vector<TFile*> fileTopBkgDown;
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtTopBkgDown.root", "OPEN"));

  vector<TFile*> fileEffBinSys;
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEffBin.root", "OPEN"));

  vector<TFile*> fileEffStatUp;
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEffStatUp.root", "OPEN"));

  vector<TFile*> fileEffStatDown;
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEffStatDown.root", "OPEN"));

  vector<TFile*> fileEffSigShapeSys;
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEffSigShape.root", "OPEN"));

  vector<TFile*> fileEffBkgShapeSys;
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtEffBkgShape.root", "OPEN"));

  vector<TFile*> fileResScaleSys;
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtResScale.root", "OPEN"));

  vector<TFile*> fileUnfoldModelSys;
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtUnfoldModel_Smoothed.root", "OPEN"));

  vector<TFile*> fileUnfoldMatrixSys;
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZPtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputPhiStarUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputZRapUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1PtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2PtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep1EtaUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLep2EtaUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepNegPtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmumu/UnfoldingOutputLepPosPtUnfoldMatrix.root", "OPEN"));

  
  // plot output file format
  const TString format("png");

   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  
   
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_Systematics.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");

  double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
  double PhiStarBins[]={0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.01,0.012,0.014,0.016,0.018,0.021,0.024,0.027,0.030,0.034,0.038,0.044,0.050,0.058,0.066,0.076,0.088,0.10,0.12,0.14,0.16,0.18,0.20,0.24,0.28,0.34,0.42,0.52,0.64,0.8,1.0,1.5,2,3};
  double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double Lep2PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150};
  double LepNegPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};
  double LepPosPtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};

  const int nBinsZPt= sizeof(ZPtBins)/sizeof(double)-1;
  const int nBinsPhiStar= sizeof(PhiStarBins)/sizeof(double)-1;
  const int nBinsLep1Pt= sizeof(Lep1PtBins)/sizeof(double)-1;
  const int nBinsLep2Pt= sizeof(Lep2PtBins)/sizeof(double)-1;
  const int nBinsLepNegPt= sizeof(LepNegPtBins)/sizeof(double)-1;
  const int nBinsLepPosPt= sizeof(LepPosPtBins)/sizeof(double)-1;
  
  // histograms
  TH2D *ZPT_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[0]->Get("hCov_ResScale"));
  TH2D *PHISTAR_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[1]->Get("hCov_ResScale"));
  TH2D *ZRAP_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[2]->Get("hCov_ResScale"));
  TH2D *LEP1PT_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[3]->Get("hCov_ResScale"));
  TH2D *LEP2PT_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[4]->Get("hCov_ResScale"));
  TH2D *LEP1ETA_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[5]->Get("hCov_ResScale"));
  TH2D *LEP2ETA_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[6]->Get("hCov_ResScale"));
  TH2D *LEPNEGPT_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[7]->Get("hCov_ResScale"));
  TH2D *LEPPOSPT_COV_MATRIX_RESSCALE=(TH2D*)(fileResScaleSys[8]->Get("hCov_ResScale"));

  TH2D *ZPT_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[0]->Get("hCov_MatrixStat"));
  TH2D *PHISTAR_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[1]->Get("hCov_MatrixStat"));
  TH2D *ZRAP_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[2]->Get("hCov_MatrixStat"));
  TH2D *LEP1PT_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[3]->Get("hCov_MatrixStat"));
  TH2D *LEP2PT_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[4]->Get("hCov_MatrixStat"));
  TH2D *LEP1ETA_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[5]->Get("hCov_MatrixStat"));
  TH2D *LEP2ETA_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[6]->Get("hCov_MatrixStat"));
  TH2D *LEPNEGPT_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[7]->Get("hCov_MatrixStat"));
  TH2D *LEPPOSPT_COV_MATRIX_MATRIXSTAT=(TH2D*)(fileUnfoldMatrixSys[8]->Get("hCov_MatrixStat"));

  TH2D *ZPT_COV_MATRIX=new TH2D((string("ZPT_COV_MATRIX")).c_str(),(string("ZPT_COV_MATRIX")).c_str(),nBinsZPt,ZPtBins,nBinsZPt,ZPtBins);
  TH2D *PHISTAR_COV_MATRIX=new TH2D((string("PHISTAR_COV_MATRIX")).c_str(),(string("PHISTAR_COV_MATRIX")).c_str(),nBinsPhiStar,PhiStarBins,nBinsPhiStar,PhiStarBins);
  TH2D *ZRAP_COV_MATRIX=new TH2D((string("ZRAP_COV_MATRIX")).c_str(),(string("ZRAP_COV_MATRIX")).c_str(),24,0,2.4,24,0,2.4);
  TH2D *LEP1PT_COV_MATRIX=new TH2D((string("LEP1PT_COV_MATRIX")).c_str(),(string("LEP1PT_COV_MATRIX")).c_str(),nBinsLep1Pt,Lep1PtBins,nBinsLep1Pt,Lep1PtBins);
  TH2D *LEP2PT_COV_MATRIX=new TH2D((string("LEP2PT_COV_MATRIX")).c_str(),(string("LEP2PT_COV_MATRIX")).c_str(),nBinsLep2Pt,Lep2PtBins,nBinsLep2Pt,Lep2PtBins);
  TH2D *LEP1ETA_COV_MATRIX=new TH2D((string("LEP1ETA_COV_MATRIX")).c_str(),(string("LEP1ETA_COV_MATRIX")).c_str(),24,0,2.4,24,0,2.4);
  TH2D *LEP2ETA_COV_MATRIX=new TH2D((string("LEP2ETA_COV_MATRIX")).c_str(),(string("LEP2ETA_COV_MATRIX")).c_str(),24,0,2.4,24,0,2.4);
  TH2D *LEPNEGPT_COV_MATRIX=new TH2D((string("LEPNEGPT_COV_MATRIX")).c_str(),(string("LEPNEGPT_COV_MATRIX")).c_str(),nBinsLepNegPt,LepNegPtBins,nBinsLepNegPt,LepNegPtBins);
  TH2D *LEPPOSPT_COV_MATRIX=new TH2D((string("LEPPOSPT_COV_MATRIX")).c_str(),(string("LEPPOSPT_COV_MATRIX")).c_str(),nBinsLepPosPt,LepPosPtBins,nBinsLepPosPt,LepPosPtBins);

  TH2D *ZPT_CORR_MATRIX=new TH2D((string("ZPT_CORR_MATRIX")).c_str(),(string("ZPT_CORR_MATRIX")).c_str(),nBinsZPt,ZPtBins,nBinsZPt,ZPtBins);
  TH2D *PHISTAR_CORR_MATRIX=new TH2D((string("PHISTAR_CORR_MATRIX")).c_str(),(string("PHISTAR_CORR_MATRIX")).c_str(),nBinsPhiStar,PhiStarBins,nBinsPhiStar,PhiStarBins);
  TH2D *ZRAP_CORR_MATRIX=new TH2D((string("ZRAP_CORR_MATRIX")).c_str(),(string("ZRAP_CORR_MATRIX")).c_str(),24,0,2.4,24,0,2.4);
  TH2D *LEP1PT_CORR_MATRIX=new TH2D((string("LEP1PT_CORR_MATRIX")).c_str(),(string("LEP1PT_CORR_MATRIX")).c_str(),nBinsLep1Pt,Lep1PtBins,nBinsLep1Pt,Lep1PtBins);
  TH2D *LEP2PT_CORR_MATRIX=new TH2D((string("LEP2PT_CORR_MATRIX")).c_str(),(string("LEP2PT_CORR_MATRIX")).c_str(),nBinsLep2Pt,Lep2PtBins,nBinsLep2Pt,Lep2PtBins);
  TH2D *LEP1ETA_CORR_MATRIX=new TH2D((string("LEP1ETA_CORR_MATRIX")).c_str(),(string("LEP1ETA_CORR_MATRIX")).c_str(),24,0,2.4,24,0,2.4);
  TH2D *LEP2ETA_CORR_MATRIX=new TH2D((string("LEP2ETA_CORR_MATRIX")).c_str(),(string("LEP2ETA_CORR_MATRIX")).c_str(),24,0,2.4,24,0,2.4);
  TH2D *LEPNEGPT_CORR_MATRIX=new TH2D((string("LEPNEGPT_CORR_MATRIX")).c_str(),(string("LEPNEGPT_CORR_MATRIX")).c_str(),nBinsLepNegPt,LepNegPtBins,nBinsLepNegPt,LepNegPtBins);
  TH2D *LEPPOSPT_CORR_MATRIX=new TH2D((string("LEPPOSPT_CORR_MATRIX")).c_str(),(string("LEPPOSPT_CORR_MATRIX")).c_str(),nBinsLepPosPt,LepPosPtBins,nBinsLepPosPt,LepPosPtBins);


  //--------------------------------------------------------------------------
  //                           Z Pt
  //--------------------------------------------------------------------------



  TH1D * hUnfoldZPt=(TH1D*)(file[0]->Get("hUnfold"));
  
  TH1D * hUnfoldZPtEWKBkgUp;
  TH1D * hUnfoldZPtEWKBkgDown;
  hUnfoldZPtEWKBkgUp=(TH1D*)(fileEWKBkgUp[0]->Get("hUnfold"));
  hUnfoldZPtEWKBkgDown=(TH1D*)(fileEWKBkgDown[0]->Get("hUnfold"));

  TH1D * hUnfoldZPtTopBkgUp;
  TH1D * hUnfoldZPtTopBkgDown;
  hUnfoldZPtTopBkgUp=(TH1D*)(fileTopBkgUp[0]->Get("hUnfold"));
  hUnfoldZPtTopBkgDown=(TH1D*)(fileTopBkgDown[0]->Get("hUnfold"));

  TH1D * hUnfoldZPtEffStatUp;
  TH1D * hUnfoldZPtEffStatDown;
  hUnfoldZPtEffStatUp=(TH1D*)(fileEffStatUp[0]->Get("hUnfold"));
  hUnfoldZPtEffStatDown=(TH1D*)(fileEffStatDown[0]->Get("hUnfold"));

  TH1D * hUnfoldZPtEffBinSysUp;
  TH1D * hUnfoldZPtEffBinSysDown;
  hUnfoldZPtEffBinSysUp=(TH1D*)(fileEffBinSys[0]->Get("hUnfold"));
  hUnfoldZPtEffBinSysDown=(TH1D*)(fileEffBinSys[0]->Get("hUnfold"));

  hUnfoldZPtEffBinSysDown->Add(hUnfoldZPt,-1.);
  hUnfoldZPtEffBinSysDown->Scale(-1.);
  hUnfoldZPtEffBinSysDown->Add(hUnfoldZPt,1.);

  TH1D * hUnfoldZPtEffSigShapeSysUp;
  TH1D * hUnfoldZPtEffSigShapeSysDown;
  hUnfoldZPtEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[0]->Get("hUnfold"));
  hUnfoldZPtEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[0]->Get("hUnfold"));

  hUnfoldZPtEffSigShapeSysDown->Add(hUnfoldZPt,-1.);
  hUnfoldZPtEffSigShapeSysDown->Scale(-1.);
  hUnfoldZPtEffSigShapeSysDown->Add(hUnfoldZPt,1.);

  TH1D * hUnfoldZPtEffBkgShapeSysUp;
  TH1D * hUnfoldZPtEffBkgShapeSysDown;
  hUnfoldZPtEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[0]->Get("hUnfold"));
  hUnfoldZPtEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[0]->Get("hUnfold"));

  hUnfoldZPtEffBkgShapeSysDown->Add(hUnfoldZPt,-1.);
  hUnfoldZPtEffBkgShapeSysDown->Scale(-1.);
  hUnfoldZPtEffBkgShapeSysDown->Add(hUnfoldZPt,1.);

  TH1D * hUnfoldZPtUnfoldModelSysUp;
  TH1D * hUnfoldZPtUnfoldModelSysDown;
  hUnfoldZPtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[0]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldZPtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[0]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldZPtUnfoldModelSysDown->Add(hUnfoldZPt,-1.);
  hUnfoldZPtUnfoldModelSysDown->Scale(-1.);
  hUnfoldZPtUnfoldModelSysDown->Add(hUnfoldZPt,1.);

  
  TGraphAsymmErrors* gUnfoldZPt=TH1TOTGraphAsymmErrors(hUnfoldZPt);

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPt);
  
  TGraphAsymmErrors* ZPT_EWKBKG_UNCERT_BAND_DATA;
  ZPT_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEWKBkgDown));

  TGraphAsymmErrors* ZPT_TOPBKG_UNCERT_BAND_DATA;
  ZPT_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZPtTopBkgDown));

  TGraphAsymmErrors* ZPT_EFFSTAT_UNCERT_BAND_DATA;
  ZPT_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEffStatDown));
  
  TGraphAsymmErrors* ZPT_EFFBIN_UNCERT_BAND_DATA;
  ZPT_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEffBinSysDown));
  
  TGraphAsymmErrors* ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA;
  ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEffSigShapeSysDown));

  TGraphAsymmErrors* ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA;
  ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEffBkgShapeSysDown));

  TGraphAsymmErrors* ZPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  ZPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_ZPT=ZPT_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_ZPT=ZPT_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_ZPT=ZPT_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_ZPT=ZPT_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_ZPT=ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_ZPT=ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_ZPT=ZPT_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_ZPt[7][nBinsZPt];
  double cov_ZPt[nBinsZPt][nBinsZPt];
  
  for(int i=0;i!=nBinsZPt;++i)
    {
      if(hUnfoldZPtEWKBkgUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[0][i]=EWKBKG_SYSTEMATIC_UNCERT_ZPT[i];
	}
      if(hUnfoldZPtTopBkgUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[1][i]=TOPBKG_SYSTEMATIC_UNCERT_ZPT[i];
	}
      if(hUnfoldZPtEffStatUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_ZPT[i];
	}
      if(hUnfoldZPtEffBinSysUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[3][i]=EFFBIN_SYSTEMATIC_UNCERT_ZPT[i];
	}
      if(hUnfoldZPtEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_ZPT[i];
	}
      if(hUnfoldZPtEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_ZPT[i];
	}
      if(hUnfoldZPtUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldZPt->GetBinContent(i+1)<1)
	{
	  uncert_ZPt[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_ZPT[i];
	}
      else
	{
	  uncert_ZPt[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_ZPT[i];
	}
    }

  for(int i=0;i!=nBinsZPt;++i)
    {
      for(int j=0;j!=nBinsZPt;++j)
	{
	  cov_ZPt[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_ZPt[i][j]+=uncert_ZPt[k][i]*uncert_ZPt[k][j];
	    }
	  ZPT_COV_MATRIX->SetBinContent(i+1,j+1,cov_ZPt[i][j]);
	  ZPT_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  ZPT_COV_MATRIX->Add(ZPT_COV_MATRIX_RESSCALE,1);
  ZPT_COV_MATRIX->Add(ZPT_COV_MATRIX_MATRIXSTAT,1);
  double corr_ZPt[nBinsZPt][nBinsZPt];
  for(int i=0;i!=nBinsZPt;++i)
    {
      for(int j=0;j!=nBinsZPt;++j)
	{
	  corr_ZPt[i][j]=ZPT_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(ZPT_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(ZPT_COV_MATRIX->GetBinContent(j+1,j+1)));
	  ZPT_CORR_MATRIX->SetBinContent(i+1,j+1,corr_ZPt[i][j]);
	  ZPT_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  //--------------------------------------------------------------------------
  //                           PhiStar
  //--------------------------------------------------------------------------

  TH1D * hUnfoldPhiStar=(TH1D*)(file[1]->Get("hUnfold"));
  
  TH1D * hUnfoldPhiStarEWKBkgUp;
  TH1D * hUnfoldPhiStarEWKBkgDown;
  hUnfoldPhiStarEWKBkgUp=(TH1D*)(fileEWKBkgUp[1]->Get("hUnfold"));
  hUnfoldPhiStarEWKBkgDown=(TH1D*)(fileEWKBkgDown[1]->Get("hUnfold"));

  TH1D * hUnfoldPhiStarTopBkgUp;
  TH1D * hUnfoldPhiStarTopBkgDown;
  hUnfoldPhiStarTopBkgUp=(TH1D*)(fileTopBkgUp[1]->Get("hUnfold"));
  hUnfoldPhiStarTopBkgDown=(TH1D*)(fileTopBkgDown[1]->Get("hUnfold"));

  TH1D * hUnfoldPhiStarEffStatUp;
  TH1D * hUnfoldPhiStarEffStatDown;
  hUnfoldPhiStarEffStatUp=(TH1D*)(fileEffStatUp[1]->Get("hUnfold"));
  hUnfoldPhiStarEffStatDown=(TH1D*)(fileEffStatDown[1]->Get("hUnfold"));

  TH1D * hUnfoldPhiStarEffBinSysUp;
  TH1D * hUnfoldPhiStarEffBinSysDown;
  hUnfoldPhiStarEffBinSysUp=(TH1D*)(fileEffBinSys[1]->Get("hUnfold"));
  hUnfoldPhiStarEffBinSysDown=(TH1D*)(fileEffBinSys[1]->Get("hUnfold"));

  hUnfoldPhiStarEffBinSysDown->Add(hUnfoldPhiStar,-1.);
  hUnfoldPhiStarEffBinSysDown->Scale(-1.);
  hUnfoldPhiStarEffBinSysDown->Add(hUnfoldPhiStar,1.);

  TH1D * hUnfoldPhiStarEffSigShapeSysUp;
  TH1D * hUnfoldPhiStarEffSigShapeSysDown;
  hUnfoldPhiStarEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[1]->Get("hUnfold"));
  hUnfoldPhiStarEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[1]->Get("hUnfold"));

  hUnfoldPhiStarEffSigShapeSysDown->Add(hUnfoldPhiStar,-1.);
  hUnfoldPhiStarEffSigShapeSysDown->Scale(-1.);
  hUnfoldPhiStarEffSigShapeSysDown->Add(hUnfoldPhiStar,1.);

  TH1D * hUnfoldPhiStarEffBkgShapeSysUp;
  TH1D * hUnfoldPhiStarEffBkgShapeSysDown;
  hUnfoldPhiStarEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[1]->Get("hUnfold"));
  hUnfoldPhiStarEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[1]->Get("hUnfold"));

  hUnfoldPhiStarEffBkgShapeSysDown->Add(hUnfoldPhiStar,-1.);
  hUnfoldPhiStarEffBkgShapeSysDown->Scale(-1.);
  hUnfoldPhiStarEffBkgShapeSysDown->Add(hUnfoldPhiStar,1.);

  TH1D * hUnfoldPhiStarUnfoldModelSysUp;
  TH1D * hUnfoldPhiStarUnfoldModelSysDown;
  hUnfoldPhiStarUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[1]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldPhiStarUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[1]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldPhiStarUnfoldModelSysDown->Add(hUnfoldPhiStar,-1.);
  hUnfoldPhiStarUnfoldModelSysDown->Scale(-1.);
  hUnfoldPhiStarUnfoldModelSysDown->Add(hUnfoldPhiStar,1.);

  
  TGraphAsymmErrors* gUnfoldPhiStar=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);
  
  TGraphAsymmErrors* PHISTAR_EWKBKG_UNCERT_BAND_DATA;
  PHISTAR_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEWKBkgDown));

  TGraphAsymmErrors* PHISTAR_TOPBKG_UNCERT_BAND_DATA;
  PHISTAR_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarTopBkgDown));

  TGraphAsymmErrors* PHISTAR_EFFSTAT_UNCERT_BAND_DATA;
  PHISTAR_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffStatDown));
  
  TGraphAsymmErrors* PHISTAR_EFFBIN_UNCERT_BAND_DATA;
  PHISTAR_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffBinSysDown));
  
  TGraphAsymmErrors* PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA;
  PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffSigShapeSysDown));

  TGraphAsymmErrors* PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA;
  PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffBkgShapeSysDown));

  TGraphAsymmErrors* PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA;
  PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_PhiStar[7][nBinsPhiStar];
  double cov_PhiStar[nBinsPhiStar][nBinsPhiStar];
  
  for(int i=0;i!=nBinsPhiStar;++i)
    {
      if(hUnfoldPhiStarEWKBkgUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[0][i]=EWKBKG_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      if(hUnfoldPhiStarTopBkgUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[1][i]=TOPBKG_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      if(hUnfoldPhiStarEffStatUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      if(hUnfoldPhiStarEffBinSysUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[3][i]=EFFBIN_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      if(hUnfoldPhiStarEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      if(hUnfoldPhiStarEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      if(hUnfoldPhiStarUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldPhiStar->GetBinContent(i+1)<1)
	{
	  uncert_PhiStar[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
      else
	{
	  uncert_PhiStar[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_PHISTAR[i];
	}
    }

  for(int i=0;i!=nBinsPhiStar;++i)
    {
      for(int j=0;j!=nBinsPhiStar;++j)
	{
	  cov_PhiStar[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_PhiStar[i][j]+=uncert_PhiStar[k][i]*uncert_PhiStar[k][j];
	    }
	  PHISTAR_COV_MATRIX->SetBinContent(i+1,j+1,cov_PhiStar[i][j]);
	  PHISTAR_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  PHISTAR_COV_MATRIX->Add(PHISTAR_COV_MATRIX_RESSCALE,1);
  PHISTAR_COV_MATRIX->Add(PHISTAR_COV_MATRIX_MATRIXSTAT,1);
  double corr_PhiStar[nBinsPhiStar][nBinsPhiStar];
  for(int i=0;i!=nBinsPhiStar;++i)
    {
      for(int j=0;j!=nBinsPhiStar;++j)
	{
	  corr_PhiStar[i][j]=PHISTAR_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(PHISTAR_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(PHISTAR_COV_MATRIX->GetBinContent(j+1,j+1)));
	  PHISTAR_CORR_MATRIX->SetBinContent(i+1,j+1,corr_PhiStar[i][j]);
	  PHISTAR_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }

  


  //--------------------------------------------------------------------------
  //                           Z Rapidity
  //--------------------------------------------------------------------------

  TH1D * hUnfoldZRap=(TH1D*)(file[2]->Get("hUnfold"));
  
  TH1D * hUnfoldZRapEWKBkgUp;
  TH1D * hUnfoldZRapEWKBkgDown;
  hUnfoldZRapEWKBkgUp=(TH1D*)(fileEWKBkgUp[2]->Get("hUnfold"));
  hUnfoldZRapEWKBkgDown=(TH1D*)(fileEWKBkgDown[2]->Get("hUnfold"));

  TH1D * hUnfoldZRapTopBkgUp;
  TH1D * hUnfoldZRapTopBkgDown;
  hUnfoldZRapTopBkgUp=(TH1D*)(fileTopBkgUp[2]->Get("hUnfold"));
  hUnfoldZRapTopBkgDown=(TH1D*)(fileTopBkgDown[2]->Get("hUnfold"));

  TH1D * hUnfoldZRapEffStatUp;
  TH1D * hUnfoldZRapEffStatDown;
  hUnfoldZRapEffStatUp=(TH1D*)(fileEffStatUp[2]->Get("hUnfold"));
  hUnfoldZRapEffStatDown=(TH1D*)(fileEffStatDown[2]->Get("hUnfold"));

  TH1D * hUnfoldZRapEffBinSysUp;
  TH1D * hUnfoldZRapEffBinSysDown;
  hUnfoldZRapEffBinSysUp=(TH1D*)(fileEffBinSys[2]->Get("hUnfold"));
  hUnfoldZRapEffBinSysDown=(TH1D*)(fileEffBinSys[2]->Get("hUnfold"));

  hUnfoldZRapEffBinSysDown->Add(hUnfoldZRap,-1.);
  hUnfoldZRapEffBinSysDown->Scale(-1.);
  hUnfoldZRapEffBinSysDown->Add(hUnfoldZRap,1.);

  TH1D * hUnfoldZRapEffSigShapeSysUp;
  TH1D * hUnfoldZRapEffSigShapeSysDown;
  hUnfoldZRapEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[2]->Get("hUnfold"));
  hUnfoldZRapEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[2]->Get("hUnfold"));

  hUnfoldZRapEffSigShapeSysDown->Add(hUnfoldZRap,-1.);
  hUnfoldZRapEffSigShapeSysDown->Scale(-1.);
  hUnfoldZRapEffSigShapeSysDown->Add(hUnfoldZRap,1.);

  TH1D * hUnfoldZRapEffBkgShapeSysUp;
  TH1D * hUnfoldZRapEffBkgShapeSysDown;
  hUnfoldZRapEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[2]->Get("hUnfold"));
  hUnfoldZRapEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[2]->Get("hUnfold"));

  hUnfoldZRapEffBkgShapeSysDown->Add(hUnfoldZRap,-1.);
  hUnfoldZRapEffBkgShapeSysDown->Scale(-1.);
  hUnfoldZRapEffBkgShapeSysDown->Add(hUnfoldZRap,1.);

  TH1D * hUnfoldZRapUnfoldModelSysUp;
  TH1D * hUnfoldZRapUnfoldModelSysDown;
  hUnfoldZRapUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[2]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldZRapUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[2]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldZRapUnfoldModelSysDown->Add(hUnfoldZRap,-1.);
  hUnfoldZRapUnfoldModelSysDown->Scale(-1.);
  hUnfoldZRapUnfoldModelSysDown->Add(hUnfoldZRap,1.);

  
  TGraphAsymmErrors* gUnfoldZRap=TH1TOTGraphAsymmErrors(hUnfoldZRap);

  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRap);
  
  TGraphAsymmErrors* ZRAP_EWKBKG_UNCERT_BAND_DATA;
  ZRAP_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEWKBkgDown));

  TGraphAsymmErrors* ZRAP_TOPBKG_UNCERT_BAND_DATA;
  ZRAP_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZRapTopBkgDown));

  TGraphAsymmErrors* ZRAP_EFFSTAT_UNCERT_BAND_DATA;
  ZRAP_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEffStatDown));
  
  TGraphAsymmErrors* ZRAP_EFFBIN_UNCERT_BAND_DATA;
  ZRAP_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEffBinSysDown));
  
  TGraphAsymmErrors* ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA;
  ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEffSigShapeSysDown));

  TGraphAsymmErrors* ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA;
  ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEffBkgShapeSysDown));

  TGraphAsymmErrors* ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA;
  ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_ZRAP=ZRAP_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_ZRAP=ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_ZRap[7][24];
  double cov_ZRap[24][24];
  
  for(int i=0;i!=24;++i)
    {
      if(hUnfoldZRapEWKBkgUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[0][i]=EWKBKG_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      if(hUnfoldZRapTopBkgUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[1][i]=TOPBKG_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      if(hUnfoldZRapEffStatUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      if(hUnfoldZRapEffBinSysUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[3][i]=EFFBIN_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      if(hUnfoldZRapEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      if(hUnfoldZRapEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      if(hUnfoldZRapUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldZRap->GetBinContent(i+1)<1)
	{
	  uncert_ZRap[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_ZRAP[i];
	}
      else
	{
	  uncert_ZRap[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_ZRAP[i];
	}
    }

  for(int i=0;i!=24;++i)
    {
      for(int j=0;j!=24;++j)
	{
	  cov_ZRap[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_ZRap[i][j]+=uncert_ZRap[k][i]*uncert_ZRap[k][j];
	    }
	  ZRAP_COV_MATRIX->SetBinContent(i+1,j+1,cov_ZRap[i][j]);
	  ZRAP_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  ZRAP_COV_MATRIX->Add(ZRAP_COV_MATRIX_RESSCALE,1);
  ZRAP_COV_MATRIX->Add(ZRAP_COV_MATRIX_MATRIXSTAT,1);
  double corr_ZRap[24][24];
  for(int i=0;i!=24;++i)
    {
      for(int j=0;j!=24;++j)
	{
	  corr_ZRap[i][j]=ZRAP_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(ZRAP_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(ZRAP_COV_MATRIX->GetBinContent(j+1,j+1)));
	  ZRAP_CORR_MATRIX->SetBinContent(i+1,j+1,corr_ZRap[i][j]);
	  ZRAP_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
 
  
  //--------------------------------------------------------------------------
  //                           Lep1 Pt
  //--------------------------------------------------------------------------

  TH1D * hUnfoldLep1Pt=(TH1D*)(file[3]->Get("hUnfold"));
  
  TH1D * hUnfoldLep1PtEWKBkgUp;
  TH1D * hUnfoldLep1PtEWKBkgDown;
  hUnfoldLep1PtEWKBkgUp=(TH1D*)(fileEWKBkgUp[3]->Get("hUnfold"));
  hUnfoldLep1PtEWKBkgDown=(TH1D*)(fileEWKBkgDown[3]->Get("hUnfold"));

  TH1D * hUnfoldLep1PtTopBkgUp;
  TH1D * hUnfoldLep1PtTopBkgDown;
  hUnfoldLep1PtTopBkgUp=(TH1D*)(fileTopBkgUp[3]->Get("hUnfold"));
  hUnfoldLep1PtTopBkgDown=(TH1D*)(fileTopBkgDown[3]->Get("hUnfold"));

  TH1D * hUnfoldLep1PtEffStatUp;
  TH1D * hUnfoldLep1PtEffStatDown;
  hUnfoldLep1PtEffStatUp=(TH1D*)(fileEffStatUp[3]->Get("hUnfold"));
  hUnfoldLep1PtEffStatDown=(TH1D*)(fileEffStatDown[3]->Get("hUnfold"));

  TH1D * hUnfoldLep1PtEffBinSysUp;
  TH1D * hUnfoldLep1PtEffBinSysDown;
  hUnfoldLep1PtEffBinSysUp=(TH1D*)(fileEffBinSys[3]->Get("hUnfold"));
  hUnfoldLep1PtEffBinSysDown=(TH1D*)(fileEffBinSys[3]->Get("hUnfold"));

  hUnfoldLep1PtEffBinSysDown->Add(hUnfoldLep1Pt,-1.);
  hUnfoldLep1PtEffBinSysDown->Scale(-1.);
  hUnfoldLep1PtEffBinSysDown->Add(hUnfoldLep1Pt,1.);

  TH1D * hUnfoldLep1PtEffSigShapeSysUp;
  TH1D * hUnfoldLep1PtEffSigShapeSysDown;
  hUnfoldLep1PtEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[3]->Get("hUnfold"));
  hUnfoldLep1PtEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[3]->Get("hUnfold"));

  hUnfoldLep1PtEffSigShapeSysDown->Add(hUnfoldLep1Pt,-1.);
  hUnfoldLep1PtEffSigShapeSysDown->Scale(-1.);
  hUnfoldLep1PtEffSigShapeSysDown->Add(hUnfoldLep1Pt,1.);

  TH1D * hUnfoldLep1PtEffBkgShapeSysUp;
  TH1D * hUnfoldLep1PtEffBkgShapeSysDown;
  hUnfoldLep1PtEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[3]->Get("hUnfold"));
  hUnfoldLep1PtEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[3]->Get("hUnfold"));

  hUnfoldLep1PtEffBkgShapeSysDown->Add(hUnfoldLep1Pt,-1.);
  hUnfoldLep1PtEffBkgShapeSysDown->Scale(-1.);
  hUnfoldLep1PtEffBkgShapeSysDown->Add(hUnfoldLep1Pt,1.);

  TH1D * hUnfoldLep1PtUnfoldModelSysUp;
  TH1D * hUnfoldLep1PtUnfoldModelSysDown;
  hUnfoldLep1PtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[3]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep1PtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[3]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep1PtUnfoldModelSysDown->Add(hUnfoldLep1Pt,-1.);
  hUnfoldLep1PtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep1PtUnfoldModelSysDown->Add(hUnfoldLep1Pt,1.);

  
  TGraphAsymmErrors* gUnfoldLep1Pt=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);

  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);
  
  TGraphAsymmErrors* LEP1PT_EWKBKG_UNCERT_BAND_DATA;
  LEP1PT_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEWKBkgDown));

  TGraphAsymmErrors* LEP1PT_TOPBKG_UNCERT_BAND_DATA;
  LEP1PT_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtTopBkgDown));

  TGraphAsymmErrors* LEP1PT_EFFSTAT_UNCERT_BAND_DATA;
  LEP1PT_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffStatDown));
  
  TGraphAsymmErrors* LEP1PT_EFFBIN_UNCERT_BAND_DATA;
  LEP1PT_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffBinSysDown));
  
  TGraphAsymmErrors* LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA;
  LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffSigShapeSysDown));

  TGraphAsymmErrors* LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA;
  LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffBkgShapeSysDown));

  TGraphAsymmErrors* LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_Lep1Pt[7][nBinsLep1Pt];
  double cov_Lep1Pt[nBinsLep1Pt][nBinsLep1Pt];
  
  for(int i=0;i!=nBinsLep1Pt;++i)
    {
      if(hUnfoldLep1PtEWKBkgUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[0][i]=EWKBKG_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      if(hUnfoldLep1PtTopBkgUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[1][i]=TOPBKG_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      if(hUnfoldLep1PtEffStatUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      if(hUnfoldLep1PtEffBinSysUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[3][i]=EFFBIN_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      if(hUnfoldLep1PtEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      if(hUnfoldLep1PtEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      if(hUnfoldLep1PtUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldLep1Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Pt[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
      else
	{
	  uncert_Lep1Pt[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP1PT[i];
	}
    }

  for(int i=0;i!=nBinsLep1Pt;++i)
    {
      for(int j=0;j!=nBinsLep1Pt;++j)
	{
	  cov_Lep1Pt[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_Lep1Pt[i][j]+=uncert_Lep1Pt[k][i]*uncert_Lep1Pt[k][j];
	    }
	  LEP1PT_COV_MATRIX->SetBinContent(i+1,j+1,cov_Lep1Pt[i][j]);
	  LEP1PT_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  LEP1PT_COV_MATRIX->Add(LEP1PT_COV_MATRIX_RESSCALE,1);
  LEP1PT_COV_MATRIX->Add(LEP1PT_COV_MATRIX_MATRIXSTAT,1);
  double corr_Lep1Pt[nBinsLep1Pt][nBinsLep1Pt];
  for(int i=0;i!=nBinsLep1Pt;++i)
    {
      for(int j=0;j!=nBinsLep1Pt;++j)
	{
	  corr_Lep1Pt[i][j]=LEP1PT_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(LEP1PT_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(LEP1PT_COV_MATRIX->GetBinContent(j+1,j+1)));
	  LEP1PT_CORR_MATRIX->SetBinContent(i+1,j+1,corr_Lep1Pt[i][j]);
	  LEP1PT_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
 
  //--------------------------------------------------------------------------
  //                           Lep2 Pt
  //--------------------------------------------------------------------------

  TH1D * hUnfoldLep2Pt=(TH1D*)(file[4]->Get("hUnfold"));
  
  TH1D * hUnfoldLep2PtEWKBkgUp;
  TH1D * hUnfoldLep2PtEWKBkgDown;
  hUnfoldLep2PtEWKBkgUp=(TH1D*)(fileEWKBkgUp[4]->Get("hUnfold"));
  hUnfoldLep2PtEWKBkgDown=(TH1D*)(fileEWKBkgDown[4]->Get("hUnfold"));

  TH1D * hUnfoldLep2PtTopBkgUp;
  TH1D * hUnfoldLep2PtTopBkgDown;
  hUnfoldLep2PtTopBkgUp=(TH1D*)(fileTopBkgUp[4]->Get("hUnfold"));
  hUnfoldLep2PtTopBkgDown=(TH1D*)(fileTopBkgDown[4]->Get("hUnfold"));

  TH1D * hUnfoldLep2PtEffStatUp;
  TH1D * hUnfoldLep2PtEffStatDown;
  hUnfoldLep2PtEffStatUp=(TH1D*)(fileEffStatUp[4]->Get("hUnfold"));
  hUnfoldLep2PtEffStatDown=(TH1D*)(fileEffStatDown[4]->Get("hUnfold"));

  TH1D * hUnfoldLep2PtEffBinSysUp;
  TH1D * hUnfoldLep2PtEffBinSysDown;
  hUnfoldLep2PtEffBinSysUp=(TH1D*)(fileEffBinSys[4]->Get("hUnfold"));
  hUnfoldLep2PtEffBinSysDown=(TH1D*)(fileEffBinSys[4]->Get("hUnfold"));

  hUnfoldLep2PtEffBinSysDown->Add(hUnfoldLep2Pt,-1.);
  hUnfoldLep2PtEffBinSysDown->Scale(-1.);
  hUnfoldLep2PtEffBinSysDown->Add(hUnfoldLep2Pt,1.);

  TH1D * hUnfoldLep2PtEffSigShapeSysUp;
  TH1D * hUnfoldLep2PtEffSigShapeSysDown;
  hUnfoldLep2PtEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[4]->Get("hUnfold"));
  hUnfoldLep2PtEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[4]->Get("hUnfold"));

  hUnfoldLep2PtEffSigShapeSysDown->Add(hUnfoldLep2Pt,-1.);
  hUnfoldLep2PtEffSigShapeSysDown->Scale(-1.);
  hUnfoldLep2PtEffSigShapeSysDown->Add(hUnfoldLep2Pt,1.);

  TH1D * hUnfoldLep2PtEffBkgShapeSysUp;
  TH1D * hUnfoldLep2PtEffBkgShapeSysDown;
  hUnfoldLep2PtEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[4]->Get("hUnfold"));
  hUnfoldLep2PtEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[4]->Get("hUnfold"));

  hUnfoldLep2PtEffBkgShapeSysDown->Add(hUnfoldLep2Pt,-1.);
  hUnfoldLep2PtEffBkgShapeSysDown->Scale(-1.);
  hUnfoldLep2PtEffBkgShapeSysDown->Add(hUnfoldLep2Pt,1.);

  TH1D * hUnfoldLep2PtUnfoldModelSysUp;
  TH1D * hUnfoldLep2PtUnfoldModelSysDown;
  hUnfoldLep2PtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[4]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep2PtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[4]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep2PtUnfoldModelSysDown->Add(hUnfoldLep2Pt,-1.);
  hUnfoldLep2PtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep2PtUnfoldModelSysDown->Add(hUnfoldLep2Pt,1.);

  
  TGraphAsymmErrors* gUnfoldLep2Pt=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);
  
  TGraphAsymmErrors* LEP2PT_EWKBKG_UNCERT_BAND_DATA;
  LEP2PT_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEWKBkgDown));

  TGraphAsymmErrors* LEP2PT_TOPBKG_UNCERT_BAND_DATA;
  LEP2PT_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtTopBkgDown));

  TGraphAsymmErrors* LEP2PT_EFFSTAT_UNCERT_BAND_DATA;
  LEP2PT_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffStatDown));
  
  TGraphAsymmErrors* LEP2PT_EFFBIN_UNCERT_BAND_DATA;
  LEP2PT_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffBinSysDown));
  
  TGraphAsymmErrors* LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA;
  LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffSigShapeSysDown));

  TGraphAsymmErrors* LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA;
  LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffBkgShapeSysDown));

  TGraphAsymmErrors* LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_Lep2Pt[7][nBinsLep2Pt];
  double cov_Lep2Pt[nBinsLep2Pt][nBinsLep2Pt];
  
  for(int i=0;i!=nBinsLep2Pt;++i)
    {
      if(hUnfoldLep2PtEWKBkgUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[0][i]=EWKBKG_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      if(hUnfoldLep2PtTopBkgUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[1][i]=TOPBKG_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      if(hUnfoldLep2PtEffStatUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      if(hUnfoldLep2PtEffBinSysUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[3][i]=EFFBIN_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      if(hUnfoldLep2PtEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      if(hUnfoldLep2PtEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      if(hUnfoldLep2PtUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldLep2Pt->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Pt[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
      else
	{
	  uncert_Lep2Pt[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP2PT[i];
	}
    }

  for(int i=0;i!=nBinsLep2Pt;++i)
    {
      for(int j=0;j!=nBinsLep2Pt;++j)
	{
	  cov_Lep2Pt[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_Lep2Pt[i][j]+=uncert_Lep2Pt[k][i]*uncert_Lep2Pt[k][j];
	    }
	  LEP2PT_COV_MATRIX->SetBinContent(i+1,j+1,cov_Lep2Pt[i][j]);
	  LEP2PT_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  LEP2PT_COV_MATRIX->Add(LEP2PT_COV_MATRIX_RESSCALE,1);
  LEP2PT_COV_MATRIX->Add(LEP2PT_COV_MATRIX_MATRIXSTAT,1);
  double corr_Lep2Pt[nBinsLep2Pt][nBinsLep2Pt];
  for(int i=0;i!=nBinsLep2Pt;++i)
    {
      for(int j=0;j!=nBinsLep2Pt;++j)
	{
	  corr_Lep2Pt[i][j]=LEP2PT_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(LEP2PT_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(LEP2PT_COV_MATRIX->GetBinContent(j+1,j+1)));
	  LEP2PT_CORR_MATRIX->SetBinContent(i+1,j+1,corr_Lep2Pt[i][j]);
	  LEP2PT_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }

  //--------------------------------------------------------------------------
  //                           Lep1 Eta
  //--------------------------------------------------------------------------

  TH1D * hUnfoldLep1Eta=(TH1D*)(file[5]->Get("hUnfold"));
  
  TH1D * hUnfoldLep1EtaEWKBkgUp;
  TH1D * hUnfoldLep1EtaEWKBkgDown;
  hUnfoldLep1EtaEWKBkgUp=(TH1D*)(fileEWKBkgUp[5]->Get("hUnfold"));
  hUnfoldLep1EtaEWKBkgDown=(TH1D*)(fileEWKBkgDown[5]->Get("hUnfold"));

  TH1D * hUnfoldLep1EtaTopBkgUp;
  TH1D * hUnfoldLep1EtaTopBkgDown;
  hUnfoldLep1EtaTopBkgUp=(TH1D*)(fileTopBkgUp[5]->Get("hUnfold"));
  hUnfoldLep1EtaTopBkgDown=(TH1D*)(fileTopBkgDown[5]->Get("hUnfold"));

  TH1D * hUnfoldLep1EtaEffStatUp;
  TH1D * hUnfoldLep1EtaEffStatDown;
  hUnfoldLep1EtaEffStatUp=(TH1D*)(fileEffStatUp[5]->Get("hUnfold"));
  hUnfoldLep1EtaEffStatDown=(TH1D*)(fileEffStatDown[5]->Get("hUnfold"));

  TH1D * hUnfoldLep1EtaEffBinSysUp;
  TH1D * hUnfoldLep1EtaEffBinSysDown;
  hUnfoldLep1EtaEffBinSysUp=(TH1D*)(fileEffBinSys[5]->Get("hUnfold"));
  hUnfoldLep1EtaEffBinSysDown=(TH1D*)(fileEffBinSys[5]->Get("hUnfold"));

  hUnfoldLep1EtaEffBinSysDown->Add(hUnfoldLep1Eta,-1.);
  hUnfoldLep1EtaEffBinSysDown->Scale(-1.);
  hUnfoldLep1EtaEffBinSysDown->Add(hUnfoldLep1Eta,1.);

  TH1D * hUnfoldLep1EtaEffSigShapeSysUp;
  TH1D * hUnfoldLep1EtaEffSigShapeSysDown;
  hUnfoldLep1EtaEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[5]->Get("hUnfold"));
  hUnfoldLep1EtaEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[5]->Get("hUnfold"));

  hUnfoldLep1EtaEffSigShapeSysDown->Add(hUnfoldLep1Eta,-1.);
  hUnfoldLep1EtaEffSigShapeSysDown->Scale(-1.);
  hUnfoldLep1EtaEffSigShapeSysDown->Add(hUnfoldLep1Eta,1.);

  TH1D * hUnfoldLep1EtaEffBkgShapeSysUp;
  TH1D * hUnfoldLep1EtaEffBkgShapeSysDown;
  hUnfoldLep1EtaEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[5]->Get("hUnfold"));
  hUnfoldLep1EtaEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[5]->Get("hUnfold"));

  hUnfoldLep1EtaEffBkgShapeSysDown->Add(hUnfoldLep1Eta,-1.);
  hUnfoldLep1EtaEffBkgShapeSysDown->Scale(-1.);
  hUnfoldLep1EtaEffBkgShapeSysDown->Add(hUnfoldLep1Eta,1.);

  TH1D * hUnfoldLep1EtaUnfoldModelSysUp;
  TH1D * hUnfoldLep1EtaUnfoldModelSysDown;
  hUnfoldLep1EtaUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[5]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep1EtaUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[5]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep1EtaUnfoldModelSysDown->Add(hUnfoldLep1Eta,-1.);
  hUnfoldLep1EtaUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep1EtaUnfoldModelSysDown->Add(hUnfoldLep1Eta,1.);

  
  TGraphAsymmErrors* gUnfoldLep1Eta=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);
  
  TGraphAsymmErrors* LEP1ETA_EWKBKG_UNCERT_BAND_DATA;
  LEP1ETA_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEWKBkgDown));

  TGraphAsymmErrors* LEP1ETA_TOPBKG_UNCERT_BAND_DATA;
  LEP1ETA_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaTopBkgDown));

  TGraphAsymmErrors* LEP1ETA_EFFSTAT_UNCERT_BAND_DATA;
  LEP1ETA_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffStatDown));
  
  TGraphAsymmErrors* LEP1ETA_EFFBIN_UNCERT_BAND_DATA;
  LEP1ETA_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffBinSysDown));
  
  TGraphAsymmErrors* LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA;
  LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffSigShapeSysDown));

  TGraphAsymmErrors* LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA;
  LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffBkgShapeSysDown));

  TGraphAsymmErrors* LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_Lep1Eta[7][24];
  double cov_Lep1Eta[24][24];
  
  for(int i=0;i!=24;++i)
    {
      if(hUnfoldLep1EtaEWKBkgUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[0][i]=EWKBKG_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      if(hUnfoldLep1EtaTopBkgUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[1][i]=TOPBKG_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      if(hUnfoldLep1EtaEffStatUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      if(hUnfoldLep1EtaEffBinSysUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[3][i]=EFFBIN_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      if(hUnfoldLep1EtaEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      if(hUnfoldLep1EtaEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      if(hUnfoldLep1EtaUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldLep1Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep1Eta[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
      else
	{
	  uncert_Lep1Eta[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP1ETA[i];
	}
    }

  for(int i=0;i!=24;++i)
    {
      for(int j=0;j!=24;++j)
	{
	  cov_Lep1Eta[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_Lep1Eta[i][j]+=uncert_Lep1Eta[k][i]*uncert_Lep1Eta[k][j];
	    }
	  LEP1ETA_COV_MATRIX->SetBinContent(i+1,j+1,cov_Lep1Eta[i][j]);
	  LEP1ETA_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  LEP1ETA_COV_MATRIX->Add(LEP1ETA_COV_MATRIX_RESSCALE,1);
  LEP1ETA_COV_MATRIX->Add(LEP1ETA_COV_MATRIX_MATRIXSTAT,1);
  double corr_Lep1Eta[24][24];
  for(int i=0;i!=24;++i)
    {
      for(int j=0;j!=24;++j)
	{
	  corr_Lep1Eta[i][j]=LEP1ETA_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(LEP1ETA_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(LEP1ETA_COV_MATRIX->GetBinContent(j+1,j+1)));
	  LEP1ETA_CORR_MATRIX->SetBinContent(i+1,j+1,corr_Lep1Eta[i][j]);
	  LEP1ETA_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
   
  //--------------------------------------------------------------------------
  //                           Lep2 Eta
  //--------------------------------------------------------------------------

  TH1D * hUnfoldLep2Eta=(TH1D*)(file[6]->Get("hUnfold"));
  
  TH1D * hUnfoldLep2EtaEWKBkgUp;
  TH1D * hUnfoldLep2EtaEWKBkgDown;
  hUnfoldLep2EtaEWKBkgUp=(TH1D*)(fileEWKBkgUp[6]->Get("hUnfold"));
  hUnfoldLep2EtaEWKBkgDown=(TH1D*)(fileEWKBkgDown[6]->Get("hUnfold"));

  TH1D * hUnfoldLep2EtaTopBkgUp;
  TH1D * hUnfoldLep2EtaTopBkgDown;
  hUnfoldLep2EtaTopBkgUp=(TH1D*)(fileTopBkgUp[6]->Get("hUnfold"));
  hUnfoldLep2EtaTopBkgDown=(TH1D*)(fileTopBkgDown[6]->Get("hUnfold"));

  TH1D * hUnfoldLep2EtaEffStatUp;
  TH1D * hUnfoldLep2EtaEffStatDown;
  hUnfoldLep2EtaEffStatUp=(TH1D*)(fileEffStatUp[6]->Get("hUnfold"));
  hUnfoldLep2EtaEffStatDown=(TH1D*)(fileEffStatDown[6]->Get("hUnfold"));

  TH1D * hUnfoldLep2EtaEffBinSysUp;
  TH1D * hUnfoldLep2EtaEffBinSysDown;
  hUnfoldLep2EtaEffBinSysUp=(TH1D*)(fileEffBinSys[6]->Get("hUnfold"));
  hUnfoldLep2EtaEffBinSysDown=(TH1D*)(fileEffBinSys[6]->Get("hUnfold"));

  hUnfoldLep2EtaEffBinSysDown->Add(hUnfoldLep2Eta,-1.);
  hUnfoldLep2EtaEffBinSysDown->Scale(-1.);
  hUnfoldLep2EtaEffBinSysDown->Add(hUnfoldLep2Eta,1.);

  TH1D * hUnfoldLep2EtaEffSigShapeSysUp;
  TH1D * hUnfoldLep2EtaEffSigShapeSysDown;
  hUnfoldLep2EtaEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[6]->Get("hUnfold"));
  hUnfoldLep2EtaEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[6]->Get("hUnfold"));

  hUnfoldLep2EtaEffSigShapeSysDown->Add(hUnfoldLep2Eta,-1.);
  hUnfoldLep2EtaEffSigShapeSysDown->Scale(-1.);
  hUnfoldLep2EtaEffSigShapeSysDown->Add(hUnfoldLep2Eta,1.);

  TH1D * hUnfoldLep2EtaEffBkgShapeSysUp;
  TH1D * hUnfoldLep2EtaEffBkgShapeSysDown;
  hUnfoldLep2EtaEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[6]->Get("hUnfold"));
  hUnfoldLep2EtaEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[6]->Get("hUnfold"));

  hUnfoldLep2EtaEffBkgShapeSysDown->Add(hUnfoldLep2Eta,-1.);
  hUnfoldLep2EtaEffBkgShapeSysDown->Scale(-1.);
  hUnfoldLep2EtaEffBkgShapeSysDown->Add(hUnfoldLep2Eta,1.);

  TH1D * hUnfoldLep2EtaUnfoldModelSysUp;
  TH1D * hUnfoldLep2EtaUnfoldModelSysDown;
  hUnfoldLep2EtaUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[6]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep2EtaUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[6]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep2EtaUnfoldModelSysDown->Add(hUnfoldLep2Eta,-1.);
  hUnfoldLep2EtaUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep2EtaUnfoldModelSysDown->Add(hUnfoldLep2Eta,1.);

  
  TGraphAsymmErrors* gUnfoldLep2Eta=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);

  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);
  
  TGraphAsymmErrors* LEP2ETA_EWKBKG_UNCERT_BAND_DATA;
  LEP2ETA_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEWKBkgDown));

  TGraphAsymmErrors* LEP2ETA_TOPBKG_UNCERT_BAND_DATA;
  LEP2ETA_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaTopBkgDown));

  TGraphAsymmErrors* LEP2ETA_EFFSTAT_UNCERT_BAND_DATA;
  LEP2ETA_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffStatDown));
  
  TGraphAsymmErrors* LEP2ETA_EFFBIN_UNCERT_BAND_DATA;
  LEP2ETA_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffBinSysDown));
  
  TGraphAsymmErrors* LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA;
  LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffSigShapeSysDown));

  TGraphAsymmErrors* LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA;
  LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffBkgShapeSysDown));

  TGraphAsymmErrors* LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_Lep2Eta[7][24];
  double cov_Lep2Eta[24][24];
  
  for(int i=0;i!=24;++i)
    {
      if(hUnfoldLep2EtaEWKBkgUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[0][i]=EWKBKG_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      if(hUnfoldLep2EtaTopBkgUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[1][i]=TOPBKG_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      if(hUnfoldLep2EtaEffStatUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      if(hUnfoldLep2EtaEffBinSysUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[3][i]=EFFBIN_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      if(hUnfoldLep2EtaEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      if(hUnfoldLep2EtaEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      if(hUnfoldLep2EtaUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldLep2Eta->GetBinContent(i+1)<1)
	{
	  uncert_Lep2Eta[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
      else
	{
	  uncert_Lep2Eta[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_LEP2ETA[i];
	}
    }

  for(int i=0;i!=24;++i)
    {
      for(int j=0;j!=24;++j)
	{
	  cov_Lep2Eta[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_Lep2Eta[i][j]+=uncert_Lep2Eta[k][i]*uncert_Lep2Eta[k][j];
	    }
	  LEP2ETA_COV_MATRIX->SetBinContent(i+1,j+1,cov_Lep2Eta[i][j]);
	  LEP2ETA_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  LEP2ETA_COV_MATRIX->Add(LEP2ETA_COV_MATRIX_RESSCALE,1);
  LEP2ETA_COV_MATRIX->Add(LEP2ETA_COV_MATRIX_MATRIXSTAT,1);
  double corr_Lep2Eta[24][24];
  for(int i=0;i!=24;++i)
    {
      for(int j=0;j!=24;++j)
	{
	  corr_Lep2Eta[i][j]=LEP2ETA_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(LEP2ETA_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(LEP2ETA_COV_MATRIX->GetBinContent(j+1,j+1)));
	  LEP2ETA_CORR_MATRIX->SetBinContent(i+1,j+1,corr_Lep2Eta[i][j]);
	  LEP2ETA_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }

  

  //--------------------------------------------------------------------------
  //                           LepNeg Pt
  //--------------------------------------------------------------------------

  TH1D * hUnfoldLepNegPt=(TH1D*)(file[7]->Get("hUnfold"));
  
  TH1D * hUnfoldLepNegPtEWKBkgUp;
  TH1D * hUnfoldLepNegPtEWKBkgDown;
  hUnfoldLepNegPtEWKBkgUp=(TH1D*)(fileEWKBkgUp[7]->Get("hUnfold"));
  hUnfoldLepNegPtEWKBkgDown=(TH1D*)(fileEWKBkgDown[7]->Get("hUnfold"));

  TH1D * hUnfoldLepNegPtTopBkgUp;
  TH1D * hUnfoldLepNegPtTopBkgDown;
  hUnfoldLepNegPtTopBkgUp=(TH1D*)(fileTopBkgUp[7]->Get("hUnfold"));
  hUnfoldLepNegPtTopBkgDown=(TH1D*)(fileTopBkgDown[7]->Get("hUnfold"));

  TH1D * hUnfoldLepNegPtEffStatUp;
  TH1D * hUnfoldLepNegPtEffStatDown;
  hUnfoldLepNegPtEffStatUp=(TH1D*)(fileEffStatUp[7]->Get("hUnfold"));
  hUnfoldLepNegPtEffStatDown=(TH1D*)(fileEffStatDown[7]->Get("hUnfold"));

  TH1D * hUnfoldLepNegPtEffBinSysUp;
  TH1D * hUnfoldLepNegPtEffBinSysDown;
  hUnfoldLepNegPtEffBinSysUp=(TH1D*)(fileEffBinSys[7]->Get("hUnfold"));
  hUnfoldLepNegPtEffBinSysDown=(TH1D*)(fileEffBinSys[7]->Get("hUnfold"));

  hUnfoldLepNegPtEffBinSysDown->Add(hUnfoldLepNegPt,-1.);
  hUnfoldLepNegPtEffBinSysDown->Scale(-1.);
  hUnfoldLepNegPtEffBinSysDown->Add(hUnfoldLepNegPt,1.);

  TH1D * hUnfoldLepNegPtEffSigShapeSysUp;
  TH1D * hUnfoldLepNegPtEffSigShapeSysDown;
  hUnfoldLepNegPtEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[7]->Get("hUnfold"));
  hUnfoldLepNegPtEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[7]->Get("hUnfold"));

  hUnfoldLepNegPtEffSigShapeSysDown->Add(hUnfoldLepNegPt,-1.);
  hUnfoldLepNegPtEffSigShapeSysDown->Scale(-1.);
  hUnfoldLepNegPtEffSigShapeSysDown->Add(hUnfoldLepNegPt,1.);

  TH1D * hUnfoldLepNegPtEffBkgShapeSysUp;
  TH1D * hUnfoldLepNegPtEffBkgShapeSysDown;
  hUnfoldLepNegPtEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[7]->Get("hUnfold"));
  hUnfoldLepNegPtEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[7]->Get("hUnfold"));

  hUnfoldLepNegPtEffBkgShapeSysDown->Add(hUnfoldLepNegPt,-1.);
  hUnfoldLepNegPtEffBkgShapeSysDown->Scale(-1.);
  hUnfoldLepNegPtEffBkgShapeSysDown->Add(hUnfoldLepNegPt,1.);

  TH1D * hUnfoldLepNegPtUnfoldModelSysUp;
  TH1D * hUnfoldLepNegPtUnfoldModelSysDown;
  hUnfoldLepNegPtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[7]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLepNegPtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[7]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLepNegPtUnfoldModelSysDown->Add(hUnfoldLepNegPt,-1.);
  hUnfoldLepNegPtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLepNegPtUnfoldModelSysDown->Add(hUnfoldLepNegPt,1.);

  
  TGraphAsymmErrors* gUnfoldLepNegPt=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);

  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);
  
  TGraphAsymmErrors* LEPNEGPT_EWKBKG_UNCERT_BAND_DATA;
  LEPNEGPT_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEWKBkgDown));

  TGraphAsymmErrors* LEPNEGPT_TOPBKG_UNCERT_BAND_DATA;
  LEPNEGPT_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtTopBkgDown));

  TGraphAsymmErrors* LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA;
  LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffStatDown));
  
  TGraphAsymmErrors* LEPNEGPT_EFFBIN_UNCERT_BAND_DATA;
  LEPNEGPT_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffBinSysDown));
  
  TGraphAsymmErrors* LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA;
  LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffSigShapeSysDown));

  TGraphAsymmErrors* LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA;
  LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffBkgShapeSysDown));

  TGraphAsymmErrors* LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_LepNegPt[7][nBinsLepNegPt];
  double cov_LepNegPt[nBinsLepNegPt][nBinsLepNegPt];
  
  for(int i=0;i!=nBinsLepNegPt;++i)
    {
      if(hUnfoldLepNegPtEWKBkgUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[0][i]=EWKBKG_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      if(hUnfoldLepNegPtTopBkgUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[1][i]=TOPBKG_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      if(hUnfoldLepNegPtEffStatUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      if(hUnfoldLepNegPtEffBinSysUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[3][i]=EFFBIN_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      if(hUnfoldLepNegPtEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      if(hUnfoldLepNegPtEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      if(hUnfoldLepNegPtUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldLepNegPt->GetBinContent(i+1)<1)
	{
	  uncert_LepNegPt[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
      else
	{
	  uncert_LepNegPt[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_LEPNEGPT[i];
	}
    }

  for(int i=0;i!=nBinsLepNegPt;++i)
    {
      for(int j=0;j!=nBinsLepNegPt;++j)
	{
	  cov_LepNegPt[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_LepNegPt[i][j]+=uncert_LepNegPt[k][i]*uncert_LepNegPt[k][j];
	    }
	  LEPNEGPT_COV_MATRIX->SetBinContent(i+1,j+1,cov_LepNegPt[i][j]);
	  LEPNEGPT_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  LEPNEGPT_COV_MATRIX->Add(LEPNEGPT_COV_MATRIX_RESSCALE,1);
  LEPNEGPT_COV_MATRIX->Add(LEPNEGPT_COV_MATRIX_MATRIXSTAT,1);
  double corr_LepNegPt[nBinsLepNegPt][nBinsLepNegPt];
  for(int i=0;i!=nBinsLepNegPt;++i)
    {
      for(int j=0;j!=nBinsLepNegPt;++j)
	{
	  corr_LepNegPt[i][j]=LEPNEGPT_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(LEPNEGPT_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(LEPNEGPT_COV_MATRIX->GetBinContent(j+1,j+1)));
	  LEPNEGPT_CORR_MATRIX->SetBinContent(i+1,j+1,corr_LepNegPt[i][j]);
	  LEPNEGPT_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }

  //--------------------------------------------------------------------------
  //                           LepPos Pt
  //--------------------------------------------------------------------------

  TH1D * hUnfoldLepPosPt=(TH1D*)(file[8]->Get("hUnfold"));
  
  TH1D * hUnfoldLepPosPtEWKBkgUp;
  TH1D * hUnfoldLepPosPtEWKBkgDown;
  hUnfoldLepPosPtEWKBkgUp=(TH1D*)(fileEWKBkgUp[8]->Get("hUnfold"));
  hUnfoldLepPosPtEWKBkgDown=(TH1D*)(fileEWKBkgDown[8]->Get("hUnfold"));

  TH1D * hUnfoldLepPosPtTopBkgUp;
  TH1D * hUnfoldLepPosPtTopBkgDown;
  hUnfoldLepPosPtTopBkgUp=(TH1D*)(fileTopBkgUp[8]->Get("hUnfold"));
  hUnfoldLepPosPtTopBkgDown=(TH1D*)(fileTopBkgDown[8]->Get("hUnfold"));

  TH1D * hUnfoldLepPosPtEffStatUp;
  TH1D * hUnfoldLepPosPtEffStatDown;
  hUnfoldLepPosPtEffStatUp=(TH1D*)(fileEffStatUp[8]->Get("hUnfold"));
  hUnfoldLepPosPtEffStatDown=(TH1D*)(fileEffStatDown[8]->Get("hUnfold"));

  TH1D * hUnfoldLepPosPtEffBinSysUp;
  TH1D * hUnfoldLepPosPtEffBinSysDown;
  hUnfoldLepPosPtEffBinSysUp=(TH1D*)(fileEffBinSys[8]->Get("hUnfold"));
  hUnfoldLepPosPtEffBinSysDown=(TH1D*)(fileEffBinSys[8]->Get("hUnfold"));

  hUnfoldLepPosPtEffBinSysDown->Add(hUnfoldLepPosPt,-1.);
  hUnfoldLepPosPtEffBinSysDown->Scale(-1.);
  hUnfoldLepPosPtEffBinSysDown->Add(hUnfoldLepPosPt,1.);

  TH1D * hUnfoldLepPosPtEffSigShapeSysUp;
  TH1D * hUnfoldLepPosPtEffSigShapeSysDown;
  hUnfoldLepPosPtEffSigShapeSysUp=(TH1D*)(fileEffSigShapeSys[8]->Get("hUnfold"));
  hUnfoldLepPosPtEffSigShapeSysDown=(TH1D*)(fileEffSigShapeSys[8]->Get("hUnfold"));

  hUnfoldLepPosPtEffSigShapeSysDown->Add(hUnfoldLepPosPt,-1.);
  hUnfoldLepPosPtEffSigShapeSysDown->Scale(-1.);
  hUnfoldLepPosPtEffSigShapeSysDown->Add(hUnfoldLepPosPt,1.);

  TH1D * hUnfoldLepPosPtEffBkgShapeSysUp;
  TH1D * hUnfoldLepPosPtEffBkgShapeSysDown;
  hUnfoldLepPosPtEffBkgShapeSysUp=(TH1D*)(fileEffBkgShapeSys[8]->Get("hUnfold"));
  hUnfoldLepPosPtEffBkgShapeSysDown=(TH1D*)(fileEffBkgShapeSys[8]->Get("hUnfold"));

  hUnfoldLepPosPtEffBkgShapeSysDown->Add(hUnfoldLepPosPt,-1.);
  hUnfoldLepPosPtEffBkgShapeSysDown->Scale(-1.);
  hUnfoldLepPosPtEffBkgShapeSysDown->Add(hUnfoldLepPosPt,1.);

  TH1D * hUnfoldLepPosPtUnfoldModelSysUp;
  TH1D * hUnfoldLepPosPtUnfoldModelSysDown;
  hUnfoldLepPosPtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[8]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLepPosPtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[8]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLepPosPtUnfoldModelSysDown->Add(hUnfoldLepPosPt,-1.);
  hUnfoldLepPosPtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLepPosPtUnfoldModelSysDown->Add(hUnfoldLepPosPt,1.);

  
  TGraphAsymmErrors* gUnfoldLepPosPt=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);
  
  TGraphAsymmErrors* LEPPOSPT_EWKBKG_UNCERT_BAND_DATA;
  LEPPOSPT_EWKBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEWKBkgDown));

  TGraphAsymmErrors* LEPPOSPT_TOPBKG_UNCERT_BAND_DATA;
  LEPPOSPT_TOPBKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtTopBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtTopBkgDown));

  TGraphAsymmErrors* LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA;
  LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffStatDown));
  
  TGraphAsymmErrors* LEPPOSPT_EFFBIN_UNCERT_BAND_DATA;
  LEPPOSPT_EFFBIN_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffBinSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffBinSysDown));
  
  TGraphAsymmErrors* LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA;
  LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffSigShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffSigShapeSysDown));

  TGraphAsymmErrors* LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA;
  LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffBkgShapeSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffBkgShapeSysDown));

  TGraphAsymmErrors* LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysDown));
  
  double* EWKBKG_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EWKBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* TOPBKG_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_TOPBKG_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSTAT_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBIN_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EFFBIN_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA->GetEYhigh();

  double* UNFOLDMODEL_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA->GetEYhigh();

  double uncert_LepPosPt[7][nBinsLepPosPt];
  double cov_LepPosPt[nBinsLepPosPt][nBinsLepPosPt];
  
  for(int i=0;i!=nBinsLepPosPt;++i)
    {
      if(hUnfoldLepPosPtEWKBkgUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[0][i]=-EWKBKG_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[0][i]=EWKBKG_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      if(hUnfoldLepPosPtTopBkgUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[1][i]=-TOPBKG_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[1][i]=TOPBKG_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      if(hUnfoldLepPosPtEffStatUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[2][i]=-EFFSTAT_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[2][i]=EFFSTAT_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      if(hUnfoldLepPosPtEffBinSysUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[3][i]=-EFFBIN_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[3][i]=EFFBIN_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      if(hUnfoldLepPosPtEffSigShapeSysUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[4][i]=-EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[4][i]=EFFSIGSHAPE_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      if(hUnfoldLepPosPtEffBkgShapeSysUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[5][i]=-EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[5][i]=EFFBKGSHAPE_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      if(hUnfoldLepPosPtUnfoldModelSysUp->GetBinContent(i+1)/hUnfoldLepPosPt->GetBinContent(i+1)<1)
	{
	  uncert_LepPosPt[6][i]=-UNFOLDMODEL_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
      else
	{
	  uncert_LepPosPt[6][i]=UNFOLDMODEL_SYSTEMATIC_UNCERT_LEPPOSPT[i];
	}
    }

  for(int i=0;i!=nBinsLepPosPt;++i)
    {
      for(int j=0;j!=nBinsLepPosPt;++j)
	{
	  cov_LepPosPt[i][j]=0;
	  for(int k=0;k!=7;++k)
	    {
	      cov_LepPosPt[i][j]+=uncert_LepPosPt[k][i]*uncert_LepPosPt[k][j];
	    }
	  LEPPOSPT_COV_MATRIX->SetBinContent(i+1,j+1,cov_LepPosPt[i][j]);
	  LEPPOSPT_COV_MATRIX->SetBinError(i+1,j+1,0);
	}
    }
  LEPPOSPT_COV_MATRIX->Add(LEPPOSPT_COV_MATRIX_RESSCALE,1);
  LEPPOSPT_COV_MATRIX->Add(LEPPOSPT_COV_MATRIX_MATRIXSTAT,1);
  double corr_LepPosPt[nBinsLepPosPt][nBinsLepPosPt];
  for(int i=0;i!=nBinsLepPosPt;++i)
    {
      for(int j=0;j!=nBinsLepPosPt;++j)
	{
	  corr_LepPosPt[i][j]=LEPPOSPT_COV_MATRIX->GetBinContent(i+1,j+1)/(sqrt(LEPPOSPT_COV_MATRIX->GetBinContent(i+1,i+1))*sqrt(LEPPOSPT_COV_MATRIX->GetBinContent(j+1,j+1)));
	  LEPPOSPT_CORR_MATRIX->SetBinContent(i+1,j+1,corr_LepPosPt[i][j]);
	  LEPPOSPT_CORR_MATRIX->SetBinError(i+1,j+1,0);
	}
    }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char xlabel[100];     // string buffer for x-axis label
  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->cd();
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.15);
  c->SetLeftMargin(0.15);  
  c->SetRightMargin(0.15);  
  c->SetTickx(1);
  c->SetTicky(1);  
  TGaxis::SetMaxDigits(3);

  gStyle->SetTitleOffset(1.4,"Y");
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("4.2f");

  
  
  //
  // ZPt
  //   

  sprintf(ylabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  ZPT_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  ZPT_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  CPlot plotZmumuPt("zmmPtCorrelations","",xlabel,ylabel);
  plotZmumuPt.AddHist2D(ZPT_CORR_MATRIX,"COLZ");
  plotZmumuPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuPt.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuPt.SetLogx();
  plotZmumuPt.SetLogy();
  plotZmumuPt.Draw(c,kTRUE,format);

  //
  // PhiStar
  //   

  sprintf(ylabel,"#phi_{#eta}*");
  sprintf(xlabel,"#phi_{#eta}*");
  PHISTAR_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  PHISTAR_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  CPlot plotZmumuPhiStar("zmmPhiStarCorrelations","",xlabel,ylabel);
  plotZmumuPhiStar.AddHist2D(PHISTAR_CORR_MATRIX,"COLZ");
  plotZmumuPhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuPhiStar.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuPhiStar.SetLogx();
  plotZmumuPhiStar.SetLogy();
  plotZmumuPhiStar.Draw(c,kTRUE,format);

  //
  // ZRap
  //   

  sprintf(ylabel,"|y^{#mu^{+}#mu^{-}}|");
  sprintf(xlabel,"|y^{#mu^{+}#mu^{-}}|");
  ZRAP_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  ZRAP_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  CPlot plotZmumuRap("zmmRapCorrelations","",xlabel,ylabel);
  plotZmumuRap.AddHist2D(ZRAP_CORR_MATRIX,"COLZ");
  plotZmumuRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuRap.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuRap.Draw(c,kTRUE,format);

  //
  // Lep1Pt
  //   

  sprintf(ylabel,"p_{T} (leading muon) [GeV]");
  sprintf(xlabel,"p_{T} (leading muon) [GeV]");
  LEP1PT_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  LEP1PT_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  CPlot plotZmumuLep1Pt("zmmLep1PtCorrelations","",xlabel,ylabel);
  plotZmumuLep1Pt.AddHist2D(LEP1PT_CORR_MATRIX,"COLZ");
  plotZmumuLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuLep1Pt.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuLep1Pt.SetLogx();
  plotZmumuLep1Pt.SetLogy();
  plotZmumuLep1Pt.Draw(c,kTRUE,format);

  //
  // Lep2Pt
  //   

  sprintf(ylabel,"p_{T} (2nd leading muon) [GeV]");
  sprintf(xlabel,"p_{T} (2nd leading muon) [GeV]");
  LEP2PT_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  LEP2PT_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  CPlot plotZmumuLep2Pt("zmmLep2PtCorrelations","",xlabel,ylabel);
  plotZmumuLep2Pt.AddHist2D(LEP2PT_CORR_MATRIX,"COLZ");
  plotZmumuLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuLep2Pt.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuLep2Pt.SetLogx();
  plotZmumuLep2Pt.SetLogy();
  plotZmumuLep2Pt.Draw(c,kTRUE,format);

  //
  // Lep1Eta
  //   

  sprintf(ylabel,"|#eta| (leading muon)");
  sprintf(xlabel,"|#eta| (leading muon)");
  LEP1ETA_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  LEP1ETA_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  CPlot plotZmumuLep1Eta("zmmLep1EtaCorrelations","",xlabel,ylabel);
  plotZmumuLep1Eta.AddHist2D(LEP1ETA_CORR_MATRIX,"COLZ");
  plotZmumuLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuLep1Eta.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuLep1Eta.Draw(c,kTRUE,format);

  //
  // Lep2Eta
  //   

  sprintf(ylabel,"|#eta| (2nd leading muon)");
  sprintf(xlabel,"|#eta| (2nd leading muon)");
  LEP2ETA_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  LEP2ETA_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  CPlot plotZmumuLep2Eta("zmmLep2EtaCorrelations","",xlabel,ylabel);
  plotZmumuLep2Eta.AddHist2D(LEP2ETA_CORR_MATRIX,"COLZ");
  plotZmumuLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuLep2Eta.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuLep2Eta.Draw(c,kTRUE,format);

  //
  // LepNegPt
  //   

  sprintf(ylabel,"p_{T}^{#mu^{-}} [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{-}} [GeV]");
  LEP1PT_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  LEP1PT_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  CPlot plotZmumuLepNegPt("zmmLepNegPtCorrelations","",xlabel,ylabel);
  plotZmumuLepNegPt.AddHist2D(LEP1PT_CORR_MATRIX,"COLZ");
  plotZmumuLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuLepNegPt.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuLepNegPt.SetLogx();
  plotZmumuLepNegPt.SetLogy();
  plotZmumuLepNegPt.Draw(c,kTRUE,format);

  //
  // LepPosPt
  //   

  sprintf(ylabel,"p_{T}^{#mu^{+}} [GeV]");
  sprintf(xlabel,"p_{T}^{#mu^{+}} [GeV]");
  LEP1PT_CORR_MATRIX->GetZaxis()->SetRangeUser(-0.5,1);
  LEP1PT_CORR_MATRIX->GetYaxis()->SetTitleOffset(1.25);
  CPlot plotZmumuLepPosPt("zmmLepPosPtCorrelations","",xlabel,ylabel);
  plotZmumuLepPosPt.AddHist2D(LEP1PT_CORR_MATRIX,"COLZ");
  plotZmumuLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.15,0.90,0.4,0.95,0);
  plotZmumuLepPosPt.AddTextBox("Correlation Matrix",0.63,0.90,0.86,0.95,0);
  plotZmumuLepPosPt.SetLogx();
  plotZmumuLepPosPt.SetLogy();
  plotZmumuLepPosPt.Draw(c,kTRUE,format);

  
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  outFile->cd();
  outFile->Write();
  outFile->Close(); 


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;
 
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     

  gBenchmark->Show("plotZmmCorrelations");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = (TH1D*)hData->Clone("hDiff");
  hDiff->SetName(name);
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff=0;
    Double_t err=0;
    if(hData->GetBinContent(ibin)!=0)
      {
	diff = hFit->GetBinContent(ibin)/hData->GetBinContent(ibin);
	err = hFit->GetBinError(ibin)/hData->GetBinContent(ibin);
      }
    hDiff->SetBinContent(ibin,diff);
    hDiff->SetBinError(ibin,err);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.55);
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

TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1){


 if (!h1) cout << "TH1TOTGraph: histogram not found !" << endl;

 TGraphAsymmErrors* g1= new TGraphAsymmErrors();

 Double_t x, y, exh,exl, eyh,eyl;
 for (Int_t i=0; i<h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i+1);
   eyl=h1->GetBinError(i+1);
   eyh=h1->GetBinError(i+1);
   x=h1->GetBinCenter(i+1);
   exl=h1->GetBinWidth(i+1)/2;
   exh=h1->GetBinWidth(i+1)/2;

   g1->SetPoint(i,x,y);
   g1->SetPointError(i,exl,exh,eyl,eyh);

 }
 return g1;
}

TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {
   TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;

  Double_t dum;
  for (Int_t i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    Double_t tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    Double_t y3=y0;
    g3->SetPoint(i,x1,y3);

    Double_t binwl=(x1-x3)/2.;
    Double_t binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,g0->GetErrorXlow(i),g0->GetErrorXhigh(i),(fabs(y3-y2)+fabs(y1-y3))/2,(fabs(y3-y2)+fabs(y1-y3))/2);

  }
  return g3;
}

void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {

  Double_t  x1=0., y1=0.,  y2=0., y0=0;
  
  if (g1->GetN()!=g2->GetN())
    cout << " graphs have not the same # of elements " <<g1->GetN()<<" "<<g2->GetN()<< endl;
  Double_t* EYhigh1 = g1-> GetEYhigh();
  Double_t* EYlow1  = g1-> GetEYlow();
  Double_t* EYhigh2 = g2-> GetEYhigh();
  Double_t* EYlow2  = g2-> GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (y1==0) y1=1;
    if (y2==0) y2=1;

    Double_t eyh1=0., eyl1=0.,eyh2=0., eyl2=0.;
  
    eyh1=EYhigh1[i];
    eyh2=EYhigh2[i];
    eyh2=sqrt(eyh1*eyh1+eyh2*eyh2);
    g2->SetPointEYhigh(i,eyh2);
    eyl1=EYlow1[i];
    eyl2=EYlow2[i];
    eyl2=sqrt(eyl1*eyl1+eyl2*eyl2);
    g2->SetPointEYlow (i,eyl2);
  }
  return;

}

TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

  const Int_t debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  Double_t   x1=0.,   y1=0., x2=0., y2=0.;
  Double_t dx1h=0., dx1l=0.;
  Double_t dy1h=0., dy1l=0.;
  Double_t dy2h=0., dy2l=0.;

  Double_t* X1 = g1->GetX();
  Double_t* Y1 = g1->GetY();
  Double_t* EXhigh1 = g1->GetEXhigh();
  Double_t* EXlow1 =  g1->GetEXlow();
  Double_t* EYhigh1 = g1->GetEYhigh();
  Double_t* EYlow1 =  g1->GetEYlow();

  Double_t* X2 = g2->GetX();
  Double_t* Y2 = g2->GetY();
  Double_t* EXhigh2 = g2->GetEXhigh();
  Double_t* EXlow2 =  g2->GetEXlow();
  Double_t* EYhigh2 = g2->GetEYhigh();
  Double_t* EYlow2 =  g2->GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
    
    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    Double_t el=0.; Double_t eh=0.;

    if (y2!=0.) el=sqrt(dy1l*dy1l)*(y1/y2);
    if (y2!=0.) eh=sqrt(dy1h*dy1h)*(y1/y2);

    g3->SetPointError(i,dx1h,dx1l,fabs(el),fabs(eh));

  }  
  return g3;

}
