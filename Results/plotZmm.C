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
#include <TGraphErrors.h>   
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

void plotZmm(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmm");
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

  TFile* fileMcAtNlo=new TFile("../TheoryUncertainty/Zmumu/zmm_PDFUnc.root", "OPEN");
  TFile* filePowheg=new TFile("../UnfoldingInput/Zmumu/zmmph_UnfoldInputs.root", "OPEN");
  TFile* fileMadgraph=new TFile("../UnfoldingInput/Zmumu/zmmmg_UnfoldInputs.root", "OPEN");



  TFile* fileFEWZ=new TFile("FEWZ_Histograms.root", "OPEN");


  // plot output file format
  const TString format("png");

   
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  

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
  TH1D * hUnfoldZPt;
  hUnfoldZPt=(TH1D*)(file[0]->Get("hUnfold"));

  TH1D * hTruthZPtMadgraph;
  TH1D * hTruthZPtPowheg;
  hTruthZPtMadgraph=(TH1D*)(fileMadgraph->Get("hZPtTruth"));
  hTruthZPtPowheg=(TH1D*)(filePowheg->Get("hZPtTruth"));

  for(int j=0;j!=nBinsZPt;++j)
    {
      hTruthZPtMadgraph->SetBinContent(j+1,hTruthZPtMadgraph->GetBinContent(j+1)/hTruthZPtMadgraph->GetBinWidth(j+1));
      hTruthZPtMadgraph->SetBinError(j+1,hTruthZPtMadgraph->GetBinError(j+1)/hTruthZPtMadgraph->GetBinWidth(j+1));
      hTruthZPtPowheg->SetBinContent(j+1,hTruthZPtPowheg->GetBinContent(j+1)/hTruthZPtPowheg->GetBinWidth(j+1));
      hTruthZPtPowheg->SetBinError(j+1,hTruthZPtPowheg->GetBinError(j+1)/hTruthZPtPowheg->GetBinWidth(j+1));
    }

  hTruthZPtMadgraph->Scale(1/lumi);
  hTruthZPtPowheg->Scale(1/lumi);



  cout<<hUnfoldZPt->Integral("width")<<endl;

  TH1D * hUnfoldZPtLumiUp;
  TH1D * hUnfoldZPtLumiDown;
  hUnfoldZPtLumiUp=(TH1D*)(fileLumiUp[0]->Get("hUnfold"));
  hUnfoldZPtLumiDown=(TH1D*)(fileLumiDown[0]->Get("hUnfold"));

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

  TH1D * hUnfoldZPtResScaleSys;
  hUnfoldZPtResScaleSys=(TH1D*)(fileResScaleSys[0]->Get("hUnfold"));

  TH1D * hUnfoldZPtUnfoldModelSysUp;
  TH1D * hUnfoldZPtUnfoldModelSysDown;
  hUnfoldZPtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[0]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldZPtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[0]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldZPtUnfoldModelSysDown->Add(hUnfoldZPt,-1.);
  hUnfoldZPtUnfoldModelSysDown->Scale(-1.);
  hUnfoldZPtUnfoldModelSysDown->Add(hUnfoldZPt,1.);

  TH1D * hUnfoldZPtUnfoldMatrixSys;
  hUnfoldZPtUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[0]->Get("hUnfold"));

  TGraphAsymmErrors* gUnfoldZPt=TH1TOTGraphAsymmErrors(hUnfoldZPt);
  TGraphAsymmErrors* gTruthZPtMadgraph=TH1TOTGraphAsymmErrors(hTruthZPtMadgraph);
  TGraphAsymmErrors* gTruthZPtPowheg=TH1TOTGraphAsymmErrors(hTruthZPtPowheg);


  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPt);
  TGraphAsymmErrors* ZPT_LUMI_UNCERT_BAND_DATA;
  ZPT_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtLumiUp),TH1TOTGraphAsymmErrors(hUnfoldZPtLumiDown));

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

  TGraphAsymmErrors* ZPT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPtResScaleSys);

  TGraphAsymmErrors* ZPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  ZPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysDown));

  TGraphAsymmErrors* ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldMatrixSys);

  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPt);

  myAddtoBand(ZPT_LUMI_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EWKBKG_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_TOPBKG_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSTAT_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFBIN_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_RESSCALE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMODEL_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_TOT_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthZPtMadgraph,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthZPtPowheg,ZPT_STAT_UNCERT_BAND_DATA);


  TH1D * hTruthZPtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hZPtTruthNominal"));
  TH1D * hTruthZPtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hZPtTruthPDFUp"));
  TH1D * hTruthZPtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hZPtTruthPDFDown"));
  TH1D * hTruthZPtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hZPtTruthScaleUp"));
  TH1D * hTruthZPtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hZPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthZPtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthZPtMcAtNlo);

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZPtMcAtNlo);
  TGraphAsymmErrors* ZPT_PDF_UNCERT_BAND_AMCATNLO;
  ZPT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloPDFDown));

  TGraphAsymmErrors* ZPT_SCALE_UNCERT_BAND_AMCATNLO;
  ZPT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloScaleDown));

  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZPtMcAtNlo);

  myAddtoBand(ZPT_PDF_UNCERT_BAND_AMCATNLO,ZPT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(ZPT_SCALE_UNCERT_BAND_AMCATNLO,ZPT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_AMCATNLO,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_TOT_UNCERT_BAND_AMCATNLO,ZPT_STAT_UNCERT_BAND_DATA);

  /*TGraphAsymmErrors* gTruthZPtFEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("Graph"));

    TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("Graph"));*/

  TGraphErrors* gTruthZPtFEWZ=(TGraphErrors*)fileFEWZ->Get("Graph");
  const int n=22;
  double* x=gTruthZPtFEWZ->GetX();
  double* y=gTruthZPtFEWZ->GetY();
  double* exl=gTruthZPtFEWZ->GetEX();
  double* exh=gTruthZPtFEWZ->GetEX();
  double* eyl=gTruthZPtFEWZ->GetEY();
  double* eyh=gTruthZPtFEWZ->GetEY();

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_FEWZ = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_FEWZ,ZPT_STAT_UNCERT_BAND_DATA);

 
  TH1D * hUnfoldPhiStar;
  hUnfoldPhiStar=(TH1D*)(file[1]->Get("hUnfold"));

  TH1D * hTruthPhiStarMadgraph;
  TH1D * hTruthPhiStarPowheg;
  hTruthPhiStarMadgraph=(TH1D*)(fileMadgraph->Get("hPhiStarTruth"));
  hTruthPhiStarPowheg=(TH1D*)(filePowheg->Get("hPhiStarTruth"));

  for(int j=0;j!=nBinsPhiStar;++j)
    {
      hTruthPhiStarMadgraph->SetBinContent(j+1,hTruthPhiStarMadgraph->GetBinContent(j+1)/hTruthPhiStarMadgraph->GetBinWidth(j+1));
      hTruthPhiStarMadgraph->SetBinError(j+1,hTruthPhiStarMadgraph->GetBinError(j+1)/hTruthPhiStarMadgraph->GetBinWidth(j+1));
      hTruthPhiStarPowheg->SetBinContent(j+1,hTruthPhiStarPowheg->GetBinContent(j+1)/hTruthPhiStarPowheg->GetBinWidth(j+1));
      hTruthPhiStarPowheg->SetBinError(j+1,hTruthPhiStarPowheg->GetBinError(j+1)/hTruthPhiStarPowheg->GetBinWidth(j+1));
    }

  hTruthPhiStarMadgraph->Scale(1/lumi);
  hTruthPhiStarPowheg->Scale(1/lumi);

  cout<<hUnfoldPhiStar->Integral("width")<<endl;

  TH1D * hUnfoldPhiStarLumiUp;
  TH1D * hUnfoldPhiStarLumiDown;
  hUnfoldPhiStarLumiUp=(TH1D*)(fileLumiUp[1]->Get("hUnfold"));
  hUnfoldPhiStarLumiDown=(TH1D*)(fileLumiDown[1]->Get("hUnfold"));

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

  TH1D * hUnfoldPhiStarResScaleSys;
  hUnfoldPhiStarResScaleSys=(TH1D*)(fileResScaleSys[1]->Get("hUnfold"));

  TH1D * hUnfoldPhiStarUnfoldModelSysUp;
  TH1D * hUnfoldPhiStarUnfoldModelSysDown;
  hUnfoldPhiStarUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[1]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldPhiStarUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[1]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldPhiStarUnfoldModelSysDown->Add(hUnfoldPhiStar,-1.);
  hUnfoldPhiStarUnfoldModelSysDown->Scale(-1.);
  hUnfoldPhiStarUnfoldModelSysDown->Add(hUnfoldPhiStar,1.);
  
  TH1D * hUnfoldPhiStarUnfoldMatrixSys;
  hUnfoldPhiStarUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[1]->Get("hUnfold"));

  TGraphAsymmErrors* gUnfoldPhiStar=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);
  TGraphAsymmErrors* gTruthPhiStarMadgraph=TH1TOTGraphAsymmErrors(hTruthPhiStarMadgraph);
  TGraphAsymmErrors* gTruthPhiStarPowheg=TH1TOTGraphAsymmErrors(hTruthPhiStarPowheg);


  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);
  TGraphAsymmErrors* PHISTAR_LUMI_UNCERT_BAND_DATA;
  PHISTAR_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarLumiUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarLumiDown));

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

  TGraphAsymmErrors* PHISTAR_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStarResScaleSys);

  TGraphAsymmErrors* PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA;
  PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysDown));

  TGraphAsymmErrors* PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldMatrixSys);

  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);

  myAddtoBand(PHISTAR_LUMI_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EWKBKG_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_TOPBKG_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSTAT_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFBIN_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_RESSCALE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_STAT_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_TOT_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthPhiStarMadgraph,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthPhiStarPowheg,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthPhiStarMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthNominal"));
  TH1D * hTruthPhiStarMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthPDFUp"));
  TH1D * hTruthPhiStarMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthPDFDown"));
  TH1D * hTruthPhiStarMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthScaleUp"));
  TH1D * hTruthPhiStarMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthScaleDown"));

  TGraphAsymmErrors* gTruthPhiStarMcAtNlo=TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNlo);

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNlo);
  TGraphAsymmErrors* PHISTAR_PDF_UNCERT_BAND_AMCATNLO;
  PHISTAR_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthPhiStarMcAtNlo,TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloPDFDown));

  TGraphAsymmErrors* PHISTAR_SCALE_UNCERT_BAND_AMCATNLO;
  PHISTAR_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthPhiStarMcAtNlo,TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloScaleDown));

  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNlo);

  myAddtoBand(PHISTAR_PDF_UNCERT_BAND_AMCATNLO,PHISTAR_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(PHISTAR_SCALE_UNCERT_BAND_AMCATNLO,PHISTAR_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_STAT_UNCERT_BAND_AMCATNLO,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_TOT_UNCERT_BAND_AMCATNLO,PHISTAR_STAT_UNCERT_BAND_DATA);



  TH1D * hUnfoldZRap;
  TH1D * hTruthZRapMadgraph;

  hUnfoldZRap=(TH1D*)(file[2]->Get("hUnfold"));
  hTruthZRapMadgraph=(TH1D*)(fileUnfoldModelSys[2]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldZRap->Integral("width")<<endl;

  TH1D * hUnfoldZRapLumiUp;
  TH1D * hUnfoldZRapLumiDown;
  hUnfoldZRapLumiUp=(TH1D*)(fileLumiUp[2]->Get("hUnfold"));
  hUnfoldZRapLumiDown=(TH1D*)(fileLumiDown[2]->Get("hUnfold"));

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

  TH1D * hUnfoldZRapResScaleSys;
  hUnfoldZRapResScaleSys=(TH1D*)(fileResScaleSys[2]->Get("hUnfold"));

  TH1D * hUnfoldZRapUnfoldModelSysUp;
  TH1D * hUnfoldZRapUnfoldModelSysDown;
  hUnfoldZRapUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[2]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldZRapUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[2]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldZRapUnfoldModelSysDown->Add(hUnfoldZRap,-1.);
  hUnfoldZRapUnfoldModelSysDown->Scale(-1.);
  hUnfoldZRapUnfoldModelSysDown->Add(hUnfoldZRap,1.);

  TH1D * hUnfoldZRapUnfoldMatrixSys;
  hUnfoldZRapUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[2]->Get("hUnfold"));
  
  TGraphAsymmErrors* gUnfoldZRap=TH1TOTGraphAsymmErrors(hUnfoldZRap);
  TGraphAsymmErrors* gTruthZRapMadgraph=TH1TOTGraphAsymmErrors(hTruthZRapMadgraph);

  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRap);
  TGraphAsymmErrors* ZRAP_LUMI_UNCERT_BAND_DATA;
  ZRAP_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapLumiUp),TH1TOTGraphAsymmErrors(hUnfoldZRapLumiDown));

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

  TGraphAsymmErrors* ZRAP_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRapResScaleSys);

  TGraphAsymmErrors* ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA;
  ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysDown));

  TGraphAsymmErrors* ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldMatrixSys);

  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRap);

  
  myAddtoBand(ZRAP_LUMI_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EWKBKG_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_TOPBKG_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSTAT_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFBIN_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_RESSCALE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_STAT_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_TOT_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthZRapMadgraph,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthZRapMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hZRapTruthNominal"));
  TH1D * hTruthZRapMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hZRapTruthPDFUp"));
  TH1D * hTruthZRapMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hZRapTruthPDFDown"));
  TH1D * hTruthZRapMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hZRapTruthScaleUp"));
  TH1D * hTruthZRapMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hZRapTruthScaleDown"));

  TGraphAsymmErrors* gTruthZRapMcAtNlo=TH1TOTGraphAsymmErrors(hTruthZRapMcAtNlo);

  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZRapMcAtNlo);
  TGraphAsymmErrors* ZRAP_PDF_UNCERT_BAND_AMCATNLO;
  ZRAP_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZRapMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloPDFDown));

  TGraphAsymmErrors* ZRAP_SCALE_UNCERT_BAND_AMCATNLO;
  ZRAP_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZRapMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloScaleDown));

  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZRapMcAtNlo);

  myAddtoBand(ZRAP_PDF_UNCERT_BAND_AMCATNLO,ZRAP_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(ZRAP_SCALE_UNCERT_BAND_AMCATNLO,ZRAP_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_STAT_UNCERT_BAND_AMCATNLO,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_TOT_UNCERT_BAND_AMCATNLO,ZRAP_STAT_UNCERT_BAND_DATA);


  //---------------------------------------------------------------------------
  //                             Lep1 Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep1Pt;
  TH1D * hTruthLep1PtMadgraph;
  
  hUnfoldLep1Pt=(TH1D*)(file[3]->Get("hUnfold"));
  hTruthLep1PtMadgraph=(TH1D*)(fileUnfoldModelSys[3]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldLep1Pt->Integral("width")<<endl;

  TH1D * hUnfoldLep1PtLumiUp;
  TH1D * hUnfoldLep1PtLumiDown;
  hUnfoldLep1PtLumiUp=(TH1D*)(fileLumiUp[3]->Get("hUnfold"));
  hUnfoldLep1PtLumiDown=(TH1D*)(fileLumiDown[3]->Get("hUnfold"));

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

  TH1D * hUnfoldLep1PtResScaleSys;
  hUnfoldLep1PtResScaleSys=(TH1D*)(fileResScaleSys[3]->Get("hUnfold"));

  TH1D * hUnfoldLep1PtUnfoldModelSysUp;
  TH1D * hUnfoldLep1PtUnfoldModelSysDown;
  hUnfoldLep1PtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[3]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep1PtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[3]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep1PtUnfoldModelSysDown->Add(hUnfoldLep1Pt,-1.);
  hUnfoldLep1PtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep1PtUnfoldModelSysDown->Add(hUnfoldLep1Pt,1.);

TH1D * hUnfoldLep1PtUnfoldMatrixSys;
  hUnfoldLep1PtUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[3]->Get("hUnfold"));
 
  TGraphAsymmErrors* gUnfoldLep1Pt=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);
  TGraphAsymmErrors* gTruthLep1PtMadgraph=TH1TOTGraphAsymmErrors(hTruthLep1PtMadgraph);

  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);
  TGraphAsymmErrors* LEP1PT_LUMI_UNCERT_BAND_DATA;
  LEP1PT_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtLumiUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtLumiDown));

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

  TGraphAsymmErrors* LEP1PT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1PtResScaleSys);

  TGraphAsymmErrors* LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysDown));

  TGraphAsymmErrors* LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldMatrixSys);

  TGraphAsymmErrors* LEP1PT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);

  myAddtoBand(LEP1PT_LUMI_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EWKBKG_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_TOPBKG_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSTAT_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFBIN_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_RESSCALE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_STAT_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_TOT_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthLep1PtMadgraph,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep1PtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthNominal"));
  TH1D * hTruthLep1PtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthPDFUp"));
  TH1D * hTruthLep1PtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthPDFDown"));
  TH1D * hTruthLep1PtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthScaleUp"));
  TH1D * hTruthLep1PtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep1PtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNlo);

  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNlo);
  TGraphAsymmErrors* LEP1PT_PDF_UNCERT_BAND_AMCATNLO;
  LEP1PT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloPDFDown));

  TGraphAsymmErrors* LEP1PT_SCALE_UNCERT_BAND_AMCATNLO;
  LEP1PT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloScaleDown));

  TGraphAsymmErrors* LEP1PT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNlo);

  myAddtoBand(LEP1PT_PDF_UNCERT_BAND_AMCATNLO,LEP1PT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP1PT_SCALE_UNCERT_BAND_AMCATNLO,LEP1PT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_STAT_UNCERT_BAND_AMCATNLO,LEP1PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_TOT_UNCERT_BAND_AMCATNLO,LEP1PT_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             Lep2 Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep2Pt;
  TH1D * hTruthLep2PtMadgraph;
  
  hUnfoldLep2Pt=(TH1D*)(file[4]->Get("hUnfold"));
  hTruthLep2PtMadgraph=(TH1D*)(fileUnfoldModelSys[4]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldLep2Pt->Integral("width")<<endl;

  TH1D * hUnfoldLep2PtLumiUp;
  TH1D * hUnfoldLep2PtLumiDown;
  hUnfoldLep2PtLumiUp=(TH1D*)(fileLumiUp[4]->Get("hUnfold"));
  hUnfoldLep2PtLumiDown=(TH1D*)(fileLumiDown[4]->Get("hUnfold"));

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

  TH1D * hUnfoldLep2PtResScaleSys;
  hUnfoldLep2PtResScaleSys=(TH1D*)(fileResScaleSys[4]->Get("hUnfold"));

  TH1D * hUnfoldLep2PtUnfoldModelSysUp;
  TH1D * hUnfoldLep2PtUnfoldModelSysDown;
  hUnfoldLep2PtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[4]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep2PtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[4]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep2PtUnfoldModelSysDown->Add(hUnfoldLep2Pt,-1.);
  hUnfoldLep2PtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep2PtUnfoldModelSysDown->Add(hUnfoldLep2Pt,1.);

  TH1D * hUnfoldLep2PtUnfoldMatrixSys;
  hUnfoldLep2PtUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[4]->Get("hUnfold"));

  TGraphAsymmErrors* gUnfoldLep2Pt=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);
  TGraphAsymmErrors* gTruthLep2PtMadgraph=TH1TOTGraphAsymmErrors(hTruthLep2PtMadgraph);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);
  TGraphAsymmErrors* LEP2PT_LUMI_UNCERT_BAND_DATA;
  LEP2PT_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtLumiUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtLumiDown));

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

  TGraphAsymmErrors* LEP2PT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2PtResScaleSys);

  TGraphAsymmErrors* LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysDown));

  TGraphAsymmErrors* LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldMatrixSys);

  TGraphAsymmErrors* LEP2PT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);

  myAddtoBand(LEP2PT_LUMI_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EWKBKG_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_TOPBKG_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSTAT_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFBIN_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_RESSCALE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_STAT_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_TOT_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthLep2PtMadgraph,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep2PtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthNominal"));
  TH1D * hTruthLep2PtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthPDFUp"));
  TH1D * hTruthLep2PtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthPDFDown"));
  TH1D * hTruthLep2PtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthScaleUp"));
  TH1D * hTruthLep2PtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep2PtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNlo);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNlo);
  TGraphAsymmErrors* LEP2PT_PDF_UNCERT_BAND_AMCATNLO;
  LEP2PT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloPDFDown));

  TGraphAsymmErrors* LEP2PT_SCALE_UNCERT_BAND_AMCATNLO;
  LEP2PT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloScaleDown));

  TGraphAsymmErrors* LEP2PT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNlo);

  myAddtoBand(LEP2PT_PDF_UNCERT_BAND_AMCATNLO,LEP2PT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP2PT_SCALE_UNCERT_BAND_AMCATNLO,LEP2PT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_STAT_UNCERT_BAND_AMCATNLO,LEP2PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_TOT_UNCERT_BAND_AMCATNLO,LEP2PT_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             Lep1 Eta
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep1Eta;
  TH1D * hTruthLep1EtaMadgraph;
  
  hUnfoldLep1Eta=(TH1D*)(file[5]->Get("hUnfold"));
  hTruthLep1EtaMadgraph=(TH1D*)(fileUnfoldModelSys[5]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldLep1Eta->Integral("width")<<endl;

  TH1D * hUnfoldLep1EtaLumiUp;
  TH1D * hUnfoldLep1EtaLumiDown;
  hUnfoldLep1EtaLumiUp=(TH1D*)(fileLumiUp[5]->Get("hUnfold"));
  hUnfoldLep1EtaLumiDown=(TH1D*)(fileLumiDown[5]->Get("hUnfold"));

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

  TH1D * hUnfoldLep1EtaResScaleSys;
  hUnfoldLep1EtaResScaleSys=(TH1D*)(fileResScaleSys[5]->Get("hUnfold"));

  TH1D * hUnfoldLep1EtaUnfoldModelSysUp;
  TH1D * hUnfoldLep1EtaUnfoldModelSysDown;
  hUnfoldLep1EtaUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[5]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep1EtaUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[5]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep1EtaUnfoldModelSysDown->Add(hUnfoldLep1Eta,-1.);
  hUnfoldLep1EtaUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep1EtaUnfoldModelSysDown->Add(hUnfoldLep1Eta,1.);

  TH1D * hUnfoldLep1EtaUnfoldMatrixSys;
  hUnfoldLep1EtaUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[5]->Get("hUnfold"));
 
  TGraphAsymmErrors* gUnfoldLep1Eta=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);
  TGraphAsymmErrors* gTruthLep1EtaMadgraph=TH1TOTGraphAsymmErrors(hTruthLep1EtaMadgraph);

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);
  TGraphAsymmErrors* LEP1ETA_LUMI_UNCERT_BAND_DATA;
  LEP1ETA_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaLumiUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaLumiDown));

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

  TGraphAsymmErrors* LEP1ETA_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1EtaResScaleSys);

  TGraphAsymmErrors* LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysDown));

  TGraphAsymmErrors* LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldMatrixSys);

  TGraphAsymmErrors* LEP1ETA_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);

  myAddtoBand(LEP1ETA_LUMI_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EWKBKG_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_TOPBKG_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSTAT_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFBIN_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_RESSCALE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_STAT_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_TOT_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthLep1EtaMadgraph,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep1EtaMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthNominal"));
  TH1D * hTruthLep1EtaMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthPDFUp"));
  TH1D * hTruthLep1EtaMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthPDFDown"));
  TH1D * hTruthLep1EtaMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthScaleUp"));
  TH1D * hTruthLep1EtaMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep1EtaMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNlo);

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNlo);
  TGraphAsymmErrors* LEP1ETA_PDF_UNCERT_BAND_AMCATNLO;
  LEP1ETA_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloPDFDown));

  TGraphAsymmErrors* LEP1ETA_SCALE_UNCERT_BAND_AMCATNLO;
  LEP1ETA_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloScaleDown));

  TGraphAsymmErrors* LEP1ETA_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNlo);

  myAddtoBand(LEP1ETA_PDF_UNCERT_BAND_AMCATNLO,LEP1ETA_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP1ETA_SCALE_UNCERT_BAND_AMCATNLO,LEP1ETA_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_STAT_UNCERT_BAND_AMCATNLO,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_TOT_UNCERT_BAND_AMCATNLO,LEP1ETA_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             Lep2 Eta
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep2Eta;
  TH1D * hTruthLep2EtaMadgraph;
  
  hUnfoldLep2Eta=(TH1D*)(file[6]->Get("hUnfold"));
  hTruthLep2EtaMadgraph=(TH1D*)(fileUnfoldModelSys[6]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldLep2Eta->Integral("width")<<endl;

  TH1D * hUnfoldLep2EtaLumiUp;
  TH1D * hUnfoldLep2EtaLumiDown;
  hUnfoldLep2EtaLumiUp=(TH1D*)(fileLumiUp[6]->Get("hUnfold"));
  hUnfoldLep2EtaLumiDown=(TH1D*)(fileLumiDown[6]->Get("hUnfold"));

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

  TH1D * hUnfoldLep2EtaResScaleSys;
  hUnfoldLep2EtaResScaleSys=(TH1D*)(fileResScaleSys[6]->Get("hUnfold"));

  TH1D * hUnfoldLep2EtaUnfoldModelSysUp;
  TH1D * hUnfoldLep2EtaUnfoldModelSysDown;
  hUnfoldLep2EtaUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[6]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLep2EtaUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[6]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLep2EtaUnfoldModelSysDown->Add(hUnfoldLep2Eta,-1.);
  hUnfoldLep2EtaUnfoldModelSysDown->Scale(-1.);
  hUnfoldLep2EtaUnfoldModelSysDown->Add(hUnfoldLep2Eta,1.);

  TH1D * hUnfoldLep2EtaUnfoldMatrixSys;
  hUnfoldLep2EtaUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[6]->Get("hUnfold"));

  TGraphAsymmErrors* gUnfoldLep2Eta=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);
  TGraphAsymmErrors* gTruthLep2EtaMadgraph=TH1TOTGraphAsymmErrors(hTruthLep2EtaMadgraph);

  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);
  TGraphAsymmErrors* LEP2ETA_LUMI_UNCERT_BAND_DATA;
  LEP2ETA_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaLumiUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaLumiDown));

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

  TGraphAsymmErrors* LEP2ETA_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2EtaResScaleSys);

  TGraphAsymmErrors* LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysDown));

  TGraphAsymmErrors* LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldMatrixSys);

  TGraphAsymmErrors* LEP2ETA_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);

  myAddtoBand(LEP2ETA_LUMI_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EWKBKG_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_TOPBKG_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSTAT_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFBIN_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_RESSCALE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_STAT_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_TOT_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthLep2EtaMadgraph,LEP2ETA_STAT_UNCERT_BAND_DATA);
  
TH1D * hTruthLep2EtaMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthNominal"));
  TH1D * hTruthLep2EtaMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthPDFUp"));
  TH1D * hTruthLep2EtaMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthPDFDown"));
  TH1D * hTruthLep2EtaMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthScaleUp"));
  TH1D * hTruthLep2EtaMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep2EtaMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNlo);

  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNlo);
  TGraphAsymmErrors* LEP2ETA_PDF_UNCERT_BAND_AMCATNLO;
  LEP2ETA_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloPDFDown));

  TGraphAsymmErrors* LEP2ETA_SCALE_UNCERT_BAND_AMCATNLO;
  LEP2ETA_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloScaleDown));

  TGraphAsymmErrors* LEP2ETA_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNlo);

  myAddtoBand(LEP2ETA_PDF_UNCERT_BAND_AMCATNLO,LEP2ETA_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP2ETA_SCALE_UNCERT_BAND_AMCATNLO,LEP2ETA_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_STAT_UNCERT_BAND_AMCATNLO,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_TOT_UNCERT_BAND_AMCATNLO,LEP2ETA_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             LepNeg Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLepNegPt;
  TH1D * hTruthLepNegPtMadgraph;
  
  hUnfoldLepNegPt=(TH1D*)(file[7]->Get("hUnfold"));
  hTruthLepNegPtMadgraph=(TH1D*)(fileUnfoldModelSys[7]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldLepNegPt->Integral("width")<<endl;

  TH1D * hUnfoldLepNegPtLumiUp;
  TH1D * hUnfoldLepNegPtLumiDown;
  hUnfoldLepNegPtLumiUp=(TH1D*)(fileLumiUp[7]->Get("hUnfold"));
  hUnfoldLepNegPtLumiDown=(TH1D*)(fileLumiDown[7]->Get("hUnfold"));

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

  TH1D * hUnfoldLepNegPtResScaleSys;
  hUnfoldLepNegPtResScaleSys=(TH1D*)(fileResScaleSys[7]->Get("hUnfold"));

  TH1D * hUnfoldLepNegPtUnfoldModelSysUp;
  TH1D * hUnfoldLepNegPtUnfoldModelSysDown;
  hUnfoldLepNegPtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[7]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLepNegPtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[7]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLepNegPtUnfoldModelSysDown->Add(hUnfoldLepNegPt,-1.);
  hUnfoldLepNegPtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLepNegPtUnfoldModelSysDown->Add(hUnfoldLepNegPt,1.);

  TH1D * hUnfoldLepNegPtUnfoldMatrixSys;
  hUnfoldLepNegPtUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[7]->Get("hUnfold"));

  TGraphAsymmErrors* gUnfoldLepNegPt=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);
  TGraphAsymmErrors* gTruthLepNegPtMadgraph=TH1TOTGraphAsymmErrors(hTruthLepNegPtMadgraph);

  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);
  TGraphAsymmErrors* LEPNEGPT_LUMI_UNCERT_BAND_DATA;
  LEPNEGPT_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtLumiUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtLumiDown));

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

  TGraphAsymmErrors* LEPNEGPT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPtResScaleSys);

  TGraphAsymmErrors* LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysDown));

  TGraphAsymmErrors* LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldMatrixSys);

  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);

  myAddtoBand(LEPNEGPT_LUMI_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EWKBKG_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_TOPBKG_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFBIN_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_RESSCALE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_STAT_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_TOT_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthLepNegPtMadgraph,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLepNegPtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthNominal"));
  TH1D * hTruthLepNegPtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthPDFUp"));
  TH1D * hTruthLepNegPtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthPDFDown"));
  TH1D * hTruthLepNegPtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthScaleUp"));
  TH1D * hTruthLepNegPtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLepNegPtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNlo);

  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNlo);
  TGraphAsymmErrors* LEPNEGPT_PDF_UNCERT_BAND_AMCATNLO;
  LEPNEGPT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepNegPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloPDFDown));

  TGraphAsymmErrors* LEPNEGPT_SCALE_UNCERT_BAND_AMCATNLO;
  LEPNEGPT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepNegPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloScaleDown));

  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNlo);

  myAddtoBand(LEPNEGPT_PDF_UNCERT_BAND_AMCATNLO,LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEPNEGPT_SCALE_UNCERT_BAND_AMCATNLO,LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             LepPos Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLepPosPt;
  TH1D * hTruthLepPosPtMadgraph;
  
  hUnfoldLepPosPt=(TH1D*)(file[8]->Get("hUnfold"));
  hTruthLepPosPtMadgraph=(TH1D*)(fileUnfoldModelSys[8]->Get("UNFOLDMODEL/hTruth"));

  cout<<hUnfoldLepPosPt->Integral("width")<<endl;

  TH1D * hUnfoldLepPosPtLumiUp;
  TH1D * hUnfoldLepPosPtLumiDown;
  hUnfoldLepPosPtLumiUp=(TH1D*)(fileLumiUp[8]->Get("hUnfold"));
  hUnfoldLepPosPtLumiDown=(TH1D*)(fileLumiDown[8]->Get("hUnfold"));

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

  TH1D * hUnfoldLepPosPtResScaleSys;
  hUnfoldLepPosPtResScaleSys=(TH1D*)(fileResScaleSys[8]->Get("hUnfold"));

  TH1D * hUnfoldLepPosPtUnfoldModelSysUp;
  TH1D * hUnfoldLepPosPtUnfoldModelSysDown;
  hUnfoldLepPosPtUnfoldModelSysUp=(TH1D*)(fileUnfoldModelSys[8]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));
  hUnfoldLepPosPtUnfoldModelSysDown=(TH1D*)(fileUnfoldModelSys[8]->Get("SMOOTH_UNFOLDMODEL/hUnfold"));

  hUnfoldLepPosPtUnfoldModelSysDown->Add(hUnfoldLepPosPt,-1.);
  hUnfoldLepPosPtUnfoldModelSysDown->Scale(-1.);
  hUnfoldLepPosPtUnfoldModelSysDown->Add(hUnfoldLepPosPt,1.);

  TH1D * hUnfoldLepPosPtUnfoldMatrixSys;
  hUnfoldLepPosPtUnfoldMatrixSys=(TH1D*)(fileUnfoldMatrixSys[8]->Get("hUnfold"));

  TGraphAsymmErrors* gUnfoldLepPosPt=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);
  TGraphAsymmErrors* gTruthLepPosPtMadgraph=TH1TOTGraphAsymmErrors(hTruthLepPosPtMadgraph);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);
  TGraphAsymmErrors* LEPPOSPT_LUMI_UNCERT_BAND_DATA;
  LEPPOSPT_LUMI_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtLumiUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtLumiDown));

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

  TGraphAsymmErrors* LEPPOSPT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPtResScaleSys);

  TGraphAsymmErrors* LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysDown));

  TGraphAsymmErrors* LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldMatrixSys);

  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);

  myAddtoBand(LEPPOSPT_LUMI_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EWKBKG_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_TOPBKG_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFBIN_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_RESSCALE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_STAT_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_TOT_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP=myTGraphErrorsDivide_noErrGraph2(gTruthLepPosPtMadgraph,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLepPosPtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthNominal"));
  TH1D * hTruthLepPosPtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthPDFUp"));
  TH1D * hTruthLepPosPtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthPDFDown"));
  TH1D * hTruthLepPosPtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthScaleUp"));
  TH1D * hTruthLepPosPtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLepPosPtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNlo);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNlo);
  TGraphAsymmErrors* LEPPOSPT_PDF_UNCERT_BAND_AMCATNLO;
  LEPPOSPT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepPosPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloPDFDown));

  TGraphAsymmErrors* LEPPOSPT_SCALE_UNCERT_BAND_AMCATNLO;
  LEPPOSPT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepPosPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloScaleDown));

  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNlo);

  myAddtoBand(LEPPOSPT_PDF_UNCERT_BAND_AMCATNLO,LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEPPOSPT_SCALE_UNCERT_BAND_AMCATNLO,LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO,LEPPOSPT_STAT_UNCERT_BAND_DATA);
 
 
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  // plot colors
  Int_t linecolorZ   = kOrange-3;
  Int_t fillcolorZ   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorTop = kGreen+2;
  Int_t fillcolorTop = kGreen-5;
  Int_t ratioColor   = kGray+2;

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
  TGaxis::SetMaxDigits(3);
  
  //
  // ZPt
  //   
  TH1D *ZPT_HIST_DUMMY = new TH1D("ZPT_HIST_DUMMY", "ZPT_HIST_DUMMY",nBinsZPt,ZPtBins);

  TH1D *hZmumuPtDiffDummy = makeDiffHist(ZPT_HIST_DUMMY,ZPT_HIST_DUMMY,"hZmumuPtDiffDummy");
  
  sprintf(ylabel,"d#sigma/dp_{T}^{#mu^{+}#mu^{-}} [pb/GeV]");
  ZPT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  ZPT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  ZPT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  ZPT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuPt("zmmPt","","",ylabel);
  plotZmumuPt.AddHist1D(ZPT_HIST_DUMMY);
  plotZmumuPt.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuPt.AddGraph(ZPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuPt.AddGraph(ZPT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuPt.AddGraph(gTruthZPtPowheg,"POWHEG","P",kGreen+3,22,1);
  plotZmumuPt.AddGraph(gTruthZPtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuPt.AddGraph(ZPT_STAT_UNCERT_BAND_FEWZ,"FEWZ","P",kBlue,24,1);
  plotZmumuPt.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuPt.AddGraph(ZPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuPt.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuPt.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuPt.SetLogx();
  plotZmumuPt.SetLogy(0);
  plotZmumuPt.SetYRange(0.01,1.2*(hUnfoldZPt->GetMaximum() + sqrt(hUnfoldZPt->GetMaximum())));
  plotZmumuPt.TransLegend(0.05,-0.05);
  plotZmumuPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuPtDiff("zmmPt","","p_{T}^{#mu^{+}#mu^{-}} [GeV]","Pred/Data");
  plotZmumuPtDiff.AddHist1D(hZmumuPtDiffDummy);
  plotZmumuPtDiff.AddGraph(ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuPtDiff.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuPtDiff.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",kGreen+3,22,1);
  plotZmumuPtDiff.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuPtDiff.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP,"0P",kBlue,24,1);
  plotZmumuPtDiff.AddGraph(ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuPtDiff.AddGraph(ZPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuPtDiff.SetLogx();
  plotZmumuPtDiff.SetYRange(0.8,1.2);

  plotZmumuPtDiff.AddLine(0, 1,1000, 1,kBlack,1);
  plotZmumuPtDiff.AddLine(0, 1.1,1000, 1.1,kBlack,3);
  plotZmumuPtDiff.AddLine(0,0.9,1000,0.9,kBlack,3);

  plotZmumuPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuPt2("zmmPtlog","","",ylabel);
  plotZmumuPt2.AddHist1D(ZPT_HIST_DUMMY);
  plotZmumuPt2.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuPt2.AddGraph(ZPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuPt2.AddGraph(ZPT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuPt2.AddGraph(gTruthZPtPowheg,"POWHEG","P",kGreen+3,22,1);
  plotZmumuPt2.AddGraph(gTruthZPtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuPt2.AddGraph(ZPT_STAT_UNCERT_BAND_FEWZ,"FEWZ","P",kBlue,24,1);
  plotZmumuPt2.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuPt2.AddGraph(ZPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuPt2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuPt2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuPt2.SetLogx();
  plotZmumuPt2.SetLogy();
  plotZmumuPt2.SetYRange(1e-6*(hUnfoldZPt->GetMaximum()),50*(hUnfoldZPt->GetMaximum()));
  plotZmumuPt2.TransLegend(0.05,-0.05);
  plotZmumuPt2.Draw(c,kTRUE,format,1);

  //
  // PhiStar
  //   
  TH1D *PHISTAR_HIST_DUMMY = new TH1D("PHISTAR_HIST_DUMMY", "PHISTAR_HIST_DUMMY",nBinsPhiStar,PhiStarBins);

  TH1D *hZmumuPhiStarDiffDummy = makeDiffHist(PHISTAR_HIST_DUMMY,PHISTAR_HIST_DUMMY,"hZmumuPhiStarDiffDummy");
  
  sprintf(ylabel,"d#sigma/d#phi_{#eta}* [pb]");
  PHISTAR_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  PHISTAR_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  PHISTAR_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  PHISTAR_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuPhiStar("zmmPhiStar","","",ylabel);
  plotZmumuPhiStar.AddHist1D(PHISTAR_HIST_DUMMY);
  plotZmumuPhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuPhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuPhiStar.AddGraph(PHISTAR_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuPhiStar.AddGraph(gTruthPhiStarPowheg,"POWHEG","P",kGreen+3,22,1);
  plotZmumuPhiStar.AddGraph(gTruthPhiStarMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuPhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuPhiStar.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuPhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuPhiStar.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuPhiStar.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuPhiStar.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuPhiStar.SetLogx();
  plotZmumuPhiStar.SetLogy(0);
  plotZmumuPhiStar.SetYRange(0.01,1.2*(hUnfoldPhiStar->GetMaximum() + sqrt(hUnfoldPhiStar->GetMaximum())));
  plotZmumuPhiStar.TransLegend(0.05,-0.05);
  plotZmumuPhiStar.Draw(c,kFALSE,format,1);

  CPlot plotZmumuPhiStarDiff("zmmPhiStar","","#phi_{#eta}*","Pred/Data");
  plotZmumuPhiStarDiff.AddHist1D(hZmumuPhiStarDiffDummy);
  plotZmumuPhiStarDiff.AddGraph(PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuPhiStarDiff.AddGraph(PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuPhiStarDiff.AddGraph(PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",kGreen+3,22,1);
  plotZmumuPhiStarDiff.AddGraph(PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuPhiStarDiff.AddGraph(PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuPhiStarDiff.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuPhiStarDiff.SetLogx();
  plotZmumuPhiStarDiff.SetYRange(0.8,1.2);

  plotZmumuPhiStarDiff.AddLine(0, 1,3, 1,kBlack,1);
  plotZmumuPhiStarDiff.AddLine(0, 1.1,3, 1.1,kBlack,3);
  plotZmumuPhiStarDiff.AddLine(0,0.9,3,0.9,kBlack,3);

  plotZmumuPhiStarDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuPhiStar2("zmmPhiStarlog","","",ylabel);
  plotZmumuPhiStar2.AddHist1D(PHISTAR_HIST_DUMMY);
  plotZmumuPhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuPhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuPhiStar2.AddGraph(PHISTAR_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuPhiStar2.AddGraph(gTruthPhiStarPowheg,"POWHEG","P",kGreen+3,22,1);
  plotZmumuPhiStar2.AddGraph(gTruthPhiStarMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuPhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuPhiStar2.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuPhiStar2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuPhiStar2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuPhiStar2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuPhiStar2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuPhiStar2.SetLogx();
  plotZmumuPhiStar2.SetLogy();
  plotZmumuPhiStar2.SetYRange(2e-4*(hUnfoldPhiStar->GetMaximum()),10*(hUnfoldPhiStar->GetMaximum()));
  plotZmumuPhiStar2.TransLegend(0.05,-0.05);
  plotZmumuPhiStar2.Draw(c,kTRUE,format,1);
  
  //
  // Z Rapidity
  //   
  TH1D *ZRAP_HIST_DUMMY = new TH1D("ZRAP_HIST_DUMMY", "ZRAP_HIST_DUMMY",24,0,2.4);

  TH1D *hZmumuRapDiffDummy = makeDiffHist(ZRAP_HIST_DUMMY,ZRAP_HIST_DUMMY,"hZmumuRapDiffDummy");
  
  sprintf(ylabel,"d#sigma/dy^{#mu^{+}#mu^{-}} [pb]");
  ZRAP_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  ZRAP_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  ZRAP_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  ZRAP_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuRap("zmmRap","","",ylabel);
  plotZmumuRap.AddHist1D(ZRAP_HIST_DUMMY);
  plotZmumuRap.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuRap.AddGraph(ZRAP_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuRap.AddGraph(ZRAP_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuRap.AddGraph(gTruthZRapMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuRap.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuRap.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuRap.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuRap.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuRap.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuRap.SetLogx(0);
  plotZmumuRap.SetLogy(0);
  plotZmumuRap.SetYRange(0.01,1.5*(hUnfoldZRap->GetMaximum() + sqrt(hUnfoldZRap->GetMaximum())));
  plotZmumuRap.TransLegend(0.05,-0.05);
  plotZmumuRap.Draw(c,kFALSE,format,1);

  CPlot plotZmumuRapDiff("zmmRap","","|y^{#mu^{+}#mu^{-}}|","Pred/Data");
  plotZmumuRapDiff.AddHist1D(hZmumuRapDiffDummy);
  plotZmumuRapDiff.AddGraph(ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuRapDiff.AddGraph(ZRAP_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuRapDiff.AddGraph(ZRAP_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuRapDiff.AddGraph(ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuRapDiff.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuRapDiff.SetLogx(0);
  plotZmumuRapDiff.SetYRange(0.8,1.2);

  plotZmumuRapDiff.AddLine(0, 1,2.4, 1,kBlack,1);
  plotZmumuRapDiff.AddLine(0, 1.1,2.4, 1.1,kBlack,3);
  plotZmumuRapDiff.AddLine(0,0.9,2.4,0.9,kBlack,3);

  plotZmumuRapDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuRap2("zmmRaplog","","",ylabel);
  plotZmumuRap2.AddHist1D(ZRAP_HIST_DUMMY);
  plotZmumuRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuRap2.AddGraph(ZRAP_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuRap2.AddGraph(gTruthZRapMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuRap2.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuRap2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuRap2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuRap2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuRap2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuRap2.SetLogx(0);
  plotZmumuRap2.SetLogy();
  plotZmumuRap2.SetYRange(4e-2*(hUnfoldZRap->GetMaximum()),5*(hUnfoldZRap->GetMaximum()));
  plotZmumuRap2.TransLegend(0.05,-0.05);
  plotZmumuRap2.Draw(c,kTRUE,format,1);
 
  //
  // Lep1 Pt
  //   
  TH1D *LEP1PT_HIST_DUMMY = new TH1D("LEP1PT_HIST_DUMMY", "LEP1PT_HIST_DUMMY",nBinsLep1Pt,Lep1PtBins);

  TH1D *hZmumuLep1PtDiffDummy = makeDiffHist(LEP1PT_HIST_DUMMY,LEP1PT_HIST_DUMMY,"hZmumuLep1PtDiffDummy");
  
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");
  LEP1PT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  LEP1PT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP1PT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP1PT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  LEP1PT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuLep1Pt("zmmLep1Pt","","",ylabel);
  plotZmumuLep1Pt.AddHist1D(LEP1PT_HIST_DUMMY);
  plotZmumuLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep1Pt.AddGraph(LEP1PT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep1Pt.AddGraph(gTruthLep1PtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep1Pt.AddGraph(LEP1PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep1Pt.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.65,0.1,0.87,0.15,0);
  plotZmumuLep1Pt.AddTextBox("PDF set: NNPDF3.0",0.65,0.05,0.87,0.1,0);
  plotZmumuLep1Pt.SetLogx();
  plotZmumuLep1Pt.SetLogy(0);
  plotZmumuLep1Pt.SetYRange(0.01,1.2*(hUnfoldLep1Pt->GetMaximum() + sqrt(hUnfoldLep1Pt->GetMaximum())));
  plotZmumuLep1Pt.TransLegend(0.05,-0.05);
  plotZmumuLep1Pt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep1PtDiff("zmmLep1Pt","","p_{T}(leading muon) [GeV]","Pred/Data");
  plotZmumuLep1PtDiff.AddHist1D(hZmumuLep1PtDiffDummy);
  plotZmumuLep1PtDiff.AddGraph(LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuLep1PtDiff.AddGraph(LEP1PT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuLep1PtDiff.AddGraph(LEP1PT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuLep1PtDiff.AddGraph(LEP1PT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuLep1PtDiff.AddGraph(LEP1PT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuLep1PtDiff.SetLogx();
  plotZmumuLep1PtDiff.SetYRange(0.8,1.2);

  plotZmumuLep1PtDiff.AddLine(0, 1,300, 1,kBlack,1);
  plotZmumuLep1PtDiff.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZmumuLep1PtDiff.AddLine(0,0.9,300,0.9,kBlack,3);

  plotZmumuLep1PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep1Pt2("zmmLep1Ptlog","","",ylabel);
  plotZmumuLep1Pt2.AddHist1D(LEP1PT_HIST_DUMMY);
  plotZmumuLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep1Pt2.AddGraph(LEP1PT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep1Pt2.AddGraph(gTruthLep1PtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep1Pt2.AddGraph(LEP1PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep1Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep1Pt2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLep1Pt2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLep1Pt2.SetLogx();
  plotZmumuLep1Pt2.SetLogy();
  plotZmumuLep1Pt2.SetYRange(5e-5*(hUnfoldLep1Pt->GetMaximum()),10*(hUnfoldLep1Pt->GetMaximum()));
  plotZmumuLep1Pt2.TransLegend(0.05,-0.05);
  plotZmumuLep1Pt2.Draw(c,kTRUE,format,1);

  //
  // Lep2 Pt
  //   
  TH1D *LEP2PT_HIST_DUMMY = new TH1D("LEP2PT_HIST_DUMMY", "LEP2PT_HIST_DUMMY",nBinsLep2Pt,Lep2PtBins);

  TH1D *hZmumuLep2PtDiffDummy = makeDiffHist(LEP2PT_HIST_DUMMY,LEP2PT_HIST_DUMMY,"hZmumuLep2PtDiffDummy");
  
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");
  LEP2PT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  LEP2PT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP2PT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP2PT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  LEP2PT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuLep2Pt("zmmLep2Pt","","",ylabel);
  plotZmumuLep2Pt.AddHist1D(LEP2PT_HIST_DUMMY);
  plotZmumuLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep2Pt.AddGraph(LEP2PT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep2Pt.AddGraph(gTruthLep2PtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep2Pt.AddGraph(LEP2PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Pt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep2Pt.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.65,0.1,0.87,0.15,0);
  plotZmumuLep2Pt.AddTextBox("PDF set: NNPDF3.0",0.65,0.05,0.87,0.1,0);
  plotZmumuLep2Pt.SetLogx();
  plotZmumuLep2Pt.SetLogy(0);
  plotZmumuLep2Pt.SetYRange(0.01,1.2*(hUnfoldLep2Pt->GetMaximum() + sqrt(hUnfoldLep2Pt->GetMaximum())));
  plotZmumuLep2Pt.TransLegend(0.05,-0.05);
  plotZmumuLep2Pt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep2PtDiff("zmmLep2Pt","","p_{T}(2nd leading muon) [GeV]","Pred/Data");
  plotZmumuLep2PtDiff.AddHist1D(hZmumuLep2PtDiffDummy);
  plotZmumuLep2PtDiff.AddGraph(LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuLep2PtDiff.AddGraph(LEP2PT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuLep2PtDiff.AddGraph(LEP2PT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuLep2PtDiff.AddGraph(LEP2PT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuLep2PtDiff.AddGraph(LEP2PT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuLep2PtDiff.SetLogx();
  plotZmumuLep2PtDiff.SetYRange(0.8,1.2);

  plotZmumuLep2PtDiff.AddLine(0, 1,200, 1,kBlack,1);
  plotZmumuLep2PtDiff.AddLine(0, 1.1,200, 1.1,kBlack,3);
  plotZmumuLep2PtDiff.AddLine(0,0.9,200,0.9,kBlack,3);

  plotZmumuLep2PtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep2Pt2("zmmLep2Ptlog","","",ylabel);
  plotZmumuLep2Pt2.AddHist1D(LEP2PT_HIST_DUMMY);
  plotZmumuLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep2Pt2.AddGraph(LEP2PT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep2Pt2.AddGraph(gTruthLep2PtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep2Pt2.AddGraph(LEP2PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep2Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Pt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep2Pt2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLep2Pt2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLep2Pt2.SetLogx();
  plotZmumuLep2Pt2.SetLogy();
  plotZmumuLep2Pt2.SetYRange(5e-5*(hUnfoldLep2Pt->GetMaximum()),10*(hUnfoldLep2Pt->GetMaximum()));
  plotZmumuLep2Pt2.TransLegend(0.05,-0.05);
  plotZmumuLep2Pt2.Draw(c,kTRUE,format,1);

  //
  // Lep1 Eta
  //   
  TH1D *LEP1ETA_HIST_DUMMY = new TH1D("LEP1ETA_HIST_DUMMY", "LEP1ETA_HIST_DUMMY",14,0,2.4);

  TH1D *hZmumuLep1EtaDiffDummy = makeDiffHist(LEP1ETA_HIST_DUMMY,LEP1ETA_HIST_DUMMY,"hZmumuLep1EtaDiffDummy");
  
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");
  LEP1ETA_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  LEP1ETA_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP1ETA_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP1ETA_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  LEP1ETA_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuLep1Eta("zmmLep1Eta","","",ylabel);
  plotZmumuLep1Eta.AddHist1D(LEP1ETA_HIST_DUMMY);
  plotZmumuLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep1Eta.AddGraph(LEP1ETA_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep1Eta.AddGraph(gTruthLep1EtaMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep1Eta.AddGraph(LEP1ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Eta.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep1Eta.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLep1Eta.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLep1Eta.SetLogx(0);
  plotZmumuLep1Eta.SetLogy(0);
  plotZmumuLep1Eta.SetYRange(70,1.2*(hUnfoldLep1Eta->GetMaximum() + sqrt(hUnfoldLep1Eta->GetMaximum())));
  plotZmumuLep1Eta.TransLegend(0.05,-0.05);
  plotZmumuLep1Eta.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep1EtaDiff("zmmLep1Eta","","|#eta(leading muon)|","Pred/Data");
  plotZmumuLep1EtaDiff.AddHist1D(hZmumuLep1EtaDiffDummy);
  plotZmumuLep1EtaDiff.AddGraph(LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuLep1EtaDiff.AddGraph(LEP1ETA_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuLep1EtaDiff.AddGraph(LEP1ETA_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuLep1EtaDiff.AddGraph(LEP1ETA_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuLep1EtaDiff.AddGraph(LEP1ETA_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuLep1EtaDiff.SetLogx(0);
  plotZmumuLep1EtaDiff.SetYRange(0.8,1.2);

  plotZmumuLep1EtaDiff.AddLine(0, 1,2.4, 1,kBlack,1);
  plotZmumuLep1EtaDiff.AddLine(0, 1.1,2.4, 1.1,kBlack,3);
  plotZmumuLep1EtaDiff.AddLine(0,0.9,2.4,0.9,kBlack,3);

  plotZmumuLep1EtaDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep1Eta2("zmmLep1Etalog","","",ylabel);
  plotZmumuLep1Eta2.AddHist1D(LEP1ETA_HIST_DUMMY);
  plotZmumuLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep1Eta2.AddGraph(LEP1ETA_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep1Eta2.AddGraph(gTruthLep1EtaMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep1Eta2.AddGraph(LEP1ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep1Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep1Eta2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep1Eta2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLep1Eta2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLep1Eta2.SetLogx(0);
  plotZmumuLep1Eta2.SetLogy();
  plotZmumuLep1Eta2.SetYRange(3e-1*(hUnfoldLep1Eta->GetMaximum()),2*(hUnfoldLep1Eta->GetMaximum()));
  plotZmumuLep1Eta2.TransLegend(0.05,-0.05);
  plotZmumuLep1Eta2.Draw(c,kTRUE,format,1);

  //
  // Lep2 Eta
  //   
  TH1D *LEP2ETA_HIST_DUMMY = new TH1D("LEP2ETA_HIST_DUMMY", "LEP2ETA_HIST_DUMMY",24,0,2.4);

  TH1D *hZmumuLep2EtaDiffDummy = makeDiffHist(LEP2ETA_HIST_DUMMY,LEP2ETA_HIST_DUMMY,"hZmumuLep2EtaDiffDummy");
  
  sprintf(ylabel,"d#sigma/d#eta [pb]");
  LEP2ETA_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  LEP2ETA_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP2ETA_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP2ETA_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  LEP2ETA_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuLep2Eta("zmmLep2Eta","","",ylabel);
  plotZmumuLep2Eta.AddHist1D(LEP2ETA_HIST_DUMMY);
  plotZmumuLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep2Eta.AddGraph(LEP2ETA_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep2Eta.AddGraph(gTruthLep2EtaMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep2Eta.AddGraph(LEP2ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Eta.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep2Eta.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLep2Eta.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLep2Eta.SetLogx(0);
  plotZmumuLep2Eta.SetLogy(0);
  plotZmumuLep2Eta.SetYRange(150,1.2*(hUnfoldLep2Eta->GetMaximum() + sqrt(hUnfoldLep2Eta->GetMaximum())));
  plotZmumuLep2Eta.TransLegend(0.05,-0.05);
  plotZmumuLep2Eta.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLep2EtaDiff("zmmLep2Eta","","|#eta(2nd leading muon)|","Pred/Data");
  plotZmumuLep2EtaDiff.AddHist1D(hZmumuLep2EtaDiffDummy);
  plotZmumuLep2EtaDiff.AddGraph(LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuLep2EtaDiff.AddGraph(LEP2ETA_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuLep2EtaDiff.AddGraph(LEP2ETA_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuLep2EtaDiff.AddGraph(LEP2ETA_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuLep2EtaDiff.AddGraph(LEP2ETA_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuLep2EtaDiff.SetLogx(0);
  plotZmumuLep2EtaDiff.SetYRange(0.8,1.2);

  plotZmumuLep2EtaDiff.AddLine(0, 1,2.4, 1,kBlack,1);
  plotZmumuLep2EtaDiff.AddLine(0, 1.1,2.4, 1.1,kBlack,3);
  plotZmumuLep2EtaDiff.AddLine(0,0.9,2.4,0.9,kBlack,3);

  plotZmumuLep2EtaDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLep2Eta2("zmmLep2Etalog","","",ylabel);
  plotZmumuLep2Eta2.AddHist1D(LEP2ETA_HIST_DUMMY);
  plotZmumuLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLep2Eta2.AddGraph(LEP2ETA_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLep2Eta2.AddGraph(gTruthLep2EtaMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLep2Eta2.AddGraph(LEP2ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLep2Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLep2Eta2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLep2Eta2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLep2Eta2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLep2Eta2.SetLogx(0);
  plotZmumuLep2Eta2.SetLogy();
  plotZmumuLep2Eta2.SetYRange(3e-1*(hUnfoldLep2Eta->GetMaximum()),1.8*(hUnfoldLep2Eta->GetMaximum()));
  plotZmumuLep2Eta2.TransLegend(0.05,-0.05);
  plotZmumuLep2Eta2.Draw(c,kTRUE,format,1);

  //
  // LepNeg Pt
  //   
  TH1D *LEPNEGPT_HIST_DUMMY = new TH1D("LEPNEGPT_HIST_DUMMY", "LEPNEGPT_HIST_DUMMY",nBinsLepNegPt,LepNegPtBins);

  TH1D *hZmumuLepNegPtDiffDummy = makeDiffHist(LEPNEGPT_HIST_DUMMY,LEPNEGPT_HIST_DUMMY,"hZmumuLepNegPtDiffDummy");
  
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");
  LEPNEGPT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  LEPNEGPT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEPNEGPT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuLepNegPt("zmmLepNegPt","","",ylabel);
  plotZmumuLepNegPt.AddHist1D(LEPNEGPT_HIST_DUMMY);
  plotZmumuLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLepNegPt.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLepNegPt.AddGraph(gTruthLepNegPtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLepNegPt.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepNegPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLepNegPt.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.65,0.1,0.87,0.15,0);
  plotZmumuLepNegPt.AddTextBox("PDF set: NNPDF3.0",0.65,0.05,0.87,0.1,0);
  plotZmumuLepNegPt.SetLogx();
  plotZmumuLepNegPt.SetLogy(0);
  plotZmumuLepNegPt.SetYRange(0.01,1.2*(hUnfoldLepNegPt->GetMaximum() + sqrt(hUnfoldLepNegPt->GetMaximum())));
  plotZmumuLepNegPt.TransLegend(0.05,-0.05);
  plotZmumuLepNegPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLepNegPtDiff("zmmLepNegPt","","p_{T}^{#mu^{-}} [GeV]","Pred/Data");
  plotZmumuLepNegPtDiff.AddHist1D(hZmumuLepNegPtDiffDummy);
  plotZmumuLepNegPtDiff.AddGraph(LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuLepNegPtDiff.AddGraph(LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuLepNegPtDiff.AddGraph(LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuLepNegPtDiff.AddGraph(LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuLepNegPtDiff.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuLepNegPtDiff.SetLogx();
  plotZmumuLepNegPtDiff.SetYRange(0.8,1.2);

  plotZmumuLepNegPtDiff.AddLine(0, 1,300, 1,kBlack,1);
  plotZmumuLepNegPtDiff.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZmumuLepNegPtDiff.AddLine(0,0.9,300,0.9,kBlack,3);

  plotZmumuLepNegPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLepNegPt2("zmmLepNegPtlog","","",ylabel);
  plotZmumuLepNegPt2.AddHist1D(LEPNEGPT_HIST_DUMMY);
  plotZmumuLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLepNegPt2.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLepNegPt2.AddGraph(gTruthLepNegPtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLepNegPt2.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLepNegPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepNegPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLepNegPt2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLepNegPt2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLepNegPt2.SetLogx();
  plotZmumuLepNegPt2.SetLogy();
  plotZmumuLepNegPt2.SetYRange(5e-5*(hUnfoldLepNegPt->GetMaximum()),10*(hUnfoldLepNegPt->GetMaximum()));
  plotZmumuLepNegPt2.TransLegend(0.05,-0.05);
  plotZmumuLepNegPt2.Draw(c,kTRUE,format,1);

  //
  // LepPos Pt
  //   
  TH1D *LEPPOSPT_HIST_DUMMY = new TH1D("LEPPOSPT_HIST_DUMMY", "LEPPOSPT_HIST_DUMMY",nBinsLepPosPt,LepPosPtBins);

  TH1D *hZmumuLepPosPtDiffDummy = makeDiffHist(LEPPOSPT_HIST_DUMMY,LEPPOSPT_HIST_DUMMY,"hZmumuLepPosPtDiffDummy");
  
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");
  LEPPOSPT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.2);
  LEPPOSPT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEPPOSPT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(kRed-10);
  LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(kRed-10);
  CPlot plotZmumuLepPosPt("zmmLepPosPt","","",ylabel);
  plotZmumuLepPosPt.AddHist1D(LEPPOSPT_HIST_DUMMY);
  plotZmumuLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLepPosPt.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLepPosPt.AddGraph(gTruthLepPosPtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLepPosPt.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepPosPt.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLepPosPt.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.65,0.1,0.87,0.15,0);
  plotZmumuLepPosPt.AddTextBox("PDF set: NNPDF3.0",0.65,0.05,0.87,0.1,0);
  plotZmumuLepPosPt.SetLogx();
  plotZmumuLepPosPt.SetLogy(0);
  plotZmumuLepPosPt.SetYRange(0.01,1.2*(hUnfoldLepPosPt->GetMaximum() + sqrt(hUnfoldLepPosPt->GetMaximum())));
  plotZmumuLepPosPt.TransLegend(0.05,-0.05);
  plotZmumuLepPosPt.Draw(c,kFALSE,format,1);

  CPlot plotZmumuLepPosPtDiff("zmmLepPosPt","","p_{T}^{#mu^{+}} [GeV]","Pred/Data");
  plotZmumuLepPosPtDiff.AddHist1D(hZmumuLepPosPtDiffDummy);
  plotZmumuLepPosPtDiff.AddGraph(LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",kRed+2,21,1);
  plotZmumuLepPosPtDiff.AddGraph(LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",kRed+2,21,1);
  plotZmumuLepPosPtDiff.AddGraph(LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_MADGRAPH_COMP,"0P",kPink,23,1);
  plotZmumuLepPosPtDiff.AddGraph(LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZmumuLepPosPtDiff.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZmumuLepPosPtDiff.SetLogx();
  plotZmumuLepPosPtDiff.SetYRange(0.8,1.2);

  plotZmumuLepPosPtDiff.AddLine(0, 1,300, 1,kBlack,1);
  plotZmumuLepPosPtDiff.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZmumuLepPosPtDiff.AddLine(0,0.9,300,0.9,kBlack,3);

  plotZmumuLepPosPtDiff.Draw(c,kTRUE,format,2);
  
  CPlot plotZmumuLepPosPt2("zmmLepPosPtlog","","",ylabel);
  plotZmumuLepPosPt2.AddHist1D(LEPPOSPT_HIST_DUMMY);
  plotZmumuLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZmumuLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",kRed+2,21,1);
  plotZmumuLepPosPt2.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO,"P",kRed+2,21,1);
  plotZmumuLepPosPt2.AddGraph(gTruthLepPosPtMadgraph,"MADGRAPH","P",kPink,23,1);
  plotZmumuLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZmumuLepPosPt2.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZmumuLepPosPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmumuLepPosPt2.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotZmumuLepPosPt2.AddTextBox("|#eta|<2.4, p_{T}>25 GeV",0.18,0.1,0.4,0.15,0);
  plotZmumuLepPosPt2.AddTextBox("PDF set: NNPDF3.0",0.18,0.05,0.4,0.1,0);
  plotZmumuLepPosPt2.SetLogx();
  plotZmumuLepPosPt2.SetLogy();
  plotZmumuLepPosPt2.SetYRange(5e-5*(hUnfoldLepPosPt->GetMaximum()),10*(hUnfoldLepPosPt->GetMaximum()));
  plotZmumuLepPosPt2.TransLegend(0.05,-0.05);
  plotZmumuLepPosPt2.Draw(c,kTRUE,format,1);
  
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  


  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  cout << endl;
 
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     

  gBenchmark->Show("plotZmm");
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
    //cout<<i<<" "<<EXhigh1[1]<<endl;
    //cout<<i<<" "<<EXlow1[1]<<endl;
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
