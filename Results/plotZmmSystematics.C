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

void plotZmmSystematics(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZmmSystematics");
  gStyle->SetTitleOffset(0.75,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  //
  // input ntuple file names
  //
 vector<TFile*> file;
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStar.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRap.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPt.root", "OPEN"));

  vector<TFile*> fileLumiUp;
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtLumiUp.root", "OPEN"));

  vector<TFile*> fileLumiDown;
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtLumiDown.root", "OPEN"));

  vector<TFile*> fileEWKBkgUp;
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEWKBkgUp.root", "OPEN"));

  vector<TFile*> fileEWKBkgDown;
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEWKBkgDown.root", "OPEN"));

  vector<TFile*> fileTopBkgUp;
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtTopBkgUp.root", "OPEN"));

  vector<TFile*> fileTopBkgDown;
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtTopBkgDown.root", "OPEN"));

  vector<TFile*> fileEffBinSys;
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEffBin.root", "OPEN"));

  vector<TFile*> fileEffStatUp;
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEffStatUp.root", "OPEN"));

  vector<TFile*> fileEffStatDown;
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEffStatDown.root", "OPEN"));

  vector<TFile*> fileEffSigShapeSys;
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEffSigShape.root", "OPEN"));

  vector<TFile*> fileEffBkgShapeSys;
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtEffBkgShape.root", "OPEN"));

  vector<TFile*> fileResScaleSys;
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtResScale.root", "OPEN"));

  vector<TFile*> fileUnfoldModelSys;
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtUnfoldModel_Smoothed.root", "OPEN"));

  vector<TFile*> fileUnfoldMatrixSys;
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZPtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputPhiStarUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputZRapUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1PtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2PtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep1EtaUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLep2EtaUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepNegPtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zmm/UnfoldingOutputLepPosPtUnfoldMatrix.root", "OPEN"));

  
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

  //--------------------------------------------------------------------------
  //                           Z Pt
  //--------------------------------------------------------------------------



  TH1D * hUnfoldZPt;
  TH1D * hTruthZPt;
  
  hUnfoldZPt=(TH1D*)(file[0]->Get("hUnfold"));
  hTruthZPt=(TH1D*)(file[0]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthZPt=TH1TOTGraphAsymmErrors(hTruthZPt);

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

  TGraphAsymmErrors* ZPT_BKG_UNCERT_BAND_DATA;
  ZPT_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEWKBkgDown));

  myAddtoBand(ZPT_TOPBKG_UNCERT_BAND_DATA,ZPT_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_EFF_UNCERT_BAND_DATA;
  ZPT_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEffStatDown));

  //myAddtoBand(ZPT_EFFBIN_UNCERT_BAND_DATA,ZPT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA,ZPT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA,ZPT_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_UNFOLD_UNCERT_BAND_DATA;
  ZPT_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysDown));

  myAddtoBand(ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA,ZPT_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_EXP_UNCERT_BAND_DATA;
  ZPT_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZPtEWKBkgDown));
 
  myAddtoBand(ZPT_TOPBKG_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSTAT_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_EFFBIN_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_RESSCALE_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMODEL_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA,ZPT_EXP_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPt);

  myAddtoBand(ZPT_LUMI_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EWKBKG_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_TOPBKG_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSTAT_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_EFFBIN_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_RESSCALE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMODEL_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
 

  TGraphAsymmErrors* ZPT_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZPT_BKG_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D *ZPT_REL_BKG_SYS = new TH1D("ZPT_REL_BKG_SYS", "ZPT_REL_BKG_SYS",nBinsZPt,ZPtBins);

  double* BKG_SYSTEMATIC_UNCERT_ZPT=ZPT_BKG_SYS_BAND->GetEYhigh();

  cout<<"ZPt Background Systematics"<<endl;

  for(int i=0;i!=ZPT_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_ZPT[i]*100<<endl;
      ZPT_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_ZPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZPT_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZPT_EFF_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D *ZPT_REL_EFF_SYS = new TH1D("ZPT_REL_EFF_SYS", "ZPT_REL_EFF_SYS",nBinsZPt,ZPtBins);

  double* EFF_SYSTEMATIC_UNCERT_ZPT=ZPT_EFF_SYS_BAND->GetEYhigh();

  cout<<"ZPt Efficiency Systematics"<<endl;

  for(int i=0;i!=ZPT_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_ZPT[i]*100<<endl;
      ZPT_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_ZPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZPT_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZPT_RESSCALE_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D *ZPT_REL_RESSCALE_SYS = new TH1D("ZPT_REL_RESSCALE_SYS", "ZPT_REL_RESSCALE_SYS",nBinsZPt,ZPtBins);

  double* RESSCALE_SYSTEMATIC_UNCERT_ZPT=ZPT_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"ZPt Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=ZPT_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_ZPT[i]*100<<endl;
      ZPT_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_ZPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZPT_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZPT_UNFOLD_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D *ZPT_REL_UNFOLD_SYS = new TH1D("ZPT_REL_UNFOLD_SYS", "ZPT_REL_UNFOLD_SYS",nBinsZPt,ZPtBins);

  double* UNFOLD_SYSTEMATIC_UNCERT_ZPT=ZPT_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"ZPt Unfolding Systematics"<<endl;

  for(int i=0;i!=ZPT_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_ZPT[i]*100<<endl;
      ZPT_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_ZPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZPT_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZPT_LUMI_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D *ZPT_REL_LUMI_SYS = new TH1D("ZPT_REL_LUMI_SYS", "ZPT_REL_LUMI_SYS",nBinsZPt,ZPtBins);

  double* LUMI_SYSTEMATIC_UNCERT_ZPT=ZPT_LUMI_SYS_BAND->GetEYhigh();

  cout<<"ZPt Luminosity Uncertainty"<<endl;

  for(int i=0;i!=ZPT_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_ZPT[i]*100<<endl;
      ZPT_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_ZPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZPT_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZPT_EXP_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D *ZPT_REL_EXP_SYS = new TH1D("ZPT_REL_EXP_SYS", "ZPT_REL_EXP_SYS",nBinsZPt,ZPtBins);

  double* EXP_SYSTEMATIC_UNCERT_ZPT=ZPT_EXP_SYS_BAND->GetEYhigh();

  cout<<"ZPt Experimental Uncertainty"<<endl;

  for(int i=0;i!=ZPT_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_ZPT[i]*100<<endl;
      ZPT_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_ZPT[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           PhiStar
  //--------------------------------------------------------------------------



  TH1D * hUnfoldPhiStar;
  TH1D * hTruthPhiStar;
  
  hUnfoldPhiStar=(TH1D*)(file[1]->Get("hUnfold"));
  hTruthPhiStar=(TH1D*)(file[1]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthPhiStar=TH1TOTGraphAsymmErrors(hTruthPhiStar);

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

  TGraphAsymmErrors* PHISTAR_BKG_UNCERT_BAND_DATA;
  PHISTAR_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEWKBkgDown));

  myAddtoBand(PHISTAR_TOPBKG_UNCERT_BAND_DATA,PHISTAR_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_EFF_UNCERT_BAND_DATA;
  PHISTAR_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEffStatDown));

  //myAddtoBand(PHISTAR_EFFBIN_UNCERT_BAND_DATA,PHISTAR_EFF_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA,PHISTAR_EFF_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA,PHISTAR_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_UNFOLD_UNCERT_BAND_DATA;
  PHISTAR_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysDown));

  myAddtoBand(PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA,PHISTAR_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_EXP_UNCERT_BAND_DATA;
  PHISTAR_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarEWKBkgDown));
 
  myAddtoBand(PHISTAR_TOPBKG_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSTAT_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_EFFBIN_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_RESSCALE_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA,PHISTAR_EXP_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);

  myAddtoBand(PHISTAR_LUMI_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EWKBKG_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_TOPBKG_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSTAT_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_EFFBIN_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_RESSCALE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(PHISTAR_BKG_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D *PHISTAR_REL_BKG_SYS = new TH1D("PHISTAR_REL_BKG_SYS", "PHISTAR_REL_BKG_SYS",nBinsPhiStar,PhiStarBins);

  double* BKG_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_BKG_SYS_BAND->GetEYhigh();

  cout<<"PhiStar Background Systematics"<<endl;

  for(int i=0;i!=PHISTAR_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_PHISTAR[i]*100<<endl;
      PHISTAR_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_PHISTAR[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* PHISTAR_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(PHISTAR_EFF_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D *PHISTAR_REL_EFF_SYS = new TH1D("PHISTAR_REL_EFF_SYS", "PHISTAR_REL_EFF_SYS",nBinsPhiStar,PhiStarBins);

  double* EFF_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EFF_SYS_BAND->GetEYhigh();

  cout<<"PhiStar Efficiency Systematics"<<endl;

  for(int i=0;i!=PHISTAR_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_PHISTAR[i]*100<<endl;
      PHISTAR_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_PHISTAR[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* PHISTAR_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(PHISTAR_RESSCALE_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D *PHISTAR_REL_RESSCALE_SYS = new TH1D("PHISTAR_REL_RESSCALE_SYS", "PHISTAR_REL_RESSCALE_SYS",nBinsPhiStar,PhiStarBins);

  double* RESSCALE_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"PhiStar Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=PHISTAR_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_PHISTAR[i]*100<<endl;
      PHISTAR_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_PHISTAR[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* PHISTAR_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(PHISTAR_UNFOLD_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D *PHISTAR_REL_UNFOLD_SYS = new TH1D("PHISTAR_REL_UNFOLD_SYS", "PHISTAR_REL_UNFOLD_SYS",nBinsPhiStar,PhiStarBins);

  double* UNFOLD_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"PhiStar Unfolding Systematics"<<endl;

  for(int i=0;i!=PHISTAR_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_PHISTAR[i]*100<<endl;
      PHISTAR_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_PHISTAR[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* PHISTAR_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(PHISTAR_LUMI_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D *PHISTAR_REL_LUMI_SYS = new TH1D("PHISTAR_REL_LUMI_SYS", "PHISTAR_REL_LUMI_SYS",nBinsPhiStar,PhiStarBins);

  double* LUMI_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_LUMI_SYS_BAND->GetEYhigh();

  cout<<"PhiStar Luminosity Uncertainty"<<endl;

  for(int i=0;i!=PHISTAR_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_PHISTAR[i]*100<<endl;
      PHISTAR_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_PHISTAR[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* PHISTAR_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(PHISTAR_EXP_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D *PHISTAR_REL_EXP_SYS = new TH1D("PHISTAR_REL_EXP_SYS", "PHISTAR_REL_EXP_SYS",nBinsPhiStar,PhiStarBins);

  double* EXP_SYSTEMATIC_UNCERT_PHISTAR=PHISTAR_EXP_SYS_BAND->GetEYhigh();

  cout<<"PhiStar Experimental Uncertainty"<<endl;

  for(int i=0;i!=PHISTAR_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_PHISTAR[i]*100<<endl;
      PHISTAR_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_PHISTAR[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           Z Rapidity
  //--------------------------------------------------------------------------


  TH1D * hUnfoldZRap;
  TH1D * hTruthZRap;
  
  hUnfoldZRap=(TH1D*)(file[2]->Get("hUnfold"));
  hTruthZRap=(TH1D*)(file[2]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthZRap=TH1TOTGraphAsymmErrors(hTruthZRap);

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

  TGraphAsymmErrors* ZRAP_BKG_UNCERT_BAND_DATA;
  ZRAP_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEWKBkgDown));

  myAddtoBand(ZRAP_TOPBKG_UNCERT_BAND_DATA,ZRAP_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_EFF_UNCERT_BAND_DATA;
  ZRAP_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEffStatDown));

  //myAddtoBand(ZRAP_EFFBIN_UNCERT_BAND_DATA,ZRAP_EFF_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA,ZRAP_EFF_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA,ZRAP_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_UNFOLD_UNCERT_BAND_DATA;
  ZRAP_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysDown));

  myAddtoBand(ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA,ZRAP_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_EXP_UNCERT_BAND_DATA;
  ZRAP_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldZRapEWKBkgDown));
  
  myAddtoBand(ZRAP_TOPBKG_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSTAT_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_EFFBIN_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_RESSCALE_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA,ZRAP_EXP_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRap);

  myAddtoBand(ZRAP_LUMI_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EWKBKG_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_TOPBKG_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSTAT_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_EFFBIN_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_RESSCALE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* ZRAP_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZRAP_BKG_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D *ZRAP_REL_BKG_SYS = new TH1D("ZRAP_REL_BKG_SYS", "ZRAP_REL_BKG_SYS",24,0,2.4);

  double* BKG_SYSTEMATIC_UNCERT_ZRAP=ZRAP_BKG_SYS_BAND->GetEYhigh();

  cout<<"ZRap Background Systematics"<<endl;

  for(int i=0;i!=ZRAP_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_ZRAP[i]*100<<endl;
      ZRAP_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_ZRAP[i]);
    }

  cout<<endl<<endl;


  TGraphAsymmErrors* ZRAP_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZRAP_EFF_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D *ZRAP_REL_EFF_SYS = new TH1D("ZRAP_REL_EFF_SYS", "ZRAP_REL_EFF_SYS",24,0,2.4);

  double* EFF_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EFF_SYS_BAND->GetEYhigh();

  cout<<"ZRap Efficiency Systematics"<<endl;

  for(int i=0;i!=ZRAP_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_ZRAP[i]*100<<endl;
      ZRAP_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_ZRAP[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZRAP_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZRAP_RESSCALE_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D *ZRAP_REL_RESSCALE_SYS = new TH1D("ZRAP_REL_RESSCALE_SYS", "ZRAP_REL_RESSCALE_SYS",24,0,2.4);

  double* RESSCALE_SYSTEMATIC_UNCERT_ZRAP=ZRAP_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"ZRap Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=ZRAP_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_ZRAP[i]*100<<endl;
      ZRAP_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_ZRAP[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZRAP_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZRAP_UNFOLD_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D *ZRAP_REL_UNFOLD_SYS = new TH1D("ZRAP_REL_UNFOLD_SYS", "ZRAP_REL_UNFOLD_SYS",24,0,2.4);

  double* UNFOLD_SYSTEMATIC_UNCERT_ZRAP=ZRAP_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"ZRap Unfolding Systematics"<<endl;

  for(int i=0;i!=ZRAP_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_ZRAP[i]*100<<endl;
      ZRAP_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_ZRAP[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZRAP_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZRAP_LUMI_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D *ZRAP_REL_LUMI_SYS = new TH1D("ZRAP_REL_LUMI_SYS", "ZRAP_REL_LUMI_SYS",24,0,2.4);

  double* LUMI_SYSTEMATIC_UNCERT_ZRAP=ZRAP_LUMI_SYS_BAND->GetEYhigh();

  cout<<"ZRap Luminosity Uncertainty"<<endl;

  for(int i=0;i!=ZRAP_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_ZRAP[i]*100<<endl;
      ZRAP_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_ZRAP[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* ZRAP_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(ZRAP_EXP_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D *ZRAP_REL_EXP_SYS = new TH1D("ZRAP_REL_EXP_SYS", "ZRAP_REL_EXP_SYS",24,0,2.4);

  double* EXP_SYSTEMATIC_UNCERT_ZRAP=ZRAP_EXP_SYS_BAND->GetEYhigh();

  cout<<"ZRap Experimental Uncertainty"<<endl;

  for(int i=0;i!=ZRAP_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_ZRAP[i]*100<<endl;
      ZRAP_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_ZRAP[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           Lep1 Pt
  //--------------------------------------------------------------------------



  TH1D * hUnfoldLep1Pt;
  TH1D * hTruthLep1Pt;
  
  hUnfoldLep1Pt=(TH1D*)(file[3]->Get("hUnfold"));
  hTruthLep1Pt=(TH1D*)(file[3]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthLep1Pt=TH1TOTGraphAsymmErrors(hTruthLep1Pt);

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

  TGraphAsymmErrors* LEP1PT_BKG_UNCERT_BAND_DATA;
  LEP1PT_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEWKBkgDown));

  myAddtoBand(LEP1PT_TOPBKG_UNCERT_BAND_DATA,LEP1PT_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_EFF_UNCERT_BAND_DATA;
  LEP1PT_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEffStatDown));

  //myAddtoBand(LEP1PT_EFFBIN_UNCERT_BAND_DATA,LEP1PT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1PT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1PT_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_UNFOLD_UNCERT_BAND_DATA;
  LEP1PT_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysDown));

  myAddtoBand(LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1PT_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_EXP_UNCERT_BAND_DATA;
  LEP1PT_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtEWKBkgDown));
  
  myAddtoBand(LEP1PT_TOPBKG_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSTAT_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_EFFBIN_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_RESSCALE_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1PT_EXP_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEP1PT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);

  myAddtoBand(LEP1PT_LUMI_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EWKBKG_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_TOPBKG_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSTAT_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_EFFBIN_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_RESSCALE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEP1PT_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1PT_BKG_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);


  TH1D *LEP1PT_REL_BKG_SYS = new TH1D("LEP1PT_REL_BKG_SYS", "LEP1PT_REL_BKG_SYS",nBinsLep1Pt,Lep1PtBins);

  double* BKG_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_BKG_SYS_BAND->GetEYhigh();

  cout<<"Lep1Pt Background Systematics"<<endl;

  for(int i=0;i!=LEP1PT_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_LEP1PT[i]*100<<endl;
      LEP1PT_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_LEP1PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1PT_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1PT_EFF_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1PT_REL_EFF_SYS = new TH1D("LEP1PT_REL_EFF_SYS", "LEP1PT_REL_EFF_SYS",nBinsLep1Pt,Lep1PtBins);

  double* EFF_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EFF_SYS_BAND->GetEYhigh();

  cout<<"Lep1Pt Efficiency Systematics"<<endl;

  for(int i=0;i!=LEP1PT_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_LEP1PT[i]*100<<endl;
      LEP1PT_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_LEP1PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1PT_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1PT_RESSCALE_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1PT_REL_RESSCALE_SYS = new TH1D("LEP1PT_REL_RESSCALE_SYS", "LEP1PT_REL_RESSCALE_SYS",nBinsLep1Pt,Lep1PtBins);

  double* RESSCALE_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"Lep1Pt Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=LEP1PT_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_LEP1PT[i]*100<<endl;
      LEP1PT_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_LEP1PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1PT_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1PT_UNFOLD_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1PT_REL_UNFOLD_SYS = new TH1D("LEP1PT_REL_UNFOLD_SYS", "LEP1PT_REL_UNFOLD_SYS",nBinsLep1Pt,Lep1PtBins);

  double* UNFOLD_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"Lep1Pt Unfolding Systematics"<<endl;

  for(int i=0;i!=LEP1PT_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_LEP1PT[i]*100<<endl;
      LEP1PT_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_LEP1PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1PT_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1PT_LUMI_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1PT_REL_LUMI_SYS = new TH1D("LEP1PT_REL_LUMI_SYS", "LEP1PT_REL_LUMI_SYS",nBinsLep1Pt,Lep1PtBins);

  double* LUMI_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_LUMI_SYS_BAND->GetEYhigh();

  cout<<"Lep1Pt Luminosity Uncertainty"<<endl;

  for(int i=0;i!=LEP1PT_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_LEP1PT[i]*100<<endl;
      LEP1PT_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_LEP1PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1PT_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1PT_EXP_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1PT_REL_EXP_SYS = new TH1D("LEP1PT_REL_EXP_SYS", "LEP1PT_REL_EXP_SYS",nBinsLep1Pt,Lep1PtBins);

  double* EXP_SYSTEMATIC_UNCERT_LEP1PT=LEP1PT_EXP_SYS_BAND->GetEYhigh();

  cout<<"Lep1Pt Experimental Uncertainty"<<endl;

  for(int i=0;i!=LEP1PT_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_LEP1PT[i]*100<<endl;
      LEP1PT_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_LEP1PT[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           Lep2 Pt
  //--------------------------------------------------------------------------



  TH1D * hUnfoldLep2Pt;
  TH1D * hTruthLep2Pt;
  
  hUnfoldLep2Pt=(TH1D*)(file[4]->Get("hUnfold"));
  hTruthLep2Pt=(TH1D*)(file[4]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthLep2Pt=TH1TOTGraphAsymmErrors(hTruthLep2Pt);

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
  
  TGraphAsymmErrors* LEP2PT_BKG_UNCERT_BAND_DATA;
  LEP2PT_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEWKBkgDown));

  myAddtoBand(LEP2PT_TOPBKG_UNCERT_BAND_DATA,LEP2PT_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_EFF_UNCERT_BAND_DATA;
  LEP2PT_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEffStatDown));

  //myAddtoBand(LEP2PT_EFFBIN_UNCERT_BAND_DATA,LEP2PT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2PT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2PT_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_UNFOLD_UNCERT_BAND_DATA;
  LEP2PT_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysDown));

  myAddtoBand(LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2PT_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_EXP_UNCERT_BAND_DATA;
  LEP2PT_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtEWKBkgDown));
  
  myAddtoBand(LEP2PT_TOPBKG_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSTAT_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_EFFBIN_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_RESSCALE_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2PT_EXP_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* LEP2PT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);

  myAddtoBand(LEP2PT_LUMI_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EWKBKG_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_TOPBKG_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSTAT_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_EFFBIN_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_RESSCALE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* LEP2PT_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2PT_BKG_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2PT_REL_BKG_SYS = new TH1D("LEP2PT_REL_BKG_SYS", "LEP2PT_REL_BKG_SYS",nBinsLep2Pt,Lep2PtBins);

  double* BKG_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_BKG_SYS_BAND->GetEYhigh();

  cout<<"Lep2Pt Background Systematics"<<endl;

  for(int i=0;i!=LEP2PT_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_LEP2PT[i]*100<<endl;
      LEP2PT_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_LEP2PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2PT_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2PT_EFF_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);


  
  TH1D *LEP2PT_REL_EFF_SYS = new TH1D("LEP2PT_REL_EFF_SYS", "LEP2PT_REL_EFF_SYS",nBinsLep2Pt,Lep2PtBins);

  double* EFF_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EFF_SYS_BAND->GetEYhigh();

  cout<<"Lep2Pt Efficiency Systematics"<<endl;

  for(int i=0;i!=LEP2PT_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_LEP2PT[i]*100<<endl;
      LEP2PT_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_LEP2PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2PT_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2PT_RESSCALE_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2PT_REL_RESSCALE_SYS = new TH1D("LEP2PT_REL_RESSCALE_SYS", "LEP2PT_REL_RESSCALE_SYS",nBinsLep2Pt,Lep2PtBins);

  double* RESSCALE_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"Lep2Pt Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=LEP2PT_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_LEP2PT[i]*100<<endl;
      LEP2PT_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_LEP2PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2PT_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2PT_UNFOLD_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2PT_REL_UNFOLD_SYS = new TH1D("LEP2PT_REL_UNFOLD_SYS", "LEP2PT_REL_UNFOLD_SYS",nBinsLep2Pt,Lep2PtBins);

  double* UNFOLD_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"Lep2Pt Unfolding Systematics"<<endl;

  for(int i=0;i!=LEP2PT_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_LEP2PT[i]*100<<endl;
      LEP2PT_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_LEP2PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2PT_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2PT_LUMI_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2PT_REL_LUMI_SYS = new TH1D("LEP2PT_REL_LUMI_SYS", "LEP2PT_REL_LUMI_SYS",nBinsLep2Pt,Lep2PtBins);

  double* LUMI_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_LUMI_SYS_BAND->GetEYhigh();

  cout<<"Lep2Pt Luminosity Uncertainty"<<endl;

  for(int i=0;i!=LEP2PT_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_LEP2PT[i]*100<<endl;
      LEP2PT_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_LEP2PT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2PT_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2PT_EXP_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2PT_REL_EXP_SYS = new TH1D("LEP2PT_REL_EXP_SYS", "LEP2PT_REL_EXP_SYS",nBinsLep2Pt,Lep2PtBins);

  double* EXP_SYSTEMATIC_UNCERT_LEP2PT=LEP2PT_EXP_SYS_BAND->GetEYhigh();

  cout<<"Lep2Pt Experimental Uncertainty"<<endl;

  for(int i=0;i!=LEP2PT_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_LEP2PT[i]*100<<endl;
      LEP2PT_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_LEP2PT[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           Lep1 Eta
  //--------------------------------------------------------------------------



  TH1D * hUnfoldLep1Eta;
  TH1D * hTruthLep1Eta;
  
  hUnfoldLep1Eta=(TH1D*)(file[5]->Get("hUnfold"));
  hTruthLep1Eta=(TH1D*)(file[5]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthLep1Eta=TH1TOTGraphAsymmErrors(hTruthLep1Eta);

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

  TGraphAsymmErrors* LEP1ETA_BKG_UNCERT_BAND_DATA;
  LEP1ETA_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEWKBkgDown));

  myAddtoBand(LEP1ETA_TOPBKG_UNCERT_BAND_DATA,LEP1ETA_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_EFF_UNCERT_BAND_DATA;
  LEP1ETA_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEffStatDown));

  //myAddtoBand(LEP1ETA_EFFBIN_UNCERT_BAND_DATA,LEP1ETA_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1ETA_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1ETA_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_UNFOLD_UNCERT_BAND_DATA;
  LEP1ETA_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysDown));

  myAddtoBand(LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1ETA_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_EXP_UNCERT_BAND_DATA;
  LEP1ETA_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaEWKBkgDown));
  
  myAddtoBand(LEP1ETA_TOPBKG_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSTAT_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_EFFBIN_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_RESSCALE_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1ETA_EXP_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* LEP1ETA_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);

  myAddtoBand(LEP1ETA_LUMI_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EWKBKG_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_TOPBKG_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSTAT_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_EFFBIN_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_RESSCALE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* LEP1ETA_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_BKG_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1ETA_REL_BKG_SYS = new TH1D("LEP1ETA_REL_BKG_SYS", "LEP1ETA_REL_BKG_SYS",24,0,2.4);

  double* BKG_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_BKG_SYS_BAND->GetEYhigh();

  cout<<"Lep1Eta Background Systematics"<<endl;

  for(int i=0;i!=LEP1ETA_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_LEP1ETA[i]*100<<endl;
      LEP1ETA_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_LEP1ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1ETA_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_EFF_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1ETA_REL_EFF_SYS = new TH1D("LEP1ETA_REL_EFF_SYS", "LEP1ETA_REL_EFF_SYS",24,0,2.4);

  double* EFF_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EFF_SYS_BAND->GetEYhigh();

  cout<<"Lep1Eta Efficiency Systematics"<<endl;

  for(int i=0;i!=LEP1ETA_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_LEP1ETA[i]*100<<endl;
      LEP1ETA_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_LEP1ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1ETA_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_RESSCALE_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1ETA_REL_RESSCALE_SYS = new TH1D("LEP1ETA_REL_RESSCALE_SYS", "LEP1ETA_REL_RESSCALE_SYS",24,0,2.4);

  double* RESSCALE_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"Lep1Eta Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=LEP1ETA_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_LEP1ETA[i]*100<<endl;
      LEP1ETA_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_LEP1ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1ETA_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_UNFOLD_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1ETA_REL_UNFOLD_SYS = new TH1D("LEP1ETA_REL_UNFOLD_SYS", "LEP1ETA_REL_UNFOLD_SYS",24,0,2.4);

  double* UNFOLD_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"Lep1Eta Unfolding Systematics"<<endl;

  for(int i=0;i!=LEP1ETA_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_LEP1ETA[i]*100<<endl;
      LEP1ETA_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_LEP1ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1ETA_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_LUMI_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1ETA_REL_LUMI_SYS = new TH1D("LEP1ETA_REL_LUMI_SYS", "LEP1ETA_REL_LUMI_SYS",24,0,2.4);

  double* LUMI_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_LUMI_SYS_BAND->GetEYhigh();

  cout<<"Lep1Eta Luminosity Uncertainty"<<endl;

  for(int i=0;i!=LEP1ETA_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_LEP1ETA[i]*100<<endl;
      LEP1ETA_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_LEP1ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP1ETA_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_EXP_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP1ETA_REL_EXP_SYS = new TH1D("LEP1ETA_REL_EXP_SYS", "LEP1ETA_REL_EXP_SYS",24,0,2.4);

  double* EXP_SYSTEMATIC_UNCERT_LEP1ETA=LEP1ETA_EXP_SYS_BAND->GetEYhigh();

  cout<<"Lep1Eta Experimental Uncertainty"<<endl;

  for(int i=0;i!=LEP1ETA_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_LEP1ETA[i]*100<<endl;
      LEP1ETA_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_LEP1ETA[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           Lep2 Eta
  //--------------------------------------------------------------------------



  TH1D * hUnfoldLep2Eta;
  TH1D * hTruthLep2Eta;
  
  hUnfoldLep2Eta=(TH1D*)(file[6]->Get("hUnfold"));
  hTruthLep2Eta=(TH1D*)(file[6]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthLep2Eta=TH1TOTGraphAsymmErrors(hTruthLep2Eta);

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

  TGraphAsymmErrors* LEP2ETA_BKG_UNCERT_BAND_DATA;
  LEP2ETA_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEWKBkgDown));

  myAddtoBand(LEP2ETA_TOPBKG_UNCERT_BAND_DATA,LEP2ETA_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_EFF_UNCERT_BAND_DATA;
  LEP2ETA_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEffStatDown));

  //myAddtoBand(LEP2ETA_EFFBIN_UNCERT_BAND_DATA,LEP2ETA_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2ETA_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2ETA_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_UNFOLD_UNCERT_BAND_DATA;
  LEP2ETA_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysDown));

  myAddtoBand(LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2ETA_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_EXP_UNCERT_BAND_DATA;
  LEP2ETA_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaEWKBkgDown));
 
  myAddtoBand(LEP2ETA_TOPBKG_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSTAT_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_EFFBIN_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_RESSCALE_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2ETA_EXP_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* LEP2ETA_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);

  myAddtoBand(LEP2ETA_LUMI_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EWKBKG_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_TOPBKG_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSTAT_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_EFFBIN_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_RESSCALE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* LEP2ETA_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_BKG_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2ETA_REL_BKG_SYS = new TH1D("LEP2ETA_REL_BKG_SYS", "LEP2ETA_REL_BKG_SYS",24,0,2.4);

  double* BKG_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_BKG_SYS_BAND->GetEYhigh();

  cout<<"Lep2Eta Background Systematics"<<endl;

  for(int i=0;i!=LEP2ETA_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_LEP2ETA[i]*100<<endl;
      LEP2ETA_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_LEP2ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2ETA_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_EFF_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2ETA_REL_EFF_SYS = new TH1D("LEP2ETA_REL_EFF_SYS", "LEP2ETA_REL_EFF_SYS",24,0,2.4);

  double* EFF_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EFF_SYS_BAND->GetEYhigh();

  cout<<"Lep2Eta Efficiency Systematics"<<endl;

  for(int i=0;i!=LEP2ETA_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_LEP2ETA[i]*100<<endl;
      LEP2ETA_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_LEP2ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2ETA_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_RESSCALE_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2ETA_REL_RESSCALE_SYS = new TH1D("LEP2ETA_REL_RESSCALE_SYS", "LEP2ETA_REL_RESSCALE_SYS",24,0,2.4);

  double* RESSCALE_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"Lep2Eta Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=LEP2ETA_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_LEP2ETA[i]*100<<endl;
      LEP2ETA_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_LEP2ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2ETA_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_UNFOLD_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2ETA_REL_UNFOLD_SYS = new TH1D("LEP2ETA_REL_UNFOLD_SYS", "LEP2ETA_REL_UNFOLD_SYS",24,0,2.4);

  double* UNFOLD_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"Lep2Eta Unfolding Systematics"<<endl;

  for(int i=0;i!=LEP2ETA_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_LEP2ETA[i]*100<<endl;
      LEP2ETA_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_LEP2ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2ETA_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_LUMI_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2ETA_REL_LUMI_SYS = new TH1D("LEP2ETA_REL_LUMI_SYS", "LEP2ETA_REL_LUMI_SYS",24,0,2.4);

  double* LUMI_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_LUMI_SYS_BAND->GetEYhigh();

  cout<<"Lep2Eta Luminosity Uncertainty"<<endl;

  for(int i=0;i!=LEP2ETA_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_LEP2ETA[i]*100<<endl;
      LEP2ETA_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_LEP2ETA[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEP2ETA_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_EXP_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D *LEP2ETA_REL_EXP_SYS = new TH1D("LEP2ETA_REL_EXP_SYS", "LEP2ETA_REL_EXP_SYS",24,0,2.4);

  double* EXP_SYSTEMATIC_UNCERT_LEP2ETA=LEP2ETA_EXP_SYS_BAND->GetEYhigh();

  cout<<"Lep2Eta Experimental Uncertainty"<<endl;

  for(int i=0;i!=LEP2ETA_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_LEP2ETA[i]*100<<endl;
      LEP2ETA_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_LEP2ETA[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           LepNeg Pt
  //--------------------------------------------------------------------------



  TH1D * hUnfoldLepNegPt;
  TH1D * hTruthLepNegPt;
  
  hUnfoldLepNegPt=(TH1D*)(file[7]->Get("hUnfold"));
  hTruthLepNegPt=(TH1D*)(file[7]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthLepNegPt=TH1TOTGraphAsymmErrors(hTruthLepNegPt);

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

  TGraphAsymmErrors* LEPNEGPT_BKG_UNCERT_BAND_DATA;
  LEPNEGPT_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEWKBkgDown));

  myAddtoBand(LEPNEGPT_TOPBKG_UNCERT_BAND_DATA,LEPNEGPT_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_EFF_UNCERT_BAND_DATA;
  LEPNEGPT_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEffStatDown));

  //myAddtoBand(LEPNEGPT_EFFBIN_UNCERT_BAND_DATA,LEPNEGPT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_UNFOLD_UNCERT_BAND_DATA;
  LEPNEGPT_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysDown));

  myAddtoBand(LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPNEGPT_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_EXP_UNCERT_BAND_DATA;
  LEPNEGPT_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtEWKBkgDown));
  
  myAddtoBand(LEPNEGPT_TOPBKG_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_EFFBIN_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_RESSCALE_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPNEGPT_EXP_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);

  myAddtoBand(LEPNEGPT_LUMI_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EWKBKG_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_TOPBKG_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_EFFBIN_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_RESSCALE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEPNEGPT_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_BKG_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPNEGPT_REL_BKG_SYS = new TH1D("LEPNEGPT_REL_BKG_SYS", "LEPNEGPT_REL_BKG_SYS",nBinsLepNegPt,LepNegPtBins);

  double* BKG_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_BKG_SYS_BAND->GetEYhigh();

  cout<<"LepNegPt Background Systematics"<<endl;

  for(int i=0;i!=LEPNEGPT_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_LEPNEGPT[i]*100<<endl;
      LEPNEGPT_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_LEPNEGPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPNEGPT_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_EFF_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPNEGPT_REL_EFF_SYS = new TH1D("LEPNEGPT_REL_EFF_SYS", "LEPNEGPT_REL_EFF_SYS",nBinsLepNegPt,LepNegPtBins);

  double* EFF_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EFF_SYS_BAND->GetEYhigh();

  cout<<"LepNegPt Efficiency Systematics"<<endl;

  for(int i=0;i!=LEPNEGPT_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_LEPNEGPT[i]*100<<endl;
      LEPNEGPT_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_LEPNEGPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPNEGPT_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_RESSCALE_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);


  TH1D *LEPNEGPT_REL_RESSCALE_SYS = new TH1D("LEPNEGPT_REL_RESSCALE_SYS", "LEPNEGPT_REL_RESSCALE_SYS",nBinsLepNegPt,LepNegPtBins);

  double* RESSCALE_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"LepNegPt Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=LEPNEGPT_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_LEPNEGPT[i]*100<<endl;
      LEPNEGPT_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_LEPNEGPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPNEGPT_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_UNFOLD_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPNEGPT_REL_UNFOLD_SYS = new TH1D("LEPNEGPT_REL_UNFOLD_SYS", "LEPNEGPT_REL_UNFOLD_SYS",nBinsLepNegPt,LepNegPtBins);

  double* UNFOLD_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"LepNegPt Unfolding Systematics"<<endl;

  for(int i=0;i!=LEPNEGPT_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_LEPNEGPT[i]*100<<endl;
      LEPNEGPT_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_LEPNEGPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPNEGPT_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_LUMI_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPNEGPT_REL_LUMI_SYS = new TH1D("LEPNEGPT_REL_LUMI_SYS", "LEPNEGPT_REL_LUMI_SYS",nBinsLepNegPt,LepNegPtBins);

  double* LUMI_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_LUMI_SYS_BAND->GetEYhigh();

  cout<<"LepNegPt Luminosity Uncertainty"<<endl;

  for(int i=0;i!=LEPNEGPT_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_LEPNEGPT[i]*100<<endl;
      LEPNEGPT_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_LEPNEGPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPNEGPT_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_EXP_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPNEGPT_REL_EXP_SYS = new TH1D("LEPNEGPT_REL_EXP_SYS", "LEPNEGPT_REL_EXP_SYS",nBinsLepNegPt,LepNegPtBins);

  double* EXP_SYSTEMATIC_UNCERT_LEPNEGPT=LEPNEGPT_EXP_SYS_BAND->GetEYhigh();

  cout<<"LepNegPt Experimental Uncertainty"<<endl;

  for(int i=0;i!=LEPNEGPT_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_LEPNEGPT[i]*100<<endl;
      LEPNEGPT_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_LEPNEGPT[i]);
    }

  cout<<endl<<endl;

  //--------------------------------------------------------------------------
  //                           LepPos Pt
  //--------------------------------------------------------------------------



  TH1D * hUnfoldLepPosPt;
  TH1D * hTruthLepPosPt;
  
  hUnfoldLepPosPt=(TH1D*)(file[8]->Get("hUnfold"));
  hTruthLepPosPt=(TH1D*)(file[8]->Get("hTruth"));

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
  TGraphAsymmErrors* gTruthLepPosPt=TH1TOTGraphAsymmErrors(hTruthLepPosPt);

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

  TGraphAsymmErrors* LEPPOSPT_BKG_UNCERT_BAND_DATA;
  LEPPOSPT_BKG_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEWKBkgDown));

  myAddtoBand(LEPPOSPT_TOPBKG_UNCERT_BAND_DATA,LEPPOSPT_BKG_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_EFF_UNCERT_BAND_DATA;
  LEPPOSPT_EFF_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffStatUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEffStatDown));

  //myAddtoBand(LEPPOSPT_EFFBIN_UNCERT_BAND_DATA,LEPPOSPT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_EFF_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_EFF_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_UNFOLD_UNCERT_BAND_DATA;
  LEPPOSPT_UNFOLD_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysDown));

  myAddtoBand(LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPPOSPT_UNFOLD_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_EXP_UNCERT_BAND_DATA;
  LEPPOSPT_EXP_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEWKBkgUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtEWKBkgDown));
  
  myAddtoBand(LEPPOSPT_TOPBKG_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_EFFBIN_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_RESSCALE_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPPOSPT_EXP_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);

  myAddtoBand(LEPPOSPT_LUMI_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EWKBKG_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_TOPBKG_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_EFFBIN_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_RESSCALE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEPPOSPT_BKG_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_BKG_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPPOSPT_REL_BKG_SYS = new TH1D("LEPPOSPT_REL_BKG_SYS", "LEPPOSPT_REL_BKG_SYS",nBinsLepPosPt,LepPosPtBins);

  double* BKG_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_BKG_SYS_BAND->GetEYhigh();

  cout<<"LepPosPt Background Systematics"<<endl;

  for(int i=0;i!=LEPPOSPT_REL_BKG_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<BKG_SYSTEMATIC_UNCERT_LEPPOSPT[i]*100<<endl;
      LEPPOSPT_REL_BKG_SYS->SetBinContent(i+1,BKG_SYSTEMATIC_UNCERT_LEPPOSPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPPOSPT_EFF_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_EFF_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPPOSPT_REL_EFF_SYS = new TH1D("LEPPOSPT_REL_EFF_SYS", "LEPPOSPT_REL_EFF_SYS",nBinsLepPosPt,LepPosPtBins);

  double* EFF_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EFF_SYS_BAND->GetEYhigh();

  cout<<"LepPosPt Efficiency Systematics"<<endl;

  for(int i=0;i!=LEPPOSPT_REL_EFF_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EFF_SYSTEMATIC_UNCERT_LEPPOSPT[i]*100<<endl;
      LEPPOSPT_REL_EFF_SYS->SetBinContent(i+1,EFF_SYSTEMATIC_UNCERT_LEPPOSPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPPOSPT_RESSCALE_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_RESSCALE_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);


  TH1D *LEPPOSPT_REL_RESSCALE_SYS = new TH1D("LEPPOSPT_REL_RESSCALE_SYS", "LEPPOSPT_REL_RESSCALE_SYS",nBinsLepPosPt,LepPosPtBins);

  double* RESSCALE_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_RESSCALE_SYS_BAND->GetEYhigh();

  cout<<"LepPosPt Scale and Resolution Systematics"<<endl;

  for(int i=0;i!=LEPPOSPT_REL_RESSCALE_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<RESSCALE_SYSTEMATIC_UNCERT_LEPPOSPT[i]*100<<endl;
      LEPPOSPT_REL_RESSCALE_SYS->SetBinContent(i+1,RESSCALE_SYSTEMATIC_UNCERT_LEPPOSPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPPOSPT_UNFOLD_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_UNFOLD_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPPOSPT_REL_UNFOLD_SYS = new TH1D("LEPPOSPT_REL_UNFOLD_SYS", "LEPPOSPT_REL_UNFOLD_SYS",nBinsLepPosPt,LepPosPtBins);

  double* UNFOLD_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_UNFOLD_SYS_BAND->GetEYhigh();

  cout<<"LepPosPt Unfolding Systematics"<<endl;

  for(int i=0;i!=LEPPOSPT_REL_UNFOLD_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<UNFOLD_SYSTEMATIC_UNCERT_LEPPOSPT[i]*100<<endl;
      LEPPOSPT_REL_UNFOLD_SYS->SetBinContent(i+1,UNFOLD_SYSTEMATIC_UNCERT_LEPPOSPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPPOSPT_LUMI_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_LUMI_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPPOSPT_REL_LUMI_SYS = new TH1D("LEPPOSPT_REL_LUMI_SYS", "LEPPOSPT_REL_LUMI_SYS",nBinsLepPosPt,LepPosPtBins);

  double* LUMI_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_LUMI_SYS_BAND->GetEYhigh();

  cout<<"LepPosPt Luminosity Uncertainty"<<endl;

  for(int i=0;i!=LEPPOSPT_REL_LUMI_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<LUMI_SYSTEMATIC_UNCERT_LEPPOSPT[i]*100<<endl;
      LEPPOSPT_REL_LUMI_SYS->SetBinContent(i+1,LUMI_SYSTEMATIC_UNCERT_LEPPOSPT[i]);
    }

  cout<<endl<<endl;

  TGraphAsymmErrors* LEPPOSPT_EXP_SYS_BAND=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_EXP_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D *LEPPOSPT_REL_EXP_SYS = new TH1D("LEPPOSPT_REL_EXP_SYS", "LEPPOSPT_REL_EXP_SYS",nBinsLepPosPt,LepPosPtBins);

  double* EXP_SYSTEMATIC_UNCERT_LEPPOSPT=LEPPOSPT_EXP_SYS_BAND->GetEYhigh();

  cout<<"LepPosPt Experimental Uncertainty"<<endl;

  for(int i=0;i!=LEPPOSPT_REL_EXP_SYS->GetNbinsX();++i)
    {
      cout<<"Bin "<<i+1<<": "<<EXP_SYSTEMATIC_UNCERT_LEPPOSPT[i]*100<<endl;
      LEPPOSPT_REL_EXP_SYS->SetBinContent(i+1,EXP_SYSTEMATIC_UNCERT_LEPPOSPT[i]);
    }

  cout<<endl<<endl;



  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char xlabel[100];     // string buffer for x-axis label
  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  // plot colors
  Int_t linestyleTotSys = 1;
  Int_t linecolorTotSys = 1;
  Int_t linestyleIDSys = 2;
  Int_t linecolorIDSys = 601;
  Int_t linestyleScaleSys = 1;
  Int_t linecolorScaleSys = 634;
  Int_t linestyleBkgSys = 3;
  Int_t linecolorBkgSys = 419;
  Int_t linestyleUnfoldSys = 4;
  Int_t linecolorUnfoldSys = kGreen;
  Int_t linestyleLumiSys = 9;
  Int_t linecolorLumiSys = 96;

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->cd();
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.15);
  c->SetLeftMargin(0.15);  
  c->SetRightMargin(0.10);  
  c->SetTickx(1);
  c->SetTicky(1);  
  TGaxis::SetMaxDigits(3);

  
  
  //
  // ZPt
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"p_{T}^{#mu^{+}#mu^{-}} [GeV]");
  ZPT_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmPt("zmmPtSys","",xlabel,ylabel);
  plotZmmPt.AddHist1D(ZPT_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmPt.AddHist1D(ZPT_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmPt.AddHist1D(ZPT_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmPt.AddHist1D(ZPT_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmPt.AddHist1D(ZPT_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmPt.AddHist1D(ZPT_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmPt.SetLogx();
  plotZmmPt.SetLogy(0);
  plotZmmPt.SetYRange(0,1.5*(ZPT_REL_EXP_SYS->GetMaximum()));
  plotZmmPt.TransLegend(-0.05,-0.01);
  plotZmmPt.Draw(c,kTRUE,format);

  //
  // PhiStar
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"#phi_{#eta}*");
  PHISTAR_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmPhiStar("zmmPhiStarSys","",xlabel,ylabel);
  plotZmmPhiStar.AddHist1D(PHISTAR_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmPhiStar.AddHist1D(PHISTAR_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmPhiStar.AddHist1D(PHISTAR_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmPhiStar.AddHist1D(PHISTAR_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmPhiStar.AddHist1D(PHISTAR_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmPhiStar.AddHist1D(PHISTAR_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmPhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmPhiStar.SetLogx();
  plotZmmPhiStar.SetLogy(0);
  plotZmmPhiStar.SetYRange(0,1.5*(PHISTAR_REL_EXP_SYS->GetMaximum()));
  plotZmmPhiStar.TransLegend(-0.05,-0.01);
  plotZmmPhiStar.Draw(c,kTRUE,format);

  //
  // ZRap
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"|y^{#mu^{+}#mu^{-}}|");
  ZRAP_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmRap("zmmRapSys","",xlabel,ylabel);
  plotZmmRap.AddHist1D(ZRAP_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmRap.AddHist1D(ZRAP_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmRap.AddHist1D(ZRAP_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmRap.AddHist1D(ZRAP_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmRap.AddHist1D(ZRAP_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmRap.AddHist1D(ZRAP_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmRap.SetLogy(0);
  //plotZmmRap.SetYRange(0,2.0*(ZRAP_REL_EXP_SYS->GetMaximum()));
  plotZmmRap.SetYRange(0,0.1);
  plotZmmRap.TransLegend(-0.05,-0.01);
  plotZmmRap.Draw(c,kTRUE,format);

  //
  // Lep1 Pt
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"p_{T}(leading muon) [GeV]");
  LEP1PT_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmLep1Pt("zmmLep1PtSys","",xlabel,ylabel);
  plotZmmLep1Pt.AddHist1D(LEP1PT_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmLep1Pt.AddHist1D(LEP1PT_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmLep1Pt.AddHist1D(LEP1PT_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmLep1Pt.AddHist1D(LEP1PT_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmLep1Pt.AddHist1D(LEP1PT_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmLep1Pt.AddHist1D(LEP1PT_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmLep1Pt.SetLogx();
  plotZmmLep1Pt.SetLogy(0);
  //plotZmmLep1Pt.SetYRange(0,2*(LEP1PT_REL_EXP_SYS->GetMaximum()));
  plotZmmLep1Pt.SetYRange(0,0.1);
  plotZmmLep1Pt.TransLegend(-0.05,-0.01);
  plotZmmLep1Pt.Draw(c,kTRUE,format);

  //
  // Lep2 Pt
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"p_{T}(2nd leading muon) [GeV]");
  LEP2PT_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmLep2Pt("zmmLep2PtSys","",xlabel,ylabel);
  plotZmmLep2Pt.AddHist1D(LEP2PT_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmLep2Pt.AddHist1D(LEP2PT_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmLep2Pt.AddHist1D(LEP2PT_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmLep2Pt.AddHist1D(LEP2PT_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmLep2Pt.AddHist1D(LEP2PT_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmLep2Pt.AddHist1D(LEP2PT_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmLep2Pt.SetLogx();
  plotZmmLep2Pt.SetLogy(0);
  plotZmmLep2Pt.SetYRange(0,1.5*(LEP2PT_REL_EXP_SYS->GetMaximum()));
  plotZmmLep2Pt.TransLegend(-0.05,-0.01);
  plotZmmLep2Pt.Draw(c,kTRUE,format);

  //
  // Lep1 Eta
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"|#eta|(leading muon)");
  LEP1ETA_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmLep1Eta("zmmLep1EtaSys","",xlabel,ylabel);
  plotZmmLep1Eta.AddHist1D(LEP1ETA_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmLep1Eta.AddHist1D(LEP1ETA_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmLep1Eta.AddHist1D(LEP1ETA_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmLep1Eta.AddHist1D(LEP1ETA_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmLep1Eta.AddHist1D(LEP1ETA_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmLep1Eta.AddHist1D(LEP1ETA_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmLep1Eta.SetLogy(0);
  //plotZmmLep1Eta.SetYRange(0,1.5*(LEP1ETA_REL_EXP_SYS->GetMaximum()));
  plotZmmLep1Eta.SetYRange(0,0.1);
  plotZmmLep1Eta.TransLegend(-0.05,-0.01);
  plotZmmLep1Eta.Draw(c,kTRUE,format);

  //
  // Lep2 Eta
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"|#eta|(2nd leading muon)");
  LEP2ETA_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmLep2Eta("zmmLep2EtaSys","",xlabel,ylabel);
  plotZmmLep2Eta.AddHist1D(LEP2ETA_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmLep2Eta.AddHist1D(LEP2ETA_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmLep2Eta.AddHist1D(LEP2ETA_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmLep2Eta.AddHist1D(LEP2ETA_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmLep2Eta.AddHist1D(LEP2ETA_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmLep2Eta.AddHist1D(LEP2ETA_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmLep2Eta.SetLogy(0);
  //plotZmmLep2Eta.SetYRange(0,1.5*(LEP2ETA_REL_EXP_SYS->GetMaximum()));
  plotZmmLep2Eta.SetYRange(0,0.1);
  plotZmmLep2Eta.TransLegend(-0.05,-0.01);
  plotZmmLep2Eta.Draw(c,kTRUE,format);

  //
  // LepNeg Pt
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"p_{T}^{#mu^{-}} [GeV]");
  LEPNEGPT_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmLepNegPt("zmmLepNegPtSys","",xlabel,ylabel);
  plotZmmLepNegPt.AddHist1D(LEPNEGPT_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmLepNegPt.AddHist1D(LEPNEGPT_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmLepNegPt.AddHist1D(LEPNEGPT_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmLepNegPt.AddHist1D(LEPNEGPT_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmLepNegPt.AddHist1D(LEPNEGPT_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmLepNegPt.AddHist1D(LEPNEGPT_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmLepNegPt.SetLogx();
  plotZmmLepNegPt.SetLogy(0);
  plotZmmLepNegPt.SetYRange(0,1.5*(LEPNEGPT_REL_EXP_SYS->GetMaximum()));
  plotZmmLepNegPt.TransLegend(-0.05,-0.01);
  plotZmmLepNegPt.Draw(c,kTRUE,format);

  //
  // LepPos Pt
  //   

  sprintf(ylabel,"Relative Systematic Uncertainty");
  sprintf(xlabel,"p_{T}^{#mu^{+}} [GeV]");
  LEPPOSPT_REL_EXP_SYS->GetYaxis()->SetTitleOffset(1.4);
  CPlot plotZmmLepPosPt("zmmLepPosPtSys","",xlabel,ylabel);
  plotZmmLepPosPt.AddHist1D(LEPPOSPT_REL_EXP_SYS,"Total Systematic","",linecolorTotSys,linestyleTotSys);
  plotZmmLepPosPt.AddHist1D(LEPPOSPT_REL_EFF_SYS,"Charge, Reco and ID","",linecolorIDSys,linestyleIDSys);
  plotZmmLepPosPt.AddHist1D(LEPPOSPT_REL_RESSCALE_SYS,"Scale and Resolution","",linecolorScaleSys,linestyleScaleSys);
  plotZmmLepPosPt.AddHist1D(LEPPOSPT_REL_BKG_SYS,"Background","",linecolorBkgSys,linestyleBkgSys);
  plotZmmLepPosPt.AddHist1D(LEPPOSPT_REL_UNFOLD_SYS,"Unfolding","",linecolorUnfoldSys,linestyleUnfoldSys);
  plotZmmLepPosPt.AddHist1D(LEPPOSPT_REL_LUMI_SYS,"Luminosity","",linecolorLumiSys,linestyleLumiSys);
  plotZmmLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.88,0);
  plotZmmLepPosPt.SetLogx();
  plotZmmLepPosPt.SetLogy(0);
  plotZmmLepPosPt.SetYRange(0,1.5*(LEPPOSPT_REL_EXP_SYS->GetMaximum()));
  plotZmmLepPosPt.TransLegend(-0.05,-0.01);
  plotZmmLepPosPt.Draw(c,kTRUE,format);


  
  
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

  gBenchmark->Show("plotZmmSystematics");
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
