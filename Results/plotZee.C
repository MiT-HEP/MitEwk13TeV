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
#include "TLatex.h"
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

#include <cmath>

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1);
TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);
TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

//=== MAIN MACRO ================================================================================================= 

void plotZee(const TString  outputDir,   // output directory
             const Double_t lumi         // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotZee");
  gStyle->SetTitleOffset(0.75,"Y");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  //
  // input ntuple file names
  //
  vector<TFile*> file;
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStar.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRap.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2Pt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2Eta.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPt.root", "OPEN"));
  file.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPt.root", "OPEN"));

  vector<TFile*> fileLumiUp;
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtLumiUp.root", "OPEN"));
  fileLumiUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtLumiUp.root", "OPEN"));

  vector<TFile*> fileLumiDown;
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtLumiDown.root", "OPEN"));
  fileLumiDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtLumiDown.root", "OPEN"));

  vector<TFile*> fileEWKBkgUp;
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEWKBkgUp.root", "OPEN"));
  fileEWKBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEWKBkgUp.root", "OPEN"));

  vector<TFile*> fileEWKBkgDown;
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEWKBkgDown.root", "OPEN"));
  fileEWKBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEWKBkgDown.root", "OPEN"));

  vector<TFile*> fileTopBkgUp;
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtTopBkgUp.root", "OPEN"));
  fileTopBkgUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtTopBkgUp.root", "OPEN"));

  vector<TFile*> fileTopBkgDown;
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtTopBkgDown.root", "OPEN"));
  fileTopBkgDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtTopBkgDown.root", "OPEN"));

  vector<TFile*> fileEffBinSys;
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEffBin.root", "OPEN"));
  fileEffBinSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEffBin.root", "OPEN"));

  vector<TFile*> fileEffStatUp;
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEffStatUp.root", "OPEN"));
  fileEffStatUp.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEffStatUp.root", "OPEN"));

  vector<TFile*> fileEffStatDown;
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEffStatDown.root", "OPEN"));
  fileEffStatDown.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEffStatDown.root", "OPEN"));

  vector<TFile*> fileEffSigShapeSys;
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEffSigShape.root", "OPEN"));
  fileEffSigShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEffSigShape.root", "OPEN"));

  vector<TFile*> fileEffBkgShapeSys;
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtEffBkgShape.root", "OPEN"));
  fileEffBkgShapeSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtEffBkgShape.root", "OPEN"));

  /*vector<TFile*> fileResScaleSys;
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtResScale.root", "OPEN"));
  fileResScaleSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtResScale.root", "OPEN"));*/

  vector<TFile*> fileUnfoldModelSys;
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtUnfoldModel_Smoothed.root", "OPEN"));
  fileUnfoldModelSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtUnfoldModel_Smoothed.root", "OPEN"));

  vector<TFile*> fileUnfoldMatrixSys;
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZPtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputPhiStarUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputZRapUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1PtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2PtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep1EtaUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLep2EtaUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepNegPtUnfoldMatrix.root", "OPEN"));
  fileUnfoldMatrixSys.push_back(new TFile("../Unfolding/Zee/UnfoldingOutputLepPosPtUnfoldMatrix.root", "OPEN"));

  TFile* fileMcAtNlo=new TFile("../TheoryUncertainty/Zee/zee_PDFUnc.root", "OPEN");
  TFile* filePowheg=new TFile("../TheoryUncertainty/Zee/zeeph_PDFUnc.root", "OPEN");
  

  TFile* fileFEWZ=new TFile("FEWZ_Histograms.root", "OPEN");

  // plot output file format, all=root+pdf+png
  //const TString format("png");
  const TString format("all");

   
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

  //TH1D * hUnfoldZPtResScaleSys;
  //hUnfoldZPtResScaleSys=(TH1D*)(fileResScaleSys[0]->Get("hUnfold"));

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

  //TGraphAsymmErrors* ZPT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPtResScaleSys);

  TGraphAsymmErrors* ZPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  ZPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZPt,TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldModelSysDown));

  TGraphAsymmErrors* ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPtUnfoldMatrixSys);

  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZPt);

  myAddtoBand(ZPT_LUMI_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_EWKBKG_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_TOPBKG_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_EFFSTAT_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_EFFBIN_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_EFFSIGSHAPE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_EFFBKGSHAPE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZPT_RESSCALE_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMODEL_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZPT_UNFOLDMATRIX_UNCERT_BAND_DATA,ZPT_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_TOT_UNCERT_BAND_DATA,ZPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthZPtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hZPtTruthNominal"));
  TH1D * hTruthZPtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hZPtTruthPDFUp"));
  TH1D * hTruthZPtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hZPtTruthPDFDown"));
  TH1D * hTruthZPtMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hZPtTruthAlphasUp"));
  TH1D * hTruthZPtMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hZPtTruthAlphasDown"));
  TH1D * hTruthZPtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hZPtTruthScaleUp"));
  TH1D * hTruthZPtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hZPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthZPtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthZPtMcAtNlo);

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZPtMcAtNlo);
  
  TGraphAsymmErrors* ZPT_PDF_UNCERT_BAND_AMCATNLO;
  ZPT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloPDFDown));

  TGraphAsymmErrors* ZPT_ALPHAS_UNCERT_BAND_AMCATNLO;
  ZPT_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloAlphasDown));

  TGraphAsymmErrors* ZPT_SCALE_UNCERT_BAND_AMCATNLO;
  ZPT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthZPtMcAtNloScaleDown));

  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZPtMcAtNlo);

  myAddtoBand(ZPT_PDF_UNCERT_BAND_AMCATNLO,ZPT_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(ZPT_ALPHAS_UNCERT_BAND_AMCATNLO,ZPT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(ZPT_SCALE_UNCERT_BAND_AMCATNLO,ZPT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_AMCATNLO,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_TOT_UNCERT_BAND_AMCATNLO,ZPT_STAT_UNCERT_BAND_DATA);

  
  TH1D * hTruthZPtPowheg=(TH1D*)(filePowheg->Get("hZPtTruthNominal"));
  TH1D * hTruthZPtPowhegPDFUp=(TH1D*)(filePowheg->Get("hZPtTruthPDFUp"));
  TH1D * hTruthZPtPowhegPDFDown=(TH1D*)(filePowheg->Get("hZPtTruthPDFDown"));
  TH1D * hTruthZPtPowhegAlphasUp=(TH1D*)(filePowheg->Get("hZPtTruthAlphasUp"));
  TH1D * hTruthZPtPowhegAlphasDown=(TH1D*)(filePowheg->Get("hZPtTruthAlphasDown"));
  TH1D * hTruthZPtPowhegScaleUp=(TH1D*)(filePowheg->Get("hZPtTruthScaleUp"));
  TH1D * hTruthZPtPowhegScaleDown=(TH1D*)(filePowheg->Get("hZPtTruthScaleDown"));
  

  TGraphAsymmErrors* gTruthZPtPowheg=TH1TOTGraphAsymmErrors(hTruthZPtPowheg);

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthZPtPowheg);

  TGraphAsymmErrors* ZPT_PDF_UNCERT_BAND_POWHEG;
  ZPT_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthZPtPowheg,TH1TOTGraphAsymmErrors(hTruthZPtPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthZPtPowhegPDFDown));

  TGraphAsymmErrors* ZPT_ALPHAS_UNCERT_BAND_POWHEG;
  ZPT_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthZPtPowheg,TH1TOTGraphAsymmErrors(hTruthZPtPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthZPtPowhegAlphasDown));
  
  TGraphAsymmErrors* ZPT_SCALE_UNCERT_BAND_POWHEG;
  ZPT_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthZPtPowheg,TH1TOTGraphAsymmErrors(hTruthZPtPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthZPtPowhegScaleDown));

  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthZPtPowheg);

  myAddtoBand(ZPT_PDF_UNCERT_BAND_POWHEG,ZPT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(ZPT_ALPHAS_UNCERT_BAND_POWHEG,ZPT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(ZPT_SCALE_UNCERT_BAND_POWHEG,ZPT_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_POWHEG,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_TOT_UNCERT_BAND_POWHEG,ZPT_STAT_UNCERT_BAND_DATA);


  TGraphAsymmErrors* gTruthZPtFEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZPt_stat"));

  TGraphAsymmErrors* ZPT_STAT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZPt_stat"));

  TGraphAsymmErrors* ZPT_PDF_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZPt_pdf"));

  TGraphAsymmErrors* ZPT_SCALE_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZPt_scale"));

  TGraphAsymmErrors* ZPT_TOT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZPt_stat"));

  myAddtoBand(ZPT_PDF_UNCERT_BAND_FEWZ,ZPT_TOT_UNCERT_BAND_FEWZ);
  myAddtoBand(ZPT_SCALE_UNCERT_BAND_FEWZ,ZPT_TOT_UNCERT_BAND_FEWZ);

  TGraphAsymmErrors* ZPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_STAT_UNCERT_BAND_FEWZ,ZPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(ZPT_TOT_UNCERT_BAND_FEWZ,ZPT_STAT_UNCERT_BAND_DATA);

 
  TH1D * hUnfoldPhiStar;
  hUnfoldPhiStar=(TH1D*)(file[1]->Get("hUnfold"));

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

  //TH1D * hUnfoldPhiStarResScaleSys;
  //hUnfoldPhiStarResScaleSys=(TH1D*)(fileResScaleSys[1]->Get("hUnfold"));

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

  //TGraphAsymmErrors* PHISTAR_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStarResScaleSys);

  TGraphAsymmErrors* PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA;
  PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldPhiStar,TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldModelSysDown));

  TGraphAsymmErrors* PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStarUnfoldMatrixSys);

  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldPhiStar);

  myAddtoBand(PHISTAR_LUMI_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_EWKBKG_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_TOPBKG_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_EFFSTAT_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_EFFBIN_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_EFFSIGSHAPE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_EFFBKGSHAPE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(PHISTAR_RESSCALE_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMODEL_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);
  myAddtoBand(PHISTAR_UNFOLDMATRIX_UNCERT_BAND_DATA,PHISTAR_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_STAT_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_TOT_UNCERT_BAND_DATA,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthPhiStarMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthNominal"));
  TH1D * hTruthPhiStarMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthPDFUp"));
  TH1D * hTruthPhiStarMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthPDFDown"));
  TH1D * hTruthPhiStarMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthAlphasUp"));
  TH1D * hTruthPhiStarMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthAlphasDown"));
  TH1D * hTruthPhiStarMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthScaleUp"));
  TH1D * hTruthPhiStarMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hPhiStarTruthScaleDown"));
  

  TGraphAsymmErrors* gTruthPhiStarMcAtNlo=TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNlo);

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNlo);
  TGraphAsymmErrors* PHISTAR_PDF_UNCERT_BAND_AMCATNLO;
  PHISTAR_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthPhiStarMcAtNlo,TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloPDFDown));

  TGraphAsymmErrors* PHISTAR_ALPHAS_UNCERT_BAND_AMCATNLO;
  PHISTAR_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthPhiStarMcAtNlo,TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloAlphasDown));

  TGraphAsymmErrors* PHISTAR_SCALE_UNCERT_BAND_AMCATNLO;
  PHISTAR_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthPhiStarMcAtNlo,TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNloScaleDown));

  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthPhiStarMcAtNlo);

  myAddtoBand(PHISTAR_PDF_UNCERT_BAND_AMCATNLO,PHISTAR_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(PHISTAR_ALPHAS_UNCERT_BAND_AMCATNLO,PHISTAR_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(PHISTAR_SCALE_UNCERT_BAND_AMCATNLO,PHISTAR_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_STAT_UNCERT_BAND_AMCATNLO,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_TOT_UNCERT_BAND_AMCATNLO,PHISTAR_STAT_UNCERT_BAND_DATA);


  TH1D * hTruthPhiStarPowheg=(TH1D*)(filePowheg->Get("hPhiStarTruthNominal"));
  TH1D * hTruthPhiStarPowhegPDFUp=(TH1D*)(filePowheg->Get("hPhiStarTruthPDFUp"));
  TH1D * hTruthPhiStarPowhegPDFDown=(TH1D*)(filePowheg->Get("hPhiStarTruthPDFDown"));
  TH1D * hTruthPhiStarPowhegAlphasUp=(TH1D*)(filePowheg->Get("hPhiStarTruthAlphasUp"));
  TH1D * hTruthPhiStarPowhegAlphasDown=(TH1D*)(filePowheg->Get("hPhiStarTruthAlphasDown"));
  TH1D * hTruthPhiStarPowhegScaleUp=(TH1D*)(filePowheg->Get("hPhiStarTruthScaleUp"));
  TH1D * hTruthPhiStarPowhegScaleDown=(TH1D*)(filePowheg->Get("hPhiStarTruthScaleDown"));
  

  TGraphAsymmErrors* gTruthPhiStarPowheg=TH1TOTGraphAsymmErrors(hTruthPhiStarPowheg);

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthPhiStarPowheg);
  TGraphAsymmErrors* PHISTAR_PDF_UNCERT_BAND_POWHEG;
  PHISTAR_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthPhiStarPowheg,TH1TOTGraphAsymmErrors(hTruthPhiStarPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthPhiStarPowhegPDFDown));

  TGraphAsymmErrors* PHISTAR_ALPHAS_UNCERT_BAND_POWHEG;
  PHISTAR_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthPhiStarPowheg,TH1TOTGraphAsymmErrors(hTruthPhiStarPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthPhiStarPowhegAlphasDown));

  TGraphAsymmErrors* PHISTAR_SCALE_UNCERT_BAND_POWHEG;
  PHISTAR_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthPhiStarPowheg,TH1TOTGraphAsymmErrors(hTruthPhiStarPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthPhiStarPowhegScaleDown));

  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthPhiStarPowheg);

  myAddtoBand(PHISTAR_PDF_UNCERT_BAND_POWHEG,PHISTAR_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(PHISTAR_ALPHAS_UNCERT_BAND_POWHEG,PHISTAR_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(PHISTAR_SCALE_UNCERT_BAND_POWHEG,PHISTAR_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_STAT_UNCERT_BAND_POWHEG,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_TOT_UNCERT_BAND_POWHEG,PHISTAR_STAT_UNCERT_BAND_DATA);


  TGraphAsymmErrors* gTruthPhiStarFEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("PhiStar_stat"));

  TGraphAsymmErrors* PHISTAR_STAT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("PhiStar_stat"));

  TGraphAsymmErrors* PHISTAR_PDF_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("PhiStar_pdf"));

  TGraphAsymmErrors* PHISTAR_SCALE_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("PhiStar_scale"));

  TGraphAsymmErrors* PHISTAR_TOT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("PhiStar_stat"));

  myAddtoBand(PHISTAR_PDF_UNCERT_BAND_FEWZ,PHISTAR_TOT_UNCERT_BAND_FEWZ);
  myAddtoBand(PHISTAR_SCALE_UNCERT_BAND_FEWZ,PHISTAR_TOT_UNCERT_BAND_FEWZ);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_STAT_UNCERT_BAND_FEWZ,PHISTAR_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(PHISTAR_TOT_UNCERT_BAND_FEWZ,PHISTAR_STAT_UNCERT_BAND_DATA);

  TH1D * hUnfoldZRap;
  hUnfoldZRap=(TH1D*)(file[2]->Get("hUnfold"));

  cout<<"Data: "<<hUnfoldZRap->Integral("width")<<endl;
  
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

  //TH1D * hUnfoldZRapResScaleSys;
  //hUnfoldZRapResScaleSys=(TH1D*)(fileResScaleSys[2]->Get("hUnfold"));

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

  //TGraphAsymmErrors* ZRAP_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRapResScaleSys);

  TGraphAsymmErrors* ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA;
  ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldZRap,TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldModelSysDown));

  TGraphAsymmErrors* ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRapUnfoldMatrixSys);

  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldZRap);

  
  myAddtoBand(ZRAP_LUMI_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_EWKBKG_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_TOPBKG_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_EFFSTAT_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_EFFBIN_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_EFFSIGSHAPE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_EFFBKGSHAPE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(ZRAP_RESSCALE_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMODEL_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  myAddtoBand(ZRAP_UNFOLDMATRIX_UNCERT_BAND_DATA,ZRAP_TOT_UNCERT_BAND_DATA);
  
  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_STAT_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_TOT_UNCERT_BAND_DATA,ZRAP_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthZRapMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hZRapTruthNominal"));
  TH1D * hTruthZRapMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hZRapTruthPDFUp"));
  TH1D * hTruthZRapMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hZRapTruthPDFDown"));
  TH1D * hTruthZRapMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hZRapTruthAlphasUp"));
  TH1D * hTruthZRapMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hZRapTruthAlphasDown"));
  TH1D * hTruthZRapMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hZRapTruthScaleUp"));
  TH1D * hTruthZRapMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hZRapTruthScaleDown"));

  cout<<"aMC@NLO: "<<hTruthZRapMcAtNlo->Integral("width")<<endl;

  TGraphAsymmErrors* gTruthZRapMcAtNlo=TH1TOTGraphAsymmErrors(hTruthZRapMcAtNlo);

  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZRapMcAtNlo);
  
  TGraphAsymmErrors* ZRAP_PDF_UNCERT_BAND_AMCATNLO;
  ZRAP_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZRapMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloPDFDown));

  TGraphAsymmErrors* ZRAP_ALPHAS_UNCERT_BAND_AMCATNLO;
  ZRAP_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZRapMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloAlphasDown));

  TGraphAsymmErrors* ZRAP_SCALE_UNCERT_BAND_AMCATNLO;
  ZRAP_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthZRapMcAtNlo,TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthZRapMcAtNloScaleDown));

  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthZRapMcAtNlo);

  myAddtoBand(ZRAP_PDF_UNCERT_BAND_AMCATNLO,ZRAP_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(ZRAP_ALPHAS_UNCERT_BAND_AMCATNLO,ZRAP_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(ZRAP_SCALE_UNCERT_BAND_AMCATNLO,ZRAP_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_STAT_UNCERT_BAND_AMCATNLO,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_TOT_UNCERT_BAND_AMCATNLO,ZRAP_STAT_UNCERT_BAND_DATA);


  TH1D * hTruthZRapPowheg=(TH1D*)(filePowheg->Get("hZRapTruthNominal"));
  TH1D * hTruthZRapPowhegPDFUp=(TH1D*)(filePowheg->Get("hZRapTruthPDFUp"));
  TH1D * hTruthZRapPowhegPDFDown=(TH1D*)(filePowheg->Get("hZRapTruthPDFDown"));
  TH1D * hTruthZRapPowhegAlphasUp=(TH1D*)(filePowheg->Get("hZRapTruthAlphasUp"));
  TH1D * hTruthZRapPowhegAlphasDown=(TH1D*)(filePowheg->Get("hZRapTruthAlphasDown"));
  TH1D * hTruthZRapPowhegScaleUp=(TH1D*)(filePowheg->Get("hZRapTruthScaleUp"));
  TH1D * hTruthZRapPowhegScaleDown=(TH1D*)(filePowheg->Get("hZRapTruthScaleDown"));

  cout<<"aMC@NLO: "<<hTruthZRapPowheg->Integral("width")<<endl;

  TGraphAsymmErrors* gTruthZRapPowheg=TH1TOTGraphAsymmErrors(hTruthZRapPowheg);

  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthZRapPowheg);
  
  TGraphAsymmErrors* ZRAP_PDF_UNCERT_BAND_POWHEG;
  ZRAP_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthZRapPowheg,TH1TOTGraphAsymmErrors(hTruthZRapPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthZRapPowhegPDFDown));

  TGraphAsymmErrors* ZRAP_ALPHAS_UNCERT_BAND_POWHEG;
  ZRAP_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthZRapPowheg,TH1TOTGraphAsymmErrors(hTruthZRapPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthZRapPowhegAlphasDown));

  TGraphAsymmErrors* ZRAP_SCALE_UNCERT_BAND_POWHEG;
  ZRAP_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthZRapPowheg,TH1TOTGraphAsymmErrors(hTruthZRapPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthZRapPowhegScaleDown));

  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthZRapPowheg);

  myAddtoBand(ZRAP_PDF_UNCERT_BAND_POWHEG,ZRAP_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(ZRAP_ALPHAS_UNCERT_BAND_POWHEG,ZRAP_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(ZRAP_SCALE_UNCERT_BAND_POWHEG,ZRAP_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_STAT_UNCERT_BAND_POWHEG,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_TOT_UNCERT_BAND_POWHEG,ZRAP_STAT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* gTruthZRapFEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZRap_stat"));

  TGraphAsymmErrors* ZRAP_STAT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZRap_stat"));

  TGraphAsymmErrors* ZRAP_PDF_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZRap_pdf"));

  TGraphAsymmErrors* ZRAP_SCALE_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZRap_scale"));

  TGraphAsymmErrors* ZRAP_TOT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("ZRap_stat"));

  myAddtoBand(ZRAP_PDF_UNCERT_BAND_FEWZ,ZRAP_TOT_UNCERT_BAND_FEWZ);
  myAddtoBand(ZRAP_SCALE_UNCERT_BAND_FEWZ,ZRAP_TOT_UNCERT_BAND_FEWZ);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_STAT_UNCERT_BAND_FEWZ,ZRAP_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(ZRAP_TOT_UNCERT_BAND_FEWZ,ZRAP_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             Lep1 Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep1Pt;
  hUnfoldLep1Pt=(TH1D*)(file[3]->Get("hUnfold"));

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

  //TH1D * hUnfoldLep1PtResScaleSys;
  //hUnfoldLep1PtResScaleSys=(TH1D*)(fileResScaleSys[3]->Get("hUnfold"));

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

  //TGraphAsymmErrors* LEP1PT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1PtResScaleSys);

  TGraphAsymmErrors* LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Pt,TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldModelSysDown));

  TGraphAsymmErrors* LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1PtUnfoldMatrixSys);

  TGraphAsymmErrors* LEP1PT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Pt);

  myAddtoBand(LEP1PT_LUMI_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_EWKBKG_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_TOPBKG_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_EFFSTAT_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_EFFBIN_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1PT_RESSCALE_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1PT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_STAT_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_TOT_UNCERT_BAND_DATA,LEP1PT_STAT_UNCERT_BAND_DATA);
  

  TH1D * hTruthLep1PtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthNominal"));
  TH1D * hTruthLep1PtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthPDFUp"));
  TH1D * hTruthLep1PtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthPDFDown"));
  TH1D * hTruthLep1PtMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthAlphasUp"));
  TH1D * hTruthLep1PtMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthAlphasDown"));
  TH1D * hTruthLep1PtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthScaleUp"));
  TH1D * hTruthLep1PtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep1PtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep1PtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNlo);

  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNlo);
  
  TGraphAsymmErrors* LEP1PT_PDF_UNCERT_BAND_AMCATNLO;
  LEP1PT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloPDFDown));
  
  TGraphAsymmErrors* LEP1PT_ALPHAS_UNCERT_BAND_AMCATNLO;
  LEP1PT_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloAlphasDown));
  
  TGraphAsymmErrors* LEP1PT_SCALE_UNCERT_BAND_AMCATNLO;
  LEP1PT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNloScaleDown));

  TGraphAsymmErrors* LEP1PT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1PtMcAtNlo);

  myAddtoBand(LEP1PT_PDF_UNCERT_BAND_AMCATNLO,LEP1PT_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(LEP1PT_ALPHAS_UNCERT_BAND_AMCATNLO,LEP1PT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP1PT_SCALE_UNCERT_BAND_AMCATNLO,LEP1PT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_STAT_UNCERT_BAND_AMCATNLO,LEP1PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_TOT_UNCERT_BAND_AMCATNLO,LEP1PT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep1PtPowheg=(TH1D*)(filePowheg->Get("hLep1PtTruthNominal"));
  TH1D * hTruthLep1PtPowhegPDFUp=(TH1D*)(filePowheg->Get("hLep1PtTruthPDFUp"));
  TH1D * hTruthLep1PtPowhegPDFDown=(TH1D*)(filePowheg->Get("hLep1PtTruthPDFDown"));
  TH1D * hTruthLep1PtPowhegAlphasUp=(TH1D*)(filePowheg->Get("hLep1PtTruthAlphasUp"));
  TH1D * hTruthLep1PtPowhegAlphasDown=(TH1D*)(filePowheg->Get("hLep1PtTruthAlphasDown"));
  TH1D * hTruthLep1PtPowhegScaleUp=(TH1D*)(filePowheg->Get("hLep1PtTruthScaleUp"));
  TH1D * hTruthLep1PtPowhegScaleDown=(TH1D*)(filePowheg->Get("hLep1PtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep1PtPowheg=TH1TOTGraphAsymmErrors(hTruthLep1PtPowheg);

  TGraphAsymmErrors* LEP1PT_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep1PtPowheg);
  
  TGraphAsymmErrors* LEP1PT_PDF_UNCERT_BAND_POWHEG;
  LEP1PT_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep1PtPowheg,TH1TOTGraphAsymmErrors(hTruthLep1PtPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthLep1PtPowhegPDFDown));
  
  TGraphAsymmErrors* LEP1PT_ALPHAS_UNCERT_BAND_POWHEG;
  LEP1PT_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep1PtPowheg,TH1TOTGraphAsymmErrors(hTruthLep1PtPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep1PtPowhegAlphasDown));
  
  TGraphAsymmErrors* LEP1PT_SCALE_UNCERT_BAND_POWHEG;
  LEP1PT_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep1PtPowheg,TH1TOTGraphAsymmErrors(hTruthLep1PtPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthLep1PtPowhegScaleDown));

  TGraphAsymmErrors* LEP1PT_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep1PtPowheg);

  myAddtoBand(LEP1PT_PDF_UNCERT_BAND_POWHEG,LEP1PT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP1PT_ALPHAS_UNCERT_BAND_POWHEG,LEP1PT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP1PT_SCALE_UNCERT_BAND_POWHEG,LEP1PT_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_STAT_UNCERT_BAND_POWHEG,LEP1PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1PT_TOT_UNCERT_BAND_POWHEG,LEP1PT_STAT_UNCERT_BAND_DATA);

  

  //---------------------------------------------------------------------------
  //                             Lep2 Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep2Pt;
  hUnfoldLep2Pt=(TH1D*)(file[4]->Get("hUnfold"));

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

  //TH1D * hUnfoldLep2PtResScaleSys;
  //hUnfoldLep2PtResScaleSys=(TH1D*)(fileResScaleSys[4]->Get("hUnfold"));

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

  //TGraphAsymmErrors* LEP2PT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2PtResScaleSys);

  TGraphAsymmErrors* LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Pt,TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldModelSysDown));

  TGraphAsymmErrors* LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2PtUnfoldMatrixSys);

  TGraphAsymmErrors* LEP2PT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Pt);

  myAddtoBand(LEP2PT_LUMI_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_EWKBKG_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_TOPBKG_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_EFFSTAT_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_EFFBIN_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2PT_RESSCALE_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2PT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2PT_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_STAT_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_TOT_UNCERT_BAND_DATA,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep2PtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthNominal"));
  TH1D * hTruthLep2PtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthPDFUp"));
  TH1D * hTruthLep2PtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthPDFDown"));
  TH1D * hTruthLep2PtMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthAlphasUp"));
  TH1D * hTruthLep2PtMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthAlphasDown"));
  TH1D * hTruthLep2PtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthScaleUp"));
  TH1D * hTruthLep2PtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep2PtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep2PtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNlo);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNlo);
  
  TGraphAsymmErrors* LEP2PT_PDF_UNCERT_BAND_AMCATNLO;
  LEP2PT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloPDFDown));
  
  TGraphAsymmErrors* LEP2PT_ALPHAS_UNCERT_BAND_AMCATNLO;
  LEP2PT_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloAlphasDown));
  
  TGraphAsymmErrors* LEP2PT_SCALE_UNCERT_BAND_AMCATNLO;
  LEP2PT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2PtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNloScaleDown));

  TGraphAsymmErrors* LEP2PT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2PtMcAtNlo);

  myAddtoBand(LEP2PT_PDF_UNCERT_BAND_AMCATNLO,LEP2PT_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(LEP2PT_ALPHAS_UNCERT_BAND_AMCATNLO,LEP2PT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP2PT_SCALE_UNCERT_BAND_AMCATNLO,LEP2PT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_STAT_UNCERT_BAND_AMCATNLO,LEP2PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_TOT_UNCERT_BAND_AMCATNLO,LEP2PT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep2PtPowheg=(TH1D*)(filePowheg->Get("hLep2PtTruthNominal"));
  TH1D * hTruthLep2PtPowhegPDFUp=(TH1D*)(filePowheg->Get("hLep2PtTruthPDFUp"));
  TH1D * hTruthLep2PtPowhegPDFDown=(TH1D*)(filePowheg->Get("hLep2PtTruthPDFDown"));
  TH1D * hTruthLep2PtPowhegAlphasUp=(TH1D*)(filePowheg->Get("hLep2PtTruthAlphasUp"));
  TH1D * hTruthLep2PtPowhegAlphasDown=(TH1D*)(filePowheg->Get("hLep2PtTruthAlphasDown"));
  TH1D * hTruthLep2PtPowhegScaleUp=(TH1D*)(filePowheg->Get("hLep2PtTruthScaleUp"));
  TH1D * hTruthLep2PtPowhegScaleDown=(TH1D*)(filePowheg->Get("hLep2PtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep2PtPowheg=TH1TOTGraphAsymmErrors(hTruthLep2PtPowheg);

  TGraphAsymmErrors* LEP2PT_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep2PtPowheg);
  
  TGraphAsymmErrors* LEP2PT_PDF_UNCERT_BAND_POWHEG;
  LEP2PT_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep2PtPowheg,TH1TOTGraphAsymmErrors(hTruthLep2PtPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthLep2PtPowhegPDFDown));
  
  TGraphAsymmErrors* LEP2PT_ALPHAS_UNCERT_BAND_POWHEG;
  LEP2PT_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep2PtPowheg,TH1TOTGraphAsymmErrors(hTruthLep2PtPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep2PtPowhegAlphasDown));
  
  TGraphAsymmErrors* LEP2PT_SCALE_UNCERT_BAND_POWHEG;
  LEP2PT_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep2PtPowheg,TH1TOTGraphAsymmErrors(hTruthLep2PtPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthLep2PtPowhegScaleDown));

  TGraphAsymmErrors* LEP2PT_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep2PtPowheg);

  myAddtoBand(LEP2PT_PDF_UNCERT_BAND_POWHEG,LEP2PT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP2PT_ALPHAS_UNCERT_BAND_POWHEG,LEP2PT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP2PT_SCALE_UNCERT_BAND_POWHEG,LEP2PT_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_STAT_UNCERT_BAND_POWHEG,LEP2PT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2PT_TOT_UNCERT_BAND_POWHEG,LEP2PT_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             Lep1 Eta
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep1Eta;
  hUnfoldLep1Eta=(TH1D*)(file[5]->Get("hUnfold"));

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

  //TH1D * hUnfoldLep1EtaResScaleSys;
  //hUnfoldLep1EtaResScaleSys=(TH1D*)(fileResScaleSys[5]->Get("hUnfold"));

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

  //TGraphAsymmErrors* LEP1ETA_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1EtaResScaleSys);

  TGraphAsymmErrors* LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep1Eta,TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldModelSysDown));

  TGraphAsymmErrors* LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1EtaUnfoldMatrixSys);

  TGraphAsymmErrors* LEP1ETA_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep1Eta);

  myAddtoBand(LEP1ETA_LUMI_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_EWKBKG_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_TOPBKG_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_EFFSTAT_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_EFFBIN_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP1ETA_RESSCALE_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP1ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP1ETA_TOT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_STAT_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_TOT_UNCERT_BAND_DATA,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep1EtaMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthNominal"));
  TH1D * hTruthLep1EtaMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthPDFUp"));
  TH1D * hTruthLep1EtaMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthPDFDown"));
  TH1D * hTruthLep1EtaMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthAlphasUp"));
  TH1D * hTruthLep1EtaMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthAlphasDown"));
  TH1D * hTruthLep1EtaMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthScaleUp"));
  TH1D * hTruthLep1EtaMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep1EtaTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep1EtaMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNlo);

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNlo);
  
  TGraphAsymmErrors* LEP1ETA_PDF_UNCERT_BAND_AMCATNLO;
  LEP1ETA_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloPDFDown));
  
  TGraphAsymmErrors* LEP1ETA_ALPHAS_UNCERT_BAND_AMCATNLO;
  LEP1ETA_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloAlphasDown));
  
  TGraphAsymmErrors* LEP1ETA_SCALE_UNCERT_BAND_AMCATNLO;
  LEP1ETA_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep1EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNloScaleDown));

  TGraphAsymmErrors* LEP1ETA_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep1EtaMcAtNlo);

  myAddtoBand(LEP1ETA_PDF_UNCERT_BAND_AMCATNLO,LEP1ETA_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(LEP1ETA_ALPHAS_UNCERT_BAND_AMCATNLO,LEP1ETA_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP1ETA_SCALE_UNCERT_BAND_AMCATNLO,LEP1ETA_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_STAT_UNCERT_BAND_AMCATNLO,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_TOT_UNCERT_BAND_AMCATNLO,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep1EtaPowheg=(TH1D*)(filePowheg->Get("hLep1EtaTruthNominal"));
  TH1D * hTruthLep1EtaPowhegPDFUp=(TH1D*)(filePowheg->Get("hLep1EtaTruthPDFUp"));
  TH1D * hTruthLep1EtaPowhegPDFDown=(TH1D*)(filePowheg->Get("hLep1EtaTruthPDFDown"));
  TH1D * hTruthLep1EtaPowhegAlphasUp=(TH1D*)(filePowheg->Get("hLep1EtaTruthAlphasUp"));
  TH1D * hTruthLep1EtaPowhegAlphasDown=(TH1D*)(filePowheg->Get("hLep1EtaTruthAlphasDown"));
  TH1D * hTruthLep1EtaPowhegScaleUp=(TH1D*)(filePowheg->Get("hLep1EtaTruthScaleUp"));
  TH1D * hTruthLep1EtaPowhegScaleDown=(TH1D*)(filePowheg->Get("hLep1EtaTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep1EtaPowheg=TH1TOTGraphAsymmErrors(hTruthLep1EtaPowheg);

  TGraphAsymmErrors* LEP1ETA_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep1EtaPowheg);
  
  TGraphAsymmErrors* LEP1ETA_PDF_UNCERT_BAND_POWHEG;
  LEP1ETA_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep1EtaPowheg,TH1TOTGraphAsymmErrors(hTruthLep1EtaPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaPowhegPDFDown));
  
  TGraphAsymmErrors* LEP1ETA_ALPHAS_UNCERT_BAND_POWHEG;
  LEP1ETA_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep1EtaPowheg,TH1TOTGraphAsymmErrors(hTruthLep1EtaPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaPowhegAlphasDown));
  
  TGraphAsymmErrors* LEP1ETA_SCALE_UNCERT_BAND_POWHEG;
  LEP1ETA_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep1EtaPowheg,TH1TOTGraphAsymmErrors(hTruthLep1EtaPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthLep1EtaPowhegScaleDown));

  TGraphAsymmErrors* LEP1ETA_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep1EtaPowheg);

  myAddtoBand(LEP1ETA_PDF_UNCERT_BAND_POWHEG,LEP1ETA_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP1ETA_ALPHAS_UNCERT_BAND_POWHEG,LEP1ETA_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP1ETA_SCALE_UNCERT_BAND_POWHEG,LEP1ETA_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_STAT_UNCERT_BAND_POWHEG,LEP1ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP1ETA_TOT_UNCERT_BAND_POWHEG,LEP1ETA_STAT_UNCERT_BAND_DATA);
  
  //---------------------------------------------------------------------------
  //                             Lep2 Eta
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLep2Eta;
  hUnfoldLep2Eta=(TH1D*)(file[6]->Get("hUnfold"));

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

  //TH1D * hUnfoldLep2EtaResScaleSys;
  //hUnfoldLep2EtaResScaleSys=(TH1D*)(fileResScaleSys[6]->Get("hUnfold"));

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

  //TGraphAsymmErrors* LEP2ETA_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2EtaResScaleSys);

  TGraphAsymmErrors* LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLep2Eta,TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldModelSysDown));

  TGraphAsymmErrors* LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2EtaUnfoldMatrixSys);

  TGraphAsymmErrors* LEP2ETA_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLep2Eta);

  myAddtoBand(LEP2ETA_LUMI_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_EWKBKG_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_TOPBKG_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_EFFSTAT_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_EFFBIN_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_EFFSIGSHAPE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_EFFBKGSHAPE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEP2ETA_RESSCALE_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMODEL_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEP2ETA_UNFOLDMATRIX_UNCERT_BAND_DATA,LEP2ETA_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_STAT_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_TOT_UNCERT_BAND_DATA,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep2EtaMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthNominal"));
  TH1D * hTruthLep2EtaMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthPDFUp"));
  TH1D * hTruthLep2EtaMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthPDFDown"));
  TH1D * hTruthLep2EtaMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthAlphasUp"));
  TH1D * hTruthLep2EtaMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthAlphasDown"));
  TH1D * hTruthLep2EtaMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthScaleUp"));
  TH1D * hTruthLep2EtaMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLep2EtaTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep2EtaMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNlo);

  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNlo);
  
  TGraphAsymmErrors* LEP2ETA_PDF_UNCERT_BAND_AMCATNLO;
  LEP2ETA_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloPDFDown));
  
  TGraphAsymmErrors* LEP2ETA_ALPHAS_UNCERT_BAND_AMCATNLO;
  LEP2ETA_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloAlphasDown));
  
  TGraphAsymmErrors* LEP2ETA_SCALE_UNCERT_BAND_AMCATNLO;
  LEP2ETA_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLep2EtaMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNloScaleDown));

  TGraphAsymmErrors* LEP2ETA_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLep2EtaMcAtNlo);

  myAddtoBand(LEP2ETA_PDF_UNCERT_BAND_AMCATNLO,LEP2ETA_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(LEP2ETA_ALPHAS_UNCERT_BAND_AMCATNLO,LEP2ETA_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEP2ETA_SCALE_UNCERT_BAND_AMCATNLO,LEP2ETA_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_STAT_UNCERT_BAND_AMCATNLO,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_TOT_UNCERT_BAND_AMCATNLO,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLep2EtaPowheg=(TH1D*)(filePowheg->Get("hLep2EtaTruthNominal"));
  TH1D * hTruthLep2EtaPowhegPDFUp=(TH1D*)(filePowheg->Get("hLep2EtaTruthPDFUp"));
  TH1D * hTruthLep2EtaPowhegPDFDown=(TH1D*)(filePowheg->Get("hLep2EtaTruthPDFDown"));
  TH1D * hTruthLep2EtaPowhegAlphasUp=(TH1D*)(filePowheg->Get("hLep2EtaTruthAlphasUp"));
  TH1D * hTruthLep2EtaPowhegAlphasDown=(TH1D*)(filePowheg->Get("hLep2EtaTruthAlphasDown"));
  TH1D * hTruthLep2EtaPowhegScaleUp=(TH1D*)(filePowheg->Get("hLep2EtaTruthScaleUp"));
  TH1D * hTruthLep2EtaPowhegScaleDown=(TH1D*)(filePowheg->Get("hLep2EtaTruthScaleDown"));

  TGraphAsymmErrors* gTruthLep2EtaPowheg=TH1TOTGraphAsymmErrors(hTruthLep2EtaPowheg);

  TGraphAsymmErrors* LEP2ETA_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep2EtaPowheg);
  
  TGraphAsymmErrors* LEP2ETA_PDF_UNCERT_BAND_POWHEG;
  LEP2ETA_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep2EtaPowheg,TH1TOTGraphAsymmErrors(hTruthLep2EtaPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaPowhegPDFDown));
  
  TGraphAsymmErrors* LEP2ETA_ALPHAS_UNCERT_BAND_POWHEG;
  LEP2ETA_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep2EtaPowheg,TH1TOTGraphAsymmErrors(hTruthLep2EtaPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaPowhegAlphasDown));
  
  TGraphAsymmErrors* LEP2ETA_SCALE_UNCERT_BAND_POWHEG;
  LEP2ETA_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLep2EtaPowheg,TH1TOTGraphAsymmErrors(hTruthLep2EtaPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthLep2EtaPowhegScaleDown));

  TGraphAsymmErrors* LEP2ETA_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLep2EtaPowheg);

  myAddtoBand(LEP2ETA_PDF_UNCERT_BAND_POWHEG,LEP2ETA_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP2ETA_ALPHAS_UNCERT_BAND_POWHEG,LEP2ETA_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEP2ETA_SCALE_UNCERT_BAND_POWHEG,LEP2ETA_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_STAT_UNCERT_BAND_POWHEG,LEP2ETA_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEP2ETA_TOT_UNCERT_BAND_POWHEG,LEP2ETA_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             LepNeg Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLepNegPt;
  hUnfoldLepNegPt=(TH1D*)(file[7]->Get("hUnfold"));

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

  //TH1D * hUnfoldLepNegPtResScaleSys;
  //hUnfoldLepNegPtResScaleSys=(TH1D*)(fileResScaleSys[7]->Get("hUnfold"));

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

  //TGraphAsymmErrors* LEPNEGPT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPtResScaleSys);

  TGraphAsymmErrors* LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepNegPt,TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldModelSysDown));

  TGraphAsymmErrors* LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPtUnfoldMatrixSys);

  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepNegPt);

  myAddtoBand(LEPNEGPT_LUMI_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_EWKBKG_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_TOPBKG_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_EFFSTAT_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_EFFBIN_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPNEGPT_RESSCALE_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPNEGPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPNEGPT_TOT_UNCERT_BAND_DATA);
 
  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_STAT_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_TOT_UNCERT_BAND_DATA,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLepNegPtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthNominal"));
  TH1D * hTruthLepNegPtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthPDFUp"));
  TH1D * hTruthLepNegPtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthPDFDown"));
  TH1D * hTruthLepNegPtMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthAlphasUp"));
  TH1D * hTruthLepNegPtMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthAlphasDown"));
  TH1D * hTruthLepNegPtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthScaleUp"));
  TH1D * hTruthLepNegPtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLepNegPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLepNegPtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNlo);

  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNlo);
  
  TGraphAsymmErrors* LEPNEGPT_PDF_UNCERT_BAND_AMCATNLO;
  LEPNEGPT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepNegPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloPDFDown));
  
  TGraphAsymmErrors* LEPNEGPT_ALPHAS_UNCERT_BAND_AMCATNLO;
  LEPNEGPT_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepNegPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloAlphasDown));
  
  TGraphAsymmErrors* LEPNEGPT_SCALE_UNCERT_BAND_AMCATNLO;
  LEPNEGPT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepNegPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNloScaleDown));

  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepNegPtMcAtNlo);

  myAddtoBand(LEPNEGPT_PDF_UNCERT_BAND_AMCATNLO,LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(LEPNEGPT_ALPHAS_UNCERT_BAND_AMCATNLO,LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEPNEGPT_SCALE_UNCERT_BAND_AMCATNLO,LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLepNegPtPowheg=(TH1D*)(filePowheg->Get("hLepNegPtTruthNominal"));
  TH1D * hTruthLepNegPtPowhegPDFUp=(TH1D*)(filePowheg->Get("hLepNegPtTruthPDFUp"));
  TH1D * hTruthLepNegPtPowhegPDFDown=(TH1D*)(filePowheg->Get("hLepNegPtTruthPDFDown"));
  TH1D * hTruthLepNegPtPowhegAlphasUp=(TH1D*)(filePowheg->Get("hLepNegPtTruthAlphasUp"));
  TH1D * hTruthLepNegPtPowhegAlphasDown=(TH1D*)(filePowheg->Get("hLepNegPtTruthAlphasDown"));
  TH1D * hTruthLepNegPtPowhegScaleUp=(TH1D*)(filePowheg->Get("hLepNegPtTruthScaleUp"));
  TH1D * hTruthLepNegPtPowhegScaleDown=(TH1D*)(filePowheg->Get("hLepNegPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLepNegPtPowheg=TH1TOTGraphAsymmErrors(hTruthLepNegPtPowheg);

  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLepNegPtPowheg);
  
  TGraphAsymmErrors* LEPNEGPT_PDF_UNCERT_BAND_POWHEG;
  LEPNEGPT_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLepNegPtPowheg,TH1TOTGraphAsymmErrors(hTruthLepNegPtPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtPowhegPDFDown));
  
  TGraphAsymmErrors* LEPNEGPT_ALPHAS_UNCERT_BAND_POWHEG;
  LEPNEGPT_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLepNegPtPowheg,TH1TOTGraphAsymmErrors(hTruthLepNegPtPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtPowhegAlphasDown));
  
  TGraphAsymmErrors* LEPNEGPT_SCALE_UNCERT_BAND_POWHEG;
  LEPNEGPT_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLepNegPtPowheg,TH1TOTGraphAsymmErrors(hTruthLepNegPtPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthLepNegPtPowhegScaleDown));

  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLepNegPtPowheg);

  myAddtoBand(LEPNEGPT_PDF_UNCERT_BAND_POWHEG,LEPNEGPT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEPNEGPT_ALPHAS_UNCERT_BAND_POWHEG,LEPNEGPT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEPNEGPT_SCALE_UNCERT_BAND_POWHEG,LEPNEGPT_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_STAT_UNCERT_BAND_POWHEG,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_TOT_UNCERT_BAND_POWHEG,LEPNEGPT_STAT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* gTruthLepNegPtFEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepNegPt_stat"));

  TGraphAsymmErrors* LEPNEGPT_STAT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepNegPt_stat"));

  TGraphAsymmErrors* LEPNEGPT_PDF_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepNegPt_pdf"));

  TGraphAsymmErrors* LEPNEGPT_SCALE_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepNegPt_scale"));

  TGraphAsymmErrors* LEPNEGPT_TOT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepNegPt_stat"));

  myAddtoBand(LEPNEGPT_PDF_UNCERT_BAND_FEWZ,LEPNEGPT_TOT_UNCERT_BAND_FEWZ);
  myAddtoBand(LEPNEGPT_SCALE_UNCERT_BAND_FEWZ,LEPNEGPT_TOT_UNCERT_BAND_FEWZ);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_STAT_UNCERT_BAND_FEWZ,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(LEPNEGPT_TOT_UNCERT_BAND_FEWZ,LEPNEGPT_STAT_UNCERT_BAND_DATA);

  //---------------------------------------------------------------------------
  //                             LepPos Pt
  //---------------------------------------------------------------------------

  TH1D * hUnfoldLepPosPt;
  hUnfoldLepPosPt=(TH1D*)(file[8]->Get("hUnfold"));

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

  //TH1D * hUnfoldLepPosPtResScaleSys;
  //hUnfoldLepPosPtResScaleSys=(TH1D*)(fileResScaleSys[8]->Get("hUnfold"));

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

  //TGraphAsymmErrors* LEPPOSPT_RESSCALE_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPtResScaleSys);

  TGraphAsymmErrors* LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA;
  LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA=myMakeBandSymmetric(gUnfoldLepPosPt,TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysUp),TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldModelSysDown));

  TGraphAsymmErrors* LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPtUnfoldMatrixSys);

  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_DATA=TH1TOTGraphAsymmErrors(hUnfoldLepPosPt);

  myAddtoBand(LEPPOSPT_LUMI_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_EWKBKG_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_TOPBKG_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_EFFSTAT_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_EFFBIN_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_EFFSIGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_EFFBKGSHAPE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  //myAddtoBand(LEPPOSPT_RESSCALE_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMODEL_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);
  myAddtoBand(LEPPOSPT_UNFOLDMATRIX_UNCERT_BAND_DATA,LEPPOSPT_TOT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_STAT_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_TOT_UNCERT_BAND_DATA,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLepPosPtMcAtNlo=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthNominal"));
  TH1D * hTruthLepPosPtMcAtNloPDFUp=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthPDFUp"));
  TH1D * hTruthLepPosPtMcAtNloPDFDown=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthPDFDown"));
  TH1D * hTruthLepPosPtMcAtNloAlphasUp=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthAlphasUp"));
  TH1D * hTruthLepPosPtMcAtNloAlphasDown=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthAlphasDown"));
  TH1D * hTruthLepPosPtMcAtNloScaleUp=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthScaleUp"));
  TH1D * hTruthLepPosPtMcAtNloScaleDown=(TH1D*)(fileMcAtNlo->Get("hLepPosPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLepPosPtMcAtNlo=TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNlo);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNlo);
  
  TGraphAsymmErrors* LEPPOSPT_PDF_UNCERT_BAND_AMCATNLO;
  LEPPOSPT_PDF_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepPosPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloPDFUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloPDFDown));
  
  TGraphAsymmErrors* LEPPOSPT_ALPHAS_UNCERT_BAND_AMCATNLO;
  LEPPOSPT_ALPHAS_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepPosPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloAlphasUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloAlphasDown));
  
  TGraphAsymmErrors* LEPPOSPT_SCALE_UNCERT_BAND_AMCATNLO;
  LEPPOSPT_SCALE_UNCERT_BAND_AMCATNLO=myMakeBandSymmetric(gTruthLepPosPtMcAtNlo,TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloScaleUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNloScaleDown));

  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO=TH1TOTGraphAsymmErrors(hTruthLepPosPtMcAtNlo);

  myAddtoBand(LEPPOSPT_PDF_UNCERT_BAND_AMCATNLO,LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO);
  //myAddtoBand(LEPPOSPT_ALPHAS_UNCERT_BAND_AMCATNLO,LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO);
  myAddtoBand(LEPPOSPT_SCALE_UNCERT_BAND_AMCATNLO,LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TH1D * hTruthLepPosPtPowheg=(TH1D*)(filePowheg->Get("hLepPosPtTruthNominal"));
  TH1D * hTruthLepPosPtPowhegPDFUp=(TH1D*)(filePowheg->Get("hLepPosPtTruthPDFUp"));
  TH1D * hTruthLepPosPtPowhegPDFDown=(TH1D*)(filePowheg->Get("hLepPosPtTruthPDFDown"));
  TH1D * hTruthLepPosPtPowhegAlphasUp=(TH1D*)(filePowheg->Get("hLepPosPtTruthAlphasUp"));
  TH1D * hTruthLepPosPtPowhegAlphasDown=(TH1D*)(filePowheg->Get("hLepPosPtTruthAlphasDown"));
  TH1D * hTruthLepPosPtPowhegScaleUp=(TH1D*)(filePowheg->Get("hLepPosPtTruthScaleUp"));
  TH1D * hTruthLepPosPtPowhegScaleDown=(TH1D*)(filePowheg->Get("hLepPosPtTruthScaleDown"));

  TGraphAsymmErrors* gTruthLepPosPtPowheg=TH1TOTGraphAsymmErrors(hTruthLepPosPtPowheg);

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLepPosPtPowheg);
  
  TGraphAsymmErrors* LEPPOSPT_PDF_UNCERT_BAND_POWHEG;
  LEPPOSPT_PDF_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLepPosPtPowheg,TH1TOTGraphAsymmErrors(hTruthLepPosPtPowhegPDFUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtPowhegPDFDown));
  
  TGraphAsymmErrors* LEPPOSPT_ALPHAS_UNCERT_BAND_POWHEG;
  LEPPOSPT_ALPHAS_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLepPosPtPowheg,TH1TOTGraphAsymmErrors(hTruthLepPosPtPowhegAlphasUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtPowhegAlphasDown));
  
  TGraphAsymmErrors* LEPPOSPT_SCALE_UNCERT_BAND_POWHEG;
  LEPPOSPT_SCALE_UNCERT_BAND_POWHEG=myMakeBandSymmetric(gTruthLepPosPtPowheg,TH1TOTGraphAsymmErrors(hTruthLepPosPtPowhegScaleUp),TH1TOTGraphAsymmErrors(hTruthLepPosPtPowhegScaleDown));

  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_POWHEG=TH1TOTGraphAsymmErrors(hTruthLepPosPtPowheg);

  myAddtoBand(LEPPOSPT_PDF_UNCERT_BAND_POWHEG,LEPPOSPT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEPPOSPT_ALPHAS_UNCERT_BAND_POWHEG,LEPPOSPT_TOT_UNCERT_BAND_POWHEG);
  myAddtoBand(LEPPOSPT_SCALE_UNCERT_BAND_POWHEG,LEPPOSPT_TOT_UNCERT_BAND_POWHEG);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_STAT_UNCERT_BAND_POWHEG,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_TOT_UNCERT_BAND_POWHEG,LEPPOSPT_STAT_UNCERT_BAND_DATA);
  

  TGraphAsymmErrors* gTruthLepPosPtFEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepPosPt_stat"));

  TGraphAsymmErrors* LEPPOSPT_STAT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepPosPt_stat"));

  TGraphAsymmErrors* LEPPOSPT_PDF_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepPosPt_pdf"));

  TGraphAsymmErrors* LEPPOSPT_SCALE_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepPosPt_scale"));

  TGraphAsymmErrors* LEPPOSPT_TOT_UNCERT_BAND_FEWZ=(TGraphAsymmErrors*)(fileFEWZ->Get("LepPosPt_stat"));

  myAddtoBand(LEPPOSPT_PDF_UNCERT_BAND_FEWZ,LEPPOSPT_TOT_UNCERT_BAND_FEWZ);
  myAddtoBand(LEPPOSPT_SCALE_UNCERT_BAND_FEWZ,LEPPOSPT_TOT_UNCERT_BAND_FEWZ);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_STAT_UNCERT_BAND_FEWZ,LEPPOSPT_STAT_UNCERT_BAND_DATA);

  TGraphAsymmErrors* LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP=myTGraphErrorsDivide_noErrGraph2(LEPPOSPT_TOT_UNCERT_BAND_FEWZ,LEPPOSPT_STAT_UNCERT_BAND_DATA);
 
 
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  char xlabel[100];     // string buffer for x-axis label
  char ylabel[100];     // string buffer for y-axis label
  
  // label for lumi
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000);  
  
  // plot colors
  Int_t linecolorPowheg = kGreen+3;
  Int_t fillcolorPowheg = kGreen-10;
  Int_t markerstylePowheg = 22;
  Int_t linecolorAMCAtNlo = kRed+2;
  Int_t fillcolorAMCAtNlo = kRed-10;
  Int_t markerstyleAMCAtNLO = 21;
  Int_t linecolorFEWZ = kBlue;
  Int_t fillcolorFEWZ = kBlue-10;
  Int_t markerstyleFEWZ = 24;

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

  TCanvas *c1 = MakeCanvas("c1","c1",800,800);
  c1->cd()->SetTopMargin(0.10);
  c1->cd()->SetBottomMargin(0.15);
  c1->cd()->SetLeftMargin(0.15);  
  c1->cd()->SetRightMargin(0.07);  
  c1->cd()->SetTickx(1);
  c1->cd()->SetTicky(1);  
  TGaxis::SetMaxDigits(3);

  // canvas for ratio, will be modified
  vector<TCanvas*> cR; // SPLIT
  cR.push_back(  MakeCanvas("c_0_ratio","c",800,800) );
  cR.push_back(  MakeCanvas("c_1_ratio","c",800,400) );
  cR.push_back(  MakeCanvas("c_2_ratio","c",800,600) );
  cR.push_back(  MakeCanvas("c_3_ratio","c",800,800) );
  cR.push_back(  MakeCanvas("c_4_ratio","c",800,1000) );
  double topmargin=0.10;
  double bottommargin=0.149;
  double firstextra=0.05; // for CMS Preliminary, etc..
  
  for(int n=2;n<=4 ;++n)
    {
      double base = (1.0 - 0.)*(1.-topmargin-bottommargin-firstextra)/n;
      
      cR[n]->Divide(1,n,0,0);
      cR[n]->Draw();
      for( int subpad=1; subpad<=n; ++subpad)
  	{
	  double ymin = bottommargin + (subpad-1)* base;
	  double ymax = bottommargin + (subpad)* base;
	  
	  cR[n]->cd(subpad)->SetTopMargin(0);
	  cR[n]->cd(subpad)->SetBottomMargin(0);
	  cR[n]->cd(subpad)->SetLeftMargin(0.15);  
	  cR[n]->cd(subpad)->SetRightMargin(0.07);  
	  cR[n]->cd(subpad)->SetTickx(1);
	  cR[n]->cd(subpad)->SetTicky(1);  
	  if (subpad==1){
	    ymin=0;
	    cR[n]->cd(subpad)->SetBottomMargin(bottommargin/(bottommargin+base));
	  }
	  if (subpad==n){
	    ymax=1;
	    cR[n]->cd(subpad)->SetTopMargin(topmargin/(topmargin+base+firstextra));
	  }
	  
	  cout <<" subpad "<<subpad <<" : "<<ymin << " -- "<<ymax<<endl;
	  cR[n]->cd(subpad)->SetPad(0,ymin,1,ymax);
	  
  	}
      
    }
  cR[0]->cd();

  int n=4;
  double base = (1.0 - 0.)*(1.-topmargin-bottommargin-firstextra)/n;
  double eps=0.01;



  //
  // ZPt
  //   

  sprintf(xlabel,"p_{T}^{e^{+}e^{-}} [GeV]");
  sprintf(ylabel,"d#sigma/dp_{T}^{e^{+}e^{-}} [pb/GeV]");

  TH1D *ZPT_HIST_DUMMY = new TH1D("ZPT_HIST_DUMMY", "ZPT_HIST_DUMMY",nBinsZPt,ZPtBins);

  ZPT_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  ZPT_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  ZPT_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  ZPT_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  ZPT_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  ZPT_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  ZPT_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  ZPT_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  ZPT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  ZPT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  ZPT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  ZPT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  ZPT_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  ZPT_TOT_UNCERT_BAND_FEWZ->SetFillColor(fillcolorFEWZ);
  ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP->SetFillColor(fillcolorFEWZ);

  //linear scale

  CPlot plotZeePt("zeePt","",xlabel,ylabel);
  plotZeePt.AddHist1D(ZPT_HIST_DUMMY);
  plotZeePt.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeePt.AddGraph(ZPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePt.AddGraph(ZPT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeePt.AddGraph(ZPT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeePt.AddGraph(ZPT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeePt.AddGraph(ZPT_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePt.AddGraph(ZPT_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePt.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeePt.AddGraph(ZPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeePt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeePt.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeePt.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeePt.SetLogx();
  plotZeePt.SetLogy(0);
  plotZeePt.SetYRange(0,1.18*(hUnfoldZPt->GetMaximum() + sqrt(hUnfoldZPt->GetMaximum())));
  plotZeePt.SetLegend(0.6,0.64,0.95,0.87);
  plotZeePt.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeePt2("zeePtlog","",xlabel,ylabel);
  plotZeePt2.AddHist1D(ZPT_HIST_DUMMY);
  plotZeePt2.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeePt2.AddGraph(ZPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePt2.AddGraph(ZPT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeePt2.AddGraph(ZPT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeePt2.AddGraph(ZPT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeePt2.AddGraph(ZPT_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePt2.AddGraph(ZPT_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePt2.AddGraph(ZPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeePt2.AddGraph(ZPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeePt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeePt2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeePt2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeePt2.SetLogx();
  plotZeePt2.SetLogy();
  plotZeePt2.SetYRange(1e-6*(hUnfoldZPt->GetMaximum()),50*(hUnfoldZPt->GetMaximum()));
  plotZeePt2.SetLegend(0.6,0.64,0.95,0.87);
  plotZeePt2.Draw(c1,kTRUE,format);
  
  // AMCATNLO POWHEG MADG FEWZ
  // draw split
  TH1D *hZeePtDiffDummySplit =  makeDiffHist(ZPT_HIST_DUMMY,ZPT_HIST_DUMMY,"hZeePtDiffDummySplit");
  hZeePtDiffDummySplit->GetYaxis()->SetTitleFont(43);
  hZeePtDiffDummySplit->GetYaxis()->SetTitleSize(27);
  hZeePtDiffDummySplit->GetYaxis()->SetTitleOffset(2);
  hZeePtDiffDummySplit->GetYaxis()->SetLabelFont(43);
  hZeePtDiffDummySplit->GetYaxis()->SetLabelSize(25);
  hZeePtDiffDummySplit->GetYaxis()->SetNdivisions(408);
  hZeePtDiffDummySplit->GetXaxis()->SetTitleFont(43);
  hZeePtDiffDummySplit->GetXaxis()->SetTitleSize(35);
  hZeePtDiffDummySplit->GetXaxis()->SetTitleOffset(3.2);
  hZeePtDiffDummySplit->GetXaxis()->SetLabelFont(43);
  hZeePtDiffDummySplit->GetXaxis()->SetLabelSize(30);
  hZeePtDiffDummySplit->GetXaxis()->SetTickSize(.1);
  
  
  CPlot plotZeePtDiffSplit_AMCATNLO("zeePtSplit","","","aMC@NLO/Data");
  plotZeePtDiffSplit_AMCATNLO.AddHist1D( hZeePtDiffDummySplit);
  plotZeePtDiffSplit_AMCATNLO.AddGraph(ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePtDiffSplit_AMCATNLO.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePtDiffSplit_AMCATNLO.AddGraph(ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeePtDiffSplit_AMCATNLO.AddGraph(ZPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeePtDiffSplit_AMCATNLO.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.51,0.465,0.66,0);
  plotZeePtDiffSplit_AMCATNLO.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.56,0.55,0.765,0.7,0);
  plotZeePtDiffSplit_AMCATNLO.AddTextBox(lumitext,0.69,0.721,0.93,0.867,0);
  
  plotZeePtDiffSplit_AMCATNLO.SetLogx();
  plotZeePtDiffSplit_AMCATNLO.SetYRange(0.8+eps,1.2 + (1.2-0.8) *(firstextra)/ base );
  plotZeePtDiffSplit_AMCATNLO.AddLine(0, 1,1000, 1,kBlack,1);
  plotZeePtDiffSplit_AMCATNLO.AddLine(0, 1.1,1000, 1.1,kBlack,3);
  plotZeePtDiffSplit_AMCATNLO.AddLine(0,0.9,1000,0.9,kBlack,3);
  
  hZeePtDiffDummySplit->GetYaxis()->SetNdivisions(405);
    
  CPlot plotZeePtDiffSplit_POWHEG("zeePtSplit","","","POWHEG/Data");
  plotZeePtDiffSplit_POWHEG.AddHist1D( hZeePtDiffDummySplit);
  plotZeePtDiffSplit_POWHEG.AddGraph(ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP,"0E2",linecolorPowheg,markerstylePowheg,1);
  plotZeePtDiffSplit_POWHEG.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",linecolorPowheg,markerstylePowheg,1);
  plotZeePtDiffSplit_POWHEG.AddGraph(ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeePtDiffSplit_POWHEG.AddGraph(ZPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeePtDiffSplit_POWHEG.SetLogx();
  plotZeePtDiffSplit_POWHEG.SetYRange(0.8+eps,1.2-eps);
  plotZeePtDiffSplit_POWHEG.AddLine(0, 1,1000, 1,kBlack,1);
  plotZeePtDiffSplit_POWHEG.AddLine(0, 1.1,1000, 1.1,kBlack,3);
  plotZeePtDiffSplit_POWHEG.AddLine(0,0.9,1000,0.9,kBlack,3);
  
  CPlot plotZeePtDiffSplit_FEWZ("zeePtSplit","",xlabel,"FEWZ/Data");
  plotZeePtDiffSplit_FEWZ.AddHist1D(hZeePtDiffDummySplit);
  plotZeePtDiffSplit_FEWZ.AddGraph(ZPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP,"0E2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePtDiffSplit_FEWZ.AddGraph(ZPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP,"0P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePtDiffSplit_FEWZ.AddGraph(ZPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeePtDiffSplit_FEWZ.AddGraph(ZPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeePtDiffSplit_FEWZ.SetLogx();
  plotZeePtDiffSplit_FEWZ.SetYRange(0.8+eps,1.2-eps);
  plotZeePtDiffSplit_FEWZ.AddLine(0, 1,1000, 1,kBlack,1);
  plotZeePtDiffSplit_FEWZ.AddLine(0, 1.1,1000, 1.1,kBlack,3);
  plotZeePtDiffSplit_FEWZ.AddLine(0,0.9,1000,0.9,kBlack,3);
  
  plotZeePtDiffSplit_FEWZ.Draw(cR[3],kFALSE,format,1);
  plotZeePtDiffSplit_POWHEG.Draw(cR[3],kFALSE,format,2);
  plotZeePtDiffSplit_AMCATNLO.Draw(cR[3],kTRUE,format,3);

  //
  // PhiStar
  //  

  sprintf(xlabel,"#phi_{#eta}*");
  sprintf(ylabel,"d#sigma/d#phi_{#eta}* [pb]");

  TH1D *PHISTAR_HIST_DUMMY = new TH1D("PHISTAR_HIST_DUMMY", "PHISTAR_HIST_DUMMY",nBinsPhiStar,PhiStarBins);

  PHISTAR_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  PHISTAR_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  PHISTAR_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  PHISTAR_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  PHISTAR_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  PHISTAR_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  PHISTAR_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  PHISTAR_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  PHISTAR_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  PHISTAR_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  PHISTAR_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  PHISTAR_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  PHISTAR_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  PHISTAR_TOT_UNCERT_BAND_FEWZ->SetFillColor(fillcolorFEWZ);
  PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP->SetFillColor(fillcolorFEWZ);

  //linear scale

  CPlot plotZeePhiStar("zeePhiStar","",xlabel,ylabel);
  plotZeePhiStar.AddHist1D(PHISTAR_HIST_DUMMY);
  plotZeePhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeePhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePhiStar.AddGraph(PHISTAR_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeePhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeePhiStar.AddGraph(PHISTAR_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeePhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePhiStar.AddGraph(PHISTAR_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePhiStar.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeePhiStar.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeePhiStar.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeePhiStar.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeePhiStar.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeePhiStar.SetLogx();
  plotZeePhiStar.SetLogy(0);
  plotZeePhiStar.SetYRange(0.01,1.15*(hUnfoldPhiStar->GetMaximum() + sqrt(hUnfoldPhiStar->GetMaximum())));
  plotZeePhiStar.SetLegend(0.6,0.64,0.95,0.87);
  plotZeePhiStar.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeePhiStar2("zeePhiStarlog","",xlabel,ylabel);
  plotZeePhiStar2.AddHist1D(PHISTAR_HIST_DUMMY);
  plotZeePhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeePhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePhiStar2.AddGraph(PHISTAR_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeePhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeePhiStar2.AddGraph(PHISTAR_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeePhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePhiStar2.AddGraph(PHISTAR_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePhiStar2.AddGraph(PHISTAR_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeePhiStar2.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeePhiStar2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeePhiStar2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeePhiStar2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeePhiStar2.SetLogx();
  plotZeePhiStar2.SetLogy();
  plotZeePhiStar2.SetYRange(2e-4*(hUnfoldPhiStar->GetMaximum()),100*(hUnfoldPhiStar->GetMaximum()));
  plotZeePhiStar2.SetLegend(0.6,0.64,0.95,0.87);
  plotZeePhiStar2.Draw(c1,kTRUE,format);
  
  // AMCATNLO POWHEG MADG FEWZ
  // draw split
  TH1D *hZeePhiStarDiffDummySplit =  makeDiffHist(PHISTAR_HIST_DUMMY,PHISTAR_HIST_DUMMY,"hZeePhiStarDiffDummySplit");
  hZeePhiStarDiffDummySplit->GetYaxis()->SetTitleFont(43);
  hZeePhiStarDiffDummySplit->GetYaxis()->SetTitleSize(27);
  hZeePhiStarDiffDummySplit->GetYaxis()->SetTitleOffset(2);
  hZeePhiStarDiffDummySplit->GetYaxis()->SetLabelFont(43);
  hZeePhiStarDiffDummySplit->GetYaxis()->SetLabelSize(25);
  hZeePhiStarDiffDummySplit->GetYaxis()->SetNdivisions(408);
  hZeePhiStarDiffDummySplit->GetXaxis()->SetTitleFont(43);
  hZeePhiStarDiffDummySplit->GetXaxis()->SetTitleSize(35);
  hZeePhiStarDiffDummySplit->GetXaxis()->SetTitleOffset(3.2);
  hZeePhiStarDiffDummySplit->GetXaxis()->SetLabelFont(43);
  hZeePhiStarDiffDummySplit->GetXaxis()->SetLabelSize(30);
  hZeePhiStarDiffDummySplit->GetXaxis()->SetTickSize(.1);
    
  CPlot plotZeePhiStarDiffSplit_AMCATNLO("zeePhiStarSplit","","","aMC@NLO/Data");
  plotZeePhiStarDiffSplit_AMCATNLO.AddHist1D( hZeePhiStarDiffDummySplit);
  plotZeePhiStarDiffSplit_AMCATNLO.AddGraph(PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePhiStarDiffSplit_AMCATNLO.AddGraph(PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeePhiStarDiffSplit_AMCATNLO.AddGraph(PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeePhiStarDiffSplit_AMCATNLO.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeePhiStarDiffSplit_AMCATNLO.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.51,0.465,0.66,0);
  plotZeePhiStarDiffSplit_AMCATNLO.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.68,0.55,0.89,0.7,0);
  plotZeePhiStarDiffSplit_AMCATNLO.AddTextBox(lumitext,0.69,0.721,0.93,0.867,0);
  
  plotZeePhiStarDiffSplit_AMCATNLO.SetLogx();
  plotZeePhiStarDiffSplit_AMCATNLO.SetYRange(0.8+eps,1.2 + (1.2-0.8) *(firstextra)/ base );
  plotZeePhiStarDiffSplit_AMCATNLO.AddLine(0, 1,3, 1,kBlack,1);
  plotZeePhiStarDiffSplit_AMCATNLO.AddLine(0, 1.1,3, 1.1,kBlack,3);
  plotZeePhiStarDiffSplit_AMCATNLO.AddLine(0,0.9,3,0.9,kBlack,3);
  
  hZeePhiStarDiffDummySplit->GetYaxis()->SetNdivisions(405);
  
  CPlot plotZeePhiStarDiffSplit_POWHEG("zeePhiStarSplit","","","POWHEG/Data");
  plotZeePhiStarDiffSplit_POWHEG.AddHist1D( hZeePhiStarDiffDummySplit);
  plotZeePhiStarDiffSplit_POWHEG.AddGraph(PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP,"0E2",linecolorPowheg,markerstylePowheg,1);
  plotZeePhiStarDiffSplit_POWHEG.AddGraph(PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",linecolorPowheg,markerstylePowheg,1);
  plotZeePhiStarDiffSplit_POWHEG.AddGraph(PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeePhiStarDiffSplit_POWHEG.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeePhiStarDiffSplit_POWHEG.SetLogx();
  plotZeePhiStarDiffSplit_POWHEG.SetYRange(0.8+eps,1.2-eps);
  plotZeePhiStarDiffSplit_POWHEG.AddLine(0, 1,3, 1,kBlack,1);
  plotZeePhiStarDiffSplit_POWHEG.AddLine(0, 1.1,3, 1.1,kBlack,3);
  plotZeePhiStarDiffSplit_POWHEG.AddLine(0,0.9,3,0.9,kBlack,3);
  
  CPlot plotZeePhiStarDiffSplit_FEWZ("zeePhiStarSplit","",xlabel,"FEWZ/Data");
  plotZeePhiStarDiffSplit_FEWZ.AddHist1D(hZeePhiStarDiffDummySplit);
  plotZeePhiStarDiffSplit_FEWZ.AddGraph(PHISTAR_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP,"0E2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePhiStarDiffSplit_FEWZ.AddGraph(PHISTAR_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP,"0P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeePhiStarDiffSplit_FEWZ.AddGraph(PHISTAR_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeePhiStarDiffSplit_FEWZ.AddGraph(PHISTAR_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeePhiStarDiffSplit_FEWZ.SetLogx();
  plotZeePhiStarDiffSplit_FEWZ.SetYRange(0.8+eps,1.2-eps);
  plotZeePhiStarDiffSplit_FEWZ.AddLine(0, 1,3, 1,kBlack,1);
  plotZeePhiStarDiffSplit_FEWZ.AddLine(0, 1.1,3, 1.1,kBlack,3);
  plotZeePhiStarDiffSplit_FEWZ.AddLine(0,0.9,3,0.9,kBlack,3);
  
  plotZeePhiStarDiffSplit_FEWZ.Draw(cR[3],kFALSE,format,1);
  plotZeePhiStarDiffSplit_POWHEG.Draw(cR[3],kFALSE,format,2);
  plotZeePhiStarDiffSplit_AMCATNLO.Draw(cR[3],kTRUE,format,3);

  //
  // Z Rapidity
  //  

  sprintf(xlabel,"|y^{e^{+}e^{-}}|");
  sprintf(ylabel,"d#sigma/dy^{e^{+}e^{-}} [pb]");

  TH1D *ZRAP_HIST_DUMMY = new TH1D("ZRAP_HIST_DUMMY", "ZRAP_HIST_DUMMY",24,0,2.4);

  ZRAP_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  ZRAP_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  ZRAP_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  ZRAP_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  ZRAP_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  ZRAP_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  ZRAP_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  ZRAP_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  ZRAP_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  ZRAP_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  ZRAP_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  ZRAP_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  ZRAP_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  ZRAP_TOT_UNCERT_BAND_FEWZ->SetFillColor(fillcolorFEWZ);
  ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP->SetFillColor(fillcolorFEWZ);

  //linear scale

  CPlot plotZeeRap("zeeRap","",xlabel,ylabel);
  plotZeeRap.AddHist1D(ZRAP_HIST_DUMMY);
  plotZeeRap.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeRap.AddGraph(ZRAP_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeRap.AddGraph(ZRAP_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeRap.AddGraph(ZRAP_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeRap.AddGraph(ZRAP_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeRap.AddGraph(ZRAP_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeRap.AddGraph(ZRAP_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeRap.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeRap.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeRap.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeRap.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeRap.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeRap.SetLogx(0);
  plotZeeRap.SetLogy(0);
  plotZeeRap.SetYRange(0.01,1.5*(hUnfoldZRap->GetMaximum() + sqrt(hUnfoldZRap->GetMaximum())));
  plotZeeRap.SetLegend(0.6,0.64,0.95,0.87);
  plotZeeRap.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeRap2("zeeRaplog","",xlabel,ylabel);
  plotZeeRap2.AddHist1D(ZRAP_HIST_DUMMY);
  plotZeeRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeRap2.AddGraph(ZRAP_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeRap2.AddGraph(ZRAP_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeRap2.AddGraph(ZRAP_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeRap2.AddGraph(ZRAP_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeRap2.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeRap2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeRap2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeRap2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeRap2.SetLogx(0);
  plotZeeRap2.SetLogy();
  plotZeeRap2.SetYRange(5e-3*(hUnfoldZRap->GetMaximum()),30*(hUnfoldZRap->GetMaximum()));
  plotZeeRap2.SetLegend(0.6,0.64,0.95,0.87);
  plotZeeRap2.Draw(c1,kTRUE,format);
  
  // AMCATNLO POWHEG MADG FEWZ
  // draw split
  TH1D *hZeeRapDiffDummySplit =  makeDiffHist(ZRAP_HIST_DUMMY,ZRAP_HIST_DUMMY,"hZeeRapDiffDummySplit");
  hZeeRapDiffDummySplit->GetYaxis()->SetTitleFont(43);
  hZeeRapDiffDummySplit->GetYaxis()->SetTitleSize(27);
  hZeeRapDiffDummySplit->GetYaxis()->SetTitleOffset(2);
  hZeeRapDiffDummySplit->GetYaxis()->SetLabelFont(43);
  hZeeRapDiffDummySplit->GetYaxis()->SetLabelSize(25);
  hZeeRapDiffDummySplit->GetYaxis()->SetNdivisions(408);
  hZeeRapDiffDummySplit->GetXaxis()->SetTitleFont(43);
  hZeeRapDiffDummySplit->GetXaxis()->SetTitleSize(35);
  hZeeRapDiffDummySplit->GetXaxis()->SetTitleOffset(3.2);
  hZeeRapDiffDummySplit->GetXaxis()->SetLabelFont(43);
  hZeeRapDiffDummySplit->GetXaxis()->SetLabelSize(30);
  hZeeRapDiffDummySplit->GetXaxis()->SetTickSize(.1);
    
  CPlot plotZeeRapDiffSplit_AMCATNLO("zeeRapSplit","","","aMC@NLO/Data");
  plotZeeRapDiffSplit_AMCATNLO.AddHist1D( hZeeRapDiffDummySplit);
  plotZeeRapDiffSplit_AMCATNLO.AddGraph(ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeRapDiffSplit_AMCATNLO.AddGraph(ZRAP_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeRapDiffSplit_AMCATNLO.AddGraph(ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeRapDiffSplit_AMCATNLO.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeRapDiffSplit_AMCATNLO.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.51,0.465,0.66,0);
  plotZeeRapDiffSplit_AMCATNLO.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.685,0.55,0.895,0.7,0);
  plotZeeRapDiffSplit_AMCATNLO.AddTextBox(lumitext,0.69,0.721,0.93,0.867,0);
  
  plotZeeRapDiffSplit_AMCATNLO.SetLogx(0);
  plotZeeRapDiffSplit_AMCATNLO.SetYRange(0.8+eps,1.2 + (1.2-0.8) *(firstextra)/ base );
  plotZeeRapDiffSplit_AMCATNLO.AddLine(0, 1,2.4, 1,kBlack,1);
  plotZeeRapDiffSplit_AMCATNLO.AddLine(0, 1.1,2.4, 1.1,kBlack,3);
  plotZeeRapDiffSplit_AMCATNLO.AddLine(0,0.9,2.4,0.9,kBlack,3);
  
  hZeeRapDiffDummySplit->GetYaxis()->SetNdivisions(405);
  
  CPlot plotZeeRapDiffSplit_POWHEG("zeeRapSplit","","","POWHEG/Data");
  plotZeeRapDiffSplit_POWHEG.AddHist1D( hZeeRapDiffDummySplit);
  plotZeeRapDiffSplit_POWHEG.AddGraph(ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP,"0E2",linecolorPowheg,markerstylePowheg,1);
  plotZeeRapDiffSplit_POWHEG.AddGraph(ZRAP_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",linecolorPowheg,markerstylePowheg,1);
  plotZeeRapDiffSplit_POWHEG.AddGraph(ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeRapDiffSplit_POWHEG.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeRapDiffSplit_POWHEG.SetLogx(0);
  plotZeeRapDiffSplit_POWHEG.SetYRange(0.8+eps,1.2-eps);
  plotZeeRapDiffSplit_POWHEG.AddLine(0, 1,2.4, 1,kBlack,1);
  plotZeeRapDiffSplit_POWHEG.AddLine(0, 1.1,2.4, 1.1,kBlack,3);
  plotZeeRapDiffSplit_POWHEG.AddLine(0,0.9,2.4,0.9,kBlack,3);
  
  CPlot plotZeeRapDiffSplit_FEWZ("zeeRapSplit","",xlabel,"FEWZ/Data");
  plotZeeRapDiffSplit_FEWZ.AddHist1D(hZeeRapDiffDummySplit);
  plotZeeRapDiffSplit_FEWZ.AddGraph(ZRAP_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP,"0E2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeRapDiffSplit_FEWZ.AddGraph(ZRAP_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP,"0P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeRapDiffSplit_FEWZ.AddGraph(ZRAP_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeRapDiffSplit_FEWZ.AddGraph(ZRAP_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeRapDiffSplit_FEWZ.SetLogx(0);
  plotZeeRapDiffSplit_FEWZ.SetYRange(0.8+eps,1.2-eps);
  plotZeeRapDiffSplit_FEWZ.AddLine(0, 1,2.4, 1,kBlack,1);
  plotZeeRapDiffSplit_FEWZ.AddLine(0, 1.1,2.4, 1.1,kBlack,3);
  plotZeeRapDiffSplit_FEWZ.AddLine(0,0.9,2.4,0.9,kBlack,3);
  
  plotZeeRapDiffSplit_FEWZ.Draw(cR[3],kFALSE,format,1);
  plotZeeRapDiffSplit_POWHEG.Draw(cR[3],kFALSE,format,2);
  plotZeeRapDiffSplit_AMCATNLO.Draw(cR[3],kTRUE,format,3);
 
  //
  // Lep1 Pt
  //   
  sprintf(xlabel,"p_{T}(leading muon) [GeV]");
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");

  TH1D *LEP1PT_HIST_DUMMY = new TH1D("LEP1PT_HIST_DUMMY", "LEP1PT_HIST_DUMMY",nBinsLep1Pt,Lep1PtBins);

  LEP1PT_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  LEP1PT_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  LEP1PT_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  LEP1PT_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  LEP1PT_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  LEP1PT_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  LEP1PT_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  LEP1PT_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  LEP1PT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  LEP1PT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP1PT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP1PT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  LEP1PT_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  LEP1PT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  LEP1PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  
  
  //linear scale

  CPlot plotZeeLep1Pt("zeeLep1Pt","",xlabel,ylabel);
  plotZeeLep1Pt.AddHist1D(LEP1PT_HIST_DUMMY);
  plotZeeLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep1Pt.AddGraph(LEP1PT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Pt.AddGraph(LEP1PT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Pt.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep1Pt.AddGraph(LEP1PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep1Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep1Pt.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.68,0.15,0.9,0.25,0);
  plotZeeLep1Pt.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep1Pt.SetLogx();
  plotZeeLep1Pt.SetLogy(0);
  plotZeeLep1Pt.SetYRange(0.01,1.2*(hUnfoldLep1Pt->GetMaximum() + sqrt(hUnfoldLep1Pt->GetMaximum())));
  plotZeeLep1Pt.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep1Pt.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeLep1Pt2("zeeLep1Ptlog","",xlabel,ylabel);
  plotZeeLep1Pt2.AddHist1D(LEP1PT_HIST_DUMMY);
  plotZeeLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep1Pt2.AddGraph(LEP1PT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Pt2.AddGraph(LEP1PT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Pt2.AddGraph(LEP1PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep1Pt2.AddGraph(LEP1PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep1Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep1Pt2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLep1Pt2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep1Pt2.SetLogx();
  plotZeeLep1Pt2.SetLogy();
  plotZeeLep1Pt2.SetYRange(5e-5*(hUnfoldLep1Pt->GetMaximum()),10*(hUnfoldLep1Pt->GetMaximum()));
  plotZeeLep1Pt2.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep1Pt2.Draw(c1,kTRUE,format);


  //
  // Lep2 Pt
  //   

  sprintf(xlabel,"p_{T}(2nd leading muon) [GeV]");
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");

  TH1D *LEP2PT_HIST_DUMMY = new TH1D("LEP2PT_HIST_DUMMY", "LEP2PT_HIST_DUMMY",nBinsLep2Pt,Lep2PtBins);

  LEP2PT_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  LEP2PT_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  LEP2PT_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  LEP2PT_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  LEP2PT_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  LEP2PT_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  LEP2PT_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  LEP2PT_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  LEP2PT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  LEP2PT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP2PT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP2PT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  LEP2PT_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  LEP2PT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  LEP2PT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  
  
  //linear scale

  CPlot plotZeeLep2Pt("zeeLep2Pt","",xlabel,ylabel);
  plotZeeLep2Pt.AddHist1D(LEP2PT_HIST_DUMMY);
  plotZeeLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep2Pt.AddGraph(LEP2PT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Pt.AddGraph(LEP2PT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Pt.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep2Pt.AddGraph(LEP2PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep2Pt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep2Pt.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.68,0.15,0.9,0.25,0);
  plotZeeLep2Pt.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep2Pt.SetLogx();
  plotZeeLep2Pt.SetLogy(0);
  plotZeeLep2Pt.SetYRange(0.01,1.2*(hUnfoldLep2Pt->GetMaximum() + sqrt(hUnfoldLep2Pt->GetMaximum())));
  plotZeeLep2Pt.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep2Pt.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeLep2Pt2("zeeLep2Ptlog","",xlabel,ylabel);
  plotZeeLep2Pt2.AddHist1D(LEP2PT_HIST_DUMMY);
  plotZeeLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep2Pt2.AddGraph(LEP2PT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Pt2.AddGraph(LEP2PT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Pt2.AddGraph(LEP2PT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep2Pt2.AddGraph(LEP2PT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep2Pt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep2Pt2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLep2Pt2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep2Pt2.SetLogx();
  plotZeeLep2Pt2.SetLogy();
  plotZeeLep2Pt2.SetYRange(5e-5*(hUnfoldLep2Pt->GetMaximum()),10*(hUnfoldLep2Pt->GetMaximum()));
  plotZeeLep2Pt2.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep2Pt2.Draw(c1,kTRUE,format);

  
  //
  // Lep1 Eta
  //   

  sprintf(xlabel,"|#eta(leading muon)|");
  sprintf(ylabel,"d#sigma/d#eta [pb");

  TH1D *LEP1ETA_HIST_DUMMY = new TH1D("LEP1ETA_HIST_DUMMY", "LEP1ETA_HIST_DUMMY",24,0,2.4);

  LEP1ETA_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  LEP1ETA_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  LEP1ETA_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  LEP1ETA_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  LEP1ETA_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  LEP1ETA_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  LEP1ETA_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  LEP1ETA_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  LEP1ETA_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  LEP1ETA_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP1ETA_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP1ETA_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  LEP1ETA_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  LEP1ETA_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  LEP1ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  
  
  //linear scale

  CPlot plotZeeLep1Eta("zeeLep1Eta","",xlabel,ylabel);
  plotZeeLep1Eta.AddHist1D(LEP1ETA_HIST_DUMMY);
  plotZeeLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep1Eta.AddGraph(LEP1ETA_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Eta.AddGraph(LEP1ETA_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Eta.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep1Eta.AddGraph(LEP1ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep1Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep1Eta.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLep1Eta.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep1Eta.SetLogx(0);
  plotZeeLep1Eta.SetLogy(0);
  plotZeeLep1Eta.SetYRange(150,1.2*(hUnfoldLep1Eta->GetMaximum() + sqrt(hUnfoldLep1Eta->GetMaximum())));
  plotZeeLep1Eta.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep1Eta.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeLep1Eta2("zeeLep1Etalog","",xlabel,ylabel);
  plotZeeLep1Eta2.AddHist1D(LEP1ETA_HIST_DUMMY);
  plotZeeLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep1Eta2.AddGraph(LEP1ETA_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Eta2.AddGraph(LEP1ETA_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep1Eta2.AddGraph(LEP1ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep1Eta2.AddGraph(LEP1ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep1Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep1Eta2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLep1Eta2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep1Eta2.SetLogx(0);
  plotZeeLep1Eta2.SetLogy();
  plotZeeLep1Eta2.SetYRange(4e-1*(hUnfoldLep1Eta->GetMaximum()),1.5*(hUnfoldLep1Eta->GetMaximum()));
  plotZeeLep1Eta2.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep1Eta2.Draw(c1,kTRUE,format);

  //
  // Lep2 Eta
  //   

  sprintf(xlabel,"|#eta(2nd leading muon)|");
  sprintf(ylabel,"d#sigma/d#eta [pb");

  TH1D *LEP2ETA_HIST_DUMMY = new TH1D("LEP2ETA_HIST_DUMMY", "LEP2ETA_HIST_DUMMY",24,0,2.4);

  LEP2ETA_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  LEP2ETA_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  LEP2ETA_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  LEP2ETA_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  LEP2ETA_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  LEP2ETA_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  LEP2ETA_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  LEP2ETA_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  LEP2ETA_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  LEP2ETA_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEP2ETA_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEP2ETA_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  LEP2ETA_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  LEP2ETA_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  LEP2ETA_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  
  //linear scale

  CPlot plotZeeLep2Eta("zeeLep2Eta","",xlabel,ylabel);
  plotZeeLep2Eta.AddHist1D(LEP2ETA_HIST_DUMMY);
  plotZeeLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep2Eta.AddGraph(LEP2ETA_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Eta.AddGraph(LEP2ETA_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Eta.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep2Eta.AddGraph(LEP2ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep2Eta.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep2Eta.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLep2Eta.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep2Eta.SetLogx(0);
  plotZeeLep2Eta.SetLogy(0);
  plotZeeLep2Eta.SetYRange(150,1.2*(hUnfoldLep2Eta->GetMaximum() + sqrt(hUnfoldLep2Eta->GetMaximum())));
  plotZeeLep2Eta.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep2Eta.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeLep2Eta2("zeeLep2Etalog","",xlabel,ylabel);
  plotZeeLep2Eta2.AddHist1D(LEP2ETA_HIST_DUMMY);
  plotZeeLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLep2Eta2.AddGraph(LEP2ETA_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Eta2.AddGraph(LEP2ETA_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLep2Eta2.AddGraph(LEP2ETA_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLep2Eta2.AddGraph(LEP2ETA_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLep2Eta2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLep2Eta2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLep2Eta2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLep2Eta2.SetLogx(0);
  plotZeeLep2Eta2.SetLogy();
  plotZeeLep2Eta2.SetYRange(4e-1*(hUnfoldLep2Eta->GetMaximum()),1.5*(hUnfoldLep2Eta->GetMaximum()));
  plotZeeLep2Eta2.SetLegend(0.6,0.685,0.95,0.87);
  plotZeeLep2Eta2.Draw(c1,kTRUE,format);

  //
  // LepNeg Pt
  //  

  sprintf(xlabel,"p_{T}^{e^{-}} [GeV]");
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");

  TH1D *LEPNEGPT_HIST_DUMMY = new TH1D("LEPNEGPT_HIST_DUMMY", "LEPNEGPT_HIST_DUMMY",nBinsLepNegPt,LepNegPtBins);

  LEPNEGPT_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  LEPNEGPT_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  LEPNEGPT_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  LEPNEGPT_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  LEPNEGPT_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  LEPNEGPT_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  LEPNEGPT_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  LEPNEGPT_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  LEPNEGPT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  LEPNEGPT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEPNEGPT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  LEPNEGPT_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  LEPNEGPT_TOT_UNCERT_BAND_FEWZ->SetFillColor(fillcolorFEWZ);
  LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP->SetFillColor(fillcolorFEWZ);

  //linear scale

  CPlot plotZeeLepNegPt("zeeLepNegPt","",xlabel,ylabel);
  plotZeeLepNegPt.AddHist1D(LEPNEGPT_HIST_DUMMY);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLepNegPt.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLepNegPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLepNegPt.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLepNegPt.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLepNegPt.SetLogx();
  plotZeeLepNegPt.SetLogy(0);
  plotZeeLepNegPt.SetYRange(0.01,1.15*(hUnfoldLepNegPt->GetMaximum() + sqrt(hUnfoldLepNegPt->GetMaximum())));
  plotZeeLepNegPt.SetLegend(0.6,0.64,0.95,0.87);
  plotZeeLepNegPt.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeLepNegPt2("zeeLepNegPtlog","",xlabel,ylabel);
  plotZeeLepNegPt2.AddHist1D(LEPNEGPT_HIST_DUMMY);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLepNegPt2.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLepNegPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLepNegPt2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLepNegPt2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLepNegPt2.SetLogx();
  plotZeeLepNegPt2.SetLogy();
  plotZeeLepNegPt2.SetYRange(5e-5*(hUnfoldLepNegPt->GetMaximum()),100*(hUnfoldLepNegPt->GetMaximum()));
  plotZeeLepNegPt2.SetLegend(0.6,0.64,0.95,0.87);
  plotZeeLepNegPt2.Draw(c1,kTRUE,format);
  
  // AMCATNLO POWHEG MADG FEWZ
  // draw split
  TH1D *hZeeLepNegPtDiffDummySplit =  makeDiffHist(LEPNEGPT_HIST_DUMMY,LEPNEGPT_HIST_DUMMY,"hZeeLepNegPtDiffDummySplit");
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetTitleFont(43);
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetTitleSize(27);
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetTitleOffset(2);
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetLabelFont(43);
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetLabelSize(25);
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetNdivisions(408);
  hZeeLepNegPtDiffDummySplit->GetXaxis()->SetTitleFont(43);
  hZeeLepNegPtDiffDummySplit->GetXaxis()->SetTitleSize(35);
  hZeeLepNegPtDiffDummySplit->GetXaxis()->SetTitleOffset(3.2);
  hZeeLepNegPtDiffDummySplit->GetXaxis()->SetLabelFont(43);
  hZeeLepNegPtDiffDummySplit->GetXaxis()->SetLabelSize(30);
  hZeeLepNegPtDiffDummySplit->GetXaxis()->SetTickSize(.1);
    
  CPlot plotZeeLepNegPtDiffSplit_AMCATNLO("zeeLepNegPtSplit","","","aMC@NLO/Data");
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddHist1D( hZeeLepNegPtDiffDummySplit);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddGraph(LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddGraph(LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddGraph(LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.51,0.465,0.66,0);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.69,0.55,0.90,0.7,0);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddTextBox(lumitext,0.69,0.721,0.93,0.867,0);
  
  plotZeeLepNegPtDiffSplit_AMCATNLO.SetLogx();
  plotZeeLepNegPtDiffSplit_AMCATNLO.SetYRange(0.8+eps,1.2 + (1.2-0.8) *(firstextra)/ base );
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddLine(0, 1,300, 1,kBlack,1);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZeeLepNegPtDiffSplit_AMCATNLO.AddLine(0,0.9,300,0.9,kBlack,3);
  
  hZeeLepNegPtDiffDummySplit->GetYaxis()->SetNdivisions(405);
  
  CPlot plotZeeLepNegPtDiffSplit_POWHEG("zeeLepNegPtSplit","","","POWHEG/Data");
  plotZeeLepNegPtDiffSplit_POWHEG.AddHist1D( hZeeLepNegPtDiffDummySplit);
  plotZeeLepNegPtDiffSplit_POWHEG.AddGraph(LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP,"0E2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepNegPtDiffSplit_POWHEG.AddGraph(LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepNegPtDiffSplit_POWHEG.AddGraph(LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeLepNegPtDiffSplit_POWHEG.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeLepNegPtDiffSplit_POWHEG.SetLogx();
  plotZeeLepNegPtDiffSplit_POWHEG.SetYRange(0.8+eps,1.2-eps);
  plotZeeLepNegPtDiffSplit_POWHEG.AddLine(0, 1,300, 1,kBlack,1);
  plotZeeLepNegPtDiffSplit_POWHEG.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZeeLepNegPtDiffSplit_POWHEG.AddLine(0,0.9,300,0.9,kBlack,3);
  
  CPlot plotZeeLepNegPtDiffSplit_FEWZ("zeeLepNegPtSplit","",xlabel,"FEWZ/Data");
  plotZeeLepNegPtDiffSplit_FEWZ.AddHist1D(hZeeLepNegPtDiffDummySplit);
  plotZeeLepNegPtDiffSplit_FEWZ.AddGraph(LEPNEGPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP,"0E2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepNegPtDiffSplit_FEWZ.AddGraph(LEPNEGPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP,"0P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepNegPtDiffSplit_FEWZ.AddGraph(LEPNEGPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeLepNegPtDiffSplit_FEWZ.AddGraph(LEPNEGPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeLepNegPtDiffSplit_FEWZ.SetLogx();
  plotZeeLepNegPtDiffSplit_FEWZ.SetYRange(0.8+eps,1.2-eps);
  plotZeeLepNegPtDiffSplit_FEWZ.AddLine(0, 1,300, 1,kBlack,1);
  plotZeeLepNegPtDiffSplit_FEWZ.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZeeLepNegPtDiffSplit_FEWZ.AddLine(0,0.9,300,0.9,kBlack,3);
  
  plotZeeLepNegPtDiffSplit_FEWZ.Draw(cR[3],kFALSE,format,1);
  plotZeeLepNegPtDiffSplit_POWHEG.Draw(cR[3],kFALSE,format,2);
  plotZeeLepNegPtDiffSplit_AMCATNLO.Draw(cR[3],kTRUE,format,3);

  //
  // LepPos Pt
  //   

  sprintf(xlabel,"p_{T}^{e^{+}} [GeV]");
  sprintf(ylabel,"d#sigma/dp_{T} [pb/GeV]");

  TH1D *LEPPOSPT_HIST_DUMMY = new TH1D("LEPPOSPT_HIST_DUMMY", "LEPPOSPT_HIST_DUMMY",nBinsLepPosPt,LepPosPtBins);

  LEPPOSPT_HIST_DUMMY->GetYaxis()->SetTitleFont(43);
  LEPPOSPT_HIST_DUMMY->GetYaxis()->SetTitleSize(35);
  LEPPOSPT_HIST_DUMMY->GetYaxis()->SetLabelFont(43);
  LEPPOSPT_HIST_DUMMY->GetYaxis()->SetLabelSize(30);
  LEPPOSPT_HIST_DUMMY->GetXaxis()->SetTitleFont(43);
  LEPPOSPT_HIST_DUMMY->GetXaxis()->SetTitleSize(35);
  LEPPOSPT_HIST_DUMMY->GetXaxis()->SetLabelFont(43);
  LEPPOSPT_HIST_DUMMY->GetXaxis()->SetLabelSize(30);
  LEPPOSPT_HIST_DUMMY->GetYaxis()->SetTitleOffset(1.4);
  LEPPOSPT_TOT_UNCERT_BAND_DATA->SetFillStyle(3554);
  LEPPOSPT_TOT_UNCERT_BAND_DATA->SetFillColor(TColor::GetColor("#828282"));
  LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO->SetFillColor(fillcolorAMCAtNlo);
  LEPPOSPT_TOT_UNCERT_BAND_POWHEG->SetFillColor(fillcolorPowheg);
  LEPPOSPT_TOT_UNCERT_BAND_FEWZ->SetFillColor(fillcolorFEWZ);
  LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP->SetFillStyle(3554);
  LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP->SetFillColor(fillcolorAMCAtNlo);
  LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP->SetFillColor(fillcolorPowheg);
  LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP->SetFillColor(fillcolorFEWZ);

  //linear scale

  CPlot plotZeeLepPosPt("zeeLepPosPt","",xlabel,ylabel);
  plotZeeLepPosPt.AddHist1D(LEPPOSPT_HIST_DUMMY);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLepPosPt.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLepPosPt.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLepPosPt.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLepPosPt.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLepPosPt.SetLogx();
  plotZeeLepPosPt.SetLogy(0);
  plotZeeLepPosPt.SetYRange(0.01,1.15*(hUnfoldLepPosPt->GetMaximum() + sqrt(hUnfoldLepPosPt->GetMaximum())));
  plotZeeLepPosPt.SetLegend(0.6,0.64,0.95,0.87);
  plotZeeLepPosPt.Draw(c1,kTRUE,format);

  //log scale

  CPlot plotZeeLepPosPt2("zeeLepPosPtlog","",xlabel,ylabel);
  plotZeeLepPosPt2.AddHist1D(LEPPOSPT_HIST_DUMMY);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"Data","PE2",1,20,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_AMCATNLO,"aMC@NLO","PE2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_AMCATNLO,"P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  //plotZeeLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_POWHEG,"POWHEG","PE2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_POWHEG,"POWHEG","P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_FEWZ,"FEWZ","PE2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_FEWZ,"P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_TOT_UNCERT_BAND_DATA,"PE2",1,20,1);
  plotZeeLepPosPt2.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA,"P",1,20,1);
  plotZeeLepPosPt2.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.80,0.465,0.85,0);
  plotZeeLepPosPt2.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.18,0.15,0.4,0.25,0);
  plotZeeLepPosPt2.AddTextBox(lumitext,0.69,0.90,0.93,0.94,0);
  plotZeeLepPosPt2.SetLogx();
  plotZeeLepPosPt2.SetLogy();
  plotZeeLepPosPt2.SetYRange(5e-5*(hUnfoldLepPosPt->GetMaximum()),100*(hUnfoldLepPosPt->GetMaximum()));
  plotZeeLepPosPt2.SetLegend(0.6,0.64,0.95,0.87);
  plotZeeLepPosPt2.Draw(c1,kTRUE,format);
  
  // AMCATNLO POWHEG MADG FEWZ
  // draw split
  TH1D *hZeeLepPosPtDiffDummySplit =  makeDiffHist(LEPPOSPT_HIST_DUMMY,LEPPOSPT_HIST_DUMMY,"hZeeLepPosPtDiffDummySplit");
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetTitleFont(43);
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetTitleSize(27);
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetTitleOffset(2);
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetLabelFont(43);
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetLabelSize(25);
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetNdivisions(408);
  hZeeLepPosPtDiffDummySplit->GetXaxis()->SetTitleFont(43);
  hZeeLepPosPtDiffDummySplit->GetXaxis()->SetTitleSize(35);
  hZeeLepPosPtDiffDummySplit->GetXaxis()->SetTitleOffset(3.2);
  hZeeLepPosPtDiffDummySplit->GetXaxis()->SetLabelFont(43);
  hZeeLepPosPtDiffDummySplit->GetXaxis()->SetLabelSize(30);
  hZeeLepPosPtDiffDummySplit->GetXaxis()->SetTickSize(.1);
    
  CPlot plotZeeLepPosPtDiffSplit_AMCATNLO("zeeLepPosPtSplit","","","aMC@NLO/Data");
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddHist1D( hZeeLepPosPtDiffDummySplit);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddGraph(LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_AMCATNLO_COMP,"0E2",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddGraph(LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_AMCATNLO_COMP,"0P",linecolorAMCAtNlo,markerstyleAMCAtNLO,1);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddGraph(LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.205,0.51,0.465,0.66,0);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddTextBox("#splitline{|#eta|<2.4, p_{T}>25 GeV}{PDF set: NNPDF3.0}",0.505,0.55,0.71,0.7,0);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddTextBox(lumitext,0.69,0.721,0.93,0.867,0);
  
  plotZeeLepPosPtDiffSplit_AMCATNLO.SetLogx();
  plotZeeLepPosPtDiffSplit_AMCATNLO.SetYRange(0.8+eps,1.2 + (1.2-0.8) *(firstextra)/ base );
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddLine(0, 1,300, 1,kBlack,1);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZeeLepPosPtDiffSplit_AMCATNLO.AddLine(0,0.9,300,0.9,kBlack,3);
  
  hZeeLepPosPtDiffDummySplit->GetYaxis()->SetNdivisions(405);
    
  CPlot plotZeeLepPosPtDiffSplit_POWHEG("zeeLepPosPtSplit","","","POWHEG/Data");
  plotZeeLepPosPtDiffSplit_POWHEG.AddHist1D( hZeeLepPosPtDiffDummySplit);
  plotZeeLepPosPtDiffSplit_POWHEG.AddGraph(LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_POWHEG_COMP,"0E2",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepPosPtDiffSplit_POWHEG.AddGraph(LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_POWHEG_COMP,"0P",linecolorPowheg,markerstylePowheg,1);
  plotZeeLepPosPtDiffSplit_POWHEG.AddGraph(LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeLepPosPtDiffSplit_POWHEG.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeLepPosPtDiffSplit_POWHEG.SetLogx();
  plotZeeLepPosPtDiffSplit_POWHEG.SetYRange(0.8+eps,1.2-eps);
  plotZeeLepPosPtDiffSplit_POWHEG.AddLine(0, 1,300, 1,kBlack,1);
  plotZeeLepPosPtDiffSplit_POWHEG.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZeeLepPosPtDiffSplit_POWHEG.AddLine(0,0.9,300,0.9,kBlack,3);
  
  CPlot plotZeeLepPosPtDiffSplit_FEWZ("zeeLepPosPtSplit","",xlabel,"FEWZ/Data");
  plotZeeLepPosPtDiffSplit_FEWZ.AddHist1D(hZeeLepPosPtDiffDummySplit);
  plotZeeLepPosPtDiffSplit_FEWZ.AddGraph(LEPPOSPT_RATIO_STAT_SYS_UNCERT_BAND_DATA_FEWZ_COMP,"0E2",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepPosPtDiffSplit_FEWZ.AddGraph(LEPPOSPT_RATIO_STAT_UNCERT_BAND_DATA_FEWZ_COMP,"0P",linecolorFEWZ,markerstyleFEWZ,1);
  plotZeeLepPosPtDiffSplit_FEWZ.AddGraph(LEPPOSPT_STAT_SYS_UNCERT_BAND_DATA_COMP,"E2",TColor::GetColor("#828282"),0,1);
  plotZeeLepPosPtDiffSplit_FEWZ.AddGraph(LEPPOSPT_STAT_UNCERT_BAND_DATA_COMP,"z",1,0,1);
  plotZeeLepPosPtDiffSplit_FEWZ.SetLogx();
  plotZeeLepPosPtDiffSplit_FEWZ.SetYRange(0.8+eps,1.2-eps);
  plotZeeLepPosPtDiffSplit_FEWZ.AddLine(0, 1,300, 1,kBlack,1);
  plotZeeLepPosPtDiffSplit_FEWZ.AddLine(0, 1.1,300, 1.1,kBlack,3);
  plotZeeLepPosPtDiffSplit_FEWZ.AddLine(0,0.9,300,0.9,kBlack,3);
  
  plotZeeLepPosPtDiffSplit_FEWZ.Draw(cR[3],kFALSE,format,1);
  plotZeeLepPosPtDiffSplit_POWHEG.Draw(cR[3],kFALSE,format,2);
  plotZeeLepPosPtDiffSplit_AMCATNLO.Draw(cR[3],kTRUE,format,3);


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

  gBenchmark->Show("plotZee");
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
