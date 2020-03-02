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

#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
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

// Bacon
#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

// Custom tools from this package
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"              // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/AppEffSF.cc"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

void fillMETs(bool doMET,TH1D** h,vector<double> met, int nMET, double wgtLum,double mtCorr);
void fillWeights(bool doMET,TH1D** h,double met, int nWeight,vector<double> wgtLum,double mtCorr);

void fillLHE(TH1D** hlhe, double met, double evtweight, vector<double> *lheweight);
void calcLHE(TH1D* hQCD, TH1D* hPDF, TH1D** hlhe, TH1D* hMain, bool isSignal);
void drawLHE(TH1D** hlhe, TH1D* hMain, TString name, TString outdir, bool isSignal);

void drawShapes(TH1D** vars, TH1D* hMain, TString outdir, TString name, vector<string> leg, int max);
void makeUncMT(vector<Double_t> &metVars, vector<Double_t> &metVarsPhi, TLorentzVector* lep);
// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

void makeDataHistPdf(string dh, string hp, TH1D* hIn, vector<RooDataHist*> &vDataHist, vector<RooHistPdf*> &vHistPdf, RooRealVar &x, int it,string sfx);

void drawWMetPlots(string plotname, TH1D *diff, RooRealVar &x, RooDataHist* dat, RooAddPdf* pdf, RooHistPdf* ewk, RooAbsPdf* qcd, RooHistPdf* wsigp, string lumitext, TH1D* hData);

void drawWMetPlotsSplit(string plotname, TH1D *diff, RooRealVar &x, RooDataHist* dat, RooAddPdf* pdf, RooHistPdf* wx, RooHistPdf* zxx, RooHistPdf* dib, RooHistPdf* ttb, RooAbsPdf* qcd, RooHistPdf* wsigp, string lumitext, TH1D* hData);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
               const Double_t ksprob, const Double_t ksprobpe);

// make webpage
void makeHTML(const TString outDir);

// check a flat Iso cut
bool pass2015Iso(double relIso, double eta){
  if(abs(eta)<1.4442){
    if(relIso >= 0.0354) return false;
    // if(relIso >= 0.02) return false;
  } else {
    if(relIso >= 0.0646) return false;
    // if(relIso >= 0.02) return false;
  }
  return true;
}

// some trash
void setRelIsoVarCR(double &relIso, double pT,double eta){
  if(abs(eta)<1.4442){
    // if(pT < 30) {
      // if(relIso > 0.0354) relIso = 1.0;
    // } else {
      // if(relIso > 0.05) relIso = 1.0;
    // }
    // relIso += 0.10-0.0354;
    // relIso*=10;
    // relIso += 0.10-0.035;
  } else {
    // if(pT < 30) {
      // if(relIso > 0.0646) relIso = 1.0;
    // } else {
      // if(relIso > 0.1) relIso = 1.0;
    // }
    // relIso += 0.10-0.0646;
    relIso -= 0.0646-0.0354;
    // relIso*=10;
    // relIso += 0.10-0.060;
  }
}

// global variables yolo
const  Int_t linecolorW   = kOrange-3;
const  Int_t fillcolorW   = kOrange-2;
const  Int_t linecolorEWK = kOrange+10;
const  Int_t fillcolorEWK = kOrange+7;
const  Int_t linecolorQCD = kViolet+2;
const  Int_t fillcolorQCD = kViolet-5;
const  Int_t ratioColor   = kGray+2;

const Int_t    NBINS   = 50;
// const Int_t    NBINS   = 75;
// const Int_t    NBINS   = 20;
// const Double_t METMIN  = 50;
// const Double_t METMIN  = 0;
// const Double_t METMAX  = 100;
const Double_t METMIN  = 40;
const Double_t METMAX  = 140;
const Int_t    nPDF = 100;
const Int_t    nQCD = 6;

//=== MAIN MACRO ================================================================================================= 

void fitWlnu(const TString  outputDir,   // output directory 
                const TString ntupleDir,
                const TString flav,
                const Double_t lumi,        // integrated luminosity (/fb)'
                const Double_t lumi2 // lumi for the anti-isolation trigger
) {
  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================     

  // file format for output plots
  const TString format("png"); 
  
  bool doMTCut = false;
  // use this to change if we use mT or MET
  bool overwriteMET = true;
  bool doTemplate = true;
  // don't change this part
  bool doMET = true;

double isoSigCut=9999;
double isoTrkCut=9999;

// double isoSigCut=0.05;
// double isoTrkCut=0.05;

  double yscale=0.5;
  
    // Control the types of uncertainties
  enum{pdf,uqcd};
  const string vaLHE[]={"pdf","qcd"};
  // const vector<string> vLHE{"pdf","qcd"};
  int nLHE = sizeof(vaLHE)/sizeof(vaLHE[0]);
  std::cout << "nlhe" << nLHE << std::endl;
    std::vector<string> vLHE;
  for(int i = 0; i < nLHE; i++){vLHE.push_back(vaLHE[i]);}
  
  // main set
  enum{no,cent,eta,keys,ru,rd,stat0,stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8,stat9};
  const string vaMET[]={"no","main","eta","keys","ru","rd","stat0","stat1","stat2","stat3","stat4","stat5","stat6","stat7","stat8","stat9"};
  int nMET = sizeof(vaMET)/sizeof(vaMET[0]);
  std::vector<string> vMET;
  for(int i = 0; i < nMET; i++){vMET.push_back(vaMET[i]);}
  // int ns=nMET-nNV;
  // front half should be nMET-nNV
  
  enum{main,mc,fsr,bkg,tagpt,effstat,pfireu,pfired};
  const string vaWeight[]={"eff","mc","fsr","bkg","tagpt","effstat","pfireu","pfired"};
  int nWeight = sizeof(vaWeight)/sizeof(vaWeight[0]);
  std::vector<string> vWeight;
  for(int i = 0; i < nWeight; i++){vWeight.push_back(vaWeight[i]);}
  
  std::cout << "size of weight array is " << nWeight << std::endl;
  std::cout << "size of met array is " << nMET << std::endl;
  
  Double_t vIsoBins[] = {0.0,0.15,0.30,0.45,0.55,0.65};
  int nIsoBins = sizeof(vIsoBins)/sizeof(vIsoBins[0])-1;
  std::cout << "size of isobin array is " << nIsoBins << std::endl;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t MT_CUT = 40.0;
  const Double_t mu_MASS = 0.1057;
  
  TH2D *hErr  = new TH2D("hErr", "",10,0,10,20,0,20);
  
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  // -----------------------------------------------------

   
  //
  // input ntuple file names
  //
  enum {eData, eWlnu, eZxx, eWx, eTtb, eDib, eQCD, eAntiData, eAntiWlnu, eAntiQCD, eAntiTtb, eAntiDib, eAntiWx, eAntiZxx };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  if(flav.CompareTo("Wenu") == 0 && lumi < 250){
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"));  typev.push_back(eData);
   // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wep_powheg_select.root"));  typev.push_back(eWlnu);
   // fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wem_powheg_select.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we0_select.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we1_select.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we2_select.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx0_select.root"));  typev.push_back(eWx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx1_select.root"));  typev.push_back(eWx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx2_select.root"));  typev.push_back(eWx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eZxx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.root"));  typev.push_back(eDib);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.root"));  typev.push_back(eDib);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.root"));  typev.push_back(eDib);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top1_select.root")); typev.push_back(eTtb);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top2_select.root")); typev.push_back(eTtb);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top3_select.root")); typev.push_back(eTtb);

   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/data_select.root")); typev.push_back(eAntiData);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx0_select.root")); typev.push_back(eAntiWx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx1_select.root")); typev.push_back(eAntiWx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx2_select.root")); typev.push_back(eAntiWx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eAntiZxx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/ww_select.root")); typev.push_back(eAntiDib);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wz_select.root")); typev.push_back(eAntiDib);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zz_select.root")); typev.push_back(eAntiDib);
   // fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wep_powheg_select.root")); typev.push_back(eAntiWlnu);
   // fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wem_powheg_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we0_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we1_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we2_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top1_select.root"));  typev.push_back(eAntiTtb);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top2_select.root"));  typev.push_back(eAntiTtb);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top3_select.root"));  typev.push_back(eAntiTtb);
  }
  // // // // 13 TEV Muon Channel
  if(flav.CompareTo("Wmunu") == 0 && lumi < 250){
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"));    typev.push_back(eData);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm0_select.raw.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm1_select.raw.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm2_select.raw.root"));  typev.push_back(eWlnu);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx0_select.raw.root"));  typev.push_back(eWx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx1_select.raw.root"));  typev.push_back(eWx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx2_select.raw.root"));  typev.push_back(eWx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.raw.root")); typev.push_back(eZxx);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.raw.root"));  typev.push_back(eDib);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.raw.root"));  typev.push_back(eDib);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.raw.root"));  typev.push_back(eDib);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top1_select.raw.root")); typev.push_back(eTtb);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top2_select.raw.root")); typev.push_back(eTtb);
   fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top3_select.raw.root")); typev.push_back(eTtb);

   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/data_select.root")); typev.push_back(eAntiData);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx0_select.root")); typev.push_back(eAntiWx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx1_select.root")); typev.push_back(eAntiWx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx2_select.root")); typev.push_back(eAntiWx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eAntiZxx);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/ww_select.root")); typev.push_back(eAntiDib);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wz_select.root")); typev.push_back(eAntiDib);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zz_select.root")); typev.push_back(eAntiDib);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm0_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm1_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm2_select.root")); typev.push_back(eAntiWlnu);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top1_select.root"));  typev.push_back(eAntiTtb);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top2_select.root"));  typev.push_back(eAntiTtb);
   fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top3_select.root"));  typev.push_back(eAntiTtb);
  }
 
  // // // For the 5 TeV
  if(flav.CompareTo("Wmunu") == 0 && lumi > 250){
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"));    typev.push_back(eData);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wm_select.raw.root"));  typev.push_back(eWlnu);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx_select.raw.root"));  typev.push_back(eWx);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.raw.root")); typev.push_back(eZxx);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.raw.root"));  typev.push_back(eDib);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.raw.root"));  typev.push_back(eDib);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.raw.root"));  typev.push_back(eDib);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top_select.raw.root")); typev.push_back(eTtb);

    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/data_select.root")); typev.push_back(eAntiData);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx_select.root")); typev.push_back(eAntiWx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eAntiZxx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/ww_select.root")); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wz_select.root")); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zz_select.root")); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wm_select.root")); typev.push_back(eAntiWlnu);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top_select.root"));  typev.push_back(eAntiTtb);
  }
    // // // For the 5 TeV
  if(flav.CompareTo("Wenu") == 0 && lumi > 250){
    cout << "reading the ntuples for 5 TeV Wenu " << endl;
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/data_select.root"));    typev.push_back(eData);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/we_select.root"));  typev.push_back(eWlnu);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wx_select.root"));  typev.push_back(eWx);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eZxx);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/zz_select.root"));  typev.push_back(eDib);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/ww_select.root"));  typev.push_back(eDib);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/wz_select.root"));  typev.push_back(eDib);
    fnamev.push_back(ntupleDir+TString("/")+flav+TString("/ntuples/top_select.root")); typev.push_back(eTtb);

    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/data_select.root")); typev.push_back(eAntiData);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wx_select.root")); typev.push_back(eAntiWx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zxx_select.root")); typev.push_back(eAntiZxx);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/ww_select.root")); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/wz_select.root")); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/zz_select.root")); typev.push_back(eAntiDib);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/we_select.root")); typev.push_back(eAntiWlnu);
    fnamev.push_back(ntupleDir+TString("/Anti")+flav+TString("/ntuples/top_select.root"));  typev.push_back(eAntiTtb);
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  

  // TH1D *hGenWeights;
  //
  // Declare MET histograms
  //
  TH1D *hDataMet   = new TH1D("hDataMet","",  NBINS,METMIN,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm  = new TH1D("hDataMetm","", NBINS,METMIN,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp  = new TH1D("hDataMetp","", NBINS,METMIN,METMAX); hDataMetp->Sumw2();
  TH1D *hWlnuMet  = new TH1D("hWlnuMet","", NBINS,METMIN,METMAX); hWlnuMet->Sumw2();
  TH1D *hWlnuMetp = new TH1D("hWlnuMetp","",NBINS,METMIN,METMAX); hWlnuMetp->Sumw2();
  TH1D *hWlnuMetm = new TH1D("hWlnuMetm","",NBINS,METMIN,METMAX); hWlnuMetm->Sumw2();
  
  TH1D *hQCDMet    = new TH1D("hQCDMet", "",  NBINS,METMIN,METMAX); hQCDMet->Sumw2();
  TH1D *hQCDMetp   = new TH1D("hQCDMetp", "", NBINS,METMIN,METMAX); hQCDMetp->Sumw2();
  TH1D *hQCDMetm   = new TH1D("hQCDMetm", "", NBINS,METMIN,METMAX); hQCDMetm->Sumw2();
  
  TH1D *hEWKMet    = new TH1D("hEWKMet", "",  NBINS,METMIN,METMAX); hEWKMet->Sumw2();
  TH1D *hEWKMetp   = new TH1D("hEWKMetp", "", NBINS,METMIN,METMAX); hEWKMetp->Sumw2();
  TH1D *hEWKMetm   = new TH1D("hEWKMetm", "", NBINS,METMIN,METMAX); hEWKMetm->Sumw2();
  
  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet", "", NBINS,METMIN,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","", NBINS,METMIN,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","", NBINS,METMIN,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWlnuMet   = new TH1D("hAntiWlnuMet", "", NBINS,METMIN,METMAX); hAntiWlnuMet->Sumw2();
  TH1D *hAntiWlnuMetp  = new TH1D("hAntiWlnuMetp","", NBINS,METMIN,METMAX); hAntiWlnuMetp->Sumw2();
  TH1D *hAntiWlnuMetm  = new TH1D("hAntiWlnuMetm","", NBINS,METMIN,METMAX); hAntiWlnuMetm->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet",  "", NBINS,METMIN,METMAX); hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "", NBINS,METMIN,METMAX); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "", NBINS,METMIN,METMAX); hAntiEWKMetm->Sumw2();
  
  TH1D *hAntiQCDMet    = new TH1D("hAntiQCDMet", "",  NBINS,METMIN,METMAX); hAntiQCDMet->Sumw2();
  TH1D *hAntiQCDMetp   = new TH1D("hAntiQCDMetp", "", NBINS,METMIN,METMAX); hAntiQCDMetp->Sumw2();
  TH1D *hAntiQCDMetm   = new TH1D("hAntiQCDMetm", "", NBINS,METMIN,METMAX); hAntiQCDMetm->Sumw2();
  
  TH1D *hMuonEtaDatap = new TH1D("hMuonEtaDatap","",25,0,2.5); hMuonEtaDatap->Sumw2();
  TH1D *hMuonEtaDatam = new TH1D("hMuonEtaDatam","",25,0,2.5); hMuonEtaDatam->Sumw2();
  TH1D *hMuonEtaMCp   = new TH1D("hMuonEtaMCp",  "",25,0,2.5); hMuonEtaMCp->Sumw2();
  TH1D *hMuonEtaMCm   = new TH1D("hMuonEtaMCm",  "",25,0,2.5); hMuonEtaMCm->Sumw2();
  TH1D *hMuonEtaAntiDatap = new TH1D("hMuonEtaAntiDatap","",25,0,2.5); hMuonEtaAntiDatap->Sumw2();
  TH1D *hMuonEtaAntiDatam = new TH1D("hMuonEtaAntiDatam","",25,0,2.5); hMuonEtaAntiDatam->Sumw2();
  
  TH1D *hDataMetpPhi   = new TH1D("hDataMetpPhi","",  100,-3.15, 6.30); hDataMetpPhi->Sumw2();
  TH1D *hDataMetmPhi   = new TH1D("hDataMetmPhi","",  100,-3.15, 6.30); hDataMetmPhi->Sumw2();
  TH1D *hWlnuMetpPhi   = new TH1D("hWlnuMetpPhi","",  100,-3.15, 6.30); hWlnuMetpPhi->Sumw2();
  TH1D *hWlnuMetmPhi   = new TH1D("hWlnuMetmPhi","",  100,-3.15, 6.30); hWlnuMetmPhi->Sumw2();
  

  // All teh uncertainties? 
  TH1D ***hWlnupMETU  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnumMETU  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  TH1D ***hWlnupMETD  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnumMETD  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  
  TH1D ***hEWKpMETU  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKmMETU  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  TH1D ***hEWKpMETD  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKmMETD  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D ***hWxpMETU  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxmMETU  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  TH1D ***hWxpMETD  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxmMETD  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  
  TH1D ***hZxxpMETU  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxmMETU  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  TH1D ***hZxxpMETD  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxmMETD  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D ***hWlnupWeightU  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnumWeightU  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  TH1D ***hWlnupWeightD  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnumWeightD  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  
  TH1D ***hEWKpWeightU  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKmWeightU  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  TH1D ***hEWKpWeightD  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKmWeightD  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D ***hWxpWeightU  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxmWeightU  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  TH1D ***hWxpWeightD  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxmWeightD  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  
  TH1D ***hZxxpWeightU  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxmWeightU  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  TH1D ***hZxxpWeightD  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxmWeightD  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  
  TH1D ***hDibpWeightU  = new TH1D**[nIsoBins];// hAntiDibMetp->Sumw2();
  TH1D ***hDibmWeightU  = new TH1D**[nIsoBins];// hAntiDibMetm->Sumw2();
  TH1D ***hDibpWeightD  = new TH1D**[nIsoBins];// hAntiDibMetp->Sumw2();
  TH1D ***hDibmWeightD  = new TH1D**[nIsoBins];// hAntiDibMetm->Sumw2();
  
  TH1D ***hTtbpWeightU  = new TH1D**[nIsoBins];// hAntiTtbMetp->Sumw2();
  TH1D ***hTtbmWeightU  = new TH1D**[nIsoBins];// hAntiTtbMetm->Sumw2();
  TH1D ***hTtbpWeightD  = new TH1D**[nIsoBins];// hAntiTtbMetp->Sumw2();
  TH1D ***hTtbmWeightD  = new TH1D**[nIsoBins];// hAntiTtbMetm->Sumw2();
  
  // shit about to get outta hand
  TH1D ***hWlnupLHE = new TH1D**[nIsoBins];
  TH1D ***hEWKpLHE  = new TH1D**[nIsoBins];
  TH1D ***hWxpLHE   = new TH1D**[nIsoBins];
  TH1D ***hZxxpLHE  = new TH1D**[nIsoBins];
  TH1D ***hDibpLHE  = new TH1D**[nIsoBins];
  TH1D ***hTtbpLHE  = new TH1D**[nIsoBins];
  
  TH1D ***hWlnumLHE = new TH1D**[nIsoBins];
  TH1D ***hEWKmLHE  = new TH1D**[nIsoBins];
  TH1D ***hWxmLHE   = new TH1D**[nIsoBins];
  TH1D ***hZxxmLHE  = new TH1D**[nIsoBins];
  TH1D ***hDibmLHE  = new TH1D**[nIsoBins];
  TH1D ***hTtbmLHE  = new TH1D**[nIsoBins];
  
  
    ////////////////////////////////////////////////////////////////////////////////////////
  TH1D ***hWlnupThyUncU  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnumThyUncU  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  TH1D ***hWlnupThyUncD  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnumThyUncD  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  
  TH1D ***hEWKpThyUncU  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKmThyUncU  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  TH1D ***hEWKpThyUncD  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKmThyUncD  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D ***hWxpThyUncU  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxmThyUncU  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  TH1D ***hWxpThyUncD  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxmThyUncD  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  
  TH1D ***hZxxpThyUncU  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxmThyUncU  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  TH1D ***hZxxpThyUncD  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxmThyUncD  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  
  TH1D ***hDibpThyUncU  = new TH1D**[nIsoBins];// hAntiDibMetp->Sumw2();
  TH1D ***hDibmThyUncU  = new TH1D**[nIsoBins];// hAntiDibMetm->Sumw2();
  TH1D ***hDibpThyUncD  = new TH1D**[nIsoBins];// hAntiDibMetp->Sumw2();
  TH1D ***hDibmThyUncD  = new TH1D**[nIsoBins];// hAntiDibMetm->Sumw2();
  
  TH1D ***hTtbpThyUncU  = new TH1D**[nIsoBins];// hAntiTtbMetp->Sumw2();
  TH1D ***hTtbmThyUncU  = new TH1D**[nIsoBins];// hAntiTtbMetm->Sumw2();
  TH1D ***hTtbpThyUncD  = new TH1D**[nIsoBins];// hAntiTtbMetp->Sumw2();
  TH1D ***hTtbmThyUncD  = new TH1D**[nIsoBins];// hAntiTtbMetm->Sumw2();
  
  // do a loop to create the second array
  for(int i=0; i < nIsoBins; ++i){
    // w signal recoil
    hWlnupMETU[i] = new TH1D*[nMET];
    hWlnumMETU[i] = new TH1D*[nMET];
    
    hWlnupMETD[i] = new TH1D*[nMET];
    hWlnumMETD[i] = new TH1D*[nMET];
    
    // w signal efficiency 
    hWlnupWeightU[i] = new TH1D*[nWeight];
    hWlnumWeightU[i] = new TH1D*[nWeight];
    
    hWlnupWeightD[i] = new TH1D*[nWeight];
    hWlnumWeightD[i] = new TH1D*[nWeight];
    
    hWlnupThyUncU[i] = new TH1D*[nLHE];
    hWlnumThyUncU[i] = new TH1D*[nLHE];
    
    hWlnupThyUncD[i] = new TH1D*[nLHE];
    hWlnumThyUncD[i] = new TH1D*[nLHE];
    
    // ewk total recoil
    hEWKpMETU[i] = new TH1D*[nMET];
    hEWKmMETU[i] = new TH1D*[nMET];
    
    hEWKpMETD[i] = new TH1D*[nMET];
    hEWKmMETD[i] = new TH1D*[nMET];
    
    // ewk total efficiency
    hEWKpWeightU[i] = new TH1D*[nWeight];
    hEWKmWeightU[i] = new TH1D*[nWeight];
    
    hEWKpWeightD[i] = new TH1D*[nWeight];
    hEWKmWeightD[i] = new TH1D*[nWeight];
    
    hEWKpThyUncU[i] = new TH1D*[nLHE];
    hEWKmThyUncU[i] = new TH1D*[nLHE];
    
    hEWKpThyUncD[i] = new TH1D*[nLHE];
    hEWKmThyUncD[i] = new TH1D*[nLHE];
    
    // wx  recoil
    hWxpMETU[i] = new TH1D*[nMET];
    hWxmMETU[i] = new TH1D*[nMET];
    
    hWxpMETD[i] = new TH1D*[nMET];
    hWxmMETD[i] = new TH1D*[nMET];
    
    // wx  efficiency
    hWxpWeightU[i] = new TH1D*[nWeight];
    hWxmWeightU[i] = new TH1D*[nWeight];
    
    hWxpWeightD[i] = new TH1D*[nWeight];
    hWxmWeightD[i] = new TH1D*[nWeight];
    
    // wx  efficiency
    hWxpThyUncU[i] = new TH1D*[nLHE];
    hWxmThyUncU[i] = new TH1D*[nLHE];
    
    hWxpThyUncD[i] = new TH1D*[nLHE];
    hWxmThyUncD[i] = new TH1D*[nLHE];
    
    
    // zxx recoil
    hZxxpMETU[i] = new TH1D*[nMET];
    hZxxmMETU[i] = new TH1D*[nMET];
    
    hZxxpMETD[i] = new TH1D*[nMET];
    hZxxmMETD[i] = new TH1D*[nMET];
    
    // zxx efficiency
    hZxxpWeightU[i] = new TH1D*[nWeight];
    hZxxmWeightU[i] = new TH1D*[nWeight];
    
    hZxxpWeightD[i] = new TH1D*[nWeight];
    hZxxmWeightD[i] = new TH1D*[nWeight];
    
    hZxxpThyUncU[i] = new TH1D*[nLHE];
    hZxxmThyUncU[i] = new TH1D*[nLHE];
    
    hZxxpThyUncD[i] = new TH1D*[nLHE];
    hZxxmThyUncD[i] = new TH1D*[nLHE];
    
    // diboson efficiency
    hDibpWeightU[i] = new TH1D*[nWeight];
    hDibmWeightU[i] = new TH1D*[nWeight];
    
    hDibpWeightD[i] = new TH1D*[nWeight];
    hDibmWeightD[i] = new TH1D*[nWeight];
    
    // ttbar efficiency
    hTtbpWeightU[i] = new TH1D*[nWeight];
    hTtbmWeightU[i] = new TH1D*[nWeight];
    
    hTtbpWeightD[i] = new TH1D*[nWeight];
    hTtbmWeightD[i] = new TH1D*[nWeight];
    
    hWlnupLHE[i] = new TH1D*[nQCD+nPDF];
    hEWKpLHE[i]  = new TH1D*[nQCD+nPDF];
    hWxpLHE[i]   = new TH1D*[nQCD+nPDF];
    hZxxpLHE[i]  = new TH1D*[nQCD+nPDF];
    hDibpLHE[i]  = new TH1D*[nQCD+nPDF];
    hTtbpLHE[i]  = new TH1D*[nQCD+nPDF];
  
    hWlnumLHE[i] = new TH1D*[nQCD+nPDF];
    hEWKmLHE[i]  = new TH1D*[nQCD+nPDF];
    hWxmLHE[i]   = new TH1D*[nQCD+nPDF];
    hZxxmLHE[i]  = new TH1D*[nQCD+nPDF];
    hDibmLHE[i]  = new TH1D*[nQCD+nPDF];
    hTtbmLHE[i]  = new TH1D*[nQCD+nPDF];
    
  }
  
  // -----------------------------------------------------------------------
  //           the main shapes  
  // -----------------------------------------------------------------------
  TH1D **hDataMetm2d  = new TH1D*[nIsoBins];// hAntiDataMetm->Sumw2();  
  TH1D **hDataMetp2d  = new TH1D*[nIsoBins];// hAntiDataMetp->Sumw2();
  
  TH1D **hWlnuMetp2d  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D **hWlnuMetm2d  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  
  TH1D **hEWKMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hEWKMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hDibMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hDibMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hTtbMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hTtbMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hWxMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hWxMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hZxxMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hZxxMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();

  TH1D **hQCDMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hQCDMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  
  TH1D **hMetpIsoValues = new TH1D*[nIsoBins];
  TH1D **hMetmIsoValues = new TH1D*[nIsoBins];
  // Create a histogram pointer in each space in the array
  for(int i = 0; i < nIsoBins; i++){
    hDataMetp2d[i]  = new TH1D(("hDataMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDataMetm2d[i]  = new TH1D(("hDataMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hWlnuMetp2d[i]  = new TH1D(("hWlnuMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWlnuMetm2d[i]  = new TH1D(("hWlnuMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hEWKMetp2d[i]      = new TH1D(("hEwkMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hEWKMetm2d[i]      = new TH1D(("hEwkMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hDibMetp2d[i]  = new TH1D(("hDibMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDibMetm2d[i]  = new TH1D(("hDibMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hTtbMetp2d[i]  = new TH1D(("hTtbMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hTtbMetm2d[i]  = new TH1D(("hTtbMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
	
    hWxMetp2d[i]      = new TH1D(("hWxMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hWxMetm2d[i]      = new TH1D(("hWxMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hZxxMetp2d[i]      = new TH1D(("hZxxMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hZxxMetm2d[i]      = new TH1D(("hZxxMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    hQCDMetp2d[i]      = new TH1D(("hQCDMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hQCDMetm2d[i]      = new TH1D(("hQCDMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
    // Create a loop over the # of uncertainty shapes to produce the uncertainty histograms for the "up" shapes
    // Do the recoil ones here
    for(int j = 0; j < nMET; ++j){
      char hname[150]; char type[50];
      if(j==ru){
        sprintf(type,"%i_rochUp",i);
      } else if (j==rd){
        sprintf(type,"%i_rochDown",i);
      } else {
        sprintf(type,"%i_%sUp",i,(vMET[j]).c_str());
      }

      // Wlnu
      sprintf(hname,"hWlnuMetpBin%s",type);
      hWlnupMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWlnuMetmBin%s",type);
      hWlnumMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      // ewk sum
      sprintf(hname,"hEWKMetpBin%s",type);
      hEWKpMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hEWKMetmBin%s",type);
      hEWKmMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      // wx
      sprintf(hname,"hWxMetpBin%s",type);
      hWxpMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWxMetmBin%s",type);
      hWxmMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      // zxx
      sprintf(hname,"hZxxMetpBin%s",type);
      hZxxpMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hZxxMetmBin%s",type);
      hZxxmMETU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      // }
    }
    
    for(int j=0; j < nWeight; ++j){
      char hname[150]; char type[50];
      //sdu,sdd,smu,smd,pfireu,pfired
      if (j==pfireu){
        sprintf(type,"%i_prefireUp",i);
      } else if (j==pfired){
        sprintf(type,"%i_prefireDown",i);
      } else {
        sprintf(type,"%i_%sUp",i,(vWeight[j]).c_str());
      }
      
      sprintf(hname,"hWlnuMetpBin%s",type);
      hWlnupWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWlnuMetmBin%s",type);
      hWlnumWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      sprintf(hname,"hEwkMetpBin%s",type);
      hEWKpWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hEwkMetmBin%s",type);
      hEWKmWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      sprintf(hname,"hWxMetpBin%s",type);
      hWxpWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hWxMetmBin%s",type);
      hWxmWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
  
      sprintf(hname,"hZxxMetpBin%s",type);
      hZxxpWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hZxxMetmBin%s",type);
      hZxxmWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
      sprintf(hname,"hDibMetpBin%s",type);
      hDibpWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hDibMetmBin%s",type);
      hDibmWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
  
      sprintf(hname,"hTtbMetpBin%s",type);
      hTtbpWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      sprintf(hname,"hTtbMetmBin%s",type);
      hTtbmWeightU[i][j] = new TH1D(hname,"",NBINS,METMIN,METMAX);
      
    }
    
    for(int j = 0; j < nQCD + nPDF; ++j){ // the 100ish LHE weight histos
      hWlnupLHE[i][j] = new TH1D(("hWlnuMetpBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      hWlnumLHE[i][j] = new TH1D(("hWlnuMetmBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      
      hEWKpLHE[i][j] = new TH1D(("hEwkMetpBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      hEWKmLHE[i][j] = new TH1D(("hEwkMetmBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      
      hWxpLHE[i][j] = new TH1D(("hWxMetpBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      hWxmLHE[i][j] = new TH1D(("hWxMetmBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
  
      hZxxpLHE[i][j] = new TH1D(("hZxxMetpBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      hZxxmLHE[i][j] = new TH1D(("hZxxMetmBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      
      hDibpLHE[i][j] = new TH1D(("hDibMetpBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      hDibmLHE[i][j] = new TH1D(("hDibMetmBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
  
      hTtbpLHE[i][j] = new TH1D(("hTtbMetpBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
      hTtbmLHE[i][j] = new TH1D(("hTtbMetmBin"+std::to_string(i)+"_lhe"+std::to_string(j)).c_str(),"",NBINS,METMIN,METMAX);
    }
    
    for(int j=0; j < nLHE; ++j){ // the final up/down shapes for the QCD and PDF uncertainty
      // w signal
      hWlnupThyUncU[i][j] = new TH1D(("hWlnuMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hWlnumThyUncU[i][j] = new TH1D(("hWlnuMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // ewk total
      hEWKpThyUncU[i][j] = new TH1D(("hEwkMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hEWKmThyUncU[i][j] = new TH1D(("hEwkMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // wx
      hWxpThyUncU[i][j] = new TH1D(("hWxMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hWxmThyUncU[i][j] = new TH1D(("hWxMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // zxx
      hZxxpThyUncU[i][j] = new TH1D(("hZxxMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hZxxmThyUncU[i][j] = new TH1D(("hZxxMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
    }
    

    hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
    hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
  }
  

  // Some test plots: 
  // pfiso/pt vs rel comb iso in barrel
  // in endcap
  TH2D *hBarrelIsoPosData = new TH2D("hBarrelIsoPosData","hBarrelIsoPosData",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoPosData = new TH2D("hEndcapIsoPosData","hEndcapIsoPosData",100,25,100,100,0,0.50);
  
  TH2D *hBarrelIsoPosWsig = new TH2D("hBarrelIsoPosWsig","hBarrelIsoPosWsig",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoPosWsig = new TH2D("hEndcapIsoPosWsig","hEndcapIsoPosWsig",100,25,100,100,0,0.50);
  
  TH2D *hBarrelIsoNegData = new TH2D("hBarrelIsoNegData","hBarrelIsoNegData",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoNegData = new TH2D("hEndcapIsoNegData","hEndcapIsoNegData",100,25,100,100,0,0.50);
  
  TH2D *hBarrelIsoNegWsig = new TH2D("hBarrelIsoNegWsig","hBarrelIsoNegWsig",100,25,100,100,0,0.50);
  TH2D *hEndcapIsoNegWsig = new TH2D("hEndcapIsoNegWsig","hEndcapIsoNegWsig",100,25,100,100,0,0.50);

  double tolerance = ROOT::Math::MinimizerOptions::DefaultTolerance();
  string algo = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
  string type = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
  int strategy= ROOT::Math::MinimizerOptions::DefaultStrategy();

  int precision= ROOT::Math::MinimizerOptions::DefaultPrecision();
  int MaxFunctionCalls= ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls();
  int MaxIterations= ROOT::Math::MinimizerOptions::DefaultMaxIterations();

  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum, npv, npu;
  Float_t genVPt, genVPhi, genVy, genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight, prefirePhoton, prefireJet;
  Float_t prefireUp, prefireDown;
  Float_t met, metPhi, sumEt, mt, u1, u2, deta;
  Int_t   q;
  TLorentzVector *lep=0, *lep_raw=0, *genV=0, *genLep=0;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso, trkIso;
  Double_t mtCorr;
  Double_t effSFweight=1, relIso;


  vector<Double_t>  *metVars=0, *metVarsPhi=0;
  vector<Double_t>  *evtWeight=0;
  vector<Double_t>  *lheweight=0;
    
  
  TFile *infile=0;
  TTree *intree=0;

  double noPrefire_Wp = 0, prefire_Wp=0, prefireJet_Wp=0, prefirePhoton_Wp=0;
  double noPrefire_Wm = 0, prefire_Wm=0, prefireJet_Wm=0, prefirePhoton_Wm=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);
    std::cout << "type = " << typev[ifile] << std::endl;

    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genVy",    &genVy);   // GEN W boson phi (signal MC)
    intree->SetBranchAddress("genLepPt",   &genLepPt);    // GEN lepton pT (signal MC)
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);   // GEN lepton phi (signal MC)
    intree->SetBranchAddress("prefireWeight",  &prefireWeight);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefireUp",      &prefireUp);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefireDown",    &prefireDown);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefirePhoton",  &prefirePhoton);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("prefireJet",    &prefireJet);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("scale1fb",      &scale1fb);  // MCwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",    &scale1fbUp);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",  &scale1fbDown);  // eventwgtLum[main] per 1/fb (MC)
    intree->SetBranchAddress("sumEt",         &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",            &mt);        // transverse mass
    intree->SetBranchAddress("mtCorr",        &mtCorr);        // transverse mass
    intree->SetBranchAddress("q",             &q);         // lepton charge
    intree->SetBranchAddress("lep",           &lep);       // lepton 4-vector
    intree->SetBranchAddress("lep_raw",       &lep_raw);       // lepton 4-vector
    intree->SetBranchAddress("genLep",        &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV",          &genV);       // lepton 4-vector
    intree->SetBranchAddress("trkIso",        &trkIso);
    intree->SetBranchAddress("pfChIso",       &pfChIso);
    intree->SetBranchAddress("pfGamIso",      &pfGamIso);
    intree->SetBranchAddress("pfNeuIso",      &pfNeuIso);
    intree->SetBranchAddress("pfCombIso",     &pfCombIso);       // lepton 4-vector
    intree->SetBranchAddress("relIso",        &relIso);       // relative isolation for the lepton
    intree->SetBranchAddress("evtWeight",     &evtWeight); // eventwgtLum[main] vector
    intree->SetBranchAddress("metVars",       &metVars);            // contains the different met variations
    intree->SetBranchAddress("metVarsPhi",    &metVarsPhi);         // met phi for variations
    intree->SetBranchAddress("lheweight",     &lheweight);         // pdf and qcdwgtLum[main]s
    // intree->SetBranchAddress("deta",        &deta);         // pdf and qcdwgtLum[main]s
  
    TH1D* hGenWeights; double totalNorm = 1.0;
    // cout << "Hello " << endl;
    // if(typev[ifile] != eData && typev[ifile] != eAntiData && !(flav.CompareTo("Wmunu") == 0)){
    if(typev[ifile] != eData && typev[ifile] != eAntiData){
      // cout << "get gen weights" << endl;
      hGenWeights = (TH1D*)infile->Get("hGenWeights");
      totalNorm = hGenWeights->Integral();
      // cout << totalNorm << endl;
    }
  
    UInt_t iterator=10;
    //
    // loop over events
    //
    // double frac=0.1;
    // if(typev[ifile]==eAntiData || typev[ifile]==eAntiWlnu ||typev[ifile]==eAntiWx|| typev[ifile]==eAntiZxx || typev[ifile]==eAntiDib || typev[ifile]==eAntiTtb) {
      // // frac = 1.0;
      // frac = 0.2;
    // }
    // if(typev[ifile]==eTtb) frac=1;
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<((UInt_t)intree->GetEntries()); ientry+=iterator) {
      intree->GetEntry(ientry);
      if(ientry%1000000==0)  cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;


      if(lep->Pt() < PT_CUT) continue;//std::cout << " pass PT " << std::endl;
      if(fabs(lep->Eta()) > ETA_CUT) continue;//std::cout << " pass eta " << std::endl;
      Double_t corr=1, corrdu=1, corrdd=1, corrmu=1, corrmd=1;
      Double_t corrFSR=1;
      Double_t corrMC=1;
      Double_t corrBkg=1;
      Double_t corrTag=1;
      vector<double> wgtLum;

      // accidentally set up prefire uncertainties in a way that a very rare event gets nan for the up/down variation
      for(int jt=0; jt < nWeight; jt++) {
        wgtLum.push_back(lumi*(isnan((*evtWeight)[jt]) ? (*evtWeight)[main] : (*evtWeight)[jt])/totalNorm);
      }
      
      if(overwriteMET && typev[ifile] != eData && typev[ifile] != eAntiData){
        makeUncMT(*metVars, *metVarsPhi, lep);
      } else if (overwriteMET) (*metVars)[no] = mtCorr;

      if(doMTCut&&(mtCorr<MT_CUT)) continue;//std::cout << " pass mt " << std::endl;
      if(typev[ifile]==eData) {
        if(relIso > isoSigCut) continue;
        if(trkIso > isoTrkCut) continue;
        
        hDataMet->Fill((*metVars)[no]);
        if(q>0) {
          doMET ? hDataMetp->Fill((*metVars)[no]) : hDataMetp->Fill(mtCorr);
          hDataMetpPhi->Fill((*metVarsPhi)[no]);
          hMuonEtaDatap->Fill(fabs(lep->Eta()));
          doMET ? hDataMetp2d[0]->Fill((*metVars)[no]) : hDataMetp2d[0]->Fill(mtCorr);
          doMET ? hQCDMetp2d[0]->Fill((*metVars)[no]) : hQCDMetp2d[0]->Fill(mtCorr);
          hMetpIsoValues[0]->Fill(relIso);
        } else {
          doMET ? hDataMetm->Fill((*metVars)[no]) : hDataMetm->Fill(mtCorr);
          hMuonEtaDatam->Fill(fabs(lep->Eta()));
          hDataMetmPhi->Fill((*metVarsPhi)[no]);
          doMET ? hDataMetm2d[0]->Fill((*metVars)[no]) : hDataMetm2d[0]->Fill(mtCorr);
          doMET ? hQCDMetm2d[0]->Fill((*metVars)[no]) : hQCDMetm2d[0]->Fill(mtCorr); 
          hMetmIsoValues[0]->Fill(relIso);
        }
      } else if(typev[ifile]==eAntiData) {
        if(abs(lep->Eta())<1.4442){
          if(q > 0)hBarrelIsoPosData->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hBarrelIsoNegData->Fill(lep->Pt(),relIso);
        } else {
          if(q > 0)hEndcapIsoPosData->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hEndcapIsoNegData->Fill(lep->Pt(),relIso);
        }
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            hAntiDataMet->Fill((*metVars)[no]);
            if(q>0) { 
              doMET ? hAntiDataMetp->Fill((*metVars)[no]) : hAntiDataMetp->Fill(mtCorr);
              hMuonEtaAntiDatap->Fill(fabs(lep->Eta()));
              doMET ? hDataMetp2d[it]->Fill((*metVars)[no]) : hDataMetp2d[it]->Fill(mtCorr);
              doMET ? hQCDMetp2d[it]->Fill((*metVars)[no]) : hQCDMetp2d[it]->Fill(mtCorr);
              hMetpIsoValues[it]->Fill(relIso);
            } else { 
              doMET ? hAntiDataMetm->Fill((*metVars)[no]) : hAntiDataMetm->Fill(mtCorr);
              hMuonEtaAntiDatap->Fill(fabs(lep->Eta()));
              doMET ? hDataMetm2d[it]->Fill((*metVars)[no]) : hDataMetm2d[it]->Fill(mtCorr);
              doMET ? hQCDMetm2d[it]->Fill((*metVars)[no]) : hQCDMetm2d[it]->Fill(mtCorr);
              hMetmIsoValues[it]->Fill(relIso);
              break;
            }
          }
        }
      } else if(typev[ifile]==eWlnu) {
        if(relIso > isoSigCut) continue;
        if(trkIso > isoTrkCut) continue;
        int bin=0;
        if(typev[ifile]==eWlnu){
          // start filling the counters
          if(q > 0){
            noPrefire_Wp+=wgtLum[main]/prefireWeight;
            prefire_Wp+=wgtLum[main];
            prefireJet_Wp+=wgtLum[main]*prefireJet/prefireWeight;
            prefirePhoton_Wp+=wgtLum[main]*prefirePhoton/prefireWeight;
          } else {
            noPrefire_Wm+=wgtLum[main]/prefireWeight;
            prefire_Wm+=wgtLum[main];
            prefireJet_Wm+=wgtLum[main]*prefireJet/prefireWeight;
            prefirePhoton_Wm+=wgtLum[main]*prefirePhoton/prefireWeight;
          }
          
        }
        
        hWlnuMet->Fill((*metVars)[cent],wgtLum[main]);
        if(q>0){
          hMuonEtaMCp->Fill(fabs(lep->Eta()),wgtLum[main]);
          hWlnuMetpPhi->Fill((*metVarsPhi)[cent]);
          doMET ? hWlnuMetp->Fill((*metVars)[cent],wgtLum[main]) : hWlnuMetp->Fill(mtCorr,wgtLum[main]);
          doMET ? hWlnuMetp2d[0]->Fill((*metVars)[cent],wgtLum[main]) : hWlnuMetp2d[0]    ->Fill(mtCorr,wgtLum[main]);
          fillMETs(doMET,hWlnupMETU[0],(*metVars),nMET,wgtLum[main],mtCorr);
          fillWeights(doMET,hWlnupWeightU[0],(*metVars)[cent],nWeight,wgtLum,mtCorr);

        } else {
          hMuonEtaMCm->Fill(fabs(lep->Eta()),wgtLum[main]);
          hWlnuMetmPhi->Fill((*metVarsPhi)[cent]);
          doMET ? hWlnuMetm->Fill((*metVars)[cent],wgtLum[main]) : hWlnuMetm->Fill(mtCorr,wgtLum[main]);
          doMET ? hWlnuMetm2d[0]->Fill((*metVars)[cent],wgtLum[main]) : hWlnuMetm2d[0]->Fill(mtCorr,wgtLum[main]);
          fillMETs(doMET,hWlnumMETU[0],(*metVars),nMET,wgtLum[main],mtCorr);
          fillWeights(doMET,hWlnumWeightU[0],(*metVars)[cent],nWeight,wgtLum,mtCorr);
        }
      } else if(typev[ifile]==eWx) {
        if(relIso > isoSigCut) continue;
        if(trkIso > isoTrkCut) continue;
        doMET ? hEWKMet->Fill((*metVars)[cent],wgtLum[main]) : hEWKMet->Fill(mtCorr,wgtLum[main]);
        if(q>0){
          doMET ? hEWKMetp->Fill((*metVars)[cent],wgtLum[main]) : hEWKMetp->Fill(mtCorr,wgtLum[main]);
          doMET ? hWxMetp2d[0]->Fill((*metVars)[cent],wgtLum[main]) : hWxMetp2d[0]->Fill(mtCorr,wgtLum[main]);
          fillMETs(doMET,hWxpMETU[0],(*metVars),nMET,wgtLum[main],mtCorr);
          fillWeights(doMET,hWxpWeightU[0],(*metVars)[cent],nWeight,wgtLum,mtCorr);
        } else {
          doMET ? hEWKMetm->Fill((*metVars)[cent],wgtLum[main]): hEWKMetm->Fill(mtCorr,wgtLum[main]);
          doMET ? hWxMetm2d[0]     ->Fill((*metVars)[cent],wgtLum[main]) : hWxMetm2d[0] ->Fill(mtCorr,wgtLum[main]);
          fillMETs(doMET,hWxmMETU[0],(*metVars),nMET,wgtLum[main],mtCorr);
          fillWeights(doMET,hWxmWeightU[0],(*metVars)[cent],nWeight,wgtLum,mtCorr);
        }
      } else if(typev[ifile]==eZxx){
        if(relIso > isoSigCut) continue;
        if(trkIso > isoTrkCut) continue;
        doMET ? hEWKMet->Fill((*metVars)[cent],wgtLum[main]) : hEWKMet->Fill(mtCorr,wgtLum[main]);
        if(q>0){
 				  doMET ? hEWKMetp->Fill((*metVars)[cent],wgtLum[main]) : hEWKMetp->Fill(mtCorr,wgtLum[main]);
          doMET ? hZxxMetp2d[0]    ->Fill((*metVars)[cent],wgtLum[main]) : hZxxMetp2d[0]    ->Fill(mtCorr,wgtLum[main]);
          fillMETs(doMET,hZxxpMETU[0],(*metVars),nMET,wgtLum[main],mtCorr);
          fillWeights(doMET,hZxxpWeightU[0],(*metVars)[cent],nWeight,wgtLum,mtCorr);
        } else {
          doMET ? hEWKMetm->Fill((*metVars)[cent],wgtLum[main]) :  hEWKMetm->Fill(mtCorr,wgtLum[main]);
          doMET ? hZxxMetm2d[0]    ->Fill((*metVars)[cent],wgtLum[main]) : hZxxMetm2d[0]    ->Fill(mtCorr,wgtLum[main]);
          fillMETs(doMET,hZxxmMETU[0],(*metVars),nMET,wgtLum[0],mtCorr);
          fillWeights(doMET,hZxxmWeightU[0],(*metVars)[cent],nWeight,wgtLum,mtCorr);
        }
      } else if(typev[ifile]==eDib) {
        if(relIso > isoSigCut) continue;
        if(trkIso > isoTrkCut) continue;
        doMET ? hEWKMet->Fill((*metVars)[no],wgtLum[main]) : hEWKMet->Fill(mtCorr,wgtLum[main]);
        if(q>0){
          doMET ? hEWKMetp->Fill((*metVars)[no],wgtLum[main]) : hEWKMetp->Fill(mtCorr,wgtLum[main]);
          doMET ? hDibMetp2d[0]->Fill((*metVars)[no],wgtLum[main]) : hDibMetp2d[0]->Fill(mtCorr,wgtLum[main]);
          fillWeights(doMET,hDibpWeightU[0],(*metVars)[no],nWeight,wgtLum,mtCorr);
        } else {
          doMET ? hEWKMetm->Fill((*metVars)[no],wgtLum[main]) : hEWKMetm->Fill(mtCorr,wgtLum[main]);
          doMET ? hDibMetm2d[0]->Fill((*metVars)[no],wgtLum[main]) : hDibMetm2d[0]->Fill(mtCorr,wgtLum[main]);
          fillWeights(doMET,hDibmWeightU[0],(*metVars)[no],nWeight,wgtLum,mtCorr);
        }
      } else if(typev[ifile]==eTtb) {
        if(relIso > isoSigCut) continue;
        if(trkIso > isoTrkCut) continue;
        doMET ? hEWKMet->Fill((*metVars)[no],wgtLum[main]) : hEWKMet->Fill(mtCorr,wgtLum[main]);
        if(q>0){
          doMET ? hEWKMetp->Fill((*metVars)[no],wgtLum[main]) : hEWKMetp->Fill(mtCorr,wgtLum[main]);
          doMET ? hTtbMetp2d[0]->Fill((*metVars)[no],wgtLum[main]) : hTtbMetp2d[0]->Fill(mtCorr,wgtLum[main]);
          fillWeights(doMET,hTtbpWeightU[0],(*metVars)[no],nWeight,wgtLum,mtCorr);
        } else {
          doMET ? hEWKMetm->Fill((*metVars)[no],wgtLum[main]) : hEWKMetm->Fill(mtCorr,wgtLum[main]);
          doMET ? hTtbMetm2d[0]->Fill((*metVars)[no],wgtLum[main]) : hTtbMetm2d[0]->Fill(mtCorr,wgtLum[main]);
          fillWeights(doMET,hTtbmWeightU[0],(*metVars)[no],nWeight,wgtLum,mtCorr);
        }
      } else if(typev[ifile]==eAntiWlnu){
        if(abs(lep->Eta())<1.4442){
          if(q > 0)hBarrelIsoPosWsig->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hBarrelIsoNegWsig->Fill(lep->Pt(),relIso);
        } else {
          if(q > 0)hEndcapIsoPosWsig->Fill(lep->Pt(),pfCombIso/lep->Pt());
          else hEndcapIsoNegWsig->Fill(lep->Pt(),relIso);
        }
        hAntiWlnuMet->Fill((*metVars)[no],wgtLum[main]);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0) {              
              doMET ? hAntiWlnuMetp->Fill((*metVars)[no],wgtLum[main]) : hAntiWlnuMetp->Fill(mtCorr,wgtLum[main]);
              doMET ? hWlnuMetp2d[it]    ->Fill((*metVars)[no],wgtLum[main]) : hWlnuMetp2d[it]    ->Fill(mtCorr,wgtLum[main]);
              fillMETs(doMET,hWlnupMETU[it],(*metVars),nMET,wgtLum[main],mtCorr);
              fillWeights(doMET,hWlnupWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            } else {
              doMET ? hAntiWlnuMetm->Fill((*metVars)[no],wgtLum[main]) : hAntiWlnuMetm->Fill(mtCorr,wgtLum[main]);
              doMET ? hWlnuMetm2d[it]    ->Fill((*metVars)[no],wgtLum[main]) : hWlnuMetm2d[it]    ->Fill(mtCorr,wgtLum[main]);
              fillMETs(doMET,hWlnumMETU[it],(*metVars),nMET,wgtLum[main],mtCorr);
              fillWeights(doMET,hWlnumWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            }
          }
        }
      } else if(typev[ifile]==eAntiWx){
        doMET ? hAntiEWKMet->Fill((*metVars)[no],wgtLum[main]) : hAntiEWKMet->Fill(mtCorr,wgtLum[main]);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              doMET ? hAntiEWKMetp->Fill((*metVars)[no],wgtLum[main]) :  hAntiEWKMetp->Fill(mtCorr,wgtLum[main]);
              doMET ? hWxMetp2d[it]     ->Fill((*metVars)[no],wgtLum[main]) : hWxMetp2d[it]     ->Fill(mtCorr,wgtLum[main]);
              fillMETs(doMET,hWxpMETU[it],(*metVars),nMET,wgtLum[main],mtCorr);
              fillWeights(doMET,hWxpWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            } else {
              doMET ? hAntiEWKMetm->Fill((*metVars)[no],wgtLum[main]) : hAntiEWKMetm->Fill(mtCorr,wgtLum[main]);
              doMET ? hWxMetm2d[it]     ->Fill((*metVars)[no],wgtLum[main]) : hWxMetm2d[it]     ->Fill(mtCorr,wgtLum[main]);
              fillMETs(doMET,hWxmMETU[it],(*metVars),nMET,wgtLum[main],mtCorr);
              fillWeights(doMET,hWxmWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            }
            break;
          }
        }
      } else if(typev[ifile]==eAntiZxx){
        doMET ? hAntiEWKMet->Fill((*metVars)[no],wgtLum[main]) : hAntiEWKMet->Fill(mtCorr,wgtLum[main]);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill((*metVars)[no],wgtLum[main]); 
              doMET ? hZxxMetp2d[it]    ->Fill((*metVars)[no],wgtLum[main]) : hZxxMetp2d[it]->Fill(mtCorr,wgtLum[main]);
              fillMETs(doMET,hZxxpMETU[it],(*metVars),nMET,wgtLum[main],mtCorr);
              fillWeights(doMET,hZxxpWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            } else {
              hAntiEWKMetm->Fill((*metVars)[no],wgtLum[main]); 
              doMET ? hZxxMetm2d[it]    ->Fill((*metVars)[no],wgtLum[main]) : hZxxMetm2d[it]    ->Fill(mtCorr,wgtLum[main]);
              fillMETs(doMET,hZxxmMETU[it],(*metVars),nMET,wgtLum[main],mtCorr);
              fillWeights(doMET,hZxxmWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            }
          }
        }
      } else if(typev[ifile]==eAntiDib){
        doMET ? hAntiEWKMet->Fill((*metVars)[no],wgtLum[main]) : hAntiEWKMet->Fill(mtCorr,wgtLum[main]);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill((*metVars)[no],wgtLum[main]); 
              doMET ? hDibMetp2d[it]->Fill((*metVars)[no],wgtLum[main]) : hDibMetp2d[it]->Fill(mtCorr,wgtLum[main]);
              fillWeights(doMET,hDibpWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            } else {
              hAntiEWKMetm->Fill((*metVars)[no],wgtLum[main]); 
              doMET ? hDibMetm2d[it]->Fill((*metVars)[no],wgtLum[main]) : hDibMetm2d[it]->Fill(mtCorr,wgtLum[main]);
              fillWeights(doMET,hDibmWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            }
          }
        }
      } else if(typev[ifile]==eAntiTtb){
        doMET ? hAntiEWKMet->Fill((*metVars)[no],wgtLum[main]) : hAntiEWKMet->Fill(mtCorr,wgtLum[main]);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill((*metVars)[no],wgtLum[main]); 
              doMET ? hTtbMetp2d[it]->Fill((*metVars)[no],wgtLum[main]) : hTtbMetp2d[it]->Fill(mtCorr,wgtLum[main]);
              fillWeights(doMET,hTtbpWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            } else {
              hAntiEWKMetm->Fill((*metVars)[no],wgtLum[main]); 
              doMET ? hTtbMetm2d[it]->Fill((*metVars)[no],wgtLum[main]) : hTtbMetm2d[it]->Fill(mtCorr,wgtLum[main]);
              fillWeights(doMET,hTtbmWeightU[it],(*metVars)[no],nWeight,wgtLum,mtCorr);
            }
          }
        }
      } else if(typev[ifile]==eAntiQCD) {
        hAntiQCDMet->Fill((*metVars)[no],wgtLum[main]);
        q>0 ? hAntiQCDMetp->Fill((*metVars)[no],wgtLum[main]) : hAntiQCDMetm->Fill((*metVars)[no],wgtLum[main]); 
      }
    }
  }
  
  
  // Print the yields in each of the isolation bins
  ofstream txtIso;
  char txtnameIso[100];
  std::cout << "Printing Yields by Isolation Bin: " << std::endl;
  sprintf(txtnameIso,"%s/isoBins_Yields.txt",CPlot::sOutDir.Data());
  txtIso.open(txtnameIso);
  assert(txtIso.is_open());
	txtIso << "W+ " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtIso << "bin " << i << "  Data# "  << hDataMetp2d[i]->Integral() << std::endl;
		txtIso << "bin " << i << "  Wsig# "  << hWlnuMetp2d[i]->Integral() << std::endl;
		txtIso << "bin " << i << "  ewk # "  << hEWKMetp2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "   wx # "  << hWxMetp2d[i]->Integral()   << std::endl;
		txtIso << "bin " << i << "  zxx # "  << hZxxMetp2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  ttb # "  << hTtbMetp2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  dib # "  << hDibMetp2d[i]->Integral()  << std::endl;
	}
	
	txtIso << "W- " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtIso << "bin " << i << "  Data# "  << hDataMetm2d[i]->Integral() << std::endl;
		txtIso << "bin " << i << "  Wsig# "  << hWlnuMetm2d[i]->Integral() << std::endl;
		txtIso << "bin " << i << "  ewk # "  << hEWKMetm2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "   wx # "  << hWxMetm2d[i]->Integral()   << std::endl;
		txtIso << "bin " << i << "  zxx # "  << hZxxMetm2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  ttb # "  << hTtbMetm2d[i]->Integral()  << std::endl;
		txtIso << "bin " << i << "  dib # "  << hDibMetm2d[i]->Integral()  << std::endl;
	}
	txtIso.close();
  
  
  char sname[50];
  for(int it = 0; it < nIsoBins; it++){
    
    // sprintf(sname,"hWlnuMetp2d%i",it);  drawLHE(hWlnupLHE[it], hWlnuMetp2d[it], sname, CPlot::sOutDir.Data(), 1);
    // sprintf(sname,"hWlnuMetm2d%i",it);  drawLHE(hWlnumLHE[it], hWlnuMetm2d[it], sname, CPlot::sOutDir.Data(), 1);
    
    // sprintf(sname,"hEWKMetp2d%i",it);  drawLHE(hEWKpLHE[it], hEWKMetp2d[it], sname, CPlot::sOutDir.Data(), 0);
    // sprintf(sname,"hEWKMetm2d%i",it);  drawLHE(hEWKmLHE[it], hEWKMetm2d[it], sname, CPlot::sOutDir.Data(), 0);
    
    // sprintf(sname,"hZxxMetp2d%i",it);  drawLHE(hZxxpLHE[it], hZxxMetp2d[it], sname, CPlot::sOutDir.Data(), 0);
    // sprintf(sname,"hZxxMetm2d%i",it);  drawLHE(hZxxmLHE[it], hZxxMetm2d[it], sname, CPlot::sOutDir.Data(), 0);
    
    // sprintf(sname,"hWxMetp2d%i",it);  drawLHE(hWxpLHE[it], hWxMetp2d[it], sname, CPlot::sOutDir.Data(), 0);
    // sprintf(sname,"hWxMetm2d%i",it);  drawLHE(hWxmLHE[it], hWxMetm2d[it], sname, CPlot::sOutDir.Data(), 0);
    
    // sprintf(sname,"hDibMetp2d%i",it);  drawLHE(hDibpLHE[it], hDibMetp2d[it], sname, CPlot::sOutDir.Data(), 0);
    // sprintf(sname,"hDibMetm2d%i",it);  drawLHE(hDibmLHE[it], hDibMetm2d[it], sname, CPlot::sOutDir.Data(), 0);
    
    // sprintf(sname,"hTtbMetp2d%i",it);  drawLHE(hTtbpLHE[it], hTtbMetp2d[it], sname, CPlot::sOutDir.Data(), 0);
    // sprintf(sname,"hTtbMetm2d%i",it);  drawLHE(hTtbmLHE[it], hTtbMetm2d[it], sname, CPlot::sOutDir.Data(), 0);
    // std::cout << "calc lhe1" << std::endl;
    // std::cout << hWlnupThyUncU[it][0]->Integral() << std::endl;
    // std::cout << hWlnupThyUncU[it][1]->Integral() << std::endl;
    // calcLHE(hWlnupThyUncU[it][uqcd], hWlnupThyUncU[it][pdf], hWlnupLHE[it], hWlnuMetp2d[it], 1);
    // calcLHE(hWlnumThyUncU[it][uqcd], hWlnumThyUncU[it][pdf], hWlnumLHE[it], hWlnuMetm2d[it], 1);
    // // std::cout << "calc lhe2" << std::endl;
    // calcLHE(hEWKpThyUncU[it][uqcd], hEWKpThyUncU[it][pdf], hEWKpLHE[it], hEWKMetp2d[it], 0);
    // calcLHE(hEWKmThyUncU[it][uqcd], hEWKmThyUncU[it][pdf], hEWKmLHE[it], hEWKMetm2d[it], 0);
    // // std::cout << "calc lhe3" << std::endl;
    // calcLHE(hWxpThyUncU[it][uqcd], hWxpThyUncU[it][pdf], hWxpLHE[it], hWxMetp2d[it], 0);
    // calcLHE(hWxmThyUncU[it][uqcd], hWxmThyUncU[it][pdf], hWxmLHE[it], hWxMetm2d[it], 0);
    // // std::cout << "calc lhe4" << std::endl;
    // calcLHE(hZxxpThyUncU[it][uqcd], hZxxpThyUncU[it][pdf], hZxxpLHE[it], hZxxMetp2d[it], 0);
    // calcLHE(hZxxmThyUncU[it][uqcd], hZxxmThyUncU[it][pdf], hZxxmLHE[it], hZxxMetm2d[it], 0);
    
    // CPlot::sOutDir.Data()
    // calcLHE(TH1D* hQCD, TH1D* hPDF, hTtbmLHE[it], hTtbMetm2d[it], 0);
    // calcLHE(TH1D* hQCD, TH1D* hPDF, hTtbmLHE[it], hTtbMetm2d[it], 0);
    
    hEWKMetp2d[it]->Add(hTtbMetp2d[it],1);
    hEWKMetp2d[it]->Add(hDibMetp2d[it],1);
    hEWKMetp2d[it]->Add(hZxxMetp2d[it],1);
    hEWKMetp2d[it]->Add(hWxMetp2d[it],1);
    
    hEWKMetm2d[it]->Add(hTtbMetm2d[it],1);
    hEWKMetm2d[it]->Add(hDibMetm2d[it],1);
    hEWKMetm2d[it]->Add(hZxxMetm2d[it],1);
    hEWKMetm2d[it]->Add(hWxMetm2d[it],1);

    hQCDMetp2d[it]->Add(hEWKMetp2d[it],-1);
    hQCDMetp2d[it]->Add(hWlnuMetp2d[it],-1);
	
    hQCDMetm2d[it]->Add(hEWKMetm2d[it],-1);
    hQCDMetm2d[it]->Add(hWlnuMetm2d[it],-1);

    for(Int_t ibin=1; ibin<=hQCDMetp2d[it]->GetNbinsX(); ibin++) {
      if(hQCDMetp2d[it]->GetBinContent(ibin)<0) {
        hQCDMetp2d[it]->SetBinContent(ibin,0);
        hQCDMetp2d[it]->SetBinError(ibin,1.8);
      }
      if(hQCDMetm2d[it]->GetBinContent(ibin)<0) {
        hQCDMetm2d[it]->SetBinContent(ibin,0);
        hQCDMetm2d[it]->SetBinError(ibin,1.8);
      }
    }
   
    
    for(int j=0; j< nWeight; j++){
      hEWKpWeightU[it][j]->Add(hTtbpWeightU[it][j],1);
      hEWKpWeightU[it][j]->Add(hDibpWeightU[it][j],1);
      hEWKpWeightU[it][j]->Add(hZxxpWeightU[it][j],1);
      hEWKpWeightU[it][j]->Add(hWxpWeightU[it][j],1);
      
      hEWKmWeightU[it][j]->Add(hTtbmWeightU[it][j],1);
      hEWKmWeightU[it][j]->Add(hDibmWeightU[it][j],1);
      hEWKmWeightU[it][j]->Add(hZxxmWeightU[it][j],1);
      hEWKmWeightU[it][j]->Add(hWxmWeightU[it][j],1);
    }
    for(int j=0; j < nMET; j++){
      hEWKpMETU[it][j]->Add(hZxxpMETU[it][j],1);
      hEWKpMETU[it][j]->Add(hWxpMETU[it][j],1);
      
      hEWKmMETU[it][j]->Add(hZxxmMETU[it][j],1);
      hEWKmMETU[it][j]->Add(hWxmMETU[it][j],1);
    }
    
  }
  // std::cout << "blahasfadsfa " << std::endl;
  for(int it=0; it<1; it++){
    
    for(int i=0; i < nWeight; i++){
      std::cout << hWlnupWeightU[it][i]->Integral() << std::endl;
    }
    for(int i=0; i < nMET; i++){
      std::cout << hWlnupMETU[it][i]->Integral() << std::endl;
    }
    
    drawShapes(hWlnupMETU[it], hWlnuMetp2d[it], CPlot::sOutDir.Data(), "WlnupMET", vMET, nMET);
    drawShapes(hWlnumMETU[it], hWlnuMetm2d[it], CPlot::sOutDir.Data(), "WlnumMET", vMET, nMET);
    
    drawShapes(hWlnupWeightU[it], hWlnuMetp2d[it], CPlot::sOutDir.Data(), "WlnupWeight", vWeight, nWeight);
    drawShapes(hWlnumWeightU[it], hWlnuMetm2d[it], CPlot::sOutDir.Data(), "WlnumWeight", vWeight, nWeight);
    
    
    drawShapes(hEWKpMETU[it], hEWKMetp2d[it], CPlot::sOutDir.Data(), "EWKpMET", vMET, nMET);
    drawShapes(hEWKmMETU[it], hEWKMetm2d[it], CPlot::sOutDir.Data(), "EWKmMET", vMET, nMET);
    
    drawShapes(hWxpMETU[it], hWxMetp2d[it], CPlot::sOutDir.Data(), "WxpMET", vMET, nMET);
    drawShapes(hWxmMETU[it], hWxMetm2d[it], CPlot::sOutDir.Data(), "WxmMET", vMET, nMET);
    
    drawShapes(hZxxpMETU[it], hZxxMetp2d[it], CPlot::sOutDir.Data(), "ZxxpMET", vMET, nMET);
    drawShapes(hZxxmMETU[it], hZxxMetm2d[it], CPlot::sOutDir.Data(), "ZxxmMET", vMET, nMET);
    
    
    drawShapes(hEWKpWeightU[it], hEWKMetp2d[it], CPlot::sOutDir.Data(), "EWKpWeight", vWeight, nWeight);
    drawShapes(hEWKmWeightU[it], hEWKMetm2d[it], CPlot::sOutDir.Data(), "EWKmWeight", vWeight, nWeight);
    
    drawShapes(hWxpWeightU[it], hWxMetp2d[it], CPlot::sOutDir.Data(), "WxpWeight", vWeight, nWeight);
    drawShapes(hWxmWeightU[it], hWxMetm2d[it], CPlot::sOutDir.Data(), "WxmWeight", vWeight, nWeight);
    
    drawShapes(hZxxpWeightU[it], hZxxMetp2d[it], CPlot::sOutDir.Data(), "ZxxpWeight", vWeight, nWeight);
    drawShapes(hZxxmWeightU[it], hZxxMetm2d[it], CPlot::sOutDir.Data(), "ZxxmWeight", vWeight, nWeight);
    
    drawShapes(hDibpWeightU[it], hDibMetp2d[it], CPlot::sOutDir.Data(), "DibpWeight", vWeight, nWeight);
    drawShapes(hDibmWeightU[it], hDibMetm2d[it], CPlot::sOutDir.Data(), "DibmWeight", vWeight, nWeight);
    
    drawShapes(hTtbpWeightU[it], hTtbMetp2d[it], CPlot::sOutDir.Data(), "TtbpWeight", vWeight, nWeight);
    drawShapes(hTtbmWeightU[it], hTtbMetm2d[it], CPlot::sOutDir.Data(), "TtbmWeight", vWeight, nWeight);
  }
  
  delete infile;
  infile=0, intree=0;   
   
  ofstream txtfilePF;
  char txtfname2[150];
  std::cout << "Printing prefiring  values" << std::endl;
  sprintf(txtfname2,"%s/prefire_yields.txt",CPlot::sOutDir.Data());
  txtfilePF.open(txtfname2);
  assert(txtfilePF.is_open());
  txtfilePF << "---- W+ ----- " << endl;
  txtfilePF << "No PF: " << noPrefire_Wp << endl;
  txtfilePF << "PF: " << prefire_Wp << "  scale factor: " << noPrefire_Wp/prefire_Wp << endl;
  txtfilePF << "PF (jet only): " << prefireJet_Wp << "  scale factor: " << noPrefire_Wp/prefireJet_Wp << endl;
  txtfilePF << "PF (photon only): " << prefirePhoton_Wp << "  scale factor: " << noPrefire_Wp/prefirePhoton_Wp << endl << endl;
  txtfilePF << "---- W- ----- " << endl;
  txtfilePF << "No PF: " << noPrefire_Wm << endl;
  txtfilePF << "PF: " << prefire_Wm << "  scale factor: " << noPrefire_Wm/prefire_Wm << endl;
  txtfilePF << "PF (jet only): " << prefireJet_Wm << "  scale factor: " << noPrefire_Wm/prefireJet_Wm << endl;
  txtfilePF << "PF (photon only): " << prefirePhoton_Wm << "  scale factor: " << noPrefire_Wm/prefirePhoton_Wm << endl << endl;
  txtfilePF.close();
   
    ofstream txtfile2;
    // char txtfname2[100];
    std::cout << "Printing We+. " << std::endl;
    sprintf(txtfname2,"%s/isoAvg.txt",CPlot::sOutDir.Data());
    txtfile2.open(txtfname2);
    assert(txtfile2.is_open());
 
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
    j==0 ? c2->SetLogy(1) : c2->SetLogy(0);
	  hMetpIsoValues[j]->Draw();
	  c2->Update();c2->SaveAs(plotname);
	  sprintf(plotname,"%s/hMetmIso_%d.png",CPlot::sOutDir.Data(),j);
    c2->Clear();
	  hMetmIsoValues[j]->Draw();c2->Update();c2->SaveAs(plotname);
  }
  
    txtfile2.close();
    sprintf(plotname,"%s/hMetpPhi_SR.png",CPlot::sOutDir.Data());
    c2->Clear();
    double norm = hDataMetpPhi->Integral()/hWlnuMetpPhi->Integral();
    hWlnuMetpPhi->Scale(norm);
    hWlnuMetpPhi->Draw();
    hDataMetpPhi->Draw("same");
	  c2->Update();c2->SaveAs(plotname);
    sprintf(plotname,"%s/hMetmPhi_SR.png",CPlot::sOutDir.Data());
    c2->Clear();
    norm = hDataMetmPhi->Integral()/hWlnuMetmPhi->Integral();
    hWlnuMetmPhi->Scale(norm);
    hWlnuMetmPhi->Draw();
    hDataMetmPhi->Draw("same");
	  c2->Update();c2->SaveAs(plotname);
    
  
  
  std::cout << " ==== nEvents === " << std::endl;
  std::cout << "sig single  = " << hWlnuMetm->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWlnuMetm2d[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetm->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetm2d[0]->Integral() << std::endl;
  
  std::cout << "sig single  = " << hWlnuMetp->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWlnuMetp2d[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetp->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetp2d[0]->Integral() << std::endl;
  
  
 
 // std::cout << "blah" << std::endl;
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig("nSig","nSig",0.7*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar nQCD("nQCD","nQCD",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar cewk("cewk","cewk",0.1,0,5) ;
  cewk.setVal(hEWKMet->Integral()/hWlnuMet->Integral());
//   cewk.setConstant(kTRUE);
  RooFormulaVar nEWK("nEWK","nEWK","cewk*nSig",RooArgList(nSig,cewk));
  RooRealVar nAntiSig("nAntiSig","nAntiSig",hAntiWlnuMet->Integral()*0.9,0,hAntiDataMet->Integral());
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",0.9*(hDataMet->Integral()),0,hDataMet->Integral());
  RooRealVar dewk("dewk","dewk",0.1,0,5) ;
  dewk.setVal(hAntiEWKMet->Integral()/hAntiWlnuMet->Integral());
  dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK("nAntiEWK","nAntiEWK","dewk*nAntiSig",RooArgList(nAntiSig,dewk));
  
  RooRealVar nSigp("nSigp","nSigp",0.7*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWlnuMetp->Integral(),0,hAntiDataMetp->Integral());
//   RooRealVar nSigp("nSigp","nSigp",90000,0,hDataMetp->Integral());
  //RooRealVar nQCDp("nQCDp","nQCDp",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
    
  RooRealVar nQCDp("nQCDp","nQCDp",hDataMetp->Integral()*0.3,0,hDataMetp->Integral());
  RooRealVar cewkp("cewkp","cewkp",0.1,0,5) ;
  cewkp.setVal(hEWKMetp->Integral()/hWlnuMetp->Integral());
//   cewkp.setConstant(kTRUE);
  RooFormulaVar nEWKp("nEWKp","nEWKp","cewkp*nSigp",RooArgList(nSigp,cewkp));
  //RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",hAntiWlnuMetp->Integral()*1.0,0,hAntiDataMetp->Integral());
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar dewkp("dewkp","dewkp",0.1,0,5) ;
  dewkp.setVal(hAntiEWKMetp->Integral()/hAntiWlnuMetp->Integral());
  dewkp.setConstant(kTRUE);
  RooFormulaVar nAntiEWKp("nAntiEWKp","nAntiEWKp","dewkp*nAntiSigp",RooArgList(nAntiSigp,dewkp));
  
  RooRealVar nSigm("nSigm","nSigm",0.7*(hDataMetm->Integral()),0,hDataMetm->Integral());
//   RooRealVar nSigm("nSigm","nSigm",75000,0,hDataMetm->Integral());
  //RooRealVar nQCDm("nQCDm","nQCDm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  RooRealVar nQCDm("nQCDm","nQCDm",hDataMetm->Integral()*0.3,0,hDataMetm->Integral());
  RooRealVar cewkm("cewkm","cewkm",0.1,0,5) ;
  cewkm.setVal(hEWKMetm->Integral()/hWlnuMetm->Integral());
//   cewkm.setConstant(kTRUE);
  RooFormulaVar nEWKm("nEWKm","nEWKm","cewkm*nSigm",RooArgList(nSigm,cewkm));  
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",hAntiWlnuMetm->Integral()*1.0,0,hAntiDataMetm->Integral());
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  RooRealVar dewkm("dewkm","dewkm",0.1,0,5) ;
  dewkm.setVal(hAntiEWKMetm->Integral()/hAntiWlnuMetm->Integral());
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
   // RooRealVar* nEWKpSR = new RooRealVar("nEWKpSR","nEWKpSR",hEWKMetp2d[0]->Integral());
   // RooRealVar* nEWKmSR = new RooRealVar("nEWKmSR","nEWKmSR",hEWKMetm2d[0]->Integral());
  char nname[50];
  char formula[50];
  for (int j = 0; j < nIsoBins; ++j){
	  // W+
	  double qcdFac = 0.9;
      if(j==0) qcdFac = 0.15;
      sprintf(nname, "nAntiSigp%d",j);
      nSigp_[j] = new RooRealVar(nname,nname,hWlnuMetp2d[j]->Integral(),0,hDataMetp2d[j]->Integral());      
      // sprintf(nname, "nAntiQCDp%d",j);
      // nQCDp_[j] = new RooRealVar(nname,nname,(qcdFac)*hDataMetp2d[j]->Integral(),0,hDataMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiTtbp%d",j);
      nTtbp_[j] = new RooRealVar(nname,nname,hTtbMetp2d[j]->Integral(),0.8*hTtbMetp2d[j]->Integral(),hDataMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiDibp%d",j);
      nDibp_[j] = new RooRealVar(nname,nname,hDibMetp2d[j]->Integral(),0.8*hDibMetp2d[j]->Integral(),hDataMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiZxxp%d",j);
      nZxxp_[j] = new RooRealVar(nname,nname,hZxxMetp2d[j]->Integral(),0.8*hZxxMetp2d[j]->Integral(),hDataMetp2d[j]->Integral());
	  
	  double qcd_init = hDataMetp2d[j]->Integral() - hWlnuMetp2d[j]->Integral() - hEWKMetp2d[j]->Integral() ;
	  if(j==0){
         sprintf(nname, "nAntiQCDp%d",j);
		 // nQCDp_[j] = new RooRealVar(nname,nname,0.4*hDataMetp2d[j]->Integral(),0.0*hDataMetp2d[j]->Integral(),0.6*hDataMetp2d[j]->Integral());
         nQCDp_[j] = new RooRealVar(nname,nname,qcd_init,0,2.0*qcd_init);
         // nQCDp_[j]->setConstant(kTRUE);
	  } else {
         sprintf(nname, "nAntiQCDp%d",j);
         // nQCDp_[j] = new RooRealVar(nname,nname,hDataMetp2d[j]->Integral(),0.95*hDataMetp2d[j]->Integral(),1.05*hDataMetp2d[j]->Integral());
         nQCDp_[j] = new RooRealVar(nname,nname,hDataMetp2d[j]->Integral(),0.95*hDataMetp2d[j]->Integral(),1.05*hDataMetp2d[j]->Integral());
	     // nQCDp_[j]->setConstant(kTRUE);
      }
	  // nQCDp_[j]->setConstant(kTRUE);
	  
	  	  // the ratio between W->enu and W->tau nu
	  sprintf(nname, "dwxp%d",j);
      dwxp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dwxp_[j]->setVal(hWxMetp2d[j]->Integral()/hWlnuMetp2d[j]->Integral());
	  
	  sprintf(nname, "nAntiWxp%d",j); sprintf(formula,"dwxp%d*nAntiSigp%d",j,j);
	  nWxp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*dwxp_[j],*nSigp_[j]));
	  
      // sprintf(nname, "dewkp%d",j);
      // dewkp_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkp_[j]->setVal(hEWKMetp2d[j]->Integral()/hWlnuMetp2d[j]->Integral());
	  
      sprintf(nname, "nAntiEWKp%d",j); sprintf(formula,"dewkp%d*nAntiSigp%d",j,j);
      // // nEWKp_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigp_[j],*dewkp_[j]));
	  nEWKp_[j] = new RooRealVar(nname,nname,hEWKMetp2d[j]->Integral(),0.0,(2.0)*hEWKMetp2d[j]->Integral());
	
	  // sprintf(nname, "nAntiDibp%d",j); 
	  // nDibp_[j] = new RooRealVar(nname,nname,hDibMetp2d[j]->Integral(),0.0,(2.0)*hDibMetp2d[j]->Integral());	  
	  
	  // sprintf(nname, "nAntiTtbp%d",j); 
	  // nTtbp_[j] = new RooRealVar(nname,nname,hTtbMetp2d[j]->Integral(),0.0,(2.0)*hTtbMetp2d[j]->Integral());
	
	 
      // W-	 
	  sprintf(nname, "nAntiSigm%d",j);
      nSigm_[j] = new RooRealVar(nname,nname,hWlnuMetm2d[j]->Integral(),0,hDataMetm2d[j]->Integral());
	  
	  	  
	  sprintf(nname, "nAntiTtbm%d",j);
      nTtbm_[j] = new RooRealVar(nname,nname,hTtbMetm2d[j]->Integral(),0.8*hTtbMetm2d[j]->Integral(),hDataMetm2d[j]->Integral());
	  
	  sprintf(nname, "nAntiDibm%d",j);
      nDibm_[j] = new RooRealVar(nname,nname,hDibMetm2d[j]->Integral(),0.8*hDibMetm2d[j]->Integral(),hDataMetm2d[j]->Integral());
	  
	  sprintf(nname, "nAntiZxxm%d",j);
      nZxxm_[j] = new RooRealVar(nname,nname,hZxxMetm2d[j]->Integral(),0.8*hZxxMetm2d[j]->Integral(),hDataMetm2d[j]->Integral());

	  qcd_init = hDataMetm2d[j]->Integral() - hWlnuMetm2d[j]->Integral() - hEWKMetm2d[j]->Integral() ;
    // j==0?qcd_init = hDataMetm2d[j]->Integral() - hWlnuMetm2d[j]->Integral() - hEWKMetm2d[j]->Integral(): hDataMetm2d[j]->Integral();
	  if(j==0){
         sprintf(nname, "nAntiQCDm%d",j);
         nQCDm_[j] = new RooRealVar(nname,nname,qcd_init,0,qcd_init*2.0);
	  } else {
		 sprintf(nname, "nAntiQCDm%d",j);
         nQCDm_[j] = new RooRealVar(nname,nname,hDataMetm2d[j]->Integral(),0.95*hDataMetm2d[j]->Integral(),1.05*hDataMetm2d[j]->Integral());
      }
	  
    // the W->tau nu / W->e nu ratio
	  sprintf(nname, "dwxm%d",j);
      dwxm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      dwxm_[j]->setVal(hWxMetm2d[j]->Integral()/hWlnuMetm2d[j]->Integral());
	  
	  sprintf(nname, "nAntiWxm%d",j); sprintf(formula,"dwxm%d*nAntiSigm%d",j,j);
	  nWxm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*dwxm_[j],*nSigm_[j]));
	  
      // sprintf(nname, "dewkm%d",j);
      // dewkm_[j] = new RooRealVar(nname,nname,0.1,0,10.0);
      // dewkm_[j]->setVal(hEWKMetm2d[j]->Integral()/hWlnuMetm2d[j]->Integral());
	  
      sprintf(nname, "nAntiEWKm%d",j); sprintf(formula,"dewkm%d*nAntiSigm%d",j,j);
      // // nEWKm_[j] = new RooFormulaVar(nname,nname,formula,RooArgList(*nSigm_[j],*dewkm_[j]));
	  nEWKm_[j] = new RooRealVar(nname,nname,hEWKMetm2d[j]->Integral(),0.0,(2.0)*hEWKMetm2d[j]->Integral());
	  
	  
	  std::cout << "blah 3." << j << std::endl;
	  
  }
  
   std::cout << "finished normalizations" << std::endl;
  
  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet("pfmet","pfmet",0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist WlnuMet ("WlnuMET", "WlnuMET", RooArgSet(pfmet),hWlnuMet);  RooHistPdf pdfWm ("wm", "wm", pfmet,WlnuMet, 1);
  RooDataHist WlnuMetp("WlnuMETp","WlnuMETp",RooArgSet(pfmet),hWlnuMetp); RooHistPdf pdfWmp("wmp","wmp",pfmet,WlnuMetp,1);
  RooDataHist WlnuMetm("WlnuMETm","WlnuMETm",RooArgSet(pfmet),hWlnuMetm); RooHistPdf pdfWmm("wmm","wmm",pfmet,WlnuMetm,1); 
  
  std::cout << "blalalala1" << std::endl;
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", pfmet,ewkMet, 1);
  std::cout << "blalalala1bbb" << std::endl;
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",pfmet,ewkMetp,1); 
  std::cout << "blalalalaaaa1" << std::endl;
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",pfmet,ewkMetm,1); 
  
  std::cout << "blalalala2" << std::endl;
//   // test using the reversed dEta, dPhi cuts as background
  RooDataHist qcdMet ("qcdMET", "qcdMET", RooArgSet(pfmet),hQCDMet);  RooHistPdf pdfQCD ("qcd", "qcd", pfmet,qcdMet, 1);
  RooDataHist qcdMetp("qcdMetp","qcdMetp",RooArgSet(pfmet),hQCDMetp); RooHistPdf pdfQCDp("qcdp","qcdp",pfmet,qcdMetp,1); 
  RooDataHist qcdMetm("qcdMetm","qcdMetm",RooArgSet(pfmet),hQCDMetm); RooHistPdf pdfQCDm("qcdm","qcdm",pfmet,qcdMetm,1); 
  
  std::cout << "blalalala3" << std::endl;
  
  vector<RooDataHist*> WlnuMetp_(nIsoBins), WlnuMetm_(nIsoBins);
  vector<RooDataHist*> ewkMetp_(nIsoBins), ewkMetm_(nIsoBins);
  vector<RooDataHist*> qcdMetp_(nIsoBins), qcdMetm_(nIsoBins);
  vector<RooDataHist*> ttbMetp_(nIsoBins), dibMetp_(nIsoBins), wxMetp_(nIsoBins), zxxMetp_(nIsoBins);
  vector<RooDataHist*> ttbMetm_(nIsoBins), dibMetm_(nIsoBins), wxMetm_(nIsoBins), zxxMetm_(nIsoBins);
  
  vector<RooDataHist*> WlnuMetpEta_(nIsoBins), WlnuMetpKeys_(nIsoBins), WlnuMetpStat_(nIsoBins);
  vector<RooDataHist*> WlnuMetmEta_(nIsoBins), WlnuMetmKeys_(nIsoBins), WlnuMetmStat_(nIsoBins);
  vector<RooDataHist*> wxMetpEta_(nIsoBins) , wxMetpKeys_(nIsoBins) , wxMetpStat_(nIsoBins);
  vector<RooDataHist*> wxMetmEta_(nIsoBins) , wxMetmKeys_(nIsoBins) , wxMetmStat_(nIsoBins);
  vector<RooDataHist*> zxxMetpEta_(nIsoBins), zxxMetpKeys_(nIsoBins), zxxMetpStat_(nIsoBins);
  vector<RooDataHist*> zxxMetmEta_(nIsoBins), zxxMetmKeys_(nIsoBins), zxxMetmStat_(nIsoBins);
  vector<RooDataHist*> ewkMetpEta_(nIsoBins), ewkMetpKeys_(nIsoBins), ewkMetpStat_(nIsoBins);
  vector<RooDataHist*> ewkMetmEta_(nIsoBins), ewkMetmKeys_(nIsoBins), ewkMetmStat_(nIsoBins);
  
  vector<RooDataHist*> WlnuMetpEtaD_(nIsoBins), WlnuMetpKeysD_(nIsoBins), WlnuMetpStatD_(nIsoBins);
  vector<RooDataHist*> WlnuMetmEtaD_(nIsoBins), WlnuMetmKeysD_(nIsoBins), WlnuMetmStatD_(nIsoBins);
  vector<RooDataHist*> wxMetpEtaD_(nIsoBins)  , wxMetpKeysD_(nIsoBins)  , wxMetpStatD_(nIsoBins);
  vector<RooDataHist*> wxMetmEtaD_(nIsoBins)  , wxMetmKeysD_(nIsoBins)  , wxMetmStatD_(nIsoBins);
  vector<RooDataHist*> zxxMetpEtaD_(nIsoBins) , zxxMetpKeysD_(nIsoBins) , zxxMetpStatD_(nIsoBins);
  vector<RooDataHist*> zxxMetmEtaD_(nIsoBins) , zxxMetmKeysD_(nIsoBins) , zxxMetmStatD_(nIsoBins);
  vector<RooDataHist*> ewkMetpEtaD_(nIsoBins) , ewkMetpKeysD_(nIsoBins) , ewkMetpStatD_(nIsoBins);
  vector<RooDataHist*> ewkMetmEtaD_(nIsoBins) , ewkMetmKeysD_(nIsoBins) , ewkMetmStatD_(nIsoBins);
  
  
  vector<RooHistPdf*> pdfWep_(nIsoBins), pdfEWKp_(nIsoBins);
  vector<RooHistPdf*> pdfQCDp_(nIsoBins), pdfQCDm_(nIsoBins);
  vector<RooHistPdf*> pdfWem_(nIsoBins), pdfEWKm_(nIsoBins);
  vector<RooHistPdf*> pdfTtbp_(nIsoBins), pdfDibp_(nIsoBins), pdfWxp_(nIsoBins), pdfZxxp_(nIsoBins);
  vector<RooHistPdf*> pdfTtbm_(nIsoBins), pdfDibm_(nIsoBins), pdfWxm_(nIsoBins), pdfZxxm_(nIsoBins);
  
  vector<RooHistPdf*> pdfWepEta_(nIsoBins) , pdfWepKeys_(nIsoBins) , pdfWepStat_(nIsoBins);
  vector<RooHistPdf*> pdfWemEta_(nIsoBins) , pdfWemKeys_(nIsoBins) , pdfWemStat_(nIsoBins);
  vector<RooHistPdf*> pdfEWKpEta_(nIsoBins), pdfEWKpKeys_(nIsoBins), pdfEWKpStat_(nIsoBins);
  vector<RooHistPdf*> pdfEWKmEta_(nIsoBins), pdfEWKmKeys_(nIsoBins), pdfEWKmStat_(nIsoBins);
  vector<RooHistPdf*> pdfWxpEta_(nIsoBins) , pdfWxpKeys_(nIsoBins) , pdfWxpStat_(nIsoBins);
  vector<RooHistPdf*> pdfWxmEta_(nIsoBins) , pdfWxmKeys_(nIsoBins) , pdfWxmStat_(nIsoBins);
  vector<RooHistPdf*> pdfZxxpEta_(nIsoBins), pdfZxxpKeys_(nIsoBins), pdfZxxpStat_(nIsoBins);
  vector<RooHistPdf*> pdfZxxmEta_(nIsoBins), pdfZxxmKeys_(nIsoBins), pdfZxxmStat_(nIsoBins);
  
  vector<RooHistPdf*> pdfWepEtaD_(nIsoBins) , pdfWepKeysD_(nIsoBins) , pdfWepStatD_(nIsoBins);
  vector<RooHistPdf*> pdfWemEtaD_(nIsoBins) , pdfWemKeysD_(nIsoBins) , pdfWemStatD_(nIsoBins);
  vector<RooHistPdf*> pdfEWKpEtaD_(nIsoBins), pdfEWKpKeysD_(nIsoBins), pdfEWKpStatD_(nIsoBins);
  vector<RooHistPdf*> pdfEWKmEtaD_(nIsoBins), pdfEWKmKeysD_(nIsoBins), pdfEWKmStatD_(nIsoBins);
  vector<RooHistPdf*> pdfWxpEtaD_(nIsoBins) , pdfWxpKeysD_(nIsoBins) , pdfWxpStatD_(nIsoBins);
  vector<RooHistPdf*> pdfWxmEtaD_(nIsoBins) , pdfWxmKeysD_(nIsoBins) , pdfWxmStatD_(nIsoBins);
  vector<RooHistPdf*> pdfZxxpEtaD_(nIsoBins), pdfZxxpKeysD_(nIsoBins), pdfZxxpStatD_(nIsoBins);
  vector<RooHistPdf*> pdfZxxmEtaD_(nIsoBins), pdfZxxmKeysD_(nIsoBins), pdfZxxmStatD_(nIsoBins);
  
  // std::cout << "blah here" << std::endl;
  
  TH1D *hIsoBinQCDp = (TH1D*) hDataMetp2d[2]->Clone("hQCDMetpBin");
  TH1D *hIsoBinQCDm = (TH1D*) hDataMetm2d[2]->Clone("hQCDMetmBin");
  hIsoBinQCDp->Add(hWlnuMetp2d[2],-1); hIsoBinQCDp->Add(hEWKMetp2d[2],-1);
  hIsoBinQCDm->Add(hWlnuMetm2d[2],-1); hIsoBinQCDm->Add(hEWKMetm2d[2],-1);
  
  c2->Clear();
  TH1D **hRatios2  = new TH1D*[nIsoBins];
  TH1D **hRatiosNext  = new TH1D*[nIsoBins];
  TH1D *hBin1Norm = (TH1D*) hIsoBinQCDp->Clone("denom");
  hBin1Norm->Scale(1/hBin1Norm->Integral());
  
  
  // doing the direct QCD shape extrapolation from the data histograms
  for(int it=1; it < nIsoBins; ++it){
    it==1?hRatiosNext[it] = (TH1D*) hIsoBinQCDp->Clone("ratio"):hRatiosNext[it] = (TH1D*) hDataMetp2d[it]->Clone("ratio");
    hRatiosNext[it]->Scale(1/hRatiosNext[it]->Integral());

    it==1?hRatios2[it] = (TH1D*) hIsoBinQCDp->Clone("ratio"):hRatios2[it] = (TH1D*) hDataMetp2d[it]->Clone("ratio");
    hRatios2[it]->Scale(1/hRatios2[it]->Integral());

    hRatios2[it]->Divide(hBin1Norm);
    hRatios2[it]->SetMarkerColor(it);
    hRatios2[it]->GetYaxis()->SetRangeUser(-1,3);
    it==1?hRatios2[it]->Draw():hRatios2[it]->Draw("same");
    // c2->Update();
  }
  
  TLine *line = new TLine(0,1,METMAX,1);
  line->SetLineColor(kBlack);
  line->Draw("same");
  sprintf(plotname,"%s/bkgRatios_vs1stBkg.png",CPlot::sOutDir.Data());
  c2->SaveAs(plotname); c2->Clear();

  
  hBin1Norm->Draw();
  for(int it = 1; it < nIsoBins; ++it){
    hRatiosNext[it]->SetMarkerColor(it);
    hRatiosNext[it]->Draw("same");
    
  } 

  c2->SaveAs(plotname); c2->Clear();
  // hRatiosNext[1]->Divide(hRatiosNext[i]);
  for(int it=1; it < nIsoBins-1; ++it){
    hRatiosNext[it]->Divide(hRatiosNext[it+1]);
    hRatiosNext[it]->SetMarkerColor(it);
    hRatiosNext[it]->GetYaxis()->SetRangeUser(0.5,1.5);
    it==1?hRatiosNext[it]->Draw():hRatiosNext[it]->Draw("same");
    c2->Update();
  }
  line->Draw("same");
  sprintf(plotname,"%s/bkgRatios_next.png",CPlot::sOutDir.Data());
  c2->SaveAs(plotname); c2->Clear();
  
  TH1D *hIsoBinHold = (TH1D*)hIsoBinQCDp->Clone("blah");
  hIsoBinHold->Draw();
  
  // hIsoBinQCDp->Multiply(hRatiosNext[1]);
  hIsoBinQCDp->SetMarkerColor(3);
  hIsoBinQCDp->Draw("same");
  sprintf(plotname,"%s/compareWithCorrection.png",CPlot::sOutDir.Data());
  c2->SaveAs(plotname); c2->Clear();
  
  // jesus i should have made 2-d vectors and a giant loop
  for(int j=0; j < nIsoBins; ++j){
    // W+ signal shape "down" uncertainty shapes
    TH1D *hh_diff;
    for(int k=0; k < nMET; ++k){
      if(k==ru||k==rd) continue;
      
      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWlnupMETU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnupMETD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnupMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWlnumMETU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnumMETD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnumMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWEwkMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hEWKpMETU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKpMETD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKpMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hEWKmMETU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKmMETD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKmMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWxpMETU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxpMETD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxpMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hWxmMETU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxmMETD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxmMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hZxxpMETU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxpMETD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxpMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vMET[k].c_str());
      hh_diff =  (TH1D*)hZxxmMETU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxmMETD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxmMETD[j][k]->Add(hh_diff,-1); //delete hh_diff;
    }
    
    for(int k=0; k < nWeight; ++k){
      if(k==pfireu||k==pfired)continue;
        // sprintf(type,"%i_prefireDown",i);

      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWlnupWeightU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnupWeightD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnupWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWlnumWeightU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnumWeightD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnumWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hEWKpWeightU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKpWeightD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hEWKmWeightU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKmWeightD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWxpWeightU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxpWeightD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hWxmWeightU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxmWeightD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hZxxpWeightU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxpWeightD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hZxxmWeightU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxmWeightD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hDibMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hDibpWeightU[j][k]->Clone("diff"); hh_diff->Add(hDibMetp2d[j],-1);
      hDibpWeightD[j][k] = (TH1D*) hDibMetp2d[j]->Clone(nname); hDibpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hDibMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hDibmWeightU[j][k]->Clone("diff"); hh_diff->Add(hDibMetm2d[j],-1);
      hDibmWeightD[j][k] = (TH1D*) hDibMetm2d[j]->Clone(nname); hDibmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hTtbMetpBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hTtbpWeightU[j][k]->Clone("diff"); hh_diff->Add(hTtbMetp2d[j],-1);
      hTtbpWeightD[j][k] = (TH1D*) hTtbMetp2d[j]->Clone(nname); hTtbpWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hTtbMetmBin%d_%sDown",j,vWeight[k].c_str());
      hh_diff =  (TH1D*)hTtbmWeightU[j][k]->Clone("diff"); hh_diff->Add(hTtbMetm2d[j],-1);
      hTtbmWeightD[j][k] = (TH1D*) hTtbMetm2d[j]->Clone(nname); hTtbmWeightD[j][k]->Add(hh_diff,-1); //delete hh_diff;
    }
    
        
    for(int k=0; k < nLHE; ++k){ // the final up/down shapes for the QCD and PDF uncertainty
      // w signal
      
      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hWlnupThyUncU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnupThyUncD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnupThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hWlnumThyUncU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnumThyUncD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnumThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      // ewk total
            
      sprintf(nname,"hEwkMetpBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hEWKpThyUncU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKpThyUncD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKpThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hEWKmThyUncU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKmThyUncD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKmThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      // hEWKpThyUncU[i][j] = new TH1D(("hEwkMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      // hEWKmThyUncU[i][j] = new TH1D(("hEwkMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // wx
            
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hWxpThyUncU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxpThyUncD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxpThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hWxmThyUncU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxmThyUncD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxmThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      // hWxpThyUncU[i][j] = new TH1D(("hWxMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      // hWxmThyUncU[i][j] = new TH1D(("hWxMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // zxx
            
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hZxxpThyUncU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxpThyUncD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxpThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vLHE[k].c_str());
      hh_diff =  (TH1D*)hZxxmThyUncU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxmThyUncD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxmThyUncD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      // hZxxpThyUncU[i][j] = new TH1D(("hZxxMetpBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      // hZxxmThyUncU[i][j] = new TH1D(("hZxxMetmBin"+std::to_string(i)+"_"+vLHE[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
    }
  }
  
   std::cout << "making pdfs from histograms" << std::endl;
  // char hname[50];
  // char pname[50];
  
  std::cout << "integrals" << std::endl;
  std::cout << "signal " << hWlnuMetp2d[0]->Integral() << std::endl;
  std::cout << "ewk " << hEWKMetp2d[0]->Integral() << std::endl;
  std::cout << "ttb " << hTtbMetp2d[0]->Integral() << std::endl;
  std::cout << "dib " << hDibMetp2d[0]->Integral() << std::endl;
  std::cout << "wx " << hWxMetp2d[0]->Integral() << std::endl;
  std::cout << "zxx " << hZxxMetp2d[0]->Integral() << std::endl;
  
  for (int j = 0; j < nIsoBins; ++j){
      
      std::cout << "start datahistPDFs for Bin " << j << std::endl;
      // function to do the 2-step pdf construction (also write it)
      makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2d[j],WlnuMetp_,pdfWep_,pfmet,j,"");
      makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2d[j],ewkMetp_,pdfEWKp_,pfmet,j,"");
      makeDataHistPdf("ttbMETp","ttbp",hTtbMetp2d[j],ttbMetp_,pdfTtbp_,pfmet,j,"");
      makeDataHistPdf("dibMETp","dibp",hDibMetp2d[j],dibMetp_,pdfDibp_,pfmet,j,"");
      makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2d[j] ,wxMetp_ ,pdfWxp_ ,pfmet,j,"");
      makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2d[j],zxxMetp_,pdfZxxp_,pfmet,j,""); 
      
      j==0?makeDataHistPdf("qcdMETp","qcdp",hIsoBinQCDp,qcdMetp_,pdfQCDp_,pfmet,j,""):makeDataHistPdf("qcdMETp","qcdp",hDataMetp2d[j],qcdMetp_,pdfQCDp_,pfmet,j,"");
	  // ----------------------------------------- W- ---------------------------------
      makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2d[j],WlnuMetm_,pdfWem_,pfmet,j,"");
      makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2d[j],ewkMetm_,pdfEWKm_,pfmet,j,"");
      makeDataHistPdf("ttbMETm","ttbm",hTtbMetm2d[j],ttbMetm_,pdfTtbm_,pfmet,j,"");
      makeDataHistPdf("dibMETm","dibm",hDibMetm2d[j],dibMetm_,pdfDibm_,pfmet,j,"");
      makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2d[j] ,wxMetm_ ,pdfWxm_ ,pfmet,j,"");
      makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2d[j],zxxMetm_,pdfZxxm_,pfmet,j,"");
      
      std::cout << "do QCD" << std::endl;
      j==0?makeDataHistPdf("qcdMetm","qcdm",hIsoBinQCDm,qcdMetm_,pdfQCDm_,pfmet,j,""):makeDataHistPdf("qcdMetm","qcdm",hDataMetm2d[j],qcdMetm_,pdfQCDm_,pfmet,j,"");
      
      std::cout << "at end of looop" << std::endl;
  }
  
  std::cout << "finished making the datahistpdf" << std::endl;
  
  // QCD Pdfs
  
  CPepeModel2 qcd("qcd",pfmet);
  CPepeModel2 qcdp("qcdp",pfmet);
  CPepeModel2 qcdm("qcdm",pfmet);
 
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
  
  
    
  
  // // Anti-Signal PDFs
  // RooDataHist aWlnuMet ("aWlnuMET", "aWlnuMET", RooArgSet(pfmet),hAntiWlnuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,aWlnuMet, 1);
  // RooDataHist aWlnuMetp("aWlnuMETp","aWlnuMETp",RooArgSet(pfmet),hAntiWlnuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,aWlnuMetp,1);
  // RooDataHist aWlnuMetm("aWlnuMETm","aWlnuMETm",RooArgSet(pfmet),hAntiWlnuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,aWlnuMetm,1); 
  
  // // Anti-EWK+top PDFs
  // RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", pfmet,aewkMet, 1);
  // RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",pfmet,aewkMetp,1); 
  // RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",pfmet,aewkMetm,1); 
  
// //   // Anti-QCD Pdfs
  // CPepeModel2 aqcd("aqcd",pfmet, qcd.a1);
  // CPepeModel2 aqcdp("aqcdp",pfmet, qcdp.a1);
  // CPepeModel2 aqcdm("aqcdm",pfmet, qcdm.a1);

//   CPepeModel2 aqcd("aqcd",pfmet);
//   CPepeModel2 aqcdp("aqcdp",pfmet);
//   CPepeModel2 aqcdm("aqcdm",pfmet);
  
  // // Anti-selection PDFs
  // RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  // RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  // RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
  
  // regular linear dependence
  vector <CPepeModel2isobinsMuons*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  qcdp_[0] = new CPepeModel2isobinsMuons("qcdp2d0",pfmet, 0.0);
  qcdm_[0] = new CPepeModel2isobinsMuons("qcdm2d0",pfmet, 0.0);  
  
  // Make some goddamn Pepe
  for (int j = 1; j < nIsoBins; ++j){
	  // // Original isolation linear dependence
      sprintf(nname, "qcdp2d%d",j);
      qcdp_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  sprintf(nname, "qcdm2d%d",j);
      qcdm_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
  
  }

   vector < RooAddPdf*> pdfMetp_(nIsoBins), pdfMetm_(nIsoBins);
  for(int j = 0; j < nIsoBins; ++j){
	 if(j==0||j==1||j==2||j==3||j==4){

     sprintf(nname,"pdfWep%d",j);

		 (doMET&&(!doTemplate)) ? pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j])) : pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfWxp_[j],*pdfZxxp_[j],*pdfDibp_[j],*pdfTtbp_[j],*pdfQCDp_[j]),RooArgList(*nSigp_[j],*nWxp_[j],*nZxxp_[j],*nDibp_[j],*nTtbp_[j],*nQCDp_[j]));
		 
		 sprintf(nname,"pdfWem%d",j);
		 (doMET&&(!doTemplate)) ? pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j])) : pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfWxm_[j],*pdfZxxm_[j],*pdfDibm_[j],*pdfTtbm_[j],*pdfQCDm_[j]),RooArgList(*nSigm_[j],*nWxm_[j],*nZxxm_[j],*nDibm_[j],*nTtbm_[j],*nQCDm_[j]));		
    
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
  
  // std::cout << "made pdfs" << std::endl;
  
  
  // // PDF for simultaneous fit  
  // RooCategory rooCat("rooCat","rooCat");
  // rooCat.defineType("Select");
  // rooCat.defineType("Anti");
  
  // RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  // pdfTotal.addPdf(pdfMet, "Select");
  // //pdfTotal.addPdf(apdfMet,"Anti");
  
  // RooSimultaneous pdfTotalp("pdfTotalp","pdfTotalp",rooCat);
  // pdfTotalp.addPdf(pdfMetp, "Select");
  // pdfTotalp.addPdf(apdfMetp,"Anti");
  
  // RooSimultaneous pdfTotalm("pdfTotalm","pdfTotalm",rooCat);
  // pdfTotalm.addPdf(pdfMetm, "Select");
  // pdfTotalm.addPdf(apdfMetm,"Anti");

  //
  // Perform fits
  //
  // RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  // RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  // RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);
 
  // RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
  // RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
  // RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);

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
	  dataMetp_[nname]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetp2d[i]);
      sprintf(nname, "isom%d",i);
	  dataMetm_[nname]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetm2d[i]);
      sprintf(nname, "isop%d",i);
      dataMetpHist_[i]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetp2d[i]);
	  sprintf(nname, "isom%d",i);
      dataMetmHist_[i]=new RooDataHist(nname,nname, RooArgSet(pfmet), hDataMetm2d[i]);
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  // set up RooCategory for data
  // test only 2 bins at first
  
  RooCategory rooCat2dTest("rooCat2dTest","rooCat2dTest");
  for (int i = 0; i < nIsoBins; ++ i){
	  sprintf(nname,"iso%d",i); rooCat2dTest.defineType(nname);
  }
  // RooDataHist sigDatap("dataTotalp","dataTotalp", RooArgList(pfmet),Import(dataMetp_));
  
  
  // RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp2d[0]);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm2d[0]);

  RooDataHist combDatap("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetp_));
  RooDataHist combDatam("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetm_));
  
  RooSimultaneous simPdfp("simPdfp","simultaneous pdf W+",rooCat2dTest) ;
  RooSimultaneous simPdfm("simPdfm","simultaneous pdf W-",rooCat2dTest) ;
  for(int i = 0; i < nIsoBins; ++i ){
	  sprintf(nname, "isop%d", i); simPdfp.addPdf(*pdfMetp_[i],nname);
	  sprintf(nname, "isom%d", i); simPdfm.addPdf(*pdfMetm_[i],nname);
  }

  
  // cout << "Starting values for Wlnu yields: " << endl;
  // cout << "Selected: " << hDataMet->Integral() << endl;
  // cout << "   sig: " << hWlnuMet->Integral() << endl;
  // cout << "   EWK: " << hEWKMet->Integral() << endl;
  // cout << "   qcd: " << hDataMet->Integral()-hWlnuMet->Integral()-hEWKMet->Integral() << endl;

  // cout << "Starting values for Wlnu_p yields: " << endl;
  // cout << "   sig: " << hWlnuMetp->Integral() << endl;
  // cout << "   EWK: " << hEWKMetp->Integral() << endl;
  // cout << "   qcd: " << hDataMetp->Integral()-hWlnuMetp->Integral()-hEWKMetp->Integral() << endl;

  // cout << "Starting values for Wlnu_m yields: " << endl;
  // cout << "   sig: " << hWlnuMetm->Integral() << endl;
  // cout << "   EWK: " << hEWKMetm->Integral() << endl;
  // cout << "   qcd: " << hDataMetm->Integral()-hWlnuMetm->Integral()-hEWKMetm->Integral() << endl;
  
  
  // cout << "Starting values for AntiWlnu yields: " << endl;
  // cout << "Selected: " << hAntiDataMet->Integral() << endl;
  // cout << "   sig: " << hAntiWlnuMet->Integral() << endl;
  // cout << "   EWK: " << hAntiEWKMet->Integral() << endl;
  // cout << "   qcd: " << hAntiDataMet->Integral()-hAntiWlnuMet->Integral()-hAntiEWKMet->Integral() << endl;

  // cout << "Starting values for AntiWlnu_p yields: " << endl;
  // cout << "   sig: " << hWlnuMetp->Integral() << endl;
  // cout << "   EWK: " << hEWKMetp->Integral() << endl;
  // cout << "   qcd: " << hAntiDataMetp->Integral()-hAntiWlnuMetp->Integral()-hAntiEWKMetp->Integral() << endl;

  // cout << "Starting values for AntiWlnu_m yields: " << endl;
  // cout << "   sig: " << hAntiWlnuMetm->Integral() << endl;
  // cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;
  // cout << "   qcd: " << hAntiDataMetm->Integral()-hAntiWlnuMetm->Integral()-hAntiEWKMetm->Integral() << endl;

// //   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());
// //   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",0.3*(hDataMet->Integral()),0,hDataMet->Integral());

  RooRealVar pepe2Pdf_qcdp_norm("pepe2Pdf_qcdp_norm","pepe2Pdf_qcdp_norm",0.3*(hDataMetp->Integral()),0,hDataMetp->Integral());
  RooRealVar pepe2Pdf_qcdm_norm("pepe2Pdf_qcdm_norm","pepe2Pdf_qcdm_norm",0.3*(hDataMetm->Integral()),0,hDataMetm->Integral());
  
  RooRealVar pepe2Pdf_aqcdp_norm("pepe2Pdf_aqcdp_norm","pepe2Pdf_aqcdp_norm",0.9*(hAntiDataMetp->Integral()),0,hAntiDataMetp->Integral());
  RooRealVar pepe2Pdf_aqcdm_norm("pepe2Pdf_aqcdm_norm","pepe2Pdf_aqcdm_norm",0.9*(hAntiDataMetm->Integral()),0,hAntiDataMetm->Integral());
  
  
  
  std::vector<RooRealVar*> pepe2Pdf_qcdp_norms_(nIsoBins), pepe2Pdf_qcdm_norms_(nIsoBins);
  for(int j = 0; j < nIsoBins; ++j){
     // sprintf(nname,"pepe2Pdf_qcdp_norm");
     if(j==0){
       sprintf(nname,"pepe2Pdf_qcdp2d%i_norm",j);
       pepe2Pdf_qcdp_norms_[j] = new RooRealVar(nname,nname,0.4*hDataMetp2d[j]->Integral(),0,0.6*hDataMetp2d[j]->Integral());
       sprintf(nname,"pepe2Pdf_qcdm2d%i_norm",j);
       pepe2Pdf_qcdm_norms_[j] = new RooRealVar(nname,nname,0.4*hDataMetm2d[j]->Integral(),0,0.6*hDataMetm2d[j]->Integral());
     } else {
       sprintf(nname,"pepe2Pdf_qcdp2d%i_norm",j);
       pepe2Pdf_qcdp_norms_[j] = new RooRealVar(nname,nname,hDataMetp2d[j]->Integral(),0.95*hDataMetp2d[j]->Integral(),1.05*hDataMetp2d[j]->Integral());
       sprintf(nname,"pepe2Pdf_qcdm2d%i_norm",j);
       pepe2Pdf_qcdm_norms_[j] = new RooRealVar(nname,nname,hDataMetm2d[j]->Integral(),0.95*hDataMetm2d[j]->Integral(),1.05*hDataMetm2d[j]->Integral());
     }
  }
  std::cout << "mark 3 " << std::endl;

  TH1 *qcdp2 = (TH1*) hQCDMetp2d[2]->Clone();   TH1 *qcdp3=(TH1*) hQCDMetp2d[3]->Clone();
  TH1 *qcdm2 = (TH1*) hQCDMetm2d[2]->Clone();   TH1 *qcdm3=(TH1*) hQCDMetm2d[3]->Clone();
  qcdp3->Scale(qcdp2->Integral()/qcdp3->Integral());
  qcdm3->Scale(qcdm2->Integral()/qcdm3->Integral());
  TH1D *hQCDMetpBin2_shapeUp   = (TH1D*) hQCDMetp2d[2]->Clone();  
  hQCDMetpBin2_shapeUp->SetTitle("hQCDMetpBin2_shapeUp");
  hQCDMetpBin2_shapeUp->SetName("hQCDMetpBin2_shapeUp");
  TH1D *hQCDMetpBin2_shapeDown = (TH1D*) hQCDMetp2d[2]->Clone();
  hQCDMetpBin2_shapeDown->SetTitle("hQCDMetpBin2_shapeDown");
  hQCDMetpBin2_shapeDown->SetName("hQCDMetpBin2_shapeDown");
  TH1D *hQCDMetmBin2_shapeUp   = (TH1D*) hQCDMetm2d[2]->Clone();
  hQCDMetmBin2_shapeUp->SetTitle("hQCDMetmBin2_shapeUp");
  hQCDMetmBin2_shapeUp->SetName("hQCDMetmBin2_shapeUp");
  TH1D *hQCDMetmBin2_shapeDown  = (TH1D*) hQCDMetm2d[2]->Clone();
  hQCDMetmBin2_shapeDown->SetTitle("hQCDMetmBin2_shapeDown");
  hQCDMetmBin2_shapeDown->SetName("hQCDMetmBin2_shapeDown");

  for(Int_t ibin=1; ibin<=qcdp2->GetNbinsX(); ibin++) 
    {
      double targetp = qcdp3->GetBinContent(ibin);
      double medianp = qcdp2->GetBinContent(ibin);
      double targetm = qcdm3->GetBinContent(ibin);
      double medianm = qcdm2->GetBinContent(ibin);
      hQCDMetpBin2_shapeUp->SetBinContent(ibin,targetp);
      hQCDMetmBin2_shapeUp->SetBinContent(ibin,targetm);
      hQCDMetpBin2_shapeDown->SetBinContent(ibin,2*medianp-targetp);
      hQCDMetmBin2_shapeDown->SetBinContent(ibin,2*medianm-targetm);
    }
  //create qcd up/down shapes

//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",40000,0,100000);
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",40000,0,100000);
  TString histfname = outputDir + TString("/Wlnu_Histograms.root");
  TFile *histFile = new TFile(histfname,"RECREATE");
  histFile->cd();
  hQCDMetpBin2_shapeUp->Write();
  hQCDMetmBin2_shapeUp->Write();
  hQCDMetpBin2_shapeDown->Write();
  hQCDMetmBin2_shapeDown->Write();
  for(int j = 0; j < nIsoBins; ++j){
    // std::cout << "writing" << std::endl;
    hDataMetp2d[j]->Write();
    hDataMetm2d[j]->Write();
    
    hWlnuMetp2d[j]->Write();
    hWlnuMetm2d[j]->Write();
    
    hEWKMetp2d[j]->Write();
    hEWKMetm2d[j]->Write();
    
    hZxxMetp2d[j]->Write();
    hZxxMetm2d[j]->Write();
    
    hWxMetp2d[j]->Write();
    hWxMetm2d[j]->Write();
    
    hTtbMetp2d[j]->Write();
    hTtbMetm2d[j]->Write();
    
    hDibMetp2d[j]->Write();
    hDibMetm2d[j]->Write();

    hQCDMetp2d[j]->Write();
    hQCDMetm2d[j]->Write();
    
    hIsoBinQCDp->Write();
    hIsoBinQCDm->Write();
    
    for(int k=0; k < nMET; k++){
    // std::cout << "writing MET" << std::endl;
      hWlnupMETU[j][k]->Write();
      hWlnumMETU[j][k]->Write();
      
      hEWKpMETU[j][k]->Write();
      hEWKmMETU[j][k]->Write();
      
      hWxpMETU[j][k]->Write();
      hWxmMETU[j][k]->Write();
      
      hZxxpMETU[j][k]->Write();
      hZxxmMETU[j][k]->Write();
      
      // std::cout << "writing MET DOWN" << std::endl;
      if(k==rd||k==ru)continue;
      hWlnupMETD[j][k]->Write();
      hWlnumMETD[j][k]->Write();
            
      hWxpMETD[j][k]->Write();
      hWxmMETD[j][k]->Write();
            
      hEWKpMETD[j][k]->Write();
      hEWKmMETD[j][k]->Write();
      
      hZxxpMETD[j][k]->Write();
      hZxxmMETD[j][k]->Write();
      
    }
    
    for(int k=0; k < nWeight; k++){
      // std::cout << "writing weight" << std::endl;
      hWlnupWeightU[j][k]->Write();
      hWlnumWeightU[j][k]->Write();

      
      hEWKpWeightU[j][k]->Write();
      hEWKmWeightU[j][k]->Write();
      

      
      hWxpWeightU[j][k]->Write();
      hWxmWeightU[j][k]->Write();

      
      hZxxpWeightU[j][k]->Write();
      hZxxmWeightU[j][k]->Write();

      hDibpWeightU[j][k]->Write();
      hDibmWeightU[j][k]->Write();

      hTtbpWeightU[j][k]->Write();
      hTtbmWeightU[j][k]->Write();
      // std::cout << "writing weight DOWN" << std::endl;
      if(k==pfireu||k==pfired)continue;
            
      hWlnupWeightD[j][k]->Write();
      hWlnumWeightD[j][k]->Write();
      
      hEWKpWeightD[j][k]->Write();
      hEWKmWeightD[j][k]->Write();
      
      hWxpWeightD[j][k]->Write();
      hWxmWeightD[j][k]->Write();
      
      hZxxpWeightD[j][k]->Write();
      hZxxmWeightD[j][k]->Write();
      
      hDibpWeightD[j][k]->Write();
      hDibmWeightD[j][k]->Write();
      
      hTtbpWeightD[j][k]->Write();
      hTtbmWeightD[j][k]->Write();
    }
    // for(int k=0; k < nLHE; k++){
      // // std::cout << "writing LHE" << std::endl;
      // hWlnupThyUncU[j][k]->Write();
      // hWlnumThyUncU[j][k]->Write();
      
      // hWlnupThyUncD[j][k]->Write();
      // hWlnumThyUncD[j][k]->Write();
      
      // hEWKpThyUncU[j][k]->Write();
      // hEWKmThyUncU[j][k]->Write();
      
      // hEWKpThyUncD[j][k]->Write();
      // hEWKmThyUncD[j][k]->Write();
      
      // hWxpThyUncU[j][k]->Write();
      // hWxmThyUncU[j][k]->Write();
      
      // hWxpThyUncD[j][k]->Write();
      // hWxmThyUncD[j][k]->Write();
      
      // hZxxpThyUncU[j][k]->Write();
      // hZxxmThyUncU[j][k]->Write();
      
      // hZxxpThyUncD[j][k]->Write();
      // hZxxmThyUncD[j][k]->Write();
    // }
    
  }
  histFile->Write();
  histFile->Close();
  
  

  
  
    RooGaussian constm("constm","constm",nEWKm,RooConst(hEWKMetm->Integral()),RooConst(0.15*hEWKMetm->Integral()));
    RooGaussian constp("constp","constp",nEWKp,RooConst(hEWKMetp->Integral()),RooConst(0.15*hEWKMetp->Integral()));
    
    RooGaussian const_wxm("const_wxm","const_wxm",*nWxm_[0],RooConst(hWxMetm2d[0]->Integral()),RooConst(0.15*hWxMetm2d[0]->Integral()));
    RooGaussian const_wxp("const_wxp","const_wxp",*nWxp_[0],RooConst(hWxMetp2d[0]->Integral()),RooConst(0.15*hWxMetp2d[0]->Integral()));
    
    RooGaussian const_zxxm("const_zxxm","const_zxxm",*nZxxm_[0],RooConst(hZxxMetm2d[0]->Integral()),RooConst(0.15*hZxxMetm2d[0]->Integral()));
    RooGaussian const_zxxp("const_zxxp","const_zxxp",*nZxxp_[0],RooConst(hZxxMetp2d[0]->Integral()),RooConst(0.15*hZxxMetp2d[0]->Integral()));
    
    RooGaussian const_dibm("const_dibm","const_dibm",*nDibm_[0],RooConst(hDibMetm2d[0]->Integral()),RooConst(0.15*hDibMetm2d[0]->Integral()));
    RooGaussian const_dibp("const_dibp","const_dibp",*nDibp_[0],RooConst(hDibMetp2d[0]->Integral()),RooConst(0.15*hDibMetp2d[0]->Integral()));
    
    RooGaussian const_ttbm("const_ttbm","const_ttbm",*nTtbm_[0],RooConst(hTtbMetm2d[0]->Integral()),RooConst(0.15*hTtbMetm2d[0]->Integral()));
    RooGaussian const_ttbp("const_ttbp","const_ttbp",*nTtbp_[0],RooConst(hTtbMetp2d[0]->Integral()),RooConst(0.15*hTtbMetp2d[0]->Integral()));
    
    RooGaussian constantim("constantim","constantim",nAntiSigm,RooConst(hAntiWlnuMetm->Integral()),RooConst(0.15*hAntiWlnuMetm->Integral()));
    RooGaussian constantip("constantip","constantip",nAntiSigp,RooConst(hAntiWlnuMetp->Integral()),RooConst(0.15*hAntiWlnuMetp->Integral()));
	


  // RooFitResult *fitResp2dCatTest = pdfMetp.fitTo(dataMetp,Extended(),Save(kTRUE),RooFit::Strategy(2),Minos(kTRUE),Minimizer("Minuit2","minimize"),PrintEvalErrors(-1));
  // RooFitResult *fitResm2dCatTest = pdfMetm.fitTo(dataMetm,Extended(),Save(kTRUE),RooFit::Strategy(2),Minos(kTRUE),Minimizer("Minuit2","minimize"),PrintEvalErrors(-1));

 
  
  RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),/*ExternalConstraints(RooArgList(const_wxp,const_zxxp,const_dibp,const_ttbp)),*/RooFit::Strategy(2),Minos(kTRUE),Minimizer("Minuit2","minimize"),PrintEvalErrors(-1));
  RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),/*ExternalConstraints(RooArgList(const_wxm,const_zxxm,const_dibm,const_ttbm)),*/RooFit::Strategy(2),Minos(kTRUE),Minimizer("Minuit2","minimize"),PrintEvalErrors(-1));
  
  // TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  // hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  // TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  // hMetDiff->SetMarkerStyle(kFullCircle); hMetDiff->SetMarkerSize(0.9);
   
  // TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  // for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWlnuMetp->GetBinError(ibin));}
  // // std::cout << nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal() << std::endl;
  // // std::cout << nSigp.getVal()+nWxp.getVal()+nZxxp.getVal()+nDibp.getVal()+nTtbp.getVal()+nQCDp.getVal() << std::endl;
  // hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  // TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  // hMetpDiff->SetMarkerStyle(kFullCircle); hMetpDiff->SetMarkerSize(0.9);
    
  // TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet)); // why did we not just clone the original histogram...
  // for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWlnuMetm->GetBinError(ibin));}
  // hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  // TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  // hMetmDiff->SetMarkerStyle(kFullCircle); hMetmDiff->SetMarkerSize(0.9);
  
  // // the diff hists for the anti-selection
  
  // TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  
  // hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  // TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  // hAntiMetDiff->SetMarkerStyle(kFullCircle);
  // hAntiMetDiff->SetMarkerSize(0.9);
   
  // TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
   // for(int ibin = 1; ibin < hPdfAntiMetp->GetNbinsX(); ++ibin){hPdfAntiMetp->SetBinError(ibin, hAntiWlnuMetp->GetBinError(ibin));}
  // hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  // TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  // hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  // hAntiMetpDiff->SetMarkerSize(0.9);
    
  // TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
   // for(int ibin = 1; ibin < hPdfAntiMetm->GetNbinsX(); ++ibin){hPdfAntiMetm->SetBinError(ibin, hAntiWlnuMetm->GetBinError(ibin));}
  // hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  // TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  // hAntiMetmDiff->SetMarkerStyle(kFullCircle);
  // hAntiMetmDiff->SetMarkerSize(0.9);
    
    // char plotname[100];
    //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
    // label for lumi
  // char lumitext[100];
  // if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  (8 TeV)",lumi*1000.);
  // else         sprintf(lumitext,"%.1f fb^{-1}  (13 TeV)",lumi/1000.);
  

  
    //
  // // Dummy histograms for TLegend
  // // (I can't figure out how to properly pass RooFit objects...)
  // //
  // TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  // hDummyData->SetMarkerStyle(kFullCircle);
  // hDummyData->SetMarkerSize(0.9);
  
  // TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  // hDummyW->SetLineColor(linecolorW);
  // hDummyW->SetFillColor(fillcolorW);
  // hDummyW->SetFillStyle(1001);
  
  // TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  // hDummyEWK->SetLineColor(linecolorEWK);
  // hDummyEWK->SetFillColor(fillcolorEWK);
  // hDummyEWK->SetFillStyle(1001);
  
  // TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  // hDummyQCD->SetLineColor(linecolorQCD);
  // hDummyQCD->SetFillColor(fillcolorQCD);
  // hDummyQCD->SetFillStyle(1001);
  
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
  
  
  
  // label for lumi
  char lumitext[100];
  
  if(lumi > 200) sprintf(lumitext,"%.1f pb^{-1}  (5 TeV)",lumi);
  else         sprintf(lumitext,"%.1f pb^{-1}  (13 TeV)",lumi);

  
    
  Double_t chi2probm, chi2ndfm, chi2probp, chi2ndfp;
  Double_t ksprobm, ksprobpem, ksprobp, ksprobpep;
  
    // ofstream txtfile4;
    // char txtfname4[100];
    // std::cout << "Printing We+. " << std::endl;
    // sprintf(txtfname4,"%s/chi2.txt",CPlot::sOutDir.Data());
    // txtfile4.open(txtfname4);
    // assert(txtfile4.is_open());
  
  for(int i = 0; i < nIsoBins; ++i){
      
    // set up the diff plot
    // turn this into its own function later?
    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetp = (TH1D*)(pdfMetp_[i]->createHistogram("hPdfMetp", pfmet));
    for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWlnuMetp2d[i]->GetBinError(ibin));}
    // std::cout << nSigp_[i]->getVal()+nEWKp_[i]->getVal()+nQCDp_[i]->getVal() << std::endl;
    std::cout << nSigp_[i]->getVal()+nWxp_[i]->getVal()+nZxxp_[i]->getVal()+nDibp_[i]->getVal()+nTtbp_[i]->getVal()+nQCDp_[i]->getVal() << std::endl;
    // hPdfMetp->Scale((nSigp_[i]->getVal()+nEWKp_[i]->getVal()+nQCDp_[i]->getVal())/hPdfMetp->Integral());
    hPdfMetp->Scale((nSigp_[i]->getVal()+nWxp_[i]->getVal()+nZxxp_[i]->getVal()+nDibp_[i]->getVal()+nTtbp_[i]->getVal()+nQCDp_[i]->getVal())/hPdfMetp->Integral());
    TH1D *hMetpDiff = makeDiffHist(hDataMetp2d[i],hPdfMetp,"hMetpDiff");
    hMetpDiff->SetMarkerStyle(kFullCircle); hMetpDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetm = (TH1D*)(pdfMetm_[i]->createHistogram("hPdfMetm", pfmet));
    for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWlnuMetm2d[i]->GetBinError(ibin));}
    std::cout << nSigm_[i]->getVal()+nWxm_[i]->getVal()+nZxxm_[i]->getVal()+nDibm_[i]->getVal()+nTtbm_[i]->getVal()+nQCDm_[i]->getVal() << std::endl;
    hPdfMetm->Scale((nSigm_[i]->getVal()+nWxm_[i]->getVal()+nZxxm_[i]->getVal()+nDibm_[i]->getVal()+nTtbm_[i]->getVal()+nQCDm_[i]->getVal())/hPdfMetm->Integral());
    TH1D *hMetmDiff = makeDiffHist(hDataMetm2d[i],hPdfMetm,"hMetmDiff");
    hMetmDiff->SetMarkerStyle(kFullCircle); hMetmDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    sprintf(nname,"isop%d",i); sprintf(plotname, "wep_fitmetp_bin%i",i);
    drawWMetPlotsSplit(plotname, hMetpDiff, pfmet, dataMetp_[nname], pdfMetp_[i], pdfWxp_[i],pdfZxxp_[i],pdfDibp_[i],pdfTtbp_[i], (RooAbsPdf*)pdfQCDp_[i], pdfWep_[i], lumitext, hDataMetp2d[i]);
    
    sprintf(nname,"isom%d",i); sprintf(plotname, "wem_fitmetm_bin%i",i);
    drawWMetPlotsSplit(plotname, hMetmDiff, pfmet, dataMetm_[nname], pdfMetm_[i], pdfWxm_[i],pdfZxxm_[i],pdfDibm_[i],pdfTtbm_[i], (RooAbsPdf*)pdfQCDm_[i], pdfWem_[i], lumitext, hDataMetm2d[i]);


  }
  
  // txtfile4.close();
  
  // std::cout << "hi" << std::endl;
  ofstream txtfile;
  char txtfname1[100];
  std::cout << "Printing We+. " << std::endl;
  sprintf(txtfname1,"%s/fitresWe2dCatTest.txt",CPlot::sOutDir.Data());
  txtfile.open(txtfname1);
  assert(txtfile.is_open());
  
  // flags = txtfile.flags();
  // txtfile << setprecision(10);
  // txtfile << " *** Yields *** " << endl;
  txtfile << "Data+: " << hDataMetp2d[0]->Integral() << endl;
  txtfile << "Data-: " << hDataMetm2d[0]->Integral() << endl;
  
  
  fitResp2dCatTest->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose); 
  txtfile << endl;
  // printCorrelations(txtfile, fitResp2dCatTest);
  txtfile << endl;
  // printChi2AndKSResults(txtfile, chi2probp, chi2ndfp, ksprobp, ksprobpep);
  txtfile<< endl;
  
  fitResm2dCatTest->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  // printCorrelations(txtfile, fitResm2dCatTest);
  txtfile << endl;
  // // printChi2AndKSResults(txtfile, chi2probm, chi2ndfm, ksprobm, ksprobpem);
  txtfile <<endl;
  
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("fitWm");
}

// TH1D *QCD


//=== FUNCTION DEFINITIONS ======================================================================================


void fillMETs(bool doMET,TH1D** h,vector<double> met, int nMET, double wgt, double mtCorr){
  // return;
  for(int k =2; k < nMET; k++){
    // doMET ? h[k] ->Fill(met[k] , wgt) : h[k] ->Fill(mtCorr,wgt);
    h[k] ->Fill(met[k] , wgt);
  }
  // std::cout << "h->int" << h[3]->Integral() << std::endl;
  return;
}

void fillWeights(bool doMET,TH1D** h,double met, int nWeight,vector<double> wgt,double mtCorr){
  // return;
  // std::cout <<  wgt[0]<<" " <<  wgt[9]<<" " <<  wgt[10]<<" " << std::endl;
  for(int k =1; k < nWeight; k++){
    // doMET ? h[k] ->Fill(met , wgt[k]) : h[k] ->Fill(mtCorr,wgt[k]);
    h[k] ->Fill(met , wgt[k]);
  }
  return;
}


void fillLHE(TH1D** hlhe, double met, double evtweight, vector<double> *lheweight){
  // return;
  // std::cout << "filling lhe weights" << std::endl;
  for(int k = 0; k < nQCD+nPDF; k ++){
    // std::cout << k << "  " <<  (*lheweight)[k] << std::endl;
    // std::cout << hlhe[k]->Integral() << std::endl;
    hlhe[k]->Fill(met, evtweight*((*lheweight)[k]));
  } 
  
  // std::cout << "done with lhe weights" << std::endl;
}

void calcLHE(TH1D* hQCD, TH1D* hPDF, TH1D** hlhe, TH1D* hMain,  bool isSignal){
  // return;
  
  // rescale the histograms if it's the signal shape
  // std::cout << "1" << std::endl;
    double sigNorm=hMain->Integral();
    // cout << "other norm " << sigNorm << endl;
  if(isSignal){
    for(int i=0; i < nQCD+nPDF; i++){
      double norm=hlhe[i]->Integral();
      // cout << "old norm " << norm << endl;
      if(norm!=0)hlhe[i]->Scale(sigNorm/norm);
      // cout << "new norm " << hlhe[i]->Integral() << endl;
    }
  }

  // std::cout << "1" << std::endl;
  // star loop through the bins
  for(int ibin=1; ibin<=hMain->GetNbinsX(); ibin++){
    // cout << "WORKING WITH BIN " << ibin << " with content " << hMain->GetBinContent(ibin) << endl;
    double qcdval=0;
    // loop through the 6 for QCD
    // std::cout <<"nqcd" << nQCD << std::endl;
    for(int k=0; k < nQCD; k++){
     // std::cout << k<< std::endl;
     // std::cout << hlhe[k]->Integral() << std::endl;
     // cout << "qcd bin contents orig: " << hMain->GetBinContent(ibin) << " qcd " <<  hlhe[k]->GetBinContent(ibin) << endl;
       fabs(hlhe[k]->GetBinContent(ibin)-hMain->GetBinContent(ibin)) > qcdval ? qcdval=fabs(hlhe[k]->GetBinContent(ibin)-hMain->GetBinContent(ibin)) : 0;
    }
    // std::cout << "blah" << std::endl;
    // std::cout << hMain->Integral() << std::endl;
    // std::cout << hQCD->Integral() << std::endl;
    // cout << "qcd val " << qcdval << endl;
    // cout << "with bin content " << hMain->GetBinContent(ibin) << "  , sum : " << hMain->GetBinContent(ibin)+qcdval << endl;
    hQCD->SetBinContent(ibin,hMain->GetBinContent(ibin)+qcdval);
    
    // std::cout << "done qcd" << std::endl;
    
    // cout << "Do the PDF calculation on bin with : " << hMain->GetBinContent(ibin) <<  endl;
    double pdfval=0;
    for(int k=nQCD; k < nQCD+nPDF; k++ ){
      // cout << "pdf bin contents orig: " << hMain->GetBinContent(ibin) << " pdf " <<  hlhe[k]->GetBinContent(ibin) << endl;
      double diff = (hlhe[k]->GetBinContent(ibin)-hMain->GetBinContent(ibin))/hMain->GetBinContent(ibin);
      pdfval+=diff*diff;
      
    // std::cout << diff << std::endl;
    }
    pdfval=sqrt(pdfval/nPDF);
    
    // cout << "pdf? " << endl;
    // std::cout << pdfval << endl;
    // cout << "make sure its empty " << hPDF->GetBinContent(ibin) << endl;;
    hPDF->SetBinContent(ibin,hMain->GetBinContent(ibin)*(1+pdfval));
    
    // cout << "orig content " << hMain->GetBinContent(ibin) <<  " new "  << hMain->GetBinContent(ibin)*(1+pdfval) << endl; 
    
  // std::cout << "2" << std::endl;
  }
  
  // std::cout << "3" << std::endl;
}

void drawLHE(TH1D** hlhe, TH1D* hMain, TString name, TString outdir, bool isSignal){
  // return;
  // std::cout << "draw lhe " << std::endl;
  // rescale the histograms if it's the signal shape
  if(isSignal){
    double sigNorm=hMain->Integral();
    for(int i=0; i < nQCD+nPDF; i++){
      double norm=hlhe[i]->Integral();
      if(norm!=0)hlhe[i]->Scale(sigNorm/norm);
    }
  }
  // make a canvas to draw onto
  TCanvas *c = new TCanvas("c","c",800,600);
  char outfile[150];
  // draw main shape 
  hMain->SetLineColor(kBlack);
  hMain->SetMarkerSize(0);
  hMain->SetLineWidth(5);
  hMain->Draw("L");
  // loop through the 6 for QCD
  for(int k=0; k < nQCD; k++){
    std::cout << hlhe[k]->Integral() << std::endl;
    hlhe[k]->SetMarkerColor(kAzure);
    hlhe[k]->SetLineColor(kAzure);
    hlhe[k]->SetMarkerSize(0);
    hlhe[k]->SetLineWidth(1);
    hlhe[k]->Draw("same");
  }
  // save and clear canvas
  sprintf(outfile,"%s/%s_QCDshapes.png",outdir.Data(),name.Data());
  c->SaveAs(outfile);
  c->Clear();
  // draw main shape again
  hMain->SetLineColor(kBlack);
  hMain->SetMarkerSize(0);
  hMain->SetLineWidth(5);
  hMain->Draw("L");
  // loop through the 100 for PDF
  for(int k=nQCD; k < nQCD+nPDF; k++ ){
    hlhe[k]->SetMarkerColor(kAzure);
    hlhe[k]->SetLineColor(kAzure);
    hlhe[k]->SetMarkerSize(0);
    hlhe[k]->SetLineWidth(1);
    hlhe[k]->Draw("same");
  }
  sprintf(outfile,"%s/%s_PDFshapes.png",outdir.Data(),name.Data());
  c->SaveAs(outfile);
  // std::cout << "blah" << std::endl;
  delete c;
}

void drawShapes(TH1D** vars, TH1D* hMain, TString outdir, TString name, vector<string> leg, int max){
   // make a canvas to draw onto
  TCanvas *c = new TCanvas("c","c",800,600);
  TLegend* legend = new TLegend(0.6, 0.5, .99, .99);
  char outfile[150];
  // draw main shape 
  hMain->SetLineColor(kBlack);
  hMain->SetMarkerSize(0);
  hMain->SetLineWidth(5);
  hMain->Draw("L");
  // loop through the 6 for QCD
  // std::cout << "here" << std::endl;
  for(int k=1; k < max; k++){
    // std::cout << vars[k]->Integral() << std::endl;
    vars[k]->SetMarkerColor(k+1);
    vars[k]->SetLineColor(k+1);
    vars[k]->SetMarkerSize(0);
    vars[k]->SetLineWidth(1);
    vars[k]->Draw("same");
    legend->AddEntry(vars[k],(leg[k]).c_str(),"l");
    // std::cout << "ljljljlk" << std::endl;
  }
  // save and clear canvas
  legend->Draw();
  sprintf(outfile,"%s/%s_all.png",outdir.Data(),name.Data());
  c->SaveAs(outfile);
  c->Clear();
  
  // std::cout << "blah" << std::endl;
  delete c;
  
}

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    if(hData->GetBinContent(ibin) == 0) diff = 0;
    
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    if(hData->GetBinContent(ibin) == 0) err = 0;
    std::cout << "ibin " << ibin << "  data " << hData->GetBinContent(ibin) << "  fit " << hFit->GetBinContent(ibin) << "  err " << err << std::endl;
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

void makeDataHistPdf(string dh, string hp, TH1D* hIn, vector<RooDataHist*> &vDataHist, vector<RooHistPdf*> &vHistPdf, RooRealVar &x, int it, string sfx){
  char nname[100];
  sprintf(nname, "%s%d%s",dh.c_str(),it,sfx.c_str()); vDataHist[it] = new RooDataHist(nname, nname, RooArgSet(x),hIn);
  sprintf(nname, "%s%d%s",hp.c_str(),it,sfx.c_str()); vHistPdf[it] = new RooHistPdf(nname,nname,x,*vDataHist[it],1);   
  return;
}

void drawWMetPlots(string plotname, TH1D *diff, RooRealVar &x, RooDataHist* dat, RooAddPdf* pdf, RooHistPdf* ewk, RooAbsPdf* qcd, RooHistPdf* wsigp, string lumitext, TH1D* hData){
  double yscale=0.05;
  
  const TString format("png"); 
  char ylabel[100];  // string buffer for y-axis label
  // plot colors
  // Int_t linecolorW   = kOrange-3;
  // Int_t fillcolorW   = kOrange-2;
  // Int_t linecolorEWK = kOrange+10;
  // Int_t fillcolorEWK = kOrange+7;
  // Int_t linecolorQCD = kViolet+2;
  // Int_t fillcolorQCD = kViolet-5;
  // Int_t ratioColor   = kGray+2;

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

  RooPlot *frame = x.frame(Bins(NBINS)); 
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetXaxis()->SetLabelOffset(2.0);
  
  // sprintf(nname,"isop%d",i);
  // sprintf(plotname, "wep_fitmetp_bin%i",i);
  dat->plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdf->plotOn(frame,FillColor(fillcolorW),DrawOption("F"));
  pdf->plotOn(frame,LineColor(linecolorW));
  
  pdf->plotOn(frame,Components(RooArgSet(*ewk,*qcd)),FillColor(fillcolorEWK),DrawOption("F"));
  pdf->plotOn(frame,Components(RooArgSet(*ewk,*qcd)),LineColor(linecolorEWK));
  pdf->plotOn(frame,Components(RooArgSet(*qcd)),LineColor(linecolorQCD));
  pdf->plotOn(frame,Components(RooArgSet(*qcd)),FillColor(fillcolorQCD),DrawOption("F"));

  pdf->plotOn(frame,Components(RooArgSet(*wsigp)),LineColor(linecolorW),LineStyle(2));
  dat->plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
 
  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plot(plotname,frame,"","",ylabel);
  plot.SetLegend(0.68,0.57,0.93,0.77);
  plot.GetLegend()->AddEntry(hDummyData,"data","PL");
  plot.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#bar{#nu}","F");
  plot.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plot.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plot.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plot.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plot.Draw(c,kFALSE,format,1);
  
  std::cout << "Draw the W plot diff" << std::endl;
  CPlot plotDiff(plotname,"","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  diff->GetYaxis()->SetTitleOffset(0.5);
  diff->GetYaxis()->SetLabelSize(0.11);
  plotDiff.SetYRange(-yscale,yscale);
  plotDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotDiff.AddHist1D(diff,"EX0",ratioColor);
  plotDiff.Draw(c,kTRUE,format,2);
  plotDiff.Draw(c,kTRUE,"pdf",2);

  std::cout << "Draw the W plot log" << std::endl;
  plot.SetName((plotname+"_log").c_str());
  plot.SetLogy();
  plot.SetYRange(1e-5*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plot.Draw(c,kTRUE,format,1);
  plot.Draw(c,kTRUE,"pdf",1);
}

void drawWMetPlotsSplit(string plotname, TH1D *diff, RooRealVar &x, RooDataHist* dat, RooAddPdf* pdf, RooHistPdf* wx,RooHistPdf* zxx,RooHistPdf* dib,RooHistPdf* ttb, RooAbsPdf* qcd, RooHistPdf* wsigp, string lumitext, TH1D* hData){
  double yscale=0.05;
  
  const TString format("png"); 
  char ylabel[100];  // string buffer for y-axis label
  // plot colors
  // Int_t linecolorW   = kOrange-3;
  // Int_t fillcolorW   = kOrange-2;
  // Int_t linecolorEWK = kOrange+10;
  // Int_t fillcolorEWK = kOrange+7;
  // Int_t linecolorQCD = kViolet+2;
  // Int_t fillcolorQCD = kViolet-5;
  // Int_t ratioColor   = kGray+2;

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

  RooPlot *frame = x.frame(Bins(NBINS)); 
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetXaxis()->SetLabelOffset(2.0);
  
  // sprintf(nname,"isop%d",i);
  // sprintf(plotname, "wep_fitmetp_bin%i",i);
  dat->plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdf->plotOn(frame,FillColor(fillcolorW),DrawOption("F"));
  pdf->plotOn(frame,LineColor(linecolorW));
  
  pdf->plotOn(frame,Components(RooArgSet(*wx,*qcd)),FillColor(fillcolorEWK),DrawOption("F"));
  pdf->plotOn(frame,Components(RooArgSet(*wx,*qcd)),LineColor(linecolorEWK));
  pdf->plotOn(frame,Components(RooArgSet(*zxx,*wx,*qcd)),FillColor(fillcolorEWK),DrawOption("F"));
  pdf->plotOn(frame,Components(RooArgSet(*zxx,*wx,*qcd)),LineColor(linecolorEWK));
  pdf->plotOn(frame,Components(RooArgSet(*ttb,*zxx,*wx,*qcd)),FillColor(fillcolorEWK),DrawOption("F"));
  pdf->plotOn(frame,Components(RooArgSet(*ttb,*zxx,*wx,*qcd)),LineColor(linecolorEWK));
  pdf->plotOn(frame,Components(RooArgSet(*dib,*ttb,*zxx,*wx,*qcd)),FillColor(fillcolorEWK),DrawOption("F"));
  pdf->plotOn(frame,Components(RooArgSet(*dib,*ttb,*zxx,*wx,*qcd)),LineColor(linecolorEWK));
  pdf->plotOn(frame,Components(RooArgSet(*qcd)),LineColor(linecolorQCD));
  pdf->plotOn(frame,Components(RooArgSet(*qcd)),FillColor(fillcolorQCD),DrawOption("F"));

  pdf->plotOn(frame,Components(RooArgSet(*wsigp)),LineColor(linecolorW),LineStyle(2));
  dat->plotOn(frame,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
 
  sprintf(ylabel,"Events / %.1f GeV",hData->GetBinWidth(1));
  CPlot plot(plotname,frame,"","",ylabel);
  plot.SetLegend(0.68,0.57,0.93,0.77);
  plot.GetLegend()->AddEntry(hDummyData,"data","PL");
  plot.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#bar{#nu}","F");
  plot.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plot.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plot.AddTextBox("#bf{CMS}",0.62,0.80,0.88,0.88,0);
  plot.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plot.Draw(c,kFALSE,format,1);
  
  std::cout << "Draw the W plot diff" << std::endl;
  CPlot plotDiff(plotname,"","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  diff->GetYaxis()->SetTitleOffset(0.5);
  diff->GetYaxis()->SetLabelSize(0.11);
  plotDiff.SetYRange(-yscale,yscale);
  plotDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotDiff.AddHist1D(diff,"EX0",ratioColor);
  plotDiff.Draw(c,kTRUE,format,2);
  plotDiff.Draw(c,kTRUE,"pdf",2);

  std::cout << "Draw the W plot log" << std::endl;
  plot.SetName((plotname+"_log").c_str());
  plot.SetLogy();
  plot.SetYRange(1e-5*(hData->GetMaximum()),10*(hData->GetMaximum()));
  plot.Draw(c,kTRUE,format,1);
  plot.Draw(c,kTRUE,"pdf",1);
}

void makeUncMT(vector<Double_t> &metVars, vector<Double_t> &metVarsPhi, TLorentzVector* lep){
  for(int i = 0; i < metVars.size(); ++i){
    double mtCorr  = sqrt(2.0*(lep->Pt()) * (metVars[i]) * (1.0-cos(toolbox::deltaPhi(lep->Phi(),metVarsPhi[i]))));
    metVars[i] = mtCorr;
  }
  return;
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
  sprintf(htmlfname,"%s/WlnuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wlnu</title></head>" << endl;
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
