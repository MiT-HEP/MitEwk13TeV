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
// #include "../Utils/RecoilCorrector_asym2.hh"
// #include "../Utils/RecoilCorrector_addJets.hh"
// #include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"

// #include "ZBackgrounds.hh"

// // // helper class to handle rochester corrections
// // // #include <rochcor2015r.h>
// // // #include <muresolution_run2r.h>
// #include <../RochesterCorr/RoccoR.cc>

// // // helper class to handle efficiency tables
// #include "../Utils/CEffUser1D.hh"
// #include "../Utils/CEffUser2D.hh"

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

void makeDataHistPdf(string dh, string hp, TH1D* hIn, vector<RooDataHist*> &vDataHist, vector<RooHistPdf*> &vHistPdf, RooRealVar &x, int it,string sfx);

void drawWMetPlots(string plotname, TH1D *diff, RooRealVar &x, RooDataHist* dat, RooAddPdf* pdf, RooHistPdf* ewk, RooAbsPdf* qcd, RooHistPdf* wsigp, string lumitext, TH1D* hData);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
               const Double_t ksprob, const Double_t ksprobpe);

// make webpage
void makeHTML(const TString outDir);

// global variables yolo
const  Int_t linecolorW   = kOrange-3;
const  Int_t fillcolorW   = kOrange-2;
const  Int_t linecolorEWK = kOrange+10;
const  Int_t fillcolorEWK = kOrange+7;
const  Int_t linecolorQCD = kViolet+2;
const  Int_t fillcolorQCD = kViolet-5;
const  Int_t ratioColor   = kGray+2;

const Int_t    NBINS   = 75;
const Double_t METMIN  = 0;
const Double_t METMAX  = 150;

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
  bool doMET = true;
  bool doTemplate = true;

  double yscale=0.5;
  
  // vector<char[50]> uncShapeNames = {};
  enum{eta,keys,stat};
  enum{mc,fsr,bkg};
  const string vUncRec[]={"eta","keys","stat"};
  const string vUncEff[]={"mc","fsr","bkg"};
  int nUncRec = sizeof(vUncRec)/sizeof(vUncRec[0]);
  int nUncEff = sizeof(vUncEff)/sizeof(vUncEff[0]);
  
  // Double_t vIsoBins[] = {0.0,0.20,0.30,0.40,0.50,0.60,0.70};
  Double_t vIsoBins[] = {0.0,0.15,0.25,0.35,0.45,0.55,0.65};
  int nIsoBins = sizeof(vIsoBins)/sizeof(vIsoBins[0])-1;
  std::cout << "size of isobin array is " << nIsoBins << std::endl;
  
  // MET histogram binning and range
  // const Int_t    NBINS   = 75*4;
  // const Int_t    NBINS   = 125;
  // const Int_t    NBINS   = 50;
  // const Int_t    NBINS   = 100;
  
  const Double_t PT_CUT  = 25;
  // const Double_t PT_CUT  = 35;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t MT_CUT = 25.0;
  
  const Double_t mu_MASS = 0.1057;

  //
  // Set up output file
  //
  TString outfilename = outputDir + TString("/") + TString("Zmm_DataBkg.root");
  TFile *outFile = new TFile(outfilename,"RECREATE");
  TH1::AddDirectory(kFALSE);

  
  // not sure i'll be needing these?
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
  
  // TFile *_rat2 = new TFile("shapeDiff/UnfoldingOutputZPt.root");
  // TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  // hh_diff = (TH1D*)_rat2->Get("hUnfold");
  // hh_diff->Scale(1/hh_diff->Integral()); // normalize
  // hh_diff->Divide(hh_mc);
   
  //
  // input ntuple file names
  //
  enum {eData, eWlnu, eZxx, eWx, eTtb, eDib, eQCD, eAntiData, eAntiWlnu, eAntiQCD, eAntiTtb, eAntiDib, eAntiWx, eAntiZxx };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
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
  
  
  // // actually now we're switching to 2-d arrays haha, actually leave to test later
  // TH1D ***hWlnuMetp2d  = new TH1D*[5];// hAntiWlnuMetp->Sumw2();
  // TH1D ***hWlnuMetm2d  = new TH1D*[5];// hAntiWlnuMetm->Sumw2();
  // // TH1D **hAntiEWKMetIsoBins    = new TH1D*[5];// hAntiEWKMet->Sumw2();
  // TH1D ***hEWKMetp2d   = new TH1D*[5];// hAntiEWKMetp->Sumw2();
  // TH1D ***hEWKMetm2d   = new TH1D*[5];// hAntiEWKMetm->Sumw2();
  
  // TH1D ***hWxMetp2d   = new TH1D*[5];// hAntiEWKMetp->Sumw2();
  // TH1D ***hWxMetm2d   = new TH1D*[5];// hAntiEWKMetm->Sumw2();
  
  // TH1D ***hZxxMetp2d   = new TH1D*[5];// hAntiEWKMetp->Sumw2();
  // TH1D ***hZxxMetm2d   = new TH1D*[5];// hAntiEWKMetm->Sumw2();
  
  
  // TH1D **hAntiDataMetIsoBins   = new TH1D*[5];// hAntiDataMet->Sumw2();
  TH1D **hDataMetm2d  = new TH1D*[nIsoBins];// hAntiDataMetm->Sumw2();  
  TH1D **hDataMetp2d  = new TH1D*[nIsoBins];// hAntiDataMetp->Sumw2();
  // TH1D **hAntiWlnuMetIsoBins   = new TH1D*[5];// hAntiWlnuMet->Sumw2();
  // these are the uncorrected shapes
  TH1D **hWlnuMetp2d  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D **hWlnuMetm2d  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // All teh uncertainties? 
  // this is an array of arrays of histograms // first index is iso bin, 2nd is uncertainty shape
  TH1D ***hWlnuMetp2dUncRecU  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnuMetm2dUncRecU  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  TH1D ***hWlnuMetp2dUncRecD  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnuMetm2dUncRecD  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  
  TH1D ***hEWKMetp2dUncRecU  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKMetm2dUncRecU  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  TH1D ***hEWKMetp2dUncRecD  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKMetm2dUncRecD  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D ***hWxMetp2dUncRecU  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxMetm2dUncRecU  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  TH1D ***hWxMetp2dUncRecD  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxMetm2dUncRecD  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  
  TH1D ***hZxxMetp2dUncRecU  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxMetm2dUncRecU  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  TH1D ***hZxxMetp2dUncRecD  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxMetm2dUncRecD  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  ////////////////////////////////////////////////////////////////////////////////////////
  TH1D ***hWlnuMetp2dUncEffU  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnuMetm2dUncEffU  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  TH1D ***hWlnuMetp2dUncEffD  = new TH1D**[nIsoBins];// hAntiWlnuMetp->Sumw2();
  TH1D ***hWlnuMetm2dUncEffD  = new TH1D**[nIsoBins];// hAntiWlnuMetm->Sumw2();
  
  TH1D ***hEWKMetp2dUncEffU  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKMetm2dUncEffU  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  TH1D ***hEWKMetp2dUncEffD  = new TH1D**[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D ***hEWKMetm2dUncEffD  = new TH1D**[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D ***hWxMetp2dUncEffU  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxMetm2dUncEffU  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  TH1D ***hWxMetp2dUncEffD  = new TH1D**[nIsoBins];// hAntiWxMetp->Sumw2();
  TH1D ***hWxMetm2dUncEffD  = new TH1D**[nIsoBins];// hAntiWxMetm->Sumw2();
  
  TH1D ***hZxxMetp2dUncEffU  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxMetm2dUncEffU  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  TH1D ***hZxxMetp2dUncEffD  = new TH1D**[nIsoBins];// hAntiZxxMetp->Sumw2();
  TH1D ***hZxxMetm2dUncEffD  = new TH1D**[nIsoBins];// hAntiZxxMetm->Sumw2();
  
  TH1D ***hDibMetp2dUncEffU  = new TH1D**[nIsoBins];// hAntiDibMetp->Sumw2();
  TH1D ***hDibMetm2dUncEffU  = new TH1D**[nIsoBins];// hAntiDibMetm->Sumw2();
  TH1D ***hDibMetp2dUncEffD  = new TH1D**[nIsoBins];// hAntiDibMetp->Sumw2();
  TH1D ***hDibMetm2dUncEffD  = new TH1D**[nIsoBins];// hAntiDibMetm->Sumw2();
  
  TH1D ***hTtbMetp2dUncEffU  = new TH1D**[nIsoBins];// hAntiTtbMetp->Sumw2();
  TH1D ***hTtbMetm2dUncEffU  = new TH1D**[nIsoBins];// hAntiTtbMetm->Sumw2();
  TH1D ***hTtbMetp2dUncEffD  = new TH1D**[nIsoBins];// hAntiTtbMetp->Sumw2();
  TH1D ***hTtbMetm2dUncEffD  = new TH1D**[nIsoBins];// hAntiTtbMetm->Sumw2();
  
  // do a loop to create the second array
  for(int i=0; i < nIsoBins; ++i){
    // w signal recoil
    hWlnuMetp2dUncRecU[i] = new TH1D*[nUncRec];
    hWlnuMetm2dUncRecU[i] = new TH1D*[nUncRec];
    
    hWlnuMetp2dUncRecD[i] = new TH1D*[nUncRec];
    hWlnuMetm2dUncRecD[i] = new TH1D*[nUncRec];
    
    // w signal efficiency 
    hWlnuMetp2dUncEffU[i] = new TH1D*[nUncEff];
    hWlnuMetm2dUncEffU[i] = new TH1D*[nUncEff];
    
    hWlnuMetp2dUncEffD[i] = new TH1D*[nUncEff];
    hWlnuMetm2dUncEffD[i] = new TH1D*[nUncEff];
    
    // ewk total recoil
    hEWKMetp2dUncRecU[i] = new TH1D*[nUncRec];
    hEWKMetm2dUncRecU[i] = new TH1D*[nUncRec];
    
    hEWKMetp2dUncRecD[i] = new TH1D*[nUncRec];
    hEWKMetm2dUncRecD[i] = new TH1D*[nUncRec];
    
    // ewk total efficiency
    hEWKMetp2dUncEffU[i] = new TH1D*[nUncEff];
    hEWKMetm2dUncEffU[i] = new TH1D*[nUncEff];
    
    hEWKMetp2dUncEffD[i] = new TH1D*[nUncEff];
    hEWKMetm2dUncEffD[i] = new TH1D*[nUncEff];
    
    // wx  recoil
    hWxMetp2dUncRecU[i] = new TH1D*[nUncRec];
    hWxMetm2dUncRecU[i] = new TH1D*[nUncRec];
    
    hWxMetp2dUncRecD[i] = new TH1D*[nUncRec];
    hWxMetm2dUncRecD[i] = new TH1D*[nUncRec];
    
    // wx  efficiency
    hWxMetp2dUncEffU[i] = new TH1D*[nUncEff];
    hWxMetm2dUncEffU[i] = new TH1D*[nUncEff];
    
    hWxMetp2dUncEffD[i] = new TH1D*[nUncEff];
    hWxMetm2dUncEffD[i] = new TH1D*[nUncEff];
    
    // zxx recoil
    hZxxMetp2dUncRecU[i] = new TH1D*[nUncRec];
    hZxxMetm2dUncRecU[i] = new TH1D*[nUncRec];
    
    hZxxMetp2dUncRecD[i] = new TH1D*[nUncRec];
    hZxxMetm2dUncRecD[i] = new TH1D*[nUncRec];
    
    // zxx efficiency
    hZxxMetp2dUncEffU[i] = new TH1D*[nUncEff];
    hZxxMetm2dUncEffU[i] = new TH1D*[nUncEff];
    
    hZxxMetp2dUncEffD[i] = new TH1D*[nUncEff];
    hZxxMetm2dUncEffD[i] = new TH1D*[nUncEff];
    
    // diboson efficiency
    hDibMetp2dUncEffU[i] = new TH1D*[nUncEff];
    hDibMetm2dUncEffU[i] = new TH1D*[nUncEff];
    
    hDibMetp2dUncEffD[i] = new TH1D*[nUncEff];
    hDibMetm2dUncEffD[i] = new TH1D*[nUncEff];
    
    // ttbar efficiency
    hTtbMetp2dUncEffU[i] = new TH1D*[nUncEff];
    hTtbMetm2dUncEffU[i] = new TH1D*[nUncEff];
    
    hTtbMetp2dUncEffD[i] = new TH1D*[nUncEff];
    hTtbMetm2dUncEffD[i] = new TH1D*[nUncEff];
  }
  // Eta-binned unc 
  // TH1D **hWlnuMetp2dEta  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  // TH1D **hWlnuMetm2dEta  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // // Keys unc
  // TH1D **hWlnuMetp2dKeys  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  // TH1D **hWlnuMetm2dKeys  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // // Stat Unc
  // TH1D **hWlnuMetp2dStat  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  // TH1D **hWlnuMetm2dStat  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
    // // Eta-binned unc 
  // TH1D **hWlnuMetp2dEtaD  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  // TH1D **hWlnuMetm2dEtaD  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // // Keys unc
  // TH1D **hWlnuMetp2dKeysD  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  // TH1D **hWlnuMetm2dKeysD  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // // Stat Unc
  // TH1D **hWlnuMetp2dStatD  = new TH1D*[nIsoBins];// hAntiWlnuMetp->Sumw2();
  // TH1D **hWlnuMetm2dStatD  = new TH1D*[nIsoBins];// hAntiWlnuMetm->Sumw2();
  // TH1D **hAntiEWKMetIsoBins    = new TH1D*[5];// hAntiEWKMet->Sumw2();
  TH1D **hEWKMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hEWKMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // eta binned
  // TH1D **hEWKMetp2dEta   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hEWKMetm2dEta   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // keys
  // TH1D **hEWKMetp2dKeys   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hEWKMetm2dKeys   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // stat unc
  // TH1D **hEWKMetp2dStat   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hEWKMetm2dStat   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // eta binned
  // TH1D **hEWKMetp2dEtaD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hEWKMetm2dEtaD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // keys
  // TH1D **hEWKMetp2dKeysD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hEWKMetm2dKeysD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // stat unc
  // TH1D **hEWKMetp2dStatD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hEWKMetm2dStatD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  
  TH1D **hDibMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hDibMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hTtbMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hTtbMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hWxMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hWxMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // eta binned
  // TH1D **hWxMetp2dEta   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hWxMetm2dEta   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // keys unc
  // TH1D **hWxMetp2dKeys   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hWxMetm2dKeys   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // stat unc
  // TH1D **hWxMetp2dStat   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hWxMetm2dStat   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  // // eta binned
  // TH1D **hWxMetp2dEtaD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hWxMetm2dEtaD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // keys unc
  // TH1D **hWxMetp2dKeysD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hWxMetm2dKeysD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // stat unc
  // TH1D **hWxMetp2dStatD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hWxMetm2dStatD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  TH1D **hZxxMetp2d   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  TH1D **hZxxMetm2d   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // eta binned unc
  // TH1D **hZxxMetp2dEta   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hZxxMetm2dEta   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // keys unc
  // TH1D **hZxxMetp2dKeys   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hZxxMetm2dKeys   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // stat unc
  // TH1D **hZxxMetp2dStat   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hZxxMetm2dStat   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
    // // eta binned unc
  // TH1D **hZxxMetp2dEtaD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hZxxMetm2dEtaD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // keys unc
  // TH1D **hZxxMetp2dKeysD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hZxxMetm2dKeysD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  // // stat unc
  // TH1D **hZxxMetp2dStatD   = new TH1D*[nIsoBins];// hAntiEWKMetp->Sumw2();
  // TH1D **hZxxMetm2dStatD   = new TH1D*[nIsoBins];// hAntiEWKMetm->Sumw2();
  
  
  TH1D **hMetpIsoValues = new TH1D*[nIsoBins];
  TH1D **hMetmIsoValues = new TH1D*[nIsoBins];
  // Create a histogram pointer in each space in the array
  for(int i = 0; i < nIsoBins; i++){
    hDataMetm2d[i]  = new TH1D(("hDataMetmBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    hDataMetp2d[i]  = new TH1D(("hDataMetpBin"+std::to_string(i)).c_str(),"",  NBINS,METMIN,METMAX); 
    
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
    
    // Create a loop over the # of uncertainty shapes to produce the uncertainty histograms for the "up" shapes
    // Do the recoil ones here
    for(int j = 0; j < nUncRec; ++j){
      // w signal
      hWlnuMetp2dUncRecU[i][j] = new TH1D(("hWlnuMetpBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hWlnuMetm2dUncRecU[i][j] = new TH1D(("hWlnuMetmBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // ewk total
      hEWKMetp2dUncRecU[i][j] = new TH1D(("hEwkMetpBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hEWKMetm2dUncRecU[i][j] = new TH1D(("hEwkMetmBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // wx
      hWxMetp2dUncRecU[i][j] = new TH1D(("hWxMetpBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hWxMetm2dUncRecU[i][j] = new TH1D(("hWxMetmBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      // zxx
      hZxxMetp2dUncRecU[i][j] = new TH1D(("hZxxMetpBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hZxxMetm2dUncRecU[i][j] = new TH1D(("hZxxMetmBin"+std::to_string(i)+"_"+vUncRec[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
    }
    
    for(int j=0; j < nUncEff; ++j){
      hWlnuMetp2dUncEffU[i][j] = new TH1D(("hWlnuMetpBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hWlnuMetm2dUncEffU[i][j] = new TH1D(("hWlnuMetmBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      hEWKMetp2dUncEffU[i][j] = new TH1D(("hEwkMetpBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hEWKMetm2dUncEffU[i][j] = new TH1D(("hEwkMetmBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      hWxMetp2dUncEffU[i][j] = new TH1D(("hWxMetpBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hWxMetm2dUncEffU[i][j] = new TH1D(("hWxMetmBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
  
      hZxxMetp2dUncEffU[i][j] = new TH1D(("hZxxMetpBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hZxxMetm2dUncEffU[i][j] = new TH1D(("hZxxMetmBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
      hDibMetp2dUncEffU[i][j] = new TH1D(("hDibMetpBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hDibMetm2dUncEffU[i][j] = new TH1D(("hDibMetmBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
  
      hTtbMetp2dUncEffU[i][j] = new TH1D(("hTtbMetpBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      hTtbMetm2dUncEffU[i][j] = new TH1D(("hTtbMetmBin"+std::to_string(i)+"_"+vUncEff[j]+"Up").c_str(),"",NBINS,METMIN,METMAX);
      
    }

	  hMetpIsoValues[i] = new TH1D(("hMetpIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
    hMetmIsoValues[i] = new TH1D(("hMetmIsoValues_"+std::to_string(i)).c_str(),"",  100,vIsoBins[i],vIsoBins[i+1]);
  }
  


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

  TFile *_rat2 = new TFile("TEST_Zmm_13TeV_incl_ZpTrw_v0/zPt_Normal13TeV.root");
  TH1D *hh_diff;// = new TH1D("hh_diff","hh_diff",75,0,150);
  hh_diff = (TH1D*)_rat2->Get("hZptRatio");  

//   
  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum, npv, npu;
  Float_t genVPt, genVPhi, genVy, genLepPt, genLepPhi;
  Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0, *lep_raw=0, *genV=0, *genLep=0;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Double_t metCorrLep, metCorrMain, metCorrEta, metCorrStat, metCorrKeys, mtCorr;
  Double_t metCorrLepPhi, metCorrMainPhi, metCorrEtaPhi, metCorrStatPhi, metCorrKeysPhi;
  Double_t totalEvtWeight=1, effSFweight=1, relIso;
  Double_t evtWeightSysFSR=1, evtWeightSysMC=1, evtWeightSysBkg=1;
  
    
  
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
    intree->SetBranchAddress("prefireWeight",  &prefireWeight);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("totalEvtWeight", &totalEvtWeight);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("evtWeightSysFSR",&evtWeightSysFSR);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("evtWeightSysMC", &evtWeightSysMC);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("evtWeightSysBkg",&evtWeightSysBkg);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb",      &scale1fb);  // MC weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp",    &scale1fbUp);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown",  &scale1fbDown);  // event weight per 1/fb (MC)
    intree->SetBranchAddress("metCorrLep",      &metCorrLep);     // MET including lepton scale/smear
    intree->SetBranchAddress("metCorrLepPhi",   &metCorrLepPhi);  // MET phi including lepton scale/smear
    intree->SetBranchAddress("metCorrMain",     &metCorrMain);     // MET including lepton scale/smear (w/ main recoil corrs)
    intree->SetBranchAddress("metCorrMainPhi",  &metCorrMainPhi);  // MET phi including lepton scale/smear (w/ main recoil corrs)
    intree->SetBranchAddress("metCorrEta",      &metCorrEta);      // MET including lepton scale/smear (w/ eta recoil corrs)
    intree->SetBranchAddress("metCorrEtaPhi",   &metCorrEtaPhi);   // MET phi including lepton scale/smear  (w/ eta recoil corrs)
    intree->SetBranchAddress("metCorrStat",     &metCorrStat);     // MET including lepton scale/smear (w/ stat unc recoil corrs)
    intree->SetBranchAddress("metCorrStatPhi",  &metCorrStatPhi);  // MET phi including lepton scale/smear (w/ stat unc recoil corrs)
    intree->SetBranchAddress("metCorrKeys",     &metCorrKeys);     // MET including lepton scale/smear (w/ keyspdf recoil corrs)
    intree->SetBranchAddress("metCorrKeysPhi",  &metCorrKeysPhi);  // MET phi including lepton scale/smear (w/ keyspdf recoil corrs)
    intree->SetBranchAddress("sumEt",        &sumEt);     // Sum ET
    intree->SetBranchAddress("mt",           &mt);        // transverse mass
    intree->SetBranchAddress("mtCorr",       &mtCorr);        // transverse mass
    intree->SetBranchAddress("q",        &q);         // lepton charge
    intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    intree->SetBranchAddress("lep_raw",      &lep_raw);       // lepton 4-vector
    intree->SetBranchAddress("genLep",      &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV",         &genV);       // lepton 4-vector
    intree->SetBranchAddress("pfChIso",      &pfChIso);
    intree->SetBranchAddress("pfGamIso",     &pfGamIso);
    intree->SetBranchAddress("pfNeuIso",     &pfNeuIso);
    intree->SetBranchAddress("pfCombIso",    &pfCombIso);       // lepton 4-vector
    intree->SetBranchAddress("relIso",       &relIso);       // relative isolation for the lepton
  
    UInt_t iterator=15;
    // UInt_t iterator=1;
    if(typev[ifile]==eData||typev[ifile]==eAntiData)iterator=1;
    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << std::endl;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<(int)(intree->GetEntries()*0.1); ientry++) {
    // for(UInt_t ientry=0; ientry<((int)intree->GetEntries()); ientry+=iterator) {
      intree->GetEntry(ientry);
      if(ientry%100000==0) cout << "Event " << ientry << ". " << (double)ientry/(double)intree->GetEntries()*100 << " % done with this file." << endl;

// figure out later what to do
        // if(typev[ifile]==eWlnu || typev[ifile]==eWx || typev[ifile]==eZxx) {
          // // what was this
            // double bin = 0;
            // for(int i = 1; i <= hh_diff->GetNbinsX();++i){
              // if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
            // }
            // double w2 = 1.0;//hh_diff->GetBinContent(bin);
      // std::cout << " pass PT " << std::endl;
      if(lep_raw->Pt() < PT_CUT) continue;//std::cout << " pass PT " << std::endl;
      if(fabs(lep->Eta()) > ETA_CUT) continue;//std::cout << " pass eta " << std::endl;
      if(doMTCut&&(mtCorr<MT_CUT)) continue;//std::cout << " pass mt " << std::endl;
          
      // set up the event weights for the MC reweighting
      Double_t weight=totalEvtWeight*lumi;
      Double_t wgtFsr=evtWeightSysFSR*lumi;
      Double_t wgtBkg=evtWeightSysBkg*lumi;
      Double_t wgtMC=evtWeightSysMC*lumi;
      
      // std::cout << weight << "  " << wgtFsr << "  " << wgtBkg << "" << wgtMC <<  " pt " << lep->Pt() << "  eta " <<  lep->Eta()<< std::endl;
      // Double_t weight=totalEvtWeight*lumi2;    
      
        // std::cout << typev[ifile]<< std::endl;
      if(typev[ifile]==eData) {
        hDataMet->Fill(metCorrLep);
        if(q>0) {
          doMET ? hDataMetp->Fill(metCorrLep) : hDataMetp->Fill(mtCorr);
          hDataMetpPhi->Fill(metCorrLepPhi);
          hMuonEtaDatap->Fill(fabs(lep->Eta()));
          doMET ? hDataMetp2d[0]->Fill(metCorrLep) : hDataMetp2d[0]->Fill(mtCorr);
          hMetpIsoValues[0]->Fill(relIso);
        } else {
          doMET ? hDataMetm->Fill(metCorrLep) : hDataMetm->Fill(mtCorr);
          hMuonEtaDatam->Fill(fabs(lep->Eta()));
          hDataMetmPhi->Fill(metCorrLepPhi);
          doMET ? hDataMetm2d[0]->Fill(metCorrLep) : hDataMetm2d[0]->Fill(mtCorr);
          hMetmIsoValues[0]->Fill(relIso);
        }
      } else if(typev[ifile]==eAntiData) {
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            hAntiDataMet->Fill(metCorrLep);
            if(q>0) { 
              doMET ? hAntiDataMetp->Fill(metCorrLep) : hAntiDataMetp->Fill(mtCorr);
              hMuonEtaAntiDatap->Fill(fabs(lep->Eta()));
              doMET ? hDataMetp2d[it]->Fill(metCorrLep) : hDataMetp2d[it]->Fill(mtCorr);
              hMetpIsoValues[it]->Fill(relIso);
            } else { 
              doMET ? hAntiDataMetm->Fill(metCorrLep) : hAntiDataMetm->Fill(mtCorr);
              hMuonEtaAntiDatap->Fill(fabs(lep->Eta()));
              doMET ? hDataMetm2d[it]->Fill(metCorrLep) : hDataMetm2d[it]->Fill(mtCorr);
              hMetmIsoValues[it]->Fill(relIso);
              break;
            }
          }
        }
      } else if(typev[ifile]==eWlnu ) {
        // std::cout << "doing signal" << std::endl;
        int bin=0;
        for(int i = 0; i <= hh_diff->GetNbinsX();++i){
          if(genVPt > hh_diff->GetBinLowEdge(i) && genVPt < hh_diff->GetBinLowEdge(i+1)){ bin = i; break; }
        }
        // double w2 =  hh_diff->GetBinContent(bin);
        // weight*=w2;
        hWlnuMet->Fill(metCorrMain,weight);
        if(q>0){
          hMuonEtaMCp->Fill(fabs(lep->Eta()),weight);
          hWlnuMetpPhi->Fill(metCorrMainPhi);
          doMET ? hWlnuMetp->Fill(metCorrMain,weight) : hWlnuMetp->Fill(mtCorr,weight);
          doMET ? hWlnuMetp2d[0]    ->Fill(metCorrMain, weight) : hWlnuMetp2d[0]    ->Fill(mtCorr, weight);
          doMET ? hWlnuMetp2dUncRecU[0][eta] ->Fill(metCorrEta ,  weight) : hWlnuMetp2dUncRecU[0][eta] ->Fill(mtCorr, weight);
          doMET ? hWlnuMetp2dUncRecU[0][keys]->Fill(metCorrKeys,  weight) : hWlnuMetp2dUncRecU[0][keys]->Fill(mtCorr, weight);
          doMET ? hWlnuMetp2dUncRecU[0][stat]->Fill(metCorrStat,  weight) : hWlnuMetp2dUncRecU[0][stat]->Fill(mtCorr, weight);
          doMET ? hWlnuMetp2dUncEffU[0][fsr] ->Fill(metCorrMain,  wgtFsr) : hWlnuMetp2dUncEffU[0][fsr] ->Fill(mtCorr, wgtFsr);
          doMET ? hWlnuMetp2dUncEffU[0][bkg] ->Fill(metCorrMain,  wgtBkg) : hWlnuMetp2dUncEffU[0][bkg] ->Fill(mtCorr, wgtBkg);
          doMET ? hWlnuMetp2dUncEffU[0][mc]  ->Fill(metCorrMain,  wgtMC ) : hWlnuMetp2dUncEffU[0][mc]  ->Fill(mtCorr, wgtMC );
          // // doMET ? hWlnuMetp2dKeys[0]->Fill(metCorrKeys, weight) : hWlnuMetp2dKeys[0]->Fill(mtCorr, weight);
          // doMET ? hWlnuMetp2dStat[0]->Fill(metCorrStat, weight) : hWlnuMetp2dStat[0]->Fill(mtCorr, weight);
        } else {
          hMuonEtaMCm->Fill(fabs(lep->Eta()),weight);
          hWlnuMetmPhi->Fill(metCorrMainPhi);
          doMET ? hWlnuMetm->Fill(metCorrMain,weight) : hWlnuMetm->Fill(mtCorr,weight);
          doMET ? hWlnuMetm2d[0]    ->Fill(metCorrMain,weight) : hWlnuMetm2d[0]    ->Fill(mtCorr,weight);
          doMET ? hWlnuMetm2dUncRecU[0][eta] ->Fill(metCorrEta , weight) : hWlnuMetm2dUncRecU[0][eta] ->Fill(mtCorr,weight);
          doMET ? hWlnuMetm2dUncRecU[0][keys]->Fill(metCorrKeys, weight) : hWlnuMetm2dUncRecU[0][keys]->Fill(mtCorr,weight);
          doMET ? hWlnuMetm2dUncRecU[0][stat]->Fill(metCorrStat, weight) : hWlnuMetm2dUncRecU[0][stat]->Fill(mtCorr,weight);
          doMET ? hWlnuMetm2dUncEffU[0][fsr] ->Fill(metCorrMain, wgtFsr) : hWlnuMetm2dUncEffU[0][fsr] ->Fill(mtCorr,wgtFsr);
          doMET ? hWlnuMetm2dUncEffU[0][bkg] ->Fill(metCorrMain, wgtBkg) : hWlnuMetm2dUncEffU[0][bkg] ->Fill(mtCorr,wgtBkg);
          doMET ? hWlnuMetm2dUncEffU[0][mc]  ->Fill(metCorrMain, wgtMC ) : hWlnuMetm2dUncEffU[0][mc]  ->Fill(mtCorr,wgtMC);
        }
      } else if(typev[ifile]==eWx) {
        // std::cout << "doing Wx" << std::endl;
        doMET ? hEWKMet->Fill(metCorrMain,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          doMET ? hEWKMetp->Fill(metCorrMain,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hWxMetp2d[0]     ->Fill(metCorrMain,weight) : hWxMetp2d[0]     ->Fill(mtCorr,weight);
          doMET ? hWxMetp2dUncRecU[0][eta] ->Fill(metCorrEta ,weight) : hWxMetp2dUncRecU[0][eta] ->Fill(mtCorr,weight);
          doMET ? hWxMetp2dUncRecU[0][keys]->Fill(metCorrKeys,weight) : hWxMetp2dUncRecU[0][keys]->Fill(mtCorr,weight);
          doMET ? hWxMetp2dUncRecU[0][stat]->Fill(metCorrStat,weight) : hWxMetp2dUncRecU[0][stat]->Fill(mtCorr,weight);
          doMET ? hWxMetp2dUncEffU[0][fsr] ->Fill(metCorrMain,wgtFsr) : hWxMetp2dUncEffU[0][fsr] ->Fill(mtCorr,wgtFsr);
          doMET ? hWxMetp2dUncEffU[0][bkg] ->Fill(metCorrMain,wgtBkg) : hWxMetp2dUncEffU[0][bkg] ->Fill(mtCorr,wgtBkg);
          doMET ? hWxMetp2dUncEffU[0][mc]  ->Fill(metCorrMain,wgtMC ) : hWxMetp2dUncEffU[0][mc]  ->Fill(mtCorr,wgtMC );
        } else {
          doMET ? hEWKMetm->Fill(metCorrMain,weight): hEWKMetm->Fill(mtCorr,weight);
          doMET ? hWxMetm2d[0]     ->Fill(metCorrMain,weight) : hWxMetm2d[0] ->Fill(mtCorr,weight);
          doMET ? hWxMetm2dUncRecU[0][eta] ->Fill(metCorrEta , weight) : hWxMetm2dUncRecU[0][eta] ->Fill(mtCorr,weight);
          doMET ? hWxMetm2dUncRecU[0][keys]->Fill(metCorrKeys, weight) : hWxMetm2dUncRecU[0][keys]->Fill(mtCorr,weight);
          doMET ? hWxMetm2dUncRecU[0][stat]->Fill(metCorrStat, weight) : hWxMetm2dUncRecU[0][stat]->Fill(mtCorr,weight);
          doMET ? hWxMetm2dUncEffU[0][fsr] ->Fill(metCorrMain, wgtFsr) : hWxMetm2dUncEffU[0][fsr] ->Fill(mtCorr,wgtFsr);
          doMET ? hWxMetm2dUncEffU[0][bkg] ->Fill(metCorrMain, wgtBkg) : hWxMetm2dUncEffU[0][bkg] ->Fill(mtCorr,wgtBkg);
          doMET ? hWxMetm2dUncEffU[0][mc]  ->Fill(metCorrMain, wgtMC ) : hWxMetm2dUncEffU[0][mc] ->Fill(mtCorr,wgtMC );
        }
      } else if(typev[ifile]==eZxx){
        // std::cout << "doing Zxx" << std::endl;
        doMET ? hEWKMet->Fill(metCorrMain,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          
 				  doMET ? hEWKMetp->Fill(metCorrMain,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hZxxMetp2d[0]    ->Fill(metCorrMain,weight) : hZxxMetp2d[0]    ->Fill(mtCorr,weight);
          doMET ? hZxxMetp2dUncRecU[0][eta] ->Fill(metCorrEta , weight) : hZxxMetp2dUncRecU[0][eta] ->Fill(mtCorr,weight);
          doMET ? hZxxMetp2dUncRecU[0][keys]->Fill(metCorrKeys, weight) : hZxxMetp2dUncRecU[0][keys]->Fill(mtCorr,weight);
          doMET ? hZxxMetp2dUncRecU[0][stat]->Fill(metCorrStat, weight) : hZxxMetp2dUncRecU[0][stat]->Fill(mtCorr,weight);
          doMET ? hZxxMetp2dUncEffU[0][fsr] ->Fill(metCorrMain, wgtFsr) : hZxxMetp2dUncEffU[0][fsr] ->Fill(mtCorr,wgtFsr);
          doMET ? hZxxMetp2dUncEffU[0][bkg] ->Fill(metCorrMain, wgtBkg) : hZxxMetp2dUncEffU[0][bkg] ->Fill(mtCorr,wgtBkg);
          doMET ? hZxxMetp2dUncEffU[0][mc]  ->Fill(metCorrMain, wgtMC ) : hZxxMetp2dUncEffU[0][mc]  ->Fill(mtCorr,wgtMC );
        } else {
          doMET ? hEWKMetm->Fill(metCorrMain,weight) :  hEWKMetm->Fill(mtCorr,weight);
          doMET ? hZxxMetm2d[0]    ->Fill(metCorrMain,weight) : hZxxMetm2d[0]    ->Fill(mtCorr,weight);
          doMET ? hZxxMetm2dUncRecU[0][eta] ->Fill(metCorrEta , weight) : hZxxMetm2dUncRecU[0][eta] ->Fill(mtCorr,weight);
          doMET ? hZxxMetm2dUncRecU[0][keys]->Fill(metCorrKeys, weight) : hZxxMetm2dUncRecU[0][keys]->Fill(mtCorr,weight);
          doMET ? hZxxMetm2dUncRecU[0][stat]->Fill(metCorrStat, weight) : hZxxMetm2dUncRecU[0][stat]->Fill(mtCorr,weight);
          doMET ? hZxxMetm2dUncEffU[0][fsr] ->Fill(metCorrMain, wgtFsr) : hZxxMetm2dUncEffU[0][fsr] ->Fill(mtCorr,wgtFsr);
          doMET ? hZxxMetm2dUncEffU[0][bkg] ->Fill(metCorrMain, wgtBkg) : hZxxMetm2dUncEffU[0][bkg] ->Fill(mtCorr,wgtBkg);
          doMET ? hZxxMetm2dUncEffU[0][mc]  ->Fill(metCorrMain, wgtMC ) : hZxxMetm2dUncEffU[0][mc]  ->Fill(mtCorr,wgtMC );
        }
      } else if(typev[ifile]==eDib) {
        // std::cout << "doing Dibosons" << std::endl;
        doMET ? hEWKMet->Fill(metCorrLep,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          // std::cout << "filling dib + " << std::endl;
          doMET ? hEWKMetp->Fill(metCorrLep,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hDibMetp2d[0]->Fill(metCorrLep,weight) : hDibMetp2d[0]->Fill(mtCorr,weight);
          doMET ? hDibMetp2dUncEffU[0][fsr]->Fill(metCorrLep,wgtFsr) : hDibMetp2dUncEffU[0][fsr]->Fill(mtCorr,wgtFsr);
          doMET ? hDibMetp2dUncEffU[0][bkg]->Fill(metCorrLep,wgtBkg) : hDibMetp2dUncEffU[0][bkg]->Fill(mtCorr,wgtBkg);
          doMET ? hDibMetp2dUncEffU[0][mc] ->Fill(metCorrLep,wgtMC ) : hDibMetp2dUncEffU[0][mc] ->Fill(mtCorr,wgtMC );
        } else {
          doMET ? hEWKMetm->Fill(metCorrLep,weight) : hEWKMetm->Fill(mtCorr,weight);
          doMET ? hDibMetm2d[0]->Fill(metCorrLep,weight) : hDibMetm2d[0]->Fill(mtCorr,weight);
          doMET ? hDibMetm2dUncEffU[0][fsr]->Fill(metCorrLep,wgtFsr) : hDibMetm2dUncEffU[0][fsr]->Fill(mtCorr,wgtFsr);
          doMET ? hDibMetm2dUncEffU[0][bkg]->Fill(metCorrLep,wgtBkg) : hDibMetm2dUncEffU[0][bkg]->Fill(mtCorr,wgtBkg);
          doMET ? hDibMetm2dUncEffU[0][mc] ->Fill(metCorrLep,wgtMC ) : hDibMetm2dUncEffU[0][mc] ->Fill(mtCorr,wgtMC );
        }
      } else if(typev[ifile]==eTtb) {
        // std::cout << "doing TTBar" << std::endl;
        doMET ? hEWKMet->Fill(metCorrLep,weight) : hEWKMet->Fill(mtCorr,weight);
        if(q>0){
          doMET ? hEWKMetp->Fill(metCorrLep,weight) : hEWKMetp->Fill(mtCorr,weight);
          doMET ? hTtbMetp2d[0]->Fill(metCorrLep,weight) : hTtbMetp2d[0]->Fill(mtCorr,weight);
          doMET ? hTtbMetp2dUncEffU[0][fsr]->Fill(metCorrLep,weight) : hTtbMetp2dUncEffU[0][fsr]->Fill(mtCorr,weight);
          doMET ? hTtbMetp2dUncEffU[0][bkg]->Fill(metCorrLep,wgtBkg) : hTtbMetp2dUncEffU[0][bkg]->Fill(mtCorr,wgtBkg);
          doMET ? hTtbMetp2dUncEffU[0][mc] ->Fill(metCorrLep,wgtMC ) : hTtbMetp2dUncEffU[0][mc] ->Fill(mtCorr,wgtMC );
        } else {
          doMET ? hEWKMetm->Fill(metCorrLep,weight) : hEWKMetm->Fill(mtCorr,weight);
          doMET ? hTtbMetm2dUncEffU[0][fsr]->Fill(metCorrLep,wgtFsr) : hTtbMetm2dUncEffU[0][fsr]->Fill(mtCorr,wgtFsr);
          doMET ? hTtbMetm2dUncEffU[0][bkg]->Fill(metCorrLep,wgtBkg) : hTtbMetm2dUncEffU[0][bkg]->Fill(mtCorr,wgtBkg);
          doMET ? hTtbMetm2dUncEffU[0][mc] ->Fill(metCorrLep,wgtMC ) : hTtbMetm2dUncEffU[0][mc] ->Fill(mtCorr,wgtMC );
        }
      } else if(typev[ifile]==eAntiWlnu){
        hAntiWlnuMet->Fill(metCorrMain,weight);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0) {              
              doMET ? hAntiWlnuMetp->Fill(metCorrMain,weight) : hAntiWlnuMetp->Fill(mtCorr,weight);
              doMET ? hWlnuMetp2d[it]    ->Fill(metCorrMain, weight) : hWlnuMetp2d[it]    ->Fill(mtCorr, weight);
              doMET ? hWlnuMetp2dUncRecU[it][eta] ->Fill(metCorrEta ,  weight) : hWlnuMetp2dUncRecU[it][eta] ->Fill(mtCorr, weight);
              doMET ? hWlnuMetp2dUncRecU[it][keys]->Fill(metCorrKeys,  weight) : hWlnuMetp2dUncRecU[it][keys]->Fill(mtCorr, weight);
              doMET ? hWlnuMetp2dUncRecU[it][stat]->Fill(metCorrStat,  weight) : hWlnuMetp2dUncRecU[it][stat]->Fill(mtCorr, weight);
              doMET ? hWlnuMetp2dUncEffU[it][fsr] ->Fill(metCorrMain,  wgtFsr) : hWlnuMetp2dUncEffU[it][fsr] ->Fill(mtCorr, wgtFsr);
              doMET ? hWlnuMetp2dUncEffU[it][bkg] ->Fill(metCorrMain,  wgtBkg) : hWlnuMetp2dUncEffU[it][bkg] ->Fill(mtCorr, wgtBkg);
              doMET ? hWlnuMetp2dUncEffU[it][mc]  ->Fill(metCorrMain,  wgtMC ) : hWlnuMetp2dUncEffU[it][mc]  ->Fill(mtCorr, wgtMC );
            } else {
              doMET ? hAntiWlnuMetm->Fill(metCorrMain,weight) : hAntiWlnuMetm->Fill(mtCorr,weight);
              doMET ? hWlnuMetm2d[it]    ->Fill(metCorrMain, weight) : hWlnuMetm2d[it]    ->Fill(mtCorr, weight);
              doMET ? hWlnuMetm2dUncRecU[it][eta] ->Fill(metCorrEta , weight) : hWlnuMetm2dUncRecU[it][eta] ->Fill(mtCorr,weight);
              doMET ? hWlnuMetm2dUncRecU[it][keys]->Fill(metCorrKeys, weight) : hWlnuMetm2dUncRecU[it][keys]->Fill(mtCorr,weight);
              doMET ? hWlnuMetm2dUncRecU[it][stat]->Fill(metCorrStat, weight) : hWlnuMetm2dUncRecU[it][stat]->Fill(mtCorr,weight);
              doMET ? hWlnuMetm2dUncEffU[it][fsr] ->Fill(metCorrMain, wgtFsr) : hWlnuMetm2dUncEffU[it][fsr] ->Fill(mtCorr,wgtFsr);
              doMET ? hWlnuMetm2dUncEffU[it][bkg] ->Fill(metCorrMain, wgtBkg) : hWlnuMetm2dUncEffU[it][bkg] ->Fill(mtCorr,wgtBkg);
              doMET ? hWlnuMetm2dUncEffU[it][mc]  ->Fill(metCorrMain, wgtMC ) : hWlnuMetm2dUncEffU[it][mc]  ->Fill(mtCorr,wgtMC );
            }
          }
        }
      } else if(typev[ifile]==eAntiWx){
        doMET ? hAntiEWKMet->Fill(metCorrMain,weight) : hAntiEWKMet->Fill(mtCorr,weight);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              doMET ? hAntiEWKMetp->Fill(metCorrMain,weight) :  hAntiEWKMetp->Fill(mtCorr,weight);
              doMET ? hWxMetp2d[it]     ->Fill(metCorrMain, weight) : hWxMetp2d[it]     ->Fill(mtCorr, weight);
              doMET ? hWxMetp2dUncRecU[it][eta] ->Fill(metCorrEta,  weight) : hWxMetp2dUncRecU[it][eta] ->Fill(mtCorr, weight);
              doMET ? hWxMetp2dUncRecU[it][keys]->Fill(metCorrKeys, weight) : hWxMetp2dUncRecU[it][keys]->Fill(mtCorr, weight);
              doMET ? hWxMetp2dUncRecU[it][stat]->Fill(metCorrStat, weight) : hWxMetp2dUncRecU[it][stat]->Fill(mtCorr, weight);
              doMET ? hWxMetp2dUncEffU[it][fsr] ->Fill(metCorrKeys, wgtFsr) : hWxMetp2dUncEffU[it][fsr] ->Fill(mtCorr, wgtFsr);
              doMET ? hWxMetp2dUncEffU[it][bkg] ->Fill(metCorrKeys, wgtBkg) : hWxMetp2dUncEffU[it][bkg] ->Fill(mtCorr, wgtBkg);
              doMET ? hWxMetp2dUncEffU[it][mc]  ->Fill(metCorrKeys, wgtMC ) : hWxMetp2dUncEffU[it][mc]  ->Fill(mtCorr, wgtMC );
            } else {
              doMET ? hAntiEWKMetm->Fill(metCorrMain,weight) : hAntiEWKMetm->Fill(mtCorr,weight);
              doMET ? hWxMetm2d[it]     ->Fill(metCorrMain, weight) : hWxMetm2d[it]     ->Fill(mtCorr, weight);
              doMET ? hWxMetm2dUncRecU[it][eta] ->Fill(metCorrEta ,  weight) : hWxMetm2dUncRecU[it][eta] ->Fill(mtCorr, weight);
              doMET ? hWxMetm2dUncRecU[it][keys]->Fill(metCorrKeys,  weight) : hWxMetm2dUncRecU[it][keys]->Fill(mtCorr, weight);
              doMET ? hWxMetm2dUncRecU[it][stat]->Fill(metCorrStat,  weight) : hWxMetm2dUncRecU[it][stat]->Fill(mtCorr, weight);
              doMET ? hWxMetm2dUncEffU[it][fsr] ->Fill(metCorrMain,  wgtFsr) : hWxMetm2dUncEffU[it][fsr] ->Fill(mtCorr, wgtFsr);
              doMET ? hWxMetm2dUncEffU[it][bkg] ->Fill(metCorrMain,  wgtBkg) : hWxMetm2dUncEffU[it][bkg] ->Fill(mtCorr, wgtBkg);
              doMET ? hWxMetm2dUncEffU[it][mc]  ->Fill(metCorrMain,  wgtMC ) : hWxMetm2dUncEffU[it][mc]  ->Fill(mtCorr, wgtMC );
            }
            break;
          }
        }
      } else if(typev[ifile]==eAntiZxx){
        doMET ? hAntiEWKMet->Fill(metCorrMain,weight) : hAntiEWKMet->Fill(mtCorr,weight);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill(metCorrMain,weight); 
              doMET ? hZxxMetp2d[it]    ->Fill(metCorrMain, weight) : hZxxMetp2d[it]->Fill(mtCorr, weight);
              doMET ? hZxxMetp2dUncRecU[it][eta] ->Fill(metCorrEta,  weight) : hZxxMetp2dUncRecU[it][eta] ->Fill(mtCorr, weight);
              doMET ? hZxxMetp2dUncRecU[it][keys]->Fill(metCorrEta,  weight) : hZxxMetp2dUncRecU[it][keys]->Fill(mtCorr, weight);
              doMET ? hZxxMetp2dUncRecU[it][stat]->Fill(metCorrEta,  weight) : hZxxMetp2dUncRecU[it][stat]->Fill(mtCorr, weight);
              doMET ? hZxxMetp2dUncEffU[it][fsr] ->Fill(metCorrEta,  wgtFsr) : hZxxMetp2dUncEffU[it][fsr] ->Fill(mtCorr, wgtFsr);
              doMET ? hZxxMetp2dUncEffU[it][bkg] ->Fill(metCorrEta,  wgtBkg) : hZxxMetp2dUncEffU[it][bkg] ->Fill(mtCorr, wgtBkg);
              doMET ? hZxxMetp2dUncEffU[it][mc]  ->Fill(metCorrEta,  wgtMC ) : hZxxMetp2dUncEffU[it][mc]  ->Fill(mtCorr, wgtMC );
            } else {
              hAntiEWKMetm->Fill(metCorrMain,weight); 
              doMET ? hZxxMetm2d[it]    ->Fill(metCorrMain, weight) : hZxxMetm2d[it]    ->Fill(mtCorr, weight);
              doMET ? hZxxMetm2dUncRecU[it][eta] ->Fill(metCorrEta ,  weight) : hZxxMetm2dUncRecU[it][eta] ->Fill(mtCorr, weight);
              doMET ? hZxxMetm2dUncRecU[it][keys]->Fill(metCorrKeys,  weight) : hZxxMetm2dUncRecU[it][keys]->Fill(mtCorr, weight);
              doMET ? hZxxMetm2dUncRecU[it][stat]->Fill(metCorrStat,  weight) : hZxxMetm2dUncRecU[it][stat]->Fill(mtCorr, weight);
              doMET ? hZxxMetm2dUncEffU[it][fsr] ->Fill(metCorrMain,  wgtFsr) : hZxxMetm2dUncEffU[it][fsr] ->Fill(mtCorr, wgtFsr);
              doMET ? hZxxMetm2dUncEffU[it][bkg] ->Fill(metCorrMain,  wgtBkg) : hZxxMetm2dUncEffU[it][bkg] ->Fill(mtCorr, wgtBkg);
              doMET ? hZxxMetm2dUncEffU[it][mc]  ->Fill(metCorrMain,  wgtMC ) : hZxxMetm2dUncEffU[it][mc]  ->Fill(mtCorr, wgtMC );
            }
          }
        }
      } else if(typev[ifile]==eAntiDib){
        doMET ? hAntiEWKMet->Fill(metCorrLep,weight) : hAntiEWKMet->Fill(mtCorr,weight);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill(metCorrLep,weight); 
              doMET ? hDibMetp2d[it]->Fill(metCorrLep, weight) : hDibMetp2d[it]->Fill(mtCorr, weight);
              doMET ? hDibMetp2dUncEffU[it][fsr]->Fill(metCorrLep, wgtFsr) : hDibMetp2dUncEffU[it][fsr]->Fill(mtCorr, wgtFsr);
              doMET ? hDibMetp2dUncEffU[it][bkg]->Fill(metCorrLep, wgtBkg) : hDibMetp2dUncEffU[it][bkg]->Fill(mtCorr, wgtBkg);
              doMET ? hDibMetp2dUncEffU[it][mc] ->Fill(metCorrLep, wgtMC ) : hDibMetp2dUncEffU[it][mc] ->Fill(mtCorr, wgtMC );
            } else {
              hAntiEWKMetm->Fill(metCorrLep,weight); 
              doMET ? hDibMetm2d[it]->Fill(metCorrLep, weight) : hDibMetm2d[it]->Fill(mtCorr, weight);
              doMET ? hDibMetm2dUncEffU[it][fsr]->Fill(metCorrLep, wgtFsr) : hDibMetm2dUncEffU[it][fsr]->Fill(mtCorr, wgtFsr);
              doMET ? hDibMetm2dUncEffU[it][bkg]->Fill(metCorrLep, wgtBkg) : hDibMetm2dUncEffU[it][bkg]->Fill(mtCorr, wgtBkg);
              doMET ? hDibMetm2dUncEffU[it][mc] ->Fill(metCorrLep, wgtMC ) : hDibMetm2dUncEffU[it][mc] ->Fill(mtCorr, wgtMC );
            }
          }
        }
      } else if(typev[ifile]==eAntiTtb){
        doMET ? hAntiEWKMet->Fill(metCorrLep,weight) : hAntiEWKMet->Fill(mtCorr,weight);
        for(int it=1; it < nIsoBins; ++it){
          if(relIso >= vIsoBins[it] && relIso < vIsoBins[it+1]) {
            if(q>0){
              hAntiEWKMetp->Fill(metCorrLep,weight); 
              doMET ? hTtbMetp2d[it]->Fill(metCorrLep, weight) : hTtbMetp2d[it]->Fill(mtCorr, weight);
              doMET ? hTtbMetp2dUncEffU[it][fsr]->Fill(metCorrLep, wgtFsr) : hTtbMetp2dUncEffU[it][fsr]->Fill(mtCorr, wgtFsr);
              doMET ? hTtbMetp2dUncEffU[it][bkg]->Fill(metCorrLep, wgtBkg) : hTtbMetp2dUncEffU[it][bkg]->Fill(mtCorr, wgtBkg);
              doMET ? hTtbMetp2dUncEffU[it][mc] ->Fill(metCorrLep, wgtMC ) : hTtbMetp2dUncEffU[it][mc] ->Fill(mtCorr, wgtMC );
            } else {
              hAntiEWKMetm->Fill(metCorrLep,weight); 
              doMET ? hTtbMetm2d[it]->Fill(metCorrLep, weight) : hTtbMetm2d[it]->Fill(mtCorr, weight);
              doMET ? hTtbMetm2dUncEffU[it][fsr]->Fill(metCorrLep, wgtFsr) : hTtbMetm2dUncEffU[it][fsr]->Fill(mtCorr, wgtFsr);
              doMET ? hTtbMetm2dUncEffU[it][bkg]->Fill(metCorrLep, wgtBkg) : hTtbMetm2dUncEffU[it][bkg]->Fill(mtCorr, wgtBkg);
              doMET ? hTtbMetm2dUncEffU[it][mc] ->Fill(metCorrLep, wgtMC ) : hTtbMetm2dUncEffU[it][mc] ->Fill(mtCorr, wgtMC );
            }
          }
        }
      } else if(typev[ifile]==eAntiQCD) {
        hAntiQCDMet->Fill(metCorrLep,weight);
        q>0 ? hAntiQCDMetp->Fill(metCorrLep,weight) : hAntiQCDMetm->Fill(metCorrLep,weight); 
      }
    }
  }
  
  for(int it = 0; it < nIsoBins; it++){
    
    hEWKMetp2d[it]->Add(hTtbMetp2d[it],1);
    hEWKMetp2d[it]->Add(hDibMetp2d[it],1);
    hEWKMetp2d[it]->Add(hZxxMetp2d[it],1);
    hEWKMetp2d[it]->Add(hWxMetp2d[it],1);
    
    hEWKMetm2d[it]->Add(hTtbMetm2d[it],1);
    hEWKMetm2d[it]->Add(hDibMetm2d[it],1);
    hEWKMetm2d[it]->Add(hZxxMetm2d[it],1);
    hEWKMetm2d[it]->Add(hWxMetm2d[it],1);
    
    for(int j=0; j< nUncEff; j++){
      hEWKMetp2dUncEffU[it][j]->Add(hTtbMetp2dUncEffU[it][j],1);
      hEWKMetp2dUncEffU[it][j]->Add(hDibMetp2dUncEffU[it][j],1);
      hEWKMetp2dUncEffU[it][j]->Add(hZxxMetp2dUncEffU[it][j],1);
      hEWKMetp2dUncEffU[it][j]->Add(hWxMetp2dUncEffU[it][j],1);
      
      hEWKMetm2dUncEffU[it][j]->Add(hTtbMetm2dUncEffU[it][j],1);
      hEWKMetm2dUncEffU[it][j]->Add(hDibMetm2dUncEffU[it][j],1);
      hEWKMetm2dUncEffU[it][j]->Add(hZxxMetm2dUncEffU[it][j],1);
      hEWKMetm2dUncEffU[it][j]->Add(hWxMetm2dUncEffU[it][j],1);
    }
    for(int j=0; j < nUncRec; j++){
      hEWKMetp2dUncRecU[it][j]->Add(hZxxMetp2dUncRecU[it][j],1);
      hEWKMetp2dUncRecU[it][j]->Add(hWxMetp2dUncRecU[it][j],1);
      
      hEWKMetm2dUncRecU[it][j]->Add(hZxxMetm2dUncRecU[it][j],1);
      hEWKMetm2dUncRecU[it][j]->Add(hWxMetm2dUncRecU[it][j],1);
    }
  }
  
  delete infile;
  infile=0, intree=0;   
   
    ofstream txtfile2;
    char txtfname2[100];
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
    
  
  ofstream txtfile3;
  char txtfname3[100];
  std::cout << "Printing We+. " << std::endl;
  sprintf(txtfname3,"%s/isoBins_Yields.txt",CPlot::sOutDir.Data());
  txtfile3.open(txtfname3);
  assert(txtfile3.is_open());
	txtfile3 << "W+ " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  Data# " << hDataMetp2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  Wsig# " << hWlnuMetp2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  ewk# "  << hEWKMetp2d[i]->Integral()  << std::endl;
	}
	
	txtfile3 << "W- " << std::endl << "data" << std::endl;
	for(int i = 0; i < nIsoBins; ++i){
		txtfile3 << "bin " << i << "  Data# " << hDataMetm2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  Wsig# " << hWlnuMetm2d[i]->Integral() << std::endl;
		txtfile3 << "bin " << i << "  ewk# "  << hEWKMetm2d[i]->Integral()  << std::endl;
	}

	txtfile3.close();
  
  
  std::cout << " ==== nEvents === " << std::endl;
  std::cout << "sig single  = " << hWlnuMetm->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWlnuMetm2d[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetm->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetm2d[0]->Integral() << std::endl;
  
  std::cout << "sig single  = " << hWlnuMetp->Integral() << std::endl;
  std::cout << "sig isobin  = " << hWlnuMetp2d[0]->Integral() << std::endl;
  std::cout << "dat single  = " << hDataMetp->Integral() << std::endl;
  std::cout << "dat isobin  = " << hDataMetp2d[0]->Integral() << std::endl;
  
  
 
 std::cout << "blah" << std::endl;
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
  
  std::cout << "blah here" << std::endl;
  
  TH1D *hIsoBinQCDp = (TH1D*) hDataMetp2d[1]->Clone("hQCDMetpBin");
  TH1D *hIsoBinQCDm = (TH1D*) hDataMetm2d[1]->Clone("hQCDMetmBin");
  hIsoBinQCDp->Add(hWlnuMetp2d[1],-1); hIsoBinQCDp->Add(hEWKMetp2d[1],-1);
  hIsoBinQCDm->Add(hWlnuMetm2d[1],-1); hIsoBinQCDm->Add(hEWKMetm2d[1],-1);
  
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
    for(int k=0; k < nUncRec; ++k){
      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetp2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnuMetp2dUncRecD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetm2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnuMetm2dUncRecD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWEwkMetpBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hEWKMetp2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKMetp2dUncRecD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hEWKMetm2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKMetm2dUncRecD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hWxMetp2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxMetp2dUncRecD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxMetp2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hWxMetm2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxMetm2dUncRecD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxMetm2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetp2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxMetp2dUncRecD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxMetp2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vUncRec[k].c_str());
      hh_diff =  (TH1D*)hZxxMetm2dUncRecU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxMetm2dUncRecD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxMetm2dUncRecD[j][k]->Add(hh_diff,-1); //delete hh_diff;
    }
    
    for(int k=0; k < nUncEff; ++k){
      sprintf(nname,"hWlnuMetpBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetp2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
      hWlnuMetp2dUncEffD[j][k] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWlnuMetmBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetm2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
      hWlnuMetm2dUncEffD[j][k] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWEwkMetpBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hEWKMetp2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
      hEWKMetp2dUncEffD[j][k] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hEwkMetmBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hEWKMetm2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
      hEWKMetm2dUncEffD[j][k] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetpBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hWxMetp2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
      hWxMetp2dUncEffD[j][k] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxMetp2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hWxMetmBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hWxMetm2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
      hWxMetm2dUncEffD[j][k] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxMetm2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetpBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetp2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
      hZxxMetp2dUncEffD[j][k] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxMetp2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hZxxMetmBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hZxxMetm2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
      hZxxMetm2dUncEffD[j][k] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxMetm2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hDibMetpBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hDibMetp2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hDibMetp2d[j],-1);
      hDibMetp2dUncEffD[j][k] = (TH1D*) hDibMetp2d[j]->Clone(nname); hDibMetp2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hDibMetmBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hDibMetm2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hDibMetm2d[j],-1);
      hDibMetm2dUncEffD[j][k] = (TH1D*) hDibMetm2d[j]->Clone(nname); hDibMetm2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hTtbMetpBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hWlnuMetp2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hTtbMetp2d[j],-1);
      hTtbMetp2dUncEffD[j][k] = (TH1D*) hTtbMetp2d[j]->Clone(nname); hTtbMetp2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
      
      sprintf(nname,"hTtbMetmBin%d_%sDown",j,vUncEff[k].c_str());
      hh_diff =  (TH1D*)hTtbMetm2dUncEffU[j][k]->Clone("diff"); hh_diff->Add(hTtbMetm2d[j],-1);
      hTtbMetm2dUncEffD[j][k] = (TH1D*) hTtbMetm2d[j]->Clone(nname); hTtbMetm2dUncEffD[j][k]->Add(hh_diff,-1); //delete hh_diff;
    }
    
    // TH1D *hh_diff =  (TH1D*)hWlnuMetp2dEta[j]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
    // hWlnuMetp2dEtaD[j] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetpBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hWlnuMetp2dKeys[j]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
    // hWlnuMetp2dKeysD[j] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetpBin%d_statDown",j);
    // hh_diff =  (TH1D*)hWlnuMetp2dStat[j]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
    // hWlnuMetp2dStatD[j] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dStatD[j]->Add(hh_diff,-1);// delete hh_diff;
    
    // sprintf(nname,"hWlnuMetpBin%d_mcDown",j);
    // // clone the up shape and subtract the main shape
    // hh_diff =  (TH1D*)hWlnuMetp2dEffMC[j]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
    // hWlnuMetp2dEffMCD[j] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dEffMCD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetpBin%d_bkgDown",j);
    // hh_diff =  (TH1D*)hWlnuMetp2dEffBkg[j]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
    // hWlnuMetp2dEffBkgD[j] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dEffBkgD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetpBin%d_fsrDown",j);
    // hh_diff =  (TH1D*)hWlnuMetp2dEffFsr[j]->Clone("diff"); hh_diff->Add(hWlnuMetp2d[j],-1);
    // hWlnuMetp2dEffFsrD[j] = (TH1D*) hWlnuMetp2d[j]->Clone(nname); hWlnuMetp2dEffFsrD[j]->Add(hh_diff,-1);// delete hh_diff;
    //////////////////////////
    //------------------------------------------------------------
    // sprintf(nname,"hEWKMetpBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hEWKMetp2dEta[j]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
    // hEWKMetp2dEtaD[j] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetpBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hEWKMetp2dKeys[j]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
    // hEWKMetp2dKeysD[j] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetpBin%d_statDown",j);
    // hh_diff =  (TH1D*)hEWKMetp2dStat[j]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
    // hEWKMetp2dStatD[j] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dStatD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetpBin%d_mcDown",j);
    // hh_diff =  (TH1D*)hEWKMetp2dEffMC[j]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
    // hEWKMetp2dEffMCD[j] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dEffMCD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetpBin%d_bkgDown",j);
    // hh_diff =  (TH1D*)hEWKMetp2dEffBkg[j]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
    // hEWKMetp2dEffBkgD[j] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dEffBkgD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetpBin%d_fsrDown",j);
    // hh_diff =  (TH1D*)hEWKMetp2dEffFsr[j]->Clone("diff"); hh_diff->Add(hEWKMetp2d[j],-1);
    // hEWKMetp2dEffFsrD[j] = (TH1D*) hEWKMetp2d[j]->Clone(nname); hEWKMetp2dEffFsrD[j]->Add(hh_diff,-1); //delete hh_diff;
    ////////////////////////////
    //------------------------------------------------------------
    // sprintf(nname,"hWxMetpBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hWxMetp2dEta[j]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
    // hWxMetp2dEtaD[j] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxMetp2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWxMetpBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hWxMetp2dKeys[j]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
    // hWxMetp2dKeysD[j] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxMetp2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWxMetpBin%d_statDown",j);
    // hh_diff =  (TH1D*)hWxMetp2dStat[j]->Clone("diff"); hh_diff->Add(hWxMetp2d[j],-1);
    // hWxMetp2dStatD[j] = (TH1D*) hWxMetp2d[j]->Clone(nname); hWxMetp2dStatD[j]->Add(hh_diff,-1); //delete hh_diff;
    // /////////////////////
    // //------------------------------------------------------------
    // sprintf(nname,"hZxxMetpBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hZxxMetp2dEta[j]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
    // hZxxMetp2dEtaD[j] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxMetp2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hZxxMetpBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hZxxMetp2dKeys[j]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
    // hZxxMetp2dKeysD[j] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxMetp2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hZxxMetpBin%d_statDown",j);
    // hh_diff =  (TH1D*)hZxxMetp2dStat[j]->Clone("diff"); hh_diff->Add(hZxxMetp2d[j],-1);
    // hZxxMetp2dStatD[j] = (TH1D*) hZxxMetp2d[j]->Clone(nname); hZxxMetp2dStatD[j]->Add(hh_diff,-1); //delete hh_diff;
    // ==============================================================================================================
    // ---------------------   W-   ------------------------------------
        // // clone the up shape and subtract the main shape
    // sprintf(nname,"hWlnuMetmBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hWlnuMetm2dEta[j]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
    // hWlnuMetm2dEtaD[j] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetmBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hWlnuMetm2dKeys[j]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
    // hWlnuMetm2dKeysD[j] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetmBin%d_statDown",j);
    // hh_diff =  (TH1D*)hWlnuMetm2dStat[j]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
    // hWlnuMetm2dStatD[j] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dStatD[j]->Add(hh_diff,-1);// delete hh_diff;
    
    // sprintf(nname,"hWlnuMetmBin%d_mcDown",j);
    // hh_diff =  (TH1D*)hWlnuMetm2dEffMC[j]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
    // hWlnuMetm2dEffMCD[j] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dEffMCD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetmBin%d_bkgDown",j);
    // hh_diff =  (TH1D*)hWlnuMetm2dEffBkg[j]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
    // hWlnuMetm2dEffBkgD[j] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dEffBkgD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWlnuMetmBin%d_fsrDown",j);
    // hh_diff =  (TH1D*)hWlnuMetm2dEffFsr[j]->Clone("diff"); hh_diff->Add(hWlnuMetm2d[j],-1);
    // hWlnuMetm2dEffFsrD[j] = (TH1D*) hWlnuMetm2d[j]->Clone(nname); hWlnuMetm2dEffFsrD[j]->Add(hh_diff,-1);// delete hh_diff;
    // //////////////////////////
    //------------------------------------------------------------
    // sprintf(nname,"hEWKMetmBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hEWKMetm2dEta[j]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
    // hEWKMetm2dEtaD[j] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetmBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hEWKMetm2dKeys[j]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
    // hEWKMetm2dKeysD[j] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetmBin%d_statDown",j);
    // hh_diff =  (TH1D*)hEWKMetm2dStat[j]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
    // hEWKMetm2dStatD[j] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dStatD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetmBin%d_mcDown",j);
    // hh_diff =  (TH1D*)hEWKMetm2dEffMC[j]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
    // hEWKMetm2dEffMCD[j] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dEffMCD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetmBin%d_bkgDown",j);
    // hh_diff =  (TH1D*)hEWKMetm2dEffBkg[j]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
    // hEWKMetm2dEffBkgD[j] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dEffBkgD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hEWKMetmBin%d_fsrDown",j);
    // hh_diff =  (TH1D*)hEWKMetm2dEffFsr[j]->Clone("diff"); hh_diff->Add(hEWKMetm2d[j],-1);
    // hEWKMetm2dEffFsrD[j] = (TH1D*) hEWKMetm2d[j]->Clone(nname); hEWKMetm2dEffFsrD[j]->Add(hh_diff,-1); //delete hh_diff;
    ////////////////////////////
    //------------------------------------------------------------
    // sprintf(nname,"hWxMetmBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hWxMetm2dEta[j]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
    // hWxMetm2dEtaD[j] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxMetm2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWxMetmBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hWxMetm2dKeys[j]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
    // hWxMetm2dKeysD[j] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxMetm2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hWxMetmBin%d_statDown",j);
    // hh_diff =  (TH1D*)hWxMetm2dStat[j]->Clone("diff"); hh_diff->Add(hWxMetm2d[j],-1);
    // hWxMetm2dStatD[j] = (TH1D*) hWxMetm2d[j]->Clone(nname); hWxMetm2dStatD[j]->Add(hh_diff,-1); //delete hh_diff;
    // /////////////////////
    // //------------------------------------------------------------
    // sprintf(nname,"hZxxMetmBin%d_etaDown",j);
    // hh_diff =  (TH1D*)hZxxMetm2dEta[j]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
    // hZxxMetm2dEtaD[j] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxMetm2dEtaD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hZxxMetmBin%d_keysDown",j);
    // hh_diff =  (TH1D*)hZxxMetm2dKeys[j]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
    // hZxxMetm2dKeysD[j] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxMetm2dKeysD[j]->Add(hh_diff,-1); //delete hh_diff;
    
    // sprintf(nname,"hZxxMetmBin%d_statDown",j);
    // hh_diff =  (TH1D*)hZxxMetm2dStat[j]->Clone("diff"); hh_diff->Add(hZxxMetm2d[j],-1);
    // hZxxMetm2dStatD[j] = (TH1D*) hZxxMetm2d[j]->Clone(nname); hZxxMetm2dStatD[j]->Add(hh_diff,-1); //delete hh_diff;
      
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
      
      //============== Recoil Up Shapes =====================
      // makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2dEta[j] ,WlnuMetpEta_ ,pdfWepEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2dKeys[j],WlnuMetpKeys_,pdfWepKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2dStat[j],WlnuMetpStat_,pdfWepStat_,pfmet,j,"_statUp");

      // makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2dEta[j] ,ewkMetpEta_ ,pdfEWKpEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2dKeys[j],ewkMetpKeys_,pdfEWKpKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2dStat[j],ewkMetpStat_,pdfEWKpStat_,pfmet,j,"_statUp");
      
      // makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2dEta[j]  ,wxMetpEta_ ,pdfWxpEta_  ,pfmet,j,"_etaUp");
      // makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2dKeys[j] ,wxMetpKeys_ ,pdfWxpKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2dStat[j] ,wxMetpStat_ ,pdfWxpStat_,pfmet,j,"_statUp");
      
      // makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2dEta[j] ,zxxMetpEta_ ,pdfZxxpEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2dKeys[j],zxxMetpKeys_,pdfZxxpKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2dStat[j],zxxMetpStat_,pdfZxxpStat_,pfmet,j,"_statUp");
      
      //============= Recoil Down Shapes ====================
      // makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2dEtaD[j] ,WlnuMetpEtaD_ ,pdfWepEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2dKeysD[j],WlnuMetpKeysD_,pdfWepKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("WlnuMETp","wep",hWlnuMetp2dStatD[j],WlnuMetpStatD_,pdfWepStatD_,pfmet,j,"_statDown");

      // makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2dEtaD[j] ,ewkMetpEtaD_ ,pdfEWKpEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2dKeysD[j],ewkMetpKeysD_,pdfEWKpKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("ewkMETp","ewkp",hEWKMetp2dStatD[j],ewkMetpStatD_,pdfEWKpStatD_,pfmet,j,"_statDown");
      
      // makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2dEtaD[j]  ,wxMetpEtaD_ ,pdfWxpEtaD_  ,pfmet,j,"_etaDown");
      // makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2dKeysD[j] ,wxMetpKeysD_ ,pdfWxpKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("wxMETp" ,"wxp" ,hWxMetp2dStatD[j] ,wxMetpStatD_ ,pdfWxpStatD_,pfmet,j,"_statDown");
      
      // makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2dEtaD[j] ,zxxMetpEtaD_ ,pdfZxxpEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2dKeysD[j],zxxMetpKeysD_,pdfZxxpKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("zxxMETp","zxxp",hZxxMetp2dStatD[j],zxxMetpStatD_,pdfZxxpStatD_,pfmet,j,"_statDown");
      
      j==0?makeDataHistPdf("qcdMETp","qcdp",hIsoBinQCDp,qcdMetp_,pdfQCDp_,pfmet,j,""):makeDataHistPdf("qcdMETp","qcdp",hDataMetp2d[j],qcdMetp_,pdfQCDp_,pfmet,j,"");
	  // ----------------------------------------- W- ---------------------------------
      makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2d[j],WlnuMetm_,pdfWem_,pfmet,j,"");
      makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2d[j],ewkMetm_,pdfEWKm_,pfmet,j,"");
      makeDataHistPdf("ttbMETm","ttbm",hTtbMetm2d[j],ttbMetm_,pdfTtbm_,pfmet,j,"");
      makeDataHistPdf("dibMETm","dibm",hDibMetm2d[j],dibMetm_,pdfDibm_,pfmet,j,"");
      makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2d[j] ,wxMetm_ ,pdfWxm_ ,pfmet,j,"");
      makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2d[j],zxxMetm_,pdfZxxm_,pfmet,j,"");
      
      //============== Recoil Up Shapes
      // makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2dEta[j] ,WlnuMetmEta_ ,pdfWemEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2dKeys[j],WlnuMetmKeys_,pdfWemKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2dStat[j],WlnuMetmStat_,pdfWemStat_,pfmet,j,"_statUp");
      
      // makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2dEta[j] ,ewkMetmEta_ ,pdfEWKmEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2dKeys[j],ewkMetmKeys_,pdfEWKmKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2dStat[j],ewkMetmStat_,pdfEWKmStat_,pfmet,j,"_statUp");
      
      // makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2dEta[j] ,wxMetmEta_ ,pdfWxmEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2dKeys[j],wxMetmKeys_,pdfWxmKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2dStat[j],wxMetmStat_,pdfWxmStat_,pfmet,j,"_statUp");
      
      // makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2dEta[j] ,zxxMetmEta_ ,pdfZxxmEta_ ,pfmet,j,"_etaUp");
      // makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2dKeys[j],zxxMetmKeys_,pdfZxxmKeys_,pfmet,j,"_keysUp");
      // makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2dStat[j],zxxMetmStat_,pdfZxxmStat_,pfmet,j,"_statUp");
      // //============== Recoil Down Shapes
      // makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2dEtaD[j] ,WlnuMetmEtaD_ ,pdfWemEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2dKeysD[j],WlnuMetmKeysD_,pdfWemKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("WlnuMETm","wem",hWlnuMetm2dStatD[j],WlnuMetmStatD_,pdfWemStatD_,pfmet,j,"_statDown");
      
      // makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2dEtaD[j] ,ewkMetmEtaD_ ,pdfEWKmEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2dKeysD[j],ewkMetmKeysD_,pdfEWKmKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("ewkMETm","ewkm",hEWKMetm2dStatD[j],ewkMetmStatD_,pdfEWKmStatD_,pfmet,j,"_statDown");
      
      // makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2dEtaD[j] ,wxMetmEtaD_ ,pdfWxmEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2dKeysD[j],wxMetmKeysD_,pdfWxmKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("wxMETm" ,"wxm" ,hWxMetm2dStatD[j],wxMetmStatD_,pdfWxmStatD_,pfmet,j,"_statDown");
      
      // makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2dEtaD[j] ,zxxMetmEtaD_ ,pdfZxxmEtaD_ ,pfmet,j,"_etaDown");
      // makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2dKeysD[j],zxxMetmKeysD_,pdfZxxmKeysD_,pfmet,j,"_keysDown");
      // makeDataHistPdf("zxxMETm","zxxm",hZxxMetm2dStatD[j],zxxMetmStatD_,pdfZxxmStatD_,pfmet,j,"_statDown");
      
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
  
  
    
  
  // Anti-Signal PDFs
  RooDataHist aWlnuMet ("aWlnuMET", "aWlnuMET", RooArgSet(pfmet),hAntiWlnuMet);  RooHistPdf apdfWm ("awm", "awm", pfmet,aWlnuMet, 1);
  RooDataHist aWlnuMetp("aWlnuMETp","aWlnuMETp",RooArgSet(pfmet),hAntiWlnuMetp); RooHistPdf apdfWmp("awmp","awmp",pfmet,aWlnuMetp,1);
  RooDataHist aWlnuMetm("aWlnuMETm","aWlnuMETm",RooArgSet(pfmet),hAntiWlnuMetm); RooHistPdf apdfWmm("awmm","awmm",pfmet,aWlnuMetm,1); 
  
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
  vector <CPepeModel2isobinsMuons*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  qcdp_[0] = new CPepeModel2isobinsMuons("qcdp2d0",pfmet, 0.0);
  qcdm_[0] = new CPepeModel2isobinsMuons("qcdm2d0",pfmet, 0.0);
  
  // quadratic dependence
  // vector <CPepeModel2isobinsQuad*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  // qcdp_[0] = new CPepeModel2isobinsQuad("qcdp2d0",pfmet, 0.075);
  // qcdm_[0] = new CPepeModel2isobinsQuad("qcdm2d0",pfmet, 0.075);
  
  // totally uncorrelated
  // vector <CPepeModel2*> qcdp_(nIsoBins), qcdm_(nIsoBins);
  // qcdp_[0] = new CPepeModel2("qcdp2d0",pfmet);
  // qcdm_[0] = new CPepeModel2("qcdm2d0",pfmet);
  std::cout << "start doing QCD" << std::endl;
  
  for (int j = 1; j < nIsoBins; ++j){
	  
      // quadratic isolation dependence
	  // sprintf(nname, "qcdp2d%d",j);
	  // qcdp_[j] = new CPepeModel2isobinsQuad(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->b1, qcdp_[0]->b2, qcdp_[0]->b3, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  // printf(nname, "qcdm2d%d",j);
	  // qcdm_[j] = new CPepeModel2isobinsQuad(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2,  qcdm_[0]->b1, qcdm_[0]->b2, qcdm_[0]->b3, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
	  
	  // // Original isolation linear dependence
      sprintf(nname, "qcdp2d%d",j);
      qcdp_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdp_[0]->c1, qcdp_[0]->c2, qcdp_[0]->c3, qcdp_[0]->d1, qcdp_[0]->d2, qcdp_[0]->d3);
	  
	  sprintf(nname, "qcdm2d%d",j);
      qcdm_[j] = new CPepeModel2isobinsMuons(nname,pfmet, (vIsoBins[j]+vIsoBins[j+1])/2, qcdm_[0]->c1, qcdm_[0]->c2, qcdm_[0]->c3, qcdm_[0]->d1, qcdm_[0]->d2, qcdm_[0]->d3);
  
      // sprintf(nname, "qcdp2d%d",j);
      // qcdp_[j] = new CPepeModel2(nname,pfmet);
	  
	  // sprintf(nname, "qcdm2d%d",j);
      // qcdm_[j] = new CPepeModel2(nname,pfmet);
  
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

		 (doMET&&(!doTemplate)) ? pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*(qcdp_[j]->model)),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j])) : pdfMetp_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWep_[j],*pdfEWKp_[j],*pdfQCDp_[j]),RooArgList(*nSigp_[j],*nEWKp_[j],*nQCDp_[j]));
		 
		 sprintf(nname,"pdfWem%d",j);
		 (doMET&&(!doTemplate)) ? pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*(qcdm_[j]->model)),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j])) : pdfMetm_[j] = new RooAddPdf(nname,nname,RooArgList(*pdfWem_[j],*pdfEWKm_[j],*pdfQCDm_[j]),RooArgList(*nSigm_[j],*nEWKm_[j],*nQCDm_[j]));		
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

  RooDataHist combDatap("dataTotalp","dataTotalp", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetp_));
  RooDataHist combDatam("dataTotalm","dataTotalm", RooArgList(pfmet), Index(rooCat2dTest),Import(dataMetm_));
  
  RooSimultaneous simPdfp("simPdfp","simultaneous pdf W+",rooCat2dTest) ;
  RooSimultaneous simPdfm("simPdfm","simultaneous pdf W-",rooCat2dTest) ;
  for(int i = 0; i < nIsoBins; ++i ){
	  sprintf(nname, "isop%d", i); simPdfp.addPdf(*pdfMetp_[i],nname);
	  sprintf(nname, "isom%d", i); simPdfm.addPdf(*pdfMetm_[i],nname);
  }

  
  cout << "Starting values for Wlnu yields: " << endl;
  cout << "Selected: " << hDataMet->Integral() << endl;
  cout << "   sig: " << hWlnuMet->Integral() << endl;
  cout << "   EWK: " << hEWKMet->Integral() << endl;
  cout << "   qcd: " << hDataMet->Integral()-hWlnuMet->Integral()-hEWKMet->Integral() << endl;

  cout << "Starting values for Wlnu_p yields: " << endl;
  cout << "   sig: " << hWlnuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hDataMetp->Integral()-hWlnuMetp->Integral()-hEWKMetp->Integral() << endl;

  cout << "Starting values for Wlnu_m yields: " << endl;
  cout << "   sig: " << hWlnuMetm->Integral() << endl;
  cout << "   EWK: " << hEWKMetm->Integral() << endl;
  cout << "   qcd: " << hDataMetm->Integral()-hWlnuMetm->Integral()-hEWKMetm->Integral() << endl;
  
  
  cout << "Starting values for AntiWlnu yields: " << endl;
  cout << "Selected: " << hAntiDataMet->Integral() << endl;
  cout << "   sig: " << hAntiWlnuMet->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMet->Integral() << endl;
  cout << "   qcd: " << hAntiDataMet->Integral()-hAntiWlnuMet->Integral()-hAntiEWKMet->Integral() << endl;

  cout << "Starting values for AntiWlnu_p yields: " << endl;
  cout << "   sig: " << hWlnuMetp->Integral() << endl;
  cout << "   EWK: " << hEWKMetp->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetp->Integral()-hAntiWlnuMetp->Integral()-hAntiEWKMetp->Integral() << endl;

  cout << "Starting values for AntiWlnu_m yields: " << endl;
  cout << "   sig: " << hAntiWlnuMetm->Integral() << endl;
  cout << "   EWK: " << hAntiEWKMetm->Integral() << endl;
  cout << "   qcd: " << hAntiDataMetm->Integral()-hAntiWlnuMetm->Integral()-hAntiEWKMetm->Integral() << endl;

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
  
//   RooRealVar pepe1Pdf_qcdp_norm("pepe1Pdf_qcdp_norm","pepe1Pdf_qcdp_norm",40000,0,100000);
//   RooRealVar pepe1Pdf_qcdm_norm("pepe1Pdf_qcdm_norm","pepe1Pdf_qcdm_norm",40000,0,100000);
  TString histfname = outputDir + TString("/Wlnu_Histograms.root");
  TFile *histFile = new TFile(histfname,"RECREATE");
  histFile->cd();
  for(int j = 0; j < nIsoBins; ++j){
    
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
    hIsoBinQCDp->Write();
    hIsoBinQCDm->Write();
    
    for(int k=0; k < nUncRec; k++){
      hWlnuMetp2dUncRecU[j][k]->Write();
      hWlnuMetm2dUncRecU[j][k]->Write();
      
      hWlnuMetp2dUncRecD[j][k]->Write();
      hWlnuMetm2dUncRecD[j][k]->Write();
      
      hEWKMetp2dUncRecU[j][k]->Write();
      hEWKMetm2dUncRecU[j][k]->Write();
      
      hEWKMetp2dUncRecD[j][k]->Write();
      hEWKMetm2dUncRecD[j][k]->Write();
      
      hWxMetp2dUncRecU[j][k]->Write();
      hWxMetm2dUncRecU[j][k]->Write();
      
      hWxMetp2dUncRecD[j][k]->Write();
      hWxMetm2dUncRecD[j][k]->Write();
      
      hZxxMetp2dUncRecU[j][k]->Write();
      hZxxMetm2dUncRecU[j][k]->Write();
      
      hZxxMetp2dUncRecD[j][k]->Write();
      hZxxMetm2dUncRecD[j][k]->Write();
      
    }
    
    for(int k=0; k < nUncEff; k++){
      hWlnuMetp2dUncEffU[j][k]->Write();
      hWlnuMetm2dUncEffU[j][k]->Write();
      
      hWlnuMetp2dUncEffD[j][k]->Write();
      hWlnuMetm2dUncEffD[j][k]->Write();
      
      hEWKMetp2dUncEffU[j][k]->Write();
      hEWKMetm2dUncEffU[j][k]->Write();
      
      hEWKMetp2dUncEffD[j][k]->Write();
      hEWKMetm2dUncEffD[j][k]->Write();
      
      hWxMetp2dUncEffU[j][k]->Write();
      hWxMetm2dUncEffU[j][k]->Write();
      
      hWxMetp2dUncEffD[j][k]->Write();
      hWxMetm2dUncEffD[j][k]->Write();
      
      hZxxMetp2dUncEffU[j][k]->Write();
      hZxxMetm2dUncEffU[j][k]->Write();
      
      hZxxMetp2dUncEffD[j][k]->Write();
      hZxxMetm2dUncEffD[j][k]->Write();
      
      hDibMetp2dUncEffU[j][k]->Write();
      hDibMetm2dUncEffU[j][k]->Write();
      
      hDibMetp2dUncEffD[j][k]->Write();
      hDibMetm2dUncEffD[j][k]->Write();
      
      hTtbMetp2dUncEffU[j][k]->Write();
      hTtbMetm2dUncEffU[j][k]->Write();
      
      hTtbMetp2dUncEffD[j][k]->Write();
      hTtbMetm2dUncEffD[j][k]->Write();
    }
    
    
    // hWlnuMetp2dEta[j]->Write();
    // hWlnuMetp2dKeys[j]->Write();
    // hWlnuMetp2dStat[j]->Write();
    
    // hWlnuMetp2dEtaD[j]->Write();
    // hWlnuMetp2dKeysD[j]->Write();
    // hWlnuMetp2dStatD[j]->Write();
    
    // hEWKMetp2dEta[j]->Write();
    // hEWKMetp2dKeys[j]->Write();
    // hEWKMetp2dStat[j]->Write();
    
    // hEWKMetp2dEtaD[j]->Write();
    // hEWKMetp2dKeysD[j]->Write();
    // hEWKMetp2dStatD[j]->Write();
    
    // hZxxMetp2dEta[j]->Write();
    // hZxxMetp2dKeys[j]->Write();
    // hZxxMetp2dStat[j]->Write();
    
    // hZxxMetp2dEtaD[j]->Write();
    // hZxxMetp2dKeysD[j]->Write();
    // hZxxMetp2dStatD[j]->Write();
    
    // hWxMetp2dEta[j]->Write();
    // hWxMetp2dKeys[j]->Write();
    // hWxMetp2dStat[j]->Write();
    
    // hWxMetp2dEtaD[j]->Write();
    // hWxMetp2dKeysD[j]->Write();
    // hWxMetp2dStatD[j]->Write();
    
    // // hWlnuMetm2dEta[j]->Write();
    // // hWlnuMetm2dKeys[j]->Write();
    // // hWlnuMetm2dStat[j]->Write();
    
    // // hWlnuMetm2dEtaD[j]->Write();
    // // hWlnuMetm2dKeysD[j]->Write();
    // // hWlnuMetm2dStatD[j]->Write();
    
    // hEWKMetm2dEta[j]->Write();
    // hEWKMetm2dKeys[j]->Write();
    // hEWKMetm2dStat[j]->Write();
    
    // hEWKMetm2dEtaD[j]->Write();
    // hEWKMetm2dKeysD[j]->Write();
    // hEWKMetm2dStatD[j]->Write();
    
    // hZxxMetm2dEta[j]->Write();
    // hZxxMetm2dKeys[j]->Write();
    // hZxxMetm2dStat[j]->Write();
    
    // hZxxMetm2dEtaD[j]->Write();
    // hZxxMetm2dKeysD[j]->Write();
    // hZxxMetm2dStatD[j]->Write();
    
    // hWxMetm2dEta[j]->Write();
    // hWxMetm2dKeys[j]->Write();
    // hWxMetm2dStat[j]->Write();
    
    // hWxMetm2dEtaD[j]->Write();
    // hWxMetm2dKeysD[j]->Write();
    // hWxMetm2dStatD[j]->Write();
    
    
  }
  histFile->Write();
  histFile->Close();


  // RooWorkspace combine_workspace("combine_workspace");
  // combine_workspace.import(dataMet);
  // combine_workspace.import(dataMetp);
  // combine_workspace.import(dataMetm);
  // combine_workspace.import(pepe2Pdf_qcdp_norm);
  // combine_workspace.import(pepe2Pdf_qcdm_norm);

  // combine_workspace.import(pdfWm);
  // combine_workspace.import(pdfWmp);
  // combine_workspace.import(pdfWmm);
  // combine_workspace.import(pdfEWK);
  // combine_workspace.import(pdfEWKp);
  // combine_workspace.import(pdfEWKm);
// //   combine_workspace.import(*(qcd.model));
// //   combine_workspace.import(*(qcdp.model));
// //   combine_workspace.import(*(qcdm.model));

  // combine_workspace.import(pdfQCD);
  // combine_workspace.import(pdfQCDp);
  // combine_workspace.import(pdfQCDm);
  
  // combine_workspace.import(antiMet);
  // combine_workspace.import(antiMetp);
  // combine_workspace.import(antiMetm);
  // combine_workspace.import(pepe2Pdf_aqcdp_norm);
  // combine_workspace.import(pepe2Pdf_aqcdm_norm);
  
  // combine_workspace.import(apdfWm);
  // combine_workspace.import(apdfWmp);
  // combine_workspace.import(apdfWmm);
  // combine_workspace.import(apdfEWK);
  // combine_workspace.import(apdfEWKp);
  // combine_workspace.import(apdfEWKm);
  // combine_workspace.import(*(aqcd.model));
  // //combine_workspace.import(qcdpn);
  // combine_workspace.import(*(aqcdp.model));
  // combine_workspace.import(*(aqcdm.model));

  // sprintf(nname,"%s/Wlnu_pdfTemplates.root",CPlot::sOutDir.Data());
  // combine_workspace.writeToFile(nname);
  
  
   // // separate file for the binned workspace just so i don't fuck up the original one
  // RooWorkspace binned_workspace("binned_workspace");
  // // loop through the number of bins, import the appropriate pdf for each one
  // for(int j = 0; j < nIsoBins; ++j){
      // // binned_workspace.import(dataMet);
      // binned_workspace.import(*dataMetpHist_[j]);
      // binned_workspace.import(*dataMetmHist_[j]);
	  // // QCD normalization RooRealVars
      // binned_workspace.import(*pepe2Pdf_qcdp_norms_[j]);
      // binned_workspace.import(*pepe2Pdf_qcdm_norms_[j]);
	  // // QCD shapes
      // // binned_workspace.import(*(qcd.model));
      // binned_workspace.import(*(qcdp_[j]->model));
      // binned_workspace.import(*(qcdm_[j]->model));
      // std::cout << "importing   " << qcdm_[j]->model->GetName() << std::endl;
	  // //The DataHist
	  // // combine_workspace.import(*wenuMet_[j]);
	  // binned_workspace.import(*WlnuMetp_[j]);
	  // binned_workspace.import(*WlnuMetm_[j]);
      // // binned_workspace.import(pdfWe);
      // // binned_workspace.import(*pdfWep_[j]);
      // // binned_workspace.import(*pdfWem_[j]);
	  // // The DataHist
	  // // combine_workspace.import(*ewkMet_[j]);
	  // binned_workspace.import(*ewkMetp_[j]);
	  // binned_workspace.import(*ewkMetm_[j]);
      
    // binned_workspace.import(*wxMetp_[j]);
	  // binned_workspace.import(*wxMetm_[j]);
      
    // binned_workspace.import(*zxxMetp_[j]);
	  // binned_workspace.import(*zxxMetm_[j]);
      
    // binned_workspace.import(*ttbMetp_[j]);
	  // binned_workspace.import(*ttbMetm_[j]);
      
    // binned_workspace.import(*dibMetp_[j]);
	  // binned_workspace.import(*dibMetm_[j]);
    
    // binned_workspace.import(*WlnuMetpEta_[j]);
	  // binned_workspace.import(*WlnuMetmEta_[j]);
    // binned_workspace.import(*WlnuMetpKeys_[j]);
	  // binned_workspace.import(*WlnuMetmKeys_[j]);
    // binned_workspace.import(*WlnuMetpStat_[j]);
	  // binned_workspace.import(*WlnuMetmStat_[j]);
    
    // binned_workspace.import(*WlnuMetpEtaD_[j]);
	  // binned_workspace.import(*WlnuMetmEtaD_[j]);
    // binned_workspace.import(*WlnuMetpKeysD_[j]);
	  // binned_workspace.import(*WlnuMetmKeysD_[j]);
    // binned_workspace.import(*WlnuMetpStatD_[j]);
	  // binned_workspace.import(*WlnuMetmStatD_[j]);
    
    // binned_workspace.import(*ewkMetpEta_[j]);
	  // binned_workspace.import(*ewkMetmEta_[j]);
    // binned_workspace.import(*ewkMetpKeys_[j]);
	  // binned_workspace.import(*ewkMetmKeys_[j]);
    // binned_workspace.import(*ewkMetpStat_[j]);
	  // binned_workspace.import(*ewkMetmStat_[j]);
    
    // binned_workspace.import(*ewkMetpEtaD_[j]);
	  // binned_workspace.import(*ewkMetmEtaD_[j]);
    // binned_workspace.import(*ewkMetpKeysD_[j]);
	  // binned_workspace.import(*ewkMetmKeysD_[j]);
    // binned_workspace.import(*ewkMetpStatD_[j]);
	  // binned_workspace.import(*ewkMetmStatD_[j]);
    
    // binned_workspace.import(*wxMetpEta_[j]);
	  // binned_workspace.import(*wxMetmEta_[j]);
    // binned_workspace.import(*wxMetpKeys_[j]);
	  // binned_workspace.import(*wxMetmKeys_[j]);
    // binned_workspace.import(*wxMetpStat_[j]);
	  // binned_workspace.import(*wxMetmStat_[j]);
    
    // binned_workspace.import(*wxMetpEtaD_[j]);
	  // binned_workspace.import(*wxMetmEtaD_[j]);
    // binned_workspace.import(*wxMetpKeysD_[j]);
	  // binned_workspace.import(*wxMetmKeysD_[j]);
    // binned_workspace.import(*wxMetpStatD_[j]);
	  // binned_workspace.import(*wxMetmStatD_[j]);
    
    // binned_workspace.import(*zxxMetpEta_[j]);
	  // binned_workspace.import(*zxxMetmEta_[j]);
    // binned_workspace.import(*zxxMetpKeys_[j]);
	  // binned_workspace.import(*zxxMetmKeys_[j]);
    // binned_workspace.import(*zxxMetpStat_[j]);
	  // binned_workspace.import(*zxxMetmStat_[j]);
    
    // binned_workspace.import(*zxxMetpEtaD_[j]);
	  // binned_workspace.import(*zxxMetmEtaD_[j]);
    // binned_workspace.import(*zxxMetpKeysD_[j]);
	  // binned_workspace.import(*zxxMetmKeysD_[j]);
    // binned_workspace.import(*zxxMetpStatD_[j]);
	  // binned_workspace.import(*zxxMetmStatD_[j]);
    
	  // //The Pdfs
      // // binned_workspace.import(pdfEWK);
      // // binned_workspace.import(*pdfEWKp_[j]);
      // // binned_workspace.import(*pdfEWKm_[j]);
	  
	  // binned_workspace.import(*hDataMetp2d[j]);
	  // binned_workspace.import(*hDataMetm2d[j]);
	  // binned_workspace.import(*hWlnuMetp2d[j]);
	  // binned_workspace.import(*hWlnuMetm2d[j]);
	  // binned_workspace.import(*hEWKMetp2d[j]);
	  // binned_workspace.import(*hEWKMetm2d[j]);
    // binned_workspace.import(*hWxMetp2d[j]);
	  // binned_workspace.import(*hWxMetm2d[j]);
    // binned_workspace.import(*hZxxMetp2d[j]);
	  // binned_workspace.import(*hZxxMetm2d[j]);
    // binned_workspace.import(*hDibMetp2d[j]);
	  // binned_workspace.import(*hDibMetm2d[j]);
    // binned_workspace.import(*hTtbMetp2d[j]);
	  // binned_workspace.import(*hTtbMetm2d[j]);
	  // // binned_workspace.import(*hDataMetp2d[j]);
	  

  // // }// end of loop, now save the workspace
  // sprintf(nname, "%s/Wlnu_pdfTemplates_binned.root",CPlot::sOutDir.Data());
  // binned_workspace.writeToFile(nname);
  
  

  
  
    RooGaussian constm("constm","constm",nEWKm,RooConst(hEWKMetm->Integral()),RooConst(0.15*hEWKMetm->Integral()));
    RooGaussian constp("constp","constp",nEWKp,RooConst(hEWKMetp->Integral()),RooConst(0.15*hEWKMetp->Integral()));
    
    RooGaussian constantim("constantim","constantim",nAntiSigm,RooConst(hAntiWlnuMetm->Integral()),RooConst(0.15*hAntiWlnuMetm->Integral()));
    RooGaussian constantip("constantip","constantip",nAntiSigp,RooConst(hAntiWlnuMetp->Integral()),RooConst(0.15*hAntiWlnuMetp->Integral()));
	

  
     // had commented out the Min2/Strat2 when running some of the other options
  // RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),RooFit::Strategy(2)/*,Minimizer("Minuit2","minimize")*/,ExternalConstraints(RooArgList(const_wxp,const_zxxp,const_dibp,const_ttbp)),PrintEvalErrors(-1));
  // RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),RooFit::Strategy(2)/*,Minimizer("Minuit2","minimize")*/,ExternalConstraints(RooArgList(const_wxm,const_zxxm,const_dibm,const_ttbm)),PrintEvalErrors(-1));
  
  RooFitResult *fitResp2dCatTest = simPdfp.fitTo(combDatap,Extended(),Save(kTRUE),ExternalConstraints(RooArgSet(constp_sr/*,constp_cr*/)),RooFit::Strategy(2),Minos(kTRUE),/*Minimizer("Minuit2","minimize"),*/PrintEvalErrors(-1));
  RooFitResult *fitResm2dCatTest = simPdfm.fitTo(combDatam,Extended(),Save(kTRUE),ExternalConstraints(RooArgSet(constm_sr/*,constm_cr*/)),RooFit::Strategy(2),Minos(kTRUE),/*Minimizer("Minuit2","minimize"),*/PrintEvalErrors(-1));
  
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle); hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWlnuMetp->GetBinError(ibin));}
  std::cout << nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal() << std::endl;
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle); hMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", pfmet)); // why did we not just clone the original histogram...
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWlnuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); hMetmDiff->SetMarkerSize(0.9);
  
  // the diff hists for the anti-selection
  
  TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram("hPdfAntiMet", pfmet));
  
  hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hPdfAntiMet,"hAntiMetDiff");
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfAntiMetp = (TH1D*)(apdfMetp.createHistogram("hPdfAntiMetp", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetp->GetNbinsX(); ++ibin){hPdfAntiMetp->SetBinError(ibin, hAntiWlnuMetp->GetBinError(ibin));}
  hPdfAntiMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hPdfAntiMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hPdfAntiMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hPdfAntiMetm = (TH1D*)(apdfMetm.createHistogram("hPdfAntiMetm", pfmet));
   for(int ibin = 1; ibin < hPdfAntiMetm->GetNbinsX(); ++ibin){hPdfAntiMetm->SetBinError(ibin, hAntiWlnuMetm->GetBinError(ibin));}
  hPdfAntiMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hPdfAntiMetm->Integral());
  TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hPdfAntiMetm,"hAntiMetmDiff");
  hAntiMetmDiff->SetMarkerStyle(kFullCircle);
  hAntiMetmDiff->SetMarkerSize(0.9);
    
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
    for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWlnuMetp2d[i]->GetBinError(ibin));}
    std::cout << nSigp_[i]->getVal()+nEWKp_[i]->getVal()+nQCDp_[i]->getVal() << std::endl;
    hPdfMetp->Scale((nSigp_[i]->getVal()+nEWKp_[i]->getVal()+nQCDp_[i]->getVal())/hPdfMetp->Integral());
    TH1D *hMetpDiff = makeDiffHist(hDataMetp2d[i],hPdfMetp,"hMetpDiff");
    hMetpDiff->SetMarkerStyle(kFullCircle); hMetpDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    std::cout << "set up diff plot #" << i << std::endl;
    TH1D *hPdfMetm = (TH1D*)(pdfMetm_[i]->createHistogram("hPdfMetm", pfmet));
    for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWlnuMetm2d[i]->GetBinError(ibin));}
    std::cout << nSigm_[i]->getVal()+nEWKm_[i]->getVal()+nQCDm_[i]->getVal() << std::endl;
    hPdfMetm->Scale((nSigm_[i]->getVal()+nEWKm_[i]->getVal()+nQCDm_[i]->getVal())/hPdfMetm->Integral());
    TH1D *hMetmDiff = makeDiffHist(hDataMetm2d[i],hPdfMetm,"hMetmDiff");
    hMetmDiff->SetMarkerStyle(kFullCircle); hMetmDiff->SetMarkerSize(0.9);
    std::cout << "did diff " <<  i << std::endl;

    sprintf(nname,"isop%d",i); sprintf(plotname, "wep_fitmetp_bin%i",i);
    drawWMetPlots(plotname, hMetpDiff, pfmet, dataMetp_[nname], pdfMetp_[i], pdfEWKp_[i], doTemplate?(RooAbsPdf*)pdfQCDp_[i]:(RooAbsPdf*)qcdp_[i]->model, pdfWep_[i], lumitext, hDataMetp2d[i]);
    
    sprintf(nname,"isom%d",i); sprintf(plotname, "wem_fitmetm_bin%i",i);
    drawWMetPlots(plotname, hMetmDiff, pfmet, dataMetm_[nname], pdfMetm_[i], pdfEWKm_[i], doTemplate?(RooAbsPdf*)pdfQCDm_[i]:(RooAbsPdf*)qcdm_[i]->model, pdfWem_[i], lumitext, hDataMetm2d[i]);


    chi2probp = hDataMetp2d[i]->Chi2Test(hPdfMetp,"PUW");
    chi2ndfp  = hDataMetp2d[i]->Chi2Test(hPdfMetp,"CHI2/NDFUW");
    ksprobp   = hDataMetp2d[i]->KolmogorovTest(hPdfMetp);
    ksprobpep = hDataMetp2d[i]->KolmogorovTest(hPdfMetp,"DX"); 
   
	
    chi2probm = hDataMetm2d[i]->Chi2Test(hPdfMetm,"PUW");
    chi2ndfm  = hDataMetm2d[i]->Chi2Test(hPdfMetm,"CHI2/NDFUW");
    ksprobm   = hDataMetm2d[i]->KolmogorovTest(hPdfMetm);
    ksprobpem = hDataMetm2d[i]->KolmogorovTest(hPdfMetm,"DX"); 
	
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
  txtfile << "Data+: " << hDataMetp2d[0]->Integral() << endl;
  txtfile << "Data-: " << hDataMetm2d[0]->Integral() << endl;
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
  double yscale=0.2;
  
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
