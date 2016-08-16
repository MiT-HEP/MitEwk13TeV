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
#include <TMath.h> // ROOT math library
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
#include "TGraphAsymmErrors.h"

#include "BaconAna/DataFormats/interface/TGenParticle.hh"  

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
// #include "../Utils/RecoilCorrector_asym2.hh"    // class to handle recoil corrections for MET
#include "../Utils/LeptonCorr.hh"         // Scale and resolution corrections

// #include "ZBackgrounds.hh"

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);
TGraphAsymmErrors* makeShapeErrorBand(TH1 *hUp, TH1 *hDown);
TGraphAsymmErrors* TH1TOTGraphAsymmErrors(TH1 *h1);
TGraphAsymmErrors* myMakeBandSymmetric(TGraphAsymmErrors* g0, TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);
TGraphAsymmErrors* myTGraphErrorsDivide_noErrGraph2(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe);

void apply_postfit_shape(TH1D *histogram, TH1D *up, TH1D *down, double shift, double shift_err);

void apply_postfit_shape(TH1D *histogram, TH1D *up1, TH1D *down1, TH1D *up2, TH1D *down2, double shift, double shift_err);

// make webpage
void makeHTML(const TString outDir);


//=== MAIN MACRO ================================================================================================= 

void postFitWm(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("postFitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS   = 75;
  const Double_t METMAX  = 150;
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;

//   TString pufname = "../Tools/pileup_weights_2015B.root";

  // file format for output plots
  const TString format("png"); 


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  // Read in the Wmunu_pdfTemplates file here and start to get pdfs.
  TString inputName = "./Wmunu_pileup_simultaneous_check/Wmunu_pdfTemplates.root";
//   TString inputName = "./Wmunu_etaBins/Wmunu_pdfTemplates_3.root";
  TFile *inWmunuShapes = new TFile(inputName); assert(inWmunuShapes);
  RooWorkspace* combine_workspace = (RooWorkspace*) inWmunuShapes->Get("combine_workspace");
  combine_workspace->Print();
  // Get the pfmet variable from the workspace
  RooRealVar *pfmet =  combine_workspace->var("pfmet");
  
  // --------------------------------------------------
  // Start getting the PDFs and turning them into histograms ? (wut)
  // Do data first
  RooAbsData *data  = combine_workspace->data("dataMet");
  RooAbsData *datap = combine_workspace->data("dataMetp");
  RooAbsData *datam = combine_workspace->data("dataMetm");
  //  Make Histograms
  TH1D *hDataMet  = (TH1D*) data->createHistogram("pfmet",NBINS);  hDataMet->Sumw2();
  TH1D *hDataMetp = (TH1D*) datap->createHistogram("pfmet",NBINS);hDataMetp->Sumw2();
  TH1D *hDataMetm = (TH1D*) datam->createHistogram("pfmet",NBINS);hDataMetm->Sumw2();
  
  std::cout << "----DATA EVENT COUNTS----- " << std::endl;
  std::cout << "TOTAL:  " << hDataMet->Integral() << std::endl;
  std::cout << "W+:  " << hDataMetp->Integral() << std::endl;
  std::cout << "W-:  " << hDataMetm->Integral() << std::endl;
  
  //  Signal MC
  RooAbsPdf *wmet  = combine_workspace->pdf("wm");
  RooAbsPdf *wmetp = combine_workspace->pdf("wmp");
  RooAbsPdf *wmetm = combine_workspace->pdf("wmm");
  // Make Histograms
  TH1D *hWmunuMet  = (TH1D*) wmet->createHistogram("pfmet",NBINS);  hWmunuMet->Sumw2();
  TH1D *hWmunuMetp = (TH1D*) wmetp->createHistogram("pfmet",NBINS);hWmunuMetp->Sumw2();
  TH1D *hWmunuMetm = (TH1D*) wmetm->createHistogram("pfmet",NBINS);hWmunuMetm->Sumw2();
 
  
  RooAbsData *wmetData = combine_workspace->embeddedData("wmunuMET");
  RooAbsData *wmetpData = combine_workspace->embeddedData("wmunuMETp");
  RooAbsData *wmetmData = combine_workspace->embeddedData("wmunuMETm");
  
  TH1D *hWmunuMetWerr  = (TH1D*) wmetData->createHistogram("pfmet",NBINS);  hWmunuMetWerr->Sumw2();
  TH1D *hWmunuMetpWerr  = (TH1D*) wmetpData->createHistogram("pfmet",NBINS);  hWmunuMetpWerr->Sumw2();
  TH1D *hWmunuMetmWerr  = (TH1D*) wmetmData->createHistogram("pfmet",NBINS);  hWmunuMetmWerr->Sumw2();
  
  // Electroweak background contributions
  RooAbsPdf *ewk  = combine_workspace->pdf("ewk");
  RooAbsPdf *ewkp = combine_workspace->pdf("ewkp");
  RooAbsPdf *ewkm = combine_workspace->pdf("ewkm");
  // Make histograms
  TH1D *hEWKMet  = (TH1D*) ewk->createHistogram("pfmet",NBINS);  hEWKMet->Sumw2();
  TH1D *hEWKMetp = (TH1D*) ewkp->createHistogram("pfmet",NBINS);hEWKMetp->Sumw2();
  TH1D *hEWKMetm = (TH1D*) ewkm->createHistogram("pfmet",NBINS);hEWKMetm->Sumw2();
  
  // Shape variation, recoil uncertainty - upwards
  RooAbsPdf *ewk_pu  = combine_workspace->pdf("ewkm_PileupUp");
  RooAbsPdf *ewkp_pu = combine_workspace->pdf("ewkp_PileupUp");
  RooAbsPdf *ewkm_pu = combine_workspace->pdf("ewkm_PileupUp");
  // Make histograms
  TH1D *hEWKMet_PileupUp  = (TH1D*) ewk_pu->createHistogram("pfmet",NBINS);  hEWKMet_PileupUp->Sumw2();
  TH1D *hEWKMetp_PileupUp = (TH1D*) ewkp_pu->createHistogram("pfmet",NBINS); hEWKMetp_PileupUp->Sumw2();
  TH1D *hEWKMetm_PileupUp = (TH1D*) ewkm_pu->createHistogram("pfmet",NBINS); hEWKMetm_PileupUp->Sumw2();
  
  // Shape variation, recoil uncertainty - downwards
  RooAbsPdf *ewk_pd  = combine_workspace->pdf("ewkm_PileupDown");
  RooAbsPdf *ewkp_pd = combine_workspace->pdf("ewkp_PileupDown");
  RooAbsPdf *ewkm_pd = combine_workspace->pdf("ewkm_PileupDown");
  // Make histograms
  TH1D *hEWKMet_PileupDown  = (TH1D*) ewk_pd->createHistogram("pfmet",NBINS);  hEWKMet_PileupDown->Sumw2();
  TH1D *hEWKMetp_PileupDown = (TH1D*) ewkp_pd->createHistogram("pfmet",NBINS); hEWKMetp_PileupDown->Sumw2();
  TH1D *hEWKMetm_PileupDown = (TH1D*) ewkm_pd->createHistogram("pfmet",NBINS); hEWKMetm_PileupDown->Sumw2();
  
  // Shape variation, recoil uncertainty - upwards
  RooAbsPdf *wmet_ru  = combine_workspace->pdf("wm_RecoilUp");
  RooAbsPdf *wmetp_ru = combine_workspace->pdf("wmp_RecoilUp");
  RooAbsPdf *wmetm_ru = combine_workspace->pdf("wmm_RecoilUp");
  // Make histograms
  TH1D *hWmunuMet_RecoilUp  = (TH1D*) wmet_ru->createHistogram("pfmet",NBINS);   hWmunuMet_RecoilUp->Sumw2();
  TH1D *hWmunuMetp_RecoilUp = (TH1D*) wmetp_ru->createHistogram("pfmet",NBINS); hWmunuMetp_RecoilUp->Sumw2();
  TH1D *hWmunuMetm_RecoilUp = (TH1D*) wmetm_ru->createHistogram("pfmet",NBINS); hWmunuMetm_RecoilUp->Sumw2();
  
  // Shape variation, recoil uncertainty - downwards
  RooAbsPdf *wmet_rd  = combine_workspace->pdf("wm_RecoilDown");
  RooAbsPdf *wmetp_rd = combine_workspace->pdf("wmp_RecoilDown");
  RooAbsPdf *wmetm_rd = combine_workspace->pdf("wmm_RecoilDown");
  // Make histograms
  TH1D *hWmunuMet_RecoilDown  = (TH1D*) wmet_rd->createHistogram("pfmet",NBINS); hWmunuMet_RecoilDown->Sumw2();
  TH1D *hWmunuMetp_RecoilDown = (TH1D*) wmetp_rd->createHistogram("pfmet",NBINS); hWmunuMetp_RecoilDown->Sumw2();
  TH1D *hWmunuMetm_RecoilDown = (TH1D*) wmetm_rd->createHistogram("pfmet",NBINS); hWmunuMetm_RecoilDown->Sumw2();
  
  // Shape variation, recoil uncertainty - upwards
  RooAbsPdf *wmet_pu  = combine_workspace->pdf("wm_PileupUp");
  RooAbsPdf *wmetp_pu = combine_workspace->pdf("wmp_PileupUp");
  RooAbsPdf *wmetm_pu = combine_workspace->pdf("wmm_PileupUp");
  // Make histograms
  TH1D *hWmunuMet_PileupUp  = (TH1D*) wmet_pu->createHistogram("pfmet",NBINS);   hWmunuMet_PileupUp->Sumw2();
  TH1D *hWmunuMetp_PileupUp = (TH1D*) wmetp_pu->createHistogram("pfmet",NBINS); hWmunuMetp_PileupUp->Sumw2();
  TH1D *hWmunuMetm_PileupUp = (TH1D*) wmetm_pu->createHistogram("pfmet",NBINS); hWmunuMetm_PileupUp->Sumw2();
  
  // Shape variation, recoil uncertainty - downwards
  RooAbsPdf *wmet_pd  = combine_workspace->pdf("wm_PileupDown");
  RooAbsPdf *wmetp_pd = combine_workspace->pdf("wmp_PileupDown");
  RooAbsPdf *wmetm_pd = combine_workspace->pdf("wmm_PileupDown");
  // Make histograms
  TH1D *hWmunuMet_PileupDown  = (TH1D*) wmet_pd->createHistogram("pfmet",NBINS); hWmunuMet_PileupDown->Sumw2();
  TH1D *hWmunuMetp_PileupDown = (TH1D*) wmetp_pd->createHistogram("pfmet",NBINS); hWmunuMetp_PileupDown->Sumw2();
  TH1D *hWmunuMetm_PileupDown = (TH1D*) wmetm_pd->createHistogram("pfmet",NBINS); hWmunuMetm_PileupDown->Sumw2();
  
  
  // ======= The Anti-selections =============
    // Do data first
  RooAbsData *anti  = combine_workspace->data("antiMet");
  RooAbsData *antip = combine_workspace->data("antiMetp");
  RooAbsData *antim = combine_workspace->data("antiMetm");
  //  Make Histograms
  TH1D *hAntiDataMet  = (TH1D*) anti->createHistogram("pfmet",NBINS);  hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetp = (TH1D*) antip->createHistogram("pfmet",NBINS); hAntiDataMetp->Sumw2();
  TH1D *hAntiDataMetm = (TH1D*) antim->createHistogram("pfmet",NBINS); hAntiDataMetm->Sumw2();
  
  //  Signal MC
  RooAbsPdf *awmet  = combine_workspace->pdf("awm");
  RooAbsPdf *awmetp = combine_workspace->pdf("awmp");
  RooAbsPdf *awmetm = combine_workspace->pdf("awmm");
  // Make Histograms
  TH1D *hAntiWmunuMet  = (TH1D*) awmet->createHistogram("pfmet",NBINS);  hAntiWmunuMet->Sumw2();
  TH1D *hAntiWmunuMetp = (TH1D*) awmetp->createHistogram("pfmet",NBINS); hAntiWmunuMetp->Sumw2();
  TH1D *hAntiWmunuMetm = (TH1D*) awmetm->createHistogram("pfmet",NBINS); hAntiWmunuMetm->Sumw2();
  
  RooAbsData *awmetData = combine_workspace->embeddedData("awmunuMET");
  RooAbsData *awmetpData = combine_workspace->embeddedData("awmunuMETp");
  RooAbsData *awmetmData = combine_workspace->embeddedData("awmunuMETm");
  
  TH1D *hAntiWmunuMetWerr  = (TH1D*) awmetData->createHistogram("pfmet",NBINS);    hAntiWmunuMetWerr->Sumw2();
  TH1D *hAntiWmunuMetpWerr  = (TH1D*) awmetpData->createHistogram("pfmet",NBINS);  hAntiWmunuMetpWerr->Sumw2();
  TH1D *hAntiWmunuMetmWerr  = (TH1D*) awmetmData->createHistogram("pfmet",NBINS);  hAntiWmunuMetmWerr->Sumw2();
  
  // Electroweak background contributions
  RooAbsPdf *aewk  = combine_workspace->pdf("aewk");
  RooAbsPdf *aewkp = combine_workspace->pdf("aewkp");
  RooAbsPdf *aewkm = combine_workspace->pdf("aewkm");
  // Make histograms
  TH1D *hAntiEWKMet  = (TH1D*) aewk->createHistogram("pfmet",NBINS);  hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp = (TH1D*) aewkp->createHistogram("pfmet",NBINS); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm = (TH1D*) aewkm->createHistogram("pfmet",NBINS); hAntiEWKMetm->Sumw2();
  
  
  // Input the  parameters for signal and background yields
  // These are hardcoded since we get them from previous steps
  RooRealVar nSig("nSig","nSig",1000);
  RooRealVar nQCD("nQCD","nQCD",1000);
  RooRealVar nEWK("nEWK","nEWK",1000);
  
  RooRealVar nAntiSig("nAntiSig","nAntiSig",1000);
  RooRealVar nAntiQCD("nAntiQCD","nAntiQCD",1000);
  RooRealVar nAntiEWK("nAntiEWK","nAntiEWK",1000);


  RooRealVar nSigp("nSigp","nSigp",1.0004*9522749.21);
  RooRealVar nQCDp("nQCDp","nQCDp",10);
  RooRealVar nEWKp("nEWKp","nEWKp",1.0004*883693.8833);
  //   
  RooRealVar nSigm("nSigm","nSigm",0.9994*7363369.521);
  RooRealVar nQCDm("nQCDm","nQCDm",10);
  RooRealVar nEWKm("nEWKm","nEWKm",0.9994*794150.8014); 
  
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",1.0910*25167.44273);
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",10);
  RooRealVar nAntiEWKp("nAntiEWKp","nAntiEWKp",1.0910*4200.392376);
  //   
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",1.0464*16618.25939);
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",10);
  RooRealVar nAntiEWKm("nAntiEWKm","nAntiEWKm",1.0464*3472.668905); 
  
//   RooRealVar nSigp("nSigp","nSigp",1.0157*9380080.387);
//   RooRealVar nQCDp("nQCDp","nQCDp",10);
//   RooRealVar nEWKp("nEWKp","nEWKp",1.0157*870454.4749);
//   //   
//   RooRealVar nSigm("nSigm","nSigm",0.9850*7441924.78);
//   RooRealVar nQCDm("nQCDm","nQCDm",10);
//   RooRealVar nEWKm("nEWKm","nEWKm",0.9850*802623.1076); 
//   
//   RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",1.0157*247914.7684);
//   RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",10);
//   RooRealVar nAntiEWKp("nAntiEWKp","nAntiEWKp",1.0157*48324.53414);
//   //   
//   RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",0.6329*224590.0848);
//   RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",10);
//   RooRealVar nAntiEWKm("nAntiEWKm","nAntiEWKm",0.6329*56698.27263); 


  double nRecoilWmp = 0.2665;
  double nRecoilWmm = 0.0254;

  double nPileupWmp = 0.7986;
  double nPileupWmm = -0.8398;
  
  
  apply_postfit_shape(hWmunuMetp, hWmunuMetp_RecoilUp, hWmunuMetp_RecoilDown,hWmunuMetp_PileupUp, hWmunuMetp_PileupDown, nRecoilWmp, nPileupWmp);
  apply_postfit_shape(hWmunuMetm, hWmunuMetm_RecoilUp, hWmunuMetm_RecoilDown,hWmunuMetm_PileupUp, hWmunuMetm_PileupDown, nRecoilWmm, nPileupWmm);
  
  apply_postfit_shape(hEWKMetp, hEWKMetp_PileupUp, hEWKMetp_PileupDown, nPileupWmp, 0);
  apply_postfit_shape(hEWKMetm, hEWKMetm_PileupUp, hEWKMetm_PileupDown, nPileupWmm, 0);
  
  
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(*pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", *pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(*pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",*pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(*pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",*pfmet,wmunuMetm,1); 
  
  RooDataHist wmunuMet_RecoilUp("wmunuMET_RecoilUp", "wmunuMET_RecoilUp", RooArgSet(*pfmet),hWmunuMet_RecoilUp);  
  RooHistPdf pdfWm_RecoilUp("wm_RecoilUp", "wm_RecoilUp", *pfmet,wmunuMet_RecoilUp, 1);
  RooDataHist wmunuMetp_RecoilUp("wmunuMETp_RecoilUp","wmunuMETp_RecoilUp",RooArgSet(*pfmet),hWmunuMetp_RecoilUp); 
  RooHistPdf pdfWmp_RecoilUp("wmp_RecoilUp","wmp_RecoilUp",*pfmet,wmunuMetp_RecoilUp,1);
  RooDataHist wmunuMetm_RecoilUp("wmunuMETm_RecoilUp","wmunuMETm_RecoilUp",RooArgSet(*pfmet),hWmunuMetm_RecoilUp); 
  RooHistPdf pdfWmm_RecoilUp("wmm_RecoilUp","wmm_RecoilUp",*pfmet,wmunuMetm_RecoilUp,1); 
  RooDataHist wmunuMet_RecoilDown("wmunuMET_RecoilDown", "wmunuMET_RecoilDown", RooArgSet(*pfmet),hWmunuMet_RecoilDown);  
  RooHistPdf pdfWm_RecoilDown("wm_RecoilDown", "wm_RecoilDown", *pfmet,wmunuMet_RecoilDown, 1);
  RooDataHist wmunuMetp_RecoilDown("wmunuMETp_RecoilDown","wmunuMETp_RecoilDown",RooArgSet(*pfmet),hWmunuMetp_RecoilDown); 
  RooHistPdf pdfWmp_RecoilDown("wmp_RecoilDown","wmp_RecoilDown",*pfmet,wmunuMetp_RecoilDown,1);
  RooDataHist wmunuMetm_RecoilDown("wmunuMETm_RecoilDown","wmunuMETm_RecoilDown",RooArgSet(*pfmet),hWmunuMetm_RecoilDown); 
  RooHistPdf pdfWmm_RecoilDown("wmm_RecoilDown","wmm_RecoilDown",*pfmet,wmunuMetm_RecoilDown,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(*pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", *pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(*pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",*pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(*pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",*pfmet,ewkMetm,1); 
  
  // QCD Pdfs
  CPepeModel2 qcd("qcd",*pfmet);
  CPepeModel2 qcdp("qcdp",*pfmet);
  CPepeModel2 qcdm("qcdm",*pfmet);//,qcdp.sigma);
  
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(*pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(*pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(*pfmet),hDataMetm);
  
    // EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(*pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", *pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(*pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",*pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(*pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",*pfmet,aewkMetm,1); 
  
  // QCD Pdfs
  CPepeModel2 aqcd("aqcd",*pfmet);
  CPepeModel2 aqcdp("aqcdp",*pfmet);
  CPepeModel2 aqcdm("aqcdm",*pfmet);//,qcdp.sigma);
  
  
  RooDataHist awmunuMet ("awmunuMET", "awmunuMET", RooArgSet(*pfmet),hAntiWmunuMet);  RooHistPdf apdfWm ("awm", "awm", *pfmet,awmunuMet, 1);
  RooDataHist awmunuMetp("awmunuMETp","awmunuMETp",RooArgSet(*pfmet),hAntiWmunuMetp); RooHistPdf apdfWmp("awmp","awmp",*pfmet,awmunuMetp,1);
  RooDataHist awmunuMetm("awmunuMETm","awmunuMETm",RooArgSet(*pfmet),hAntiWmunuMetm); RooHistPdf apdfWmm("awmm","awmm",*pfmet,awmunuMetm,1); 
  
  // Signal + Background PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWm,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWmp,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWmm,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));

  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(*pfmet),hAntiDataMet);
  RooDataHist antiMetp("antiMetp","antiMetp",RooArgSet(*pfmet),hAntiDataMetp);
  RooDataHist antiMetm("antiMetm","antiMetm",RooArgSet(*pfmet),hAntiDataMetm);
  

  std::cout  << "Selectedp: " << hDataMetp->Integral() << endl;
  std::cout  << "Selectedm: " << hDataMetm->Integral() << endl;

// // Fixed Tail on QCD
//   qcdp.a1->setVal(0.1439);
//   qcdm.a1->setVal(0.1427);
//   qcdp.sigma->setVal(15.1711);
//   qcdm.sigma->setVal(15.2956);
//   nQCDp.setVal(1057310.7764);
//   nQCDm.setVal(1050856.9819);
//   
  qcdp.a1->setVal(0.0378);
  qcdp.a2->setVal(3.3264);
  qcdp.a3->setVal(243.3918);
  qcdm.a1->setVal(0.0386);
  qcdm.a2->setVal(3.3491);
  qcdm.a3->setVal(246.4277);
  nQCDp.setVal(1055089.3766);
  nQCDm.setVal(1052941.8978);
  
  aqcdp.a1->setVal(0.0378);
  aqcdp.a2->setVal(6.6133);
  aqcdp.a3->setVal(278.1532);
  aqcdm.a1->setVal(0.0386);
  aqcdm.a2->setVal(6.5470);
  aqcdm.a3->setVal(283.5787);
  nAntiQCDp.setVal(974322.6621);
  nAntiQCDm.setVal(975626.1161);

//   qcdp.a1->setVal(0.1352);
//   qcdm.a1->setVal(0.1274);
//   qcdp.sigma->setVal(16.3930);
//   qcdm.sigma->setVal(16.0603);
//   nQCDp.setVal(1160976.1229);
//   nQCDm.setVal(1061211.7311);
//   
//   aqcdp.a1->setVal(0.1352);
//   aqcdm.a1->setVal(0.1274);
//   aqcdp.sigma->setVal(18.5697);
//   aqcdm.sigma->setVal(18.9457);
//   nAntiQCDp.setVal(734720.0639);
//   nAntiQCDm.setVal(737052.2311);



  
  std::cout  << "EWKp: " << nEWKp.getVal() << endl;
  std::cout  << "EWKm: " << nEWKm.getVal() << endl;


  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  
  for(int ibin = 1; ibin < hWmunuMet->GetNbinsX(); ++ibin){
//     std::cout << "bin error " << hWmunuMetp->GetBinError(ibin) << std::endl;
    hWmunuMet->SetBinError(ibin, hWmunuMetWerr->GetBinError(ibin));
    hWmunuMetp->SetBinError(ibin, hWmunuMetpWerr->GetBinError(ibin));
    hWmunuMetm->SetBinError(ibin, hWmunuMetmWerr->GetBinError(ibin));
    
    hAntiWmunuMetp->SetBinError(ibin, hAntiWmunuMetpWerr->GetBinError(ibin));
    hAntiWmunuMetm->SetBinError(ibin, hAntiWmunuMetmWerr->GetBinError(ibin));
  }
  
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", *pfmet));
  for(int ibin = 1; ibin < hPdfMet->GetNbinsX(); ++ibin){hPdfMet->SetBinError(ibin, hWmunuMet->GetBinError(ibin));}
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", *pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);

  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", *pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
  
     
  TH1D *hAntiPdfMetp = (TH1D*)(apdfMetp.createHistogram("hAntiPdfMetp", *pfmet));
  for(int ibin = 1; ibin < hAntiPdfMetp->GetNbinsX(); ++ibin){hAntiPdfMetp->SetBinError(ibin, hAntiWmunuMetp->GetBinError(ibin));}
  hAntiPdfMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hAntiPdfMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hAntiPdfMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);

  TH1D *hAntiPdfMetm = (TH1D*)(apdfMetm.createHistogram("hAntiPdfMetm", *pfmet));
  for(int ibin = 1; ibin < hAntiPdfMetm->GetNbinsX(); ++ibin){hAntiPdfMetm->SetBinError(ibin, hAntiWmunuMetm->GetBinError(ibin));}
  hAntiPdfMetm->Scale((nAntiSigm.getVal()+nAntiEWKm.getVal()+nAntiQCDm.getVal())/hAntiPdfMetm->Integral());
  TH1D *hAntiMetmDiff = makeDiffHist(hAntiDataMetm,hAntiPdfMetm,"hAntiMetmDiff");
  hAntiMetmDiff->SetMarkerStyle(kFullCircle); 
  hAntiMetmDiff->SetMarkerSize(0.9);

  
 
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  
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
  else         sprintf(lumitext,"%.0f pb^{-1}  (13 TeV)",lumi);
  
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
   
  //
  // W MET plot
  //
  RooPlot *wmframe = pfmet->frame(Bins(NBINS)); 
  wmframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(wmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(wmframe,LineColor(linecolorW));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
  pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfWm)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("fitmet",wmframe,"","",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox("#bf{CMS}",0.64,0.80,0.72,0.88,0);
  plotMet.AddTextBox("#it{Preliminary}",0.72,0.80,0.88,0.86,0);
  plotMet.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kFALSE,format,1);

  CPlot plotMetDiff("fitmet","","#slash{E}_{T} [GeV]","#chi");
  //CPlot plotMetDiff("fitmet","","mT [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-8,8);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  
  plotMet.SetName("fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
       
  //
  // W+ MET plot
  //
  RooPlot *wmpframe = pfmet->frame(Bins(NBINS));
  wmpframe->GetYaxis()->SetNdivisions(505);
  wmpframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wmpframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,LineColor(linecolorW));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("fitmetp",wmpframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
//   plotMetpDiff.AddGraph(gMetpDiff,"Data","PE2",3,1,3);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.2,0.2);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  plotMetp.Draw(c,kTRUE,"pdf",1);
  
  // Anti- WM plot
  RooPlot *awmpframe = pfmet->frame(Bins(NBINS));
  awmpframe->GetYaxis()->SetNdivisions(505);
  awmpframe->GetXaxis()->SetLabelOffset(2.0);
  antiMetp.plotOn(awmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetp.plotOn(awmpframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetp.plotOn(awmpframe,LineColor(linecolorW));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),LineColor(linecolorEWK));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(*(aqcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(*(aqcdp.model))),LineColor(linecolorQCD));
  apdfMetp.plotOn(awmpframe,Components(RooArgSet(apdfWmp)),LineColor(linecolorW),LineStyle(2));
  antiMetp.plotOn(awmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
  CPlot plotAntiMetp("fitantimetp",awmpframe,"","",ylabel);
  plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrow#mu^{+}#nu","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotAntiMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetp.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetpDiff("fitantimetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetpDiff->GetYaxis()->SetLabelSize(0.11);
//   plotMetpDiff.AddGraph(gMetpDiff,"Data","PE2",3,1,3);
  plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
  plotAntiMetpDiff.SetYRange(-0.2,0.2);
  plotAntiMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetpDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotAntiMetpDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
  plotAntiMetpDiff.Draw(c,kTRUE,format,2);
  plotAntiMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotAntiMetp.SetName("fitantimetplog");
  plotAntiMetp.SetLogy();
  plotAntiMetp.SetYRange(1e-3*(hAntiDataMetp->GetMaximum()),10*(hAntiDataMetp->GetMaximum()));
  plotAntiMetp.Draw(c,kTRUE,format,1);
  plotAntiMetp.Draw(c,kTRUE,"pdf",1);
  
  //
  // W- MET plot
  //
  RooPlot *wmmframe = pfmet->frame(Bins(NBINS)); 
  wmmframe->GetYaxis()->SetNdivisions(505);
  wmmframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wmmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,LineColor(linecolorW));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
  pdfMetm.plotOn(wmmframe,Components(RooArgSet(pdfWmm)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("fitmetm",wmmframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetm.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetm.Draw(c,kFALSE,format,1);

  CPlot plotMetmDiff("fitmetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetmDiff.SetYRange(-0.2,0.2);
  plotMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
//   plotMetmDiff.AddGraph(gMetmDiff,"Err","FPE2",kGreen,1,kGreen);
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
  plotMetm.Draw(c,kTRUE,"pdf",1);
  
  //
  // W- MET plot anti-selection
  //
  RooPlot *awmmframe = pfmet->frame(Bins(NBINS)); 
  awmmframe->GetYaxis()->SetNdivisions(505);
  awmmframe->GetXaxis()->SetLabelOffset(2.0);
  antiMetm.plotOn(awmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetm.plotOn(awmmframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetm.plotOn(awmmframe,LineColor(linecolorW));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),LineColor(linecolorEWK));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(*(aqcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(*(aqcdm.model))),LineColor(linecolorQCD));
  apdfMetm.plotOn(awmmframe,Components(RooArgSet(apdfWmm)),LineColor(linecolorW),LineStyle(2));
  antiMetm.plotOn(awmmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetm->GetBinWidth(1));
  CPlot plotAntiMetm("fitantimetm",awmmframe,"","",ylabel);
  plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrow#mu^{-}#bar{#nu}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetm.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetm.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotAntiMetm.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetm.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetmDiff("fitantimetm","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetmDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetmDiff->GetYaxis()->SetLabelSize(0.11);
  plotAntiMetmDiff.SetYRange(-0.2,0.2);
  plotAntiMetmDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotAntiMetmDiff.AddLine(0, 0.10,METMAX, 0.10,kBlack,3);
  plotAntiMetmDiff.AddLine(0,-0.10,METMAX,-0.10,kBlack,3);
//   plotMetmDiff.AddGraph(gMetmDiff,"Err","FPE2",kGreen,1,kGreen);
  plotAntiMetmDiff.AddHist1D(hAntiMetmDiff,"EX0",ratioColor);
  plotAntiMetmDiff.Draw(c,kTRUE,format,2);
  plotAntiMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotAntiMetm.SetName("fitantimetmlog");
  plotAntiMetm.SetLogy();
  plotAntiMetm.SetYRange(1e-3*(hAntiDataMetm->GetMaximum()),10*(hAntiDataMetm->GetMaximum()));
  plotAntiMetm.Draw(c,kTRUE,format,1);
  plotAntiMetm.Draw(c,kTRUE,"pdf",1);
    
  double chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  double chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  std::cout << chi2prob << " " << chi2ndf << std::endl;
  chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  std::cout << chi2prob << " " << chi2ndf << std::endl;

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;  
  
  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  
  gBenchmark->Show("postFitWm");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff0 = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    Double_t diff = diff0/hData->GetBinContent(ibin);
    std::cout << "bin # " << ibin << std::endl;
    std::cout << "fit bin content " << hFit->GetBinContent(ibin) << std::endl;
    std::cout << "fit bin err " << hFit->GetBinError(ibin) << std::endl;
    std::cout << "data bin content " << hData->GetBinContent(ibin) << std::endl;
    std::cout << "data bin err " << hData->GetBinError(ibin) << std::endl;
    std::cout << "diff =  " << diff << std::endl;
//     Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    //Double_t err = sqrt(hData->GetBinContent(ibin));
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    std::cout << "err = " << err << std::endl;
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
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}


void apply_postfit_shape(TH1D *histogram, TH1D *up, TH1D *down, double shift, double shift_err)
{
  for(int bin=1; bin<=histogram->GetNbinsX(); bin++)
    {
      double upper = up->GetBinContent(bin);
      double lower = down->GetBinContent(bin);
      double central = histogram->GetBinContent(bin);
      double value = 1;
      std::cout << "bin = " << bin << std::endl;
      std::cout << "central = " << central << std::endl;
      std::cout << "lower = " << lower << std::endl;
      std::cout << "upper = " << upper << std::endl;
      if (down->GetBinContent(bin) < 0.0001)
	lower = 0;
      if(shift > 0 && central > 0 && upper != central)
	value = shift * (upper / central - 1) + 1;
      else if(shift < 0 && central > 0 && lower != central){
	value = abs(shift) * (lower / central - 1) + 1;
    std::cout << "shift = " << abs(shift) << std::endl;
    std::cout << "low/c = " << lower / central << std::endl;
    std::cout << "low/c -1 = " << lower / central -1 << std::endl;
      }
      if(value != 1)
	{
      std::cout << "value = " << value << std::endl;
	  central *= value;
      std::cout << "new central = " << central << std::endl;
	  histogram->SetBinContent(bin, central);
	}
      double uncertainty = 1;
      if (lower && central)
	{
	  uncertainty = TMath::Max(uncertainty,central/lower);
	  uncertainty = TMath::Max(uncertainty,lower/central);
	}
      if (upper && central)
	{
	  uncertainty = TMath::Max(uncertainty,upper / central);
	  uncertainty = TMath::Max(uncertainty,central / upper);
	}
      double temp = uncertainty-1;
      uncertainty = shift_err * TMath::Min(2.0,temp);
      double error = sqrt(histogram->GetBinError(bin)*histogram->GetBinError(bin)+ (histogram->GetBinContent(bin) * uncertainty)*(histogram->GetBinContent(bin) * uncertainty));
      histogram->SetBinError(bin, error);
      std::cout << "histo bin error " << histogram->GetBinError(bin) << std::endl;
      std::cout << "error =" << error << std::endl;
    }
}
void apply_postfit_shape(TH1D *histogram, TH1D *up1, TH1D *down1,TH1D *up2, TH1D *down2, double shift1, double shift2)
{
  for(int bin=1; bin<=histogram->GetNbinsX(); bin++)
    {
      double upper1 = up1->GetBinContent(bin);
      double lower1 = down1->GetBinContent(bin);
      double upper2 = up2->GetBinContent(bin);
      double lower2 = down2->GetBinContent(bin);
      double central = histogram->GetBinContent(bin);
      double value1 = 1, value2=1;
      if (down1->GetBinContent(bin) < 0.00001)
	lower1 = 0;
      if (down2->GetBinContent(bin) < 0.00001)
	lower2 = 0;
      if(shift1 > 0 && central > 0 && upper1 != central)
	value1 = shift1 * (upper1 / central - 1) + 1;
      else if(shift1 < 0 && central > 0 && lower1 != central)
	value1 = abs(shift1) * (lower1 / central - 1) + 1;
      if(shift2 > 0 && central > 0 && upper2 != central)
	value2 = shift2 * (upper2 / central - 1) + 1;
      else if(shift2 < 0 && central > 0 && lower2 != central)
	value2 = abs(shift2) * (lower2 / central - 1) + 1;
      if(value1 != 1 && value2 != 1)
	{
	  //central *= value1;
	  central = central*value1+central*value2-central;
	  histogram->SetBinContent(bin, central);
	}
      /*  double uncertainty = 1;
      if (lower && central)
	{
	  uncertainty = TMath::Max(uncertainty,central/lower);
	  uncertainty = TMath::Max(uncertainty,lower/central);
	}
      if (upper && central)
	{
	  uncertainty = TMath::Max(uncertainty,upper / central);
	  uncertainty = TMath::Max(uncertainty,central / upper);
	}
      double temp = uncertainty-1;
      uncertainty = shift_err * TMath::Min(2.0,temp);
      double error = sqrt(histogram->GetBinError(bin)*histogram->GetBinError(bin)+ (histogram->GetBinContent(bin) * uncertainty)*(histogram->GetBinContent(bin) * uncertainty));
      histogram->SetBinError(bin, error);*/
    }
}


TGraphAsymmErrors* makeShapeErrorBand(TH1 *hUp, TH1 *hDown){
  if (!hUp || !hDown) cout << "TH1TOTGraph: histogram not found !" << endl;
  TGraphAsymmErrors* g1= new TGraphAsymmErrors();
  
//   if (hUp->GetN()!=g2->GetN())
//   cout << " graphs have not the same # of elements " <<g1->GetN()<<" "<<g2->GetN()<< endl;
    
  Double_t x, y, exh,exl, eyh,eyl;
  for (Int_t i=0; i<hUp->GetNbinsX(); i++) {
    y=0;
    if(hDown->GetBinContent(i+1) < 0){
     eyl=fabs(hDown->GetBinContent(i+1));
     eyh=fabs(hUp->GetBinContent(i+1));
    } else {
     eyl=fabs(hUp->GetBinContent(i+1));
     eyh=fabs(hDown->GetBinContent(i+1));
    }
    x=hUp->GetBinCenter(i+1);
    exl=hUp->GetBinWidth(i+1)/2;
    exh=hUp->GetBinWidth(i+1)/2;
    
//     std::cout << "low = "

    g1->SetPoint(i,x,y);
    g1->SetPointError(i,exl,exh,eyl,eyh);

  }
  return g1;
  
}