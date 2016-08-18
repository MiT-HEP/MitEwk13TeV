//================================================================================================
//
// Perform fit to extract W->enu signal
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
#include "RooAbsData.h"
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

void postFitWe(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("postFitWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS  = 75;
  const Double_t METMAX = 150;

  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.5;

  TString pufname = "../Tools/pileup_weights_2015B.root";

  // file format for output plots
  const TString format("png"); 


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  // Read in the Wenu_pdfTemplates file here and start to get pdfs.
//   TString inputName = "Wenu_wShapes/Wenu_pdfTemplates.root";
  TString inputName = "Wenu_pileup_simultaneous_fixLumi/Wenu_pdfTemplates.root";
  TFile *inWenuShapes = new TFile(inputName); assert(inWenuShapes);
  RooWorkspace* combine_workspace = (RooWorkspace*) inWenuShapes->Get("combine_workspace");
  std::cout << "what" << std::endl;
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
 
  
  //  Signal MC
  RooAbsPdf *wmet  = combine_workspace->pdf("we");
  RooAbsPdf *wmetp = combine_workspace->pdf("wep");
  RooAbsPdf *wmetm = combine_workspace->pdf("wem");
  // Make Histograms
  TH1D *hWenuMet  = (TH1D*) wmet->createHistogram("pfmet",NBINS);  hWenuMet->Sumw2();
  TH1D *hWenuMetp = (TH1D*) wmetp->createHistogram("pfmet",NBINS);hWenuMetp->Sumw2();
  TH1D *hWenuMetm = (TH1D*) wmetm->createHistogram("pfmet",NBINS);hWenuMetm->Sumw2();
  
  TH1D *hWenuMetp_UncUp   = (TH1D*) hWenuMetp->Clone("hWenuMetp_UncUp");
  TH1D *hWenuMetp_UncDown = (TH1D*) hWenuMetp->Clone("hWenuMetp_UncDown");
  
  TH1D *hWenuMetm_UncUp   = (TH1D*) hWenuMetm->Clone("hWenuMetm_UncUp");
  TH1D *hWenuMetm_UncDown = (TH1D*) hWenuMetm->Clone("hWenuMetm_UncDown");
  
  // Electroweak background contributions
  RooAbsPdf *ewk  = combine_workspace->pdf("ewk");
  RooAbsPdf *ewkp = combine_workspace->pdf("ewkp");
  RooAbsPdf *ewkm = combine_workspace->pdf("ewkm");
  // Make histograms
  TH1D *hEWKMet  = (TH1D*) ewk->createHistogram("pfmet",NBINS);  hEWKMet->Sumw2();
  TH1D *hEWKMetp = (TH1D*) ewkp->createHistogram("pfmet",NBINS);hEWKMetp->Sumw2();
  TH1D *hEWKMetm = (TH1D*) ewkm->createHistogram("pfmet",NBINS);hEWKMetm->Sumw2();
  
  // Shape variation, recoil uncertainty - upwards
  RooAbsPdf *wmet_ru  = combine_workspace->pdf("we_RecoilUp");
  RooAbsPdf *wmetp_ru = combine_workspace->pdf("wep_RecoilUp");
  RooAbsPdf *wmetm_ru = combine_workspace->pdf("wem_RecoilUp");
  // Make histograms
  TH1D *hWenuMet_RecoilUp  = (TH1D*) wmet_ru->createHistogram("pfmet",NBINS);   hWenuMet_RecoilUp->Sumw2();
  TH1D *hWenuMetp_RecoilUp = (TH1D*) wmetp_ru->createHistogram("pfmet",NBINS); hWenuMetp_RecoilUp->Sumw2();
  TH1D *hWenuMetm_RecoilUp = (TH1D*) wmetm_ru->createHistogram("pfmet",NBINS); hWenuMetm_RecoilUp->Sumw2();
  
  // Shape variation, recoil uncertainty - downwards
  RooAbsPdf *wmet_rd  = combine_workspace->pdf("we_RecoilDown");
  RooAbsPdf *wmetp_rd = combine_workspace->pdf("wep_RecoilDown");
  RooAbsPdf *wmetm_rd = combine_workspace->pdf("wem_RecoilDown");
  
  TH1D *hWenuMet_RecoilDown  = (TH1D*) wmet_rd->createHistogram("pfmet",NBINS);   hWenuMet_RecoilDown->Sumw2();
  TH1D *hWenuMetp_RecoilDown = (TH1D*) wmetp_rd->createHistogram("pfmet",NBINS); hWenuMetp_RecoilDown->Sumw2();
  TH1D *hWenuMetm_RecoilDown = (TH1D*) wmetm_rd->createHistogram("pfmet",NBINS); hWenuMetm_RecoilDown->Sumw2();
  
    // Shape variation, Scale uncertainty - upwards
  RooAbsPdf *wmet_su  = combine_workspace->pdf("we_ScaleUp");
  RooAbsPdf *wmetp_su = combine_workspace->pdf("wep_ScaleUp");
  RooAbsPdf *wmetm_su = combine_workspace->pdf("wem_ScaleUp");
  // Make histograms
  TH1D *hWenuMet_ScaleUp  = (TH1D*) wmet_su->createHistogram("pfmet",NBINS);   hWenuMet_ScaleUp->Sumw2();
  TH1D *hWenuMetp_ScaleUp = (TH1D*) wmetp_su->createHistogram("pfmet",NBINS); hWenuMetp_ScaleUp->Sumw2();
  TH1D *hWenuMetm_ScaleUp = (TH1D*) wmetm_su->createHistogram("pfmet",NBINS); hWenuMetm_ScaleUp->Sumw2();
  
  // Shape variation, Scale uncertainty - downwards
  RooAbsPdf *wmet_sd  = combine_workspace->pdf("we_ScaleDown");
  RooAbsPdf *wmetp_sd = combine_workspace->pdf("wep_ScaleDown");
  RooAbsPdf *wmetm_sd = combine_workspace->pdf("wem_ScaleDown");
  
  TH1D *hWenuMet_ScaleDown  = (TH1D*) wmet_sd->createHistogram("pfmet",NBINS);   hWenuMet_ScaleDown->Sumw2();
  TH1D *hWenuMetp_ScaleDown = (TH1D*) wmetp_sd->createHistogram("pfmet",NBINS); hWenuMetp_ScaleDown->Sumw2();
  TH1D *hWenuMetm_ScaleDown = (TH1D*) wmetm_sd->createHistogram("pfmet",NBINS); hWenuMetm_ScaleDown->Sumw2();
  
  
  // Shape variation, pileup uncertainty - upwards
  RooAbsPdf *wmet_pu  = combine_workspace->pdf("we_PileupUp");
  RooAbsPdf *wmetp_pu = combine_workspace->pdf("wep_PileupUp");
  RooAbsPdf *wmetm_pu = combine_workspace->pdf("wem_PileupUp");
  // Make histograms
  TH1D *hWenuMet_PileupUp  = (TH1D*) wmet_pu->createHistogram("pfmet",NBINS);  hWenuMet_PileupUp->Sumw2();
  TH1D *hWenuMetp_PileupUp = (TH1D*) wmetp_pu->createHistogram("pfmet",NBINS); hWenuMetp_PileupUp->Sumw2();
  TH1D *hWenuMetm_PileupUp = (TH1D*) wmetm_pu->createHistogram("pfmet",NBINS); hWenuMetm_PileupUp->Sumw2();
  
  // Shape variation, pileup uncertainty - downwards
  RooAbsPdf *wmet_pd  = combine_workspace->pdf("we_PileupDown");
  RooAbsPdf *wmetp_pd = combine_workspace->pdf("wep_PileupDown");
  RooAbsPdf *wmetm_pd = combine_workspace->pdf("wem_PileupDown");
  
  TH1D *hWenuMet_PileupDown  = (TH1D*) wmet_pd->createHistogram("pfmet",NBINS);   hWenuMet_PileupDown->Sumw2();
  TH1D *hWenuMetp_PileupDown = (TH1D*) wmetp_pd->createHistogram("pfmet",NBINS); hWenuMetp_PileupDown->Sumw2();
  TH1D *hWenuMetm_PileupDown = (TH1D*) wmetm_pd->createHistogram("pfmet",NBINS); hWenuMetm_PileupDown->Sumw2();
  
    // pileup uncertainty on the ewk background
  RooAbsPdf *ewk_pu  = combine_workspace->pdf("ewk_PileupUp");
  RooAbsPdf *ewkp_pu = combine_workspace->pdf("ewkp_PileupUp");
  RooAbsPdf *ewkm_pu = combine_workspace->pdf("ewkm_PileupUp");
  // Make histograms
  TH1D *hEWKMet_PileupUp  = (TH1D*) ewk_pu->createHistogram("pfmet",NBINS);  hEWKMet_PileupUp->Sumw2();
  TH1D *hEWKMetp_PileupUp = (TH1D*) ewkp_pu->createHistogram("pfmet",NBINS); hEWKMetp_PileupUp->Sumw2();
  TH1D *hEWKMetm_PileupUp = (TH1D*) ewkm_pu->createHistogram("pfmet",NBINS); hEWKMetm_PileupUp->Sumw2();
  
  // Shape variation, pileup uncertainty - downwards
  RooAbsPdf *ewk_pd  = combine_workspace->pdf("ewk_PileupDown");
  RooAbsPdf *ewkp_pd = combine_workspace->pdf("ewkp_PileupDown");
  RooAbsPdf *ewkm_pd = combine_workspace->pdf("ewkm_PileupDown");
  
  TH1D *hEWKMet_PileupDown  = (TH1D*) ewk_pd->createHistogram("pfmet",NBINS);  hEWKMet_PileupDown->Sumw2();
  TH1D *hEWKMetp_PileupDown = (TH1D*) ewkp_pd->createHistogram("pfmet",NBINS); hEWKMetp_PileupDown->Sumw2();
  TH1D *hEWKMetm_PileupDown = (TH1D*) ewkm_pd->createHistogram("pfmet",NBINS); hEWKMetm_PileupDown->Sumw2();
  
  
  RooAbsData *wmetData = combine_workspace->embeddedData("wenuMET");
  RooAbsData *wmetpData = combine_workspace->embeddedData("wenuMETp");
  RooAbsData *wmetmData = combine_workspace->embeddedData("wenuMETm");
  // Make histograms
  TH1D *hWenuMetWerr  = (TH1D*) wmetData->createHistogram("pfmet",NBINS);  hWenuMetWerr->Sumw2();
  TH1D *hWenuMetpWerr  = (TH1D*) wmetpData->createHistogram("pfmet",NBINS);  hWenuMetpWerr->Sumw2();
  TH1D *hWenuMetmWerr  = (TH1D*) wmetmData->createHistogram("pfmet",NBINS);  hWenuMetmWerr->Sumw2();
  
  // ----------------------------------------------------------------------------------------------
  // pick up some of the anti-isolation selection for the background fits
  RooAbsData *anti  = combine_workspace->data("antiMet");
  RooAbsData *antip = combine_workspace->data("antiMetp");
  RooAbsData *antim = combine_workspace->data("antiMetm");
  //  Make Histograms
  TH1D *hAntiDataMet  = (TH1D*) anti->createHistogram("pfmet",NBINS);  hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetp = (TH1D*) antip->createHistogram("pfmet",NBINS); hAntiDataMetp->Sumw2();
  TH1D *hAntiDataMetm = (TH1D*) antim->createHistogram("pfmet",NBINS); hAntiDataMetm->Sumw2();
  
  //  Signal MC
  RooAbsPdf *awmet  = combine_workspace->pdf("awe");
  RooAbsPdf *awmetp = combine_workspace->pdf("awep");
  RooAbsPdf *awmetm = combine_workspace->pdf("awem");
  // Make Histograms
  TH1D *hAntiWenuMet  = (TH1D*) awmet->createHistogram("pfmet",NBINS);  hAntiWenuMet->Sumw2();
  TH1D *hAntiWenuMetp = (TH1D*) awmetp->createHistogram("pfmet",NBINS); hAntiWenuMetp->Sumw2();
  TH1D *hAntiWenuMetm = (TH1D*) awmetm->createHistogram("pfmet",NBINS); hAntiWenuMetm->Sumw2();
  
  // Electroweak background contributions
  RooAbsPdf *aewk  = combine_workspace->pdf("aewk");
  RooAbsPdf *aewkp = combine_workspace->pdf("aewkp");
  RooAbsPdf *aewkm = combine_workspace->pdf("aewkm");
  // Make histograms
  TH1D *hAntiEWKMet  = (TH1D*) aewk->createHistogram("pfmet",NBINS);  hAntiEWKMet->Sumw2();
  TH1D *hAntiEWKMetp = (TH1D*) aewkp->createHistogram("pfmet",NBINS); hAntiEWKMetp->Sumw2();
  TH1D *hAntiEWKMetm = (TH1D*) aewkm->createHistogram("pfmet",NBINS); hAntiEWKMetm->Sumw2();
  
  RooAbsData *wmetAntiData = combine_workspace->embeddedData("awenuMET");
  RooAbsData *wmetpAntiData = combine_workspace->embeddedData("awenuMETp");
  RooAbsData *wmetmAntiData = combine_workspace->embeddedData("awenuMETm");
  // Make histograms
  TH1D *hAntiWenuMetWerr  = (TH1D*) wmetAntiData->createHistogram("pfmet",NBINS);    hAntiWenuMetWerr->Sumw2();
  TH1D *hAntiWenuMetpWerr  = (TH1D*) wmetpAntiData->createHistogram("pfmet",NBINS);  hAntiWenuMetpWerr->Sumw2();
  TH1D *hAntiWenuMetmWerr  = (TH1D*) wmetmAntiData->createHistogram("pfmet",NBINS);  hAntiWenuMetmWerr->Sumw2();
  
  // Input the  parameters for signal and background yields
  // These are hardcoded since we get them from previous steps
  RooRealVar nSig("nSig","nSig",1000);
  RooRealVar nQCD("nQCD","nQCD",1000);
  RooRealVar nEWK("nEWK","nEWK",1000);
  
  RooRealVar nAntiSig("nSig","nSig",1000);
  RooRealVar nAntiQCD("nQCD","nQCD",1000);
  RooRealVar nAntiEWK("nEWK","nEWK",1000);
  
// // //   // === The combine results for simultaneous a1 fit === 
//   RooRealVar nSigp("nSigp","nSigp",0.9977*6524675.797);
//   RooRealVar nQCDp("nQCDp","nQCDp",10);
//   RooRealVar nEWKp("nEWKp","nEWKp",0.9977*736755.8992);
//   //   
//   RooRealVar nSigm("nSigm","nSigm",0.9960*5020548.647);
//   RooRealVar nQCDm("nQCDm","nQCDm",10);
//   RooRealVar nEWKm("nEWKm","nEWKm",0.9960*669877.1517); 
//   
//   // Anti
//   RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",1.0125*287843.1745);
//   RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",10);
//   RooRealVar nAntiEWKp("nAntiEWKp","nAntiEWKp",1.0125*42818.06689);
//   //   
//   RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",1.0178*231372.4399);
//   RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",10);
//   RooRealVar nAntiEWKm("nAntiEWKm","nAntiEWKm",1.0178*43406.11125); 
  
  // //   // === The combine results for new simultaneous a1 fit === 
  RooRealVar nSigp("nSigp","nSigp",0.9977*6524675.797);
  RooRealVar nQCDp("nQCDp","nQCDp",10);
  RooRealVar nEWKp("nEWKp","nEWKp",0.9977*736755.8992);
  //   
  RooRealVar nSigm("nSigm","nSigm",0.9965*5020548.647);
  RooRealVar nQCDm("nQCDm","nQCDm",10);
  RooRealVar nEWKm("nEWKm","nEWKm",0.9965*669877.1517); 
  
  // Anti
  RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",1.0123*287843.1745);
  RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",10);
  RooRealVar nAntiEWKp("nAntiEWKp","nAntiEWKp",1.0123*42818.06689);
  //   
  RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",1.0293*231372.4399);
  RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",10);
  RooRealVar nAntiEWKm("nAntiEWKm","nAntiEWKm",1.0293*43406.11125); 

//   // === The combine results for free a1 fit === 
//   RooRealVar nSigp("nSigp","nSigp",0.9963*6524675.797);
//   RooRealVar nQCDp("nQCDp","nQCDp",10);
//   RooRealVar nEWKp("nEWKp","nEWKp",0.9963*736755.8992);
//   //   
//   RooRealVar nSigm("nSigm","nSigm",0.9888*5020548.647);
//   RooRealVar nQCDm("nQCDm","nQCDm",10);
//   RooRealVar nEWKm("nEWKm","nEWKm",0.9888*669877.1517); 
//   
//   // Anti
//   RooRealVar nAntiSigp("nAntiSigp","nAntiSigp",0.9963*287843.1745);
//   RooRealVar nAntiQCDp("nAntiQCDp","nAntiQCDp",10);
//   RooRealVar nAntiEWKp("nAntiEWKp","nAntiEWKp",0.9963*42818.06689);
//   //   
//   RooRealVar nAntiSigm("nAntiSigm","nAntiSigm",0.9888*231372.4399);
//   RooRealVar nAntiQCDm("nAntiQCDm","nAntiQCDm",10);
//   RooRealVar nAntiEWKm("nAntiEWKm","nAntiEWKm",0.9888*43406.11125); 
  
//   return;
  
  // Signal PDFs
  // Only using the Recoil shape for now, no Scale variations
  
// //   // For the simultaneous fit
//   double nRecoilWep = 0.1448;
//   double nRecoilWem = 0.1605;
// 
//   double nPileupWep = 0.6824;
//   double nPileupWem = 0.6347;

//   // For the simultaneous fit
  double nRecoilWep = 0.1449;
  double nRecoilWem = 0.1451;

  double nPileupWep = 0.6823;
  double nPileupWem = -0.6972;

//   // For the free fit
//   double nRecoilWep = 0.0144;
//   double nRecoilWem = 0.0396;
// 
//   double nPileupWep = 0.5639;
//   double nPileupWem = 0.5762;
  
  
  apply_postfit_shape(hWenuMetp, hWenuMetp_RecoilUp, hWenuMetp_RecoilDown,hWenuMetp_PileupUp, hWenuMetp_PileupDown, nRecoilWep, nPileupWep);
  apply_postfit_shape(hWenuMetm, hWenuMetm_RecoilUp, hWenuMetm_RecoilDown,hWenuMetm_PileupUp, hWenuMetm_PileupDown, nRecoilWem, nPileupWem);
  
  apply_postfit_shape(hEWKMetp, hEWKMetp_PileupUp, hEWKMetp_PileupDown, nPileupWep, 0);
  apply_postfit_shape(hEWKMetm, hEWKMetm_PileupUp, hEWKMetm_PileupDown, nPileupWem, 0);

  
  RooDataHist wenuMet ("wenuMET", "wenuMET", RooArgSet(*pfmet),hWenuMet);  RooHistPdf pdfWe ("we", "we", *pfmet,wenuMet, 1);
  RooDataHist wenuMetp("wenuMETp","wenuMETp",RooArgSet(*pfmet),hWenuMetp); RooHistPdf pdfWep("wep","wep",*pfmet,wenuMetp,1);
  RooDataHist wenuMetm("wenuMETm","wenuMETm",RooArgSet(*pfmet),hWenuMetm); RooHistPdf pdfWem("wem","wem",*pfmet,wenuMetm,1); 
  
  RooDataHist wenuMetp_UncUp("wenuMETp_UncUp","wenuMETp_UncUp",RooArgSet(*pfmet),hWenuMetp_UncUp); RooHistPdf pdfWep_UncUp("wep_UncUp","wep_UncUp",*pfmet,wenuMetp_UncUp,1);
  RooDataHist wenuMetm_UncUp("wenuMETm_UncUp","wenuMETm_UncUp",RooArgSet(*pfmet),hWenuMetm_UncUp); RooHistPdf pdfWem_UncUp("wem_UncUp","wem_UncUp",*pfmet,wenuMetm_UncUp,1); 
  RooDataHist wenuMetp_UncDown("wenuMETp_UncDown","wenuMETp_UncDown",RooArgSet(*pfmet),hWenuMetp_UncDown); RooHistPdf pdfWep_UncDown("wep_UncDown","wep_UncDown",*pfmet,wenuMetp_UncDown,1);
  RooDataHist wenuMetm_UncDown("wenuMETm_UncDown","wenuMETm_UncDown",RooArgSet(*pfmet),hWenuMetm_UncDown); RooHistPdf pdfWem_UncDown("wem_UncDown","wem_UncDown",*pfmet,wenuMetm_UncDown,1); 
  
  RooDataHist wenuMet_RecoilUp("wenuMET_RecoilUp", "wenuMET_RecoilUp", RooArgSet(*pfmet),hWenuMet_RecoilUp);  
  RooHistPdf pdfWe_RecoilUp("we_RecoilUp", "we_RecoilUp", *pfmet,wenuMet_RecoilUp, 1);
  RooDataHist wenuMetp_RecoilUp("wenuMETp_RecoilUp","wenuMETp_RecoilUp",RooArgSet(*pfmet),hWenuMetp_RecoilUp); 
  RooHistPdf pdfWep_RecoilUp("wep_RecoilUp","wep_RecoilUp",*pfmet,wenuMetp_RecoilUp,1);
  RooDataHist wenuMetm_RecoilUp("wenuMETm_RecoilUp","wenuMETm_RecoilUp",RooArgSet(*pfmet),hWenuMetm_RecoilUp); 
  RooHistPdf pdfWem_RecoilUp("wem_RecoilUp","wem_RecoilUp",*pfmet,wenuMetm_RecoilUp,1); 
  RooDataHist wenuMet_RecoilDown("wenuMET_RecoilDown", "wenuMET_RecoilDown", RooArgSet(*pfmet),hWenuMet_RecoilDown);  
  RooHistPdf pdfWe_RecoilDown("we_RecoilDown", "we_RecoilDown", *pfmet,wenuMet_RecoilDown, 1);
  RooDataHist wenuMetp_RecoilDown("wenuMETp_RecoilDown","wenuMETp_RecoilDown",RooArgSet(*pfmet),hWenuMetp_RecoilDown); 
  RooHistPdf pdfWep_RecoilDown("wep_RecoilDown","wep_RecoilDown",*pfmet,wenuMetp_RecoilDown,1);
  RooDataHist wenuMetm_RecoilDown("wenuMETm_RecoilDown","wenuMETm_RecoilDown",RooArgSet(*pfmet),hWenuMetm_RecoilDown); 
  RooHistPdf pdfWem_RecoilDown("wem_RecoilDown","wem_RecoilDown",*pfmet,wenuMetm_RecoilDown,1); 
  
  // EWK+top PDFs
  RooDataHist ewkMet ("ewkMET", "ewkMET", RooArgSet(*pfmet),hEWKMet);  RooHistPdf pdfEWK ("ewk", "ewk", *pfmet,ewkMet, 1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(*pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",*pfmet,ewkMetp,1); 
  RooDataHist ewkMetm("ewkMETm","ewkMETm",RooArgSet(*pfmet),hEWKMetm); RooHistPdf pdfEWKm("ewkm","ewkm",*pfmet,ewkMetm,1); 
  
  // QCD Pdfs
  CPepeModel2 qcd("qcd",*pfmet);
  CPepeModel2 qcdp("qcdp",*pfmet);
  CPepeModel2 qcdm("qcdm",*pfmet);//,qcdp.sigma);
  
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWe,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWep,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWem,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));
  
  RooAddPdf pdfMetp_UncUp("pdfMetp_UncUp","pdfMetp_UncUp",RooArgList(pdfWep_UncUp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm_UncUp("pdfMetm_UncUp","pdfMetm_UncUp",RooArgList(pdfWem_UncUp,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));
  
  RooAddPdf pdfMetp_UncDown("pdfMetp_UncDown","pdfMetp_UncDown",RooArgList(pdfWep_UncDown,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm_UncDown("pdfMetm_UncDown","pdfMetm_UncDown",RooArgList(pdfWem_UncDown,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(*pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(*pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(*pfmet),hDataMetm);

  std::cout  << "Selectedp: " << hDataMetp->Integral() << endl;
  std::cout  << "Selectedm: " << hDataMetm->Integral() << endl;
  
  
  // ========== Anti-selection stuff  ================
  RooDataHist awenuMet ("awenuMET", "awenuMET", RooArgSet(*pfmet),hAntiWenuMet);  RooHistPdf apdfWe ("awe", "awe", *pfmet,awenuMet, 1);
  RooDataHist awenuMetp("awenuMETp","awenuMETp",RooArgSet(*pfmet),hAntiWenuMetp); RooHistPdf apdfWep("awep","awep",*pfmet,awenuMetp,1);
  RooDataHist awenuMetm("awenuMETm","awenuMETm",RooArgSet(*pfmet),hAntiWenuMetm); RooHistPdf apdfWem("awem","awem",*pfmet,awenuMetm,1); 
  
  // EWK+top PDFs
  RooDataHist aewkMet ("aewkMET", "aewkMET", RooArgSet(*pfmet),hAntiEWKMet);  RooHistPdf apdfEWK ("aewk", "aewk", *pfmet,aewkMet, 1);
  RooDataHist aewkMetp("aewkMETp","aewkMETp",RooArgSet(*pfmet),hAntiEWKMetp); RooHistPdf apdfEWKp("aewkp","aewkp",*pfmet,aewkMetp,1); 
  RooDataHist aewkMetm("aewkMETm","aewkMETm",RooArgSet(*pfmet),hAntiEWKMetm); RooHistPdf apdfEWKm("aewkm","aewkm",*pfmet,aewkMetm,1); 
  
  // QCD Pdfs
  CPepeModel2 aqcd("aqcd",*pfmet);
  CPepeModel2 aqcdp("aqcdp",*pfmet);
  CPepeModel2 aqcdm("aqcdm",*pfmet);//,qcdp.sigma);
  
  // Signal + Background PDFs
  RooAddPdf apdfMet ("apdfMet", "apdfMet", RooArgList(apdfWe,apdfEWK,*(aqcd.model)),   RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  RooAddPdf apdfMetp("apdfMetp","apdfMetp",RooArgList(apdfWep,apdfEWKp,*(aqcdp.model)),RooArgList(nAntiSigp,nAntiEWKp,nAntiQCDp));
  RooAddPdf apdfMetm("apdfMetm","apdfMetm",RooArgList(apdfWem,apdfEWKm,*(aqcdm.model)),RooArgList(nAntiSigm,nAntiEWKm,nAntiQCDm));
  
  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(*pfmet),hAntiDataMet);
  RooDataHist antiMetp("antiMetp","antiMetp",RooArgSet(*pfmet),hAntiDataMetp);
  RooDataHist antiMetm("antiMetm","antiMetm",RooArgSet(*pfmet),hAntiDataMetm);
  
// // // Simultaneous a1 background fit
//   qcdp.a1->setVal(0.0125);
//   qcdp.a2->setVal(6.4466);
//   qcdp.a3->setVal(156.2003);
//   qcdm.a1->setVal(0.0144);
//   qcdm.a2->setVal(6.5348);
//   qcdm.a3->setVal(157.5949);
//   nQCDp.setVal(3998995.2872);
//   nQCDm.setVal(4082089.9516);
//   
//   aqcdp.a1->setVal(0.0125);
//   aqcdp.a2->setVal(7.4424);
//   aqcdp.a3->setVal(181.1628);
//   aqcdm.a1->setVal(0.0144);
//   aqcdm.a2->setVal(7.1129);
//   aqcdm.a3->setVal(187.7958);
//   nAntiQCDp.setVal(3800916.7512);
//   nAntiQCDm.setVal(3635825.0707);

// // // Simultaneous a1 background fit
  qcdp.a1->setVal(0.0125);
  qcdp.a2->setVal(6.4468);
  qcdp.a3->setVal(156.1972);
  qcdm.a1->setVal(0.0154);
  qcdm.a2->setVal(6.4399);
  qcdm.a3->setVal(159.9687);
  nQCDp.setVal(3999006.3775);
  nQCDm.setVal(4079359.2314);
  
  aqcdp.a1->setVal(0.0125);
  aqcdp.a2->setVal(7.4428);
  aqcdp.a3->setVal(181.1576);
  aqcdm.a1->setVal(0.0154);
  aqcdm.a2->setVal(7.0094);
  aqcdm.a3->setVal(190.0784);
  nAntiQCDp.setVal(3800981.0354);
  nAntiQCDm.setVal(3632642.4333);


// // // // free fit
//   qcdp.a1->setVal(-0.0260);
//   qcdp.a2->setVal(9.0455);
//   qcdp.a3->setVal(112.5038);
//   qcdm.a1->setVal(-0.0120);
//   qcdm.a2->setVal(8.5756);
//   qcdm.a3->setVal(121.4293);
//   nQCDp.setVal(4003412.9099);
//   nQCDm.setVal(4125308.9488);
//   
//   aqcdp.a1->setVal(0.0125);
//   aqcdp.a2->setVal(7.4426);
//   aqcdp.a3->setVal(181.1589);
//   aqcdm.a1->setVal(0.0144);
//   aqcdm.a2->setVal(7.1130);
//   aqcdm.a3->setVal(187.7957);
//   nAntiQCDp.setVal(3800916.6120);
//   nAntiQCDm.setVal(3635839.0451);

  
  std::cout  << "EWKp: " << nEWKp.getVal() << endl;
  std::cout  << "EWKm: " << nEWKm.getVal() << endl;


  
  for(int ibin = 1; ibin < hWenuMet->GetNbinsX(); ++ibin){
//     std::cout << "bin error " << hWmunuMetp->GetBinError(ibin) << std::endl;
    hWenuMet->SetBinError(ibin, hWenuMetWerr->GetBinError(ibin));
    hWenuMetp->SetBinError(ibin, hWenuMetpWerr->GetBinError(ibin));
    hWenuMetm->SetBinError(ibin, hWenuMetmWerr->GetBinError(ibin));
    hAntiWenuMet->SetBinError(ibin, hAntiWenuMetWerr->GetBinError(ibin));
    hAntiWenuMetp->SetBinError(ibin, hAntiWenuMetpWerr->GetBinError(ibin));
    hAntiWenuMetm->SetBinError(ibin, hAntiWenuMetmWerr->GetBinError(ibin));
  }
  
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram("hPdfMet", *pfmet));
  for(int ibin = 1; ibin < hPdfMet->GetNbinsX(); ++ibin){hPdfMet->SetBinError(ibin, hWenuMet->GetBinError(ibin));}
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(hDataMet,hPdfMet,"hMetDiff");
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
   
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", *pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWenuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
  
  TH1D *hPdfMetp_UncUp = (TH1D*)(pdfMetp_UncUp.createHistogram("hPdfMetp_UncUp", *pfmet));
  hPdfMetp_UncUp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp_UncUp->Integral());
  TH1D *hMetpDiff_UncUp = makeDiffHist(hDataMetp,hPdfMetp_UncUp,"hMetpDiff_UncUp");
  hMetpDiff_UncUp->SetLineColor(kBlue);
//   hMetpDiff_UncUp->SetMarkerSize(0.9);
  
  TH1D *hPdfMetp_UncDown = (TH1D*)(pdfMetp_UncDown.createHistogram("hPdfMetp_UncDown", *pfmet));
  hPdfMetp_UncDown->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp_UncDown->Integral());
  TH1D *hMetpDiff_UncDown = makeDiffHist(hDataMetp,hPdfMetp_UncDown,"hMetpDiff_UncDown");
  hMetpDiff_UncDown->SetLineColor(kRed);
//   hMetpDiff_UncDown->SetMarkerSize(0.9);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", *pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWenuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
  
  // The Anti-selecitons
  TH1D *hAntiPdfMet = (TH1D*)(apdfMet.createHistogram("hAntiPdfMet", *pfmet));
  for(int ibin = 1; ibin < hAntiPdfMet->GetNbinsX(); ++ibin){hAntiPdfMet->SetBinError(ibin, hAntiWenuMet->GetBinError(ibin));}
  hAntiPdfMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hAntiPdfMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(hAntiDataMet,hAntiPdfMet,"hAntiMetDiff");
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
   
  TH1D *hAntiPdfMetp = (TH1D*)(apdfMetp.createHistogram("hAntiPdfMetp", *pfmet));
  for(int ibin = 1; ibin < hAntiPdfMetp->GetNbinsX(); ++ibin){hAntiPdfMetp->SetBinError(ibin, hAntiWenuMetp->GetBinError(ibin));}
  hAntiPdfMetp->Scale((nAntiSigp.getVal()+nAntiEWKp.getVal()+nAntiQCDp.getVal())/hAntiPdfMetp->Integral());
  TH1D *hAntiMetpDiff = makeDiffHist(hAntiDataMetp,hAntiPdfMetp,"hAntiMetpDiff");
  hAntiMetpDiff->SetMarkerStyle(kFullCircle);
  hAntiMetpDiff->SetMarkerSize(0.9);
    
  TH1D *hAntiPdfMetm = (TH1D*)(apdfMetm.createHistogram("hAntiPdfMetm", *pfmet));
  for(int ibin = 1; ibin < hAntiPdfMetm->GetNbinsX(); ++ibin){hAntiPdfMetm->SetBinError(ibin, hAntiWenuMetm->GetBinError(ibin));}
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
  RooPlot *weframe = pfmet->frame(Bins(NBINS));
  weframe->GetYaxis()->SetNdivisions(505);
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(weframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(weframe,LineColor(linecolorW));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfEWK,*(qcd.model))),LineColor(linecolorEWK));
  pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMet.plotOn(weframe,Components(RooArgSet(*(qcd.model))),LineColor(linecolorQCD));
  pdfMet.plotOn(weframe,Components(RooArgSet(pdfWe)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(weframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMet->GetBinWidth(1));
  CPlot plotMet("fitmet",weframe,"","mT [GeV]",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrowe#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  //plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  //plotMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plotMet.AddTextBox("CMS Preliminary",0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox(lumitext,0.63,0.92,0.95,0.99,0);
//  plotMet.SetYRange(0.1,1.1*(hDataMet->GetMaximum()));
//plotMet.SetYRange(0.1,1e4);
  plotMet.Draw(c,kTRUE,format,1);

  CPlot plotMetDiff("fitmet","","#slash{E}_{T} [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-0.1,0.1);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 0.05,METMAX, 0.05,kBlack,3);
  plotMetDiff.AddLine(0,-0.05,METMAX,-0.05,kBlack,3);
  plotMetDiff.Draw(c,kTRUE,format,2);
  
  plotMet.SetName("fitmetlog");
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(hDataMet->GetMaximum()),10*(hDataMet->GetMaximum()));
  plotMet.Draw(c,kTRUE,format,1);
    
  //
  // W+ MET plot
  //
  RooPlot *wepframe = pfmet->frame(Bins(NBINS));    
  wepframe->GetYaxis()->SetNdivisions(505);
  wepframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wepframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wepframe,LineColor(linecolorW));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfEWKp,*(qcdp.model))),LineColor(linecolorEWK));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(*(qcdp.model))),LineColor(linecolorQCD));
  pdfMetp.plotOn(wepframe,Components(RooArgSet(pdfWep)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("fitmetp",wepframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrowe^{+}#nu","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-0.20,0.20);
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
  
  // Anti-selection W+ background fits
  RooPlot *awepframe = pfmet->frame(Bins(NBINS));    
  awepframe->GetYaxis()->SetNdivisions(505);
  awepframe->GetXaxis()->SetLabelOffset(2.0);
  antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetp.plotOn(awepframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetp.plotOn(awepframe,LineColor(linecolorW));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfEWKp,*(aqcdp.model))),LineColor(linecolorEWK));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(*(aqcdp.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(*(aqcdp.model))),LineColor(linecolorQCD));
  apdfMetp.plotOn(awepframe,Components(RooArgSet(apdfWep)),LineColor(linecolorW),LineStyle(2));
  antiMetp.plotOn(awepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetp->GetBinWidth(1));
  CPlot plotAntiMetp("fitantimetp",awepframe,"","",ylabel);
  plotAntiMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetp.GetLegend()->AddEntry(hDummyW,"W^{+}#rightarrowe^{+}#nu","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotAntiMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotAntiMetp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotAntiMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotAntiMetp.Draw(c,kFALSE,format,1);

  CPlot plotAntiMetpDiff("fitantimetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hAntiMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hAntiMetpDiff->GetYaxis()->SetLabelSize(0.11);
  plotAntiMetpDiff.AddHist1D(hAntiMetpDiff,"EX0",ratioColor);
  plotAntiMetpDiff.SetYRange(-0.20,0.20);
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
  RooPlot *wemframe = pfmet->frame(Bins(NBINS)); 
  wemframe->GetYaxis()->SetNdivisions(505);
  wemframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetm.plotOn(wemframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetm.plotOn(wemframe,LineColor(linecolorW));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfEWKm,*(qcdm.model))),LineColor(linecolorEWK));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(*(qcdm.model))),LineColor(linecolorQCD));
  pdfMetm.plotOn(wemframe,Components(RooArgSet(pdfWem)),LineColor(linecolorW),LineStyle(2));
  dataMetm.plotOn(wemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetm->GetBinWidth(1));
  CPlot plotMetm("fitmetm",wemframe,"","",ylabel);
  plotMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrowe^{-}#bar{#nu}","F");
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
  plotMetmDiff.AddHist1D(hMetmDiff,"EX0",ratioColor);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
  plotMetm.Draw(c,kTRUE,"pdf",1);
  
  // Anti-selection W- background fits
  // fix these
  RooPlot *awemframe = pfmet->frame(Bins(NBINS)); 
  awemframe->GetYaxis()->SetNdivisions(505);
  awemframe->GetXaxis()->SetLabelOffset(2.0);
  antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  apdfMetm.plotOn(awemframe,FillColor(fillcolorW),DrawOption("F"));
  apdfMetm.plotOn(awemframe,LineColor(linecolorW));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),FillColor(fillcolorEWK),DrawOption("F"));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfEWKm,*(aqcdm.model))),LineColor(linecolorEWK));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(*(aqcdm.model))),FillColor(fillcolorQCD),DrawOption("F"));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(*(aqcdm.model))),LineColor(linecolorQCD));
  apdfMetm.plotOn(awemframe,Components(RooArgSet(apdfWem)),LineColor(linecolorW),LineStyle(2));
  antiMetm.plotOn(awemframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  
  sprintf(ylabel,"Events / %.1f GeV",hAntiDataMetm->GetBinWidth(1));
  CPlot plotAntiMetm("fitantimetm",awemframe,"","",ylabel);
  plotAntiMetm.SetLegend(0.68,0.57,0.93,0.77);
  plotAntiMetm.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotAntiMetm.GetLegend()->AddEntry(hDummyW,"W^{-}#rightarrowe^{-}#bar{#nu}","F");
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
  
  gBenchmark->Show("postFitWe");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
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
  sprintf(htmlfname,"%s/WenuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wenu</title></head>" << endl;
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
      if (down1->GetBinContent(bin) < 0.0001)
	lower1 = 0;
      if (down2->GetBinContent(bin) < 0.0001)
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
