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
  TString inputName = "./Wmunu_reweight2/Wmunu_pdfTemplates.root";
//   TString inputName = "./Wmunu_etaBins/Wmunu_pdfTemplates_3.root";
  TFile *inWmunuShapes = new TFile(inputName); assert(inWmunuShapes);
  RooWorkspace* combine_workspace = (RooWorkspace*) inWmunuShapes->Get("combine_workspace");
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
  RooAbsPdf *wmet  = combine_workspace->pdf("wm");
  RooAbsPdf *wmetp = combine_workspace->pdf("wmp");
  RooAbsPdf *wmetm = combine_workspace->pdf("wmm");
  // Make Histograms
  TH1D *hWmunuMet  = (TH1D*) wmet->createHistogram("pfmet",NBINS);  hWmunuMet->Sumw2();
  TH1D *hWmunuMetp = (TH1D*) wmetp->createHistogram("pfmet",NBINS);hWmunuMetp->Sumw2();
  TH1D *hWmunuMetm = (TH1D*) wmetm->createHistogram("pfmet",NBINS);hWmunuMetm->Sumw2();
  
  TH1D *hWmunuMetp_UncUp   = (TH1D*) hWmunuMetp->Clone("hWmunuMetp_UncUp");
  TH1D *hWmunuMetp_UncDown = (TH1D*) hWmunuMetp->Clone("hWmunuMetp_UncDown");
  
  TH1D *hWmunuMetm_UncUp   = (TH1D*) hWmunuMetm->Clone("hWmunuMetm_UncUp");
  TH1D *hWmunuMetm_UncDown = (TH1D*) hWmunuMetm->Clone("hWmunuMetm_UncDown");
  
  RooAbsData *wmetData = combine_workspace->embeddedData("wmunuMET");
  RooAbsData *wmetpData = combine_workspace->embeddedData("wmunuMETp");
  RooAbsData *wmetmData = combine_workspace->embeddedData("wmunuMETm");
  
  TH1D *hWmunuMetWerr  = (TH1D*) wmetData->createHistogram("pfmet",NBINS);  hWmunuMetWerr->Sumw2();
  TH1D *hWmunuMetpWerr  = (TH1D*) wmetpData->createHistogram("pfmet",NBINS);  hWmunuMetpWerr->Sumw2();
  TH1D *hWmunuMetmWerr  = (TH1D*) wmetmData->createHistogram("pfmet",NBINS);  hWmunuMetmWerr->Sumw2();
  
//   for(int ibin = 1; ibin < hWmunuMet->GetNbinsX(); ++ibin){
//     hWmunuMet->SetBinError(ibin, hWmunuMetWerr->GetBinError(ibin));
//     hWmunuMetp->SetBinError(ibin, hWmunuMetpWerr->GetBinError(ibin));
//     hWmunuMetm->SetBinError(ibin, hWmunuMetmWerr->GetBinError(ibin));
//   }
  
  // Electroweak background contributions
  RooAbsPdf *ewk  = combine_workspace->pdf("ewk");
  RooAbsPdf *ewkp = combine_workspace->pdf("ewkp");
  RooAbsPdf *ewkm = combine_workspace->pdf("ewkm");
  // Make histograms
  TH1D *hEWKMet  = (TH1D*) ewk->createHistogram("pfmet",NBINS);  hEWKMet->Sumw2();
  TH1D *hEWKMetp = (TH1D*) ewkp->createHistogram("pfmet",NBINS);hEWKMetp->Sumw2();
  TH1D *hEWKMetm = (TH1D*) ewkm->createHistogram("pfmet",NBINS);hEWKMetm->Sumw2();
  
  // Shape variation, recoil uncertainty - upwards
  RooAbsPdf *wmet_u  = combine_workspace->pdf("wm_RecoilUp");
  RooAbsPdf *wmetp_u = combine_workspace->pdf("wmp_RecoilUp");
  RooAbsPdf *wmetm_u = combine_workspace->pdf("wmm_RecoilUp");
  // Make histograms
  TH1D *hWmunuMet_RecoilUp  = (TH1D*) wmet_u->createHistogram("pfmet",NBINS);   hWmunuMet_RecoilUp->Sumw2();
  TH1D *hWmunuMetp_RecoilUp = (TH1D*) wmetp_u->createHistogram("pfmet",NBINS); hWmunuMetp_RecoilUp->Sumw2();
  TH1D *hWmunuMetm_RecoilUp = (TH1D*) wmetm_u->createHistogram("pfmet",NBINS); hWmunuMetm_RecoilUp->Sumw2();
  
  // Shape variation, recoil uncertainty - downwards
  RooAbsPdf *wmet_d  = combine_workspace->pdf("wm_RecoilDown");
  RooAbsPdf *wmetp_d = combine_workspace->pdf("wmp_RecoilDown");
  RooAbsPdf *wmetm_d = combine_workspace->pdf("wmm_RecoilDown");
  // Make histograms
  TH1D *hWmunuMet_RecoilDown  = (TH1D*) wmet_d->createHistogram("pfmet",NBINS); hWmunuMet_RecoilDown->Sumw2();
  TH1D *hWmunuMetp_RecoilDown = (TH1D*) wmetp_d->createHistogram("pfmet",NBINS); hWmunuMetp_RecoilDown->Sumw2();
  TH1D *hWmunuMetm_RecoilDown = (TH1D*) wmetm_d->createHistogram("pfmet",NBINS); hWmunuMetm_RecoilDown->Sumw2();
  
  
  // Input the  parameters for signal and background yields
  // These are hardcoded since we get them from previous steps
  RooRealVar nSig("nSig","nSig",1000);
  RooRealVar nQCD("nQCD","nQCD",1000);
  RooRealVar nEWK("nEWK","nEWK",1000);

  // set some values for W+
//   RooRealVar nSigp("nSigp","nSigp",1.0004*9526797.471);
//   RooRealVar nQCDp("nQCDp","nQCDp",10);
//   RooRealVar nEWKp("nEWKp","nEWKp",1.0004*871008.6457);
//   //   
//   RooRealVar nSigm("nSigm","nSigm",1.0000*7370960.276);
//   RooRealVar nQCDm("nQCDm","nQCDm",10);
//   RooRealVar nEWKm("nEWKm","nEWKm",1.0000*781686.9035); 

// pf met
//   RooRealVar nSigp("nSigp","nSigp",0.9905*9202459.033);
//   RooRealVar nQCDp("nQCDp","nQCDp",10);
//   RooRealVar nEWKp("nEWKp","nEWKp",0.9905*841315.7742);
//   //   
//   RooRealVar nSigm("nSigm","nSigm",0.9022*7917395.403);
//   RooRealVar nQCDm("nQCDm","nQCDm",10);
//   RooRealVar nEWKm("nEWKm","nEWKm",0.9022*839599.0043); 
  
//     RooRealVar nSigp("nSigp","nSigp",1.0042*850815.7225);
//   RooRealVar nQCDp("nQCDp","nQCDp",10);
//   RooRealVar nEWKp("nEWKp","nEWKp",1.0042*80458.9007);
//   //   
//   RooRealVar nSigm("nSigm","nSigm",0.9999*665063.0581);
//   RooRealVar nQCDm("nQCDm","nQCDm",10);
//   RooRealVar nEWKm("nEWKm","nEWKm",0.9999*71977.76818); 

  RooRealVar nSigp("nSigp","nSigp",0.9999*9326252.746);
  RooRealVar nQCDp("nQCDp","nQCDp",10);
  RooRealVar nEWKp("nEWKp","nEWKp",0.9999*866619.6083);
  //   
  RooRealVar nSigm("nSigm","nSigm",1.0021*7145364.342);
  RooRealVar nQCDm("nQCDm","nQCDm",10);
  RooRealVar nEWKm("nEWKm","nEWKm",1.0021*770341.466); 

  
//   return;
  
  // Signal PDFs
  
//   double shape_p = 0.2143;
//   double shape_m = -0.0078;
//   double shape_p_err = 0.0102;
//   double shape_m_err =0.0110;

  double shape_p = 0.0220;
  double shape_m = -0.1634;
  double shape_p_err = 0.0210;
  double shape_m_err = 0.0199;

  std::cout << "Nominal Shape W+" << std::endl;
  apply_postfit_shape(hWmunuMetp, hWmunuMetp_RecoilUp, hWmunuMetp_RecoilDown, shape_p, shape_p_err);
  std::cout << "Nominal Shape W-" << std::endl;
  apply_postfit_shape(hWmunuMetm, hWmunuMetm_RecoilUp, hWmunuMetm_RecoilDown, shape_m, shape_m_err);
  
  std::cout << "hWmunumetp error = " << hWmunuMetp->GetBinError(60) << std::endl;
  
  std::cout << "Up Shape W+" << std::endl;
  apply_postfit_shape(hWmunuMetp_UncUp, hWmunuMetp_RecoilUp, hWmunuMetp_RecoilDown, 1.0, shape_p_err);
  std::cout << "Down Shape W+" << std::endl;
  apply_postfit_shape(hWmunuMetp_UncDown, hWmunuMetp_RecoilUp, hWmunuMetp_RecoilDown, -1.0, shape_p_err);
  
  std::cout << "Up Shape W-" << std::endl;
  apply_postfit_shape(hWmunuMetm_UncUp, hWmunuMetm_RecoilUp, hWmunuMetm_RecoilDown, 1.0, shape_m_err);
  std::cout << "Down Shape W-" << std::endl;
  apply_postfit_shape(hWmunuMetm_UncDown, hWmunuMetm_RecoilUp, hWmunuMetm_RecoilDown, -1.0, shape_m_err);
  
  RooDataHist wmunuMet ("wmunuMET", "wmunuMET", RooArgSet(*pfmet),hWmunuMet);  RooHistPdf pdfWm ("wm", "wm", *pfmet,wmunuMet, 1);
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(*pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",*pfmet,wmunuMetp,1);
  RooDataHist wmunuMetm("wmunuMETm","wmunuMETm",RooArgSet(*pfmet),hWmunuMetm); RooHistPdf pdfWmm("wmm","wmm",*pfmet,wmunuMetm,1); 
  
  RooDataHist wmunuMetp_UncUp("wmunuMETp_UncUp","wmunuMETp_UncUp",RooArgSet(*pfmet),hWmunuMetp_UncUp); RooHistPdf pdfWmp_UncUp("wmp_UncUp","wmp_UncUp",*pfmet,wmunuMetp_UncUp,1);
  RooDataHist wmunuMetm_UncUp("wmunuMETm_UncUp","wmunuMETm_UncUp",RooArgSet(*pfmet),hWmunuMetm_UncUp); RooHistPdf pdfWmm_UncUp("wmm_UncUp","wmm_UncUp",*pfmet,wmunuMetm_UncUp,1); 
  RooDataHist wmunuMetp_UncDown("wmunuMETp_UncDown","wmunuMETp_UncDown",RooArgSet(*pfmet),hWmunuMetp_UncDown); RooHistPdf pdfWmp_UncDown("wmp_UncDown","wmp_UncDown",*pfmet,wmunuMetp_UncDown,1);
  RooDataHist wmunuMetm_UncDown("wmunuMETm_UncDown","wmunuMETm_UncDown",RooArgSet(*pfmet),hWmunuMetm_UncDown); RooHistPdf pdfWmm_UncDown("wmm_UncDown","wmm_UncDown",*pfmet,wmunuMetm_UncDown,1); 
  
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
  CPepeModel1 qcd("qcd",*pfmet);
  CPepeModel1 qcdp("qcdp",*pfmet);
  CPepeModel1 qcdm("qcdm",*pfmet);//,qcdp.sigma);
  
  // Signal + Background PDFs
  RooAddPdf pdfMet ("pdfMet", "pdfMet", RooArgList(pdfWm,pdfEWK,*(qcd.model)),   RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm("pdfMetm","pdfMetm",RooArgList(pdfWmm,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));
  
  RooAddPdf pdfMetp_UncUp("pdfMetp_UncUp","pdfMetp_UncUp",RooArgList(pdfWmp_UncUp,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm_UncUp("pdfMetm_UncUp","pdfMetm_UncUp",RooArgList(pdfWmm_UncUp,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));
  
  RooAddPdf pdfMetp_UncDown("pdfMetp_UncDown","pdfMetp_UncDown",RooArgList(pdfWmp_UncDown,pdfEWKp,*(qcdp.model)),RooArgList(nSigp,nEWKp,nQCDp));
  RooAddPdf pdfMetm_UncDown("pdfMetm_UncDown","pdfMetm_UncDown",RooArgList(pdfWmm_UncDown,pdfEWKm,*(qcdm.model)),RooArgList(nSigm,nEWKm,nQCDm));

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(*pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(*pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(*pfmet),hDataMetm);

  std::cout  << "Selectedp: " << hDataMetp->Integral() << endl;
  std::cout  << "Selectedm: " << hDataMetm->Integral() << endl;

// // Fixed Tail on QCD
//   qcdp.a1->setVal(0.1421);
//   qcdm.a1->setVal(0.1423);
//   qcdp.sigma->setVal(15.2281);
//   qcdm.sigma->setVal(15.2638);
//   nQCDp.setVal(1063991.4498);
//   nQCDm.setVal(1052542.2145);
  
  qcdp.a1->setVal(0.2116);
  qcdm.a1->setVal(0.2112);
  qcdp.sigma->setVal(13.8811);
  qcdm.sigma->setVal(14.0460);
  nQCDp.setVal(1274213.1519);
  nQCDm.setVal(1272907.1591);



  
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
  
  TH1D *hPdfMetp_UncUp = (TH1D*)(pdfMetp_UncUp.createHistogram("hPdfMetp_UncUp", *pfmet));
  hPdfMetp_UncUp->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp_UncUp->Integral());
  TH1D *hMetpDiff_UncUp = makeDiffHist(hDataMetp,hPdfMetp_UncUp,"hMetpDiff_UncUp");
  hMetpDiff_UncUp->Add(hMetpDiff,-1);
//   hMetpDiff_UncUp->SetLineColor(kBlue);
//   hMetpDiff_UncUp->SetMarkerSize(0.5);
  
  TH1D *hPdfMetp_UncDown = (TH1D*)(pdfMetp_UncDown.createHistogram("hPdfMetp_UncDown", *pfmet));
  hPdfMetp_UncDown->Scale((nSigp.getVal()+nEWKp.getVal()+nQCDp.getVal())/hPdfMetp_UncDown->Integral());
  TH1D *hMetpDiff_UncDown = makeDiffHist(hDataMetp,hPdfMetp_UncDown,"hMetpDiff_UncDown");
  hMetpDiff_UncDown->Add(hMetpDiff,-1);
//   hMetpDiff_UncDown->SetLineColor(kRed);
//   hMetpDiff_UncDown->SetMarkerSize(0.5);

  TGraphAsymmErrors* gMetpDiff=makeShapeErrorBand(hMetpDiff_UncUp, hMetpDiff_UncDown);
    
  TH1D *hPdfMetm = (TH1D*)(pdfMetm.createHistogram("hPdfMetm", *pfmet));
  for(int ibin = 1; ibin < hPdfMetm->GetNbinsX(); ++ibin){hPdfMetm->SetBinError(ibin, hWmunuMetm->GetBinError(ibin));}
  hPdfMetm->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm->Integral());
  TH1D *hMetmDiff = makeDiffHist(hDataMetm,hPdfMetm,"hMetmDiff");
  hMetmDiff->SetMarkerStyle(kFullCircle); 
  hMetmDiff->SetMarkerSize(0.9);
  
  TH1D *hPdfMetm_UncUp = (TH1D*)(pdfMetm_UncUp.createHistogram("hPdfMetm_UncUp", *pfmet));
  hPdfMetm_UncUp->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm_UncUp->Integral());
  TH1D *hMetmDiff_UncUp = makeDiffHist(hDataMetm,hPdfMetm_UncUp,"hMetmDiff_UncUp");
  hMetmDiff_UncUp->Add(hMetmDiff,-1);
  
  TH1D *hPdfMetm_UncDown = (TH1D*)(pdfMetm_UncDown.createHistogram("hPdfMetm_UncDown", *pfmet));
  hPdfMetm_UncDown->Scale((nSigm.getVal()+nEWKm.getVal()+nQCDm.getVal())/hPdfMetm_UncDown->Integral());
  TH1D *hMetmDiff_UncDown = makeDiffHist(hDataMetm,hPdfMetm_UncDown,"hMetmDiff_UncDown");
  hMetmDiff_UncDown->Add(hMetmDiff,-1);
//   hMetmDiff_UncDown->Draw();
//   return;
  
  TGraphAsymmErrors* gMetmDiff=makeShapeErrorBand(hMetmDiff_UncUp, hMetmDiff_UncDown);
  
 
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
//   plotMetmDiff.AddHist1D(hMetmDiff_UncUp,"F",kBlue);
//   plotMetmDiff.AddHist1D(hMetmDiff_UncDown,"F",kRed);
  plotMetmDiff.Draw(c,kTRUE,format,2);
  plotMetmDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetm.SetName("fitmetmlog");
  plotMetm.SetLogy();
  plotMetm.SetYRange(1e-3*(hDataMetm->GetMaximum()),10*(hDataMetm->GetMaximum()));
  plotMetm.Draw(c,kTRUE,format,1);
  plotMetm.Draw(c,kTRUE,"pdf",1);
    
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
      if (down1->GetBinContent(bin) < 0.01)
	lower1 = 0;
      if (down2->GetBinContent(bin) < 0.01)
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