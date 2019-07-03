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

void plotResiduals(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
  gBenchmark->Start("postFitWm");
  bool noRecoil = true;

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
  
  TH1D *errorBandRec = new TH1D("errorBandRec","errorBandRec",NBINS,0,METMAX);
  TH1D *errorBandMet = new TH1D("errorBandMet","errorBandMet",NBINS,0,METMAX);
  
  TH1D *envRec = new TH1D("envRec","envRec",NBINS,0,METMAX);
  TH1D *envMet = new TH1D("envMet","envMet",NBINS,0,METMAX);
  
	for(int i =0 ; i < NBINS ; ++i){
		errorBandRec->SetBinContent(i,0);
		errorBandMet->SetBinContent(i,0);
        envRec->SetBinContent(i,0);
		envMet->SetBinContent(i,0);
        errorBandRec->SetBinError(i,0);
		errorBandMet->SetBinError(i,0);
        envRec->SetBinError(i,0);
		envMet->SetBinError(i,0);
	}
    
    
  TString lumiSec = "";
  
  // MAIN FILE
  // Read in the Wmunu_pdfTemplates file here and start to get pdfs.
  // TString inputName = "./Zmm_plot_full/Zmumu_pdfTemplates.root";
  // TString inputName = "./Zmm_plot_eta_new/Zmumu_pdfTemplates.root";
  TString inputName = "./TEST_Zm"+lumiSec+"_Main_2G_ZptReweight/Zmumu_pdfTemplates.root";
  // TString inputName = "./Wmunu_etaBins/Wmunu_pdfTemplates_3.root";
  TFile *inWmunuShapes = new TFile(inputName); assert(inWmunuShapes);
  RooWorkspace* combine_workspace = (RooWorkspace*) inWmunuShapes->Get("combine_workspace");
  combine_workspace->Print();
  // Get the pfmet variable from the workspace
  RooRealVar *pfmet =  combine_workspace->var("pfmet");
  
    // --------------------------------------------------
  // Start getting the PDFs and turning them into histograms ? (wut)
  // Do data first
  RooAbsData *datap = combine_workspace->data("dataMetp");
  TH1D *hDataMetp = (TH1D*) datap->createHistogram("pfmet",NBINS);hDataMetp->Sumw2();
  
  RooAbsPdf *wmetp = combine_workspace->pdf("wmp");
  TH1D *hWmunuMetp = (TH1D*) wmetp->createHistogram("pfmet",NBINS);hWmunuMetp->Sumw2();
  
  RooAbsData *wmetpData = combine_workspace->embeddedData("wmunuMETp");
  TH1D *hWmunuMetpWerr  = (TH1D*) wmetpData->createHistogram("pfmet",NBINS);  hWmunuMetpWerr->Sumw2();
  
  RooAbsPdf *ewkp = combine_workspace->pdf("ewkp");
  TH1D *hEWKMetp = (TH1D*) ewkp->createHistogram("pfmet",NBINS);hEWKMetp->Sumw2();
  
  RooAbsData *datarecoilp; RooAbsPdf *wrecoilp; RooAbsData *wrecoilpData;RooAbsPdf *ewkrecoilp;
  TH1D *hDataRecoilp;TH1D *hWmunuRecoilp; TH1D *hWmunuRecoilpWerr; TH1D *hEWKRecoilp;
  if(noRecoil){
      datarecoilp = combine_workspace->data("dataMetp");
      hDataRecoilp = (TH1D*) datarecoilp->createHistogram("pfmet",NBINS);hDataRecoilp->Sumw2();
  
      wrecoilp = combine_workspace->pdf("wmp");
      hWmunuRecoilp = (TH1D*) wrecoilp->createHistogram("pfmet",NBINS);hWmunuRecoilp->Sumw2();
  
      wrecoilpData = combine_workspace->embeddedData("wmunuMETp");
      hWmunuRecoilpWerr  = (TH1D*) wrecoilpData->createHistogram("pfmet",NBINS);  hWmunuRecoilpWerr->Sumw2();
  
      ewkrecoilp = combine_workspace->pdf("ewkp");
      hEWKRecoilp = (TH1D*) ewkrecoilp->createHistogram("pfmet",NBINS);hEWKRecoilp->Sumw2();
    }else {
      datarecoilp = combine_workspace->data("dataRecoilp");
      hDataRecoilp = (TH1D*) datarecoilp->createHistogram("pfmet",NBINS);hDataRecoilp->Sumw2();
   
      wrecoilp = combine_workspace->pdf("wmprec");
      hWmunuRecoilp = (TH1D*) wrecoilp->createHistogram("pfmet",NBINS);hWmunuRecoilp->Sumw2();
  
      wrecoilpData = combine_workspace->embeddedData("wmunuRecoilp");
      hWmunuRecoilpWerr  = (TH1D*) wrecoilpData->createHistogram("pfmet",NBINS);  hWmunuRecoilpWerr->Sumw2();
  
      ewkrecoilp = combine_workspace->pdf("ewkprec");
      hEWKRecoilp = (TH1D*) ewkrecoilp->createHistogram("pfmet",NBINS);hEWKRecoilp->Sumw2();
    }
  
    // Input the  parameters for signal and background yields
  // These are hardcoded since we get them from previous steps
  RooRealVar nSigp("nSigp","nSigp",hDataRecoilp->Integral()*0.995);
  RooRealVar nEWKp("nEWKp","nEWKp",hDataRecoilp->Integral()*0.005);
  
    // RooRealVar nSigp("nSigp","nSigp",8.0096e+05);
  // RooRealVar nEWKp("nEWKp","nEWKp",8.0096e+05*4.4803e-03);
  
  // list of names: 
  // Int_t nShapes = 3;
  std::vector<TString> fnamev;
    // inputName = "./Zmm_plot_eta_diag_new/Zmumu_pdfTemplates.root";
  // fnamev.push_back(inputName);  
  // inputName = "./Zmm_plot_eta/Zmumu_pdfTemplates.root";
  // fnamev.push_back(inputName);  
  // inputName = "./Zmm_plot_keys/Zmumu_pdfTemplates.root";
  // fnamev.push_back(inputName);  
  // inputName = "./Zmm_plot_diag2/Zmumu_pdfTemplates.root";
  // fnamev.push_back(inputName);  
  inputName = "./TEST_Zm"+lumiSec+"_Keys_2G_ZptReweight/Zmumu_pdfTemplates.root";
  fnamev.push_back(inputName);
  inputName = "./TEST_Zm"+lumiSec+"_Eta_2G_ZptReweight/Zmumu_pdfTemplates.root";
  fnamev.push_back(inputName);
  inputName = "./TEST_Zm"+lumiSec+"m_StatFixed_PlotCheck/Zmumu_pdfTemplates.root";
  fnamev.push_back(inputName);
  
  int nShapes = fnamev.size();
  // loop through the alternate shapes:
  for(int i = 0; i < nShapes; ++i){
	  std::cout << "blah" << std::endl;
	  inputName = (TString)fnamev.at(i);
	  std::cout << inputName.Data() << std::endl;
	  TFile *inWmunuShapes_alt = new TFile(inputName); assert(inWmunuShapes_alt);
	  RooWorkspace* combine_wksp_alt = (RooWorkspace*) inWmunuShapes_alt->Get("combine_workspace");
	  combine_wksp_alt->Print();
	  RooRealVar *pfmet_alt =  combine_wksp_alt->var("pfmet");

	  
	  RooAbsPdf *wmetp_a = combine_wksp_alt->pdf("wmp");
      TH1D *hWmunuMetp_alt = (TH1D*) wmetp_a->createHistogram("pfmet",NBINS);hWmunuMetp_alt->Sumw2();
	  
      RooAbsPdf *wrecoilp_a; TH1D *hWmunuRecoilp_alt;
      if(noRecoil){
	      wrecoilp_a = combine_wksp_alt->pdf("wmp");
          hWmunuRecoilp_alt = (TH1D*) wrecoilp_a->createHistogram("pfmet",NBINS);hWmunuRecoilp_alt->Sumw2();
      } else {
          wrecoilp_a = combine_wksp_alt->pdf("wmprec");
          hWmunuRecoilp_alt = (TH1D*) wrecoilp_a->createHistogram("pfmet",NBINS);hWmunuRecoilp_alt->Sumw2();
      }
	  
	  // loop through the bins, take the relative difference between the main and alternate
	  // add in quadrature to existing bin contents
	  for(int i = 0; i <= NBINS; ++i){
		  
		  double current = (errorBandRec->GetBinError(i))*(errorBandRec->GetBinError(i));
		  double addon = (hWmunuRecoilp_alt->GetBinContent(i)-hWmunuRecoilp->GetBinContent(i))/hWmunuRecoilp->GetBinContent(i);
          if(hWmunuRecoilp->GetBinContent(i)==0) addon=0;
		  std::cout << "addon = " << addon << std::endl;
		  current+=addon*addon;
		  errorBandRec->SetBinError(i,sqrt(current));
		  std::cout << "current Recoil " << errorBandRec->GetBinContent(i) << std::endl;
		  current=0;addon=0;
          

		  
		  current = (errorBandMet->GetBinError(i))*(errorBandMet->GetBinError(i));
		  addon = (hWmunuMetp_alt->GetBinContent(i)-hWmunuMetp->GetBinContent(i))/hWmunuMetp->GetBinContent(i);
          if(hWmunuMetp->GetBinContent(i)==0) addon=0;
		  std::cout << "addon = " << addon << std::endl;
		  current+=addon*addon;
		  errorBandMet->SetBinError(i,sqrt(current));
          current=0;addon=0;
          
          
          current = envRec->GetBinError(i);
		  addon = (hWmunuRecoilp_alt->GetBinContent(i)-hWmunuRecoilp->GetBinContent(i))/hWmunuRecoilp->GetBinContent(i);
          if(hWmunuRecoilp->GetBinContent(i)==0) addon=0;
		  std::cout << "addon = " << addon << std::endl;
		  current+=addon;
		  envRec->SetBinError(i,current);
		  std::cout << "env Recoil " << envRec->GetBinContent(i) << std::endl;
		  current=0;addon=0;
          
          current = envMet->GetBinError(i);
          std::cout << "current " << current << std::endl;
		  addon = (hWmunuMetp_alt->GetBinContent(i)-hWmunuMetp->GetBinContent(i))/hWmunuMetp->GetBinContent(i);
          if(hWmunuMetp->GetBinContent(i)==0) addon=0;
		  std::cout << "addon = " << addon << std::endl;
		  current+=addon;
		  envMet->SetBinError(i,current);
          current=0;addon=0;
		  
		  std::cout << "env MET " << envMet->GetBinContent(i) << std::endl;
		  // for(int ibin = 1; ibin < errorBandMet->GetNbinsX(); ++ibin){errorBandMet->SetBinError(ibin, hWmunuRecoilp->GetBinError(i)*errorBandMet->GetBinContent(i));}
		  
	  }
	  delete inWmunuShapes_alt;
	  delete combine_wksp_alt;
  
  }
  
   TH1D *errorBandNegMet = (TH1D*) errorBandMet->Clone("errNegMet");
   errorBandNegMet->Scale(-1);
   TH1D *errorBandNegRec = (TH1D*) errorBandRec->Clone("errNegRec");
   errorBandNegRec->Scale(-1);
   
   TH1D *envNegMet = (TH1D*) envMet->Clone("envNegMet");
   envNegMet->Scale(-1);
   TH1D *envNegRec = (TH1D*) envRec->Clone("envNegRec");
   // envNegRec->Scale(-1);
  


  
  RooDataHist wmunuMetp("wmunuMETp","wmunuMETp",RooArgSet(*pfmet),hWmunuMetp); RooHistPdf pdfWmp("wmp","wmp",*pfmet,wmunuMetp,1);
  RooDataHist ewkMetp("ewkMETp","ewkMETp",RooArgSet(*pfmet),hEWKMetp); RooHistPdf pdfEWKp("ewkp","ewkp",*pfmet,ewkMetp,1); 
  RooAddPdf pdfMetp("pdfMetp","pdfMetp",RooArgList(pdfWmp,pdfEWKp),RooArgList(nSigp,nEWKp));
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(*pfmet),hDataMetp);
  
  RooDataHist wmunuRecoilp("wmunuRecoilp","wmunuRecoilp",RooArgSet(*pfmet),hWmunuRecoilp); RooHistPdf pdfWmprec("wmprec","wmprec",*pfmet,wmunuRecoilp,1);
  RooDataHist ewkRecoilp("ewkRecoilp","ewkRecoilp",RooArgSet(*pfmet),hEWKRecoilp); RooHistPdf pdfEWKprec("ewkprec","ewkprec",*pfmet,ewkRecoilp,1); 
  RooAddPdf pdfRecoilp("pdfRecoilp","pdfRecoilp",RooArgList(pdfWmprec,pdfEWKprec),RooArgList(nSigp,nEWKp));
  RooDataHist dataRecoilp("dataRecoilp","dataRecoilp",RooArgSet(*pfmet),hDataRecoilp);

  std::cout  << "Selectedp: " << hDataMetp->Integral() << endl;
  std::cout  << "EWKp: " << nEWKp.getVal() << endl;
  
  for(int ibin = 1; ibin < hWmunuMetp->GetNbinsX(); ++ibin){
    hWmunuMetp->SetBinError(ibin, hWmunuMetpWerr->GetBinError(ibin));
    hWmunuRecoilp->SetBinError(ibin, hWmunuRecoilpWerr->GetBinError(ibin));
  }
  TH1D *hPdfMetp = (TH1D*)(pdfMetp.createHistogram("hPdfMetp", *pfmet));
  for(int ibin = 1; ibin < hPdfMetp->GetNbinsX(); ++ibin){hPdfMetp->SetBinError(ibin, hWmunuMetp->GetBinError(ibin));}
  hPdfMetp->Scale((nSigp.getVal()+nEWKp.getVal())/hPdfMetp->Integral());
  TH1D *hMetpDiff = makeDiffHist(hDataMetp,hPdfMetp,"hMetpDiff");
  hMetpDiff->SetMarkerStyle(kFullCircle);
  hMetpDiff->SetMarkerSize(0.9);
  
  TH1D *hPdfRecoilp = (TH1D*)(pdfRecoilp.createHistogram("hPdfRecoilp", *pfmet));
  for(int ibin = 1; ibin < hPdfRecoilp->GetNbinsX(); ++ibin){hPdfRecoilp->SetBinError(ibin, hWmunuRecoilp->GetBinError(ibin));}
  hPdfRecoilp->Scale((nSigp.getVal()+nEWKp.getVal())/hPdfRecoilp->Integral());
  TH1D *hRecoilpDiff = makeDiffHist(hDataRecoilp,hPdfRecoilp,"hRecoilpDiff");
  hRecoilpDiff->SetMarkerStyle(kFullCircle);
  hRecoilpDiff->SetMarkerSize(0.9);
 
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
  
  
  // double ymax = 0.1;
  double ymax = 0.5;
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
  // W+ MET plot
  //
  RooPlot *wmpframe = pfmet->frame(Bins(NBINS));
  wmpframe->GetYaxis()->SetNdivisions(505);
  wmpframe->GetXaxis()->SetLabelOffset(2.0);
  dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMetp.plotOn(wmpframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,LineColor(linecolorW));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfEWKp)),LineColor(linecolorEWK));
  pdfMetp.plotOn(wmpframe,Components(RooArgSet(pdfWmp)),LineColor(linecolorW),LineStyle(2));
  dataMetp.plotOn(wmpframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataMetp->GetBinWidth(1));
  CPlot plotMetp("fitmetp",wmpframe,"","",ylabel);
  plotMetp.SetLegend(0.68,0.57,0.93,0.77);
  plotMetp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMetp.GetLegend()->AddEntry(hDummyW,"Signal","F");
  plotMetp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  // plotMetp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMetp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotMetp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotMetp.Draw(c,kFALSE,format,1);

  CPlot plotMetpDiff("fitmetp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
//   plotMetpDiff.AddGraph(gMetpDiff,"Data","PE2",3,1,3);
  plotMetpDiff.AddHist1D(errorBandMet,"E3",kRed,1,1);
  // plotMetpDiff.AddHist1D(errorBandNegMet,"C",kRed);
  plotMetpDiff.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpDiff.SetYRange(-ymax,ymax);
  plotMetpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpDiff.AddLine(0, ymax*0.5,METMAX, ymax*0.5,kBlack,3);
  plotMetpDiff.AddLine(0,-ymax*0.5,METMAX,-ymax*0.5,kBlack,3);
  plotMetpDiff.Draw(c,kTRUE,format,2);
  plotMetpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("fitmetplog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  plotMetp.Draw(c,kTRUE,"pdf",1);
  
  plotMetp.SetLogy(kFALSE);
  // c->Clear();
    // plotMetp.Draw(c,kFALSE,format,1);
  CPlot plotMetpEnv("fitmetpenv","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hMetpDiff->GetYaxis()->SetTitleOffset(0.5);
  hMetpDiff->GetYaxis()->SetLabelSize(0.11);
//   plotMetpEnv.AddGraph(gMetpDiff,"Data","PE2",3,1,3);
  plotMetpEnv.AddHist1D(envMet,"E3",kRed,1,1);
  // plotMetpEnv.AddHist1D(envNegMet,"C",kRed);
  plotMetpEnv.AddHist1D(hMetpDiff,"EX0",ratioColor);
  plotMetpEnv.SetYRange(-ymax,ymax);
  plotMetpEnv.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetpEnv.AddLine(0, ymax*0.5,METMAX, ymax*0.5,kBlack,3);
  plotMetpEnv.AddLine(0,-ymax*0.5,METMAX,-ymax*0.5,kBlack,3);
  plotMetpEnv.Draw(c,kTRUE,format,2);
  plotMetpEnv.Draw(c,kTRUE,"pdf",2);
  
  plotMetp.SetName("fitmetpenvlog");
  plotMetp.SetLogy();
  plotMetp.SetYRange(1e-3*(hDataMetp->GetMaximum()),10*(hDataMetp->GetMaximum()));
  plotMetp.Draw(c,kTRUE,format,1);
  plotMetp.Draw(c,kTRUE,"pdf",1);
  
  RooPlot *wmpframerec = pfmet->frame(Bins(NBINS));
  wmpframerec->GetYaxis()->SetNdivisions(505);
  wmpframerec->GetXaxis()->SetLabelOffset(2.0);
  dataRecoilp.plotOn(wmpframerec,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfRecoilp.plotOn(wmpframerec,FillColor(fillcolorW),DrawOption("F"));
  pdfRecoilp.plotOn(wmpframerec,LineColor(linecolorW));
  pdfRecoilp.plotOn(wmpframerec,Components(RooArgSet(pdfEWKprec)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfRecoilp.plotOn(wmpframerec,Components(RooArgSet(pdfEWKprec)),LineColor(linecolorEWK));
  pdfRecoilp.plotOn(wmpframerec,Components(RooArgSet(pdfWmprec)),LineColor(linecolorW),LineStyle(2));
  dataRecoilp.plotOn(wmpframerec,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",hDataRecoilp->GetBinWidth(1));
  CPlot plotRecoilp("fitRecoilp",wmpframerec,"","",ylabel);
  plotRecoilp.SetLegend(0.68,0.57,0.93,0.77);
  plotRecoilp.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotRecoilp.GetLegend()->AddEntry(hDummyW,"Signal","F");
  plotRecoilp.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  // plotRecoilp.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotRecoilp.AddTextBox("#bf{CMS} #scale[0.75]{#it{Preliminary}}",0.62,0.80,0.88,0.88,0);
  plotRecoilp.AddTextBox(lumitext,0.66,0.91,0.95,0.96,0);
  plotRecoilp.Draw(c,kFALSE,format,1);

  CPlot plotRecoilpDiff("fitRecoilp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hRecoilpDiff->GetYaxis()->SetTitleOffset(0.5);
  hRecoilpDiff->GetYaxis()->SetLabelSize(0.11);
//   plotRecoilpDiff.AddGraph(gRecoilpDiff,"Data","PE2",3,1,3);

  plotRecoilpDiff.AddHist1D(errorBandRec,"E3",kRed,1,1);
  // plotRecoilpDiff.AddHist1D(errorBandNegRec,"C",kRed);
  plotRecoilpDiff.AddHist1D(hRecoilpDiff,"EX0",ratioColor);
  plotRecoilpDiff.SetYRange(-ymax,ymax);
  plotRecoilpDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotRecoilpDiff.AddLine(0, ymax*0.5,METMAX, ymax*0.5,kBlack,3);
  plotRecoilpDiff.AddLine(0,-ymax*0.5,METMAX,-ymax*0.5,kBlack,3);
  plotRecoilpDiff.Draw(c,kTRUE,format,2);
  plotRecoilpDiff.Draw(c,kTRUE,"pdf",2);
  
  plotRecoilp.SetName("fitRecoilplog");
  plotRecoilp.SetLogy();
  plotRecoilp.SetYRange(1e-3*(hDataRecoilp->GetMaximum()),10*(hDataRecoilp->GetMaximum()));
  plotRecoilp.Draw(c,kTRUE,format,1);
  plotRecoilp.Draw(c,kTRUE,"pdf",1);
  
  
  plotRecoilp.SetLogy(kFALSE);
  CPlot plotRecoilpEnv("fitRecoilenvp","","#slash{E}_{T} [GeV]","#frac{Data-Pred}{Data}");
  hRecoilpDiff->GetYaxis()->SetTitleOffset(0.5);
  hRecoilpDiff->GetYaxis()->SetLabelSize(0.11);
//   plotRecoilpDiff.AddGraph(gRecoilpDiff,"Data","PE2",3,1,3);
  
  plotRecoilpEnv.AddHist1D(envRec,"E3",kRed,1,1);
  // plotRecoilpEnv.AddHist1D(envNegRec,"C",kRed);
  plotRecoilpEnv.AddHist1D(hRecoilpDiff,"EX0",ratioColor);
  plotRecoilpEnv.SetYRange(-ymax,ymax);
  plotRecoilpEnv.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotRecoilpEnv.AddLine(0, ymax*0.5,METMAX, ymax*0.5,kBlack,3);
  plotRecoilpEnv.AddLine(0,-ymax*0.5,METMAX,-ymax*0.5,kBlack,3);
  plotRecoilpEnv.Draw(c,kTRUE,format,2);
  plotRecoilpEnv.Draw(c,kTRUE,"pdf",2);
  
  plotRecoilp.SetName("fitRecoilenvplog");
  plotRecoilp.SetLogy();
  plotRecoilp.SetYRange(1e-3*(hDataRecoilp->GetMaximum()),10*(hDataRecoilp->GetMaximum()));
  plotRecoilp.Draw(c,kTRUE,format,1);
  plotRecoilp.Draw(c,kTRUE,"pdf",1);
  
  
    
  double chi2prob = hDataMetp->Chi2Test(hPdfMetp,"PUW");
  double chi2ndf  = hDataMetp->Chi2Test(hPdfMetp,"CHI2/NDFUW");
  // std::cout << chi2prob << " " << chi2ndf << std::endl;
  // chi2prob = hDataMetm->Chi2Test(hPdfMetm,"PUW");
  // chi2ndf  = hDataMetm->Chi2Test(hPdfMetm,"CHI2/NDFUW");
  // std::cout << chi2prob << " " << chi2ndf << std::endl;

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

    if(hData->GetBinContent(ibin) == 0) diff = 0;
    // std::cout << "data " << hData->GetBinContent(ibin) << std::endl;
    // std::cout << "fits " << hFit->GetBinContent(ibin) << std::endl;
//     Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((1.0/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    Double_t err = (hFit->GetBinContent(ibin)/hData->GetBinContent(ibin))*sqrt((hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))*(hFit->GetBinError(ibin)/hFit->GetBinContent(ibin))+(1.0/hData->GetBinContent(ibin)));
    if(hData->GetBinContent(ibin) == 0) err = 0;
    //if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    //if(err>0) hDiff->SetBinContent(ibin,diff/err);
    //else      hDiff->SetBinContent(ibin,0);
    std::cout << "bin # " << ibin << std::endl;
    std::cout << "fit bin content " << hFit->GetBinContent(ibin) << std::endl;
    std::cout << "fit bin err " << hFit->GetBinError(ibin) << std::endl;
    std::cout << "data bin content " << hData->GetBinContent(ibin) << std::endl;
    std::cout << "data bin err " << hData->GetBinError(ibin) << std::endl;
    std::cout << "diff =  " << diff << std::endl;
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