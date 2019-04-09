#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "CPlot.hh"                 // helper class for plots
#include "MitStyleRemix.hh"         // style settings for drawing
#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h" 
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooMCStudy.h"
#include "RooWorkspace.h"

#define BIN_SIZE_PASS 2
#define BIN_SIZE_FAIL 2

void makePseudoData(const TString sigDir,
		    const TString bkgDir,
		    const TString binName ="eta_0",
		    const TString outputDir=".",
		    const Int_t siglabel=0,
		    const Int_t bkglabel=0,
		    const Int_t nPsExp=1000
){
  gSystem->mkdir(outputDir,kTRUE);
  // RETRIEVE TRUTH PDFS
  TFile *fsig = new TFile(sigDir+binName+".root");
  TFile *fbkg = new TFile(bkgDir+binName+".root");
  RooWorkspace *wsig = (RooWorkspace*) fsig->Get("w");
  // std::cout << "Signal Workspace" << std::endl;
  // wsig->Print();
  RooWorkspace *wbkg = (RooWorkspace*) fbkg->Get("w");
  // std::cout << "BKG Workspace" << std::endl;
  // wbkg->Print();
  // std::cout << "blah" << std::endl;
  RooCategory sample("sample", "");
  sample.defineType("Pass", 1);
  sample.defineType("Fail", 2);

  RooRealVar m("m","mass",60, 120);
  m.setBins(10000);
  Double_t nsigPass = wsig->var("eff")->getVal() * wsig->var("Nsig")->getVal();
  Double_t nsigFail = (1. - wsig->var("eff")->getVal()) * wsig->var("Nsig")->getVal();
  Double_t nbkgPass = wbkg->var("NbkgPass")->getVal();
  Double_t nbkgFail = wbkg->var("NbkgFail")->getVal();

  std::cout << "pass eff " << wsig->var("eff")->getVal() << std::endl;
  RooRealVar NsigPass("NsigPass", "", nsigPass);
  RooRealVar NsigFail("NsigFail", "", nsigFail);
  RooRealVar NbkgPass("NbkgPass", "", nbkgPass);
  RooRealVar NbkgFail("NbkgFail", "", nbkgFail);

/*
  RooRealVar *NsigPass = wsig->var("NsigPass");
  RooRealVar *NsigFail = wsig->var("NsigFail");
  RooRealVar *NbkgPass = wbkg->var("NbkgPass");
  RooRealVar *NbkgFail = wbkg->var("NbkgFail");
*/

  RooAbsPdf *sigPass;
  RooAbsPdf *sigFail;
  if(siglabel==-1){
      // std::cout <<"wat"<<std::endl;
    sigPass = wsig->pdf("signalPass");
    sigFail = wsig->pdf("signalFail");
  }else{
      // std:cout << "looking for signal model" << std::endl;
    char sigpassname[20];
    char sigfailname[20];
    sprintf(sigpassname,"signalPass_%d",siglabel);
    sprintf(sigfailname,"signalFail_%d",siglabel);
    // std::cout << "names " << sigpassname << "  " << sigfailname << std::endl;
    sigPass = wsig->pdf(sigpassname);
    sigFail = wsig->pdf(sigfailname);
  }
  sigPass->Print();
  sigFail->Print();
  // std::cout << "look for bkg" << std::endl;

  char bkgpassname[20];
  char bkgfailname[20];
  sprintf(bkgpassname,"backgroundPass_%d",bkglabel);
  sprintf(bkgfailname,"backgroundFail_%d",bkglabel);
  RooAbsPdf *bkgPass = wbkg->pdf(bkgpassname);
  RooAbsPdf *bkgFail = wbkg->pdf(bkgfailname);
  // bkgPass->Print();
  // bkgFail->Print();
  // std::cout << "lbbbbbbbbb" << std::endl;

  RooAddPdf *modelPass, *modelFail;
  modelPass = new RooAddPdf("modelPass","Model for PASS sample", RooArgList(*sigPass, *bkgPass), RooArgList(NsigPass,NbkgPass));
  modelFail = new RooAddPdf("modelFail","Model for FAIL sample", RooArgList(*sigFail, *bkgFail), RooArgList(NsigFail,NbkgFail));
// std::cout << "aaaaaaaaaaaaaa" << std::endl;
  RooSimultaneous totalPdfGen("totalPdfGen","totalPdfGen",sample);
  totalPdfGen.addPdf(*modelPass,"Pass");  
  totalPdfGen.addPdf(*modelFail,"Fail");
// std::cout << "cccccccccccc" << std::endl;
  TTree *intree = (TTree*)fsig->Get("Bin");

  UInt_t nEvents;
  intree->SetBranchAddress("nEvents",  &nEvents);
  intree->GetEntry(0);
  // std::cout << "gggggggggggg" << std::endl;
  
  // RooRandom::randomGenerator().SetSeed(11111);
  RooMCStudy mgr(totalPdfGen, totalPdfGen, RooArgList(m, sample), "", "2r");
  mgr.generate(nPsExp, nEvents, kFALSE, outputDir+binName+"_%d.dat");
}
