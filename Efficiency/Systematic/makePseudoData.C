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

#include "EffData.hh"
#include "CEffUser1D.hh"            // class for handling efficiency graphs
#include "CEffUser2D.hh"            // class for handling efficiency tables
#include "ZSignals.hh"
#include "ZBackgrounds.hh"

#include "BinInfo.hh"

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

void generateHistTemplates(const TString infilename, const Float_t ptLo, const Float_t ptHi, const Float_t etaLo, const Float_t etaHi,
			   const Float_t phiLo, const Float_t phiHi, const Float_t npvLo, const Float_t npvHi, const Int_t absEta,
			   const Int_t charge, const Int_t iBin);

void makePseudoData(const TString inputDir="CB_MuSelEff/analysis/plots/", const TString binName ="etapt_0", const Int_t nPsExp=50, const TString outputDir="/scratch/klawhorn/EffSysStore/")
{

  // RETRIEVE TRUTH PDFS

  TFile *f = new TFile(inputDir+binName+".root");
  RooWorkspace *w = (RooWorkspace*) f->Get("w");

  RooCategory sample("sample", "");
  sample.defineType("Pass", 1);
  sample.defineType("Fail", 2);

  RooRealVar m("m","mass",60, 120);
  m.setBins(10000);

  //RooArgSet allVariables = w->allVars();



  RooAbsPdf *modelFailGen = w->pdf("modelFail");
  RooAbsPdf *modelPassGen = w->pdf("modelPass");

  RooSimultaneous totalPdfGen("totalPdfGen","totalPdfGen",sample);
  totalPdfGen.addPdf(*modelPassGen,"Pass");  
  totalPdfGen.addPdf(*modelFailGen,"Fail");

  TTree *intree = (TTree*)f->Get("Bin");
  BinInfo bin;
  intree->SetBranchAddress("Bin",&bin);
  intree->GetEntry(0);

  /*  TString histTempName="/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probes.root";

  generateHistTemplates(histTempName, bin.ptLo, bin.ptHi, bin.etaLo, bin.etaHi, bin.phiLo, bin.phiHi, bin.npvLo, bin.npvHi, bin.absEta, 0, bin.iBin);
    
  TFile *histfile = 0;
  histfile = new TFile("histTemplates.root");
  assert(histfile);

  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;

  char hname[50];
  sprintf(hname,"pass_%i",bin.iBin);
  TH1D *h = (TH1D*)histfile->Get(hname);
  assert(h);
  sigPass = new CMCTemplateConvGaussian(m,h,kTRUE);

  bkgPass = new CExponential(m,kTRUE);

  char hname2[50];
  sprintf(hname2,"fail_%i",bin.iBin);
  h = (TH1D*)histfile->Get(hname2);
  assert(h);
  sigFail = new CMCTemplateConvGaussian(m,h,kFALSE);

  bkgFail = new CExponential(m,kFALSE);

  Double_t NsigMax = bin.nEvents;
  Double_t NbkgFailMax = bin.nBkgFail;
  Double_t NbkgPassMax = bin.nBkgPass;

  RooRealVar Nsignal("Nsignal","Signal Yield",0.80*NsigMax,0,NsigMax);
  RooRealVar effCalc("effCalc","Efficiency",0.8,0,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",50,0,NbkgPassMax);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",100,0.01,NbkgFailMax);

  RooFormulaVar NsigPass("NsigPass","effCalc*Nsignal",RooArgList(effCalc,Nsignal));
  RooFormulaVar NsigFail("NsigFail","(1.0-effCalc)*Nsignal",RooArgList(effCalc,Nsignal));

  RooAddPdf *modelPass, *modelFail;
  modelPass = new RooAddPdf("modelPass", "Model for fitting pass sample", RooArgList(*(sigPass->model), *(bkgPass->model)), RooArgList(NsigPass,NbkgPass));
  modelFail = new RooAddPdf("modelFail", "Model for fitting fail sample", RooArgList(*(sigFail->model), *(bkgFail->model)), RooArgList(NsigFail,NbkgFail));

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");
  totalPdf.addPdf(*modelFail,"Fail");
  */  
  // SET UP AND RUN PSEUDO EXPERIMENTS

  RooMCStudy mgr(totalPdfGen, totalPdfGen, RooArgList(m, sample), "", "2r");

  mgr.generate(nPsExp, bin.nEvents, kFALSE, outputDir+binName+"_%04d.dat");

  //mgr.generateAndFit(1, bin.nEvents);

  //mgr.plotParam(effCalc);

  //for (Int_t i=0; i<nPsExp; i++) {
    // mgr.fitResult(i)->Print("v");
  //}
  
}

void generateHistTemplates(const TString infilename, const Float_t ptLo, const Float_t ptHi, const Float_t etaLo, const Float_t etaHi,
			   const Float_t phiLo, const Float_t phiHi, const Float_t npvLo, const Float_t npvHi, const Int_t absEta,
			   const Int_t charge, const Int_t iBin)
{

  cout << "Creating histogram templates..."; cout.flush();

  Int_t fitMassLo = 60;
  Int_t fitMassHi = 120;

  TH1D* passEtaPt;  
  TH1D* failEtaPt;

  char hname[50];
  sprintf(hname,"pass_%i",iBin);
  passEtaPt = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
  passEtaPt->SetDirectory(0);
  sprintf(hname,"fail_%i",iBin);
  failEtaPt = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
  failEtaPt->SetDirectory(0);

  TFile infile(infilename);
  TTree *intree = (TTree*)infile.Get("Events");
  EffData data;
  intree->SetBranchAddress("Events",&data);  

  for(UInt_t ientry=0;ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    if((data.q)*charge < 0) continue;

    Int_t ipt=-1;
    if (ptLo == 999) 
      ipt = 1;
    else if( (data.pt >= ptLo) && (data.pt < ptHi) )
      ipt = 1;
    if(ipt<0) continue;

    Int_t ieta=-1;
    if (etaLo == 999)
      ieta = 1;
    else if(absEta == 1) {
      assert(etaLo>=0);
      if((fabs(data.eta) >= etaLo) && (fabs(data.eta) < etaHi) )
	ieta = 1; 
    } else {
      if( (data.eta >= etaLo) && (data.eta < etaHi) )
	ieta = 1; 
    }
    if(ieta<0) continue;

    Int_t iphi=-1;
    if (phiLo == 999)
      iphi = 1;
    else if( (data.phi >= phiLo) && (data.phi < phiHi) )
      iphi = 1; 
    if(iphi<0) continue;

    Int_t inpv=-1;
    if (npvLo == 999)
      iphi = 1;
    else if( (data.npv >= npvLo) && (data.npv < npvHi) )
      inpv = 1; 
    if(inpv<0) continue;

    if(data.pass) passEtaPt->Fill(data.mass,1);
    else failEtaPt->Fill(data.mass,1);
  }
  infile.Close();

  TFile outfile("histTemplates.root", "RECREATE");
  passEtaPt->Write();
  failEtaPt->Write();
  delete passEtaPt;
  delete failEtaPt;

  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}
