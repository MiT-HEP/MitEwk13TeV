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

#define BIN_SIZE_PASS 2
#define BIN_SIZE_FAIL 2

void generateHistTemplates(const TString infilename, const Double_t ptLo, const Double_t ptHi, const Double_t etaHi, const Double_t etaLo,
			   const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const Int_t count);

void runJob(Int_t count, Double_t massFail, Double_t massPass, Double_t widthFail, Double_t widthPass,
	    Double_t NbkgFail, Double_t NbkgPass, Double_t Nsig, Double_t alphaFail, Double_t alphaPass,
	    Double_t eff, Double_t meanFail, Double_t meanPass, Double_t nFail, Double_t nPass, Double_t sigmaFail,
	    Double_t sigmaPass, Double_t tFail, Double_t tPass, Double_t passing, Double_t failing, Double_t etaLo, Double_t etaHi,
	    Double_t ptLo, Double_t ptHi, Bool_t doAbsEta) 
{

  // DILEPTON MASS -- VARIABLE

  const Double_t fitMassLo=60;
  const Double_t fitMassHi=120;

  RooRealVar m("m","mass",fitMassLo,fitMassHi);
  m.setBins(10000);

  // SETUP GENERATOR / TRUTH PDFS

  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);

  CSignalModel     *sigPassGen = 0;
  CBackgroundModel *bkgPassGen = 0;
  CSignalModel     *sigFailGen = 0;
  CBackgroundModel *bkgFailGen = 0;  

  sigPassGen = new CBreitWignerConvCrystalBall(m,kTRUE,massPass,widthPass,meanPass,sigmaPass,alphaPass,nPass);
  bkgPassGen = new CExponential(m,kTRUE,tPass);
  sigFailGen = new CBreitWignerConvCrystalBall(m,kTRUE,massFail,widthFail,meanFail,sigmaFail,alphaFail,nFail);
  bkgFailGen = new CExponential(m,kTRUE,tFail);

  Double_t sampleSize = passing+failing;
  Double_t failSize = failing;
  Double_t passSize = passing;

  RooRealVar NsigGen("NsigGen","Signal Yield",Nsig,0,sampleSize);
  
  RooRealVar effGen("effGen","Efficiency",eff,0,1.0);  

  RooRealVar NbkgPassGen("NbkgPassGen","Background count in PASS sample",NbkgPass,0,failSize);
  RooRealVar NbkgFailGen("NbkgFailGen","Background count in FAIL sample",NbkgFail,1,passSize);

  RooFormulaVar NsigPassGen("NsigPassGen","effGen*NsigGen",RooArgList(effGen,NsigGen));
  RooFormulaVar NsigFailGen("NsigFailGen","(1.0-effGen)*NsigGen",RooArgList(effGen,NsigGen));

  RooAddPdf *modelPassGen, *modelFailGen;
  modelPassGen = new RooAddPdf("modelPassGen", "Model for generating pass sample", RooArgList(*(sigPassGen->model), *(bkgPassGen->model)), RooArgList(NsigPassGen, NbkgPassGen));
  modelFailGen = new RooAddPdf("modelFailGen", "Model for generating fail sample", RooArgList(*(sigFailGen->model), *(bkgFailGen->model)), RooArgList(NsigFailGen, NbkgFailGen));

  RooSimultaneous totalPdfGen("totalPdfGen","totalPdfGen",sample);
  totalPdfGen.addPdf(*modelPassGen,"Pass");
  totalPdfGen.addPdf(*modelFailGen,"Fail");

  // SETUP FIT PDFS
  
  TString histTempName="/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probes.root";

  generateHistTemplates(histTempName, ptLo, ptHi, etaHi, etaLo, fitMassLo, fitMassHi, doAbsEta, 0, count);
  
  TFile *histfile = 0;
  histfile = new TFile("histTemplates.root");
  assert(histfile);

  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;

  char hname[50];
  sprintf(hname,"pass%s_%i","etapt",count);
  TH1D *h = (TH1D*)histfile->Get(hname);
  assert(h);
  sigPass = new CMCTemplateConvGaussian(m,h,kTRUE);

  bkgPass = new CExponential(m,kTRUE);

  char hname2[50];
  sprintf(hname2,"fail%s_%i","etapt",count);
  h = (TH1D*)histfile->Get(hname2);
  assert(h);
  sigFail = new CMCTemplateConvGaussian(m,h,kFALSE);

  bkgFail = new CExponential(m,kFALSE);

  Double_t NsigMax = passing+failing;
  Double_t NbkgFailMax = failing;
  Double_t NbkgPassMax = passing;

  RooRealVar Nsignal("Nsignal","Signal Yield",0.80*NsigMax,0,NsigMax);
  RooRealVar effCalc("effCalc","Efficiency",0.8,0,1.0);
  RooRealVar NbackgdPass("NbackgdPass","Background count in PASS sample",50,0,NbkgPassMax);
  RooRealVar NbackgdFail("NbackgdFail","Background count in FAIL sample",0.1*NbkgFailMax,0.01,NbkgFailMax);  

  RooFormulaVar NsigPass("NsigPass","effCalc*Nsignal",RooArgList(effCalc,Nsignal));
  RooFormulaVar NsigFail("NsigFail","(1.0-effCalc)*Nsignal",RooArgList(effCalc,Nsignal));

  RooAddPdf *modelPass, *modelFail;
  modelPass = new RooAddPdf("modelPass", "Model for fitting pass sample", RooArgList(*(sigPass->model), *(bkgPass->model)), RooArgList(NsigPass,NbackgdPass));
  modelFail = new RooAddPdf("modelFail", "Model for fitting fail sample", RooArgList(*(sigFail->model), *(bkgFail->model)), RooArgList(NsigFail,NbackgdFail));

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");
  totalPdf.addPdf(*modelFail,"Fail");

  // SET UP AND RUN PSEUDO EXPERIMENTS

  RooMCStudy mgr(totalPdfGen, totalPdf, RooArgList(m, sample, NsigPass, NsigFail), "", "0r");

  mgr.generateAndFit(10, passing+failing);

  mgr.plotParam(effCalc);

  for (Int_t i=0; i<10; i++) {
    mgr.fitResult(i)->Print("v");
  }
  
}

void generateHistTemplates(const TString infilename, const Double_t ptLo, const Double_t ptHi, const Double_t etaHi, const Double_t etaLo,
			   const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const Int_t count)
{

  cout << "Creating histogram templates..."; cout.flush();

  TH1D* passEtaPt;  
  TH1D* failEtaPt;

  char hname[50];
  sprintf(hname,"passetapt_%i",count);
  passEtaPt = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
  passEtaPt->SetDirectory(0);
  sprintf(hname,"failetapt_%i",count);
  failEtaPt = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
  failEtaPt->SetDirectory(0);

  TFile infile(infilename);
  TTree *intree = (TTree*)infile.Get("Events");
  EffData data;
  intree->SetBranchAddress("Events",&data);  

  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    if((data.q)*charge < 0) continue;

    Int_t ipt=-1;
    if((data.pt >= ptLo) && (data.pt < ptHi))
      ipt = 1;
    if(ipt<0) continue;

    Int_t ieta=-1;
    if(doAbsEta) {
      assert(etaLo>=0);
        if((fabs(data.eta) >= etaLo) && (fabs(data.eta) < etaHi))
          ieta = 1;
    } else {
      if((data.eta >= etaLo) && (data.eta < etaHi))
	ieta = 1;
    }
    if(ieta<0) continue;

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
