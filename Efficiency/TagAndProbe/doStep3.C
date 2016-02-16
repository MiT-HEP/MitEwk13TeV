//================================================================================================
//
// Signal Extraction
//-------------------
//  0: probe counting
//  1: Breit-Wigner convolved with Crystal Ball function
//  2: MC template convolved with Gaussian
//  3: Phil's Crystal Ball based "Voigtian" shape
//  4: Unbinned MC data convolved with Gaussian
//
// Background Model
//------------------
//  0: no background
//  1: exponential model
//  2: erfc*exp model
//  3: double exponential model
//  4: linear*exp model
//  5: quadratic*exp model
//
//________________________________________________________________________________________________

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
#endif

// RooFit headers
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"

// bin size constants
#define BIN_SIZE_PASS 1 
#define BIN_SIZE_FAIL 1

//=== MAIN MACRO ================================================================================================= 

void doStep3(const TString infilename,      // input file
		  const TString binname,         // name
		  const TString binfile,
		  const TString resultfolder,
		  const TString resultfilename)         // file with bin info
{
  gBenchmark->Start("doStep3");
  ofstream myfile;

  TString resultfile = resultfolder+"/"+resultfilename;
  myfile.open(resultfile,std::ios_base::app);


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // signal extraction mass region
  const Double_t massLo    = 60;
  const Double_t massHi    = 120;
  
  // fit mass region
  const Double_t fitMassLo = massLo;
  const Double_t fitMassHi = massHi;
  
  TH1D* histPass = new TH1D("histPass","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS,fitMassLo,fitMassHi);
  TH1D* histFail = new TH1D("histFail","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL,fitMassLo,fitMassHi);

  TString readInFile = infilename+"/"+binname;

  ifstream ifs;
  ifs.open(readInFile.Data());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') continue; 

    Double_t mass;
    Int_t state;
    stringstream ss(line);
    ss >> mass >> state;

    if (state == 1) histPass->Fill(mass);
    else if (state == 2) histFail->Fill(mass);
  }
  ifs.close();

/*
  TCanvas *c = new TCanvas("c","c");
  histPass->Draw();
  c->SaveAs("pass.png");
  histFail->Draw();
  c->SaveAs("fail.png");

  cout << "histPass has " << histPass->GetEntries() << endl;
  cout << "histFail has " << histFail->GetEntries() << endl;
  cout << "total we have " << histPass->GetEntries()+histFail->GetEntries() << endl;
*/

  RooRealVar m("m","mass",fitMassLo,fitMassHi);
  m.setBins(10000);

  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);

  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  RooAbsData *dataCombined=0;

  dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),histPass);
  dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),histFail);

  dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                 RooFit::Index(sample),
                                 RooFit::Import("Pass",*((RooDataHist*)dataPass)),
                                 RooFit::Import("Fail",*((RooDataHist*)dataFail)));

  TFile *f = new TFile(binfile);  
  RooWorkspace *w = (RooWorkspace*) f->Get("w");

  RooAbsPdf *modelFail = w->pdf("modelFail");
  RooAbsPdf *modelPass = w->pdf("modelPass");

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");
  totalPdf.addPdf(*modelFail,"Fail");

  RooFitResult *fitResult=0;
  RooMsgService::instance().setSilentMode(kTRUE);
  fitResult = totalPdf.fitTo(*dataCombined,
			     RooFit::PrintEvalErrors(-1),
                             RooFit::Extended(),
                             RooFit::Strategy(2),
			     //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  RooRealVar* eff = (RooRealVar*)fitResult->floatParsFinal().find("eff");
  RooRealVar* eff0 = (RooRealVar*)fitResult->floatParsInit().find("eff");

  if((fabs(eff->getErrorLo())<5e-5) || (eff->getErrorHi()<5e-5))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

  myfile<<"EFFICIENCY: "<<eff->getVal()<<" "<<fabs(eff->getErrorLo())<<endl;
  myfile.close();
}  
