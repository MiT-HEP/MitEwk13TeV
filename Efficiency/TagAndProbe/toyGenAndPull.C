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

void toyGenAndPull(const TString sigDir, // toy sig template
		    const TString bkgDir, // toy bkg template
		    const TString masDir, // templates used for the fit
		    const TString binName ="eta_0",
		    const TString outputDir=".",
            const TString outputName="pull_0",
		    const Int_t siglabel=0,
		    const Int_t bkglabel=0,
		    const Int_t nPsExp=1000
){
    RooTrace::active(1);
  gSystem->mkdir(outputDir,kTRUE);
  // RETRIEVE TRUTH PDFS
  TFile *fsig = new TFile(sigDir+binName+".root");
  TFile *fbkg = new TFile(bkgDir+binName+".root");
  TFile *fmas = new TFile(masDir+binName+".root");
  RooWorkspace *wsig = (RooWorkspace*) fsig->Get("w");
  // std::cout << "Signal Workspace" << std::endl;
  // wsig->Print();
  RooWorkspace *wbkg = (RooWorkspace*) fbkg->Get("w");
  RooWorkspace *wmas = (RooWorkspace*) fmas->Get("w");
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
  
  Double_t effMaster = wmas->var("eff")->getVal();

  RooRealVar NsigPass("NsigPass", "", nsigPass);
  RooRealVar NsigFail("NsigFail", "", nsigFail);
  RooRealVar NbkgPass("NbkgPass", "", nbkgPass);
  RooRealVar NbkgFail("NbkgFail", "", nbkgFail);


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
  char outHistName[20];
  sprintf(outHistName,"pullHist_%s",binName.Data());
  TH1D *h = new TH1D(outHistName,outHistName,100,-1.0,1.0);
  // wat
  // UInt_t nEvents;
  // intree->SetBranchAddress("nEvents",  &nEvents);
  // intree->GetEntry(0);
  // TFile *f = new TFile(binfile);  
  // RooWorkspace *w = (RooWorkspace*) f->Get("w");

  RooAbsPdf *fitmodelFail = wmas->pdf("modelFail");
  RooAbsPdf *fitmodelPass = wmas->pdf("modelPass");

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*fitmodelPass,"Pass");
  totalPdf.addPdf(*fitmodelFail,"Fail");
  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  RooAbsData *dataCombined=0;
  
  RooFitResult *fitResult=0;
  RooRealVar* eff =0;
  RooRealVar* eff0 =0;
  for(int i =0; i < nPsExp; ++i){
      // std::cout << "iteration " << i << std::endl;
      // RooCategory sample("sample","");
      // sample.defineType("Pass",1);
      // sample.defineType("Fail",2);
      
      

      // dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),histPass);
      // dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),histFail);  
      dataPass = modelPass->generateBinned(RooArgSet(m),NsigPass.getVal()+NbkgPass.getVal());
      dataFail = modelFail->generateBinned(RooArgSet(m),NsigFail.getVal()+NbkgFail.getVal());
      // Throw some toys from the totalPdfGen Model Pass and Fail
      // std::cout << "gggggggggggg" << std::endl;
      
      // remove this garbage shit and replace it with the fit in here




      dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                     RooFit::Index(sample),
                                     RooFit::Import("Pass",*((RooDataHist*)dataPass)),
                                     RooFit::Import("Fail",*((RooDataHist*)dataFail)));



      RooMsgService::instance().setSilentMode(kTRUE);
      fitResult = totalPdf.fitTo(*dataCombined,
                     RooFit::PrintEvalErrors(-1),
                                 RooFit::Extended(),
                                 RooFit::Strategy(2),
                     //RooFit::Minos(RooArgSet(eff)),
                                 RooFit::Save());

      eff = (RooRealVar*)fitResult->floatParsFinal().find("eff");
      eff0 = (RooRealVar*)fitResult->floatParsInit().find("eff");

      if((fabs(eff->getErrorLo())<5e-5) || (eff->getErrorHi()<5e-5))
        fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

      std::cout <<"EFFICIENCY: "<<eff->getVal()<<" "<<fabs(eff->getErrorLo())<<endl;
      h->Fill((eff->getVal()-effMaster)/eff->getErrorLo());
      
      dataPass=0;
      dataFail=0;
      dataCombined=0;
      eff=0;
      eff0=0;
      fitResult=0;
      // std::cout << "total end" << std::endl;
  }
    
    // std::cout << "blah" << std::endl;
    h->Fit("gaus");
	float meanhist = abs(h->GetMean());
	float meangaus = abs(h->GetFunction("gaus")->GetParameter(1));
	// float mean = meangaus < meanhist ? meangaus : meanhist; 
	// // float mean = meantemp > 0? meantemp : -1.*meantemp;
	float mean = meangaus;
	float sigma = h->GetFunction("gaus")->GetParameter(2);

	char outputfile[100];
	sprintf(outputfile,"%s/%s.txt",outputDir.Data(), outputName.Data());
	ofstream meanfile;
	meanfile.open(outputfile);
	meanfile<<mean<<" "<<sigma<<endl;
	meanfile.close();
    
    sprintf(outputfile, "%s/%s.png",outputDir.Data(), outputName.Data());
        TCanvas *c = new TCanvas("c","c");
        gStyle->SetOptFit();
	h->Draw();
	c->SaveAs(outputfile);
	c->Clear();
    
    RooTrace::dump();
  
  // RooMCStudy mgr(totalPdfGen, totalPdfGen, RooArgList(m, sample), "", "2r");
  // mgr.generate(nPsExp, nEvents, kFALSE, outputDir+binName+"_%d.dat");
}
