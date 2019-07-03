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
#include "../../Utils/CPlot.hh"     // helper class for plots
#include "../../Utils/MitStyleRemix.hh"// style settings for drawing
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
#include "RooHist.h"
#include "RooGenericPdf.h"

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
    
    const int NBINS=60;
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooTrace::active(1);
  gSystem->mkdir(outputDir,kTRUE);
  
  TCanvas *c = new TCanvas("c","c",800,800);
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

  // RooRealVar m("m","mass",60, 120);
  // m.setBins(NBINS);
  // m.setBins(10000);
  RooRealVar *m =  wmas->var("m");

  // Double_t nsigPass = wsig->var("eff")->getVal() * wsig->var("Nsig")->getVal();
  // Double_t nsigFail = (1. - wsig->var("eff")->getVal()) * wsig->var("Nsig")->getVal();
  // Double_t nbkgPass = wbkg->var("NbkgPass")->getVal();
  // Double_t nbkgFail = wbkg->var("NbkgFail")->getVal();
  RooRealVar *Nsig = wsig->var("Nsig");
  RooRealVar *effGen = wsig->var("eff");
  // effGen->setConstant(kTRUE);
  std::cout << "wsig_eff " << effGen->getVal()<<std::endl;
  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(*effGen,*Nsig));
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(*effGen,*Nsig));
  RooRealVar NbkgPass = *(wbkg->var("NbkgPass"));
  RooRealVar NbkgFail = *(wbkg->var("NbkgFail"));
  // NbkgPass.setConstant(kTRUE);
  // NbkgFail.setConstant(kTRUE);

  
  Double_t effMaster = wsig->var("eff")->getVal();


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
  bkgPass->Print();
  bkgFail->Print();
  // std::cout << "lbbbbbbbbb" << std::endl;

  RooAddPdf *modelPass, *modelFail;
  modelPass = new RooAddPdf("modelPass","Model for PASS sample", RooArgList(*sigPass, *bkgPass), RooArgList(NsigPass,NbkgPass));
  modelFail = new RooAddPdf("modelFail","Model for FAIL sample", RooArgList(*sigFail, *bkgFail), RooArgList(NsigFail,NbkgFail));
// std::cout << "aaaaaaaaaaaaaa" << std::endl;
  RooSimultaneous totalPdfGen("totalPdfGen","totalPdfGen",sample);
  totalPdfGen.addPdf(*modelPass,"Pass");  
  totalPdfGen.addPdf(*modelFail,"Fail");
// std::cout << "cccccccccccc" << std::endl;
  // TTree *intree = (TTree*)fsig->Get("Bin");
  char outHistName[20];
  sprintf(outHistName,"pullHist_%s",binName.Data());
  TH1D *h = new TH1D(outHistName,outHistName,100,-5.0,5.0);
  // TH1D *h2 = new TH1D(outHistName,outHistName,100,-10.0,10.0);
  // wat
  // UInt_t nEvents;
  // intree->SetBranchAddress("nEvents",  &nEvents);
  // intree->GetEntry(0);
  // TFile *f = new TFile(binfile);  
  // RooWorkspace *w = (RooWorkspace*) f->Get("w");

  TTree *intree = (TTree*)fsig->Get("Bin");

  UInt_t nEvents;
  intree->SetBranchAddress("nEvents",  &nEvents);
  intree->GetEntry(0);
  
  RooAbsPdf *fitmodelFail = wmas->pdf("modelFail");
  RooAbsPdf *fitmodelPass = wmas->pdf("modelPass");
  
  double sigInit=fabs(wmas->var("eff")->getErrorLo());
  double meanInit=wmas->var("eff")->getVal();
  std::cout << "err lo " << fabs(wmas->var("eff")->getErrorLo()) << "  err hi " << fabs(wmas->var("eff")->getErrorHi()) << std::endl;

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*fitmodelPass,"Pass");
  totalPdf.addPdf(*fitmodelFail,"Fail");
  // RooAbsData *dataPass=0;
  // RooAbsData *dataFail=0;
  // RooAbsData *dataCombined=0;
  
  // RooFitResult *fitResult=0;
  RooRealVar* eff =wmas->var("eff");
  std::cout << "wmas_eff " << eff->getVal()<<std::endl;
  
  totalPdf.Print();
  totalPdfGen.Print();
  // RooRealVar* eff0 =0;
  
  RooMCStudy mgr(totalPdfGen,  RooArgList(*m, sample),FitModel(totalPdf),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(-1)));
  // RooMCStudy mgr(totalPdfGen,  RooArgList(*m, sample),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(-1)));
  // RooMCStudy* mcstudy = new RooMCStudy(model,x,Binned(kTRUE),Silence(),Extended(),
				       // FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
                       
   mgr.generateAndFit(nPsExp, nEvents);
   std::cout << "done here" << std::endl;
   
    mgr.fitParDataSet().get(3)->Print() ;
    mgr.fitParDataSet().get(6)->Print() ;

     	char outputfile[100];   
   
    sprintf(outputfile, "%s/%s.png",outputDir.Data(), "effplot");
    RooPlot* mframe = eff->frame(0.5,1.05) ;
    mgr.plotParamOn(mframe) ;
    mframe->Draw() ;
    c->SaveAs(outputfile);
    // std::cout << "blahblha" << std::endl;
   RooPlot* mpframe = mgr.plotPull(*eff,-5,5,100,kTRUE);
   // std::cout << "plotted" << std::endl;
   mpframe->Draw();


    sprintf(outputfile, "%s/%s.png",outputDir.Data(), outputName.Data());
    // TCanvas *c = new TCanvas("c","c");
    gStyle->SetOptFit();
	// h->Draw();
	c->SaveAs(outputfile);
    
    // std::cout << "The official one" << std::endl;
	c->Clear();
     mpframe->Print();
     RooHist *pullHist =mpframe->getHist("h_fitParData_totalPdf_totalPdfGen");
     // RooHist *pullHist =mpframe->getHist("h_fitParData_totalPdfGen");
     for(int i=0; i < pullHist->GetN();++i){
         double x=0; double y=0;
         pullHist->GetPoint(i,x,y);
         h->Fill(x,y);
     }
     
    RooRealVar pvar("effpull","effpull",-5,5) ;
    pvar.setBins(100) ;
    
    RooDataHist pullDataHist("pullDataEff","pullDataEff",RooArgList(pvar),h);
    
    RooRealVar pullMean("pullMean","Mean of pull",0,-5,5) ;
    RooRealVar pullSigma("pullSigma","Width of pull",1,0.1,5) ;
    RooGenericPdf pullGauss("pullGauss","Gaussian of pull","exp(-0.5*(@0-@1)*(@0-@1)/(@2*@2))",RooArgSet(pvar,pullMean,pullSigma)) ;
    pullGauss.fitTo(pullDataHist,RooFit::Minos(0),RooFit::PrintLevel(-1)) ;
    RooPlot* frame = pvar.frame(Name("effpull_2"),Title("Pull"),Bins(pvar.getBins())) ;
    pullDataHist.plotOn(frame);
    pullGauss.plotOn(frame) ;
    pullGauss.paramOn(frame,&pullDataHist) ;
    
    frame->Draw();
    sprintf(outputfile, "%s/%s_v2.png",outputDir.Data(), outputName.Data());
	c->SaveAs(outputfile);
	c->Clear();

    // double mean=pullMean.getVal(); double sigma=pullSigma.getVal();
    // double meanErr = pullMean.getError(); double sigmaErr=pullSigma.getError();
    double mean=h->GetMean(); double sigma = h->GetRMS();
    double meanErr = h->GetMeanError(); double sigmaErr = h->GetRMSError();
	sprintf(outputfile,"%s/%s.txt",outputDir.Data(), outputName.Data());
	ofstream meanfile;
	meanfile.open(outputfile);
	meanfile<<mean<<" "<<sigma<<endl;
    // meanfile<<pullMean.getVal() << " " << pullSigma.getVal() << std::endl;
	meanfile.close();
    
    sprintf(outputfile,"%s/err_%s.txt",outputDir.Data(), outputName.Data());
	ofstream errfile;
	errfile.open(outputfile);
	errfile<<meanErr<<" "<<sigmaErr<<endl;
    errfile.close();
    
	sprintf(outputfile,"%s/Sig_%s.txt",outputDir.Data(), outputName.Data());
	ofstream finalfile;
	finalfile.open(outputfile);
	finalfile<<mean*sigInit/meanInit<<endl;
	finalfile.close();
   
    

    

  
  // comment for now
  // for(int i =0; i < nPsExp; ++i){
      // // std::cout << "iteration " << i << std::endl;
      // // RooCategory sample("sample","");
      // // sample.defineType("Pass",1);
      // // sample.defineType("Fail",2);
      
      

      // // dataCombined = totalPdfGen.generate(RooArgList(m,sample),nEvents,Extended());
      // dataCombined = totalPdfGen.generate(RooArgList(*m,sample),nEvents,Extended());
                                     

      // RooMsgService::instance().setSilentMode(kTRUE);
      // RooFitResult *fitResult = totalPdf.fitTo(*dataCombined,
                     // RooFit::PrintEvalErrors(-1),
                                 // RooFit::Extended(),
                                 // RooFit::Strategy(2),
                     // // RooFit::Minos(RooArgSet(*eff)),
                                 // RooFit::Save());

      // eff = (RooRealVar*)fitResult->floatParsFinal().find("eff");
      // eff0 = (RooRealVar*)fitResult->floatParsInit().find("eff");

      // if((fabs(eff->getErrorLo())<5e-5) || (eff->getErrorHi()<5e-5) || fitResult->status()!=0){
        // fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());
      // }
      // std::cout << "RooFit Status " << fitResult->status() << std::endl;
      // std::cout <<"EFFICIENCY: "<< eff->getVal()<<" "<<eff0->getVal()<< " " << effMaster <<  std::endl;
      // std::cout <<"ERROR " << fabs(eff->getErrorLo()) << " "<<fabs(eff->getErrorHi())<< "  "<< eff->getPropagatedError(*fitResult) << std::endl;//<< " pull " << (eff->getVal()-effMaster)/eff->getErrorLo()<< endl;
      // std::cout <<"PULL: "<<(eff->getVal()-eff0->getVal())/fabs(eff->getErrorLo()) << "  b  " <<(eff->getVal()-effMaster)/fabs(eff->getErrorLo()) <<std::endl;//<< " pull " << (eff->getVal()-effMaster)/eff->getErrorLo()<< endl;
      // h->Fill((eff->getVal()-effMaster)/fabs(eff->getErrorLo()));
      // // h2->Fill((eff->getVal()-eff0->getVal())/fabs(eff->getErrorLo()));
      // // std::cout << "eff master " << effMaster << "  eff0 " << eff0->getVal() << "  eff final " << eff->getVal() << endl;
      // TCanvas *c = new TCanvas("c","c",800,800);
      // char filename[100];
	  // sprintf(filename,"Pass_%d",i);
      // // RooPlot* framePass = m.frame(Name(filename),Title("Plot and Fit Pass"),Bins(m.)) ;
      // RooPlot* framePass = m->frame(Name(filename),Title("Plot and Fit Pass"),Bins(m->getBins())) ;
      // dataCombined->plotOn(framePass,Cut("sample==sample::Pass"));
      // dataCombined->statOn(framePass,Layout(0.55,0.99,0.99),Cut("sample==sample::Pass")) ;
      // totalPdf.plotOn(framePass,Slice(sample,"Pass"),ProjWData(sample,*dataCombined));
      // c->cd(1);
      // framePass->Draw();
      // sprintf(filename,"%s/Pass_%d.png",outputDir.Data(),i);
      
      // c->SaveAs(filename);
      
      // RooPlot* frameFail = m->frame(Name(filename),Title("Plot and Fit Fail"),Bins(m->getBins())) ;
      // dataCombined->plotOn(frameFail,Cut("sample==sample::Fail"));
      // dataCombined->statOn(frameFail,Layout(0.55,0.99,0.99),Cut("sample==sample::Fail")) ;
      // totalPdf.plotOn(frameFail,Slice(sample,"Fail"),ProjWData(sample,*dataCombined));
      // c->cd(1);
      // frameFail->Draw();
      // sprintf(filename,"%s/Fail_%d.png",outputDir.Data(),i);
      
      
      // c->SaveAs(filename);
        
      // ofstream txtfile;
      // char txtfname[100];    
      // sprintf(txtfname,"%s/fitres_%i.txt",outputDir.Data(),i);
      // txtfile.open(txtfname);
      // assert(txtfile.is_open());
      // fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
      // txtfile << endl;
      // // printCorrelations(txtfile, fitResult);
      // txtfile.close();
      

      
      // // dataPass=0;
      // // dataFail=0;
      // dataCombined=0;
      // eff->setVal(effMaster);
      // eff0->setVal(effMaster);
      // fitResult=0;
      // // std::cout << "total end" << std::endl;
  // }
  // comment for now
    
    // std::cout << "blah" << std::endl;
    
    // std::cout << "sig " << sigInit  << " mean " << meanInit << std::endl;
    
    // TF1 *f = new TF1("f","gaus",0,10);
    // f->SetParameters(nPsExp,0,1);
    // h->Fit("f");
	// float meanhist = h->GetMean();
	// // float meangaus = abs(f->GetParameter(1));
	// float meangaus = f->GetParameter(1);
	// // float mean = meangaus < meanhist ? meangaus : meanhist; 
	// // // float mean = meantemp > 0? meantemp : -1.*meantemp;
	// float mean = meangaus;
	// float sigma = f->GetParameter(2);

	// char outputfile[100];
	// sprintf(outputfile,"%s/%s.txt",outputDir.Data(), outputName.Data());
	// ofstream meanfile;
	// meanfile.open(outputfile);
	// meanfile<<mean<<" "<<sigma<<endl;
	// meanfile.close();
    
	// sprintf(outputfile,"%s/Sig_%s.txt",outputDir.Data(), outputName.Data());
	// ofstream finalfile;
	// finalfile.open(outputfile);
	// finalfile<<mean*sigInit/meanInit<<endl;
	// finalfile.close();
    
    // sprintf(outputfile, "%s/%s.png",outputDir.Data(), outputName.Data());
        // TCanvas *c = new TCanvas("c","c");
        // gStyle->SetOptFit();
	// h->Draw();
	// c->SaveAs(outputfile);
	// c->Clear();
    
    // h2->Fit("f");
	// meanhist = abs(h->GetMean());
	// meangaus = abs(f->GetParameter(1));
	// // float mean = meangaus < meanhist ? meangaus : meanhist; 
	// // // float mean = meantemp > 0? meantemp : -1.*meantemp;
	// mean = meangaus;
	// sigma = f->GetParameter(2);

	// // char outputfile[100];
	// sprintf(outputfile,"%s/%s_2.txt",outputDir.Data(), outputName.Data());
	// // ofstream meanfile;
	// meanfile.open(outputfile);
	// meanfile<<mean<<" "<<sigma<<endl;
	// meanfile.close();
    
    // sprintf(outputfile, "%s/%s_2.png",outputDir.Data(), outputName.Data());
        // // TCanvas *c = new TCanvas("c","c");
        // gStyle->SetOptFit();
	// h2->Draw();
	// c->SaveAs(outputfile);
	// c->Clear();
    
    RooTrace::dump();
  
  // RooMCStudy mgr(totalPdfGen, totalPdfGen, RooArgList(m, sample), "", "2r");
  // mgr.generate(nPsExp, nEvents, kFALSE, outputDir+binName+"_%d.dat");
}
