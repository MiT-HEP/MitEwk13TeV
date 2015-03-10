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

#include "CPlot.hh"	            // helper class for plots
#include "MitStyleRemix.hh"         // style settings for drawing
#include "CEffUser1D.hh"            // class for handling efficiency graphs
#include "CEffUser2D.hh"            // class for handling efficiency tables

// structure for input ntuple
#include "EffData.hh"

#include "ZSignals.hh"
#include "ZBackgrounds.hh"

#include "BinInfo.hh"
#endif

// RooFit headers
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
#define BIN_SIZE_PASS 2
#define BIN_SIZE_FAIL 2

//=== FUNCTION DECLARATIONS ======================================================================================

void generateHistTemplates(const TString infilename, const Float_t ptLo, const Float_t ptHi, const Float_t etaLo, const Float_t etaHi,
                           const Float_t phiLo, const Float_t phiHi, const Float_t npvLo, const Float_t npvHi, const Int_t absEta,
                           const Int_t charge, const Int_t iBin);

//=== MAIN MACRO ================================================================================================= 

void doPseudoFits(const TString infilename,      // input file
		  const TString binname,         // name
		  const TString binfile,         // file with bin info
		  const Int_t   sigModPass,      // signal extraction method for PASS sample
		  const Int_t   bkgModPass,      // background model for PASS sample
		  const Int_t   sigModFail,      // signal extraction method for FAIL sample	     
		  const Int_t   bkgModFail,      // background model for FAIL sample
		  const TString outputDir,       // output directory
		  const Int_t   charge,          // 0 (no charge requirement), -1, +1
		  const TString mcfilename="")   // ROOT file containing MC events to generate templates from
{
  gBenchmark->Start("plotEff");


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // signal extraction mass region
  const Double_t massLo    = 60;
  const Double_t massHi    = 120;
  
  // fit mass region
  const Double_t fitMassLo = massLo;
  const Double_t fitMassHi = massHi;
  
  gSystem->mkdir(outputDir,kTRUE);  

  //
  // parse pseudodata file
  //
  TH1D* histPass = new TH1D("histPass","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS,fitMassLo,fitMassHi);
  TH1D* histFail = new TH1D("histFail","",Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL,fitMassLo,fitMassHi);

  TString readInFile = infilename+"/"+binname;

  cout << readInFile << endl;

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

  cout << "histPass has " << histPass->GetEntries() << endl;
  cout << "histFail has " << histFail->GetEntries() << endl;
  cout << "total we have " << histPass->GetEntries()+histFail->GetEntries() << endl;
  
  //
  // get binning info
  //
  TFile *f = new TFile(binfile);  
  TTree *intree = (TTree*)f->Get("Bin");
  BinInfo bin;
  intree->SetBranchAddress("Bin",&bin);
  intree->GetEntry(0);

  cout << "we should have " << bin.nEvents << endl;

  cout << bin.ptLo << " " << bin.ptHi << " " << bin.etaLo << " " << bin.etaHi << " " << bin.phiLo << " " << bin.phiHi << " " << bin.npvLo << " " << bin.npvHi << " " << bin.absEta << endl;
  
  //
  // Generate histogram templates from MC if necessary
  //
  if(sigModPass==2 || sigModFail==2) {
    generateHistTemplates(mcfilename, bin.ptLo, bin.ptHi, bin.etaLo, bin.etaHi, bin.phiLo, bin.phiHi, bin.npvLo, bin.npvHi, bin.absEta, 0, bin.iBin);
  }
  
  RooRealVar m("m","mass",fitMassLo,fitMassHi);
  m.setBins(10000);

  Int_t nflpass=0, nflfail=0;

  TFile *histfile = 0;
  if(sigModPass==2 || sigModFail==2) {
    histfile = new TFile("histTemplates.root");
    assert(histfile);
  }

  
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
  
  // Define signal and background models                                                                                                                                                                    
  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;

  char hname[50];
  sprintf(hname,"pass_%i",bin.iBin);
  TH1D *h = (TH1D*)histfile->Get(hname);
  assert(h);
  sigPass = new CMCTemplateConvGaussian(m,h,kTRUE);
  //((CMCTemplateConvGaussian*)sigPass)->mean->setVal(-0.1);
  //((CMCTemplateConvGaussian*)sigPass)->sigma->setVal(0.015);
  //((CMCTemplateConvGaussian*)sigPass)->sigma->setMax(0.2);
  nflpass += 2;
  if (bkgModPass == 1) {
    bkgPass = new CExponential(m,kTRUE);
    nflpass += 1;
  }
  else if (bkgModPass == 2) {
    bkgPass = new CErfExpo(m,kTRUE); 
    nflfail += 3;
    //((CErfExpo*)bkgPass)->alfa->setVal(64.);
    //((CErfExpo*)bkgPass)->alfa->setMax(80.);
    //((CErfExpo*)bkgPass)->beta->setVal(1.0);
    //((CErfExpo*)bkgPass)->beta->setMax(3.0);
    //((CErfExpo*)bkgPass)->gamma->setVal(1.0);
    //((CErfExpo*)bkgPass)->gamma->setMax(3.0);
  }
  else {
    cout << "trying to use bkg model that's not implemented!!" << endl;
  }
  char hname2[50];
  sprintf(hname2,"fail_%i",bin.iBin);
  TH1D *h2 = (TH1D*)histfile->Get(hname2);
  assert(h2);
  sigFail = new CMCTemplateConvGaussian(m,h2,kFALSE);
  //((CMCTemplateConvGaussian*)sigFail)->mean->setVal(-0.28);
  //((CMCTemplateConvGaussian*)sigFail)->sigma->setVal(0.2);
  //((CMCTemplateConvGaussian*)sigFail)->sigma->setMax(2.0);
  nflfail += 2;

  if (bkgModFail == 1) {
    bkgFail = new CExponential(m,kFALSE);
    nflfail += 1;
  }
  else if (bkgModFail == 2) {
    bkgFail = new CErfExpo(m,kFALSE);
    nflfail += 3;
    //((CErfExpo*)bkgFail)->alfa->setVal(60.);
    //((CErfExpo*)bkgFail)->alfa->setMax(80.);
    //((CErfExpo*)bkgFail)->beta->setVal(0.07);
    //((CErfExpo*)bkgFail)->beta->setMax(0.5);
    //((CErfExpo*)bkgFail)->gamma->setVal(0.02);
    //((CErfExpo*)bkgFail)->gamma->setMax(1.0);
  }
  else {
    cout << "trying to use bkg model that's not implemented!!" << endl;
  }

  // Define free parameters
  Double_t NsigMax     = histPass->Integral()+histFail->Integral();
  cout << "NsigMax " << NsigMax << endl;
  Double_t NbkgFailMax = histFail->Integral();
  cout << "NbkgFailMax " << NbkgFailMax << endl;
  Double_t NbkgPassMax = histPass->Integral();
  RooRealVar Nsig("Nsig","Signal Yield",0.80*NsigMax,0,NsigMax);
  RooRealVar eff("eff","Efficiency",0.9,0,1.0);
  //cout << "got here" << endl;
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",10,0,NbkgPassMax);
  //cout << "chicken" << endl;
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",10,0.01,NbkgFailMax);
  //cout << "frog" << endl;

  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
  RooAddPdf *modelPass=0, *modelFail=0;
  RooExtendPdf *esignalPass=0, *ebackgroundPass=0, *esignalFail=0, *ebackgroundFail=0;

  if(massLo!=fitMassLo || massHi!=fitMassHi) {
    m.setRange("signalRange",massLo,massHi);

    esignalPass     = new RooExtendPdf("esignalPass","esignalPass",*(sigPass->model),NsigPass,"signalRange");
    ebackgroundPass = new RooExtendPdf("ebackgroundPass","ebackgroundPass",(bkgModPass>0) ? *(bkgPass->model) : *(sigPass->model),NbkgPass,"signalRange");
    modelPass       = new RooAddPdf("modelPass","Model for PASS sample",(bkgModPass>0) ? RooArgList(*esignalPass,*ebackgroundPass) : RooArgList(*esignalPass));

    esignalFail     = new RooExtendPdf("esignalFail","esignalFail",*(sigFail->model),NsigFail,"signalRange");
    ebackgroundFail = new RooExtendPdf("ebackgroundFail","ebackgroundFail",*(bkgFail->model),NbkgFail,"signalRange");
    modelFail       = new RooAddPdf("modelFail","Model for FAIL sample", (bkgModFail>0) ? RooArgList(*esignalFail,*ebackgroundFail) : RooArgList(*esignalFail));

  } else {
    modelPass = new RooAddPdf("modelPass","Model for PASS sample",
                              RooArgList(*(sigPass->model),*(bkgPass->model)),
                              RooArgList(NsigPass,NbkgPass));

    modelFail = new RooAddPdf("modelFail","Model for FAIL sample",
			      RooArgList(*(sigFail->model),*(bkgFail->model)),
			      RooArgList(NsigFail,NbkgFail));
  }
  cout << "whale?" << endl;
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");
  totalPdf.addPdf(*modelFail,"Fail");
  
  RooFitResult *fitResult=0;
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::Extended(),
                             RooFit::Strategy(1),
                             //RooFit::Minos(RooArgSet(eff)),
			     RooFit::Save());

  // Refit w/o MINOS if MINOS errors are strange...                                                                                                                                                         
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());
  
  cout << eff.getVal() << " " << fabs(eff.getErrorLo()) << " " << eff.getErrorHi() << endl;
  /*          
  RooPlot *mframePass = m.frame(Bins(Int_t(fitMassHi-fitMassLo)/BIN_SIZE_PASS));
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
  modelPass->plotOn(mframePass);
  modelPass->plotOn(mframePass,Components("backgroundPass"),LineStyle(kDashed),LineColor(kRed));
  
  RooPlot *mframeFail = m.frame(Bins(Int_t(fitMassHi-fitMassLo)/BIN_SIZE_FAIL));
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail);
  modelFail->plotOn(mframeFail,Components("backgroundFail"),LineStyle(kDashed),LineColor(kRed));

  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);

  char ylabel[50];

  //
  // Plot passing probes
  //
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_PASS);
  CPlot plotPass("pass",mframePass,"Passing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotPass.Draw(cpass,kTRUE,"png");
 
  //
  // Plot failing probes
  //
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",(Double_t)BIN_SIZE_FAIL);
  CPlot plotFail("fail",mframeFail,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotFail.Draw(cfail,kTRUE,"png"); 
  */
  //
  // Write fit results
  //

  TObjString* temp=(TObjString*)readInFile.Tokenize("/.")->At(4);

  ofstream txtfile;
  char txtfname[100];
  sprintf(txtfname,"%s/%s.output",outputDir.Data(),temp->GetString().Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile.close();
  
}  


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void generateHistTemplates(const TString infilename, const Float_t ptLo, const Float_t ptHi, const Float_t etaLo, const Float_t etaHi,
                           const Float_t phiLo, const Float_t phiHi, const Float_t npvLo, const Float_t npvHi, const Int_t absEta,
                           const Int_t charge, const Int_t iBin)
{

  cout << "Creating histogram templates..."; cout.flush();

  cout << " mcfilename " << infilename; cout.flush();
  
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

  //cout << endl; cout.flush();
  for(UInt_t ientry=0;ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    //cout << endl; cout.flush();
    //cout << " start loop"; cout.flush();

    if((data.q)*charge < 0) continue;
    //cout << " charge "; cout.flush();
 
    Int_t ipt=-1;
    if (ptLo == 999) 
      ipt = 1;
    else if( (data.pt >= ptLo) && (data.pt < ptHi) )
      ipt = 1;
    if(ipt<0) continue;
    //cout << " pt "; cout.flush();

    Int_t ieta=-1;
    if (etaLo == 999) {
      ieta = 1; //cout << " no eta cut "; cout.flush();
    } else if(absEta == 1) {
      //cout << " absEta "; cout.flush();
      assert(etaLo>=0);
      if((fabs(data.eta) >= etaLo) && (fabs(data.eta) < etaHi) )
        ieta = 1; 
    } else {
      //cout << " notAbsEta " << data.eta; cout.flush();
      if( (data.eta >= etaLo) && (data.eta < etaHi) )
        ieta = 1; 
    }
    if(ieta<0) continue;
    //cout << " eta "; cout.flush();
    Int_t iphi=-1;
    if (phiLo == 999)
      iphi = 1;
    else if( (data.phi >= phiLo) && (data.phi < phiHi) )
      iphi = 1; 
    if(iphi<0) continue;
    //cout << " phi "; cout.flush();
    Int_t inpv=-1;
    if (npvLo == 999) {
      //cout << " no npv cut "; cout.flush();
      inpv = 1;
    } else if( (data.npv >= npvLo) && (data.npv < npvHi) )
      inpv = 1; 
    if(inpv<0) continue;

    //cout << "end loop" << endl; cout.flush();

    if(data.pass) passEtaPt->Fill(data.mass,1);
    else failEtaPt->Fill(data.mass,1);
  }
  infile.Close();
  //cout << endl; cout.flush();
  
  TFile outfile("histTemplates.root", "RECREATE");
  passEtaPt->Write();
  failEtaPt->Write();
  delete passEtaPt;
  delete failEtaPt;

  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
  
}
