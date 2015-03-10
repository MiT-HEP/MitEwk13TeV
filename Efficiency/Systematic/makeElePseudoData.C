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

void generateHistTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
			   const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const TH1D* puWeights);

void makeElePseudoData() {

  // signal extraction mass region
  const Double_t massLo    = 60;
  const Double_t massHi    = 120;

  // fit mass region
  const Double_t fitMassLo = massLo;
  const Double_t fitMassHi = massHi;

  const Bool_t doAbsEta = kFALSE;
  const Int_t charge = 0;

  const TString inputDir="EleGsfSelEff/plots/";
  const TString outputDir="/scratch1/klawhorn/EleGsfSelEff_Up/";
  const Int_t nPsExp=1000;

  // bin edges for kinematic variables
  vector<Double_t> ptBinEdgesv;
  vector<Double_t> etaBinEdgesv;
  vector<Double_t> phiBinEdgesv;
  vector<Double_t> npvBinEdgesv;

  const TString conf="elgsfsel.bins";
  //
  // parse binning file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  Int_t opts[6];
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    Double_t edge;
    stringstream ss(line);
    if(state==0) {
      ss >> opts[0] >> opts[1] >> opts[2] >> opts[3] >> opts[4] >> opts[5];
    } else {
      ss >> edge;
      if(state==1) { ptBinEdgesv.push_back(edge); }
      else if(state==2) { etaBinEdgesv.push_back(edge); }
      else if(state==3) { phiBinEdgesv.push_back(edge); }
      else if(state==4) { npvBinEdgesv.push_back(edge); }
    }
  }
  ifs.close();

  //  
  // Set up binning in kinematic variables
  //        
  const UInt_t ptNbins = ptBinEdgesv.size()-1;
  Double_t ptEdges[ptBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<ptBinEdgesv.size(); iedge++)
    ptEdges[iedge] = ptBinEdgesv[iedge];
    
  const UInt_t etaNbins = etaBinEdgesv.size()-1;
  Double_t etaEdges[etaBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<etaBinEdgesv.size(); iedge++)
    etaEdges[iedge] = etaBinEdgesv[iedge];
  
  const UInt_t phiNbins = phiBinEdgesv.size()-1;
  Double_t phiEdges[phiBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<phiBinEdgesv.size(); iedge++)
    phiEdges[iedge] = phiBinEdgesv[iedge];
  
  const UInt_t npvNbins = npvBinEdgesv.size()-1;
  Double_t npvEdges[npvBinEdgesv.size()];
  for(UInt_t iedge=0; iedge<npvBinEdgesv.size(); iedge++)
    npvEdges[iedge] = npvBinEdgesv[iedge];
  
  TFile *pufile    = 0;
  TH1D  *puWeights = 0;

  //generateHistTemplates("/scratch/klawhorn/EffSysStore/test/down_probes.root", ptBinEdgesv,etaBinEdgesv,phiBinEdgesv,npvBinEdgesv,fitMassLo,fitMassHi,doAbsEta,charge,puWeights);
  generateHistTemplates("/scratch/klawhorn/EffSysStore/test/up_probes.root", ptBinEdgesv,etaBinEdgesv,phiBinEdgesv,npvBinEdgesv,fitMassLo,fitMassHi,doAbsEta,charge,puWeights);
  
  TFile *histfile = 0;
  histfile = new TFile("histTemplates.root");
  assert(histfile);

  char filename[100];
  char binname[100];
  TFile *f = 0;
  TTree *t = 0;
  RooWorkspace *w = 0;
  RooAbsPdf *backgroundFail, *backgroundPass;
  RooAbsPdf *gausFail, *gausPass;
  RooRealVar *NbkgPass, *NbkgFail, *eff, *Nsig;
  RooFormulaVar *NsigPass, *NsigFail;
  TH1D *hFail, *hPass;

  RooDataHist *dataHistFail, *dataHistPass;
  RooHistPdf *histPdfFail, *histPdfPass;
  RooFFTConvPdf *signalFail, *signalPass;
  RooAddPdf *modelFail, *modelPass;

  //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  for (Int_t ibin=0; ibin<10; ibin++) {

    cout << " generating bin " << ibin << endl;

    // RETRIEVE TRUTH PDFS
    sprintf(binname, "etapt_%i", ibin);
    sprintf(filename,"etapt_%i.root",ibin);
    f = new TFile(inputDir+filename);
    w = (RooWorkspace*) f->Get("w");

    t = (TTree*)f->Get("Bin");
    BinInfo bin;
    t->SetBranchAddress("Bin",&bin);
    t->GetEntry(0);

    RooCategory sample("sample", "");
    sample.defineType("Pass", 1);
    sample.defineType("Fail", 2);

    RooRealVar m("m","mass",60, 120);
    m.setBins(10000);
    //RooPlot *frame1 = m.frame(Name("frame"),Title("does this look right???"), Bins(30));

    backgroundFail = w->pdf("backgroundFail");
    backgroundPass = w->pdf("backgroundPass");
    gausFail = w->pdf("gausFail");
    gausPass = w->pdf("gausPass");
    NbkgFail = w->var("NbkgFail");
    NbkgPass = w->var("NbkgPass");

    eff = w->var("eff");
    Nsig = w->var("Nsig");
    NsigPass = new RooFormulaVar("NsigPass", "Nsig*eff", RooArgList(*Nsig, *eff));
    NsigFail = new RooFormulaVar("NsigFail", "Nsig*(1.0-eff)", RooArgList(*Nsig, *eff));

    char hname[50];
    sprintf(hname,"passetapt_%i",ibin);
    hPass = (TH1D*)histfile->Get(hname);
    assert(hPass);

    char vname[50];
    sprintf(vname, "dataHistPass"); dataHistPass = new RooDataHist(vname, vname, RooArgSet(m), hPass);
    sprintf(vname, "histPdfPass"); histPdfPass = new RooHistPdf(vname, vname, m, *dataHistPass, 1);
    sprintf(vname, "signalPass"); signalPass = new RooFFTConvPdf(vname, vname, m, *histPdfPass, *gausPass);

    sprintf(hname,"failetapt_%i", ibin);
    hFail = (TH1D*) histfile->Get(hname);
    assert(hFail);

    sprintf(vname, "dataHistFail"); dataHistFail = new RooDataHist(vname, vname, RooArgSet(m), hFail);
    sprintf(vname, "histPdfFail"); histPdfFail = new RooHistPdf(vname, vname, m, *dataHistFail, 1);
    //histPdfFail->plotOn(frame1);    
    sprintf(vname, "signalFail"); signalFail = new RooFFTConvPdf(vname, vname, m, *histPdfFail, *gausFail);

    modelPass = new RooAddPdf("modelPass", "modelPass", RooArgList(*signalPass, *backgroundPass), RooArgList(*NsigPass, *NbkgPass));
    modelFail = new RooAddPdf("modelFail", "modelFail", RooArgList(*signalFail, *backgroundFail), RooArgList(*NsigFail, *NbkgFail));

    RooSimultaneous totalPdf("totalPdf", "totalPdf", sample);
    totalPdf.addPdf(*modelPass, "Pass");
    totalPdf.addPdf(*modelFail, "Fail");

    //frame1->Draw();

    RooMCStudy mgr(totalPdf, totalPdf, RooArgList(m, sample), "", "2r");
    mgr.generate(nPsExp, bin.nEvents, kFALSE, outputDir+binname+"_%04d.dat");

    f->Close();
  }
  
}

void generateHistTemplates(const TString infilename,
                           const vector<Double_t> &ptEdgesv, const vector<Double_t> &etaEdgesv, const vector<Double_t> &phiEdgesv, const vector<Double_t> &npvEdgesv,
			   const Double_t fitMassLo, const Double_t fitMassHi, const Bool_t doAbsEta, const Int_t charge, const TH1D* puWeights)
{
  cout << "Creating histogram templates ... "; cout.flush();
  char hname[50];
  
  const UInt_t ptNbins  = ptEdgesv.size()-1;
  const UInt_t etaNbins = etaEdgesv.size()-1;
  const UInt_t phiNbins = phiEdgesv.size()-1;
  const UInt_t npvNbins = npvEdgesv.size()-1;
  
  TH1D* passPt[ptNbins];
  TH1D* failPt[ptNbins];
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    sprintf(hname,"passpt_%i",ibin);
    passPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passPt[ibin]->SetDirectory(0);
    sprintf(hname,"failpt_%i",ibin);
    failPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failPt[ibin]->SetDirectory(0);
  }
  
  TH1D* passEta[etaNbins];
  TH1D* failEta[etaNbins];
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
    sprintf(hname,"passeta_%i",ibin);
    passEta[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passEta[ibin]->SetDirectory(0);
    sprintf(hname,"faileta_%i",ibin);
    failEta[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failEta[ibin]->SetDirectory(0);
  }
  
  TH1D* passPhi[phiNbins];  
  TH1D* failPhi[phiNbins];
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    sprintf(hname,"passphi_%i",ibin);
    passPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passPhi[ibin]->SetDirectory(0);
    sprintf(hname,"failphi_%i",ibin);
    failPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failPhi[ibin]->SetDirectory(0);
  }
  
  TH1D* passEtaPt[etaNbins*ptNbins];  
  TH1D* failEtaPt[etaNbins*ptNbins];
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    sprintf(hname,"passetapt_%i",ibin);
    passEtaPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passEtaPt[ibin]->SetDirectory(0);
    sprintf(hname,"failetapt_%i",ibin);
    failEtaPt[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failEtaPt[ibin]->SetDirectory(0);
  }
  
  TH1D* passEtaPhi[etaNbins*phiNbins];
  TH1D* failEtaPhi[etaNbins*phiNbins];
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    sprintf(hname,"passetaphi_%i",ibin); 
    passEtaPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passEtaPhi[ibin]->SetDirectory(0);
    sprintf(hname,"failetaphi_%i",ibin);
    failEtaPhi[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failEtaPhi[ibin]->SetDirectory(0);
  }

  TH1D* passNPV[npvNbins];  
  TH1D* failNPV[npvNbins];
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    sprintf(hname,"passnpv_%i",ibin);
    passNPV[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_PASS),fitMassLo,fitMassHi);
    passNPV[ibin]->SetDirectory(0);
    sprintf(hname,"failnpv_%i",ibin);
    failNPV[ibin] = new TH1D(hname,"",Int_t((fitMassHi-fitMassLo)/BIN_SIZE_FAIL),fitMassLo,fitMassHi);
    failNPV[ibin]->SetDirectory(0);
  }
    
  TFile infile(infilename);
  TTree *intree = (TTree*)infile.Get("Events");
  EffData data;
  intree->SetBranchAddress("Events",&data);
  
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    Float_t puWgt=1;
    //if(puWeights)
    //puWgt = puWeights->GetBinContent(data.npu+1);
    puWgt = data.weight;
    
    if((data.q)*charge < 0) continue;
    
    Int_t ipt=-1;
    for(UInt_t ibin=0; ibin<ptNbins; ibin++)
      if((data.pt >= ptEdgesv[ibin]) && (data.pt < ptEdgesv[ibin+1]))
        ipt = ibin;
    if(ipt<0) continue;
    
    Int_t ieta=-1;
    for(UInt_t ibin=0; ibin<etaNbins; ibin++) {
      if(doAbsEta) {
        assert(etaEdgesv[ibin]>=0);
        if((fabs(data.eta) >= etaEdgesv[ibin]) && (fabs(data.eta) < etaEdgesv[ibin+1]))
          ieta = ibin;
      } else {
        if((data.eta >= etaEdgesv[ibin]) && (data.eta < etaEdgesv[ibin+1]))
          ieta = ibin;
      }
    }
    if(ieta<0) continue;
    
    Int_t iphi=-1;
    for(UInt_t ibin=0; ibin<phiNbins; ibin++)
      if((data.phi >= phiEdgesv[ibin]) && (data.phi < phiEdgesv[ibin+1]))
        iphi = ibin;
    if(iphi<0) continue;

    Int_t inpv=-1;
    for(UInt_t ibin=0; ibin<npvNbins; ibin++)
      if((data.npv >= npvEdgesv[ibin]) && (data.npv < npvEdgesv[ibin+1]))
        inpv = ibin;
    if(inpv<0) continue;
        
    if(data.pass) {
      passPt[ipt]->Fill(data.mass,puWgt);
      passEta[ieta]->Fill(data.mass,puWgt);
      passPhi[iphi]->Fill(data.mass,puWgt);
      passEtaPt[ipt*etaNbins + ieta]->Fill(data.mass,puWgt);
      passEtaPhi[iphi*etaNbins + ieta]->Fill(data.mass,puWgt);
      passNPV[inpv]->Fill(data.mass,puWgt);
    } else {
      failPt[ipt]->Fill(data.mass,puWgt);
      failEta[ieta]->Fill(data.mass,puWgt);
      failPhi[iphi]->Fill(data.mass,puWgt);
      failEtaPt[ipt*etaNbins + ieta]->Fill(data.mass,puWgt);
      failEtaPhi[iphi*etaNbins + ieta]->Fill(data.mass,puWgt);
      failNPV[inpv]->Fill(data.mass,puWgt);
    }    
  }
  infile.Close();

  char outname[50];
  sprintf(outname, "histTemplates.root");
 
  TFile outfile(outname, "RECREATE");
  for(UInt_t ibin=0; ibin<ptNbins; ibin++) {
    passPt[ibin]->Write();
    failPt[ibin]->Write();
    delete passPt[ibin];
    delete failPt[ibin];
  }
  for(UInt_t ibin=0; ibin<etaNbins; ibin++) { 
    passEta[ibin]->Write();
    failEta[ibin]->Write();
    delete passEta[ibin];
    delete failEta[ibin];
  }
  for(UInt_t ibin=0; ibin<phiNbins; ibin++) {
    passPhi[ibin]->Write();
    failPhi[ibin]->Write();
    delete passPhi[ibin];
    delete failPhi[ibin];
  }
  for(UInt_t ibin=0; ibin<(etaNbins*ptNbins); ibin++) {
    passEtaPt[ibin]->Write();
    failEtaPt[ibin]->Write();
    delete passEtaPt[ibin];
    delete failEtaPt[ibin];
  }
  for(UInt_t ibin=0; ibin<(etaNbins*phiNbins); ibin++) {
    passEtaPhi[ibin]->Write();
    failEtaPhi[ibin]->Write();
    delete passEtaPhi[ibin];
    delete failEtaPhi[ibin];
  }
  for(UInt_t ibin=0; ibin<npvNbins; ibin++) {
    passNPV[ibin]->Write();
    failNPV[ibin]->Write();
    delete passNPV[ibin];
    delete failNPV[ibin];
  }
  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}

