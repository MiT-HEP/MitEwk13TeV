#if !defined(__CINT__) //|| defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include <TStyle.h>
#include "Math/LorentzVector.h"     // 4-vector class

#include "../Utils/MitStyleRemix.hh"
#include "pdfUncertCT10.hh"

using namespace std;

#endif

void doInclusiveCalcs(TString incDir);
void doDifferentialPlots(TString difDir, TString chan);

void getCT10uncertainties(TString incDir="/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-acc",
			  TString difDir="/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes") {

  doInclusiveCalcs(incDir);
  doDifferentialPlots(difDir, "wme");
  doDifferentialPlots(difDir, "wpe");
  doDifferentialPlots(difDir, "wmm");
  doDifferentialPlots(difDir, "wpm");
  doDifferentialPlots(difDir, "zee");
  doDifferentialPlots(difDir, "zmm");

}

void doDifferentialPlots(TString difDir, TString chan) {

  TFile *f = new TFile(difDir+"/"+chan+"_CT10nlo_all.root");

  TH1D* hEtaNom = (TH1D*) f->Get("dEta_CT10nlo_0");
  TH1D* hEtaEnvU = new TH1D("hEtaEnvU","",20,-5,5);
  TH1D* hEtaEnvL = new TH1D("hEtaEnvL","",20,-5,5);

  TH1D* hPtNom =  (TH1D*) f->Get("dPt_CT10nlo_0");
  TH1D* hPtEnvU;
  if (chan=="zee" || chan=="zmm") hPtEnvU = new TH1D("hPtEnvU","",25,0,100);
  else hPtEnvU = new TH1D("hPtEnvU","",25,25,100);
  TH1D* hPtEnvL;
  if (chan=="zee" || chan=="zmm") hPtEnvL = new TH1D("hPtEnvL","",25,0,100);
  else hPtEnvL = new TH1D("hPtEnvL","",25,25,100);

  for (Int_t j=0; j<20+1; j++) { hEtaEnvU->SetBinContent(j,1.0); }
  for (Int_t j=0; j<20+1; j++) { hEtaEnvL->SetBinContent(j,1.0); }

  for (Int_t j=0; j<25+1; j++) { hPtEnvU->SetBinContent(j,1.0); }
  for (Int_t j=0; j<25+1; j++) { hPtEnvL->SetBinContent(j,1.0); }

  char hname[100];
  for (Int_t i=1; i<53; i++) {
    sprintf(hname,"dEta_CT10nlo_%i",i);
    TH1D* hEta =(TH1D*)f->Get(hname);
    hEta->Divide(hEtaNom);
    for (Int_t j=0; j<20+1; j++) { 
      if ( hEta->GetBinContent(j) > hEtaEnvU->GetBinContent(j) ) hEtaEnvU->SetBinContent(j,hEta->GetBinContent(j));
      else if ( hEta->GetBinContent(j) < hEtaEnvL->GetBinContent(j) ) hEtaEnvL->SetBinContent(j,hEta->GetBinContent(j));
    }
    sprintf(hname,"dPt_CT10nlo_%i",i);
    TH1D* hPt =(TH1D*)f->Get(hname);
    hPt->Divide(hPtNom);
    for (Int_t j=0; j<25+1; j++) { 
      if ( hPt->GetBinContent(j) > hPtEnvU->GetBinContent(j) ) hPtEnvU->SetBinContent(j,hPt->GetBinContent(j));
      else if ( hPt->GetBinContent(j) < hPtEnvL->GetBinContent(j) ) hPtEnvL->SetBinContent(j,hPt->GetBinContent(j));
    }
  }

  for (Int_t i=115; i<124; i++) {
    sprintf(hname,"dEta_CT10nlo_as_0%i_0",i);
    TH1D* hEta =(TH1D*)f->Get(hname);
    hEta->Divide(hEtaNom);
    for (Int_t j=0; j<20+1; j++) { 
      if ( hEta->GetBinContent(j) > hEtaEnvU->GetBinContent(j) ) hEtaEnvU->SetBinContent(j,hEta->GetBinContent(j));
      else if ( hEta->GetBinContent(j) < hEtaEnvL->GetBinContent(j) ) hEtaEnvL->SetBinContent(j,hEta->GetBinContent(j));
    }
    sprintf(hname,"dPt_CT10nlo_as_0%i_0",i);
    TH1D* hPt =(TH1D*)f->Get(hname);
    hPt->Divide(hPtNom);
    for (Int_t j=0; j<25+1; j++) { 
      if ( hPt->GetBinContent(j) > hPtEnvU->GetBinContent(j) ) hPtEnvU->SetBinContent(j,hPt->GetBinContent(j));
      else if ( hPt->GetBinContent(j) < hPtEnvL->GetBinContent(j) ) hPtEnvL->SetBinContent(j,hPt->GetBinContent(j));
    }
  }

  hEtaNom->Divide(hEtaNom);
  hPtNom->Divide(hPtNom);
  
  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 600);

  if (chan=="zee" || chan=="zmm") InitHist(hEtaNom, "Z #eta", "PDF uncertainties");
  else InitHist(hEtaNom, "lepton #eta", "PDF uncertainties");
  hEtaNom->SetLineWidth(2);
  InitHist(hEtaEnvU, "", ""); 
  hEtaEnvU->SetFillStyle(1001);
  hEtaEnvU->SetFillColor(kAzure+6);
  hEtaEnvU->SetLineColor(kAzure+6);
  InitHist(hEtaEnvL, "", ""); 
  hEtaEnvL->SetFillStyle(1001);
  hEtaEnvL->SetFillColor(10);
  hEtaEnvL->SetLineColor(10);

  hEtaNom->GetYaxis()->SetRangeUser(0.9*hEtaEnvL->GetMinimum(),1.1*hEtaEnvU->GetMaximum());
  hEtaNom->Draw("hist");
  hEtaEnvU->Draw("histsame");
  hEtaEnvL->Draw("histsame");
  hEtaNom->Draw("histsame");

  c1->RedrawAxis();
  c1->SaveAs(difDir+"/"+chan+"_CT10nlo_eta_diff.png");

  if (chan=="zee" || chan=="zmm") InitHist(hPtNom, "Z p_{T}", "PDF uncertainties");
  else InitHist(hPtNom, "lepton p_{T}", "PDF uncertainties");
  hPtNom->SetLineWidth(2);
  InitHist(hPtEnvU, "", ""); 
  hPtEnvU->SetFillStyle(1001);
  hPtEnvU->SetFillColor(kAzure+6);
  hPtEnvU->SetLineColor(kAzure+6);
  InitHist(hPtEnvL, "", ""); 
  hPtEnvL->SetFillStyle(1001);
  hPtEnvL->SetFillColor(10);
  hPtEnvL->SetLineColor(10);

  hPtNom->GetYaxis()->SetRangeUser(0.9*hPtEnvL->GetMinimum(),1.1*hPtEnvU->GetMaximum());
  hPtNom->Draw("hist");
  hPtEnvU->Draw("histsame");
  hPtEnvL->Draw("histsame");
  hPtNom->Draw("histsame");

  c1->RedrawAxis();
  c1->SaveAs(difDir+"/"+chan+"_CT10nlo_pt_diff.png");

  f->Close();

  delete f;
  delete c1;

}

void doInclusiveCalcs(TString incDir) {

  Double_t W_p_mu = xsecUncert(incDir+"/wpm_parse_ct10nlo.txt",
			       incDir+"/wpm_parse_ct10nlo_as.txt");

  Double_t W_m_mu = xsecUncert(incDir+"/wmm_parse_ct10nlo.txt",
			       incDir+"/wmm_parse_ct10nlo_as.txt");

  cout << "wpm" << endl;  

  Double_t W_pOverM_mu = chRatioUncert(incDir+"/wpm_parse_ct10nlo.txt",
				       incDir+"/wpm_parse_ct10nlo_as.txt",
				       incDir+"/wmm_parse_ct10nlo.txt",
				       incDir+"/wmm_parse_ct10nlo_as.txt");

  cout << "wm" << endl;  

  Double_t W_mu = chSumUncert(incDir+"/wpm_parse_ct10nlo.txt",
			      incDir+"/wpm_parse_ct10nlo_as.txt",
			      incDir+"/wmm_parse_ct10nlo.txt",
			      incDir+"/wmm_parse_ct10nlo_as.txt",
			      5.82288, 
			      3.94813);

  Double_t Z_mumu = xsecUncert(incDir+"/zmm_parse_ct10nlo.txt",
			       incDir+"/zmm_parse_ct10nlo_as.txt");

  cout << "zwm" << endl;

  Double_t Z_W_mu = wzRatioUncert(incDir+"/wpm_parse_ct10nlo.txt",
				  incDir+"/wpm_parse_ct10nlo_as.txt",
				  incDir+"/wmm_parse_ct10nlo.txt",
				  incDir+"/wmm_parse_ct10nlo_as.txt",
				  incDir+"/zmm_parse_ct10nlo.txt",
				  incDir+"/zmm_parse_ct10nlo_as.txt",
				  5.82288, 
				  3.94813);

  Double_t W_p_el = xsecUncert(incDir+"/wpe_parse_ct10nlo.txt",
			       incDir+"/wpe_parse_ct10nlo_as.txt");

  Double_t W_m_el = xsecUncert(incDir+"/wme_parse_ct10nlo.txt",
			       incDir+"/wme_parse_ct10nlo_as.txt");

  cout << "wpme" << endl;

  Double_t W_pOverM_el = chRatioUncert(incDir+"/wpe_parse_ct10nlo.txt",
				       incDir+"/wpe_parse_ct10nlo_as.txt",
				       incDir+"/wme_parse_ct10nlo.txt",
				       incDir+"/wme_parse_ct10nlo_as.txt");

  Double_t W_el = chSumUncert(incDir+"/wpe_parse_ct10nlo.txt",
			      incDir+"/wpe_parse_ct10nlo_as.txt",
			      incDir+"/wme_parse_ct10nlo.txt",
			      incDir+"/wme_parse_ct10nlo_as.txt",
			      5.82288, 
			      3.94813);

  Double_t Z_elel = xsecUncert(incDir+"/zee_parse_ct10nlo.txt",
			       incDir+"/zee_parse_ct10nlo_as.txt");

  cout << "zwe" << endl;

  Double_t Z_W_el = wzRatioUncert(incDir+"/wpe_parse_ct10nlo.txt",
				  incDir+"/wpe_parse_ct10nlo_as.txt",
				  incDir+"/wme_parse_ct10nlo.txt",
				  incDir+"/wme_parse_ct10nlo_as.txt",
				  incDir+"/zee_parse_ct10nlo.txt",
				  incDir+"/zee_parse_ct10nlo_as.txt",
				  5.82288, 
				  3.94813);
  
  cout << "W+    m " << W_p_mu << endl;
  cout << "W-    m " << W_m_mu << endl;
  cout << "W+/W- m " << W_pOverM_mu << endl;
  cout << "W     m " << W_mu << endl;
  cout << "Z     m " << Z_mumu << endl;
  cout << "Z/W   m " << Z_W_mu << endl;
  cout << endl;
  cout << "W+    e " << W_p_el << endl;
  cout << "W-    e " << W_m_el << endl;
  cout << "W+/W- e " << W_pOverM_el << endl;
  cout << "W     e " << W_el << endl;
  cout << "Z     e " << Z_elel << endl;
  cout << "Z/W   e " << Z_W_el << endl;
  cout << endl;

}
