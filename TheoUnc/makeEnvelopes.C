#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <TH1.h>

#include "compare.hh"

using namespace std;

#endif

void foo();
void makeFile(TString chan);
void getEnv(TFile *f, TString var, TString id, TH1D* &nom, TH1D* &hi, TH1D* &lo);
void getEnvN(TFile *f, TString var, TString id, TH1D* &nom, TH1D* &hi, TH1D* &lo);

void makeEnvelopes() {

  makeFile("wpe");
  makeFile("wpm");
  makeFile("wmm");
  makeFile("zmm");
  makeFile("zee");

}
void makeFile(TString chan) {

  TFile *fN = new TFile("/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes-new/"+chan+"_NNPDF30_nlo_as_0118.root");
  TFile *fC = new TFile("/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes-new/"+chan+"_CT10nlo.root");
  TFile *fM = new TFile("/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes-new/"+chan+"_MMHT2014nlo68cl.root");

  TFile *fOut = new TFile(chan+"_pdf.root", "recreate");
  TString var="Eta";
  TH1D* nomN=0, *hiN=0, *loN=0;
  getEnvN(fN, var, "nnpdf", nomN, hiN, loN);
  TH1D* nomC=0, *hiC=0, *loC=0;
  getEnv(fC, var, "ct", nomC, hiC, loC);
  TH1D* nomM=0, *hiM=0, *loM=0;
  getEnv(fM, var, "mmht", nomM, hiM, loM);

  var="Pt";
  TH1D* nomN2=0, *hiN2=0, *loN2=0;
  getEnvN(fN, var, "nnpdf", nomN2, hiN2, loN2);
  TH1D* nomC2=0, *hiC2=0, *loC2=0;
  getEnv(fC, var, "ct", nomC2, hiC2, loC2);
  TH1D* nomM2=0, *hiM2=0, *loM2=0;
  getEnv(fM, var, "mmht", nomM2, hiM2, loM2);

  fOut->Write();
  fOut->Close();

  fN->Close();
  fC->Close();
  fM->Close();

}
void foo() {
  return;
  /*
  hiN->Scale(1.0/nomN->Integral()); loN->Scale(1.0/nomN->Integral()); nomN->Scale(1.0/nomN->Integral());
  hiC->Scale(1.0/nomC->Integral()); loC->Scale(1.0/nomC->Integral()); nomC->Scale(1.0/nomC->Integral());
  hiM->Scale(1.0/nomM->Integral()); loM->Scale(1.0/nomM->Integral()); nomM->Scale(1.0/nomM->Integral());

  nomN->SetLineWidth(2); nomN->SetLineColor(kBlue);
  hiN->SetLineWidth(2);  hiN->SetLineColor(kBlue);
  loN->SetLineWidth(2);  loN->SetLineColor(kBlue);
  nomC->SetLineWidth(2); nomC->SetLineColor(kRed);
  hiC->SetLineWidth(2);  hiC->SetLineColor(kRed);
  loC->SetLineWidth(2);  loC->SetLineColor(kRed);
  nomM->SetLineWidth(2); nomM->SetLineColor(kGreen);
  hiM->SetLineWidth(2);  hiM->SetLineColor(kGreen);
  loM->SetLineWidth(2);  loM->SetLineColor(kGreen);

  TH1D* hNR = returnRelDiff(hiN, nomN, "hNR");
  TH1D* lNR = returnRelDiff(loN, nomN, "lNR");
  TH1D* hCR = returnRelDiff(hiC, nomN, "hCR");
  TH1D* nCR = returnRelDiff(nomC, nomN, "nCR");
  TH1D* lCR = returnRelDiff(loC, nomN, "lCR");
  TH1D* hMR = returnRelDiff(hiM, nomN, "hMR");
  TH1D* nMR = returnRelDiff(nomM, nomN, "nMR");
  TH1D* lMR = returnRelDiff(loM, nomN, "lMR");

  hNR->GetYaxis()->SetRangeUser(-0.1, 0.1);
  hNR->Draw("hist"); lNR->Draw("histsame");
  hCR->Draw("histsame"); nCR->Draw("histsame"); lCR->Draw("histsame");
  hMR->Draw("histsame"); nMR->Draw("histsame"); lMR->Draw("histsame");

    nomN->Draw("hist");
  hiN->Draw("histsame");
  loN->Draw("histsame");

  nomC->Draw("histsame");
  hiC->Draw("histsame");
  loC->Draw("histsame");

  nomM->Draw("histsame");
  hiM->Draw("histsame");
  loM->Draw("histsame");
  */
  //nomN->GetYaxis()->SetRangeUser(0, 1.3*max(hiN->GetMaximum(), hiC->GetMaximum()));
  //nomN->GetYaxis()->SetRangeUser(0, 1.3*hiN->GetMaximum());

}

void getEnvN(TFile *f, TString var, TString id, TH1D* &nom, TH1D* &hi, TH1D* &lo) {
  
  TH1D* temp;
  char name[50];

  vector<Double_t> vHi, nHi;
  vector<Double_t> vLo, nLo; 

  for (Int_t i=0; i<f->GetNkeys(); i++) {
    temp = ((TH1D*)((TKey*)f->GetListOfKeys()->At(i))->ReadObj());
    sprintf(name, "%s",temp->GetName());
    if (string(name).find(var.Data())<string(name).size()) {
      if (nom==0) {
	nom = (TH1D*) ((TKey*)f->GetListOfKeys()->At(i))->ReadObj()->Clone("central_"+id+"_"+var);
	hi = (TH1D*) ((TKey*)f->GetListOfKeys()->At(i))->ReadObj()->Clone("higher_"+id+"_"+var);
	lo = (TH1D*) ((TKey*)f->GetListOfKeys()->At(i))->ReadObj()->Clone("lower_"+id+"_"+var);
	for (Int_t j=0; j<nom->GetNbinsX(); j++) {
	  vHi.push_back(0); vLo.push_back(0);
	  nHi.push_back(0); nLo.push_back(0);
	}
      }
      if (nom!=0) {
	Double_t cont;
	for (Int_t j=0; j<nom->GetNbinsX(); j++) {
	  cont = temp->GetBinContent(j);
	  if (cont-nom->GetBinContent(j) > 0) { 
	    vHi[j]+=(cont-nom->GetBinContent(j))*(cont-nom->GetBinContent(j));
	    nHi[j]++;
	  }
          else if (cont-nom->GetBinContent(j) < 0) {
	    vLo[j]+=(cont-nom->GetBinContent(j))*(cont-nom->GetBinContent(j));
	    nLo[j]++;
	  }
	}
      }
    }
  }

  for (Int_t j=0; j<nom->GetNbinsX(); j++) {
    hi->SetBinContent(j, nom->GetBinContent(j)+sqrt(vHi[j]/(nHi[j]-1)));
    lo->SetBinContent(j, nom->GetBinContent(j)-sqrt(vLo[j]/(nLo[j]-1)));
  }

}

void getEnv(TFile *f, TString var, TString id, TH1D* &nom, TH1D* &hi, TH1D* &lo) {

  TH1D* temp;
  char name[50];

  for (Int_t i=0; i<f->GetNkeys(); i++) {
    temp = ((TH1D*)((TKey*)f->GetListOfKeys()->At(i))->ReadObj());
    sprintf(name, "%s",temp->GetName());
    if (string(name).find(var.Data())<string(name).size()) {
      if (nom==0) {
	nom = (TH1D*) ((TKey*)f->GetListOfKeys()->At(i))->ReadObj()->Clone("central_"+id+"_"+var);
	hi = (TH1D*) ((TKey*)f->GetListOfKeys()->At(i))->ReadObj()->Clone("higher_"+id+"_"+var);
	lo = (TH1D*) ((TKey*)f->GetListOfKeys()->At(i))->ReadObj()->Clone("lower_"+id+"_"+var);
      }
      if (nom!=0) {
	Double_t cont;
	for (Int_t j=0; j<hi->GetNbinsX(); j++) {
	  cont = temp->GetBinContent(j);
	  if (cont > hi->GetBinContent(j)) hi->SetBinContent(j, cont);
	  else if (cont < lo->GetBinContent(j)) lo->SetBinContent(j, cont);
	}
      }
    }
  }
}
