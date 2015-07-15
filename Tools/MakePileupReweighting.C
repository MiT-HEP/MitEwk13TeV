#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void MakePileupReweighting(TString datafile = "/data/blue/Bacon/Run2/wz_bacon/SingleMuon.root",
			   TString mcfile   = "/data/blue/Bacon/Run2/wz_bacon/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
			   TString outfile  = "pileup_weights_2015B.root") {

  TFile *f_data = new TFile(datafile, "read");
  TTree *t_data = (TTree*) f_data->Get("Events");
  TH1D *h_data = new TH1D("npv_rw", "npv_rw", 40, 0, 40); h_data->Sumw2();

  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  t_data->SetBranchAddress("PV", &vertexArr); TBranch *vertexBr = t_data->GetBranch("PV");

  for (Int_t i=0; i<t_data->GetEntries()/1000; i++) {
    vertexArr->Clear();
    vertexBr->GetEntry(i);
    h_data->Fill(vertexArr->GetEntries());
  }

  h_data->Scale(1.0/h_data->Integral());

  TFile *f_mc = new TFile(mcfile, "read");
  TTree *t_mc = (TTree*) f_mc->Get("Events");
  TH1D *h_mc = new TH1D("npv_mc", "npv_mc", 40, 0, 40); h_mc->Sumw2();

  t_mc->SetBranchAddress("PV", &vertexArr); vertexBr = t_mc->GetBranch("PV");

  for (Int_t i=0; i<t_data->GetEntries(); i++) {
    vertexArr->Clear();
    vertexBr->GetEntry(i);
    h_mc->Fill(vertexArr->GetEntries());
  }

  h_mc->Scale(1.0/h_mc->Integral());

  /*  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  h_data->Draw("hist");
  h_mc->SetLineColor(kRed);
  h_mc->Draw("histsame");

  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  l->AddEntry(h_data,"Data","l");
  l->AddEntry(h_mc,"MC","l");
  l->Draw();
  
  c1->SaveAs("npv_data_mc.png");*/

  TFile *f_out = new TFile(outfile, "recreate");

  TH1D *h_scale = (TH1D*) h_data->Clone();

  h_scale->Divide(h_mc);

  f_out->Write();
  f_out->Close();

}
