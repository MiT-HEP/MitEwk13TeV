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

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void testPUReweighting(TString datafile = "root://eoscms//store/user/jlawhorn/Run2/wz_bacon/SingleMuon.root",
		       TString mcfile   = "root://eoscms//store/user/jlawhorn/Run2/wz_bacon//DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
		       TString certfile = "../Selection/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt",
		       TString infile  = "pileup_weights_2015B.root") {

  TFile *f_rw = TFile::Open(infile, "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("npv_rw");

  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  baconhep::TEventInfo *info = new baconhep::TEventInfo();

  TFile *f_mc = TFile::Open(mcfile, "read");
  TTree *t_mc = (TTree*) f_mc->Get("Events");
  TH1D *h_mc = new TH1D("npv_mc", "npv_mc", 40, 0, 40); h_mc->Sumw2();
  TH1D *h_mc2 = new TH1D("npv_mc2", "npv_mc2", 40, 0, 40); h_mc2->Sumw2();

  t_mc->SetBranchAddress("GenEvtInfo", &gen); TBranch *genBr   = t_mc->GetBranch("GenEvtInfo");
  t_mc->SetBranchAddress("Info", &info); TBranch *infoBr   = t_mc->GetBranch("Info");

  for (Int_t i=0; i<t_mc->GetEntries(); i++) {
    //for (Int_t i=0; i<10000; i++) {
    infoBr->GetEntry(i);
    genBr->GetEntry(i);
    h_mc->Fill(info->nPU, gen->weight);
    h_mc2->Fill(info->nPU,h_rw->GetBinContent(info->nPU+1)*gen->weight);
  }

  cout << h_mc->Integral() << ", " << h_mc2->Integral() << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  h_mc->SetLineColor(kRed);
  h_mc->Draw("hist");
  h_mc2->Draw("histsame");

  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  l->AddEntry(h_mc,"MC","l");
  l->AddEntry(h_mc2,"MC reweighted","l");
  l->Draw();
  

}
