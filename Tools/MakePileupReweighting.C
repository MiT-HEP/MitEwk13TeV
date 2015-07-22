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

#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void MakePileupReweighting(TString datafile = "root://eoscms//store/user/jlawhorn/Run2/wz_bacon/SingleMuon.root",
			   TString mcfile   = "root://eoscms//store/user/jlawhorn/Run2/wz_bacon//DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",
			   TString certfile = "../Selection/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt",
			   TString outfile  = "pileup_weights_2015B.root") {

  TFile *f_data = TFile::Open(datafile, "read");
  TTree *t_data = (TTree*) f_data->Get("Events");
  TH1D *h_data = new TH1D("npv_rw", "npv_rw", 40, 0, 40); h_data->Sumw2();

  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  baconhep::TEventInfo *info = new baconhep::TEventInfo();
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");

  t_data->SetBranchAddress("Info", &info);    TBranch *infoBr   = t_data->GetBranch("Info");
  t_data->SetBranchAddress("PV", &vertexArr); TBranch *vertexBr = t_data->GetBranch("PV");

  baconhep::RunLumiRangeMap rlrm;
  rlrm.addJSONFile(certfile.Data());

  for (Int_t i=0; i<t_data->GetEntries(); i++) {
    vertexArr->Clear();
    vertexBr->GetEntry(i);
    infoBr->GetEntry(i);

    baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec); 
    if (!(rlrm.hasRunLumi(rl))) continue;
    if (vertexArr->GetEntries()==0) continue;
    h_data->Fill(vertexArr->GetEntries());
  }

  h_data->Scale(1.0/h_data->Integral());

  TFile *f_mc = TFile::Open(mcfile, "read");
  TTree *t_mc = (TTree*) f_mc->Get("Events");
  TH1D *h_mc = new TH1D("npv_mc", "npv_mc", 40, 0, 40); h_mc->Sumw2();

  t_mc->SetBranchAddress("Info", &info); infoBr = t_mc->GetBranch("Info");
  t_mc->SetBranchAddress("GenEvtInfo", &gen); TBranch *genBr = t_mc->GetBranch("GenEvtInfo");

  for (Int_t i=0; i<t_data->GetEntries(); i++) {
    infoBr->GetEntry(i);
    genBr->GetEntry(i);
    h_mc->Fill(info->nPU, gen->weight);
  }

  h_mc->Scale(1.0/h_mc->Integral());
  /*
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  h_data->Draw("hist");
  h_mc->SetLineColor(kRed);
  h_mc->Draw("histsame");

  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  l->AddEntry(h_data,"Data","l");
  l->AddEntry(h_mc,"MC","l");
  l->Draw();
  
  c1->SaveAs("npv_data_mc.png");
  */
  TFile *f_out = new TFile(outfile, "recreate");

  TH1D *h_scale = (TH1D*) h_data->Clone();

  h_scale->Divide(h_mc);

  f_out->Write();
  f_out->Close();

}
