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

#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void MakePileupReweighting(TString datafile = "/data/t3home000/sabrandt/13TeVLowPU/SingleMuon.root",
			   TString mcfile   = "/data/t3home000/sabrandt/LowPUBacons/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root",
			   TString certfile = "../Selection/2017H_lowPU.json",
			   TString outfile  = "pileup_weights_2017H_v2.root") {

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

  Int_t nZero=0, nGood=0;

  std::cout << "Data File " << std::endl;
  for (Int_t i=0; i<t_data->GetEntries(); i++) {
    if(i%100000==0) std::cout << "On Entry " << i << " out of " << t_data->GetEntries() << "..." << 100*i/t_data->GetEntries() << std::endl;
    vertexArr->Clear();
    vertexBr->GetEntry(i);
    infoBr->GetEntry(i);

    baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec); 
    if (!rlrm.hasRunLumi(rl)) continue;

    nGood++;

    if (vertexArr->GetEntries()==0) {
      nZero++;
      continue;
    }
    
    h_data->Fill(vertexArr->GetEntries());

  }

  cout << nZero << " events with no reconstructed vertex out of " << nGood << endl;

  h_data->Scale(1.0/h_data->Integral());

  TFile *f_mc = TFile::Open(mcfile, "read");
  TTree *t_mc = (TTree*) f_mc->Get("Events");
  TH1D *h_mc = new TH1D("npv_mc", "npv_mc", 40, 0, 40); h_mc->Sumw2();

  t_mc->SetBranchAddress("Info", &info); infoBr = t_mc->GetBranch("Info");
  t_mc->SetBranchAddress("GenEvtInfo", &gen); TBranch *genBr = t_mc->GetBranch("GenEvtInfo");
  t_mc->SetBranchAddress("PV", &vertexArr); vertexBr = t_mc->GetBranch("PV");

  Double_t nZeroMC=0, nGoodMC=0;
  std::cout << "MC File " << std::endl;
  for (Int_t i=0; i<t_mc->GetEntries(); i++) {
     if(i%100000==0) std::cout << "On Entry " << i << " out of " << t_mc->GetEntries() << "..." << 100*i/t_mc->GetEntries() << std::endl;
    infoBr->GetEntry(i);
    genBr->GetEntry(i);
    vertexArr->Clear();
    vertexBr->GetEntry(i);
    nGoodMC+=gen->weight;
    if (vertexArr->GetEntries()==0) {
      nZeroMC+=gen->weight;
      continue;
    }
    h_mc->Fill(vertexArr->GetEntries(),gen->weight);
    
  }

  cout << nZeroMC << " events with no reconstructed vertex out of " << nGoodMC << endl;

  h_mc->Scale(1.0/h_mc->Integral());

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 600);
  h_mc->SetFillColor(798); h_mc->SetLineColor(797);
  h_mc->SetLineWidth(3); h_data->SetLineWidth(3);
  h_mc->SetFillStyle(1001);
  h_mc->SetTitle("");
  h_mc->GetXaxis()->SetTitle("n_{PV}");
  h_mc->GetYaxis()->SetTitle("A.U.");
  h_mc->GetYaxis()->SetRangeUser(0,1.2*TMath::Max(h_mc->GetMaximum(), h_data->GetMaximum()));
  h_mc->Draw("hist");
  h_data->Draw("histsame");

  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  l->SetShadowColor(0);
  l->SetLineColor(0);
  l->AddEntry(h_data,"Full Dataset","l");
  l->AddEntry(h_mc,"MC","lf");
  l->Draw();
  
  c1->SaveAs("npv_data_mc.png");

  TFile *f_out = new TFile(outfile, "recreate");

  TH1D *h_scale = (TH1D*) h_data->Clone();

  h_scale->Divide(h_mc);

  f_out->Write();
  f_out->Close();

}
