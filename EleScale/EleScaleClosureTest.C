#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting style
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TCanvas.h>                  // class for drawing
#include <TPaveText.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TH1D.h>                     // 1D histograms
#include <TGraphErrors.h>             // graphs
#include <vector>                     // STL vector class
#include <utility>                    // For STL pair class
#include <map>                        // STL map class
#include <iostream>                   // standard I/O
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <sstream>                    // class for parsing strings
#include "TLorentzVector.h"       // 4-vector class
#include "TRandom3.h"

#include "../Utils/CPlot.hh"          // helper class for plots
#include "../Utils/MitStyleRemix.hh"  // style settings for drawing
#include "../Utils/LeptonCorr.hh"
#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooLinearVar.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"

TH1D* returnRelDiff(TH1D* h, TH1D* b, TString name);
//=== MAIN MACRO ================================================================================================= 

void EleScaleClosureTest() {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // event category enumeration
  enum { eEleEle2HLT=1, eEleEle1HLT1L1, eEleEle1HLT, eEleEleNoSel, eEleSC };
  
  TString outputDir = "test";
  TString pufname = "../Tools/pileup_weights_2015B.root";
  
  vector<TString> infilenamev;
  infilenamev.push_back("/data/blue/jlawhorn/Zee/ntuples/data_select.raw.root");  // data
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat_07_23/Zee/ntuples/zee_select.root");  // MC
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat_07_23/Zee/ntuples/zee_select.root");  // MC2

  Float_t lumi=40.0;
  
  const Int_t    NBINS     = 40;
  const Double_t MASS_LOW  = 80;
  const Double_t MASS_HIGH = 100;
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS  = 0.000511;  
  
  vector<pair<Double_t,Double_t> > scEta_limits;
  /*  scEta_limits.push_back(make_pair(0.0,0.5));
  scEta_limits.push_back(make_pair(0.5,1.0));
  scEta_limits.push_back(make_pair(1.0,1.4442));
  scEta_limits.push_back(make_pair(1.566,2.5));
  scEta_limits.push_back(make_pair(2.0,2.5));*/
  scEta_limits.push_back(make_pair(0.0,0.4));
  scEta_limits.push_back(make_pair(0.4,0.8));
  scEta_limits.push_back(make_pair(0.8,1.4442));
  scEta_limits.push_back(make_pair(1.566,2.5));
  //scEta_limits.push_back(make_pair(2.0,2.5));

  CPlot::sOutDir = outputDir;
  
  const TString format("png");

  TRandom3 *rnd = new TRandom3();  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  enum { eData=0, eMC, eMC2 };
  
  TFile *pufile = new TFile(pufname); assert(pufile);
  TH1D  *h_rw = (TH1D*)pufile->Get("npv_rw");

  TH1D* hMC_Tot = new TH1D("hMC_Tot", "", NBINS, MASS_LOW, MASS_HIGH); hMC_Tot->Sumw2();
  TH1D* hData_Tot = new TH1D("hData_Tot", "", NBINS, MASS_LOW, MASS_HIGH); hData_Tot->Sumw2();
  TH1D* hData2_Tot = new TH1D("hData2_Tot", "", NBINS, MASS_LOW, MASS_HIGH); hData2_Tot->Sumw2();
  
  char hname[100];
  vector<TH1D*> hMCv, hDatav, hDatav2;  
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      sprintf(hname,"mc_%i_%i",ibin,jbin);
      hMCv.push_back(new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH));
      hMCv.back()->Sumw2();
      
      sprintf(hname,"data_%i_%i",ibin,jbin);
      hDatav.push_back(new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH));
      hDatav.back()->Sumw2();

      sprintf(hname,"data2_%i_%i",ibin,jbin);
      hDatav2.push_back(new TH1D(hname,"",NBINS,MASS_LOW,MASS_HIGH));
      hDatav2.back()->Sumw2();
    }
  }
  
  //
  // Declare output ntuple variables
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Int_t   q1, q2;
  Float_t scale1fb, puWeight;
  TLorentzVector *dilep=0, *lep1=0, *lep2=0;
  ///// electron specific /////
  TLorentzVector *sc1=0, *sc2=0;
  
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    cout << "Processing " << infilenamev[ifile] << "..." << endl;
    TFile *infile = TFile::Open(infilenamev[ifile]); assert(infile);
    TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  
    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight
    intree->SetBranchAddress("puWeight", &puWeight);  // pileup reweighting
    intree->SetBranchAddress("matchGen", &matchGen);  // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category", &category);  // dilepton category
    intree->SetBranchAddress("npv",      &npv);	      // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);	      // number of in-time PU events (MC)
    intree->SetBranchAddress("q1",       &q1);	      // charge of lead lepton
    intree->SetBranchAddress("q2",       &q2);	      // charge of trail lepton
    intree->SetBranchAddress("dilep",    &dilep);     // dilepton 4-vector
    intree->SetBranchAddress("lep1",     &lep1);      // lead lepton 4-vector
    intree->SetBranchAddress("lep2",     &lep2);      // trail lepton 4-vector
    intree->SetBranchAddress("sc1",      &sc1);	      // lead Supercluster 4-vector
    intree->SetBranchAddress("sc2",      &sc2);	      // trail Supercluster 4-vector 
  
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      Double_t weight = 1;
      if(ifile==eMC || ifile==eMC2) {
	//if(!matchGen) continue;
        weight=scale1fb*lumi*h_rw->GetBinContent(npv+1);
      }
      
      if((category!=eEleEle2HLT) && (category!=eEleEle1HLT) && (category!=eEleEle1HLT1L1)) continue;
      if(q1 == q2) continue;
      if(dilep->M()	  < MASS_LOW)  continue;
      if(dilep->M()	  > MASS_HIGH) continue;
      if(sc1->Pt()	  < PT_CUT)    continue;
      if(sc2->Pt()	  < PT_CUT)    continue;
      if(fabs(sc1->Eta()) > ETA_CUT)   continue;      
      if(fabs(sc2->Eta()) > ETA_CUT)   continue;

      TLorentzVector vLep1(0,0,0,0); 
      TLorentzVector vLep2(0,0,0,0); 
      if (ifile==eData) {
	vLep1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), ELE_MASS);
	vLep2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), ELE_MASS);
      }
      else if (ifile==eMC) {
	vLep1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), ELE_MASS);
	vLep2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), ELE_MASS);
      }
      else {
	vLep1.SetPtEtaPhiM(gRandom->Gaus(lep1->Pt()*getEleScaleCorr(lep1->Eta(),0),getEleResCorr(lep1->Eta(),0)), lep1->Eta(), lep1->Phi(), ELE_MASS);
	vLep2.SetPtEtaPhiM(gRandom->Gaus(lep2->Pt()*getEleScaleCorr(lep2->Eta(),0),getEleResCorr(lep2->Eta(),0)), lep2->Eta(), lep2->Phi(), ELE_MASS);
      }
      TLorentzVector vDilep = vLep1 + vLep2;
    
      Int_t bin1=-1, bin2=-1;
      for(UInt_t i=0; i<scEta_limits.size(); i++) {
        Double_t etalow  = scEta_limits.at(i).first;
        Double_t etahigh = scEta_limits.at(i).second;
        if(fabs(sc1->Eta())>=etalow && fabs(sc1->Eta())<=etahigh) bin1=i;
        if(fabs(sc2->Eta())>=etalow && fabs(sc2->Eta())<=etahigh) bin2=i;
      }
      assert(bin1>=0);
      assert(bin2>=0);
      Int_t ibin= (bin1<=bin2) ? bin1 : bin2;
      Int_t jbin= (bin1<=bin2) ? bin2 : bin1;

      if (ifile==eData) hData_Tot->Fill(vDilep.M(),weight);
      else if (ifile==eMC) hMC_Tot->Fill(vDilep.M(),weight);
      else if (ifile==eMC2) hData2_Tot->Fill(vDilep.M(),weight);

      UInt_t n=jbin-ibin;
      for(Int_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);
      
      if(ifile==eData) {
	hDatav[n]->Fill(vDilep.M(),weight);
      }
      else if(ifile==eMC)   {
	hMCv[n]->Fill(vDilep.M(),weight);
      }
      else if(ifile==eMC2) {
	hDatav2[n]->Fill(vDilep.M(),weight);
      }
    }  
    delete infile;
    infile=0, intree=0;
  }

  TCanvas *c1 = MakeCanvas("c1", "", 800, 800);
  char pname[100];

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.05);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(0);

  TGaxis::SetMaxDigits(3);

  for (UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      UInt_t n=jbin-ibin;
      for(UInt_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);

      cout << hMCv[n]->Integral() << ", " << hDatav[n]->Integral() << ", " << hDatav2[n]->Integral() << endl;

      //hMCv[n]   ->Scale(1.0/hMCv[n]->Integral());
      //hDatav[n] ->Scale(1.0/hDatav[n]->Integral());
      //hDatav2[n]->Scale(1.0/hDatav2[n]->Integral());

      c1->cd(1);

      hMCv[n]->SetLineColor(kRed);
      hMCv[n]->GetYaxis()->SetTitleOffset(1.100);
      hMCv[n]->GetYaxis()->SetTitle("Events");
      hMCv[n]->GetYaxis()->SetRangeUser(0.01, 1.3*TMath::Max(hMCv[n]->GetMaximum(),hDatav[n]->GetMaximum()));
      hMCv[n]->Draw("hist");
      hDatav[n]->Draw("EX0 same");
      hDatav2[n]->SetLineColor(kBlue);
      hDatav2[n]->Draw("histsame");

      c1->cd(2);

      TH1D* hDiffMC = returnRelDiff(hMCv[n],hDatav[n],"foo");
      TH1D* hDiffMC2 = returnRelDiff(hDatav2[n],hDatav[n],"foo2");

      hDiffMC->GetYaxis()->SetRangeUser(-1.0,1.0);
      hDiffMC->GetXaxis()->SetTitle("m_{ee} [GeV]");
      hDiffMC->GetYaxis()->SetTitle("#chi");

      hDiffMC->GetYaxis()->SetTitleOffset(0.42);
      hDiffMC->GetYaxis()->SetTitleSize(0.13);
      hDiffMC->GetXaxis()->SetTitleSize(0.13);
      hDiffMC->GetYaxis()->SetLabelSize(0.12);
      hDiffMC->GetXaxis()->SetLabelSize(0.12);
      hDiffMC->GetYaxis()->SetNdivisions(102);
      hDiffMC->GetYaxis()->CenterTitle();
      hDiffMC->GetXaxis()->SetTitleOffset(1.2);
      hDiffMC->GetXaxis()->CenterTitle();

      hDiffMC->Draw("hist");
      TLine l(80,0.0,100,0.0);
      l.Draw();
      hDiffMC->Draw("histsame");

      hDiffMC2->SetMarkerColor(kBlue);
      hDiffMC2->SetMarkerSize(1);
      hDiffMC2->SetLineColor(kBlue);
      hDiffMC2->Draw("EX0 same");

      c1->cd(1);

      TLegend *leg = new TLegend(0.65, 0.55, 0.90, 0.80);
      leg->SetShadowColor(0); leg->SetLineColor(0);
      leg->AddEntry(hMCv[n],"Raw MC","l");
      leg->AddEntry(hDatav[n],"Data","l");
      leg->AddEntry(hDatav2[n],"Corr. MC","l");
      leg->Draw();

      // CMS label
      TPaveText tb1(0.65,0.92,0.95,0.99,"NDC");
      tb1.SetFillStyle(0);
      tb1.SetBorderSize(0);
      tb1.SetTextAlign(32);
      tb1.AddText("CMS Preliminary");
      tb1.Draw();

      char buffer[200];
      // lumi label
      sprintf(buffer,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
      TPaveText tb2(0.55,0.82,0.90,0.90,"NDC");
      tb2.SetFillStyle(0);
      tb2.SetBorderSize(0);
      tb2.SetTextAlign(32);
      tb2.AddText(buffer);
      tb2.Draw();

      char str1[200],str2[200];
      sprintf(str1,"[%.1f, %.1f]",scEta_limits.at(ibin).first,scEta_limits.at(ibin).second);
      sprintf(str2,"[%.1f, %.1f]",scEta_limits.at(jbin).first,scEta_limits.at(jbin).second);
      TPaveText *a = new TPaveText(0.16,0.75,0.40,0.82,"NDC");
      a->SetFillColor(0); a->SetShadowColor(0); a->SetLineColor(0);
      a->AddText(str1);
      TPaveText *b = new TPaveText(0.16,0.68,0.40,0.75,"NDC");
      b->SetFillColor(0); b->SetShadowColor(0); b->SetLineColor(0);
      b->AddText(str2);

      a->Draw();
      b->Draw();

      sprintf(pname,"ele_comp_%i_%i.png",ibin,jbin); 
      c1->SaveAs(outputDir+"/"+pname);

      delete hDiffMC; delete hDiffMC2;
    }
  }

  cout << endl;
  cout << hMC_Tot->Integral() << ", " << hData_Tot->Integral() << ", " << hData2_Tot->Integral() << endl;

  c1->cd(1);

  hMC_Tot->SetLineColor(kRed); hMC_Tot->SetMarkerColor(kRed);
  hMC_Tot->GetYaxis()->SetTitleOffset(1.100);
  hMC_Tot->GetYaxis()->SetTitle("Events");
  hMC_Tot->GetYaxis()->SetRangeUser(0.01, 1.2*hMC_Tot->GetMaximum());
  hMC_Tot->Draw("hist");
  hData_Tot->Draw("EX0 same");
  hData2_Tot->SetLineColor(kBlue); hData2_Tot->SetMarkerColor(kBlue);
  hData2_Tot->Draw("hist same");

  c1->cd(2);
  
  TH1D* hDiffMC = returnRelDiff(hMC_Tot,hData_Tot,"foo");
  TH1D* hDiffMC2 = returnRelDiff(hData2_Tot,hData_Tot,"foo2");
  
  hDiffMC->GetYaxis()->SetRangeUser(-1.0,1.0);
  hDiffMC->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hDiffMC->GetYaxis()->SetTitle("#chi");

  hDiffMC->GetYaxis()->SetTitleOffset(0.42);
  hDiffMC->GetYaxis()->SetTitleSize(0.13);
  hDiffMC->GetXaxis()->SetTitleSize(0.13);
  hDiffMC->GetYaxis()->SetLabelSize(0.12);
  hDiffMC->GetXaxis()->SetLabelSize(0.12);
  hDiffMC->GetYaxis()->SetNdivisions(102);
  hDiffMC->GetYaxis()->CenterTitle();
  hDiffMC->GetXaxis()->SetTitleOffset(1.2);
  hDiffMC->GetXaxis()->CenterTitle();

  hDiffMC->SetMarkerColor(kRed);
  hDiffMC->SetMarkerSize(1);
  hDiffMC->SetLineColor(kRed);
  hDiffMC->Draw("hist");
  TLine l(80,0.0,100,0.0);
  l.Draw();
  //hDiffMC->Draw("EX0 same");

  hDiffMC2->SetMarkerColor(kBlue);
  hDiffMC2->SetMarkerSize(1);
  hDiffMC2->SetLineColor(kBlue);
  hDiffMC2->Draw("EX0 same");

  c1->cd(1);

  TLegend *leg = new TLegend(0.65, 0.55, 0.90, 0.80);
  leg->SetShadowColor(0); leg->SetLineColor(0);
  leg->AddEntry(hMC_Tot,"Raw MC","l");
  leg->AddEntry(hData_Tot,"Data","l");
  leg->AddEntry(hData2_Tot,"Corr. MC","l");
  leg->Draw();

  // CMS label
  TPaveText tb1(0.65,0.92,0.95,0.99,"NDC");
  tb1.SetFillStyle(0);
  tb1.SetBorderSize(0);
  tb1.SetTextAlign(32);
  tb1.AddText("CMS Preliminary");
  tb1.Draw();

  char buffer[200];
  // lumi label
  sprintf(buffer,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
  TPaveText tb2(0.55,0.82,0.90,0.90,"NDC");
  tb2.SetFillStyle(0);
  tb2.SetBorderSize(0);
  tb2.SetTextAlign(32);
  tb2.AddText(buffer);
  tb2.Draw();

  char str1[200],str2[200];
  sprintf(str1,"[%.1f, %.1f]",scEta_limits.at(0).first,scEta_limits.at(scEta_limits.size()-1).second);
  sprintf(str2,"[%.1f, %.1f]",scEta_limits.at(0).first,scEta_limits.at(scEta_limits.size()-1).second);
  TPaveText *a = new TPaveText(0.16,0.75,0.40,0.82,"NDC");
  a->SetFillColor(0); a->SetShadowColor(0); a->SetLineColor(0);
  a->AddText(str1);
  TPaveText *b = new TPaveText(0.16,0.68,0.40,0.75,"NDC");
  b->SetFillColor(0); b->SetShadowColor(0); b->SetLineColor(0);
  b->AddText(str2);
  
  a->Draw();
  b->Draw();

  sprintf(pname,"comp_tot.png");
  c1->SaveAs(outputDir+"/"+pname);


}

TH1D* returnRelDiff(TH1D* h, TH1D* b, TString name) {
  TH1D* hRelDiff = new TH1D(name, "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  hRelDiff->SetLineColor(h->GetLineColor());
  hRelDiff->SetLineStyle(h->GetLineStyle());
  hRelDiff->SetLineWidth(h->GetLineWidth());

  hRelDiff->GetYaxis()->SetTitleOffset(0.42);
  hRelDiff->GetYaxis()->SetTitleSize(0.13);
  hRelDiff->GetYaxis()->SetLabelSize(0.10);
  hRelDiff->GetXaxis()->SetTitleOffset(1.2);
  hRelDiff->GetXaxis()->SetTitleSize(0.13);
  hRelDiff->GetXaxis()->SetLabelSize(0.12);
  //hRelDiff->GetXaxis()->CenterTitle();                                                                                   
  hRelDiff->GetYaxis()->CenterTitle();
  hRelDiff->GetYaxis()->SetNdivisions(303,kTRUE);

  for (Int_t i=1; i<h->GetNbinsX()+1; i++) {
    Double_t y=b->GetBinContent(i);
    Double_t val = h->GetBinContent(i) - y;
    if (y!=0) hRelDiff->SetBinContent(i, val/y);
    else hRelDiff->SetBinContent(i, 0);
  }
  return hRelDiff;
}
