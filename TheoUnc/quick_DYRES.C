#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

//#include "MitStyleRemix.hh"
#include "compare.hh"

using namespace std;

#endif


void quick_DYRES() {

  //TFile *resYk = new TFile("/afs/cern.ch/work/j/jlawhorn/resbos/dy_1Q_yk.root");
  TFile *resY = new TFile("/afs/cern.ch/work/j/jlawhorn/resbos/dy_1Q_yk.root");

  //TTree* tresYk = (TTree*) resYk->Get("h10");
  TTree* tresY = (TTree*) resY->Get("h10");

  TH1D* hresYk = new TH1D("hresYk", "", 20, -20, 20); hresYk->Sumw2();
  TH1D* hresY = new TH1D("hresY", "", 20, -20, 20); hresY->Sumw2();

  //tresY->Draw("1>>hresYk", "WT00*(M_B>60 && M_B<120)*(pT_d1>25 && pT_d2>25 && abs(y_d1)<2.4 && abs(y_d2)<2.4)*(abs(y_d1)>1.566 || abs(y_d1)<1.4442)*(abs(y_d2)>1.566 || abs(y_d2)<1.4442)");
  tresY->Draw("y_d1>>hresYk", "WT00*(M_B>60 && M_B<120)*(pT_d1>25 && pT_d2>25 && abs(y_d1)<2.1 && abs(y_d2)<2.1)");
  tresY->Draw("y_d1>>hresY", "WT00*(M_B>60 && M_B<120)");

  cout << "+NNLO, zee" << endl;

  cout << hresYk->Integral()/hresY->Integral() << endl;
  /*
  TGraphAsymmErrors *fewz = new TGraphAsymmErrors();

  fewz->SetPoint(0 ,   45.00 ,          71.4001 ); fewz->SetPointError(0 , 0, 0,        1.22505,       0.254294 );   
  fewz->SetPoint(1 ,   55.00 ,          46.0603 ); fewz->SetPointError(1 , 0, 0,       0.758499,       0.121831 );   
  fewz->SetPoint(2 ,   68.00 ,          42.0004 ); fewz->SetPointError(2 , 0, 0,       0.663454,      0.0946322 );   
  fewz->SetPoint(3 ,   81.00 ,          16.6212 ); fewz->SetPointError(3 , 0, 0,       0.250576,       0.068595 );   
  fewz->SetPoint(4 ,   91.00 ,          11.4263 ); fewz->SetPointError(4 , 0, 0,       0.167815,      0.0605012 );   
  fewz->SetPoint(5 ,  101.00 ,          8.40696 ); fewz->SetPointError(5 , 0, 0,       0.119422,      0.0501429 );   
  fewz->SetPoint(6 ,  113.00 ,          7.96726 ); fewz->SetPointError(6 , 0, 0,       0.109629,      0.0422739 );   
  fewz->SetPoint(7 ,  135.00 ,          9.05065 ); fewz->SetPointError(7 , 0, 0,       0.117449,      0.0347146 );   
  fewz->SetPoint(8 ,  175.00 ,          5.93631 ); fewz->SetPointError(8 , 0, 0,      0.0697825,      0.0266941 );   
  fewz->SetPoint(9 ,  400.00 ,          2.77051 ); fewz->SetPointError(9 , 0, 0,      0.0267556,      0.0267637 );   

  TGraphErrors *dy = new TGraphErrors();
  dy->SetPoint( 8, 42.5000, 4124.52); dy->SetPointError( 8, 2.5, 122.929);
  dy->SetPoint( 9, 47.5000, 3052.78); dy->SetPointError( 9, 2.5, 84.9985);
  dy->SetPoint(10, 52.5000, 2423.06); dy->SetPointError(10, 2.5, 80.1229); 
  dy->SetPoint(11, 57.5000, 2029.33); dy->SetPointError(11, 2.5, 68.0173); 
  dy->SetPoint(12, 62.5000, 2249.86); dy->SetPointError(12, 2.5, 428.429); 
  dy->SetPoint(13, 67.5000, 1369.66); dy->SetPointError(13, 2.5, 109.028); 
  dy->SetPoint(14, 72.5000, 1215.69); dy->SetPointError(14, 2.5, 83.8946); 
  dy->SetPoint(15, 77.5000, 877.062); dy->SetPointError(15, 2.5, 44.7429); 
  dy->SetPoint(16, 82.5000, 933.481); dy->SetPointError(16, 2.5, 166.487); 
  dy->SetPoint(17, 87.5000, 617.726); dy->SetPointError(17, 2.5, 44.7782); 
  dy->SetPoint(18, 92.5000, 610.916); dy->SetPointError(18, 2.5, 69.8936); 
  dy->SetPoint(19, 97.5000, 513.827); dy->SetPointError(19, 2.5, 78.9036);
  */
   /*150.000        148.952        3.41388
     250.000        17.3767       0.120172
     350.000        3.93049       0.375338E-01
     450.000        1.13661       0.148523E-01
     550.000       0.405330       0.625933E-02*/

  /*  dy->SetFillColor(kGreen+1);
  dy->SetMarkerColor(kGreen+1);
  dy->SetLineColor(kGreen+1);*/
  //dy->SetFillStyle(3004);
  /*
  TFile *f = new TFile("/afs/cern.ch/work/j/jlawhorn/flat-05-14/zmmLOEWK.root");
  TFile *f2 = new TFile("/afs/cern.ch/work/j/jlawhorn/flat-05-14/DYJets.root");
  TTree *tr = (TTree*) f->Get("Events");
  TTree *tr2 = (TTree*) f2->Get("Events");
  //TH1D *h = new TH1D("h", "", 20, 0, 100); h->Sumw2();
  //TH1D *h2 = new TH1D("h2", "", 20, 0, 100); h2->Sumw2();
  TH1D *h = new TH1D("h", "", 25, 0, 25); h->Sumw2();
  TH1D *h2 = new TH1D("h2", "", 25, 0, 25); h2->Sumw2();
  tr->Draw("genVf_pt>>h", "weight*(genV_m>60 && genV_m<120)*(genL1_pt>20 && genL2_pt>20 && abs(genL1_eta)<2.1 && abs(genL2_eta)<2.1)");
  tr2->Draw("genVf_pt>>h2", "weight*(genV_m>60 && genV_m<120)*(genL1_pt>20 && genL2_pt>20 && abs(genL1_eta)<2.1 && abs(genL2_eta)<2.1)");

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01); //0.01
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetRightMargin(0.07);
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);

  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.45);//0.25
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.100,"Y");

  c1->cd(1);

  Float_t scale=0;
  for (Int_t i=0; i<25; i++) { 
    Double_t x, y; 
    dy->GetPoint(i,x,y);
    dy->SetPoint(i,x,y/1000);
    dy->SetPointError(i,0.5,dy->GetErrorY(i)/1000);
    scale+=y/1000; //5.0=bin size
  }
  dy->GetXaxis()->SetRangeUser(0,25);
  dy->SetTitle("");
  dy->GetYaxis()->SetTitle("Scaled to DYRES [pb]");
  dy->GetYaxis()->SetNdivisions(310,kTRUE);
  dy->Draw("ape");

  h2->SetLineColor(kRed);
  h->SetLineColor(kBlue);
  h->Scale(scale/h->Integral());
  h2->Scale(scale/h2->Integral());

  h->Draw("histsame");
  h2->Draw("histsame");

  TH1D* cPow = returnRelDiff(h, dy, "cPow");
  TH1D* cAmc = returnRelDiff(h2, dy, "cAmc");

  hresYk->SetLineColor(kGreen+1);
  hresYk->SetMarkerColor(kGreen+1);
  hresY->SetLineColor(kGreen+1);
  hresY->SetMarkerColor(kGreen+1);

  hresYk->Scale(scale/hresYk->Integral());
  hresY->Scale(scale/hresY->Integral());
  //hresYk->Draw("same e1");
  hresY->Draw("histsame");
  
  dy->Draw("pesame");

  TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.85);
  l->SetFillColor(0); l->SetShadowColor(0); l->SetLineColor(0);
  l->AddEntry(dy, "DYRES", "pe");
  l->AddEntry(h, "Powheg", "l");
  l->AddEntry(h2, "aMC@NLO", "l");
  //l->AddEntry(hresYk, "ResBos, w/yk grid", "lp");
  l->AddEntry(hresY, "ResBos", "l");
  l->Draw();

  //TH1D* cResYk = returnRelDiff(hresYk, dy, "cResYk");
  TH1D* cResY  = returnRelDiff(hresY,  dy, "cResY");

  //dy->GetYaxis()->SetRangeUser(0, 1.2*TMath::Max(h2->GetMaximum(), h->GetMaximum()));
  //dy->GetXaxis()->SetTitle("Z p_{T}");
  //dy->GetYaxis()->SetTitle("Normalized to DYRES");

  c1->cd(2);
  cPow->GetYaxis()->SetRangeUser(-0.5,0.5);
  cPow->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  cPow->GetYaxis()->SetTitle("#chi");
  cPow->Draw("hist");
  cAmc->Draw("histsame");
  //cResYk->Draw("histsame");
  cResY->Draw("histsame");
*/

  /*
  Float_t scale=0;
  for (Int_t i=0; i<20; i++) { 
    Double_t x, y; 
    dy->GetPoint(i,x,y);
    //cout << y << ", ";
    scale+=y; //5.0=bin size
  }

  //cout << scale << endl;
  //cout << hres->Integral() << endl;
  //hres->Scale(scale/hres->Integral());
  //cout << hres->Integral() << endl;
  //hres->Draw("same e");

  //dy->GetYaxis()->SetRangeUser(1e2, 1e5);
  dy->GetYaxis()->SetRangeUser(0, 1.2*hres->GetMaximum());
  */
}
