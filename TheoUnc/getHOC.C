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


void getHOC() {

  TFile *f_DY = new TFile("/afs/cern.ch/work/j/jlawhorn/public/resbos-files-06-16/dy_1q_yk.root");
  TTree* t_DY = (TTree*) f_DY->Get("h10");

  TH1D* dy_All = new TH1D("dy_All", "", 20, -20, 20); dy_All->Sumw2();
  TH1D* dy_Mu = new TH1D("dy_Mu", "", 20, -20, 20); dy_Mu->Sumw2();
  TH1D* dy_El = new TH1D("dy_El", "", 20, -20, 20); dy_El->Sumw2();

  t_DY->Draw("1>>dy_El", "WT00*(M_B>60 && M_B<120)*(pT_d1>25 && pT_d2>25 && abs(y_d1)<2.5 && abs(y_d2)<2.5)*(abs(y_d1)>1.566 || abs(y_d1)<1.4442)*(abs(y_d2)>1.566 || abs(y_d2)<1.4442)");
  t_DY->Draw("1>>dy_Mu", "WT00*(M_B>60 && M_B<120)*(pT_d1>25 && pT_d2>25 && abs(y_d1)<2.4 && abs(y_d2)<2.4)");
  t_DY->Draw("1>>dy_All", "WT00*(M_B>60 && M_B<120)");

  TFile *f_WP = new TFile("/afs/cern.ch/work/j/jlawhorn/public/resbos-files-06-16/wp_1q_yk.root");
  TTree* t_WP = (TTree*) f_WP->Get("h10");

  TH1D* wp_All = new TH1D("wp_All", "", 20, -20, 20); wp_All->Sumw2();
  TH1D* wp_Mu  = new TH1D("wp_Mu",  "", 20, -20, 20); wp_Mu->Sumw2();
  TH1D* wp_El  = new TH1D("wp_El",  "", 20, -20, 20); wp_El->Sumw2();

  t_WP->Draw("1>>wp_El", "WT00*(pT_d1>25 && abs(y_d1)<2.5)*(abs(y_d1)>1.566 || abs(y_d1)<1.4442)");
  t_WP->Draw("1>>wp_Mu", "WT00*(pT_d1>25 && abs(y_d1)<2.4)");
  t_WP->Draw("1>>wp_All", "WT00");

  TFile *f_WM = new TFile("/afs/cern.ch/work/j/jlawhorn/public/resbos-files-06-16/wm_1q_yk.root");
  TTree* t_WM = (TTree*) f_WM->Get("h10");

  TH1D* wm_All = new TH1D("wm_All", "", 20, -20, 20); wm_All->Sumw2();
  TH1D* wm_Mu  = new TH1D("wm_Mu",  "", 20, -20, 20); wm_Mu->Sumw2();
  TH1D* wm_El  = new TH1D("wm_El",  "", 20, -20, 20); wm_El->Sumw2();

  t_WM->Draw("1>>wm_El", "WT00*(pT_d1>25 && abs(y_d1)<2.5)*(abs(y_d1)>1.566 || abs(y_d1)<1.4442)");
  t_WM->Draw("1>>wm_Mu", "WT00*(pT_d1>25 && abs(y_d1)<2.4)");
  t_WM->Draw("1>>wm_All", "WT00");

  cout << "Zmm: " << dy_Mu->Integral()/dy_All->Integral() << endl;
  cout << "Wpm: " << wp_Mu->Integral()/wp_All->Integral() << endl;
  cout << "Wmm: " << wm_Mu->Integral()/wm_All->Integral() << endl;

  cout << "Zee: " << dy_El->Integral()/dy_All->Integral() << endl;
  cout << "Wpe: " << wp_El->Integral()/wp_All->Integral() << endl;
  cout << "Wme: " << wm_El->Integral()/wm_All->Integral() << endl;

}
