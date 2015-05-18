#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
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
#include "Math/LorentzVector.h"     // 4-vector class

#include "../Utils/MitStyleRemix.hh"

using namespace std;

#endif

enum PROC1 {wme=0, wpe, wmm, wpm};
enum PROC2 {zee=0, zmm};

void calcAcceptance(TFile *f, TString chan, Int_t n, vector<Double_t> &scale);
void calcAcceptance(TFile *f, TString chan, Int_t n, Int_t as, vector<Double_t> &scale);
void makeScaleVector(TString inDir, TString chan, Double_t &nom, Double_t &nomAs, vector<Double_t> &vScales, vector<Double_t> &alphaS);
Double_t uncert(vector<Double_t> &vScale, vector<Double_t> &vScaleAs, Double_t nom, Double_t nomAs);

void getCT10uncertainties() {

  Double_t wmXsec = 8.35;
  Double_t wpXsec = 11.20;
  Double_t zXsec  = 1.92;

  TString inDir = "/afs/cern.ch/work/j/jlawhorn/public/wz-13tev-envelopes-new/";

  cout << "CT10 acceptances" << endl;

  Double_t wmeNom, wmeNomAs;
  vector<Double_t> wmeScales, wmeScaleAs;
  makeScaleVector(inDir, "wme", wmeNom, wmeNomAs, wmeScales, wmeScaleAs);

  Double_t wpeNom, wpeNomAs;
  vector<Double_t> wpeScales, wpeScaleAs;
  makeScaleVector(inDir, "wpe", wpeNom, wpeNomAs, wpeScales, wpeScaleAs);

  Double_t wmmNom, wmmNomAs; 
  vector<Double_t> wmmScales, wmmScaleAs;
  makeScaleVector(inDir, "wmm", wmmNom, wmmNomAs, wmmScales, wmmScaleAs);

  Double_t wpmNom, wpmNomAs;
  vector<Double_t> wpmScales, wpmScaleAs;
  makeScaleVector(inDir, "wpm", wpmNom, wpmNomAs, wpmScales, wpmScaleAs);

  Double_t zeeNom, zeeNomAs;
  vector<Double_t> zeeScales, zeeScaleAs;
  makeScaleVector(inDir, "zee", zeeNom, zeeNomAs, zeeScales, zeeScaleAs);

  Double_t zmmNom, zmmNomAs;
  vector<Double_t> zmmScales, zmmScaleAs;
  makeScaleVector(inDir, "zmm", zmmNom, zmmNomAs, zmmScales, zmmScaleAs);

  if (wmeScales.size()!=wpeScales.size()) cout << "SOMETHING BAD HAPPENED TO FILES" << endl;
  if (wmeScales.size()!=wpmScales.size()) cout << "SOMETHING BAD HAPPENED TO FILES" << endl;
  if (wmmScales.size()!=wpmScales.size()) cout << "SOMETHING BAD HAPPENED TO FILES" << endl;
  if (zmmScales.size()!=wpmScales.size()) cout << "SOMETHING BAD HAPPENED TO FILES" << endl;
  if (zmmScales.size()!=zeeScales.size()) cout << "SOMETHING BAD HAPPENED TO FILES" << endl;

  Double_t wpmeNom=wpeNom/wmeNom; 
  vector<Double_t> wpmeScales;
  for (UInt_t i=0; i<wmeScales.size(); i++) { wpmeScales.push_back(wpeScales[i]/wmeScales[i]); }

  Double_t wpmmNom=wpmNom/wmmNom; 
  vector<Double_t> wpmmScales;
  for (UInt_t i=0; i<wmmScales.size(); i++) { wpmmScales.push_back(wpmScales[i]/wmmScales[i]); }

  Double_t weNom=(wpeNom*wpXsec+wmeNom*wmXsec)/(wpXsec+wmXsec); 
  vector<Double_t> weScales;
  for (UInt_t i=0; i<wpeScales.size(); i++) { weScales.push_back((wpeScales[i]*wpXsec+wmeScales[i]*wmXsec)/(wpXsec+wmXsec)); }

  Double_t wmNom=(wpmNom*wpXsec+wmmNom*wmXsec)/(wpXsec+wmXsec); 
  vector<Double_t> wmScales;
  for (UInt_t i=0; i<wpmScales.size(); i++) { wmScales.push_back((wpmScales[i]*wpXsec+wmmScales[i]*wmXsec)/(wpXsec+wmXsec)); }

  Double_t wzeNom=(wpeNom*wpXsec+wmeNom*wmXsec)/(zeeNom*zXsec)*((zXsec)/(wpXsec+wmXsec)); 
  vector<Double_t> wzeScales;
  for (UInt_t i=0; i<wpeScales.size(); i++) { wzeScales.push_back(((wpeScales[i]*wpXsec+wmeScales[i]*wmXsec)/(zeeScales[i]*zXsec)*((zXsec)/(wpXsec+wmXsec)))); }

  Double_t wzmNom=(wpmNom*wpXsec+wmmNom*wmXsec)/(zmmNom*zXsec)*((zXsec)/(wpXsec+wmXsec)); 
  vector<Double_t> wzmScales;
  for (UInt_t i=0; i<wpmScales.size(); i++) { wzmScales.push_back(((wpmScales[i]*wpXsec+wmmScales[i]*wmXsec)/(zmmScales[i]*zXsec)*((zXsec)/(wpXsec+wmXsec)))); }

  Double_t wpmeNomAs=wpeNomAs/wmeNomAs; 
  vector<Double_t> wpmeScaleAs;
  for (UInt_t i=0; i<wmeScaleAs.size(); i++) { wpmeScaleAs.push_back(wpeScaleAs[i]/wmeScaleAs[i]); }

  Double_t wpmmNomAs=wpmNomAs/wmmNomAs; 
  vector<Double_t> wpmmScaleAs;
  for (UInt_t i=0; i<wmmScaleAs.size(); i++) { wpmmScaleAs.push_back(wpmScaleAs[i]/wmmScaleAs[i]); }

  Double_t weNomAs=(wpeNomAs*wpXsec+wmeNomAs*wmXsec)/(wpXsec+wmXsec); 
  vector<Double_t> weScaleAs;
  for (UInt_t i=0; i<wpeScaleAs.size(); i++) { weScaleAs.push_back((wpeScaleAs[i]*wpXsec+wmeScaleAs[i]*wmXsec)/(wpXsec+wmXsec)); }

  Double_t wmNomAs=(wpmNomAs*wpXsec+wmmNomAs*wmXsec)/(wpXsec+wmXsec); 
  vector<Double_t> wmScaleAs;
  for (UInt_t i=0; i<wpmScaleAs.size(); i++) { wmScaleAs.push_back((wpmScaleAs[i]*wpXsec+wmmScaleAs[i]*wmXsec)/(wpXsec+wmXsec)); }

  Double_t wzeNomAs=(wpeNomAs*wpXsec+wmeNomAs*wmXsec)/(zeeNomAs*zXsec)*((zXsec)/(wpXsec+wmXsec)); 
  vector<Double_t> wzeScaleAs;
  for (UInt_t i=0; i<wpeScaleAs.size(); i++) { wzeScaleAs.push_back(((wpeScaleAs[i]*wpXsec+wmeScaleAs[i]*wmXsec)/(zeeScaleAs[i]*zXsec)*((zXsec)/(wpXsec+wmXsec)))); }

  Double_t wzmNomAs=(wpmNomAs*wpXsec+wmmNomAs*wmXsec)/(zmmNomAs*zXsec)*((zXsec)/(wpXsec+wmXsec)); 
  vector<Double_t> wzmScaleAs;
  for (UInt_t i=0; i<wpmScaleAs.size(); i++) { wzmScaleAs.push_back(((wpmScaleAs[i]*wpXsec+wmmScaleAs[i]*wmXsec)/(zmmScaleAs[i]*zXsec)*((zXsec)/(wpXsec+wmXsec)))); }

  cout << "uncertainties" << endl;
  cout << "W+m:      " << setprecision(4) << uncert(wpmScales, wpmScaleAs, wpmNom, wpmNomAs) << " %" << endl;
  cout << "W-m:      " << setprecision(4) << uncert(wmmScales, wmmScaleAs, wmmNom, wmmNomAs) << " %" << endl;
  cout << "W+e:      " << setprecision(4) << uncert(wpeScales, wpeScaleAs, wpeNom, wpeNomAs) << " %" << endl;
  cout << "W-e:      " << setprecision(4) << uncert(wmeScales, wmeScaleAs, wmeNom, wmeNomAs) << " %" << endl;
  cout << "Zmm:      " << setprecision(4) << uncert(zmmScales, zmmScaleAs, zmmNom, zmmNomAs) << " %" << endl;
  cout << "Zee:      " << setprecision(4) << uncert(zeeScales, zeeScaleAs, zeeNom, zeeNomAs) << " %" << endl;

  cout << "W+/W-(m): " << setprecision(4) << uncert(wpmmScales, wpmmScaleAs, wpmmNom, wpmmNomAs) << " %" << endl;
  cout << "W+/W-(e): " << setprecision(4) << uncert(wpmeScales, wpmeScaleAs, wpmeNom, wpmeNomAs) << " %" << endl;
  cout << "W(m):     " << setprecision(4) << uncert(wmScales,   wmScaleAs,   wmNom,   wmNomAs) << " %" << endl;
  cout << "W(e):     " << setprecision(4) << uncert(weScales,   weScaleAs,   weNom,   weNomAs) << " %" << endl;
  cout << "W/Z(m):   " << setprecision(4) << uncert(wzmScales,  wzmScaleAs,  wzmNom,  wzmNomAs) << " %" << endl;
  cout << "W/Z(e):   " << setprecision(4) << uncert(wzeScales,  wzeScaleAs,  wzeNom,  wzeNomAs) << " %" << endl;

}
void makeScaleVector(TString inDir, TString chan, Double_t &nom, Double_t &nomAs, vector<Double_t> &vScales, vector<Double_t> &alphaS) {

  //cout << chan << endl;

  TString inFile = inDir+chan+"_CT10nlo.root";

  TFile *f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 52, vScales);
  delete f;

  nom=vScales[0];
  /*
  inFile = inDir+chan+"_CT10nlo_as_0115.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 15, alphaS);
  delete f;
  */
  inFile = inDir+chan+"_CT10nlo_as_0116.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 16, alphaS);
  delete f;
  /*
  inFile = inDir+chan+"_CT10nlo_as_0117.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 17, alphaS);
  delete f;
  */
  inFile = inDir+chan+"_CT10nlo_as_0118.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 18, alphaS);
  delete f;

  nomAs=alphaS[alphaS.size()-1];
  alphaS.pop_back();
  /*
  inFile = inDir+chan+"_CT10nlo_as_0119.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 19, alphaS);
  delete f;
  */
  inFile = inDir+chan+"_CT10nlo_as_0120.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 20, alphaS);
  delete f;
  /*
  inFile = inDir+chan+"_CT10nlo_as_0121.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 21, alphaS);
  delete f;

  inFile = inDir+chan+"_CT10nlo_as_0122.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 22, alphaS);
  delete f;

  inFile = inDir+chan+"_CT10nlo_as_0123.root";
  f = new TFile(inFile, "read");
  calcAcceptance(f, chan, 1, 23, alphaS);
  delete f;
  */

}

void calcAcceptance(TFile *f, TString chan, Int_t n, Int_t as, vector<Double_t> &scale) {
  char hname[100];
  
  for (Int_t i=0; i<n; i++) {
    //cout << "CT10nlo, as = 0.1" << as << ", i = " << i << endl;

    sprintf(hname, "dTot_CT10nlo_as_01%i_%i", as, i);
    //cout << hname << endl;
    TH1D *tot = (TH1D*) f->Get(hname);
    Double_t s = tot->Integral()/tot->GetEntries();

    Double_t accept;
    if (chan=="wmm" || chan=="wme" || chan=="wpe" || chan=="wpm") {
      sprintf(hname, "dPostB_CT10nlo_as_01%i_%i", as, i);
      TH1D *bar = (TH1D*) f->Get(hname);
      accept=bar->Integral()/(tot->Integral());
      //cout << "Barrel Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostE_CT10nlo_as_01%i_%i", as, i);
      TH1D *end = (TH1D*) f->Get(hname);
      accept=end->Integral()/(tot->Integral());
      //cout << "Endcap Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      accept=(bar->Integral()+end->Integral())/(tot->Integral());
    }
    else if (chan=="zmm" || chan=="zee") {
      sprintf(hname, "dPostBB_CT10nlo_as_01%i_%i", as, i);
      TH1D *bb = (TH1D*) f->Get(hname);
      accept=bb->Integral()/(tot->Integral());
      //cout << "Bar-bar Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostBE_CT10nlo_as_01%i_%i", as, i);
      TH1D *be = (TH1D*) f->Get(hname);
      accept=be->Integral()/(tot->Integral());
      //cout << "Bar-end Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostEE_CT10nlo_as_01%i_%i", as, i);
      TH1D *ee = (TH1D*) f->Get(hname);
      accept=ee->Integral()/(tot->Integral());
      //cout << "End-end Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      accept=(bb->Integral()+be->Integral()+ee->Integral())/(tot->Integral());
    }
    //cout << "Total Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;
    //cout << "Scale: " << tot->Integral()/tot->GetEntries() << endl;

    //scale.push_back(s);
    scale.push_back(accept);
  }
}

void calcAcceptance(TFile *f, TString chan, Int_t n, vector<Double_t> &scale) {
  char hname[100];
  
  for (Int_t i=0; i<n; i++) {
    if (i==0) cout << "wme ";

    sprintf(hname, "dTot_CT10nlo_%i", i);
    TH1D *tot = (TH1D*) f->Get(hname);
    Double_t s = tot->Integral()/tot->GetEntries();
    //Double_t s = tot->GetBinContent(1)/tot->GetEntries();

    Double_t accept;
    if (chan=="wmm" || chan=="wme" || chan=="wpe" || chan=="wpm") {
      sprintf(hname, "dPostB_CT10nlo_%i", i);
      TH1D *bar = (TH1D*) f->Get(hname);
      accept=bar->Integral()/(tot->Integral());
      //if (i==0) cout << "Barrel Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostE_CT10nlo_%i", i);
      TH1D *end = (TH1D*) f->Get(hname);
      accept=end->Integral()/(tot->Integral());
      //if (i==0) cout << "Endcap Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      accept=(bar->Integral()+end->Integral())/(tot->Integral());
    }
    else if (chan=="zmm" || chan=="zee") {
      sprintf(hname, "dPostBB_CT10nlo_%i", i);
      TH1D *bb = (TH1D*) f->Get(hname);
      accept=bb->Integral()/(tot->Integral());
      //if (i==0) cout << "Bar-bar Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostBE_CT10nlo_%i", i);
      TH1D *be = (TH1D*) f->Get(hname);
      accept=be->Integral()/(tot->Integral());
      //if (i==0) cout << "Bar-end Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostEE_CT10nlo_%i", i);
      TH1D *ee = (TH1D*) f->Get(hname);
      accept=ee->Integral()/(tot->Integral());
      //if (i==0) cout << "End-end Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      accept=(bb->Integral()+be->Integral()+ee->Integral())/(tot->Integral());
    }
    if (i==0) cout << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;
    //if (i==0) cout << "Scale: " << tot->Integral()/tot->GetEntries() << endl;

    //scale.push_back(s);
    scale.push_back(accept);
  }
}

Double_t uncert(vector<Double_t> &vScale, vector<Double_t> &vScaleAs, Double_t nom, Double_t nomAs) {

  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScale.size(); i+=2) {
    prevDiff=nom-vScale[i-1];
    thisDiff=nom-vScale[i];

    if (prevDiff*thisDiff>0) {
      if (fabs(prevDiff)>thisDiff) thisDiff=prevDiff;
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else negUnc+=thisDiff*thisDiff;
    }
    else {
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else if (thisDiff<0) negUnc+=thisDiff*thisDiff;

      if (prevDiff>0) posUnc+=prevDiff*prevDiff;
      if (prevDiff<0) negUnc+=prevDiff*prevDiff;
    }
  }

  posUnc/=(1.645*1.645);
  negUnc/=(1.645*1.645);

  Double_t negUncAs=0, posUncAs=0;
  prevDiff=0;
  for (UInt_t i=0; i<vScaleAs.size(); i++) {
    thisDiff=nomAs-vScaleAs[i];
    if (thisDiff>0) posUncAs+=thisDiff*thisDiff;
    else negUncAs+=thisDiff*thisDiff;
  }

  //posUncAs/=(25.0/36);
  posUncAs/=(1.645*1.645);
  //negUncAs/=(25.0/36);
  negUncAs/=(1.645*1.645);

  return max(sqrt(posUnc/nom/nom+posUncAs/nomAs/nomAs)*100,
             sqrt(negUnc/nom/nom+negUncAs/nomAs/nomAs)*100);

}
