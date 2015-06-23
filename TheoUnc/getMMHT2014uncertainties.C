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

void calcAcceptance(TFile *f, TString chan, Int_t n, vector<Double_t> &scale, TString sname);
void makeScaleVector(TString inDir, TString chan, Double_t &nom, vector<Double_t> &vScales, TString sname);
Double_t uncert(vector<Double_t> &vScale, Double_t nom);
Double_t uncert(vector<Double_t> &vScale, vector<Double_t> &vScaleAs);

void getMMHT2014uncertainties() {

  Double_t wmXsec = 8.35;
  Double_t wpXsec = 11.20;
  Double_t zXsec  = 1.92;

  TString inDir = "/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/MMHT2014_amc/";

  //cout << "MMHT2014 acceptances" << endl;

  Double_t wmeNom, wmeNomA;
  vector<Double_t> wmeScales, wmeScaleA;
  makeScaleVector(inDir, "wme", wmeNom, wmeScales, "MMHT2014nlo68cl");
  makeScaleVector(inDir, "wme", wmeNomA, wmeScaleA, "MMHT2014nlo_asmzsmallrange");

  Double_t wpeNom, wpeNomA;
  vector<Double_t> wpeScales, wpeScaleA;
  makeScaleVector(inDir, "wpe", wpeNom, wpeScales, "MMHT2014nlo68cl");
  makeScaleVector(inDir, "wpe", wpeNomA, wpeScaleA, "MMHT2014nlo_asmzsmallrange");

  Double_t wmmNom, wmmNomA; 
  vector<Double_t> wmmScales, wmmScaleA;
  makeScaleVector(inDir, "wmm", wmmNom, wmmScales, "MMHT2014nlo68cl");
  makeScaleVector(inDir, "wmm", wmmNomA, wmmScaleA, "MMHT2014nlo_asmzsmallrange");

  Double_t wpmNom, wpmNomA;
  vector<Double_t> wpmScales, wpmScaleA;
  makeScaleVector(inDir, "wpm", wpmNom, wpmScales, "MMHT2014nlo68cl");
  makeScaleVector(inDir, "wpm", wpmNomA, wpmScaleA, "MMHT2014nlo_asmzsmallrange");
  
  Double_t zeeNom, zeeNomA;
  vector<Double_t> zeeScales, zeeScaleA;
  makeScaleVector(inDir, "zee", zeeNom, zeeScales, "MMHT2014nlo68cl");
  makeScaleVector(inDir, "zee", zeeNomA, zeeScaleA, "MMHT2014nlo_asmzsmallrange");

  Double_t zmmNom, zmmNomA;
  vector<Double_t> zmmScales, zmmScaleA;
  makeScaleVector(inDir, "zmm", zmmNom, zmmScales, "MMHT2014nlo68cl");
  makeScaleVector(inDir, "zmm", zmmNomA, zmmScaleA, "MMHT2014nlo_asmzsmallrange");

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

  Double_t wpmeNomA=wpeNomA/wmeNomA; 
  vector<Double_t> wpmeScaleA;
  for (UInt_t i=0; i<wmeScaleA.size(); i++) { wpmeScaleA.push_back(wpeScaleA[i]/wmeScaleA[i]); }

  Double_t wpmmNomA=wpmNomA/wmmNomA; 
  vector<Double_t> wpmmScaleA;
  for (UInt_t i=0; i<wmmScaleA.size(); i++) { wpmmScaleA.push_back(wpmScaleA[i]/wmmScaleA[i]); }

  Double_t weNomA=(wpeNomA*wpXsec+wmeNomA*wmXsec)/(wpXsec+wmXsec); 
  vector<Double_t> weScaleA;
  for (UInt_t i=0; i<wpeScaleA.size(); i++) { weScaleA.push_back((wpeScaleA[i]*wpXsec+wmeScaleA[i]*wmXsec)/(wpXsec+wmXsec)); }

  Double_t wmNomA=(wpmNomA*wpXsec+wmmNomA*wmXsec)/(wpXsec+wmXsec); 
  vector<Double_t> wmScaleA;
  for (UInt_t i=0; i<wpmScaleA.size(); i++) { wmScaleA.push_back((wpmScaleA[i]*wpXsec+wmmScaleA[i]*wmXsec)/(wpXsec+wmXsec)); }

  Double_t wzeNomA=(wpeNomA*wpXsec+wmeNomA*wmXsec)/(zeeNomA*zXsec)*((zXsec)/(wpXsec+wmXsec)); 
  vector<Double_t> wzeScaleA;
  for (UInt_t i=0; i<wpeScaleA.size(); i++) { wzeScaleA.push_back(((wpeScaleA[i]*wpXsec+wmeScaleA[i]*wmXsec)/(zeeScaleA[i]*zXsec)*((zXsec)/(wpXsec+wmXsec)))); }

  Double_t wzmNomA=(wpmNomA*wpXsec+wmmNomA*wmXsec)/(zmmNomA*zXsec)*((zXsec)/(wpXsec+wmXsec)); 
  vector<Double_t> wzmScaleA;
  for (UInt_t i=0; i<wpmScaleA.size(); i++) { wzmScaleA.push_back(((wpmScaleA[i]*wpXsec+wmmScaleA[i]*wmXsec)/(zmmScaleA[i]*zXsec)*((zXsec)/(wpXsec+wmXsec)))); }

  cout << "PDF: MMHT2014 " << endl;
  cout << "W+m:      " << setprecision(4) <<  uncert(wpmScales,   wpmScaleA) << "\\\% " << endl;
  cout << "W-m:      " << setprecision(4) <<  uncert(wmmScales,   wmmScaleA) << "\\\% " << endl;
  cout << "W(m):     " << setprecision(4) <<  uncert(wmScales,     wmScaleA) << "\\\% " << endl;
  cout << "W+/W-(m): " << setprecision(4) <<  uncert(wpmmScales, wpmmScaleA) << "\\\% " << endl;
  cout << "Zmm:      " << setprecision(4) <<  uncert(zmmScales,   zmmScaleA) << "\\\% " << endl;
  cout << "W/Z(m):   " << setprecision(4) <<  uncert(wzmScales,   wzmScaleA) << "\\\% " << endl;

  cout << "W+e:      " << setprecision(4) <<  uncert(wpeScales,   wpeScaleA) << "\\\% " << endl;
  cout << "W-e:      " << setprecision(4) <<  uncert(wmeScales,   wmeScaleA) << "\\\% " << endl;
  cout << "W(e):     " << setprecision(4) <<  uncert(weScales,     weScaleA) << "\\\% " << endl;
  cout << "W+/W-(e): " << setprecision(4) <<  uncert(wpmeScales, wpmeScaleA) << "\\\% " << endl;
  cout << "Zee:      " << setprecision(4) <<  uncert(zeeScales,   zeeScaleA) << "\\\% " << endl;
  cout << "W/Z(e):   " << setprecision(4) <<  uncert(wzeScales,   wzeScaleA) << "\\\% " << endl;

}
void makeScaleVector(TString inDir, TString chan, Double_t &nom, vector<Double_t> &vScales, TString sname) {

  TString inFile = inDir+chan+"_"+sname+".root";

  TFile *f = new TFile(inFile, "read");
  if (sname == "MMHT2014nlo68cl") calcAcceptance(f, chan, 51, vScales, sname);
  else calcAcceptance(f, chan, 5, vScales, sname);
  delete f;

  nom=vScales[0];

}

void calcAcceptance(TFile *f, TString chan, Int_t n, vector<Double_t> &scale, TString sname) {
  char hname[100];
  
  for (Int_t i=0; i<n; i++) {
    Bool_t prnt = kTRUE;
    if (sname!="MMHT2014nlo68cl") prnt=kFALSE;
    if (sname=="MMHT2014nlo68cl" && i==0) cout << "Nominal ";
    if (prnt) cout << "A (" << chan << "): ";

    sprintf(hname, "dTot_%s_%i", sname.Data(), i);
    TH1D *tot = (TH1D*) f->Get(hname);

    Double_t accept;
    if (chan=="wmm" || chan=="wme" || chan=="wpe" || chan=="wpm") {
      sprintf(hname, "dPostB_%s_%i", sname.Data(), i);
      TH1D *bar = (TH1D*) f->Get(hname);
      accept=bar->Integral()/(tot->Integral());
      //if (prnt) cout << "Barrel Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostE_%s_%i", sname.Data(), i);
      TH1D *end = (TH1D*) f->Get(hname);
      accept=end->Integral()/(tot->Integral());
      //if (prnt) cout << "Endcap Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;
      accept=(bar->Integral()+end->Integral())/(tot->Integral());
    }
    else if (chan=="zmm" || chan=="zee") {
      sprintf(hname, "dPostBB_%s_%i", sname.Data(), i);
      TH1D *bb = (TH1D*) f->Get(hname);
      accept=bb->Integral()/(tot->Integral());
      //if (prnt) cout << "Bar-bar Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostBE_%s_%i", sname.Data(), i);
      TH1D *be = (TH1D*) f->Get(hname);
      accept=be->Integral()/(tot->Integral());
      //if (prnt) cout << "Bar-end Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

      sprintf(hname, "dPostEE_%s_%i", sname.Data(), i);
      TH1D *ee = (TH1D*) f->Get(hname);
      accept=ee->Integral()/(tot->Integral());
      //if (prnt)cout << "End-end Acceptance: " << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;
      accept=(bb->Integral()+be->Integral()+ee->Integral())/(tot->Integral());
    }
    if (prnt) cout << accept << " +/- " << sqrt(accept*(1-accept)/tot->GetEntries()) << endl;

    scale.push_back(accept);
  }
}

Double_t uncert(vector<Double_t> &vScale, Double_t nom) {

  Double_t nomScale=vScale[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScale.size(); i+=2) {
    prevDiff=nomScale-vScale[i-1];
    thisDiff=nomScale-vScale[i];

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

  Double_t centralDiff=fabs(nom-nomScale);

  return max((sqrt(posUnc)+centralDiff)/nomScale*100,
	     (sqrt(negUnc)+centralDiff)/nomScale*100);

}

Double_t uncert(vector<Double_t> &vScale, vector<Double_t> &vScaleAs) {

  Double_t nomScale=vScale[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScale.size(); i+=2) {
    prevDiff=nomScale-vScale[i-1];
    thisDiff=nomScale-vScale[i];

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

  Double_t posUncAs=abs(vScaleAs[0]-vScaleAs[3])*(vScaleAs[0]-vScaleAs[3]);
  Double_t negUncAs=(vScaleAs[0]-vScaleAs[4])*(vScaleAs[0]-vScaleAs[4]);

  return max((sqrt(posUnc+posUncAs))/nomScale*100,
	     (sqrt(negUnc+negUncAs))/nomScale*100);

}
