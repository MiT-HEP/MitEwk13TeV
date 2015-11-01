#ifndef COMPARE_HH
#define COMPARE_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TGraph.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "../Utils/MitStyleRemix.hh"

using namespace std;

class Config {
public:
  Config(TString fname, TString sname, TString lab) {
    filename   = fname;
    samplename = sname;
    label = lab;
    file = new TFile(fname, "READ");
  };
  ~Config() {};

  TString getSampleName() { return samplename; };
  TString getLabel() {return label; };

  TFile *file;

protected:
  TString filename;
  TString samplename;
  TString label;
};

enum Channel { wmm=0, wme, wpm, wpe, zmm, zee};

TString chan_name[6]= {"wmm", "wme", "wpm", "wpe", "zmm", "zee"};
void drawFourPt(Channel chan, TString var, Config fB, Config f1, Config f2, Config f3, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt);
void drawThreePt(Channel chan, TString var, Config fB, Config f1, Config f2, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt);
void drawTwoPt(Channel chan, TString var, Config fB, Config f1, Int_t nbinsPt, Int_t xminPt, Int_t xmaxPt);
void drawThreeEta(Channel chan, TString var, Config fB, Config f1, Config f2, Int_t nbinsEta, Int_t xminEta, Int_t xmaxEta);
void drawTwoEta(Channel chan, TString var, Config fB, Config f1, Int_t nbinsEta, Int_t xminEta, Int_t xmaxEta);
Double_t calcAcc(Channel chan, Config f);
TH1D* returnPlot(Channel chan, Config f, Int_t nbins, Double_t xmin, Double_t xmax, TString var);
TH1D* returnRelDiff(TH1D* h, TH1D* b, TString name);
TH1D* returnRelDiff(TH1D* h, TGraph* b, TString name);

#endif
