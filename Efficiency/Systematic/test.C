#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "CPlot.hh"                 // helper class for plots
#include "MitStyleRemix.hh"         // style settings for drawing
#include "CEffUser1D.hh"            // class for handling efficiency graphs
#include "CEffUser2D.hh"            // class for handling efficiency tables

#include "ZSignals.hh"
#include "ZBackgrounds.hh"

// structure for output ntuple
#include "BinInfo.hh"

#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"

void test() {

  TCanvas *c1 = new TCanvas("c1", "", 800, 600);

  for (Int_t x=0; x<18; x++) {

    char filename[100];
    sprintf(filename, "CB_MuSelEff/analysis/plots/etapt_%i.root",x);

    TFile f(filename);

    if (! f.IsOpen() ) {
      cout << "file " << x << " doesn't exist" << endl;
      continue;
    }
    RooWorkspace* w=(RooWorkspace*)f.Get("w");
    RooRealVar* trueEff = w->var("eff");
    Float_t truthEff = trueEff->getVal();
    
    TTree *t = (TTree*) f.Get("Bin");
    BinInfo bin;
    t->SetBranchAddress("Bin",&bin);
    t->GetEntry(0);

    Int_t nevents = bin.nEvents;
    
    f.Close();
    
    cout << "truthEff " << truthEff << endl;
    
    TH1F* effDist = new TH1F("effDist", "", 500,-0.05,0.05);
    
    ifstream ifs;
    char txtfilename[50];
    sprintf(txtfilename, "CB_MuSelUncert/etapt_%i.dat", x);
    ifs.open(txtfilename);
    if (!ifs.is_open()) {
      cout << "file " << x << " doesn't exist" << endl;
      continue;
    }
    assert(ifs.is_open());
    string line;
    Float_t eff=0;
    Float_t diff=0;
    while(getline(ifs,line)) {
      stringstream ss(line);
      ss >> eff;
      diff = eff-truthEff;
      effDist->Fill(diff);

    }

    Float_t binUncert = effDist->GetMean();
    cout << binUncert << endl;

    effDist->GetXaxis()->SetTitle("FitEff-TruthEff");
    
    effDist->Draw();

    char outfilename[50];
    sprintf(outfilename, "CB_MuSelUncert/etapt_%i.png", x);
    c1->SaveAs(outfilename);

    delete(effDist);
    
  }

}
