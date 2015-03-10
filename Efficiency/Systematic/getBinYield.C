//================================================================================================
//
// Perform fit to extract W->munu signal
//
//  * outputs plots and fit results summary
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "Math/LorentzVector.h"           // 4-vector class

#include "../../Utils/MyTools.hh"	          // various helper functions
#include "../../Utils/CPlot.hh"	          // helper class for plots
#include "../../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET

#include "BinInfo.hh"

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

//=== MAIN MACRO ================================================================================================= 

void getBinYield(const TString effType="BK_EleEff",
		 const Bool_t absEta=kFALSE,
		 const Int_t nbins=10) {

  gBenchmark->Start("fitWm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // Selection cuts
  const Double_t METMAX  = 100;
  const Double_t PT_CUT  = 30;
  const Double_t ETA_CUT = 2.1;

  // Set up binning from efficiency

  vector<Float_t> etaLoV;
  vector<Float_t> etaHiV;
  vector<Float_t> ptLoV;
  vector<Float_t> ptHiV;

  Int_t binYields[nbins];
  Float_t truthEff[nbins];
  Float_t binUncert[nbins];
  //Bool_t absEta = kFALSE;

  for (Int_t ibin=0; ibin<nbins; ibin++) {

    binYields[ibin]= 0;

    char binfile[50];
    sprintf(binfile, "%s/analysis/plots/etapt_%i.root",effType.Data(),ibin);

    TFile *f = new TFile(binfile);  
    TTree *intree = (TTree*)f->Get("Bin");
    BinInfo bin;
    intree->SetBranchAddress("Bin",&bin);
    intree->GetEntry(0);

    etaLoV.push_back(bin.etaLo);
    etaHiV.push_back(bin.etaHi);
    ptLoV.push_back(bin.ptLo);
    ptHiV.push_back(bin.ptHi);

    RooWorkspace* w=(RooWorkspace*)f->Get("w");
    RooRealVar* trueEff = w->var("eff");
    truthEff[ibin] = trueEff->getVal();

    delete(intree);
    delete(f);
  }

  //
  // Declare variables to read in ntuple
  //
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi;
  Float_t scale1fb;
  Float_t met, metPhi, sumEt, mt, u1, u2;
  Int_t   q;
  LorentzVector *lep=0;
  Float_t pfChIso, pfGamIso, pfNeuIso;
    
  TFile *infile=0;
  TTree *intree=0;

  // Read input file and get the TTrees
  infile = new TFile("/scratch/klawhorn/EWKAnaR12a/Selection/Wmunu/ntuples/data_select.root");	  assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);
  intree->SetBranchAddress("met",      &met);       // MET
  intree->SetBranchAddress("metPhi",   &metPhi);    // phi(MET)
  intree->SetBranchAddress("runNum",   &runNum);    // event run number
  intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
  intree->SetBranchAddress("evtNum",   &evtNum);    // event number
  intree->SetBranchAddress("npv",      &npv);       // number of primary vertices
  intree->SetBranchAddress("npu",      &npu);       // number of in-time PU events (MC)
  intree->SetBranchAddress("genVPt",   &genVPt);    // GEN W boson pT (signal MC)
  intree->SetBranchAddress("genVPhi",  &genVPhi);   // GEN W boson phi (signal MC)   
  intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight per 1/fb (MC)
  intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
  intree->SetBranchAddress("mt",       &mt);        // transverse mass
  intree->SetBranchAddress("u1",       &u1);        // parallel component of recoil
  intree->SetBranchAddress("u2",       &u2);        // perpendicular component of recoil
  intree->SetBranchAddress("q",        &q);	      // lepton charge
  intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
  intree->SetBranchAddress("pfChIso",  &pfChIso);   // nearby charged hadrons
  intree->SetBranchAddress("pfGamIso", &pfGamIso);  // nearby photons
  intree->SetBranchAddress("pfNeuIso", &pfNeuIso);  // nearby neutral hadrons
  
  //
  // loop over events
  //
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    if(lep->Pt()        < PT_CUT)  continue;	
    if(fabs(lep->Eta()) > ETA_CUT) continue;
    
    for (Int_t x=0; x<nbins; x++) {

      if (absEta==kFALSE) {
	//vcout << "absEta false " << etaLoV[x] << " " << etaHiV[x] << endl;
	if ( (lep->Pt() > ptLoV[x]) && (lep->Pt() < ptHiV[x]) && (lep->Eta() > etaLoV[x]) && (lep->Eta() < etaHiV[x]) && (met > 0) && (met < METMAX) ) {
	  binYields[x]++;
	  continue;
	}
      }
      else {
	if ( (lep->Pt() > ptLoV[x]) && (lep->Pt() < ptHiV[x]) && (fabs(lep->Eta()) > etaLoV[x]) && (fabs(lep->Eta()) < etaHiV[x]) && (met > 0) && (met < METMAX) ) {
	  binYields[x]++;
	  continue;
	}
      }
    }
  }
  
  delete infile;
  infile=0, intree=0;   

  // count things!
  Int_t totalYield=0;
  Float_t totalUncert=0;

  // get bin uncertainties

  for (Int_t x=0; x<nbins; x++) {
    
    TH1F* effDist = new TH1F("effDist", "", 500,-0.05,0.05);
    
    ifstream ifs;
    char txtfilename[50];
    sprintf(txtfilename, "%s_BinUncert/etapt_%i.dat", effType.Data(),x);
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
      diff = eff-truthEff[x];
      effDist->Fill(diff);

    }
    cout << "Bin " << x << " is " << etaLoV[x] << " < eta < " << etaHiV[x] << " and " << ptLoV[x] << " < pt < " << ptHiV[x] << endl;
    cout << "Bin yield for bin " << x << " is " << binYields[x] << endl;
    totalYield+=binYields[x];

    binUncert[x] = effDist->GetMean();
    cout << "Bin uncertainty is " << binUncert[x] << endl;
    cout << "Bin uncertainty time bin yield is " << binUncert[x]*binYields[x] << endl;

    totalUncert+=( binUncert[x] * binUncert[x] * binYields[x] * binYields[x] );

    delete(effDist);

  }

  cout << "---------SUMMARY----------" << endl;
  cout <<  " total yield is " << totalYield << endl;
  cout << effType << " overall uncertainty is " << TMath::Sqrt(totalUncert) << endl;
  cout << endl;
  cout << " relative uncertainty is " << TMath::Sqrt(totalUncert) * 100 / totalYield << "%" << endl;

}
