#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting style
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TCanvas.h>                  // class for drawing
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
/*
  Data->MC scale correction
$0 < |\eta| < 1.4442$ & $0.997542$ \pm $0.0011134$ \\
$1.566 < |\eta| < 2.5$ & $1.01507$ \pm $0.00126721$ \\

  MC->Data resolution correction [GeV]
0 < |\eta| < 1.4442 & $0.959767$ \pm $0.126489$ \\
1.566 < |\eta| < 2.5 & $1.28455$ \pm $0.277164$ \\
*/

Double_t getScaleCorr(const Double_t eta)
{
  if(fabs(eta) < 1.4442) { return 1.0/0.997542; }
  else                   { return 1.0/1.01507;  }
}

Double_t getResCorr(const Double_t eta)
{
  if(fabs(eta) < 1.4442) { return 0.959767; }
  else                   { return 1.28455;  }
}


//=== MAIN MACRO ================================================================================================= 

void EleScaleClosureTest() {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // event category enumeration
  enum { eEleEle2HLT=1, eEleEle1HLT1L1, eEleEle1HLT, eEleEleNoSel, eEleSC };
  
  TString outputDir = "Results";
  TString pufname = ""; 
  
  vector<TString> infilenamev;
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat/Zee/ntuples/data_select.raw.root");  // data
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat/Zee/ntuples/zee_select.root");  // MC
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat/Zee/ntuples/zee_select.root");  // MC2
  
  const Double_t MASS_LOW  = 80;
  const Double_t MASS_HIGH = 100;
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS  = 0.000511;  
  
  vector<pair<Double_t,Double_t> > scEta_limits;
  scEta_limits.push_back(make_pair(0.0,1.4442));
  scEta_limits.push_back(make_pair(1.566,2.5));

  CPlot::sOutDir = outputDir;
  
  const TString format("png");

  TRandom3 *rnd = new TRandom3();  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  enum { eData=0, eMC, eMC2 };
  
//  TFile *pufile = new TFile(pufname); assert(pufile);
//  TH1D  *puWeights = (TH1D*)pufile->Get("puWeights");
  
  char hname[100];
  vector<TH1D*> hMCv, hDatav, hDatav2;  
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      sprintf(hname,"mc_%i_%i",ibin,jbin);
      hMCv.push_back(new TH1D(hname,"",20,MASS_LOW,MASS_HIGH));
      hMCv.back()->Sumw2();
      
      sprintf(hname,"data_%i_%i",ibin,jbin);
      hDatav.push_back(new TH1D(hname,"",20,MASS_LOW,MASS_HIGH));
      hDatav.back()->Sumw2();

      sprintf(hname,"data2_%i_%i",ibin,jbin);
      hDatav2.push_back(new TH1D(hname,"",20,MASS_LOW,MASS_HIGH));
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
  Float_t scale1fb;
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
        weight=scale1fb;
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
	vLep1.SetPtEtaPhiM(gRandom->Gaus(lep1->Pt()*getScaleCorr(lep1->Eta()),getResCorr(lep1->Eta())), lep1->Eta(), lep1->Phi(), ELE_MASS);
	vLep2.SetPtEtaPhiM(gRandom->Gaus(lep2->Pt()*getScaleCorr(lep2->Eta()),getResCorr(lep2->Eta())), lep2->Eta(), lep2->Phi(), ELE_MASS);
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
      
      UInt_t n=jbin-ibin;
      for(Int_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);
      
      if(ifile==eData) hDatav[n]->Fill(vDilep.M(),weight);
      if(ifile==eMC)   hMCv[n]->Fill(vDilep.M(),weight);
      if(ifile==eMC2)   hDatav2[n]->Fill(vDilep.M(),weight);
    }  
    delete infile;
    infile=0, intree=0;
  }

  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  char pname[100];

  for (UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      UInt_t n=jbin-ibin;
      for(UInt_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);

      hMCv[n]   ->Scale(1.0/hMCv[n]->Integral());
      hDatav[n] ->Scale(1.0/hDatav[n]->Integral());
      hDatav2[n]->Scale(1.0/hDatav2[n]->Integral());

      hMCv[n]->SetLineColor(kRed);
      hMCv[n]->GetXaxis()->SetTitle("m_{ee}");
      hMCv[n]->GetYaxis()->SetRangeUser(0, 1.2*hMCv[n]->GetMaximum());
      hMCv[n]->Draw("hist");
      hDatav[n]->Draw("histsame");
      hDatav2[n]->SetLineColor(kGreen+1);
      hDatav2[n]->Draw("histsame");

      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry(hMCv[n],"MC, pre-correction","l");
      leg->AddEntry(hDatav[n],"Data","l");
      leg->AddEntry(hDatav2[n],"MC, post-correction","l");
      leg->Draw();

      sprintf(pname,"comp_%i_%i.png",ibin,jbin); 
      c1->SaveAs(outputDir+"/"+pname);
    }
  }

}
