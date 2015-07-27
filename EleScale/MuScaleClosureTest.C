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

/*
  Data->MC scale correction
$0 < |\eta| < 0.8$ & $0.998766$ \pm $0.00123151$ \\
$0.8 < |\eta| < 1.2$ & $0.998886$ \pm $0.00157559$ \\
$1.6 < |\eta| < 2.4$ & $0.999522$ \pm $0.00209911$ \\

  MC->Data resolution correction [GeV]
0 < |\eta| < 0.8 & $0.668464$ \pm $0.169995$ \\
0.8 < |\eta| < 1.2 & $0.518957$ \pm $0.276301$ \\
1.6 < |\eta| < 2.4 & $0.869929$ \pm $0.244774$ \\
*/

//=== MAIN MACRO ================================================================================================= 

void MuScaleClosureTest() {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // event category enumeration
  enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum                                                                                      
  
  TString outputDir = "muon";
  TString pufname = "../Tools/pileup_weights_2015B.root";
  
  vector<TString> infilenamev;
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat_07_23/Zmumu/ntuples/data_select.root");  // data
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat_07_23/Zmumu/ntuples/zmm_select.root");  // MC
  infilenamev.push_back("/data/blue/Bacon/Run2/wz_flat_07_23/Zmumu/ntuples/zmm_select.root");  // MC2
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;
  const Double_t PT_CUT    = 25;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MU_MASS   = 0.105658369;
  
  vector<pair<Double_t,Double_t> > scEta_limits;
  scEta_limits.push_back(make_pair(0.0,0.8));
  scEta_limits.push_back(make_pair(0.8,1.6));
  scEta_limits.push_back(make_pair(1.6,2.4));

  CPlot::sOutDir = outputDir;
  
  const TString format("png");

  TRandom3 *rnd = new TRandom3();  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  enum { eData=0, eMC, eMC2 };
  
  TFile *pufile = new TFile(pufname); assert(pufile);
  TH1D  *puWeights = (TH1D*)pufile->Get("npv_rw");

  TH1D* hMC_Tot = new TH1D("hMC_Tot", "", 30, MASS_LOW, MASS_HIGH);
  TH1D* hData_Tot = new TH1D("hData_Tot", "", 30, MASS_LOW, MASS_HIGH);
  TH1D* hData2_Tot = new TH1D("hData2_Tot", "", 30, MASS_LOW, MASS_HIGH);
  
  char hname[100];
  vector<TH1D*> hMCv, hDatav, hDatav2;  
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=ibin; jbin<scEta_limits.size(); jbin++) {
      sprintf(hname,"mc_%i_%i",ibin,jbin);
      hMCv.push_back(new TH1D(hname,"",30,MASS_LOW,MASS_HIGH));
      hMCv.back()->Sumw2();
      
      sprintf(hname,"data_%i_%i",ibin,jbin);
      hDatav.push_back(new TH1D(hname,"",30,MASS_LOW,MASS_HIGH));
      hDatav.back()->Sumw2();

      sprintf(hname,"data2_%i_%i",ibin,jbin);
      hDatav2.push_back(new TH1D(hname,"",30,MASS_LOW,MASS_HIGH));
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
  
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    cout << "Processing " << infilenamev[ifile] << "..." << endl;
    TFile *infile = TFile::Open(infilenamev[ifile]); assert(infile);
    TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  
    intree->SetBranchAddress("runNum",   &runNum);    // event run number
    intree->SetBranchAddress("lumiSec",  &lumiSec);   // event lumi section
    intree->SetBranchAddress("evtNum",   &evtNum);    // event number
    intree->SetBranchAddress("scale1fb", &scale1fb);  // event weight
    intree->SetBranchAddress("puWeight", &puWeight);  // event weight
    intree->SetBranchAddress("matchGen", &matchGen);  // event has both leptons matched to MC Z->ll
    intree->SetBranchAddress("category", &category);  // dilepton category
    intree->SetBranchAddress("npv",      &npv);	      // number of primary vertices
    intree->SetBranchAddress("npu",      &npu);	      // number of in-time PU events (MC)
    intree->SetBranchAddress("q1",       &q1);	      // charge of lead lepton
    intree->SetBranchAddress("q2",       &q2);	      // charge of trail lepton
    intree->SetBranchAddress("dilep",    &dilep);     // dilepton 4-vector
    intree->SetBranchAddress("lep1",     &lep1);      // lead lepton 4-vector
    intree->SetBranchAddress("lep2",     &lep2);      // trail lepton 4-vector
  
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      Double_t weight = 1;
      if(ifile==eMC || ifile==eMC2) {
	//if(!matchGen) continue;
        weight=scale1fb*7.3*puWeight;
      }
      
      if((category!=eMuMu2HLT) && (category!=eMuMu1HLT) && (category!=eMuMu1HLT1L1)) continue;
      if(q1 == q2) continue;
      if(dilep->M()	  < MASS_LOW)  continue;
      if(dilep->M()	  > MASS_HIGH) continue;
      if(lep1->Pt()	  < PT_CUT)    continue;
      if(lep2->Pt()	  < PT_CUT)    continue;
      if(fabs(lep1->Eta()) > ETA_CUT)   continue;      
      if(fabs(lep2->Eta()) > ETA_CUT)   continue;

      TLorentzVector vLep1(0,0,0,0); 
      TLorentzVector vLep2(0,0,0,0); 
      if (ifile==eData) {
	vLep1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), MU_MASS);
	vLep2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), MU_MASS);
      }
      else if (ifile==eMC) {
	vLep1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), MU_MASS);
	vLep2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), MU_MASS);
      }
      else {
	vLep1.SetPtEtaPhiM(gRandom->Gaus(lep1->Pt()*getMuScaleCorr(lep1->Eta(),0),getMuResCorr(lep1->Eta(),0)), lep1->Eta(), lep1->Phi(), MU_MASS);
	vLep2.SetPtEtaPhiM(gRandom->Gaus(lep2->Pt()*getMuScaleCorr(lep2->Eta(),0),getMuResCorr(lep2->Eta(),0)), lep2->Eta(), lep2->Phi(), MU_MASS);
      }
      TLorentzVector vDilep = vLep1 + vLep2;
    
      Int_t bin1=-1, bin2=-1;
      for(UInt_t i=0; i<scEta_limits.size(); i++) {
        Double_t etalow  = scEta_limits.at(i).first;
        Double_t etahigh = scEta_limits.at(i).second;
        if(fabs(lep1->Eta())>=etalow && fabs(lep1->Eta())<=etahigh) bin1=i;
        if(fabs(lep2->Eta())>=etalow && fabs(lep2->Eta())<=etahigh) bin2=i;
      }
      assert(bin1>=0);
      assert(bin2>=0);
      Int_t ibin= (bin1<=bin2) ? bin1 : bin2;
      Int_t jbin= (bin1<=bin2) ? bin2 : bin1;

      if (ifile==eData) hData_Tot->Fill(vDilep.M(),1);
      else if (ifile==eMC) hMC_Tot->Fill(vDilep.M(),weight);
      else if (ifile==eMC2) hData2_Tot->Fill(vDilep.M(),weight);
     
      UInt_t n=jbin-ibin;
      for(Int_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);

      if(ifile==eData) hDatav[n]->Fill(vDilep.M(),1);
      else if(ifile==eMC) hMCv[n]->Fill(vDilep.M(),weight);
      else if(ifile==eMC2) hDatav2[n]->Fill(vDilep.M(),weight);
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

      hMCv[n]->SetLineColor(kRed);
      hMCv[n]->GetXaxis()->SetTitle("m_{ee}");
      hMCv[n]->GetYaxis()->SetRangeUser(0, 1.2*hMCv[n]->GetMaximum());
      hMCv[n]->Draw("hist");
      hDatav[n]->Draw("epsame");
      hDatav2[n]->SetLineColor(kGreen+1);
      hDatav2[n]->Draw("histsame");

      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry(hMCv[n],"uncorrected MC","l");
      leg->AddEntry(hDatav[n],"Data","l");
      leg->AddEntry(hDatav2[n],"corrected MC","l");
      leg->Draw();

      sprintf(pname,"comp_%i_%i.png",ibin,jbin); 
      c1->SaveAs(outputDir+"/"+pname);
    }
  }

  cout << endl;
  cout << hMC_Tot->Integral() << ", " << hData_Tot->Integral() << ", " << hData2_Tot->Integral() << endl;

  hMC_Tot->SetLineColor(kRed);
  hMC_Tot->GetXaxis()->SetTitle("m_{#mu#mu}");
  hMC_Tot->GetYaxis()->SetRangeUser(0, 1.2*hMC_Tot->GetMaximum());
  hMC_Tot->Draw("hist");
  hData_Tot->Draw("epsame");
  hData2_Tot->SetLineColor(kGreen+1);
  hData2_Tot->Draw("histsame");

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetShadowColor(0);
  leg->AddEntry(hMC_Tot,"uncorrected MC","l");
  leg->AddEntry(hData_Tot,"Data","l");
  leg->AddEntry(hData2_Tot,"corrected MC","l");
  leg->Draw();

  sprintf(pname,"comp_tot.png");
  c1->SaveAs(outputDir+"/"+pname);

}
