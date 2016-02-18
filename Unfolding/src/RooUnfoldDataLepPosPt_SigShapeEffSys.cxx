#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>
#include <sstream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"


#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

#endif


//Include Headers
#include "/home/kfiekas/CMSSW_7_4_14/src/BootStrap/interface/BootStrap.hpp"

#ifdef __CINT__
gSystem->Load("/home/kfiekas/CMSSW_7_4_14/src/BootStrap/bin/libBootStrap.so")
#endif

string int2string(int i) {
  stringstream ss;
  string ret;
  ss << i;
  ss >> ret;
  return ret;
}

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;


TH1D *hReco, *hTruth, *hData, *hTop, *hEWK, *hMeas, *hUnfold;
TH2D *hMatrix, *hMatrix_hilf;

TGraphAsymmErrors *gUnfold;

//==============================================================================
// Train unfolding algorithm
//==============================================================================

Int_t Train ()
{
  cout <<"================ TRAIN ======================="<<endl;
  TFile *file1 = new TFile("../UnfoldingInput/Zmumu/zmm_UnfoldInputs.root");
  if (file1 == NULL) cout <<"File does not exists"<<endl;
  hTruth = (TH1D*)file1->Get("hLepPosPtTruth");
  hReco = (TH1D*)file1->Get("hLepPosPtReco_EffSigShape");
  hMatrix_hilf = (TH2D*)file1->Get("hLepPosPtMatrix_EffSigShape");
 
  if (hTruth==NULL) cout<< "hTruth does not exist"<<endl;
  if (hReco==NULL) cout<< "hReco does not exist"<<endl;
  if (hMatrix_hilf==NULL) cout<< "hMatrix does not exist"<<endl;

  hMatrix= (TH2D*)hMatrix_hilf->Clone("hMatrix");
  for(int j=0;j!=hMatrix_hilf->GetNbinsX();++j)
    {
      for(int k=0;k!=hMatrix_hilf->GetNbinsY();++k)
	{
	  hMatrix->SetBinContent(j+1,k+1,hMatrix_hilf->GetBinContent(k+1,j+1));
	  hMatrix->SetBinError(j+1,k+1,hMatrix_hilf->GetBinError(k+1,j+1));
	}
    }


  cout<<"Training finished"<<endl;
  return 1;
}

//==============================================================================
// Test distribution
//==============================================================================

Int_t Test ()
{
 cout <<"==================== TEST ===================="<<endl;
 TFile *f=new TFile("../SignalExtraction/Zmm/Zmm_DataBkg.root");
 
  if (f == NULL) cout<<"file does not exists"<<endl;
  hData = (TH1D*)f->Get("hDataLepPosPt");
  hTop = (TH1D*)f->Get("hTopLepPosPt_EffSigShape");
  hEWK = (TH1D*)f->Get("hEWKLepPosPt_EffSigShape");

  hMeas = (TH1D*)f->Get("hDataLepPosPt");
  hMeas->Sumw2();
  hMeas->Add(hTop,-1);
  hMeas->Add(hEWK,-1);


  cout<<"Test finished"<<endl;
  return 1;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfoldExample(const TString  outputDir,    // output directory
		      const Double_t lumi,         // integrated luminosity (/pb
		      const Int_t iterations)
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
 

  TH1::AddDirectory(kFALSE);

  if(!mkdir(outputDir,0755))
    cout<<"Created output directory"<<endl;

  double LUMI=lumi;
 
  Train();
  Test();
  
  cout << "==================================== UNFOLD ===================================" << endl;

  TFile* his=new TFile(outputDir+(string("/UnfoldingOutputLepPosPtEffSigShape.root")).c_str(), "recreate");

  gRandom->SetSeed(1234);

  RooUnfoldResponse response(hReco,hTruth,hMatrix);

  response.UseOverflow();


  RooUnfoldBayes   unfold (&response, hMeas, iterations);

  BootStrap b;
  b.SetUnfoldType(BootStrap::kBayes);
  b.SetRegParam(iterations);
  
  b.SetUMatrix(hReco,hTruth,hMatrix);
  b.SetData( (TH1D*)hMeas->Clone("bootstrap_data") );
  b.SetToyType(BootStrap::kBootstrap);
  b.run();

  TGraphAsymmErrors *g = b.result(BootStrap::kMin,.68);
  g->SetName("gBootStrap");

  unfold.SetNToys(1000);
  
  hUnfold=(TH1D*) unfold.Hreco(RooUnfold::kCovToy);
  
  for(int j=0;j!=hUnfold->GetNbinsX();++j)
    {
      hUnfold->SetBinContent(j+1,hUnfold->GetBinContent(j+1)/hUnfold->GetBinWidth(j+1));
      hUnfold->SetBinError(j+1,hUnfold->GetBinError(j+1)/hUnfold->GetBinWidth(j+1));
    }

  for(int j=0;j!=hTruth->GetNbinsX();++j)
    {
      hTruth->SetBinContent(j+1,hTruth->GetBinContent(j+1)/hTruth->GetBinWidth(j+1));
      hTruth->SetBinError(j+1,hTruth->GetBinError(j+1)/hTruth->GetBinWidth(j+1));
    }

  hUnfold->Scale(1./LUMI);
  hTruth->Scale(1./LUMI);

  his->cd();
  hTruth->Write("hTruth");
  hUnfold->Write("hUnfold");
  g->Write();
  his->Close();
}

#ifndef __CINT__
int main (int argc,char **argv) {RooUnfoldExample(argv[1],atof(argv[2]),atof(argv[3])); return 0; }  // Main program when run stand-alone
#endif
