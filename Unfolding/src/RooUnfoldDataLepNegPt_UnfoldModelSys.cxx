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
TH1D *hReco_Nominal, *hTruth_Nominal, *hData_Nominal, *hTop_Nominal, *hEWK_Nominal, *hMeas_Nominal, *hUnfold_Nominal;
TH2D *hMatrix_Nominal, *hMatrix_hilf_Nominal;

TGraphAsymmErrors *gUnfold;

//==============================================================================
// Train unfolding algorithm
//==============================================================================

Int_t TrainNominal ()
{
  cout <<"================ TRAIN ======================="<<endl;
  TFile *file1 = new TFile("../UnfoldingInput/Zmumu/zmm_UnfoldInputs.root");
  if (file1 == NULL) cout <<"File does not exists"<<endl;
  hTruth_Nominal = (TH1D*)file1->Get("hLepNegPtTruth");
  hReco_Nominal = (TH1D*)file1->Get("hLepNegPtReco");
  hMatrix_hilf_Nominal = (TH2D*)file1->Get("hLepNegPtMatrix");
 
  if (hTruth_Nominal==NULL) cout<< "hTruth does not exist"<<endl;
  if (hReco_Nominal==NULL) cout<< "hReco does not exist"<<endl;
  if (hMatrix_hilf_Nominal==NULL) cout<< "hMatrix does not exist"<<endl;

  hMatrix_Nominal= (TH2D*)hMatrix_hilf_Nominal->Clone("hMatrix");
  for(int j=0;j!=hMatrix_hilf_Nominal->GetNbinsX();++j)
    {
      for(int k=0;k!=hMatrix_hilf_Nominal->GetNbinsY();++k)
	{
	  hMatrix_Nominal->SetBinContent(j+1,k+1,hMatrix_hilf_Nominal->GetBinContent(k+1,j+1));
	  hMatrix_Nominal->SetBinError(j+1,k+1,hMatrix_hilf_Nominal->GetBinError(k+1,j+1));
	}
    }


  cout<<"Training finished"<<endl;
  return 1;
}

Int_t Train ()
{
  cout <<"================ TRAIN ======================="<<endl;
  TFile *file1 = new TFile("../UnfoldingInput/Zmumu/zmmmg_UnfoldInputs.root");
  if (file1 == NULL) cout <<"File does not exists"<<endl;
  hTruth = (TH1D*)file1->Get("hLepNegPtTruth");
  hReco = (TH1D*)file1->Get("hLepNegPtReco");
  hMatrix_hilf = (TH2D*)file1->Get("hLepNegPtMatrix");
 
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

Int_t TestNominal ()
{
 cout <<"==================== TEST ===================="<<endl;
 TFile *f=new TFile("../SignalExtraction/Zmumu/Zmm_DataBkg.root");
 
  if (f == NULL) cout<<"file does not exists"<<endl;
  hData_Nominal = (TH1D*)f->Get("hDataLepNegPt");
  hTop_Nominal = (TH1D*)f->Get("hTopLepNegPt");
  hEWK_Nominal = (TH1D*)f->Get("hEWKLepNegPt");

  hMeas_Nominal = (TH1D*)f->Get("hDataLepNegPt");
  hMeas_Nominal->Sumw2();
  hMeas_Nominal->Add(hTop_Nominal,-1);
  hMeas_Nominal->Add(hEWK_Nominal,-1);


  cout<<"Test finished"<<endl;
  return 1;
}

Int_t Test ()
{
 cout <<"==================== TEST ===================="<<endl;
 TFile *f=new TFile("../SignalExtraction/Zmumu/Zmm_DataBkg.root");
 
  if (f == NULL) cout<<"file does not exists"<<endl;
  hData = (TH1D*)f->Get("hDataLepNegPt");
  hTop = (TH1D*)f->Get("hTopLepNegPt");
  hEWK = (TH1D*)f->Get("hEWKLepNegPt");

  hMeas = (TH1D*)f->Get("hDataLepNegPt");
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
 
  
  TrainNominal();
  TestNominal();

  
  cout << "==================================== UNFOLD ===================================" << endl;

  TFile* his=new TFile(outputDir+(string("/UnfoldingOutputLepNegPtUnfoldModel.root")).c_str(), "recreate");

  gRandom->SetSeed(1234);

  RooUnfoldResponse responseNominal(hReco_Nominal,hTruth_Nominal,hMatrix_Nominal);

  responseNominal.UseOverflow();


  RooUnfoldBayes   unfoldNominal (&responseNominal, hMeas_Nominal, iterations);

  BootStrap bNominal;
  bNominal.SetUnfoldType(BootStrap::kBayes);
  bNominal.SetRegParam(iterations);
  
  bNominal.SetUMatrix(hReco_Nominal,hTruth_Nominal,hMatrix_Nominal);
  bNominal.SetData( (TH1D*)hMeas_Nominal->Clone("bootstrap_data") );
  bNominal.SetToyType(BootStrap::kBootstrap);
  bNominal.run();

  TGraphAsymmErrors *gNominal = bNominal.result(BootStrap::kMin,.68);
  gNominal->SetName("gBootStrap");

  unfoldNominal.SetNToys(1000);
  
  hUnfold_Nominal=(TH1D*) unfoldNominal.Hreco(RooUnfold::kCovToy);
  hUnfold_Nominal->SetXTitle("p_{T}^{#mu^{-}} [GeV]");
  
  for(int j=0;j!=hUnfold_Nominal->GetNbinsX();++j)
    {
      hUnfold_Nominal->SetBinContent(j+1,hUnfold_Nominal->GetBinContent(j+1)/hUnfold_Nominal->GetBinWidth(j+1));
      hUnfold_Nominal->SetBinError(j+1,hUnfold_Nominal->GetBinError(j+1)/hUnfold_Nominal->GetBinWidth(j+1));
    }

  for(int j=0;j!=hTruth_Nominal->GetNbinsX();++j)
    {
      hTruth_Nominal->SetBinContent(j+1,hTruth_Nominal->GetBinContent(j+1)/hTruth_Nominal->GetBinWidth(j+1));
      hTruth_Nominal->SetBinError(j+1,hTruth_Nominal->GetBinError(j+1)/hTruth_Nominal->GetBinWidth(j+1));
    }

  hUnfold_Nominal->Scale(1./LUMI);
  hTruth_Nominal->Scale(1./LUMI);

  his->cd();
  hTruth_Nominal->SetName("hTruth");
  hTruth_Nominal->Write("hTruth");
  hUnfold_Nominal->SetName("hUnfold");
  hUnfold_Nominal->Write("hUnfold");
  gNominal->Write();

  Train();
  Test();
  TDirectory *cdUNFOLDMODEL = his->mkdir("UNFOLDMODEL");
  cdUNFOLDMODEL->cd();

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
  hUnfold->SetXTitle("p_{T}^{#mu^{-}} [GeV]");
  
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

  hTruth->SetName("hTruth");
  hTruth->Write("hTruth");
  hUnfold->SetName("hUnfold");
  hUnfold->Write("hUnfold");
  g->Write();

  his->Write();
  his->Close();
}

#ifndef __CINT__
int main (int argc,char **argv) {RooUnfoldExample(argv[1],atof(argv[2]),atof(argv[3])); return 0; }  // Main program when run stand-alone
#endif
