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
float mean(vector<float> &a )
{
  float S=0;
  for(int i=0;i<int(a.size());i++) S+=a[i];
  S/=a.size();
  return S;
}
float rms(vector<float> &a)
{
  float S=0;
  float m=mean(a);
  for(int i=0;i<int(a.size());i++) S+=(a[i]-m)*(a[i]-m);
  S/=(a.size()-1);
  return sqrt(S);
}

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;


TH1D *hReco[100], *hTruth[100], *hData[100], *hTop[100], *hEWK[100], *hMeas[100], *hUnfold[100];
TH2D *hMatrix[100], *hMatrix_hilf[100];
TH1D *hReco_Nominal, *hTruth_Nominal, *hData_Nominal, *hTop_Nominal, *hEWK_Nominal, *hMeas_Nominal, *hUnfold_Nominal;
TH2D *hMatrix_Nominal, *hMatrix_hilf_Nominal;

vector<float> vrms;

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
  TFile *file1= new TFile("../UnfoldingInput/Zmumu/zmm_UnfoldInputs_ResScale.root");

  if (file1 == NULL) cout <<"File does not exists"<<endl;
  for(int i=0;i!=100;++i)
    {
      string si=int2string(i);
      hTruth[i] = (TH1D*)file1->Get((string("hLepNegPtTruth_")+si).c_str());
      hReco[i] = (TH1D*)file1->Get((string("hLepNegPtReco_")+si).c_str());
      hMatrix_hilf[i] = (TH2D*)file1->Get((string("hLepNegPtMatrix_")+si).c_str());
      
      if (hTruth[i]==NULL) cout<< "hTruth does not exist"<<endl;
      if (hReco[i]==NULL) cout<< "hReco does not exist"<<endl;
      if (hMatrix_hilf[i]==NULL) cout<< "hMatrix does not exist"<<endl;
      
      hMatrix[i]= (TH2D*)hMatrix_hilf[i]->Clone((string("hMatrix_")+si).c_str());
      for(int j=0;j!=hMatrix_hilf[i]->GetNbinsX();++j)
	{
	  for(int k=0;k!=hMatrix_hilf[i]->GetNbinsY();++k)
	    {
	      hMatrix[i]->SetBinContent(j+1,k+1,hMatrix_hilf[i]->GetBinContent(k+1,j+1));
	      hMatrix[i]->SetBinError(j+1,k+1,hMatrix_hilf[i]->GetBinError(k+1,j+1));
	    }
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
  TFile *f=new TFile("../SignalExtraction/Zmumu/Zmm_DataBkg_ResScale.root");
  
  if (f == NULL) cout<<"file does not exists"<<endl;
  for(int i=0;i!=100;++i)
    {
      string si=int2string(i);
      hData[i] = (TH1D*)f->Get((string("hDataLepNegPt_")+si).c_str());
      hTop[i] = (TH1D*)f->Get((string("hTopLepNegPt_")+si).c_str());
      hEWK[i] = (TH1D*)f->Get((string("hEWKLepNegPt_")+si).c_str());
      
      hMeas[i] = (TH1D*)f->Get((string("hDataLepNegPt_")+si).c_str());
      hMeas[i]->Sumw2();
      hMeas[i]->Add(hTop[i],-1);
      hMeas[i]->Add(hEWK[i],-1);
    }
      
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
 
  Train();
  Test();
  
  cout << "==================================== UNFOLD ===================================" << endl;

  TFile* his=new TFile(outputDir+(string("/UnfoldingOutputLepNegPtResScale.root")).c_str(), "recreate");

  RooUnfoldResponse responseNominal(hReco_Nominal,hTruth_Nominal,hMatrix_Nominal);

  responseNominal.UseOverflow();

  RooUnfoldBayes   unfoldNominal (&responseNominal, hMeas_Nominal, iterations);

  hUnfold_Nominal=(TH1D*) unfoldNominal.Hreco(RooUnfold::kNoError);
  
  for(int j=0;j!=hUnfold_Nominal->GetNbinsX();++j)
    {
      hUnfold_Nominal->SetBinContent(j+1,hUnfold_Nominal->GetBinContent(j+1)/hUnfold_Nominal->GetBinWidth(j+1));
    }

  for(int j=0;j!=hTruth_Nominal->GetNbinsX();++j)
    {
      hTruth_Nominal->SetBinContent(j+1,hTruth_Nominal->GetBinContent(j+1)/hTruth_Nominal->GetBinWidth(j+1));
    }

  hUnfold_Nominal->Scale(1./LUMI);
  hTruth_Nominal->Scale(1./LUMI);

  for(int i=0;i!=100;++i)
    {
      RooUnfoldResponse response(hReco[i],hTruth[i],hMatrix[i]);

      response.UseOverflow();
      

      RooUnfoldBayes   unfold (&response, hMeas[i], iterations);
      
      hUnfold[i]=(TH1D*) unfold.Hreco(RooUnfold::kNoError);
      
      for(int j=0;j!=hUnfold[i]->GetNbinsX();++j)
	{
	  hUnfold[i]->SetBinContent(j+1,hUnfold[i]->GetBinContent(j+1)/hUnfold[i]->GetBinWidth(j+1));
	}
      hUnfold[i]->Scale(1./LUMI);
    }

  int Nbins=hUnfold[0]->GetNbinsX();
  vector<float> vtoyval;
  for(int i=0;i!=Nbins;++i)
    {
      vtoyval.clear();
      for(int j=0;j!=100;++j)
	{
	  vtoyval.push_back(hUnfold[j]->GetBinContent(i+1));
	}
      vrms.push_back(rms(vtoyval));
      cout<<"RMS: "<<vrms[i]<<endl;
    }

  for(int j=0;j!=hUnfold_Nominal->GetNbinsX();++j)
    {
      hUnfold_Nominal->SetBinError(j+1,vrms[j]);
    }

  his->cd();
  hTruth_Nominal->Write("hTruth");
  hUnfold_Nominal->Write("hUnfold");
  his->Close();
}

#ifndef __CINT__
int main (int argc,char **argv) {RooUnfoldExample(argv[1],atof(argv[2]),atof(argv[3])); return 0; }  // Main program when run stand-alone
#endif
