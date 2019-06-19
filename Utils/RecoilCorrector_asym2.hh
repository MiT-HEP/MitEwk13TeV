#include <vector>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TGraphErrors.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooGenericPdf.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooKeysPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include <TFitResult.h>
#include "PdfDiagonalizer.h"

#include <stdio.h>
#include <math.h>
//
// ** apply recoil corrections **
// 
// usage: 
//    double met=rawMetValue;
//    double metphi=rawMetPhiValue;
//    RecoilCorrector corrector;
//    corrector->CorrectType*(met,metphi,GenZPt,GenZPhi,leptonPt,leptonPhi,u1,u2,sigma_mean,sigma_rms);
//    printf("corrected met: %10.2f%10.2f\n",met,metphi);
//
// where leptonPt, leptonPhi are dilepton kinematics for z->ll and single lepton kinematics for w->lnu
//

using namespace std;
using namespace RooFit;

class RecoilCorrector
{
  
public:
  RecoilCorrector(string iNameZDat, int iSeed=0xDEADBEEF);
  RecoilCorrector(string iNameZDat1, string iPrefix, int iSeed=0xDEADBEEF);
    
  void loadRooWorkspacesMCtoCorrectKeys(string iNameFile);
  void loadRooWorkspacesMCtoCorrect(string iNameFile);
  void loadRooWorkspacesMC(string iNameFile);
  void loadRooWorkspacesDiagMCtoCorrect(string iNameFile, int sigma);
  void loadRooWorkspacesDiagMC(string iNameFile, int sigma);
  void loadRooWorkspacesDiagData(string iNameFile, int sigma);
  void loadRooWorkspacesData(string iNameFile);
  void loadFileRatio(string iNameFile);
  
  void CorrectType0(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectType2(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectInvCdf(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0, bool dokeys=false, bool doDiago=false);
  void CorrectShitty(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2, double &evtWeight);
  void CorrectFromToys(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void addDataFile(std::string iNameDat);
  void addMCFile  (std::string iNameMC);
  void addFileWithGraph  (std::string iNameMC);
  Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResult  *fs);
  Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResult *fs);
  TFitResult *fitresPFu1mean, *fitresPFu1sigma1, *fitresPFu1sigma2,  *fitresPFu1sigma0;
  TFitResult *fitresPFu2mean, *fitresPFu2sigma1, *fitresPFu2sigma2,  *fitresPFu2sigma0;
  void runDiago(RooWorkspace *w, RooFitResult *result,int i, RooAbsReal *&pdfUiCdf, int sigma);
  void statUnc50nsStyle(RooWorkspace *w,  int i, RooAbsReal *&pdfUiCdf, int sigma);
  
protected:
  enum Recoil { 
    PFU1,
    PFU2,
    PFMSU1,
    PFMSU2,
    PFS1U1,
    PFS2U1,
    PFS1U2,
    PFS2U2,
    TKU1,
    TKU2,
    TKMSU1,
    TKMSU2,
    TKS1U1,
    TKS2U1,
    TKS1U2,
    TKS2U2
  };
  

  
  void readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,
          std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit,
          std::string iFName,std::string iPrefix, int type); 
  
  void readRecoilMC(std::vector<TGraphErrors*> &iU1Fit,std::vector<TGraphErrors*> &iU1MRMSFit,std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,
          std::vector<TGraphErrors*> &iU2Fit,std::vector<TGraphErrors*> &iU2MRMSFit,std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,
          std::string iFName,std::string iPrefix); 

  void readRecoilMC2(std::vector<TGraphErrors*> &iU1Fit,std::vector<TGraphErrors*> &iU1Fit2,std::vector<TGraphErrors*> &iU1MRMSFit,std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,
          std::vector<TGraphErrors*> &iU2Fit,std::vector<TGraphErrors*> &iU2MRMSFit,std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,
          std::string iFName,std::string iPrefix); 
          
                
  void metDistributionType0(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                TGraphErrors *iU1Default,
                TGraphErrors *iU1RZDatFit, 
                TGraphErrors *iU1RZDatFit2, 
                TGraphErrors *iU1MSZDatFit,  
                TGraphErrors *iU1S1ZDatFit,  
                TGraphErrors *iU1S2ZDatFit, 
                TGraphErrors *iU2RZDatFit,
                TGraphErrors *iU2MSZDatFit, 
                TGraphErrors *iU2S1ZDatFit,       
                TGraphErrors *iU2S2ZDatFit,      
                double &iU1, double &iU2,double iFluc=0,double iScale=0);
                
                     
  void metDistributionType2(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                TGraphErrors *iU1Default,
                TF1 *iU1RZDatFit,  TGraphErrors *iU1RZMCFit,
                TF1 *iU1MSZDatFit, TGraphErrors *iU1MSZMCFit, 
                TF1 *iU1S1ZDatFit, TGraphErrors *iU1S1ZMCFit, 
                TF1 *iU1S2ZDatFit, TGraphErrors *iU1S2ZMCFit, 
                TF1 *iU2MSZDatFit, TGraphErrors *iU2MSZMCFit,
                TF1 *iU2S1ZDatFit, TGraphErrors *iU2S1ZMCFit,      
                TF1 *iU2S2ZDatFit, TGraphErrors *iU2S2ZMCFit,      
                double &iU1, double &iU2,double iFluc=0,double iScale=0);
                
  void metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                TGraphErrors *iU1Default,
                double &iU1, double &iU2,double iFluc=0,double iScale=0);
                
  void metDistributionShitty(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                double &iU1, double &iU2, double &evtWeight);
  
  void metDistributionFromToys(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                double &iU1, double &iU2,double iFluc=0,double iScale=0);            
                
                
  Double_t calcErrorGraph(const TGraphErrors *graph, const Double_t x, const Double_t bins[], const Int_t bin);

  double diGausPVal    (double iVal, double iFrac,double iSimga1,double iSigma2);
  double triGausPVal    (double iVal, double iFrac2,double iFrac3,double iSimga1,double iSigma2,double iSigma3);
  double diGausPInverse(double iPVal,double iFrac,double iSigma1,double iSigma2);
  
  double triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *varDat, RooRealVar *varMC, int bin,double max);
  
  double calculate(int iMet,double iEPt,double iEPhi,double iWPhi,double iU1,double iU2);
  double getCorError2(double iVal,TF1 *iFit);
  double mag(double iV0,double iV1,double iV2,double iV3);
  double correlatedSeed(double iVal, double iCorr1,double iCorr2,double iCorr3,double iSeed0,double iSeed1,double iSeed2,double iSeed3);
  double deCorrelate   (double iVal, double iCorr1,double iCorr2,double iCorr3,double iSeed0,double iSeed1,double iSeed2,double iSeed3);
  TF1*   getFunc(bool iMC, Recoil iType);
  double CorrVal(double iPt,double iVal,Recoil iType);
  
  double getError(double iVal,TF1 *iFit,Recoil iType);
  double getError2(double iVal,TF1 *iFit);
  double getErrorSigma(double iVal,TF1 *iFit,Recoil iType);
  double getErrorMean(double iVal,TF1 *iFit,Recoil iType);

  
  TCanvas *iC = new TCanvas("C1","C1",800,600); 
  
  TRandom1 *fRandom; 
  vector<TGraphErrors*> fF1U1Fit; vector<TGraphErrors*> fF1U1RMSSMFit; vector<TGraphErrors*> fF1U1RMS1Fit; vector<TGraphErrors*> fF1U1RMS2Fit; 
  vector<TGraphErrors*> fF1U2Fit; vector<TGraphErrors*> fF1U2RMSSMFit; vector<TGraphErrors*> fF1U2RMS1Fit; vector<TGraphErrors*> fF1U2RMS2Fit; 
  vector<TGraphErrors*> fF2U1Fit; vector<TGraphErrors*> fF2U1RMSSMFit; vector<TGraphErrors*> fF2U1RMS1Fit; vector<TGraphErrors*> fF2U1RMS2Fit; 
  vector<TGraphErrors*> fF2U2Fit; vector<TGraphErrors*> fF2U2RMSSMFit; vector<TGraphErrors*> fF2U2RMS1Fit; vector<TGraphErrors*> fF2U2RMS2Fit; 
  // means 2 and 3, sigma 3
  vector<TGraphErrors*> fF1U1Fit2; 

  vector<TF1*> fD1U1Fit; vector<TF1*> fD1U1RMSSMFit; vector<TF1*> fD1U1RMS1Fit; vector<TF1*> fD1U1RMS2Fit; 
  vector<TF1*> fD1U2Fit; vector<TF1*> fD1U2RMSSMFit; vector<TF1*> fD1U2RMS1Fit; vector<TF1*> fD1U2RMS2Fit; 
  vector<TF1*> fD2U1Fit; vector<TF1*> fD2U1RMSSMFit; vector<TF1*> fD2U1RMS1Fit; vector<TF1*> fD2U1RMS2Fit; 
  vector<TF1*> fD2U2Fit; vector<TF1*> fD2U2RMSSMFit; vector<TF1*> fD2U2RMS1Fit; vector<TF1*> fD2U2RMS2Fit; 
  
  vector<TF1*> fD1U1FitFun; vector<TF1*> fD1U1RMSSMFitFun; vector<TF1*> fD1U1RMS1FitFun; vector<TF1*> fD1U1RMS2FitFun; 
  vector<TF1*> fD1U2FitFun; vector<TF1*> fD1U2RMSSMFitFun; vector<TF1*> fD1U2RMS1FitFun; vector<TF1*> fD1U2RMS2FitFun; 
  vector<TF1*> fD2U1FitFun; vector<TF1*> fD2U1RMSSMFitFun; vector<TF1*> fD2U1RMS1FitFun; vector<TF1*> fD2U1RMS2FitFun; 
  vector<TF1*> fD2U2FitFun; vector<TF1*> fD2U2RMSSMFitFun; vector<TF1*> fD2U2RMS1FitFun; vector<TF1*> fD2U2RMS2FitFun; 

  
  vector<TF1*> fFunM1U1Fit; vector<TF1*> fFunM1U1RMSSMFit; vector<TF1*> fFunM1U1RMS1Fit; vector<TF1*> fFunM1U1RMS2Fit; 
  vector<TF1*> fFunM1U2Fit; vector<TF1*> fFunM1U2RMSSMFit; vector<TF1*> fFunM1U2RMS1Fit; vector<TF1*> fFunM1U2RMS2Fit; 
  vector<TF1*> fFunM2U1Fit; vector<TF1*> fFunM2U1RMSSMFit; vector<TF1*> fFunM2U1RMS1Fit; vector<TF1*> fFunM2U1RMS2Fit; 
  vector<TF1*> fFunM2U2Fit; vector<TF1*> fFunM2U2RMSSMFit; vector<TF1*> fFunM2U2RMS1Fit; vector<TF1*> fFunM2U2RMS2Fit; 
  
  vector<TGraphErrors*> fM1U1Fit; vector<TGraphErrors*> fM1U1RMSSMFit; vector<TGraphErrors*> fM1U1RMS1Fit; vector<TGraphErrors*> fM1U1RMS2Fit; 
  vector<TGraphErrors*> fM1U2Fit; vector<TGraphErrors*> fM1U2RMSSMFit; vector<TGraphErrors*> fM1U2RMS1Fit; vector<TGraphErrors*> fM1U2RMS2Fit; 
  // means 2 and 3, and sigma 3
  vector<TGraphErrors*> fM1U1Fit2; 
  
  
  vector<TGraphErrors*> fM2U1Fit; vector<TGraphErrors*> fM2U1RMSSMFit; vector<TGraphErrors*> fM2U1RMS1Fit; vector<TGraphErrors*> fM2U1RMS2Fit; 
  vector<TGraphErrors*> fM2U2Fit; vector<TGraphErrors*> fM2U2RMSSMFit; vector<TGraphErrors*> fM2U2RMS1Fit; vector<TGraphErrors*> fM2U2RMS2Fit; 

  vector<TGraphErrors*> fMT1U1Fit; vector<TGraphErrors*> fMT1U1RMSSMFit; vector<TGraphErrors*> fMT1U1RMS1Fit; vector<TGraphErrors*> fMT1U1RMS2Fit; 
  vector<TGraphErrors*> fMT1U2Fit; vector<TGraphErrors*> fMT1U2RMSSMFit; vector<TGraphErrors*> fMT1U2RMS1Fit; vector<TGraphErrors*> fMT1U2RMS2Fit; 
  vector<TGraphErrors*> fMT2U1Fit; vector<TGraphErrors*> fMT2U1RMSSMFit; vector<TGraphErrors*> fMT2U1RMS1Fit; vector<TGraphErrors*> fMT2U1RMS2Fit; 
  vector<TGraphErrors*> fMT2U2Fit; vector<TGraphErrors*> fMT2U2RMSSMFit; vector<TGraphErrors*> fMT2U2RMS1Fit; vector<TGraphErrors*> fMT2U2RMS2Fit;
  vector<TGraphErrors*> fMT1U1Fit2; 

  vector<TGraphErrors*> fF1U1U2Corr;     vector<TGraphErrors*> fF2U1U2Corr;
  vector<TGraphErrors*> fF1F2U1Corr;     vector<TGraphErrors*> fF1F2U2Corr;
  vector<TGraphErrors*> fF1F2U1U2Corr;   vector<TGraphErrors*> fF1F2U2U1Corr;

  vector<TGraphErrors*> fM1U1U2Corr;     vector<TGraphErrors*> fM2U1U2Corr;
  vector<TGraphErrors*> fM1M2U1Corr;     vector<TGraphErrors*> fM1M2U2Corr;
  vector<TGraphErrors*> fM1M2U1U2Corr;   vector<TGraphErrors*> fM1M2U2U1Corr;
  
  RooWorkspace* rooWData[2];
  RooWorkspace* rooWMC[2];
  RooWorkspace* rooWMCtoCorr[2];
  RooWorkspace* rooWDataDiag[2];
  RooWorkspace* rooWMCDiag[2];
  RooWorkspace* rooWMCtoCorrDiag[2];
  RooWorkspace* pdfsU1zData, pdfsU2zData;
  RooWorkspace* pdfsU1zMC, pdfsU2zMC;
  RooWorkspace* pdfsU1sigMC, pdfsU2sigMC;
  int fId; int fJet;
  bool dokeys; bool doDiago;
  
  RooWorkspace rooWksDataU1;
  RooWorkspace rooWksMCU1;
  RooWorkspace rooWksDataU2;
  RooWorkspace rooWksMCU2;
  
  // oct2 binning
  //  std::vector<double> vZPtBins ={0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};

  // oct7 binning
  //  std::vector<double> vZPtBins ={0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};

  // may22 binning
  // bins from 2015 recoil corrections
  // std::vector<double> vZPtBins = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};
  // bins for the low pileup corrections
  std::vector<double> vZPtBins = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,85,90,95,100,120,140,160,180,200,220,250,300}; // the regular 13 tev 2017 one
  
  // std::vector<double> vZPtBins = {0,2.5,5.0,10,20,30,40,50,60,80,100,125,150,200,250,300}; // rebin to check effect of bin size on statitiscal unc (2017, 13 TeV)
             // Double_t ptbins[] = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};
 // std::vector<double> vZPtBins = {0,300};

  int nBins = vZPtBins.size()-1;
  TH1D **hRatiosU1 = new TH1D*[nBins];
  TH1D **hRatiosU2 = new TH1D*[nBins];
 
//   Double_t vZPtBins[] = {0,1,2.5,5.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100};
// int nZPtBins = sizeof(vZPtBins)/sizeof(Double_t)-1;

};

//-----------------------------------------------------------------------------------------------------------------------------------------
RecoilCorrector::RecoilCorrector(string iNameZDat,std::string iPrefix, int iSeed) {

  fRandom = new TRandom1(iSeed);
  readRecoilMC2(fF1U1Fit,fF1U1Fit2,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,iNameZDat,iPrefix);
  fId = 0; fJet = 0;
}

RecoilCorrector::RecoilCorrector(string iNameZ, int iSeed) {

  fRandom = new TRandom1(iSeed);
  // get fits for Z data
  readRecoilMC(fF1U1Fit,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,iNameZ,"PF");
  fId = 0; fJet = 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void RecoilCorrector::loadFileRatio(std::string iFName){
  
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName).c_str());
  // assert(lFile);
  int nBins = vZPtBins.size()-1;
  
  for(uint i = 0; i < nBins; ++i){
    char name[100];
    sprintf(name,"hRecoilU1_bin%d",i+1);
    hRatiosU1[i] = (TH1D*) lFile->Get(name);
    sprintf(name,"hRecoilU2_bin%d",i);
    hRatiosU2[i] = (TH1D*) lFile->Get(name);
  }
}

void RecoilCorrector::loadRooWorkspacesData(std::string iFName){
  
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWData[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  lFile->Delete();
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWData[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  lFile2->Delete();
  
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    TString name;
    name = Form("sig_%i",i);
    RooAbsPdf* pdf1 = rooWData[0]->pdf(name);
    RooAbsPdf* pdf2 = rooWData[1]->pdf(name);
    name = Form("u_%i",i);
    RooRealVar* myX1 = (RooRealVar*) rooWData[0]->var(name);
    RooRealVar* myX2 = (RooRealVar*) rooWData[1]->var(name);
    RooAbsReal *cdfU1 = pdf1->createCdf(*myX1);
    rooWData[0]->import(*cdfU1, RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWData[1]->import(*cdfU2, RooFit::Silence());
    
  }
  

  
  std::cout << "Loaded Workspaces...DATA "<< std::endl;
}


void RecoilCorrector::loadRooWorkspacesDiagMCtoCorrect(std::string iFName, int sigma){
  std::cout << iFName << std::endl;
  std::cout << "aaaaaaa" << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMCtoCorrDiag[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMCtoCorrDiag[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    RooFitResult* fitresultU1; //= (RooFitResult*)lFile->FindObject(Form("fitResultU1_%d",i));
    lFile->GetObject(Form("fitResultU1_%d",i),fitresultU1);
    RooFitResult* fitresultU2;//    = (RooFitResult*)lFile2->FindObject(Form("fitResultU2_%d",i));
    lFile2->GetObject(Form("fitResultU2_%d",i),fitresultU2);
    TString name;
    name = Form("sig_%d",i);
    RooAbsReal *cdfU1, *cdfU2;

    runDiago(rooWMCtoCorrDiag[0],fitresultU1,i,cdfU1,sigma);
    runDiago(rooWMCtoCorrDiag[1],fitresultU2,i,cdfU2,sigma);
  }

  lFile->Delete();
  lFile2->Delete();
  std::cout << "Loaded Workspaces Source MC - Stat Unc "<< std::endl;
}

void RecoilCorrector::loadRooWorkspacesDiagMC(std::string iFName,int sigma){
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMCDiag[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMCDiag[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    RooFitResult* fitresultU1; //= (RooFitResult*)lFile->FindObject(Form("fitResultU1_%d",i));
    lFile->GetObject(Form("fitResultU1_%d",i),fitresultU1);
    RooFitResult* fitresultU2;//    = (RooFitResult*)lFile2->FindObject(Form("fitResultU2_%d",i));
    lFile2->GetObject(Form("fitResultU2_%d",i),fitresultU2);
    TString name;
    name = Form("sig_%d",i);
    RooAbsReal *cdfU1, *cdfU2;
	  runDiago(rooWMCDiag[0],fitresultU1,i,cdfU1,sigma);
	  runDiago(rooWMCDiag[1],fitresultU2,i,cdfU2,sigma);
  }
  lFile->Delete();
  lFile2->Delete();
  std::cout << "Loaded Workspaces Z MC - Stat Unc"<< std::endl;
}

void RecoilCorrector::loadRooWorkspacesDiagData(std::string iFName,int sigma){
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWDataDiag[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWDataDiag[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    RooFitResult* fitresultU1; //= (RooFitResult*)lFile->FindObject(Form("fitResultU1_%d",i));
    lFile->GetObject(Form("fitResultU1_%d",i),fitresultU1);
    RooFitResult* fitresultU2;//    = (RooFitResult*)lFile2->FindObject(Form("fitResultU2_%d",i));
    lFile2->GetObject(Form("fitResultU2_%d",i),fitresultU2);
    TString name;
    name = Form("sig_%d",i);
    RooAbsReal *cdfU1, *cdfU2;
    runDiago(rooWDataDiag[0],fitresultU1,i,cdfU1,sigma);
    runDiago(rooWDataDiag[1],fitresultU2,i,cdfU2,sigma);
  }
  lFile->Delete();
  lFile2->Delete();
  std::cout << "Loaded Workspaces Data - Stat Unc "<< std::endl;
}

void RecoilCorrector::loadRooWorkspacesMC(std::string iFName){
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMC[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  lFile->Delete();
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMC[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  lFile2->Delete();
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    std::stringstream name;
    name << "sig_" << i;
    RooAbsPdf* pdf1 = rooWMC[0]->pdf(name.str().c_str());
    RooAbsPdf* pdf2 = rooWMC[1]->pdf(name.str().c_str());
    name.str(""); name << "u_" << i;
    RooRealVar* myX1 = (RooRealVar*) rooWMC[0]->var(name.str().c_str());
    RooRealVar* myX2 = (RooRealVar*) rooWMC[1]->var(name.str().c_str());
    RooAbsReal *cdfU1 = pdf1->createCdf(*myX1);
    rooWMC[0]->import(*cdfU1, RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWMC[1]->import(*cdfU2, RooFit::Silence());
  }
  std::cout << "Loaded Workspaces...Z MC "<< std::endl;
}

void RecoilCorrector::loadRooWorkspacesMCtoCorrect(std::string iFName){

  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMCtoCorr[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  lFile->Delete();
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMCtoCorr[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  lFile2->Delete();
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    TString name;
    name =Form("sig_%i", i);
    RooAbsPdf* pdf1 = rooWMCtoCorr[0]->pdf(name);
    RooAbsPdf* pdf2 = rooWMCtoCorr[1]->pdf(name);
    name =Form( "u_%i", i);
    RooRealVar* myX1 = (RooRealVar*) rooWMCtoCorr[0]->var(name);
    RooRealVar* myX2 = (RooRealVar*) rooWMCtoCorr[1]->var(name);
    RooAbsReal *cdfU1 = pdf1->createCdf(*myX1);
    rooWMCtoCorr[0]->import(*cdfU1, RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWMCtoCorr[1]->import(*cdfU2, RooFit::Silence());
  }
  std::cout << "Loaded Workspaces...Source MC "<< std::endl;
}


void RecoilCorrector::loadRooWorkspacesMCtoCorrectKeys(std::string iFName){

  //  RooKeysPdf::key_44

  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMCtoCorr[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  lFile->Delete();
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMCtoCorr[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  lFile2->Delete();
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    std::stringstream name;
    name << "key_" << i ;
    RooAbsPdf* pdf1 = (RooKeysPdf*) rooWMCtoCorr[0]->pdf(name.str().c_str());
    RooAbsPdf* pdf2 = (RooKeysPdf*) rooWMCtoCorr[1]->pdf(name.str().c_str());
    name.str(""); name << "u_" << i;
    RooRealVar* myX1 = (RooRealVar*) rooWMCtoCorr[0]->var(name.str().c_str());
    RooRealVar* myX2 = (RooRealVar*) rooWMCtoCorr[1]->var(name.str().c_str());
    RooAbsReal *cdfU1 = pdf1->createCdf(*myX1);
    rooWMCtoCorr[0]->import(*cdfU1, RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWMCtoCorr[1]->import(*cdfU2, RooFit::Silence());
  }
  std::cout << "Loaded Workspaces Source MC - Keys"<< std::endl;
}



void RecoilCorrector::addDataFile(std::string iNameData) {
  readRecoil(fD1U1Fit,fD1U1RMSSMFit,fD1U1RMS1Fit,fD1U1RMS2Fit,fD1U2Fit,fD1U2RMSSMFit,fD1U2RMS1Fit,fD1U2RMS2Fit,iNameData,"fcnPF",0);
  fId++;   
}
void RecoilCorrector::addMCFile  (std::string iNameMC) {
  fId++;
  readRecoilMC(fM1U1Fit,fM1U1RMSSMFit,fM1U1RMS1Fit,fM1U1RMS2Fit,fM1U2Fit,fM1U2RMSSMFit,fM1U2RMS1Fit,fM1U2RMS2Fit,iNameMC,"fcnPF");
}

void RecoilCorrector::addFileWithGraph  (std::string iNameMC) { 
  fId++;
  readRecoilMC(fMT1U1Fit,fMT1U1RMSSMFit,fMT1U1RMS1Fit,fMT1U1RMS2Fit,fMT1U2Fit,fMT1U2RMSSMFit,fMT1U2RMS1Fit,fMT1U2RMS2Fit,iNameMC,"fcnPF");
}


void RecoilCorrector::CorrectType0(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet) {
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0;
  metDistributionType0(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi,fF1U1Fit[fJet],
               fMT1U1Fit     [fJet],
               fMT1U1Fit2    [fJet],
               fMT1U1RMSSMFit[fJet],
               fMT1U1RMS1Fit [fJet],
               fMT1U1RMS2Fit [fJet],
               fMT1U2Fit     [fJet],
               fMT1U2RMSSMFit[fJet],
               fMT1U2RMS1Fit [fJet],
               fMT1U2RMS2Fit [fJet], 
               iU1,iU2,iFluc,iScale);
}

void RecoilCorrector::CorrectType2(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet) {  
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0;
  
  metDistributionType2(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi,fF1U1Fit[fJet],
               fD1U1Fit     [fJet],fM1U1Fit     [fJet],
               fD1U1RMSSMFit[fJet],fM1U1RMSSMFit[fJet],
               fD1U1RMS1Fit [fJet],fM1U1RMS1Fit [fJet],
               fD1U1RMS2Fit [fJet],fM1U1RMS2Fit [fJet],
               fD1U2RMSSMFit[fJet],fM1U2RMSSMFit[fJet],
               fD1U2RMS1Fit [fJet],fM1U2RMS1Fit [fJet],
               fD1U2RMS2Fit [fJet],fM1U2RMS2Fit [fJet], 
               iU1,iU2,iFluc,iScale);
}

void RecoilCorrector::CorrectInvCdf(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet, bool useKeys, bool diago) {
  dokeys=useKeys;
  doDiago=diago;
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0; 
  metDistributionInvCdf(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi, fF1U1Fit[fJet], iU1,iU2,iFluc,iScale);
}

void RecoilCorrector::CorrectShitty(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2, double &evtWeight) {
  metDistributionShitty(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi, iU1,iU2, evtWeight);
}

void RecoilCorrector::CorrectFromToys(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet) {  
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0; 

  metDistributionFromToys(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi, iU1,iU2,iFluc,iScale);
}

TF1* RecoilCorrector::getFunc(bool iMC, Recoil iType) { 
  if(fId == 0 || fId == 1) return 0;
  switch(iType) {
  case PFU1   : return fD1U1Fit     [fJet];
  case PFMSU1 : return fD1U1RMSSMFit[fJet]; 
  case PFS1U1 : return fD1U1RMS1Fit [fJet]; 
  case PFS2U1 : return fD1U1RMS2Fit [fJet]; 
  case PFU2   : return 0;
  case PFMSU2 : return fD1U2RMSSMFit[fJet]; 
  case PFS1U2 : return fD1U2RMS1Fit [fJet]; 
  case PFS2U2 : return fD1U2RMS2Fit [fJet]; 
  case TKU1   : return fD2U1Fit     [fJet];
  case TKMSU1 : return fD2U1RMSSMFit[fJet];
  case TKS1U1 : return fD2U1RMS1Fit [fJet];
  case TKS2U1 : return fD2U1RMS2Fit [fJet];
  case TKU2   : return 0;
  case TKMSU2 : return fD2U2RMSSMFit[fJet];
  case TKS1U2 : return fD2U2RMS1Fit [fJet];
  case TKS2U2 : return fD2U2RMS2Fit [fJet];
  }
  return 0;
}


//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,
                         std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit,
                         std::string iFName,std::string iPrefix, int type) {
  if(type==0){
  TFile *lFile  = new TFile(iFName.c_str());
  int lNJet = 0;
  std::stringstream lSS; //lSS << iPrefix;
  lSS << iPrefix << "u1mean"; 
  iU1Fit.push_back    ( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1sigma0";
  iU1MRMSFit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma1"; 
  iU1RMS1Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma2"; iU1RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u2mean"; iU2Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma0"; 
  iU2MRMSFit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma1"; iU2RMS1Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma2"; iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  fitresPFu1mean  = (TFitResult*)(lFile->Get("fitresPFu1mean")->Clone("fitresPFu1mean"));    
  fitresPFu2mean  = (TFitResult*)(lFile->Get("fitresPFu2mean")->Clone("fitresPFu2mean"));
  fitresPFu1sigma0  = (TFitResult*)(lFile->Get("fitresPFu1sigma0")->Clone("fitresPFu1sigma0"));
  fitresPFu2sigma0  = (TFitResult*)(lFile->Get("fitresPFu2sigma0")->Clone("fitresPFu2sigma0"));
  fitresPFu1sigma1  = (TFitResult*)(lFile->Get("fitresPFu1sigma1")->Clone("fitresPFu1sigma1"));
  fitresPFu2sigma1  = (TFitResult*)(lFile->Get("fitresPFu2sigma1")->Clone("fitresPFu2sigma1"));
  fitresPFu1sigma2  = (TFitResult*)(lFile->Get("fitresPFu1sigma2")->Clone("fitresPFu1sigma2"));
  fitresPFu2sigma2  = (TFitResult*)(lFile->Get("fitresPFu2sigma2")->Clone("fitresPFu2sigma2"));

  lFile->Close();
  }
  
  if(type==1){
    TFile *lFile  = new TFile(iFName.c_str());
    iU1Fit.push_back    ( (TF1*) lFile->FindObjectAny((iPrefix+"u1Mean_0").c_str()));
    iU1MRMSFit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1MeanRMS_0").c_str()));
    iU1RMS1Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1RMS1_0").c_str()));
    iU1RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1RMS2_0").c_str()));
    iU2Fit    .push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2Mean_0").c_str()));
    iU2MRMSFit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2MeanRMS_0").c_str()));
    iU2RMS1Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS1_0").c_str()));
    iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS2_0").c_str()));
    lFile->Close(); 
  }
  
}


// This is from the really old 2015 analysis, maybe it can be deleted
//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::readRecoilMC(std::vector<TGraphErrors*> &iU1Fit,    std::vector<TGraphErrors*> &iU1MRMSFit,
                                   std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,
                                   std::vector<TGraphErrors*> &iU2Fit,    std::vector<TGraphErrors*> &iU2MRMSFit,
                                   std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,
                                   std::string iFName,std::string iPrefix) {
  TFile *lFile  = new TFile(iFName.c_str());
  int lNJet = 0;
  std::stringstream lSS; //lSS << iPrefix;
  iPrefix = "grPF"; 
  lSS << iPrefix << "u1mean"; 
  iU1Fit.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1sigma0";
  iU1MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma1"; 
  iU1RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma2"; iU1RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u2mean"    << lNJet; iU2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma0"; 
  iU2MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma1"; iU2RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma2"; iU2RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lFile->Close();
}

// This is also from the really old 2015 analysis, consider deleting
void RecoilCorrector::readRecoilMC2(std::vector<TGraphErrors*> &iU1Fit,std::vector<TGraphErrors*> &iU1Fit2,std::vector<TGraphErrors*> &iU1MRMSFit,std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,
                         std::vector<TGraphErrors*> &iU2Fit,std::vector<TGraphErrors*> &iU2MRMSFit,std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,
                         std::string iFName,std::string iPrefix) {
  TFile *lFile  = new TFile(iFName.c_str());
  int lNJet = 0;
  std::stringstream lSS; //lSS << iPrefix;
  iPrefix = "grPF"; 
  lSS << iPrefix << "u1mean"; 
  iU1Fit.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1mean2"; 
  iU1Fit2.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1sigma0";
  iU1MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma1"; 
  iU1RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma2"; iU1RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u2mean"    << lNJet; iU2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma0"; 
  iU2MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma1"; iU2RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma2"; iU2RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lFile->Close(); 
  
  
}

// consider renaming this function to something less stupid
double RecoilCorrector::triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *myXd, RooRealVar *myXm, int bin, double max) {
  if(iPVal< myXm->getMin()) return iPVal;
  if(iPVal> myXm->getMax()) return iPVal;
  myXm->setVal(iPVal);
  myXd->setVal(iPVal);
  pdfDATAcdf->getVal(); // do not delete this line, for some reason is necessary when using the diagonalized cdfs for the statistical uncertainty
  double pVal=pdfDATAcdf->findRoot(*myXd,myXd->getMin(),myXd->getMax(),pdfMCcdf->getVal());
  myXd->setVal(pVal);
  return pVal;
}

double RecoilCorrector::diGausPVal(double iVal,double iFrac,double iSigma1,double iSigma2) {
  return iFrac*TMath::Erf(iVal/iSigma1) + (1-iFrac)*TMath::Erf(iVal/iSigma2);
}

double RecoilCorrector::triGausPVal(double iVal,double iFrac2,double iFrac3,double iSigma1,double iSigma2,double iSigma3) {
  return (1-iFrac2-iFrac3)*TMath::Erf(iVal/iSigma1) + iFrac2*TMath::Erf(iVal/iSigma2)+ iFrac3*TMath::Erf(iVal/iSigma3);
}

double RecoilCorrector::diGausPInverse(double iPVal,double iFrac,double iSigma1,double iSigma2) {
  double lVal = TMath::ErfInverse(iPVal);
  double lMin = lVal * ((iSigma1 < iSigma2) ? iSigma1 : iSigma2); // lVal * sigma1
  double lMax = lVal * ((iSigma1 < iSigma2) ? iSigma2 : iSigma1); // lVal * sigma2
  double lDiff = (lMax-lMin); // lVal * (sigma2 - sigma1)
  //Iterative procedure to invert a double gaussian given a PVal
  int lId = 0; int lN1 = 4;  int lN2 = 10;  // Fewer toys
  //   int lId = 0; int lN1 = 10;  int lN2 = 100;  // More toys - takes longer
  for(int i0 = 0; i0 < lN1; i0++) {
    if(i0 != 0) lMin = lMin + (lId-1)*lDiff/lN2;
    if(i0 != 0) lDiff/=lN2;
    for(int i1 = 0; i1 < lN2; i1++) { 
      double pVal = lMin + lDiff/lN2*i1;
      pVal = diGausPVal(pVal,iFrac,iSigma1,iSigma2);
      if(pVal > iPVal) {lId = i1; break;}
    }
  }
  return (lMin + (lId-0.5)*lDiff/lN2);
}

// consider deleting this b/c it's super old
void RecoilCorrector::metDistributionType0(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
        double iLepPt,double iLepPhi,
        TGraphErrors *iU1Default,
        TGraphErrors *iU1RZDatFit, 
        TGraphErrors *iU1RZDatFit2, 
        TGraphErrors *iU1MSZDatFit,  
        TGraphErrors *iU1S1ZDatFit,  
        TGraphErrors *iU1S2ZDatFit, 
        TGraphErrors *iU2RZDatFit,
        TGraphErrors *iU2MSZDatFit, 
        TGraphErrors *iU2S1ZDatFit,                                              
        TGraphErrors *iU2S2ZDatFit, 
        double &iU1,double &iU2,double iFluc,double iScale){
  
  double pfu1mean   = iU1RZDatFit->Eval(iGenPt);
  double pfu1mean2   = iU1RZDatFit->Eval(iGenPt);
  double pfu1sigma1 = iU1S1ZDatFit->Eval(iGenPt);
  double pfu1sigma2 = iU1S2ZDatFit->Eval(iGenPt);
  double pfu1sigma0 = iU1MSZDatFit->Eval(iGenPt);
  
  // why
  if(iScale!=0 || iFluc!=0) {
//     pfu1mean   += iScale*dMean(iU1RZDatFit,iGenPt,fitresPFu1mean);
//     pfu1sigma1 += iFluc*dSigma(iU1S1ZDatFit,iGenPt,fitresPFu1sigma1);
//     pfu1sigma2 += iFluc*dSigma(iU1S2ZDatFit,iGenPt,fitresPFu1sigma2);
//     pfu1sigma0 += iFluc*dSigma(iU1MSZDatFit,iGenPt,fitresPFu1sigma0);
  }
  
  double pfu2mean   = iU2RZDatFit->Eval(iGenPt);
  double pfu2sigma1 = iU2S1ZDatFit->Eval(iGenPt);
  double pfu2sigma2 = iU2S2ZDatFit->Eval(iGenPt);
  double pfu2sigma0 = iU2MSZDatFit->Eval(iGenPt);
  
  if(iScale!=0 || iFluc!=0) {
//     pfu2mean   += iScale*dMean(iU2RZDatFit,iGenPt,fitresPFu2mean);
//     pfu2sigma1 += iFluc*dSigma(iU2S1ZDatFit,iGenPt,fitresPFu2sigma1);
//     pfu2sigma2 += iFluc*dSigma(iU2S2ZDatFit,iGenPt,fitresPFu2sigma2);
//     pfu2sigma0 += iFluc*dSigma(iU2MSZDatFit,iGenPt,fitresPFu2sigma0);
  }
  
  double pfu1frac2  = (pfu1sigma0 - pfu1sigma1)/(pfu1sigma2 - pfu1sigma1);
  double pfu2frac2  = (pfu2sigma0 - pfu2sigma1)/(pfu2sigma2 - pfu2sigma1);
  
// //   TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
  
  double z1 = gRandom->Gaus(0,1);
  double z2 = gRandom->Gaus(0,1);

  double pfu1 = (gRandom->Uniform(0,1) < pfu1frac2) ? z1*pfu1sigma2+pfu1mean2 : z1*pfu1sigma1+pfu1mean;
  double pfu2 = (gRandom->Uniform(0,1) < pfu2frac2) ? z2*pfu2sigma2+pfu2mean : z2*pfu2sigma1+pfu2mean;

  iU1 = pfu1;
  iU2 = pfu2;
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,iU1,iU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,iU1,iU2);
  
}

// Also consider deleting because it's super old
void RecoilCorrector::metDistributionType2(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       TGraphErrors *iU1Default,
                       TF1 *iU1RZDatFit,  TGraphErrors *iU1RZMCFit,
                       TF1 *iU1MSZDatFit, TGraphErrors *iU1MSZMCFit, 
                       TF1 *iU1S1ZDatFit, TGraphErrors *iU1S1ZMCFit, 
                       TF1 *iU1S2ZDatFit, TGraphErrors *iU1S2ZMCFit, 
                       TF1 *iU2MSZDatFit, TGraphErrors *iU2MSZMCFit,
                       TF1 *iU2S1ZDatFit, TGraphErrors *iU2S1ZMCFit,                                              
                       TF1 *iU2S2ZDatFit, TGraphErrors *iU2S2ZMCFit, 
                       double &iU1,double &iU2,double iFluc,double iScale) {
  double pDefU1    = iU1Default->Eval(iGenPt);
  double lRescale  = sqrt((TMath::Pi())/2.);             
  double pDU1       = iU1RZDatFit ->Eval(iGenPt);

  double pDFrac1    = iU1MSZDatFit->Eval(iGenPt);
  double pDSigma1_1 = iU1S1ZDatFit->Eval(iGenPt);
  double pDSigma1_2 = iU1S2ZDatFit->Eval(iGenPt);
  double pDFrac2    = iU2MSZDatFit->Eval(iGenPt);
  double pDSigma2_1 = iU2S1ZDatFit->Eval(iGenPt);
  double pDSigma2_2 = iU2S2ZDatFit->Eval(iGenPt);

  double pMU1       = iU1RZMCFit  ->Eval(iGenPt);
  double pMU2       = 0; 

  double pMFrac1    = iU1MSZMCFit ->Eval(iGenPt);
  double pMSigma1_1 = iU1S1ZMCFit ->Eval(iGenPt);
  double pMSigma1_2 = iU1S2ZMCFit ->Eval(iGenPt);
  double pMFrac2    = iU2MSZMCFit ->Eval(iGenPt);
  double pMSigma2_1 = iU2S1ZMCFit ->Eval(iGenPt);
  double pMSigma2_2 = iU2S2ZMCFit ->Eval(iGenPt);

  
  //Uncertainty propagation
  if(iFluc != 0 || iScale != 0) {
    fId = 0;
    double lEUR1    = dMean(iU1RZDatFit,iGenPt,fitresPFu1mean);
    double lEUS1_1  = dSigma(iU1S1ZDatFit,iGenPt,fitresPFu1sigma1);
    double lEUS1_2  = dSigma(iU1S2ZDatFit,iGenPt,fitresPFu1sigma2);
    double lEU1Frac = dSigma(iU1MSZDatFit,iGenPt,fitresPFu1sigma0);
    double lEUS2_1  = dSigma(iU2S1ZDatFit,iGenPt,fitresPFu2sigma1);
    double lEUS2_2  = dSigma(iU2S2ZDatFit,iGenPt,fitresPFu2sigma2);
    double lEU2Frac = dSigma(iU2MSZDatFit,iGenPt,fitresPFu2sigma0);
        
    //Modify all the different parameters the choice of signs makes it maximal
    pDU1       = pDU1       + iScale*lEUR1;            //Recoil
    pDFrac1    = pDFrac1    + iFluc*(lEU1Frac);        //Mean RMS 
    pDSigma1_1 = pDSigma1_1 + iFluc*lEUS1_1;           //Sigma 1 smalles sigma
    pDSigma1_2 = pDSigma1_2 + iFluc*lEUS1_2;           //Sigma 2 (Maximal when oppsite sigma 1)
    pDFrac2    = pDFrac2    + iFluc*(lEU2Frac);        //Mean RMS for U2
    pDSigma2_1 = pDSigma2_1 + iFluc*lEUS2_1;           //Sigma 1 U2
    pDSigma2_2 = pDSigma2_2 + iFluc*lEUS2_2;           //Sigma 2 (Maximal when oppsite sigma 1)
  }
  pDFrac1           = (pDFrac1-pDSigma1_2)/(pDSigma1_1-pDSigma1_2);
  pDFrac2           = (pDFrac2-pDSigma2_2)/(pDSigma2_1-pDSigma2_2);
  pMFrac1           = (pMFrac1-pMSigma1_2)/(pMSigma1_1-pMSigma1_2);
  pMFrac2           = (pMFrac2-pMSigma2_2)/(pMSigma2_1-pMSigma2_2);

  double pUX   = iMet*cos(iMPhi) + iLepPt*cos(iLepPhi);
  double pUY   = iMet*sin(iMPhi) + iLepPt*sin(iLepPhi);
  double pU    = sqrt(pUX*pUX+pUY*pUY);
  double pCos  = - (pUX*cos(iGenPhi) + pUY*sin(iGenPhi))/pU;
  double pSin  =   (pUX*sin(iGenPhi) - pUY*cos(iGenPhi))/pU;
  double pU1   = pU*pCos;
  double pU2   = pU*pSin;
  double pU1Diff  = pU1-pDefU1;
  double pU2Diff  = pU2;

  double p1Charge        = pU1Diff/fabs(pU1Diff);
  double p2Charge        = pU2Diff/fabs(pU2Diff);
  double pTU1Diff        = pU1Diff;
  double pU1ValM         = diGausPVal(fabs(pU1Diff),pMFrac1,pMSigma1_1,pMSigma1_2);
  double pU2ValM         = diGausPVal(fabs(pU2Diff),pMFrac2,pMSigma2_1,pMSigma2_2);
  double pU1ValD         = diGausPInverse(pU1ValM  ,pDFrac1,pDSigma1_1,pDSigma1_2);
  double pU2ValD         = diGausPInverse(pU2ValM  ,pDFrac2,pDSigma2_1,pDSigma2_2);
  
  pU1ValD*=p1Charge;
  pU2ValD*=p2Charge;
  pDefU1 *= (pDU1/pMU1);
  
  pU1   = pDefU1             + pU1ValD;
  pU2   =                      pU2ValD;
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iU1   = pU1; 
  iU2   = pU2;
  return;
}

// The one we actually use, clean it up a bit
void RecoilCorrector::metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       TGraphErrors *iU1Default,
                       double &iU1,double &iU2,double iFluc,double iScale) {
  double iGenPt2 = 0; // literally never used, also delete this variable
  Int_t nbinsPt = vZPtBins.size()-1;
  int iBin = -1;
  for(int i = 0; i < nbinsPt-1; ++i){
    if(iGenPt > vZPtBins[nbinsPt]){
      iBin = nbinsPt-1;
      iGenPt2 = (vZPtBins[nbinsPt-1]+vZPtBins[nbinsPt-2])*0.5;
      break;
    }
    if(vZPtBins[i+1] < iGenPt) continue;
    if(vZPtBins[i] > iGenPt ) continue;
    if(iGenPt < vZPtBins[i+1] && vZPtBins[i] < iGenPt){
      iBin = i; 
      iGenPt2 = (vZPtBins[i+1]+vZPtBins[i])*0.5;
      break;
    }
  }
  if(iBin<0) return;
  
  double pUX  = iMet*cos(iMPhi) + iLepPt*cos(iLepPhi);
  double pUY  = iMet*sin(iMPhi) + iLepPt*sin(iLepPhi);
  double pU   = sqrt(pUX*pUX+pUY*pUY);
  double pCos = - (pUX*cos(iGenPhi) + pUY*sin(iGenPhi))/pU;
  double pSin =   (pUX*sin(iGenPhi) - pUY*cos(iGenPhi))/pU;
  double pU1  = pU*pCos; // U1 in sample to Correct (WMC or ZMC)
  double pU2  = pU*pSin; // U2 in sample to Correct (WMC or ZMC)
  
  RooAbsPdf *thisPdfMCU1toCorr;
  RooAbsReal *thisCdfMCU1toCorr;
  RooAbsPdf *thisPdfMCU2toCorr;
  RooAbsReal *thisCdfMCU2toCorr;
  RooAbsPdf *thisPdfDataU1; RooAbsPdf *thisPdfMCU1; 
  RooAbsReal *thisCdfDataU1; RooAbsReal *thisCdfMCU1;
  RooAbsPdf *thisPdfDataU2; RooAbsPdf *thisPdfMCU2;
  RooAbsReal *thisCdfDataU2; RooAbsReal *thisCdfMCU2;
  
  std::stringstream name;
  
  if(!doDiago){
    // std::cout << "not diagonali" << std::endl;
    name << "sig_" << iBin;
    thisPdfDataU1 = rooWData[0]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin;
    thisPdfMCU1 = rooWMC[0]->pdf(name.str().c_str()); name.str("");

    // std::cout << "get cdfs" << std::endl;
    name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfDataU1 = rooWData[0]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU1 = rooWMC[0]->function(name.str().c_str()); name.str("");

    name << "sig_" << iBin;
    thisPdfDataU2 = rooWData[1]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin;
    thisPdfMCU2 = rooWMC[1]->pdf(name.str().c_str()); name.str("");

    name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfDataU2 = rooWData[1]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU2 = rooWMC[1]->function(name.str().c_str()); name.str("");
  } else {
    // std::cout << "diago fuck" << std::endl;
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfDataU1 = rooWDataDiag[0]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfMCU1 = rooWMCDiag[0]->pdf(name.str().c_str()); name.str("");

    // std::cout << "get cdfs" << std::endl;
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfDataU1 = rooWDataDiag[0]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU1 = rooWMCDiag[0]->function(name.str().c_str()); name.str("");
    



    // ---------- Diagonalized PDFs
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfDataU2 = rooWDataDiag[1]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfMCU2 = rooWMCDiag[1]->pdf(name.str().c_str()); name.str("");

    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfDataU2 = rooWDataDiag[1]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU2 = rooWMCDiag[1]->function(name.str().c_str()); name.str("");
  }


  if(!dokeys && !doDiago) { // central values
    name << "sig_" << iBin;
    thisPdfMCU1toCorr = rooWMCtoCorr[0]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU1toCorr = rooWMCtoCorr[0]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin;
    thisPdfMCU2toCorr = rooWMCtoCorr[1]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU2toCorr = rooWMCtoCorr[1]->function(name.str().c_str()); name.str("");
  } else if (!dokeys) { // Doing the Diagonalized recoil corrections
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfMCU1toCorr = rooWMCtoCorrDiag[0]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU1toCorr = rooWMCtoCorrDiag[0]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfMCU2toCorr = rooWMCtoCorrDiag[1]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU2toCorr = rooWMCtoCorrDiag[1]->function(name.str().c_str()); name.str("");
  }else { // RooKeys recoils
    name << "key_" << iBin;
    thisPdfMCU1toCorr = rooWMCtoCorr[0]->pdf(name.str().c_str()); name.str("");
    name << "key_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU1toCorr = rooWMCtoCorr[0]->function(name.str().c_str()); name.str("");
    name << "key_" << iBin;
    thisPdfMCU2toCorr = rooWMCtoCorr[1]->pdf(name.str().c_str()); name.str("");
    name << "key_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU2toCorr = rooWMCtoCorr[1]->function(name.str().c_str()); name.str("");

  }

  std::stringstream varName;
  varName.str("");varName << "u_"<<iBin;
  RooRealVar *myXdU1,  *myXmU1, *myXmcU1, *myXdU2, *myXmU2, *myXmcU2;
  if(!doDiago) { // normal corrections
     myXdU1 =  (RooRealVar*) rooWData[0]->var(varName.str().c_str());
     myXmU1 =  (RooRealVar*) rooWMC[0]->var(varName.str().c_str());
     myXmcU1 =  (RooRealVar*) rooWMCtoCorr[0]->var(varName.str().c_str());
     myXdU2 =  (RooRealVar*) rooWData[1]->var(varName.str().c_str());
     myXmU2 =  (RooRealVar*) rooWMC[1]->var(varName.str().c_str());
     myXmcU2 =  (RooRealVar*) rooWMCtoCorr[1]->var(varName.str().c_str());
  } else {// diagonalized
     myXdU1 =  (RooRealVar*) rooWDataDiag[0]->var(varName.str().c_str());
     myXmU1 =  (RooRealVar*) rooWMCDiag[0]->var(varName.str().c_str());
     myXmcU1 =  (RooRealVar*) rooWMCtoCorrDiag[0]->var(varName.str().c_str());
     myXdU2 =  (RooRealVar*) rooWDataDiag[1]->var(varName.str().c_str());
     myXmU2 =  (RooRealVar*) rooWMCDiag[1]->var(varName.str().c_str());
     myXmcU2 =  (RooRealVar*) rooWMCtoCorrDiag[1]->var(varName.str().c_str());
  }

  // for the closure on Z events: this step should give pU1ValMzlike=pU1
  double pU1ValMzlike = triGausInvGraphPDF(pU1,iGenPt,thisCdfMCU1toCorr,thisCdfMCU1,thisPdfMCU1toCorr,thisPdfMCU1,myXmU1,myXmcU1,iBin,0);
  double pU2ValMzlike = triGausInvGraphPDF(pU2,iGenPt,thisCdfMCU2toCorr,thisCdfMCU2,thisPdfMCU2toCorr,thisPdfMCU2,myXmU2,myXmcU2,iBin,0);

  // invert the target MC (Z) to the (ZDATA)
  double pU2ValDzlike = triGausInvGraphPDF(pU2ValMzlike,iGenPt,thisCdfMCU2,thisCdfDataU2,thisPdfMCU2,thisPdfDataU2,myXdU2,myXmU2,iBin,0);
  double pU1ValDzlike = triGausInvGraphPDF(pU1ValMzlike,iGenPt,thisCdfMCU1,thisCdfDataU1,thisPdfMCU1,thisPdfDataU1,myXdU1,myXmU1,iBin,0);

  // have the newW recoil as WrecoilMC + Difference in Zdata/MC
  pU1   = pU1 + ( pU1ValDzlike - pU1ValMzlike);
  pU2   = pU2 + ( pU2ValDzlike - pU2ValMzlike);
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  
  iU1   = pU1; 
  iU2   = pU2;

  return;
}


// What the hell is this? Delete it //////////////////////
void RecoilCorrector::metDistributionShitty(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       double &iU1,double &iU2,double &evtWeight) {
  
  //  double pDefU1    = iU1Default->Eval(iGenPt);
  // std::cout << "hi 2" << std::endl;
  double iGenPt2 = 0;
  Int_t nbinsPt = vZPtBins.size()-1;
  int iBin = -1;
  for(int i = 0; i < nbinsPt-1; ++i){
    if(iGenPt > vZPtBins[nbinsPt]){
      iBin = nbinsPt-1;
      iGenPt2 = (vZPtBins[nbinsPt-1]+vZPtBins[nbinsPt-2])*0.5;
      break;
    }
    if(vZPtBins[i+1] < iGenPt) continue;
    if(vZPtBins[i] > iGenPt ) continue;
    if(iGenPt < vZPtBins[i+1] && vZPtBins[i] < iGenPt){
      iBin = i; 
      iGenPt2 = (vZPtBins[i+1]+vZPtBins[i])*0.5;
      break;
    }
  }
  if(iBin<0) return;
  // std::cout << std::endl;
      // std::cout << "ibin = " << iBin <<  " pt = " <<  iGenPt <<  "  u1 = " << iU1 << std::endl;
  double pUX  = iMet*cos(iMPhi) + iLepPt*cos(iLepPhi);
  double pUY  = iMet*sin(iMPhi) + iLepPt*sin(iLepPhi);
  double pU   = sqrt(pUX*pUX+pUY*pUY);
  double pCos = - (pUX*cos(iGenPhi) + pUY*sin(iGenPhi))/pU;
  double pSin =   (pUX*sin(iGenPhi) - pUY*cos(iGenPhi))/pU;
  double pU1  = pU*pCos; // U1 in sample to Correct (WMC or ZMC)
  double pU2  = pU*pSin; // U2 in sample to Correct (WMC or ZMC)

  double iU1Bin = 0;
  if(iBin > 51) return;
    for(int i = 0; i <= hRatiosU1[iBin]->GetNbinsX();++i){
        if(iU1 > hRatiosU1[iBin]->GetBinLowEdge(i) && iU1 < hRatiosU1[iBin]->GetBinLowEdge(i+1)){ iU1Bin = i; break; }
	}
    // std::cout << "u1Bin = " << iU1Bin << std::endl;
  evtWeight=hRatiosU1[iBin]->GetBinContent(iU1Bin);
  // pU1*=hRecoilU2[iBin]->GetBinContent(iU2Bin);
                       if(evtWeight > 2.0 || evtWeight < 0) {evtWeight=1;}
  
  // std::cout << "ibin = " << iBin << "  u1Bin = " << iU1Bin << "  evtWeight = " << evtWeight << std::endl;

  // pU1   = pU1 + ( pU1ValDzlike - pU1ValMzlike);
  // pU2   = pU2 + ( pU2ValDzlike - pU2ValMzlike);
  iMet  = iMet;//calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = iMPhi;//calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  
  // std::cout << "u1fin = " << pU1 << ", u2fin = " << pU2 << ", met " << iMet << ", metPhi " << iMPhi  << std::endl;
  iU1   = iU1; 
  iU2   = iU2;

  return;
} ////////////////////// Delete ^


// Also not used ever? Delte
void RecoilCorrector::metDistributionFromToys(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       double &iU1,double &iU2,double iFluc,double iScale) {
  
  double iGenPt2 = 0;
  Int_t nbinsPt = vZPtBins.size()-1;
  int iBin = 0;
  for(int i = 0; i < nbinsPt-1; ++i){
    if(iGenPt > vZPtBins[nbinsPt]){
      iBin = nbinsPt-1;
      iGenPt2 = (vZPtBins[nbinsPt-1]+vZPtBins[nbinsPt-2])*0.5;
      break;
    }
    if(vZPtBins[i+1] < iGenPt) continue;
    if(vZPtBins[i] > iGenPt ) continue;
    if(iGenPt < vZPtBins[i+1] && vZPtBins[i] < iGenPt){
      iBin = i; 
      iGenPt2 = (vZPtBins[i+1]+vZPtBins[i])*0.5;
      break;
    }
  }
  
  double pU1 = 0; 
  double pU2 = 0;


  std::stringstream varName;
  varName.str("");varName << "u_"<<iBin;
  RooRealVar* myXdU1 =  (RooRealVar*) rooWData[0]->var(varName.str().c_str());
  RooRealVar* myXmU1 =  (RooRealVar*) rooWMC[0]->var(varName.str().c_str());
  RooRealVar* myXdU2 =  (RooRealVar*) rooWData[1]->var(varName.str().c_str());
  RooRealVar* myXmU2 =  (RooRealVar*) rooWMC[1]->var(varName.str().c_str());
  // std::cout << "blah " << std::endl;
  
  std::stringstream name;
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfDataU1 = rooWData[0]->function(name.str().c_str()); name.str("");

  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfDataU2 = rooWData[1]->function(name.str().c_str()); name.str("");

  double r1 = ((double) rand() / (RAND_MAX));
  double r2 = ((double) rand() / (RAND_MAX));
  pU1=thisCdfDataU1->findRoot(*myXdU1,myXdU1->getMin(),myXdU1->getMax(),r1);
  pU2=thisCdfDataU2->findRoot(*myXdU2,myXdU2->getMin(),myXdU2->getMax(),r2);

  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  
  iU1   = pU1; 
  iU2   = pU2;
  
  return;
}
// ----------------------------------------------------------------------------------------------------------------------------------------
// -------------- Utility functions to do basic calculations ---------------
//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::calculate(int iMet,double iEPt,double iEPhi,double iWPhi,double iU1,double iU2) { 
  double lMX = -iEPt*cos(iEPhi) - iU1*cos(iWPhi) + iU2*sin(iWPhi);
  double lMY = -iEPt*sin(iEPhi) - iU1*sin(iWPhi) - iU2*cos(iWPhi);
  if(iMet == 0) return sqrt(lMX*lMX + lMY*lMY);
  if(iMet == 1) {if(lMX > 0) {return atan(lMY/lMX);} return (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX); } 
  if(iMet == 2) return lMX;
  if(iMet == 3) return lMY;
  return lMY;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::getCorError2(double iVal,TF1 *iFit) { 
  double lE = sqrt(iFit->GetParError(0))  + iVal*sqrt(iFit->GetParError(2));
  if(fabs(iFit->GetParError(4)) > 0) lE += iVal*iVal*sqrt(iFit->GetParError(4));
  return lE*lE;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::mag(double iV0,double iV1,double iV2,double iV3) { 
  return sqrt(iV0+iV1*iV1+2*iV1*0.88 + iV2*iV2+2.*iV2*0.88+ iV3*iV3+2.*iV3*0.88);//
}

//-----------------------------------------------------------------------------------------------------------------------------------------
double RecoilCorrector::correlatedSeed(double iVal, double iCorr1,double iCorr2,double iCorr3,double iSeed0,double iSeed1,double iSeed2,double iSeed3) { 
  double lMag = mag(1.,iCorr1,iCorr2,iCorr3); 
  double lVal = ((1./lMag) + (iCorr1/lMag)*(iSeed1) + (iCorr2/lMag)*fabs(iSeed2) + (iCorr3/lMag)*fabs(iSeed3))*iVal;
  lVal*=iSeed0;
  return lVal;
}
double RecoilCorrector::deCorrelate(double iVal, double iCorr1,double iCorr2,double iCorr3,double iSeed0,double iSeed1,double iSeed2,double iSeed3) { 
  double lMag = mag(1.,iCorr1,iCorr2,iCorr3); 
  double lVal =  (1 - iCorr1*fabs(iSeed1) - iCorr2*fabs(iSeed2) - iCorr3*fabs(iSeed3))*lMag;
  return lVal*iVal*iSeed0;
}

Double_t RecoilCorrector::dSigma(const TF1 *fcn, const Double_t x, const TFitResult *fs) {
  Double_t df[3];
  Double_t a  = fcn->GetParameter(0);
  Double_t b  = fcn->GetParameter(1);
  Double_t c  = fcn->GetParameter(2);
  
  df[0] = x*x;
  df[1] = x;
  df[2] = 1;
  
  Double_t err2=0;
  for(Int_t i=0; i<3; i++) {
    err2 += df[i]*df[i]*(fs->GetCovarianceMatrix()[i][i]);
    for(Int_t j=i+1; j<3; j++) {
      err2 += 2.0*df[i]*df[j]*(fs->GetCovarianceMatrix()[i][j]);
    }
  }
  assert(err2>=0);
  return sqrt(err2);
}

Double_t RecoilCorrector::dMean(const TF1 *fcn, const Double_t x, const TFitResult *fs) {
  Double_t df[2];
  df[0] = 1;
  df[1] = x;
  Double_t err2 = df[0]*df[0]*(fs->GetCovarianceMatrix()[0][0]) 
                  + df[1]*df[1]*(fs->GetCovarianceMatrix()[1][1]) 
          + 2.0*df[0]*df[1]*(fs->GetCovarianceMatrix()[0][1]);
  assert(err2>=0);
  return sqrt(err2);
}

Double_t RecoilCorrector::calcErrorGraph(const TGraphErrors *graph, const Double_t x, const Double_t bins[], const Int_t bin) {
  Double_t diff = (graph->GetErrorY(bin)*(bins[bin+1]-x) + graph->GetErrorY(bin+1)*(x - bins[bin]));
  Double_t avg = diff/((bins[bin+1]-bins[bin]));
  return sqrt(avg);
}


double RecoilCorrector::getError2(double iVal,TF1 *iFit) { 
  return iFit->GetParError(0);
  double lE2 = iFit->GetParError(0) + iVal*iFit->GetParError(1) + iVal*iVal*iFit->GetParError(2);
  if(fabs(iFit->GetParError(3)) > 0) lE2 += iVal*iVal*iVal*     iFit->GetParError(3);
  if(fabs(iFit->GetParError(4)) > 0) lE2 += iVal*iVal*iVal*iVal*iFit->GetParError(4);
  if(fabs(iFit->GetParError(5)) > 0 && iFit->GetParameter(3) == 0) lE2 += iVal*iVal*               iFit->GetParError(5);
  if(fabs(iFit->GetParError(5)) > 0 && iFit->GetParameter(3) != 0) lE2 += iVal*iVal*iVal*iVal*iVal*iFit->GetParError(5);
  if(fabs(iFit->GetParError(6)) > 0) lE2 += iVal*iVal*iVal*iVal*iVal*iVal*iFit->GetParError(6);
  return lE2;
}

double RecoilCorrector::getError(double iVal,TF1 *iFit,Recoil iType) {
  if(fId == 0) return sqrt(getError2(iVal,iFit));
  if(fId != 2) return sqrt(getError2(iVal,iFit)); 
  double lEW2  = getError2(iVal,iFit);
  double lEZD2 = getError2(iVal,getFunc(true ,iType));
  double lEZM2 = getError2(iVal,getFunc(false,iType));
  double lZDat = getFunc(true ,iType)->Eval(iVal);
  double lZMC  = getFunc(false,iType)->Eval(iVal);
  double lWMC  = iFit                ->Eval(iVal);
  double lR    = lZDat/lZMC;
  double lER   = lR*lR/lZDat/lZDat*lEZD2 + lR*lR/lZMC/lZMC*lEZM2;
  double lVal  = lR*lR*lEW2 + lWMC*lWMC*lER;
  return sqrt(lVal);
}

double RecoilCorrector::getErrorMean(double iVal,TF1 *iFit,Recoil iType) {
  Double_t df[2];
  df[0] = iVal;
  df[1] = 1;
  
  double lE2 = df[0]*df[0]*iFit->GetParError(0) + 2*df[0]*df[1]*iFit->GetParError(1) + df[1]*df[1]*iFit->GetParError(2);
  return sqrt(lE2);
}

double RecoilCorrector::getErrorSigma(double iVal,TF1 *iFit,Recoil iType) {
  Double_t df[3];
  df[0] = 1;
  df[1] = iVal;
  df[2] = iVal*iVal;
  
  double diag = df[0]*df[0]*iFit->GetParError(0) + df[1]*df[1]*iFit->GetParError(2) + df[2]*df[2]*iFit->GetParError(4);
  double offD = 2*(df[0]*df[1]*iFit->GetParError(1) + df[0]*df[2]*iFit->GetParError(5) + df[1]*df[2]*iFit->GetParError(3));
  return sqrt(diag+offD);
}

// ---------------------------------------------------------------------------------------------------------------------------
// -------------- Setup for the diagonalization steps for the statistical uncertainty ----------------------------------------
void RecoilCorrector::runDiago(RooWorkspace *w, RooFitResult *result, int i, RooAbsReal *&pdfUiCdf, int sigma) {
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    // std::cout << "run diago" << std::endl; 
    // std::cout << "result info " << std::endl;
    result->Print();
    std::cout << "covariance quality = " << result->covQual() << std::endl;
    std::cout << "fit status = " << result->status() << std::endl;
    // TCanvas *c = new TCanvas("c","c",800,800);
    
  char name[50];
  
  sprintf(name,"u_%i",i);
  RooRealVar* myX1 = (RooRealVar*) w->var(name);
  
  
  sprintf(name,"eig_%i",i);
  // cout << "w= " << w << " result= " << result << " fit= " << fit << " pdfUiCdf= " << pdfUiCdf << endl;
  
  // cout << "BEFORE DIAGO 1" << " fit= " << fit<< endl;
  // w->Print();
  // cout << "AFTER CALLING w->Print()"<< " fit= " << fit << endl;
  // cout << "CALLING result->Print(\"V\")" << endl;
  // result->Print("V");
  // std::cout << "start runDiago" << std::endl;
  // result->Print();
  // std::cout << "blah" << std::endl;
  PdfDiagonalizer *diago = new PdfDiagonalizer(name, w, *result);
  // std::cout << "test 6 " << std::endl;
  // char name[50];
  sprintf(name,"sig_%i",i);
  // std::cout << "name signal " << name << std::endl;
  // cout << "AFTER CALLING DIAGO 1, searching for pdf " << fit << " fit= " << fit<< endl;
  RooAddPdf* pdf_temp = (RooAddPdf*) w->pdf(name);
  // pdf_temp->Print();
  // std::cout << "first step diago" << std::endl;
  // pdf_temp->Print();
  // cout << "AFTER CALLING pdf_temp"<< " fit= " << fit << endl;
  RooAbsPdf *newpdf = diago->diagonalize(*pdf_temp);
  
  // TString nameC;
  // RooPlot *testframea = myX1->frame(); 
  // std::cout << "plot 1" << std::endl;
  // // RooPlot *testframe2 = myX2->frame(); 
  // newpdf->plotOn(testframea);
      // // mainWep->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  // nameC = Form("TEST_JUSTDIAGPDF_bin%d.png",i);
  // testframea->Draw();
  // c->SaveAs(nameC);
  // c->Clear();
  // std::cout << "diago 1 " << std::endl;
  RooAbsPdf *varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,sigma);
  // std::cout << "diago w eigenVars " << std::endl;
  // cout << "AFTER CALLING diagonalize" << endl;
  // std::cout << "new pdf  " ;
  // varpdf->Print();
  // cout << "AFTER print after CALLING diagonalize" << endl;
  
  // std::cout << "test 9 " << std::endl;
    // RooPlot *testframeb = myX1->frame(); 
    
  // std::cout << "plot 2" << std::endl;
    // // RooPlot *testframe2 = myX2->frame(); 
  // newpdf->plotOn(testframea,LineColor(kBlue));
  // varpdf->plotOn(testframeb,LineColor(kRed));
  // // mainWep->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  // nameC = Form("TEST_varPDF_bin%d.png",i);
  // testframeb->Draw();
  // c->SaveAs(nameC);
  // c->Clear();
  
  // name << ""; name << "u_" << i;
  // w->Print();
  // std::cout << "var name " << name << std::endl;
  // RooRealVar* myX1 = (RooRealVar*) rooWMCtoCorr[0]->var(name.str().c_str());
  // loop over the values on the x-axis and check if the PDF is positive
  // if it's not, replace the diagonalized PDF with the original one
  
  // std::cout << "check norm " << varpdf->getNorm() << std::endl;
  // std::cout << "check the values? " << std::endl;
      // std::cout << "logval = " << varpdf->getLogVal() << "  isfinite " << std::isfinite(varpdf->getLogVal()) << std::endl;
  for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=1){
	  myX1->setVal(i);
	  double pdfval=varpdf->getVal();
      // std::cout << "logval = " << varpdf->getLogVal() << "  isfinite " << std::isfinite(varpdf->getLogVal()) << std::endl;
	  // std::cout << pdfval << std::endl;
	  if(pdfval < 0 || pdfval > 1||!std::isfinite(pdfval) || std::isnan(pdfval) ||  !std::isfinite(varpdf->getLogVal()) || std::isnan(varpdf->getLogVal())) {
		  std::cout << "pdf is messed up" << std::endl;
		  varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,0);
		  // varpdf = pdf_temp;
          break;
	  }
  }
  pdf_temp->Print();
  varpdf->Print();
  // RooRealVar* myX1=w->var("XVar");
  // pdfUiCdf = newpdf->createCdf(*myX1,RooFit::ScanAllCdf());
  // std::cout << "make a CDF" << std::endl;

  
  pdfUiCdf = varpdf->createCdf(*myX1);
  pdfUiCdf->Print();
  
          // RooPlot *testframe1 = myX1->frame(); 
        // // RooPlot *testframe2 = myX2->frame(); 
      // pdfUiCdf->plotOn(testframe1);
      // // mainWep->plotOn(wepframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
      // nameC = Form("TEST_CDF_bin%d.png",i);
      // testframe1->Draw();
      // c->SaveAs(nameC);
      // c->Clear();
  
    for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=1){
	  myX1->setVal(i);
	  double pdfval=pdfUiCdf->getVal();
      // std::cout << "logval = " << pdfUiCdf->getLogVal() << "  isfinite " << std::isfinite(pdfUiCdf->getLogVal()) << std::endl;
	  // std::cout << pdfval << std::endl;
	  if(pdfval < 0 || pdfval > 1  || std::isnan(pdfval) || !std::isfinite(pdfval)) {
		  std::cout << "cdf is messed up" << std::endl;
		  varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,0);
          pdfUiCdf = varpdf->createCdf(*myX1);
          varpdf->Print();
          pdfUiCdf->Print();
		  break;
	  }
  }
  pdfUiCdf->Print();
  // std::cout << "import shit to wksp" << std::endl;
  w->import(*varpdf, RooFit::RecycleConflictNodes(),RooFit::Silence());
  w->import(*pdfUiCdf, RooFit::RecycleConflictNodes(),RooFit::Silence());
  
  // w->Print();
  std::cout << "done now " << std::endl;

  // delete c;
  return;
}


// Not sure this even works, I think I might just give up on developing it, also delete possibly
// do shape template by using the PDF values & resetting them to be +/-1 sigma given the roofit results
void RecoilCorrector::statUnc50nsStyle(RooWorkspace *w, int i, RooAbsReal *&pdfUiCdf, int sigma) {
  char name[50];
  sprintf(name,"eig_%i",i);
  sprintf(name,"sig_%i",i);
  RooAddPdf* pdf_temp = (RooAddPdf*) w->pdf(name);
  sprintf(name,"var_%i",i);
  RooAddPdf *varpdf = new RooAddPdf(*pdf_temp,name);

  varpdf->Print();
  RooArgList pdfList = varpdf->pdfList();
  RooArgList coefList = varpdf->coefList();
  
  pdfList.Print();
  coefList.Print();
  std::cout << "first" << std::endl;
  std::cout << "val 3 " << ((RooRealVar*)coefList.at(0))->getVal() << "  val2 " << ((RooRealVar*)coefList.at(1))->getVal() << std::endl;

  RooGaussian *gaus1 = (RooGaussian*) pdfList.at(2);
  RooGaussian *gaus2 = (RooGaussian*) pdfList.at(1);
  RooGaussian *gaus3 = (RooGaussian*) pdfList.at(0);
  gaus1->Print();
  // sprintf(name,"");
  RooRealVar mean1("mean1","mean1",0); RooRealVar sigma1("sigma1","sigma1",0);
  RooArgSet tempSet1(mean1,sigma1);
  
  RooRealVar mean2("mean2","mean2",0); RooRealVar sigma2("sigma2","sigma2",0);
  RooArgSet tempSet2(mean2,sigma2);
  
  RooRealVar mean3("mean3","mean3",0); RooRealVar sigma3("sigma3","sigma3",0);
  RooArgSet tempSet3(mean3,sigma3);
  // RooRealVar mean("mean1","mean",0);
  // RooRealVar mean("mean","mean",0);
  RooArgSet* gaus1Set = (RooArgSet*) gaus1->getParameters(tempSet1); RooArgList gaus1List(*gaus1Set);
  RooArgSet* gaus2Set = (RooArgSet*) gaus2->getParameters(tempSet2); RooArgList gaus2List(*gaus2Set);
  RooArgSet* gaus3Set = (RooArgSet*) gaus3->getParameters(tempSet3); RooArgList gaus3List(*gaus3Set);
  gaus1->Print();
  gaus2->Print();
  gaus3->Print();
  ((RooRealVar*)gaus1List.at(0))->setVal(((RooRealVar*)gaus1List.at(0))->getVal()+sigma*((RooRealVar*)gaus1List.at(0))->getError());
  ((RooRealVar*)gaus1List.at(1))->setVal(((RooRealVar*)gaus1List.at(1))->getVal()+sigma*((RooRealVar*)gaus1List.at(1))->getError());
  ((RooRealVar*)gaus2List.at(0))->setVal(((RooRealVar*)gaus2List.at(0))->getVal()+sigma*((RooRealVar*)gaus2List.at(0))->getError());
  ((RooRealVar*)gaus2List.at(1))->setVal(((RooRealVar*)gaus2List.at(1))->getVal()+sigma*((RooRealVar*)gaus2List.at(1))->getError());
  ((RooRealVar*)gaus3List.at(0))->setVal(((RooRealVar*)gaus3List.at(0))->getVal()+sigma*((RooRealVar*)gaus3List.at(0))->getError());
  ((RooRealVar*)gaus3List.at(1))->setVal(((RooRealVar*)gaus3List.at(1))->getVal()+sigma*((RooRealVar*)gaus3List.at(1))->getError());
  std::cout << "sigma 1 valu3 " << ((RooRealVar*)gaus1List.at(1))->getVal() << std::endl;
  std::cout << "sigma 1 err up " << ((RooRealVar*)gaus1List.at(1))->getError() << std::endl;
  std::cout << "sigma 1 err up " << ((RooRealVar*)gaus1List.at(1))->getAsymErrorHi() << std::endl;
  std::cout << "sigma 1 err lo " << ((RooRealVar*)gaus1List.at(1))->getAsymErrorLo() << std::endl;
  gaus1->Print();
  gaus2->Print();
  gaus3->Print();
  std::cout << "new pdf  " ;
  varpdf->Print();
  
  std::cout << "orig pdf " << std::endl;
  pdf_temp->Print();
  sprintf(name,"u_%i",i);
  RooRealVar* myX1 = (RooRealVar*) w->var(name);
  // loop over the values on the x-axis and check if the PDF is positive
  // if it's not, replace the diagonalized PDF with the original one
  
      // std::cout << "logval = " << varpdf->getLogVal() << "  isfinite " << std::isfinite(varpdf->getLogVal()) << std::endl;
  for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=5){
	  myX1->setVal(i);
	  double pdfval=varpdf->getVal();
      // std::cout << "logval = " << varpdf->getLogVal() << "  isfinite " << std::isfinite(varpdf->getLogVal()) << std::endl;
	  // std::cout << pdfval << std::endl;
	  if(pdfval <= 0 || pdfval < -100 /*|| !std::isfinite(varpdf->getLogVal()) */|| std::isnan(varpdf->getLogVal())) {
		  std::cout << "pdf is messed up" << std::endl;
		  // varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,0);
		  break;
	  }
  }
  // std::cout << "nearly done "<< std::endl;
  // pdf_temp->Print();
  // varpdf->Print();
  // RooRealVar* myX1=w->var("XVar");
  // pdfUiCdf = newpdf->createCdf(*myX1,RooFit::ScanAllCdf());
  pdfUiCdf = varpdf->createCdf(*myX1);
  pdfUiCdf->Print();
  
    for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=5){
	  myX1->setVal(i);
	  double pdfval=pdfUiCdf->getVal();
      // std::cout << "logval = " << pdfUiCdf->getLogVal() << "  isfinite " << std::isfinite(pdfUiCdf->getLogVal()) << std::endl;
	  // std::cout << pdfval << std::endl;
	  if(pdfval <= 0 || pdfval < -100 ) {
		  std::cout << "cdf is messed up" << std::endl;
		  // varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,0);
		  break;
	  }
  }
  
  w->import(*varpdf, RooFit::RecycleConflictNodes(),RooFit::Silence());
  w->import(*pdfUiCdf, RooFit::RecycleConflictNodes(),RooFit::Silence());
  
  return;

}
