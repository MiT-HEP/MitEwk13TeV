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
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include <TFitResult.h>

#include "TStopwatch.h"

//
// ** apply phil's recoil corrections **
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

TStopwatch sw;

class RecoilCorrector
{
  
public:
  RecoilCorrector(string iNameZDat, int iSeed=0xDEADBEEF);
  RecoilCorrector(string iNameZDat1, string iPrefix, int iSeed=0xDEADBEEF);
    
  void loadRooWorkspacesMC(string iNameFile);
  void loadRooWorkspacesData(string iNameFile);
  
  void CorrectAll(double &met, double &metphi, double iGenPt, double iGenPhi, double iLepPt, double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectType0(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectType2(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectType2FromGraph(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectInvCdf(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void CorrectFromToys(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0);
  void addDataFile(std::string iNameDat);
  void addMCFile  (std::string iNameMC);
  void addFileWithGraph  (std::string iNameMC);
  void addFileWithGraphAsym  (std::string iNameMC);
  void addMCFileAsym(std::string iNameDat);
  Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResult  *fs);
  Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResult *fs);
  TFitResult *fitresPFu1mean, *fitresPFu1sigma1, *fitresPFu1sigma2,  *fitresPFu1sigma0;
  TFitResult *fitresPFu2mean, *fitresPFu2sigma1, *fitresPFu2sigma2,  *fitresPFu2sigma0;
  void MakeCDFs(Int_t nGaus = 2, bool fixedMeans = false, int nJets = 0);
  void MakeCDFsU2(Int_t nGaus = 2, bool fixedMeans = false, int nJets = 0);
  
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
          
  void readRecoilMC3(std::vector<TGraphErrors*> &iU1Fit,std::vector<TGraphErrors*> &iU1Fit2,std::vector<TGraphErrors*> &iU1Fit3,
                     std::vector<TGraphErrors*> &iU1MRMSFit,std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,std::vector<TGraphErrors*> &iU1RMS3Fit, 
                     std::vector<TGraphErrors*> &iU1Frac2Fit,  std::vector<TGraphErrors*> &iU1Frac3Fit,
                     std::vector<TGraphErrors*> &iU2Fit,
                     std::vector<TGraphErrors*> &iU2MRMSFit,std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,std::vector<TGraphErrors*> &iU2RMS3Fit,
                     std::vector<TGraphErrors*> &iU2Frac2Fit,  std::vector<TGraphErrors*> &iU2Frac3Fit,
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

  void metDistributionType2FromGraph(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                TGraphErrors *iU1Default,TGraphErrors *iU1Default2,
                TGraphErrors *iU1RZDatFit,  TGraphErrors *iU1RZMCFit,
                TGraphErrors *iU1RZDatFit2,  TGraphErrors *iU1RZMCFit2,
                TGraphErrors *iU1MSZDatFit, TGraphErrors *iU1MSZMCFit, 
                TGraphErrors *iU1S1ZDatFit, TGraphErrors *iU1S1ZMCFit, 
                TGraphErrors *iU1S2ZDatFit, TGraphErrors *iU1S2ZMCFit, 
                TGraphErrors *iU2MSZDatFit, TGraphErrors *iU2MSZMCFit,
                TGraphErrors *iU2S1ZDatFit, TGraphErrors *iU2S1ZMCFit,     
                TGraphErrors *iU2S2ZDatFit, TGraphErrors *iU2S2ZMCFit,     
                double &iU1, double &iU2,double iFluc=0,double iScale=0);
                
  void metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                TGraphErrors *iU1Default,TGraphErrors *iU1Default2,
                double &iU1, double &iU2,double iFluc=0,double iScale=0);
  
  void metDistributionFromToys(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                TGraphErrors *iU1Default,TGraphErrors *iU1Default2,
                double &iU1, double &iU2,double iFluc=0,double iScale=0);            
                
                
  Double_t calcErrorGraph(const TGraphErrors *graph, const Double_t x, const Double_t bins[], const Int_t bin);

  double diGausPVal    (double iVal, double iFrac,double iSimga1,double iSigma2);
  double triGausPVal    (double iVal, double iFrac2,double iFrac3,double iSimga1,double iSigma2,double iSigma3);
//   double diGausPValAsym(double iVal,double diff1, double iFrac,double iSimga1,double iSigma2);
  double diGausPInverse(double iPVal,double iFrac,double iSigma1,double iSigma2);
  double triGausPInverse(double iPVal,double iFrac2,double iFrac3,double iSimga1,double iSigma2,double iSigma3);
double triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *varDat, RooRealVar *varMC, int bin,double max);
  
  double prepCdfs(RooWorkspace *work[]);
  
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
  vector<TGraphErrors*> fF1U1Fit2; vector<TGraphErrors*> fF1U1Fit3;
  vector<TGraphErrors*> fF1U1RMS3Fit; 
  vector<TGraphErrors*> fF1U1Frac2; vector<TGraphErrors*> fF1U1Frac3; 
  vector<TGraphErrors*> fF1U2RMS3Fit; 
  vector<TGraphErrors*> fF1U2Frac2; vector<TGraphErrors*> fF1U2Frac3;

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
  vector<TGraphErrors*> fM1U1Fit2; vector<TGraphErrors*> fM1U1Fit3;
  vector<TGraphErrors*> fM1U1RMS3Fit; 
  vector<TGraphErrors*> fM1U1Frac2; vector<TGraphErrors*> fM1U1Frac3;
  // sigma 3 and fracs for U2
  vector<TGraphErrors*> fM1U2RMS3Fit; 
  vector<TGraphErrors*> fM1U2Frac2; vector<TGraphErrors*> fM1U2Frac3;
  
  
  vector<TGraphErrors*> fM2U1Fit; vector<TGraphErrors*> fM2U1RMSSMFit; vector<TGraphErrors*> fM2U1RMS1Fit; vector<TGraphErrors*> fM2U1RMS2Fit; 
  vector<TGraphErrors*> fM2U2Fit; vector<TGraphErrors*> fM2U2RMSSMFit; vector<TGraphErrors*> fM2U2RMS1Fit; vector<TGraphErrors*> fM2U2RMS2Fit; 

  vector<TGraphErrors*> fMT1U1Fit; vector<TGraphErrors*> fMT1U1RMSSMFit; vector<TGraphErrors*> fMT1U1RMS1Fit; vector<TGraphErrors*> fMT1U1RMS2Fit; 
  vector<TGraphErrors*> fMT1U2Fit; vector<TGraphErrors*> fMT1U2RMSSMFit; vector<TGraphErrors*> fMT1U2RMS1Fit; vector<TGraphErrors*> fMT1U2RMS2Fit; 
  vector<TGraphErrors*> fMT2U1Fit; vector<TGraphErrors*> fMT2U1RMSSMFit; vector<TGraphErrors*> fMT2U1RMS1Fit; vector<TGraphErrors*> fMT2U1RMS2Fit; 
  vector<TGraphErrors*> fMT2U2Fit; vector<TGraphErrors*> fMT2U2RMSSMFit; vector<TGraphErrors*> fMT2U2RMS1Fit; vector<TGraphErrors*> fMT2U2RMS2Fit;
  // means 2 and 3
  vector<TGraphErrors*> fMT1U1Fit2; vector<TGraphErrors*> fMT1U1Fit3;
  vector<TGraphErrors*> fMT1U1RMS3Fit;
  vector<TGraphErrors*> fMT1U1Frac2; vector<TGraphErrors*> fMT1U1Frac3;
  vector<TGraphErrors*> fMT1U2RMS3Fit;
  vector<TGraphErrors*> fMT1U2Frac2; vector<TGraphErrors*> fMT1U2Frac3;

  vector<TGraphErrors*> fF1U1U2Corr;     vector<TGraphErrors*> fF2U1U2Corr;
  vector<TGraphErrors*> fF1F2U1Corr;     vector<TGraphErrors*> fF1F2U2Corr;
  vector<TGraphErrors*> fF1F2U1U2Corr;   vector<TGraphErrors*> fF1F2U2U1Corr;

  vector<TGraphErrors*> fM1U1U2Corr;     vector<TGraphErrors*> fM2U1U2Corr;
  vector<TGraphErrors*> fM1M2U1Corr;     vector<TGraphErrors*> fM1M2U2Corr;
  vector<TGraphErrors*> fM1M2U1U2Corr;   vector<TGraphErrors*> fM1M2U2U1Corr;
  
  RooWorkspace* rooWData[2];
  RooWorkspace* rooWMC[2];
  RooWorkspace* pdfsU1zData, pdfsU2zData;
  RooWorkspace* pdfsU1zMC, pdfsU2zMC;
  RooWorkspace* pdfsU1sigMC, pdfsU2sigMC;
  int fId; int fJet;
  
  RooWorkspace rooWksDataU1;
  RooWorkspace rooWksMCU1;
  RooWorkspace rooWksDataU2;
  RooWorkspace rooWksMCU2;
  
  std::vector<double> vZPtBins;
//   Double_t vZPtBins[] = {0,1,2.5,5.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100};
// int nZPtBins = sizeof(vZPtBins)/sizeof(Double_t)-1;

};

//-----------------------------------------------------------------------------------------------------------------------------------------
RecoilCorrector::RecoilCorrector(string iNameZDat,std::string iPrefix, int iSeed) {

  fRandom = new TRandom1(iSeed);
//   readRecoilMC2(fF1U1Fit,fF1U1Fit2,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,iNameZDat,iPrefix);
  readRecoilMC3(fF1U1Fit,fF1U1Fit2,fF1U1Fit3,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U1RMS3Fit,fF1U1Frac2,fF1U1Frac3,
  fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,fF1U2RMS3Fit,fF1U2Frac2,fF1U2Frac3,iNameZDat,iPrefix);
  fId = 0; fJet = 0;
}

RecoilCorrector::RecoilCorrector(string iNameZ, int iSeed) {

  fRandom = new TRandom1(iSeed);
  // get fits for Z data
  readRecoilMC(fF1U1Fit,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit,iNameZ,"PF");
  fId = 0; fJet = 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void RecoilCorrector::loadRooWorkspacesData(std::string iFName){
  vZPtBins ={0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};
  
//   std::cout << "n Bins = " << std::endl;

//   vZPtBins ={0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300};
  
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWData[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  lFile->Delete();
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWData[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  lFile2->Delete();
//   std::cout << "vZPtBins.size() = " << vZPtBins.size() << std::endl;
  for(int i = 0; i < vZPtBins.size()-1; ++i){
//     std::cout << "iBin = " << i << std::endl;
//     std::cout << "lower,upper bound = [" << vZPtBins[i] << ", " << vZPtBins[i+1] << "]" << std::endl;
    std::stringstream name;
    name << "sig_" << i;
    RooAbsPdf* pdf1 = rooWData[0]->pdf(name.str().c_str());
    RooAbsPdf* pdf2 = rooWData[1]->pdf(name.str().c_str());
    name.str(""); name << "u_" << i;
    RooRealVar* myX1 = (RooRealVar*) rooWData[0]->var(name.str().c_str());
    RooRealVar* myX2 = (RooRealVar*) rooWData[1]->var(name.str().c_str());
    RooAbsReal *cdfU1 = pdf1->createCdf(*myX1);
    rooWData[0]->import(*cdfU1, RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWData[1]->import(*cdfU2, RooFit::Silence());
    
//     myX1->Print("v");
    
//     RooPlot* xframe = myX1->frame(); 
//     pdf1->plotOn(xframe); 
//     xframe->Draw(); 
//     iC->SaveAs(("PDF_"+name.str()+".png").c_str());
  }
//   rooWData[0]->Print();
}
void RecoilCorrector::loadRooWorkspacesMC(std::string iFName){
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMC[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  lFile->Delete();
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMC[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  lFile2->Delete();
  for(int i = 0; i < vZPtBins.size()-1; ++i){
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
//   rooWMC[0]->Print();
}

void RecoilCorrector::addDataFile(std::string iNameData) {
  readRecoil(fD1U1Fit,fD1U1RMSSMFit,fD1U1RMS1Fit,fD1U1RMS2Fit,fD1U2Fit,fD1U2RMSSMFit,fD1U2RMS1Fit,fD1U2RMS2Fit,iNameData,"fcnPF",0);
  fId++;   
}
void RecoilCorrector::addMCFile  (std::string iNameMC) {
  fId++;
  readRecoilMC(fM1U1Fit,fM1U1RMSSMFit,fM1U1RMS1Fit,fM1U1RMS2Fit,fM1U2Fit,fM1U2RMSSMFit,fM1U2RMS1Fit,fM1U2RMS2Fit,iNameMC,"fcnPF");
}
void RecoilCorrector::addMCFileAsym  (std::string iNameMC) {
  fId++;
//   readRecoilMC2(fM1U1Fit,fM1U1Fit2,fM1U1RMSSMFit,fM1U1RMS1Fit,fM1U1RMS2Fit,fM1U2Fit,fM1U2RMSSMFit,fM1U2RMS1Fit,fM1U2RMS2Fit,iNameMC,"fcnPF");
  readRecoilMC3(fM1U1Fit,fM1U1Fit2,fM1U1Fit3,fM1U1RMSSMFit,fM1U1RMS1Fit,fM1U1RMS2Fit,fM1U1RMS3Fit,fM1U1Frac2,fM1U1Frac3,fM1U2Fit,fM1U2RMSSMFit,fM1U2RMS1Fit,fM1U2RMS2Fit,fM1U2RMS3Fit,fM1U2Frac2,fM1U2Frac3,iNameMC,"fcnPF");
  

  
}

void RecoilCorrector::addFileWithGraph  (std::string iNameMC) { 
  fId++;
  readRecoilMC(fMT1U1Fit,fMT1U1RMSSMFit,fMT1U1RMS1Fit,fMT1U1RMS2Fit,fMT1U2Fit,fMT1U2RMSSMFit,fMT1U2RMS1Fit,fMT1U2RMS2Fit,iNameMC,"fcnPF");
}

void RecoilCorrector::addFileWithGraphAsym  (std::string iNameMC) {
  fId++;
//   readRecoilMC2(fMT1U1Fit,fMT1U1Fit2,fMT1U1RMSSMFit,fMT1U1RMS1Fit,fMT1U1RMS2Fit,fMT1U2Fit,fMT1U2RMSSMFit,fMT1U2RMS1Fit,fMT1U2RMS2Fit,iNameMC,"fcnPF");
  readRecoilMC3(fMT1U1Fit,fMT1U1Fit2,fMT1U1Fit3,fMT1U1RMSSMFit,fMT1U1RMS1Fit,fMT1U1RMS2Fit,fMT1U1RMS3Fit,fMT1U1Frac2,fMT1U1Frac3,fMT1U2Fit,fMT1U2RMSSMFit,fMT1U2RMS1Fit,fMT1U2RMS2Fit,fMT1U2RMS3Fit,fMT1U2Frac2,fMT1U2Frac3,iNameMC,"fcnPF");
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
//   std::cout << "hello" << std::endl;
  
  
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

void RecoilCorrector::CorrectType2FromGraph(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet) {  
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0; 

  metDistributionType2FromGraph(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi,
               fF1U1Fit[fJet],fF1U1Fit2[fJet],//fF1U1Fit3[fJet],
               fMT1U1Fit     [fJet],fM1U1Fit     [fJet],
               fMT1U1Fit2    [fJet],fM1U1Fit2    [fJet],
//                fMT1U1Fit3    [fJet],fM1U1Fit3    [fJet],
               fMT1U1RMSSMFit[fJet],fM1U1RMSSMFit[fJet],
               fMT1U1RMS1Fit [fJet],fM1U1RMS1Fit [fJet],
               fMT1U1RMS2Fit [fJet],fM1U1RMS2Fit [fJet],
//                fMT1U1RMS3Fit [fJet],fM1U1RMS3Fit [fJet],
               fMT1U2RMSSMFit[fJet],fM1U2RMSSMFit[fJet],
               fMT1U2RMS1Fit [fJet],fM1U2RMS1Fit [fJet],
               fMT1U2RMS2Fit [fJet],fM1U2RMS2Fit [fJet], 
               iU1,iU2,iFluc,iScale);
}

void RecoilCorrector::CorrectInvCdf(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet) {  
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0; 

  metDistributionInvCdf(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi,
               fF1U1Fit[fJet],fF1U1Fit2[fJet],
               iU1,iU2,iFluc,iScale);
}

void RecoilCorrector::CorrectFromToys(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet) {  
  fJet = njet; if(njet > 2) fJet = 2;
  if(fJet >= int(fF1U1Fit.size())) fJet = 0; 

  metDistributionFromToys(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi,
               fF1U1Fit[fJet],fF1U1Fit2[fJet],
               iU1,iU2,iFluc,iScale);
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

// void RecoilCorrector::readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,
//            std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit,
//            std::string iFName,std::string iPrefix) {
//   TFile *lFile  = new TFile(iFName.c_str());
//   iU1Fit.push_back    ( (TF1*) lFile->FindObjectAny((iPrefix+"u1Mean_0").c_str()));
//   iU1MRMSFit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1MeanRMS_0").c_str()));
//   iU1RMS1Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1RMS1_0").c_str()));
//   iU1RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u1RMS2_0").c_str()));
//   iU2Fit    .push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2Mean_0").c_str()));
//   iU2MRMSFit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2MeanRMS_0").c_str()));
//   iU2RMS1Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS1_0").c_str()));
//   iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny((iPrefix+"u2RMS2_0").c_str()));
//   lFile->Close();
// //   return 1;
// }

//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::readRecoilMC(std::vector<TGraphErrors*> &iU1Fit,std::vector<TGraphErrors*> &iU1MRMSFit,std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,
                         std::vector<TGraphErrors*> &iU2Fit,std::vector<TGraphErrors*> &iU2MRMSFit,std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,
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
  
  // new
//   TFile *lFile  = new TFile(iFName.c_str());
//   std::cout << "file " << iFName.c_str() << std::endl;
//   int lNJet = 0;
//   std::stringstream lSS; //lSS << iPrefix;
//   iPrefix = "gr"; 
//   lSS << iPrefix << "MeanU1_0"; 
//   std::cout << "reading " << lSS.str().c_str();
//   iU1Fit.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
//   std::cout << "reading " << lSS.str().c_str();
//   lSS << iPrefix << "RMS0U1_0";
//   iU1MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
//   lSS << iPrefix << "RMS1U1_0"; 
//   iU1RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
//   lSS << iPrefix << "RMS2U1_0"; iU1RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
//   lSS << iPrefix << "MeanU2_0"    << lNJet; iU2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
//   lSS << iPrefix << "RMS0U2_0"; 
//   iU2MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
//   lSS << iPrefix << "RMS1U2_0"; iU2RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
//   lSS << iPrefix << "RMS2U2_0"; iU2RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
//   lFile->Close();
//   
}

void RecoilCorrector::readRecoilMC3(std::vector<TGraphErrors*> &iU1Fit,std::vector<TGraphErrors*> &iU1Fit2,std::vector<TGraphErrors*> &iU1Fit3,
                     std::vector<TGraphErrors*> &iU1MRMSFit,std::vector<TGraphErrors*> &iU1RMS1Fit,std::vector<TGraphErrors*> &iU1RMS2Fit,std::vector<TGraphErrors*> &iU1RMS3Fit, 
                     std::vector<TGraphErrors*> &iU1Frac2Fit,  std::vector<TGraphErrors*> &iU1Frac3Fit,
                     std::vector<TGraphErrors*> &iU2Fit,
                     std::vector<TGraphErrors*> &iU2MRMSFit,std::vector<TGraphErrors*> &iU2RMS1Fit,std::vector<TGraphErrors*> &iU2RMS2Fit,std::vector<TGraphErrors*> &iU2RMS3Fit, 
                     std::vector<TGraphErrors*> &iU2Frac2Fit,  std::vector<TGraphErrors*> &iU2Frac3Fit,
                     std::string iFName,std::string iPrefix) {
  std::cout << "reading stuff..." << std::endl;
  TFile *lFile  = new TFile(iFName.c_str());
  int lNJet = 0;
  std::stringstream lSS; //lSS << iPrefix;
  iPrefix = "grPF"; 
  lSS << iPrefix << "u1mean"; 
  iU1Fit.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1mean2"; 
  iU1Fit2.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1mean3"; 
  iU1Fit3.push_back    ( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u1sigma0";
  iU1MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma1"; 
  iU1RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma2"; iU1RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1sigma3"; iU1RMS3Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1frac2"; iU1Frac2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u1frac3"; iU1Frac3Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  
  lSS << iPrefix << "u2mean"; iU2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma0"; 
  iU2MRMSFit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma1"; iU2RMS1Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma2"; iU2RMS2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2sigma3"; iU2RMS3Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lSS << iPrefix << "u2frac2"; iU2Frac2Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
  lSS << iPrefix << "u2frac3"; iU2Frac3Fit.push_back( (TGraphErrors*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
  lFile->Close(); 
  std::cout << "done reading " << std::endl;
  
  
  
}


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


// double RecoilCorrector::triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooWorkspace *wMC, RooWorkspace *wDATA, double max) {
// double RecoilCorrector::triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, int bin, double max) {
double RecoilCorrector::triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *myXd, RooRealVar *myXm, int bin, double max) {
// double RecoilCorrector::triGausInvGraphPDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *myXd, double max) {

// std::cout << "-------" << std::endl;
// std::cout << "ipval " << iPVal << std::endl;
// std::cout << "max " << max << std::endl;
  if(TMath::Abs(iPVal-max)>=400) return iPVal;
  myXm->setVal(iPVal);
  double pVal=pdfDATAcdf->findRoot(*myXd,myXd->getMin(),myXd->getMax(),pdfMCcdf->getVal());
//   std::cout << "pVal " << pVal << std::endl;
  //if(TMath::Abs(pVal)>=max) pVal=iPVal;
  return pVal;

}

double RecoilCorrector::diGausPVal(double iVal,double iFrac,double iSigma1,double iSigma2) {
  return iFrac*TMath::Erf(iVal/iSigma1) + (1-iFrac)*TMath::Erf(iVal/iSigma2);
}

double RecoilCorrector::triGausPVal(double iVal,double iFrac2,double iFrac3,double iSigma1,double iSigma2,double iSigma3) {
  return (1-iFrac2-iFrac3)*TMath::Erf(iVal/iSigma1) + iFrac2*TMath::Erf(iVal/iSigma2)+ iFrac3*TMath::Erf(iVal/iSigma3);
//    return iFrac*TMath::Erf(iVal/iSigma1) + (1-iFrac)*TMath::Erf(iVal/iSigma2);
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

double RecoilCorrector::triGausPInverse(double iPVal,double iFrac2,double iFrac3,double iSigma1,double iSigma2,double iSigma3) {
  double lVal = TMath::ErfInverse(iPVal);
  double lMin = lVal*std::min(std::min(iSigma1,iSigma2),std::min(iSigma1,iSigma3));
  double lMax = lVal*std::max(std::max(iSigma1,iSigma2),std::max(iSigma1,iSigma3));
//   double lMin = lVal * ((iSigma1 < iSigma2) ? iSigma1 : iSigma2); // lVal * sigma1
//   double lMax = lVal * ((iSigma1 < iSigma2) ? iSigma2 : iSigma1); // lVal * sigma2
  double lDiff = (lMax-lMin); // lVal * (sigma2 - sigma1)
  //Iterative procedure to invert a double gaussian given a PVal
  int lId = 0; int lN1 = 4;  int lN2 = 10;  // Fewer toys
  //   int lId = 0; int lN1 = 10;  int lN2 = 100;  // More toys - takes longer
  for(int i0 = 0; i0 < lN1; i0++) {
    if(i0 != 0) lMin = lMin + (lId-1)*lDiff/lN2;
    if(i0 != 0) lDiff/=lN2;
    for(int i1 = 0; i1 < lN2; i1++) { 
      double pVal = lMin + lDiff/lN2*i1;
      pVal = triGausPVal(pVal,iFrac2,iFrac3,iSigma1,iSigma2,iSigma3);
      if(pVal > iPVal) {lId = i1; break;}
    }
  }
  return (lMin + (lId-0.5)*lDiff/lN2);
}


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
//   std::cout << "wtf " << pDefU1 << std::endl;
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

void RecoilCorrector::metDistributionType2FromGraph(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       TGraphErrors *iU1Default, TGraphErrors *iU1Default2,//TGraphErrors *iU1Default3,
                       TGraphErrors *iU1RZDatFit,  TGraphErrors *iU1RZMCFit,
                       TGraphErrors *iU1R2ZDatFit, TGraphErrors *iU1R2ZMCFit,
//                        TGraphErrors *iU1R3ZDatFit, TGraphErrors *iU1R3ZMCFit,
                       TGraphErrors *iU1MSZDatFit, TGraphErrors *iU1MSZMCFit, 
                       TGraphErrors *iU1S1ZDatFit, TGraphErrors *iU1S1ZMCFit, 
                       TGraphErrors *iU1S2ZDatFit, TGraphErrors *iU1S2ZMCFit, 
//                        TGraphErrors *iU1S3ZDatFit, TGraphErrors *iU1S3ZMCFit,
                       TGraphErrors *iU2MSZDatFit, TGraphErrors *iU2MSZMCFit,
                       TGraphErrors *iU2S1ZDatFit, TGraphErrors *iU2S1ZMCFit,                                             
                       TGraphErrors *iU2S2ZDatFit, TGraphErrors *iU2S2ZMCFit, 
                       double &iU1,double &iU2,double iFluc,double iScale) {
  
  // std::cout << "hello" << std::endl;
  double pDefU1    = iU1Default->Eval(iGenPt);
  double pDefU1_2    = iU1Default2->Eval(iGenPt);
//   std::cout << "pdfe1 " << pDefU1 << std::endl;
//   std::cout << "pdef2 " << pDefU1_2 << std::endl;
  double lRescale  = sqrt((TMath::Pi())/2.);             
  double pDU1       = iU1RZDatFit ->Eval(iGenPt);

  double pDFrac1    = iU1MSZDatFit->Eval(iGenPt);//*lRescale;
  double pDSigma1_1 = iU1S1ZDatFit->Eval(iGenPt);//*pDFrac1;
  double pDSigma1_2 = iU1S2ZDatFit->Eval(iGenPt);//*pDFrac1;
//   double pDFrac2    = iU2MSZDatFit->Eval(iGenPt);//*lRescale;
  double pDSigma2_1 = iU2S1ZDatFit->Eval(iGenPt);//*pDFrac2;
  double pDSigma2_2 = iU2S2ZDatFit->Eval(iGenPt);//*pDFrac2;
  
  double pDSigma2_3 = fMT1U2RMS3Fit[0]->Eval(iGenPt);
  double pDFrac2 = fMT1U2Frac2[0]->Eval(iGenPt);
  double pDFrac3 = fMT1U2Frac3[0]->Eval(iGenPt);

  double pMU1       = iU1RZMCFit  ->Eval(iGenPt);
  double pMU2       = 0; 

  double pMFrac1    = iU1MSZMCFit ->Eval(iGenPt);//*lRescale;
  double pMSigma1_1 = iU1S1ZMCFit ->Eval(iGenPt);//*pMFrac1;
  double pMSigma1_2 = iU1S2ZMCFit ->Eval(iGenPt);//*pMFrac1;
//   double pMFrac2    = iU2MSZMCFit ->Eval(iGenPt);//*lRescale;
  double pMSigma2_1 = iU2S1ZMCFit ->Eval(iGenPt);//*pMFrac2;
  double pMSigma2_2 = iU2S2ZMCFit ->Eval(iGenPt);//*pMFrac2;
  
  double pMSigma2_3 = fM1U2RMS3Fit[0]->Eval(iGenPt);
  double pMFrac2 = fM1U2Frac2[0]->Eval(iGenPt);
  double pMFrac3 = fM1U2Frac3[0]->Eval(iGenPt);
  
  double iGenPt2 = 0;
  Int_t nbinsPt = vZPtBins.size();
  int iBin = 0;
  for(int i = 0; i < nbinsPt-1; ++i){
    if(iGenPt > vZPtBins[nbinsPt-1]){
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
   // loop through Double_t vZPtBins[] to get bin of 
  }
  
//   
  if(iFluc != 0 || iScale != 0) {
    
    // get the bins for the tgraph
    Double_t ptbins[] = {0,2.5,5.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100};
//     Double_t ptbins[] ={0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100};
    
    Int_t nbins = sizeof(ptbins)/sizeof(Double_t)-1;
    
    double lEUR1    = 0;
    double lEUS1_1  = 0;
    double lEUS1_2  = 0;
    double lEU1Frac = 0;
    double lEUS2_1  = 0;
    double lEUS2_2  = 0;
    double lEU2Frac = 0;
    
//     iBin+=2*iFluc/fabs(iFluc);
//     if(iBin < 0) iBin = 0;
    
    // lookup the bin which the iGenPt belongs in
    // linear interpolate this bin and the next bin?
    for(int jBin = 0; jBin < nbins; ++jBin){
      if(iGenPt < ptbins[jBin]) continue; // skip if we are not in the right bin
      if(iGenPt > ptbins[jBin+1]) continue;
      if(iGenPt > ptbins[jBin] && iGenPt < ptbins[jBin+1]){
        lEUR1    = calcErrorGraph(iU1RZDatFit, iGenPt, ptbins, jBin); 
        lEUS1_1  = calcErrorGraph(iU1S1ZDatFit,iGenPt, ptbins, jBin); 
        lEUS1_2  = calcErrorGraph(iU1S2ZDatFit,iGenPt, ptbins, jBin);
        lEU1Frac = calcErrorGraph(iU1MSZDatFit,iGenPt, ptbins, jBin);
        lEUS2_1  = calcErrorGraph(iU2S1ZDatFit,iGenPt, ptbins, jBin);
        lEUS2_2  = calcErrorGraph(iU2S2ZDatFit,iGenPt, ptbins, jBin);
        lEU2Frac = calcErrorGraph(iU2MSZDatFit,iGenPt, ptbins, jBin);
        break; 
      }
    }
       
    pDU1       = pDU1       + iScale*lEUR1;             //Recoil
    pDFrac1    = pDFrac1    + iFluc*(lEU1Frac);        //Mean RMS 
    pDSigma1_1 = pDSigma1_1 + iFluc*lEUS1_1;//lEUS1_1*pDFrac1;    //Sigma 1 smalles sigma
    pDSigma1_2 = pDSigma1_2 + iFluc*lEUS1_2;//lEUS1_2*pDFrac1;    //Sigma 2 (Maximal when oppsite sigma 1)
    pDFrac2    = pDFrac2    + iFluc*(lEU2Frac);        //Mean RMS for U2
    pDSigma2_1 = pDSigma2_1 + iFluc*lEUS2_1;//lEUS2_1*pDFrac2;    //Sigma 1 U2
    pDSigma2_2 = pDSigma2_2 + iFluc*lEUS2_2;//(lEUS2_2)*pDFrac2;
  }
  
  pDFrac1           = (pDFrac1-pDSigma1_2)/(pDSigma1_1-pDSigma1_2);
//   pDFrac2           = (pDFrac2-pDSigma2_2)/(pDSigma2_1-pDSigma2_2);
  pMFrac1           = (pMFrac1-pMSigma1_2)/(pMSigma1_1-pMSigma1_2);
//   pMFrac2           = (pMFrac2-pMSigma2_2)/(pMSigma2_1-pMSigma2_2);

// std::cout << "Met, Phi " << iMet << " " << iMPhi << std::endl;
// std::cout << "Lep, Phi " << iLepPt << " " << iLepPhi << std::endl;
//std::cout << 
  
// calculate u1 and u2 again
  double pUX  = iMet*cos(iMPhi) + iLepPt*cos(iLepPhi);
  double pUY  = iMet*sin(iMPhi) + iLepPt*sin(iLepPhi);
  double pU   = sqrt(pUX*pUX+pUY*pUY);
  double pCos = - (pUX*cos(iGenPhi) + pUY*sin(iGenPhi))/pU;
  double pSin =   (pUX*sin(iGenPhi) - pUY*cos(iGenPhi))/pU;
  double pU1   = pU*pCos; // U1 in data
  double pU2   = pU*pSin; // U2 in data
  double pU1Diff  = pU1-pDefU1; // subtract the mean1 from MC?
  double pU1Diff2  = pU1-pDefU1_2; // subtract the mean2 also from MC?
  double pU1MeanDiff = pDefU1-pDefU1_2;
  double pU2Diff  = pU2; // don't care because expect mean to be ~0
  double p1Charge        = pU1Diff/fabs(pU1Diff);
  double p2Charge        = pU2Diff/fabs(pU2Diff);
  double pTU1Diff        = pU1Diff;
//   std::cout << "pU1 first = " << pU1 << std::endl;
//   std::cout << "pU2 first = " << pU2 << std::endl;
//   std::cout << "begin of inverting" << std::endl;
//   double pU1ValM         = diGausPValAsym(pU1Diff,pDefU1_2,pMFrac1,pMSigma1_1,pMSigma1_2);  // diffs, f1_mc, sigma1 mc, sigma2 mc
//   double pU2ValM         = diGausPVal(fabs(pU2Diff),pMFrac2,pMSigma2_1,pMSigma2_2);
//   double pU2ValM         = triGausPVal(fabs(pU2Diff),pMFrac2,pMFrac3,pMSigma2_1,pMSigma2_2,pMSigma2_3);
// //   double pU1ValD         = diGausPInverse(pU1ValM  ,pDFrac1,pDSigma1_1,pDSigma1_2); // pval, f1_data, sigma1 data, sigma2 data
//   double pU1ValD         = diGausPInverseAsym(pU1ValM  ,pU1MeanDiff, pDFrac1,pDSigma1_1,pDSigma1_2); // pval, f1_data, sigma1 data, sigma2 data, sigma1 mean? sigma2 mean??
//   double pU2ValD         = triGausPInverse(pU2ValM  ,pDFrac2,pDFrac3,pDSigma2_1,pDSigma2_2,pDSigma2_3);

// iBin = 1; // test

//   std::cout <<" hey hey hey" << iBin <<std::endl;
//   std::stringstream name;
//   name << "rooDatGausAdd_" << iBin;
//   RooAbsPdf *thisPdfDataU1 = rooWksDataU1.pdf(name.str().c_str()); name.str("");
// 
//   name << "rooMCGausAdd_" << iBin;
//   RooAbsPdf *thisPdfMCU1 = rooWksMCU1.pdf(name.str().c_str()); name.str("");
// 
// //   std::cout << "got the pdf MC" << std::endl;
//   
//   name << "rooDatGausAdd_" << iBin <<"_cdf_Int[rooU1_"<< iBin<< "_prime|CDF]_Norm[rooU1_"<< iBin<< "_prime]";
//   RooAbsReal *thisCdfDataU1 = rooWksDataU1.function(name.str().c_str()); name.str("");
//   
//   name << "rooMCGausAdd_" << iBin <<"_cdf_Int[rooU1_"<< iBin<< "_prime|CDF]_Norm[rooU1_"<< iBin<< "_prime]";
//   RooAbsReal *thisCdfMCU1 = rooWksMCU1.function(name.str().c_str()); name.str("");
// //   std::cout << "did the U1s " << std::endl;
//   name << "rooDatGausAdd_" << iBin;
//   RooAbsPdf *thisPdfDataU2 = rooWksDataU2.pdf(name.str().c_str()); name.str("");
// //   std::cout << "got the pdf data" << std::endl;
// 
//   name << "rooMCGausAdd_" << iBin;
//   RooAbsPdf *thisPdfMCU2 = rooWksMCU2.pdf(name.str().c_str()); name.str("");
//   
//   name << "rooDatGausAdd_" << iBin <<"_cdf_Int[rooU1_"<< iBin<< "_prime|CDF]_Norm[rooU1_"<< iBin<< "_prime]";
//   RooAbsReal *thisCdfDataU2 = rooWksDataU2.function(name.str().c_str()); name.str("");
//   name << "rooMCGausAdd_" << iBin <<"_cdf_Int[rooU1_"<< iBin<< "_prime|CDF]_Norm[rooU1_"<< iBin<< "_prime]";
//   RooAbsReal *thisCdfMCU2 = rooWksMCU2.function(name.str().c_str()); name.str("");
//   
//   std::stringstream varName; varName << "rooU1_"<<iBin; // not the memory leak
// 
//   RooRealVar* myXdU1 =  (RooRealVar*) rooWksDataU1.var(varName.str().c_str());
//   RooRealVar* myXmU1 =  (RooRealVar*) rooWksMCU1.var(varName.str().c_str());
// //   
//   RooRealVar* myXdU2 =  (RooRealVar*) rooWksDataU2.var(varName.str().c_str());
//   RooRealVar* myXmU2 =  (RooRealVar*) rooWksMCU2.var(varName.str().c_str());
  
//   myXdU1->setVal(pU1);
//   std::cout <<  "data cdf graph1 = "  <<thisCdfDataU1->getVal() << std::endl;
//   std::cout <<  "data pdf graph1 = "  <<thisPdfDataU1->getVal() << std::endl;
//   myXdU2->setVal(pU2);
//   std::cout <<  "data cdf graph2 = "  <<thisCdfDataU2->getVal() << std::endl;
//   std::cout <<  "data pdf graph3 = "  <<thisPdfDataU2->getVal() << std::endl;
//   myXmU1->setVal(pU1);
//   std::cout <<  "MC cdf graph1 = "  <<thisCdfMCU1->getVal() << std::endl;
//   std::cout <<  "MC pdf graph1 = "  <<thisPdfMCU1->getVal() << std::endl;
//   myXmU2->setVal(pU2);
//   std::cout <<  "MC cdf graph2 = "  <<thisCdfMCU2->getVal() << std::endl;
//   std::cout <<  "MC pdf graph2 = "  <<thisPdfMCU2->getVal() << std::endl;
//   
//   RooPlot *rooFrame = myXdU1->frame(Title(TString("Plot of MC and Data: Double-Gaussians" )));
//   thisPdfDataU1->plotOn(rooFrame,RooFit::LineColor(kRed));
// //   thisCdfDataU1->plotOn(rooFrame,RooFit::LineColor(kGreen));
//   iC->cd(); rooFrame->Draw();
//   iC->SaveAs((varName.str()+"_g.png").c_str());
//   
  std::stringstream name;
  name << "sig_" << iBin;
  RooAbsPdf *thisPdfDataU1_2 = rooWData[0]->pdf(name.str().c_str()); name.str("");

  name << "sig_" << iBin;
  RooAbsPdf *thisPdfMCU1_2 = rooWMC[0]->pdf(name.str().c_str()); name.str("");

//   std::cout << "got the pdf MC" << std::endl;
  
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfDataU1_2 = rooWData[0]->function(name.str().c_str()); name.str("");
  
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfMCU1_2 = rooWMC[0]->function(name.str().c_str()); name.str("");
//   std::cout << "did the U1s " << std::endl;
  name << "sig_" << iBin;
  RooAbsPdf *thisPdfDataU2_2 = rooWData[1]->pdf(name.str().c_str()); name.str("");
//   std::cout << "got the pdf data" << std::endl;

  name << "sig_" << iBin;
  RooAbsPdf *thisPdfMCU2_2 = rooWMC[1]->pdf(name.str().c_str()); name.str("");
  
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfDataU2_2 = rooWData[1]->function(name.str().c_str()); name.str("");
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfMCU2_2 = rooWMC[1]->function(name.str().c_str()); name.str("");
  
  std::stringstream varName;
  varName.str("");varName << "u_"<<iBin; // not the memory leak

  RooRealVar* myXdU1_2 =  (RooRealVar*) rooWData[0]->var(varName.str().c_str());
  RooRealVar* myXmU1_2 =  (RooRealVar*) rooWMC[0]->var(varName.str().c_str());
//   
  RooRealVar* myXdU2_2 =  (RooRealVar*) rooWData[1]->var(varName.str().c_str());
  RooRealVar* myXmU2_2 =  (RooRealVar*) rooWMC[1]->var(varName.str().c_str());
  
  double pU1ValD = triGausInvGraphPDF(pU1,iGenPt,thisCdfMCU1_2,thisCdfDataU1_2,thisPdfMCU1_2,thisPdfDataU1_2,myXdU1_2,myXmU1_2,iBin,pDefU1);
  double pU2ValD = triGausInvGraphPDF(pU2,iGenPt,thisCdfMCU2_2,thisCdfDataU2_2,thisPdfMCU2_2,thisPdfDataU2_2,myXdU2_2,myXmU2_2,iBin,0);
  
  
  pDefU1 *= (pDU1/pMU1);
  
  pU1   = /*pDefU1             +*/ pU1ValD;
  pU2   =                      pU2ValD;
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  

   iU1   = pU1; 
  iU2   = pU2;
 
  return;
}

void RecoilCorrector::metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       TGraphErrors *iU1Default, TGraphErrors *iU1Default2,
                       double &iU1,double &iU2,double iFluc,double iScale) {
  
  double pDefU1    = iU1Default->Eval(iGenPt);
  double pDefU1_2    = iU1Default2->Eval(iGenPt);
//   std::cout << "---------------------------" << std::endl;
//     std::cout << "original u1 = " << iU1 << std::endl;

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
  
  // calculate u1 and u2 again
  double pUX  = iMet*cos(iMPhi) + iLepPt*cos(iLepPhi);
  double pUY  = iMet*sin(iMPhi) + iLepPt*sin(iLepPhi);
  double pU   = sqrt(pUX*pUX+pUY*pUY);
  double pCos = - (pUX*cos(iGenPhi) + pUY*sin(iGenPhi))/pU;
  double pSin =   (pUX*sin(iGenPhi) - pUY*cos(iGenPhi))/pU;
  double pU1   = pU*pCos; // U1 in data
  double pU2   = pU*pSin; // U2 in data
  double pU1Diff  = pU1-pDefU1; // subtract the mean1 from MC?
  double pU1Diff2  = pU1-pDefU1_2; // subtract the mean2 also from MC?
  double pU1MeanDiff = pDefU1-pDefU1_2;
  double pU2Diff  = pU2; // don't care because expect mean to be ~0
  double p1Charge        = pU1Diff/fabs(pU1Diff);
  double p2Charge        = pU2Diff/fabs(pU2Diff);
  double pTU1Diff        = pU1Diff;

  
  
  std::stringstream name;
  name << "sig_" << iBin;
  RooAbsPdf *thisPdfDataU1 = rooWData[0]->pdf(name.str().c_str()); name.str("");
  name << "sig_" << iBin;
  RooAbsPdf *thisPdfMCU1 = rooWMC[0]->pdf(name.str().c_str()); name.str("");
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfDataU1 = rooWData[0]->function(name.str().c_str()); name.str("");
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfMCU1 = rooWMC[0]->function(name.str().c_str()); name.str("");
  name << "sig_" << iBin;
  RooAbsPdf *thisPdfDataU2 = rooWData[1]->pdf(name.str().c_str()); name.str("");
  name << "sig_" << iBin;
  RooAbsPdf *thisPdfMCU2 = rooWMC[1]->pdf(name.str().c_str()); name.str("");
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfDataU2 = rooWData[1]->function(name.str().c_str()); name.str("");
  name << "sig_" << iBin <<"_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
  RooAbsReal *thisCdfMCU2 = rooWMC[1]->function(name.str().c_str()); name.str("");
  
  std::stringstream varName;
  varName.str("");varName << "u_"<<iBin;
  RooRealVar* myXdU1 =  (RooRealVar*) rooWData[0]->var(varName.str().c_str());
  RooRealVar* myXmU1 =  (RooRealVar*) rooWMC[0]->var(varName.str().c_str());
  RooRealVar* myXdU2 =  (RooRealVar*) rooWData[1]->var(varName.str().c_str());
  RooRealVar* myXmU2 =  (RooRealVar*) rooWMC[1]->var(varName.str().c_str());
  
//   sw.Start();
  double pU1ValD = triGausInvGraphPDF(pU1,iGenPt,thisCdfMCU1,thisCdfDataU1,thisPdfMCU1,thisPdfDataU1,myXdU1,myXmU1,iBin,pDefU1);
//   sw.Stop();
//   std::cout << "Real time, inversion = " << sw.RealTime() << std::endl;
//   sw.Reset();
  double pU2ValD = triGausInvGraphPDF(pU2,iGenPt,thisCdfMCU2,thisCdfDataU2,thisPdfMCU2,thisPdfDataU2,myXdU2,myXmU2,iBin,0);

  pU1   = /*pDefU1             +*/ pU1ValD;
  pU2   =                      pU2ValD;
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
//   std::cout << "u1 = " << pU1 << std::endl;
//   std::cout << "met  = " << iMet << std::endl;
  
  iU1   = pU1; 
  iU2   = pU2;
 
  return;
}

void RecoilCorrector::metDistributionFromToys(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
                       TGraphErrors *iU1Default, TGraphErrors *iU1Default2,
                       double &iU1,double &iU2,double iFluc,double iScale) {
  
  double pDefU1    = iU1Default->Eval(iGenPt);
  double pDefU1_2    = iU1Default2->Eval(iGenPt);

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
//   double lE2 = df[0]*df[0]*iFit->GetParError(0) + df[0]*df[1]*iFit->GetParError(1) + df[1]*df[1]*iFit->GetParError(2);
  return sqrt(diag+offD);
}
