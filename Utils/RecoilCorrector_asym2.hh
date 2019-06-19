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
    // these are used
  void loadRooWorkspacesMCtoCorrectKeys(string iNameFile);
  void loadRooWorkspacesMCtoCorrect(string iNameFile);
  void loadRooWorkspacesMC(string iNameFile);
  void loadRooWorkspacesDiagMCtoCorrect(string iNameFile, int sigma);
  void loadRooWorkspacesDiagMC(string iNameFile, int sigma);
  void loadRooWorkspacesDiagData(string iNameFile, int sigma);
  void loadRooWorkspacesData(string iNameFile);
  // what's this? 
  void loadFileRatio(string iNameFile);
  
  void CorrectInvCdf(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0, bool dokeys=false, bool doDiago=false);

  // deletable?
  // delte?
  void addDataFile(std::string iNameDat);
  void addMCFile  (std::string iNameMC);
  void addFileWithGraph  (std::string iNameMC);
  Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResult  *fs);
  Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResult *fs);
  TFitResult *fitresPFu1mean, *fitresPFu1sigma1, *fitresPFu1sigma2,  *fitresPFu1sigma0;
  TFitResult *fitresPFu2mean, *fitresPFu2sigma1, *fitresPFu2sigma2,  *fitresPFu2sigma0;
  
  // Set up the diagonalized PDF
  void runDiago(RooWorkspace *w, RooFitResult *result,int i, RooAbsReal *&pdfUiCdf, int sigma);
  
  // delete?
  void statUnc50nsStyle(RooWorkspace *w,  int i, RooAbsReal *&pdfUiCdf, int sigma);
  
protected:
  enum Recoil {  // delete?
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
  
  void metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                double &iU1, double &iU2,double iFluc=0,double iScale=0);
                
  Double_t calcErrorGraph(const TGraphErrors *graph, const Double_t x, const Double_t bins[], const Int_t bin);

  double diGausPVal    (double iVal, double iFrac,double iSimga1,double iSigma2);
  double triGausPVal    (double iVal, double iFrac2,double iFrac3,double iSimga1,double iSigma2,double iSigma3);
  double diGausPInverse(double iPVal,double iFrac,double iSigma1,double iSigma2);
  
  // rename this b/c it's a stupid name
  double invertCDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *varDat, RooRealVar *varMC, int bin,double max);
  
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
  
  // bins for the low pileup corrections
  std::vector<double> vZPtBins = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,85,90,95,100,120,140,160,180,200,220,250,300}; // the regular 13 tev 2017 one
  
  // std::vector<double> vZPtBins = {0,2.5,5.0,10,20,30,40,50,60,80,100,125,150,200,250,300}; // rebin to check effect of bin size on statitiscal unc (2017, 13 TeV)

  int nBins = vZPtBins.size()-1;
  // TH1D **hRatiosU1 = new TH1D*[nBins];
  // TH1D **hRatiosU2 = new TH1D*[nBins];

};

//-----------------------------------------------------------------------------------------------------------------------------------------
// constructors, but some of the parts are pointless lol

RecoilCorrector::RecoilCorrector(string iNameZDat,std::string iPrefix, int iSeed) {
}

RecoilCorrector::RecoilCorrector(string iNameZ, int iSeed) {
}

//-----------------------------------------------------------------------------------------------------------------------------------------

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

void RecoilCorrector::CorrectInvCdf(double &met, double &metphi, double lGenPt, double lGenPhi, double lepPt, double lepPhi,double &iU1,double &iU2,double iFluc,double iScale,int njet, bool useKeys, bool diago) {
  dokeys=useKeys;
  doDiago=diago;
  fJet = njet; if(njet > 2) fJet = 2;
  metDistributionInvCdf(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi, iU1,iU2,iFluc,iScale);
}


double RecoilCorrector::invertCDF(double iPVal, double Zpt, RooAbsReal *pdfMCcdf, RooAbsReal *pdfDATAcdf, RooAbsPdf *wMC, RooAbsPdf *wDATA, RooRealVar *myXd, RooRealVar *myXm, int bin, double max) {
  if(iPVal< myXm->getMin()) return iPVal;
  if(iPVal> myXm->getMax()) return iPVal;
  myXm->setVal(iPVal);
  myXd->setVal(iPVal);
  pdfDATAcdf->getVal(); // do not delete this line, for some reason is necessary when using the diagonalized cdfs for the statistical uncertainty
  double pVal=pdfDATAcdf->findRoot(*myXd,myXd->getMin(),myXd->getMax(),pdfMCcdf->getVal());
  myXd->setVal(pVal);
  return pVal;
}

// The one we actually use, clean it up a bit
void RecoilCorrector::metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                       double iLepPt,double iLepPhi,
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
  double pU1ValMzlike = invertCDF(pU1,iGenPt,thisCdfMCU1toCorr,thisCdfMCU1,thisPdfMCU1toCorr,thisPdfMCU1,myXmU1,myXmcU1,iBin,0);
  double pU2ValMzlike = invertCDF(pU2,iGenPt,thisCdfMCU2toCorr,thisCdfMCU2,thisPdfMCU2toCorr,thisPdfMCU2,myXmU2,myXmcU2,iBin,0);

  // invert the target MC (Z) to the (ZDATA)
  double pU2ValDzlike = invertCDF(pU2ValMzlike,iGenPt,thisCdfMCU2,thisCdfDataU2,thisPdfMCU2,thisPdfDataU2,myXdU2,myXmU2,iBin,0);
  double pU1ValDzlike = invertCDF(pU1ValMzlike,iGenPt,thisCdfMCU1,thisCdfDataU1,thisPdfMCU1,thisPdfDataU1,myXdU1,myXmU1,iBin,0);

  // have the newW recoil as WrecoilMC + Difference in Zdata/MC
  pU1   = pU1 + ( pU1ValDzlike - pU1ValMzlike);
  pU2   = pU2 + ( pU2ValDzlike - pU2ValMzlike);
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  
  iU1   = pU1; 
  iU2   = pU2;

  return;
}

// ----------------------------------------------------------------------------------------------------------------------------------------
// -------------- Utility functions to do basic calculations ---------------
//-----------------------------------------------------------------------------------------------------------------------------------------
// A couple of these can also be deleted, need to check
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
  result->Print();
  std::cout << "Covariance Matrix quality (3 is Good) = " << result->covQual() << std::endl;
  char name[50];
  
  sprintf(name,"u_%i",i);
  RooRealVar* myX1 = (RooRealVar*) w->var(name);
  
  
  sprintf(name,"eig_%i",i);
  PdfDiagonalizer *diago = new PdfDiagonalizer(name, w, *result);
  sprintf(name,"sig_%i",i);
  RooAddPdf* pdf_temp = (RooAddPdf*) w->pdf(name);
  RooAbsPdf *newpdf = diago->diagonalize(*pdf_temp);

  RooAbsPdf *varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,sigma);
  for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=1){
	  myX1->setVal(i);
	  double pdfval=varpdf->getVal();
	  if(pdfval < 0 || pdfval > 1||!std::isfinite(pdfval) || std::isnan(pdfval) ||  !std::isfinite(varpdf->getLogVal()) || std::isnan(varpdf->getLogVal())) {
		  std::cout << "PDF us unphysical. Replacing with original shape." << std::endl;
		  varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,0);
      break;
	  }
  }
  pdfUiCdf = varpdf->createCdf(*myX1);
  
  for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=1){
	  myX1->setVal(i);
	  double pdfval=pdfUiCdf->getVal();
	  if(pdfval < 0 || pdfval > 1  || std::isnan(pdfval) || !std::isfinite(pdfval)) {
		  std::cout << "CDF is unphysical. Using the original shape." << std::endl;
		  varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,0,0);
      pdfUiCdf = varpdf->createCdf(*myX1);
		  break;
	  }
  }
  w->import(*varpdf, RooFit::RecycleConflictNodes(),RooFit::Silence());
  w->import(*pdfUiCdf, RooFit::RecycleConflictNodes(),RooFit::Silence());

  return;
}