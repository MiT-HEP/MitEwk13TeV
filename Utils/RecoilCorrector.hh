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

#include "RooRealIntegral.h"
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
    // can these be condensed into something less awful
  void loadRooWorkspacesMCtoCorrectKeys(string iNameFile);
  void loadRooWorkspacesMCtoCorrect(string iNameFile);
  void loadRooWorkspacesMC(string iNameFile);
  void loadRooWorkspacesData(string iNameFile);
  void loadRooWorkspacesDiagMCtoCorrect(string iNameFile,int iPar, int sigma);
  void loadRooWorkspacesDiagMC(string iNameFile,int iPar, int sigma);
  void loadRooWorkspacesDiagData(string iNameFile,int iPar, int sigma);
  
  void CorrectInvCdf(double &pfmet, double &pfmetphi,double iGenPt,double iGenPhi,double iLepPt,double iLepPhi,double &iU1,double &iU2,double iFluc,double iScale=0,int njet=0, bool dokeys=false, bool doDiago=false);

  // Double_t dMean(const TF1 *fcn, const Double_t x, const TFitResult  *fs);
  // Double_t dSigma(const TF1 *fcn, const Double_t x, const TFitResult *fs);
  
  // Set up the diagonalized PDF
  void runDiago(RooWorkspace *w, RooFitResult *result,int i, RooAbsReal *&pdfUiCdf, int iPar, int sigma);
  
  // // delete?
  // // void statUnc50nsStyle(RooWorkspace *w,  int i, RooAbsReal *&pdfUiCdf, int sigma);
  
protected:

  void metDistributionInvCdf(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
                double iLepPt,double iLepPhi,
                double &iU1, double &iU2,double iFluc=0,double iScale=0);

  double invertCDF(double p, RooAbsReal *cdfMC, RooAbsReal *cdfData, RooRealVar *xd, RooRealVar *xm);
  
  double calculate(int iMet,double iEPt,double iEPhi,double iWPhi,double iU1,double iU2);

  
  TCanvas *iC = new TCanvas("C1","C1",800,600); 
  
  // RooWorkspace* rooWData[2],  rooWMC[2], rooWMCtoCorr[2], rooWDataDiag[2], rooWMCDiag[2], rooWMCtoCorrDiag[2];
  RooWorkspace* rooWData[2];
  RooWorkspace* rooWMC[2];
  RooWorkspace* rooWMCtoCorr[2];
  RooWorkspace* rooWDataDiag[2];
  RooWorkspace* rooWMCDiag[2];
  RooWorkspace* rooWMCtoCorrDiag[2];
  bool dokeys; bool doDiago;
  
  // Binning for the 2017G and 2017H runs
  // preserving the fine binning at low pT but the higher-pT bins (>75 GeV have been adjusted to be slightly wider)
   std::vector<double> vZPtBins ={0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000}; 

  int nBins = vZPtBins.size()-1;

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
    rooWData[0]->import(*cdfU1, RooFit::RecycleConflictNodes() ,RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWData[1]->import(*cdfU2, RooFit::RecycleConflictNodes() ,RooFit::Silence() );
    
  }
  
  std::cout << "Loaded Workspaces...DATA "<< std::endl;
}




void RecoilCorrector::loadRooWorkspacesDiagMCtoCorrect(std::string iFName,int iPar, int sigma){
  std::cout << iFName << std::endl;
  // std::cout << "aaaaaaa" << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMCtoCorrDiag[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMCtoCorrDiag[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  int iPar1=0, iPar2=0, sigma1=0, sigma2=0;
  if(iPar<6){
    iPar1=iPar;
    sigma1=sigma;
  } else {
    iPar2=iPar-6;
    sigma2=sigma;
  }
  
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    RooFitResult* fitresultU1; //= (RooFitResult*)lFile->FindObject(Form("fitResultU1_%d",i));
    lFile->GetObject(Form("fitResultU1_%d",i),fitresultU1);
    RooFitResult* fitresultU2;//    = (RooFitResult*)lFile2->FindObject(Form("fitResultU2_%d",i));
    lFile2->GetObject(Form("fitResultU2_%d",i),fitresultU2);
    TString name;
    name = Form("sig_%d",i);
    RooAbsReal *cdfU1, *cdfU2;
    runDiago(rooWMCtoCorrDiag[0],fitresultU1,i,cdfU1,iPar1,sigma1);
    runDiago(rooWMCtoCorrDiag[1],fitresultU2,i,cdfU2,iPar2,sigma2);
    
    delete fitresultU1;
    delete fitresultU2;
  }

  lFile->Delete();
  lFile2->Delete();
  std::cout << "Loaded Workspaces Source MC - Stat Unc "<< std::endl;
}

void RecoilCorrector::loadRooWorkspacesDiagMC(std::string iFName,int iPar,int sigma){
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWMCDiag[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWMCDiag[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  int iPar1=0, iPar2=0, sigma1=0, sigma2=0;
  if(iPar<6){
    iPar1=iPar;
    sigma1=sigma;
  } else {
    iPar2=iPar-6;
    sigma2=sigma;
  }
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    RooFitResult* fitresultU1; //= (RooFitResult*)lFile->FindObject(Form("fitResultU1_%d",i));
    lFile->GetObject(Form("fitResultU1_%d",i),fitresultU1);
    RooFitResult* fitresultU2;//    = (RooFitResult*)lFile2->FindObject(Form("fitResultU2_%d",i));
    lFile2->GetObject(Form("fitResultU2_%d",i),fitresultU2);
    TString name;
    name = Form("sig_%d",i);
    RooAbsReal *cdfU1, *cdfU2;
	  runDiago(rooWMCDiag[0],fitresultU1,i,cdfU1,iPar1,sigma1);
	  runDiago(rooWMCDiag[1],fitresultU2,i,cdfU2,iPar2,sigma2);
    
    delete fitresultU1;
    delete fitresultU2;
  }
  
  lFile->Delete();
  lFile2->Delete();
  std::cout << "Loaded Workspaces Z MC - Stat Unc"<< std::endl;
}

void RecoilCorrector::loadRooWorkspacesDiagData(std::string iFName,int iPar,int sigma){
  std::cout << iFName << std::endl;
  TFile *lFile  = new TFile((iFName+"pdfsU1.root").c_str());
  rooWDataDiag[0] = (RooWorkspace*) lFile->Get("pdfsU1");
  TFile *lFile2  = new TFile((iFName+"pdfsU2.root").c_str());
  rooWDataDiag[1] = (RooWorkspace*) lFile2->Get("pdfsU2");
  int iPar1=0, iPar2=0, sigma1=0, sigma2=0;
  if(iPar<6){
    iPar1=iPar;
    sigma1=sigma;
  } else {
    iPar2=iPar-6;
    sigma2=sigma;
  }
  std::cout <<iPar << " " << iPar1 << "  " << sigma1 << " " <<  iPar2 << "  " << sigma2 << std::endl;
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    RooFitResult* fitresultU1; //= (RooFitResult*)lFile->FindObject(Form("fitResultU1_%d",i));
    lFile->GetObject(Form("fitResultU1_%d",i),fitresultU1);
    RooFitResult* fitresultU2;//    = (RooFitResult*)lFile2->FindObject(Form("fitResultU2_%d",i));
    lFile2->GetObject(Form("fitResultU2_%d",i),fitresultU2);
    TString name;
    name = Form("sig_%d",i);
    RooAbsReal *cdfU1, *cdfU2;
    runDiago(rooWDataDiag[0],fitresultU1,i,cdfU1,iPar1,sigma1);
    runDiago(rooWDataDiag[1],fitresultU2,i,cdfU2,iPar2,sigma2);
    
    delete fitresultU1;
    delete fitresultU2;
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
  // cout << "what" << endl;
  for(uint i = 0; i < vZPtBins.size()-1; ++i){
    std::stringstream name;
    name << "sig_" << i;
    RooAbsPdf* pdf1 = rooWMC[0]->pdf(name.str().c_str());
    RooAbsPdf* pdf2 = rooWMC[1]->pdf(name.str().c_str());
    name.str(""); name << "u_" << i;
    RooRealVar* myX1 = (RooRealVar*) rooWMC[0]->var(name.str().c_str());
    RooRealVar* myX2 = (RooRealVar*) rooWMC[1]->var(name.str().c_str());
    RooAbsReal *cdfU1 = pdf1->createCdf(*myX1);
    rooWMC[0]->import(*cdfU1, RooFit::RecycleConflictNodes() ,RooFit::Silence() );
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWMC[1]->import(*cdfU2, RooFit::RecycleConflictNodes() ,RooFit::Silence() );
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
    rooWMCtoCorr[0]->import(*cdfU1, RooFit::RecycleConflictNodes() ,RooFit::Silence());
    RooAbsReal *cdfU2 = pdf2->createCdf(*myX2);
    rooWMCtoCorr[1]->import(*cdfU2, RooFit::RecycleConflictNodes() ,RooFit::Silence());
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
  metDistributionInvCdf(met,metphi,lGenPt,lGenPhi,lepPt,lepPhi, iU1,iU2,iFluc,iScale);
}

// this could be cleaned up a bit... i think some of the PDFs are unnecessary
double RecoilCorrector::invertCDF(double p, RooAbsReal *cdfMC, RooAbsReal *cdfData, RooRealVar *xd, RooRealVar *xm) {

  if(p < xm->getMin()) return p;
  if(p > xm->getMax()) return p;
  xm->setVal(p);
  xd->setVal(p);
  cdfData->getVal();
  double pVal=cdfData->findRoot(*xd,xd->getMin(),xd->getMax(),cdfMC->getVal());
  xd->setVal(pVal);
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

    // name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    // std::cout << "get cdfs" << std::endl;
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfDataU1 = rooWDataDiag[0]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU1 = rooWMCDiag[0]->function(name.str().c_str()); name.str("");

     // std::cout << "diago fuck2" << std::endl;
    // ---------- Diagonalized PDFs
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfDataU2 = rooWDataDiag[1]->pdf(name.str().c_str()); name.str("");
    name << "sig_" << iBin << "_eig_" << iBin;
    thisPdfMCU2 = rooWMCDiag[1]->pdf(name.str().c_str()); name.str("");
   // std::cout << "get cdfs2" << std::endl;
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfDataU2 = rooWDataDiag[1]->function(name.str().c_str()); name.str("");
    name << "sig_" << iBin <<"_eig_" << iBin << "_cdf_Int[u_"<< iBin<< "_prime|CDF]_Norm[u_"<< iBin<< "_prime]";
    thisCdfMCU2 = rooWMCDiag[1]->function(name.str().c_str()); name.str("");
  }

// cout << "bllkasdfkj" << endl;

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



  

// cout << "getVar " << endl;
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
  
  
  // cout << "plot" << endl;
  // TCanvas *c = new TCanvas("c","c",800,800);
  // RooPlot *plot3 = myXdU1->frame();
  // // thisPdfMCU1toCorr->plotOn(plot3,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),RooFit::Name("a"));
  // thisPdfDataU1->plotOn(plot3,LineColor(kBlue));
  // thisPdfMCU1->plotOn(plot3,LineColor(kGreen));
  // plot3->Draw();
  // name<<"test_u1_"<<iBin<<".png";
  // c->SaveAs(name.str().c_str());  name.str(""); 

  // thisPdfDataU2->plotOn(plot3,LineColor(kRed));
  // thisPdfMCU2->plotOn(plot3,LineColor(kBlack));
  // plot3->Draw();
  // name<<"test_"<<iBin<<".png";
  // c->SaveAs(name.str().c_str());name.str(""); 
  // cout << "done" << endl;
  
// thisPdfMCU1toCorr->Print();
// thisPdfDataU1->Print();
// thisPdfMCU1->Print();

// // myXmU1->setVal(1);
// cout << "--- " << myXmU1->getVal() << endl;;
// cout << myXmU1->getVal() << endl;;

// // cout << "get val " << thisCdfMCU1->getVal() << endl;
// // cout << "get val2 " << thisCdfMCU1toCorr->getVal() << endl;

// thisCdfMCU1toCorr->Print();
// thisCdfMCU1->Print();
// thisCdfDataU1->Print();

// cout <<" inversion " << endl;
  // for the closure on Z events: this step should give pU1ValMzlike=pU1
  // invertCDF(double p, RooAbsReal *cdfMC, RooAbsReal *cdfData, RooRealVar *xd, RooRealVar *xm)
  double pU1ValMzlike = invertCDF(pU1,thisCdfMCU1toCorr,thisCdfMCU1,myXmU1,myXmcU1);
  double pU2ValMzlike = invertCDF(pU2,thisCdfMCU2toCorr,thisCdfMCU2,myXmU2,myXmcU2);


  // invert the target MC (Z) to the (ZDATA)
  double pU2ValDzlike = invertCDF(pU2ValMzlike,thisCdfMCU2,thisCdfDataU2,myXdU2,myXmU2);
  double pU1ValDzlike = invertCDF(pU1ValMzlike,thisCdfMCU1,thisCdfDataU1,myXdU1,myXmU1);
  
// cout << "done both inversion" << endl;
  // have the newW recoil as WrecoilMC + Difference in Zdata/MC
  pU1   = pU1 + ( pU1ValDzlike - pU1ValMzlike); // should this be addition or multiplication....
  pU2   = pU2 + ( pU2ValDzlike - pU2ValMzlike);
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  // cout << "done" << endl;
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

// -------------- Setup for the diagonalization steps for the statistical uncertainty ----------------------------------------
void RecoilCorrector::runDiago(RooWorkspace *w, RooFitResult *result, int i, RooAbsReal *&pdfUiCdf, int iPar, int sigma) {
  // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  // result->Print();
  // std::cout << "Covariance Matrix quality (3 is Good) = " << result->covQual() << std::endl;
  char name[50];
  
  sprintf(name,"u_%i",i);
  RooRealVar* myX1 = (RooRealVar*) w->var(name);
  
  
  sprintf(name,"eig_%i",i);
  PdfDiagonalizer *diago = new PdfDiagonalizer(name, w, *result);
  sprintf(name,"sig_%i",i);
  RooAddPdf* pdf_temp = (RooAddPdf*) w->pdf(name);
  RooAbsPdf *newpdf = diago->diagonalize(*pdf_temp);

  RooAbsPdf *varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,iPar,sigma);
  for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=1){
	  myX1->setVal(i);
	  double pdfval=varpdf->getVal();
	  if(pdfval < 0 || pdfval > 1||!std::isfinite(pdfval) || std::isnan(pdfval) ||  !std::isfinite(varpdf->getLogVal()) || std::isnan(varpdf->getLogVal())) {
		  // std::cout << "PDF us unphysical. Replacing with original shape." << std::endl;
		  varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,iPar,0);
      break;
	  }
  }
  pdfUiCdf = varpdf->createCdf(*myX1);
  
  for(int i = floor(myX1->getMin())+1; i < floor(myX1->getMax());i+=1){
	  myX1->setVal(i);
	  double pdfval=pdfUiCdf->getVal();
	  if(pdfval < 0 || pdfval > 1  || std::isnan(pdfval) || !std::isfinite(pdfval)) {
		  // std::cout << "CDF is unphysical. Using the original shape." << std::endl;
		  varpdf = diago->diagonalizeWithEigenVariations(*newpdf,*result,iPar,0);
      pdfUiCdf = varpdf->createCdf(*myX1);
		  break;
	  }
  }
  w->import(*varpdf, RooFit::RecycleConflictNodes(),RooFit::Silence());
  w->import(*pdfUiCdf, RooFit::RecycleConflictNodes(),RooFit::Silence());

  return;
}