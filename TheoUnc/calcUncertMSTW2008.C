#if !defined(__CINT__) //|| defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include "Math/LorentzVector.h"     // 4-vector class

using namespace std;

#endif

//void mstw(double *iA0,double *iA1,double *iA2,double *iA3,double *iA4);
//void mstw(double *iVal, int iN,double *iA,double iC0 = 0);

void getErr(TString filename, Double_t *errH, Double_t *errL);

void calcUncertMSTW2008(TString file0="test0.txt", TString file1="test1.txt", TString file2="test2.txt", TString file3="test3.txt", TString file4="test4.txt") {

  Double_t *vErrH = new Double_t[5];
  Double_t *vErrL = new Double_t[5];

  getErr(file0, &vErrH[0], &vErrL[0]);
  getErr(file1, &vErrH[1], &vErrL[1]);
  getErr(file2, &vErrH[2], &vErrL[2]);
  getErr(file3, &vErrH[3], &vErrL[3]);
  getErr(file4, &vErrH[4], &vErrL[4]);

  //cout << vErrH[0] << " " << vErrH[1] << " " << vErrH[2] << " " << vErrH[3] << " " << vErrH[4] << endl;
  //cout << vErrL[0] << " " << vErrL[1] << " " << vErrL[2] << " " << vErrL[3] << " " << vErrL[4] << endl;

  cout << "MSTW: + " << *max_element(vErrH, vErrH+5) << " - " << *max_element(vErrL, vErrL+5) << endl;

}

void getErr(TString filename, Double_t *errH, Double_t *errL) {

  ifstream ifs;
  ifs.open(filename);
  assert(ifs.is_open());
  string line;
  
  vector<Double_t>  vScale;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t scale;
    ss >> scale;
    vScale.push_back(scale);
  }
  ifs.close();

  Double_t nomScale=vScale[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (Int_t i=2; i<vScale.size()+1; i+=2) {
    prevDiff=nomScale-vScale[i-1];
    thisDiff=nomScale-vScale[i];

    if (prevDiff*thisDiff>0) {
      if (fabs(prevDiff)>thisDiff) {
	thisDiff=prevDiff;
      }
      if (thisDiff>0) {
	posUnc+=thisDiff*thisDiff;
      }
      else {
	negUnc+=thisDiff*thisDiff;
      }
    }
    else {
      if (thisDiff>0) {
	posUnc+=thisDiff*thisDiff;
      }
      if (thisDiff<0) {
	negUnc+=thisDiff*thisDiff;
      }
      
      if (prevDiff>0) {
	posUnc+=prevDiff*prevDiff;
      }
      if (prevDiff<0) {
	negUnc+=prevDiff*prevDiff;
      }
    }
  }

  errH[0]=sqrt(posUnc/nomScale/nomScale)*100;
  errL[0]=sqrt(negUnc/nomScale/nomScale)*100;

}
/*
void mstw(double *iA0,double *iA1,double *iA2,double *iA3,double *iA4){
  double *lE0 = new double[2]; mstw(lE0,41,iA0,iA0[0]);
  double *lE1 = new double[2]; mstw(lE1,41,iA1,iA0[0]);
  double *lE2 = new double[2]; mstw(lE2,41,iA2,iA0[0]);
  double *lE3 = new double[2]; mstw(lE3,41,iA3,iA0[0]);
  double *lE4 = new double[2]; mstw(lE4,41,iA4,iA0[0]);
  double lMax = lE0[0]; double lMin = lE0[1];
  //cout << "MSTW 0 +" << sqrt(lE0[0]) << "-" << sqrt(lE0[1]) << endl;                                                                                   
  // cout << "MSTW 1 +" << sqrt(lE1[0]) << "-" << sqrt(lE1[1]) << endl;                                                                                  
  //cout << "MSTW 2 +" << sqrt(lE2[0]) << "-" << sqrt(lE2[1]) << endl;                                                                                   
  //cout << "MSTW 3 +" << sqrt(lE3[0]) << "-" << sqrt(lE3[1]) << endl;                                                                                   
  //cout << "MSTW 4 +" << sqrt(lE4[0]) << "-" << sqrt(lE4[1]) << endl;                                                                                   
  if(lE1[0] > lMax) lMax = lE1[0];
  if(lE2[0] > lMax) lMax = lE2[0];
  if(lE3[0] > lMax) lMax = lE3[0];
  if(lE4[0] > lMax) lMax = lE4[0];
  if(lE1[1] > lMin) lMin = lE1[1];
  if(lE2[1] > lMin) lMin = lE2[1];
  if(lE3[1] > lMin) lMin = lE3[1];
  if(lE4[1] > lMin) lMin = lE4[1];
  cout << " MSTW Error => +" << sqrt(lMax)*100 << "-" << sqrt(lMin)*100 << endl;
}

void mstw(double *iVal, int iN,double *iA,double iC0 = 0) {
  double lA = iA[0];
  if(iN == 11) lA = iA[5];
  double lVal = lA;
  double lDiffP = 0;
  double lDiffM = 0;
  double lODiff = 0;
  for(int i0 = 1; i0 < iN; i0++) {
    if(iN == 11 && (i0 < 2  || i0 > 8)) continue;
    double pDiff = (lVal - iA[i0]);
    //cout << " ===> " << i0 << " -- " << pDiff << " -- " << iA[i0] << " -- " << lA << " -- " << sqrt(lDiffP)/lVal << " == " << sqrt(lDiffM)/lVal << end\
    l;                                                                                                                                                       
    if(i0 % 2 == 1 && iN != 11) {lODiff = pDiff; continue;}
    if(pDiff*lODiff > 0 && iN != 6) {
      if(fabs(lODiff) > pDiff) pDiff = lODiff;
      if(pDiff > 0) lDiffP += pDiff*pDiff;
      if(pDiff < 0) lDiffM += pDiff*pDiff;
      continue;
    }
    if(pDiff > 0)lDiffP  += pDiff*pDiff;
    if(pDiff < 0)lDiffM  += pDiff*pDiff;
    if(lODiff > 0)lDiffP += lODiff*lODiff;
    if(lODiff < 0)lDiffM += lODiff*lODiff;
  }
  //cout << " Error => +" << sqrt(lDiffP)*100/lVal << "-" << sqrt(lDiffM)*100/lVal << endl;                                                              
  //iVal[0] = lDiffP/lVal; iVal[1] = lDiffM/lVal;                                                                                                        
  double lCentral = fabs(iC0-lVal); if(iC0 == 0) lCentral = 0;
  iVal[0] = (sqrt(lDiffP)+lCentral)*(sqrt(lDiffP)+lCentral)/lVal/lVal;
  iVal[1] = (sqrt(lDiffM)+lCentral)*(sqrt(lDiffM)+lCentral)/lVal/lVal;
}
*/
