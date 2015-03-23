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

//void nnpdf(double *iA);

void calcUncertNNPDF23(TString input="test2.txt") {

  ifstream ifs;
  ifs.open(input);
  assert(ifs.is_open());
  string line;

  vector<Double_t>  vAcc;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t acc;
    ss >> acc;
    vAcc.push_back(acc);
  }

  Double_t stdP=0, stdM=0;
  Int_t nP=0, nM=0;
  for (Int_t i=1; i<vAcc.size(); i++) {
    Double_t diff=vAcc[i]-vAcc[0];
    if (diff>0) { 
      stdP+=diff*diff;
      nP++;
    }
    else if (diff<0) {
      stdM+=diff*diff;
      nM++;
    }
  }
  Double_t errP=sqrt(stdP/(nP-1.0));
  Double_t errM=sqrt(stdM/(nM-1.0));
  cout << "NNPDF err: + " << errP/vAcc[0]*100 << " - " << errM/vAcc[0]*100 << endl;

  //nnpdf(&vAcc[0]);
  
}
/*
void nnpdf(double *iA) {
  double lA = iA[0];
  double lVal = lA;
  //double lVal = 0;                                                                                                                                     
  //for(int i0 = 0; i0 < 101; i0++) {                                                                                                                    
  //  double lA = iA[i0]; double lB=iB[i0];                                                                                                              
  //  lVal += lA/lB;                                                                                                                                     
  // }                                                                                                                                                   
  //lVal/=101;                                                                                                                                           
  int lN = 305;
  double lDiffP = 0; int lNP = 0;
  double lDiffM = 0; int lNM = 0;
  for(int i0 = 1; i0 < lN; i0++) {
    double pDiff = (iA[i0] - lVal);
    double pVal = (iA[i0]);
    //cout << "===> " << iA[i0] << " -- " << iB[i0] << " -- " << i0 << endl;                                                                             
    if(pDiff > 0) {lDiffP += pDiff*pDiff; lNP++;}
    if(pDiff < 0) {lDiffM += pDiff*pDiff; lNM++;}
  }
  double lMVal = lDiffM/lNM;
  double lEMinus = sqrt(lNM/(lNM-1.)*lDiffM/lNM);//sqrt(lNM/(lNM-1.)*(lDiffM/lNM - lVal*lVal));                                                          
  double lEPlus  = sqrt(lNP/(lNP-1.)*lDiffP/lNP);
  cout << "NNPDF Error => +" <<  (lEPlus/lVal)*100  << " -" << (lEMinus/lVal)*100 << endl;
}
*/
