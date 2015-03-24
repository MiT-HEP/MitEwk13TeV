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
#include <TLine.h>
#include "Math/LorentzVector.h"     // 4-vector class

using namespace std;

#endif

void nnpdf(double *iA);

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

  //TH1D* hTest = new TH1D("hTest", "", 30, 0.32, 0.35);
  TH1D* hTest = new TH1D("hTest", "", 30, 0.46, 0.48);

  Double_t avg=0, stdv=0;
  Double_t stdP=0, stdM=0;
  Double_t stdP2=0, stdM2=0;

  Int_t nP=0, nM=0;
  Int_t nP2=0, nM2=0;
  for (Int_t i=0; i<vAcc.size(); i++) {
    hTest->Fill(vAcc[i]);
    avg+=vAcc[i];
  }
  avg/=double(vAcc.size());

  for (Int_t i=0; i<vAcc.size(); i++) {
    stdv+=(vAcc[i]-avg)*(vAcc[i]-avg);
    if ((vAcc[i]-avg)>0) {
      stdP+=(vAcc[i]-avg)*(vAcc[i]-avg);
      nP++;
    }
    else if ((vAcc[i]-avg)<0) {
      stdM+=(vAcc[i]-avg)*(vAcc[i]-avg);
      nM++;
    }

    if (i>0 && (vAcc[i]-vAcc[0])>0) {
      stdP2+=(vAcc[i]-vAcc[0])*(vAcc[i]-vAcc[0]);
      nP2++;
    }
    else if (i>0 && (vAcc[i]-vAcc[0])<0) {
      stdM2+=(vAcc[i]-vAcc[0])*(vAcc[i]-vAcc[0]);
      nM2++;
    }
  }
  stdv=sqrt(stdv/sqrt(vAcc.size()-1));
  stdP=sqrt(stdP/(nP-1));
  stdM=sqrt(stdM/(nM-1));
  stdP2=sqrt(stdP2/(nP2-1));
  stdM2=sqrt(stdM2/(nM2-1));
  cout << avg << " + " << stdP << " - " << stdM << endl;
  cout << stdP/avg*100 << ", " << stdM/avg*100 << " " << nP << " " << nM << endl;
  cout << vAcc[0] << " + " << stdP2 << " - " << stdM2 << endl;
  cout << stdP2/vAcc[0]*100 << ", " << stdM2/vAcc[0]*100 << " " << nP2 << " " << nM2 << endl;

  nnpdf(&vAcc[0]);

    hTest->Draw();
  TLine *l = new TLine(vAcc[0],0,vAcc[0],50);
  l->SetLineColor(kRed);
  l->Draw();
  
}

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
  cout << lNP << " " << lNM << endl;
}
