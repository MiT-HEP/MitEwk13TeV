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

void ctq(double *iA0,double *iA1);
void mstw(double *iVal, int iN,double *iA);

void calcUncertCT10(TString ev_errors="test_f.txt", TString as_errors="test_as_f.txt") {
//void calcUncertCT10(TString ev_errors="test_d.txt", TString as_errors="test_as_d.txt") {
//void calcUncertCT10(TString ev_errors="test_ori.txt", TString as_errors="test_ori2.txt") {
//void calcUncertCT10(TString ev_errors="phillllll2.txt", TString as_errors="phillllll.txt") {

  ifstream ifs;
  ifs.open(ev_errors);
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

    //cout << prevDiff << " " << thisDiff << " ";

    if (prevDiff*thisDiff>0) {
      //if (fabs(prevDiff)>fabs(thisDiff)) thisDiff=prevDiff;
      if (fabs(prevDiff)>thisDiff) {
	//cout << "a";
	thisDiff=prevDiff;
      }
      if (thisDiff>0) {
	//cout << "b";
	posUnc+=thisDiff*thisDiff;
      }
      else {
	//cout << "c";
	negUnc+=thisDiff*thisDiff;
      }
    }
    else {
      if (thisDiff>0) {
	//cout << "d";
	posUnc+=thisDiff*thisDiff;
      }
      if (thisDiff<0) {
	//cout << "e";
	negUnc+=thisDiff*thisDiff;
      }
      
      if (prevDiff>0) {
	//cout << "f";
	posUnc+=prevDiff*prevDiff;
      }
      if (prevDiff<0) {
	//cout << "g";
	negUnc+=prevDiff*prevDiff;
      }
    }
    //cout << endl;
  }

  ifs.open(as_errors);
  assert(ifs.is_open());

  vector<Double_t>  vScaleAs;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t scale;
    ss >> scale;
    vScaleAs.push_back(scale);
  }
  ifs.close();

  Double_t nomScaleAs=vScaleAs[4];
  Double_t negUncAs=0, posUncAs=0;
  prevDiff=0;
  for (Int_t i=0; i<vScaleAs.size(); i++) {
    if (i==4) continue;

    thisDiff=nomScaleAs-vScaleAs[i];

    //cout << prevDiff << " " << thisDiff << " ";

    if (thisDiff>0) {
      //cout << "a";
      posUncAs+=thisDiff*thisDiff;
    }
    else {
      //cout << "b";
      negUncAs+=thisDiff*thisDiff;
    }
    //cout << endl;

  }
  //cout << "from err_ev: + " << sqrt(posUnc/nomScale/nomScale)*100/1.645 << " - " << sqrt(negUnc/nomScale/nomScale)*100/1.645 << endl;
  //cout << "from as_unc: + " << sqrt(posUncAs/nomScaleAs/nomScaleAs)*100/1.645 << " - " << sqrt(negUncAs/nomScaleAs/nomScaleAs)*100/1.645 << endl;
  cout << "CT10nlo:    + " << sqrt(posUnc/nomScale/nomScale+posUncAs/nomScaleAs/nomScaleAs)*100/1.645 << " - " << sqrt(negUnc/nomScale/nomScale+negUncAs/nomScaleAs/nomScaleAs)*100/1.645 << endl;
  //cout << "uncert:        " << max(sqrt(posUnc/nomScale/nomScale+posUncAs/nomScaleAs/nomScaleAs)*100/1.645,sqrt(negUnc/nomScale/nomScale+negUncAs/nomScaleAs/nomScaleAs)*100/1.645) << endl;

  //ctq(&vScale[0],&vScaleAs[0]);  

}

/*void ctq(double *iA0,double *iA1) {
  double *lE0 = new double[2]; mstw(lE0,53,iA0);
  double *lE1 = new double[2]; mstw(lE1,11,iA1); //lE1[0] = 0.; lE1[1] = 0.;
  cout << "from err_ev: + " << sqrt(lE0[0])*100/1.645 << " - " << sqrt(lE0[1])*100/1.645 << endl;
  cout << "from as_unc: + " << sqrt(lE1[0])*100/1.645 << " - " << sqrt(lE1[1])*100/1.645 << endl;
  cout << " CTEQ Error => + " << sqrt(lE0[0]+lE1[0])*100/1.645 << " - " << sqrt(lE0[1]+lE1[1])*100/1.645 << endl;
  }*/

void ctq(double *iA0,double *iA1) {
  double *lE0 = new double[2]; mstw(lE0,53,iA0);
  //cout << "from err_ev: + " << sqrt(lE0[0])*100/1.645 << " - " << sqrt(lE0[1])*100/1.645 << endl;
  double *lE1 = new double[2]; mstw(lE1,11,iA1); //lE1[0] = 0.; lE1[1] = 0.;                     
  //cout << "from as_unc: + " << sqrt(lE1[0])*100/1.645 << " - " << sqrt(lE1[1])*100/1.645 << endl;                                                        
  //cout << " CTEQ Error => +" << sqrt(lE0[0]+lE1[0])*100/1.645 << "-" << sqrt(lE0[1]+lE1[1])*100/1.645 << endl;
}

void mstw(double *iVal, int iN,double *iA) {
  double iC0=0;
  double lA = iA[0];
  if(iN == 11) lA = iA[4];
  //if(iN == 11) lA = iA[5];
  double lVal = lA;
  double lDiffP = 0;
  double lDiffM = 0;
  double lODiff = 0;
  /**/
  double start=1;
  if (iN==11) start=0;
  /**/
  for(int i0 = start; i0 < iN; i0++) {
    //for(int i0 = 1; i0 < iN; i0++) {
    if(iN == 11 && i0 > 8) continue;
    //if(iN == 11 && (i0 < 2  || i0 > 8)) continue;
    double pDiff = (lVal - iA[i0]);

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


/*void mstw(double *iVal, int iN,double *iA) {
  double iC0=0;
  double lA = iA[0];
  //if(iN == 11) lA = iA[5];
  if(iN == 11) lA = iA[4];
  double lVal = lA;
  double lDiffP = 0;
  double lDiffM = 0;
  double lODiff = 0;
  //for(int i0 = 1; i0 < iN; i0++) {
  for(int i0 = 0; i0 < iN; i0++) {
    //if(iN == 11 && (i0 < 2  || i0 > 8)) continue;
    if(iN == 11 && i0 > 8) continue;
    double pDiff = (lVal - iA[i0]);
    if(i0 % 2 == 1 && iN != 11) {
      lODiff = pDiff; 
      continue;
    }

    cout << lODiff << " " << pDiff << " ";

    if(pDiff*lODiff > 0 && iN != 11) {
      if(fabs(lODiff) > pDiff) {
	//cout << "a";
	pDiff = lODiff;
      }
      if(pDiff > 0) {
	//cout << "b";
	lDiffP += pDiff*pDiff;
      }
      if(pDiff < 0) {
	//cout << "c";
	lDiffM += pDiff*pDiff;
      }
      continue;
    }
    if(pDiff > 0) {
      cout << "d";
      lDiffP  += pDiff*pDiff;
    }
    if(pDiff < 0) {
      cout << "e";
      lDiffM  += pDiff*pDiff;
    }
    if(lODiff > 0) {
      cout << "f";
      lDiffP += lODiff*lODiff;
    }
    if(lODiff < 0) {
      cout << "g";
      lDiffM += lODiff*lODiff;
    }
    cout << endl;
  }
  double lCentral = fabs(iC0-lVal); if(iC0 == 0) lCentral = 0;
  iVal[0] = (sqrt(lDiffP)+lCentral)*(sqrt(lDiffP)+lCentral)/lVal/lVal;
  iVal[1] = (sqrt(lDiffM)+lCentral)*(sqrt(lDiffM)+lCentral)/lVal/lVal;
}
*/
