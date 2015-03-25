#ifndef PDF_UNCERT_NNPDF23_HH
#define PDF_UNCERT_NNPDF23_HH

#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O

using namespace std;

void getVals(TString filename, vector<Double_t> &values);

Double_t xsecUncert(TString replicas);

Double_t chRatioUncert(TString replicas_p, 
		       TString replicas_m);

Double_t chSumUncert(TString replicas_p, 
		     TString replicas_m, 
		     Double_t xsec_p, // check that order is right!!
		     Double_t xsec_m);

Double_t wzRatioUncert(TString replicas_p, 
		       TString replicas_m, 
		       TString replicas_z, 
		       Double_t xsec_p, // check that order is right!!
		       Double_t xsec_m);

//-------------------------------------------------------------------//
//-------------------------------------------------------------------//

void getVals(TString filename, vector<Double_t> &values) {

  ifstream ifs;
  ifs.open(filename);
  assert(ifs.is_open());
  string line;

  while (getline(ifs,line)) {
    stringstream ss(line);
    Double_t scale;
    ss >> scale;
    values.push_back(scale);
  }
  ifs.close();
}

Double_t xsecUncert(TString replicas) {

  vector<Double_t>  vScale;
  getVals(replicas, vScale);
  
  Double_t nom=vScale[0];
  Double_t stdU=0, stdL=0;
  Int_t nP=0, nM=0;

  for (UInt_t i=1; i<vScale.size(); i++) {
    if ((vScale[i]-nom)>0) {
      stdU+=(vScale[i]-nom)*(vScale[i]-nom);
      nP++;
    }
    else if ((vScale[i]-nom)<0) {
      stdL+=(vScale[i]-nom)*(vScale[i]-nom);
      nM++;
    }
  }

  stdU=sqrt(stdU/(nP-1));
  stdL=sqrt(stdL/(nM-1));
  
  return max(stdU/nom*100, stdL/nom*100);

}

Double_t chRatioUncert(TString replicas_p, 
		       TString replicas_m) {

  vector<Double_t>  vScaleP; getVals(replicas_p, vScaleP);
  vector<Double_t>  vScaleM; getVals(replicas_m, vScaleM);
  
  Double_t nom=vScaleP[0]/vScaleM[0];
  Double_t stdU=0, stdL=0;
  Int_t nP=0, nM=0;

  if ( vScaleP.size() != vScaleM.size() ) return 999;

  for (UInt_t i=1; i<vScaleP.size(); i++) {
    if ((vScaleP[i]/vScaleM[i]-nom)>0) {
      stdU+=(vScaleP[i]/vScaleM[i]-nom)*(vScaleP[i]/vScaleM[i]-nom);
      nP++;
    }
    else if ((vScaleP[i]/vScaleM[i]-nom)<0) {
      stdL+=(vScaleP[i]/vScaleM[i]-nom)*(vScaleP[i]/vScaleM[i]-nom);
      nM++;
    }
  }

  stdU=sqrt(stdU/(nP-1));
  stdL=sqrt(stdL/(nM-1));
  
  return max(stdU/nom*100, stdL/nom*100);

}

Double_t chSumUncert(TString replicas_p, 
		     TString replicas_m,
		     Double_t xsec_p, // check that order is right!!
		     Double_t xsec_m) {

  vector<Double_t>  vScaleP; getVals(replicas_p, vScaleP);
  vector<Double_t>  vScaleM; getVals(replicas_m, vScaleM);
  
  Double_t nom=(vScaleP[0]*xsec_p+vScaleM[0]*xsec_m)/(xsec_p+xsec_m);
  Double_t stdU=0, stdL=0;
  Int_t nP=0, nM=0;

  if ( vScaleP.size() != vScaleM.size() ) return 999;

  for (UInt_t i=1; i<vScaleP.size(); i++) {
    Float_t sum=(vScaleP[i]*xsec_p+vScaleM[i]*xsec_m)/(xsec_p+xsec_m);
    if ((sum-nom)>0) {
      stdU+=(sum-nom)*(sum-nom);
      nP++;
    }
    else if ((sum-nom)<0) {
      stdL+=(sum-nom)*(sum-nom);
      nM++;
    }
  }

  stdU=sqrt(stdU/(nP-1));
  stdL=sqrt(stdL/(nM-1));
  
  return max(stdU/nom*100, stdL/nom*100);

}

Double_t wzRatioUncert(TString replicas_p, 
		       TString replicas_m,
		       TString replicas_z,
		       Double_t xsec_p, // check that order is right!!
		       Double_t xsec_m) {

  vector<Double_t>  vScaleP; getVals(replicas_p, vScaleP);
  vector<Double_t>  vScaleM; getVals(replicas_m, vScaleM);
  vector<Double_t>  vScaleZ; getVals(replicas_z, vScaleZ);
  
  Double_t nom=(vScaleP[0]*xsec_p+vScaleM[0]*xsec_m)/(xsec_p+xsec_m)/vScaleZ[0];
  Double_t stdU=0, stdL=0;
  Int_t nP=0, nM=0;

  if ( vScaleP.size() != vScaleM.size() || vScaleP.size() != vScaleZ.size() ) return 999;

  for (UInt_t i=1; i<vScaleP.size(); i++) {
    Float_t sum=(vScaleP[i]*xsec_p+vScaleM[i]*xsec_m)/(xsec_p+xsec_m)/vScaleZ[i];
    if ((sum-nom)>0) {
      stdU+=(sum-nom)*(sum-nom);
      nP++;
    }
    else if ((sum-nom)<0) {
      stdL+=(sum-nom)*(sum-nom);
      nM++;
    }
  }

  stdU=sqrt(stdU/(nP-1));
  stdL=sqrt(stdL/(nM-1));
  
  return max(stdU/nom*100, stdL/nom*100);

}

#endif
