#ifndef PDF_UNCERT_CT10_HH
#define PDF_UNCERT_CT10_HH

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

Double_t xsecUncert(TString rep_sets, 
		    TString as_var);

Double_t chRatioUncert(TString rep_sets_p, 
		       TString as_var_p,
		       TString rep_sets_m, 
		       TString as_var_m);

Double_t chSumUncert(TString rep_sets_p, 
		     TString as_var_p,
		     TString rep_sets_m, 
		     TString as_var_m,
		     Double_t xsec_p, // check that order is right!!
		     Double_t xsec_m);

Double_t wzRatioUncert(TString rep_sets_p, 
		       TString as_var_p,
		       TString rep_sets_m, 
		       TString as_var_m,
		       TString rep_sets_z, 
		       TString as_var_z,
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

Double_t xsecUncert(TString rep_sets, TString as_var) {

  vector<Double_t>  vScale;
  getVals(rep_sets, vScale);
  
  Double_t nomScale=vScale[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScale.size()+1; i+=2) {
    prevDiff=nomScale-vScale[i-1];
    thisDiff=nomScale-vScale[i];

    if (prevDiff*thisDiff>0) {
      if (fabs(prevDiff)>thisDiff) thisDiff=prevDiff;
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else negUnc+=thisDiff*thisDiff;
    }
    else {
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else if (thisDiff<0) negUnc+=thisDiff*thisDiff;
      
      if (prevDiff>0) posUnc+=prevDiff*prevDiff;
      if (prevDiff<0) negUnc+=prevDiff*prevDiff;
    }
  }

  vector<Double_t>  vScaleAs;
  getVals(as_var, vScaleAs);

  Double_t nomScaleAs=vScaleAs[4]; // should not have this hardcoded
  Double_t negUncAs=0, posUncAs=0;
  prevDiff=0;
  for (UInt_t i=0; i<vScaleAs.size(); i++) {
    if (i==4) continue; 
    thisDiff=nomScaleAs-vScaleAs[i];
    if (thisDiff>0) posUncAs+=thisDiff*thisDiff;
    else negUncAs+=thisDiff*thisDiff;
   }

  return max(sqrt(posUnc/nomScale/nomScale+posUncAs/nomScaleAs/nomScaleAs)*100/1.645, 
	     sqrt(negUnc/nomScale/nomScale+negUncAs/nomScaleAs/nomScaleAs)*100/1.645);

}

Double_t chRatioUncert(TString rep_sets_p, 
		       TString as_var_p,
		       TString rep_sets_m, 
		       TString as_var_m) {

  vector<Double_t>  vScaleP;
  getVals(rep_sets_p, vScaleP);

  vector<Double_t>  vScaleM;
  getVals(rep_sets_m, vScaleM);

  if ( vScaleM.size() != vScaleP.size() ) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomRatio=vScaleP[0]/vScaleM[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScaleP.size()+1; i+=2) {
    prevDiff=nomRatio-vScaleP[i-1]/vScaleM[i-1];
    thisDiff=nomRatio-vScaleP[i]/vScaleM[i];

    if (prevDiff*thisDiff>0) {
      if (fabs(prevDiff)>thisDiff) thisDiff=prevDiff;
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else negUnc+=thisDiff*thisDiff;
    }
    else {
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else if (thisDiff<0) negUnc+=thisDiff*thisDiff;
      
      if (prevDiff>0) posUnc+=prevDiff*prevDiff;
      if (prevDiff<0) negUnc+=prevDiff*prevDiff;
    }
  }

  vector<Double_t>  vScaleAsP;
  getVals(as_var_p, vScaleAsP);

  vector<Double_t>  vScaleAsM;
  getVals(as_var_m, vScaleAsM);

  if ( vScaleAsM.size() != vScaleAsP.size() ) {
    cout << "A_s arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomRatioAs=vScaleAsP[4]/vScaleAsM[4]; //should not have this hard coded
  Double_t negUncAs=0, posUncAs=0;
  prevDiff=0;
  for (UInt_t i=0; i<vScaleAsP.size(); i++) {
    if (i==4) continue;
    thisDiff=nomRatioAs-vScaleAsP[i]/vScaleAsM[i];
    if (thisDiff>0) posUncAs+=thisDiff*thisDiff;
    else negUncAs+=thisDiff*thisDiff;
  }
  
  return max(sqrt(posUnc/nomRatio/nomRatio+posUncAs/nomRatioAs/nomRatioAs)*100/1.645, 
	     sqrt(negUnc/nomRatio/nomRatio+negUncAs/nomRatioAs/nomRatioAs)*100/1.645);
  
}

Double_t chSumUncert(TString rep_sets_p, 
		     TString as_var_p,
		     TString rep_sets_m, 
		     TString as_var_m,
		     Double_t xsec_p, 
		     Double_t xsec_m) {

  vector<Double_t>  vScaleP;
  getVals(rep_sets_p, vScaleP);

  vector<Double_t>  vScaleM;
  getVals(rep_sets_m, vScaleM);

  if ( vScaleM.size() != vScaleP.size() ) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomSum = (vScaleP[0]*xsec_p + vScaleM[0]*xsec_m)/(xsec_p + xsec_m);
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScaleP.size()+1; i+=2) {
    prevDiff=nomSum-(vScaleP[i-1]*xsec_p + vScaleM[i-1]*xsec_m)/(xsec_p + xsec_m);
    thisDiff=nomSum-(vScaleP[i]*xsec_p + vScaleM[i]*xsec_m)/(xsec_p + xsec_m);

    if (prevDiff*thisDiff>0) {
      if (fabs(prevDiff)>thisDiff) thisDiff=prevDiff;
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else negUnc+=thisDiff*thisDiff;
    }
    else {
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else if (thisDiff<0) negUnc+=thisDiff*thisDiff;
      
      if (prevDiff>0) posUnc+=prevDiff*prevDiff;
      if (prevDiff<0) negUnc+=prevDiff*prevDiff;
    }
  }

  vector<Double_t>  vScaleAsP;
  getVals(as_var_p, vScaleAsP);

  vector<Double_t>  vScaleAsM;
  getVals(as_var_m, vScaleAsM);

  if ( vScaleAsM.size() != vScaleAsP.size() ) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomSumAs=(vScaleAsP[4]*xsec_p + vScaleAsM[4]*xsec_m)/(xsec_p + xsec_m); // shouldn't have this hard-coded
  Double_t negUncAs=0, posUncAs=0;
  prevDiff=0;
  for (UInt_t i=0; i<vScaleAsP.size(); i++) {
    if (i==4) continue;
    thisDiff=nomSumAs-(vScaleAsP[i]*xsec_p + vScaleAsM[i]*xsec_m)/(xsec_p + xsec_m);
    if (thisDiff>0) posUncAs+=thisDiff*thisDiff;
    else negUncAs+=thisDiff*thisDiff;
   }

  return max(sqrt(posUnc/nomSum/nomSum+posUncAs/nomSumAs/nomSumAs)*100/1.645, 
	     sqrt(negUnc/nomSum/nomSum+negUncAs/nomSumAs/nomSumAs)*100/1.645);

}

Double_t wzRatioUncert(TString rep_sets_p, 
		       TString as_var_p,
		       TString rep_sets_m, 
		       TString as_var_m,
		       TString rep_sets_z, 
		       TString as_var_z,
		       Double_t xsec_p, // check that order is right!!
		       Double_t xsec_m) {


  vector<Double_t>  vScaleP;
  getVals(rep_sets_p, vScaleP);

  vector<Double_t>  vScaleM;
  getVals(rep_sets_m, vScaleM);

  vector<Double_t>  vScaleZ;
  getVals(rep_sets_z, vScaleZ);

  if ( vScaleM.size() != vScaleP.size() || vScaleM.size() != vScaleZ.size()) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomRatio = (vScaleP[0]*xsec_p + vScaleM[0]*xsec_m)/(xsec_p + xsec_m)/vScaleZ[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<vScaleP.size()+1; i+=2) {
    prevDiff=nomRatio-(vScaleP[i-1]*xsec_p + vScaleM[i-1]*xsec_m)/(xsec_p + xsec_m)/vScaleZ[i-1];
    thisDiff=nomRatio-(vScaleP[i]*xsec_p + vScaleM[i]*xsec_m)/(xsec_p + xsec_m)/vScaleZ[i];

    if (prevDiff*thisDiff>0) {
      if (fabs(prevDiff)>thisDiff) thisDiff=prevDiff;
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else negUnc+=thisDiff*thisDiff;
    }
    else {
      if (thisDiff>0) posUnc+=thisDiff*thisDiff;
      else if (thisDiff<0) negUnc+=thisDiff*thisDiff;
      
      if (prevDiff>0) posUnc+=prevDiff*prevDiff;
      if (prevDiff<0) negUnc+=prevDiff*prevDiff;
    }
  }

  vector<Double_t>  vScaleAsP;
  getVals(as_var_p, vScaleAsP);

  vector<Double_t>  vScaleAsM;
  getVals(as_var_m, vScaleAsM);

  vector<Double_t>  vScaleAsZ;
  getVals(as_var_z, vScaleAsZ);

  if ( vScaleAsM.size() != vScaleAsP.size() || vScaleAsM.size() != vScaleAsZ.size()) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomRatioAs=(vScaleAsP[4]*xsec_p + vScaleAsM[4]*xsec_m)/(xsec_p + xsec_m)/vScaleAsZ[4]; // shouldn't have this hard-coded
  Double_t negUncAs=0, posUncAs=0;

  prevDiff=0;
  for (UInt_t i=0; i<vScaleAsP.size(); i++) {
    if (i==4) continue;
    thisDiff=nomRatioAs-(vScaleAsP[i]*xsec_p + vScaleAsM[i]*xsec_m)/(xsec_p + xsec_m)/vScaleAsZ[i];
    if (thisDiff>0) posUncAs+=thisDiff*thisDiff;
    else negUncAs+=thisDiff*thisDiff;
  }
  
  return max(sqrt(posUnc/nomRatio/nomRatio+posUncAs/nomRatioAs/nomRatioAs)*100/1.645, 
	     sqrt(negUnc/nomRatio/nomRatio+negUncAs/nomRatioAs/nomRatioAs)*100/1.645);

}

#endif
