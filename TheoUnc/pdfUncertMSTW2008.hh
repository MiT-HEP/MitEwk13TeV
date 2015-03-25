#ifndef PDF_UNCERT_MSTW2008_HH
#define PDF_UNCERT_MSTW2008_HH

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

Double_t xsecUncert(TString rep_set_0, 
		    TString rep_set_1,
		    TString rep_set_2,
		    TString rep_set_3,
		    TString rep_set_4);

Double_t chRatioUncert(TString rep_set_p_0, 
		       TString rep_set_p_1,
		       TString rep_set_p_2,
		       TString rep_set_p_3,
		       TString rep_set_p_4,
		       TString rep_set_m_0, 
		       TString rep_set_m_1,
		       TString rep_set_m_2,
		       TString rep_set_m_3,
		       TString rep_set_m_4);

Double_t chSumUncert(TString rep_set_p_0, 
		     TString rep_set_p_1,
		     TString rep_set_p_2,
		     TString rep_set_p_3,
		     TString rep_set_p_4,
		     TString rep_set_m_0, 
		     TString rep_set_m_1,
		     TString rep_set_m_2,
		     TString rep_set_m_3,
		     TString rep_set_m_4,
		     Double_t xsec_p, // check that order is right!!                                                             
		     Double_t xsec_m);

Double_t wzRatioUncert(TString rep_set_p_0, 
		       TString rep_set_p_1,
		       TString rep_set_p_2,
		       TString rep_set_p_3,
		       TString rep_set_p_4,
		       TString rep_set_m_0, 
		       TString rep_set_m_1,
		       TString rep_set_m_2,
		       TString rep_set_m_3,
		       TString rep_set_m_4,
		       TString rep_set_z_0, 
		       TString rep_set_z_1,
		       TString rep_set_z_2,
		       TString rep_set_z_3,
		       TString rep_set_z_4,
		       Double_t xsec_p, // check that order is right!!                                                             
		       Double_t xsec_m);

Double_t getXsecUncert(vector<Double_t> &values, 
		       Double_t nomVal);

Double_t getChRatioUncert(vector<Double_t> &values_p, 
			  vector<Double_t> &values_m, 
			  Double_t nomVal_p,
			  Double_t nomVal_m);

Double_t getChSumUncert(vector<Double_t> &values_p, 
			vector<Double_t> &values_m, 
			Double_t nomVal_p,
			Double_t nomVal_m,
			Double_t xsec_p, // check that order is right!!                                                             
			Double_t xsec_m);

Double_t getWzRatioUncert(vector<Double_t> &values_p, 
			  vector<Double_t> &values_m, 
			  vector<Double_t> &values_z, 
			  Double_t nomVal_p,
			  Double_t nomVal_m,
			  Double_t nomVal_z,
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

//-------------------------------------------------------------------//

Double_t xsecUncert(TString rep_set_0, 
		    TString rep_set_1,
		    TString rep_set_2,
		    TString rep_set_3,
		    TString rep_set_4) {
  
  vector<Double_t>  vScale0;
  getVals(rep_set_0, vScale0);

  Double_t *vErr = new Double_t[5];
  vErr[0]=getXsecUncert(vScale0, vScale0[0]);

  vector<Double_t>  vScale1;
  getVals(rep_set_1, vScale1);
  vErr[1]=getXsecUncert(vScale1, vScale0[0]);

  vector<Double_t>  vScale2;
  getVals(rep_set_2, vScale2);
  vErr[2]=getXsecUncert(vScale2, vScale0[0]);

  vector<Double_t>  vScale3;
  getVals(rep_set_3, vScale3);
  vErr[3]=getXsecUncert(vScale3, vScale0[0]);

  vector<Double_t>  vScale4;
  getVals(rep_set_4, vScale4);
  vErr[4]=getXsecUncert(vScale4, vScale0[0]);

  return *max_element(vErr, vErr+5);
}

//-------------------------------------------------------------------//

Double_t chRatioUncert(TString rep_set_p_0, 
		       TString rep_set_p_1,
		       TString rep_set_p_2,
		       TString rep_set_p_3,
		       TString rep_set_p_4,
		       TString rep_set_m_0, 
		       TString rep_set_m_1,
		       TString rep_set_m_2,
		       TString rep_set_m_3,
		       TString rep_set_m_4) {

  Double_t *vErr = new Double_t[5];

  vector<Double_t>  vScaleP0; getVals(rep_set_p_0, vScaleP0);
  vector<Double_t>  vScaleM0; getVals(rep_set_m_0, vScaleM0);
  vErr[0]=getChRatioUncert(vScaleP0, vScaleM0, vScaleP0[0], vScaleM0[0]);

  vector<Double_t>  vScaleP1; getVals(rep_set_p_1, vScaleP1);
  vector<Double_t>  vScaleM1; getVals(rep_set_m_1, vScaleM1);
  vErr[1]=getChRatioUncert(vScaleP1, vScaleM1, vScaleP0[0], vScaleM0[0]);

  vector<Double_t>  vScaleP2; getVals(rep_set_p_2, vScaleP2);
  vector<Double_t>  vScaleM2; getVals(rep_set_m_2, vScaleM2);
  vErr[2]=getChRatioUncert(vScaleP2, vScaleM2, vScaleP0[0], vScaleM0[0]);

  vector<Double_t>  vScaleP3; getVals(rep_set_p_3, vScaleP3);
  vector<Double_t>  vScaleM3; getVals(rep_set_m_3, vScaleM3);
  vErr[3]=getChRatioUncert(vScaleP3, vScaleM3, vScaleP0[0], vScaleM0[0]);

  vector<Double_t>  vScaleP4; getVals(rep_set_p_4, vScaleP4);
  vector<Double_t>  vScaleM4; getVals(rep_set_m_4, vScaleM4);
  vErr[4]=getChRatioUncert(vScaleP4, vScaleM4, vScaleP0[0], vScaleM0[0]);

  return *max_element(vErr, vErr+5);

}

//-------------------------------------------------------------------//

Double_t chSumUncert(TString rep_set_p_0, 
		     TString rep_set_p_1,
		     TString rep_set_p_2,
		     TString rep_set_p_3,
		     TString rep_set_p_4,
		     TString rep_set_m_0, 
		     TString rep_set_m_1,
		     TString rep_set_m_2,
		     TString rep_set_m_3,
		     TString rep_set_m_4,
		     Double_t xsec_p, // check that order is right!!                                                             
                     Double_t xsec_m) {

  Double_t *vErr = new Double_t[5];

  vector<Double_t>  vScaleP0; getVals(rep_set_p_0, vScaleP0);
  vector<Double_t>  vScaleM0; getVals(rep_set_m_0, vScaleM0);
  vErr[0]=getChSumUncert(vScaleP0, vScaleM0, vScaleP0[0], vScaleM0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP1; getVals(rep_set_p_1, vScaleP1);
  vector<Double_t>  vScaleM1; getVals(rep_set_m_1, vScaleM1);
  vErr[1]=getChSumUncert(vScaleP1, vScaleM1, vScaleP0[0], vScaleM0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP2; getVals(rep_set_p_2, vScaleP2);
  vector<Double_t>  vScaleM2; getVals(rep_set_m_2, vScaleM2);
  vErr[2]=getChSumUncert(vScaleP2, vScaleM2, vScaleP0[0], vScaleM0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP3; getVals(rep_set_p_3, vScaleP3);
  vector<Double_t>  vScaleM3; getVals(rep_set_m_3, vScaleM3);
  vErr[3]=getChSumUncert(vScaleP3, vScaleM3, vScaleP0[0], vScaleM0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP4; getVals(rep_set_p_4, vScaleP4);
  vector<Double_t>  vScaleM4; getVals(rep_set_m_4, vScaleM4);
  vErr[4]=getChSumUncert(vScaleP4, vScaleM4, vScaleP0[0], vScaleM0[0], xsec_p, xsec_m);

  return *max_element(vErr, vErr+5);

}

//-------------------------------------------------------------------//

Double_t wzRatioUncert(TString rep_set_p_0, 
		       TString rep_set_p_1,
		       TString rep_set_p_2,
		       TString rep_set_p_3,
		       TString rep_set_p_4,
		       TString rep_set_m_0, 
		       TString rep_set_m_1,
		       TString rep_set_m_2,
		       TString rep_set_m_3,
		       TString rep_set_m_4,
		       TString rep_set_z_0, 
		       TString rep_set_z_1,
		       TString rep_set_z_2,
		       TString rep_set_z_3,
		       TString rep_set_z_4,
		       Double_t xsec_p, // check that order is right!!                                                             
		       Double_t xsec_m) {

  Double_t *vErr = new Double_t[5];

  vector<Double_t>  vScaleP0; getVals(rep_set_p_0, vScaleP0);
  vector<Double_t>  vScaleM0; getVals(rep_set_m_0, vScaleM0);
  vector<Double_t>  vScaleZ0; getVals(rep_set_z_0, vScaleZ0);
  vErr[0]=getWzRatioUncert(vScaleP0, vScaleM0, vScaleZ0, vScaleP0[0], vScaleM0[0], vScaleZ0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP1; getVals(rep_set_p_1, vScaleP1);
  vector<Double_t>  vScaleM1; getVals(rep_set_m_1, vScaleM1);
  vector<Double_t>  vScaleZ1; getVals(rep_set_z_1, vScaleZ1);
  vErr[1]=getWzRatioUncert(vScaleP1, vScaleM1, vScaleZ1, vScaleP0[0], vScaleM0[0], vScaleZ0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP2; getVals(rep_set_p_2, vScaleP2);
  vector<Double_t>  vScaleM2; getVals(rep_set_m_2, vScaleM2);
  vector<Double_t>  vScaleZ2; getVals(rep_set_z_2, vScaleZ2);
  vErr[2]=getWzRatioUncert(vScaleP2, vScaleM2, vScaleZ2, vScaleP0[0], vScaleM0[0], vScaleZ0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP3; getVals(rep_set_p_3, vScaleP3);
  vector<Double_t>  vScaleM3; getVals(rep_set_m_3, vScaleM3);
  vector<Double_t>  vScaleZ3; getVals(rep_set_z_3, vScaleZ3);
  vErr[3]=getWzRatioUncert(vScaleP3, vScaleM3, vScaleZ3, vScaleP0[0], vScaleM0[0], vScaleZ0[0], xsec_p, xsec_m);

  vector<Double_t>  vScaleP4; getVals(rep_set_p_4, vScaleP4);
  vector<Double_t>  vScaleM4; getVals(rep_set_m_4, vScaleM4);
  vector<Double_t>  vScaleZ4; getVals(rep_set_z_4, vScaleZ4);
  vErr[4]=getWzRatioUncert(vScaleP4, vScaleM4, vScaleZ4, vScaleP0[0], vScaleM0[0], vScaleZ0[0], xsec_p, xsec_m);

  return *max_element(vErr, vErr+5);

}

//-------------------------------------------------------------------//

Double_t getXsecUncert(vector<Double_t> &values, Double_t nomVal) {

  Double_t nomScale=values[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<values.size()+1; i+=2) {
    prevDiff=nomScale-values[i-1];
    thisDiff=nomScale-values[i];

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

  Double_t centralDiff=fabs(nomVal-nomScale);

  return max((sqrt(posUnc)+centralDiff)/nomScale*100,
	     (sqrt(negUnc)+centralDiff)/nomScale*100);

}

//-------------------------------------------------------------------//

Double_t getChRatioUncert(vector<Double_t> &values_p, 
			  vector<Double_t> &values_m, 
			  Double_t nomVal_p,
			  Double_t nomVal_m) {

  if ( values_p.size() != values_m.size() ) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomRatio=values_p[0]/values_m[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<values_p.size()+1; i+=2) {
    prevDiff=nomRatio-values_p[i-1]/values_m[i-1];
    thisDiff=nomRatio-values_p[i]/values_m[i];

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

  Double_t centralDiff=fabs(nomVal_p/nomVal_m-nomRatio);

  return max((sqrt(posUnc)+centralDiff)/nomRatio*100,
	     (sqrt(negUnc)+centralDiff)/nomRatio*100);

}

//-------------------------------------------------------------------//

Double_t getChSumUncert(vector<Double_t> &values_p, 
			vector<Double_t> &values_m, 
			Double_t nomVal_p,
			Double_t nomVal_m,
			Double_t xsec_p, // check that order is right!!                                                             
			Double_t xsec_m) {

  if ( values_p.size() != values_m.size() ) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomSum= (values_p[0]*xsec_p + values_m[0]*xsec_m)/(xsec_p + xsec_m);
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<values_p.size()+1; i+=2) {
    prevDiff=nomSum-(values_p[i-1]*xsec_p + values_m[i-1]*xsec_m)/(xsec_p + xsec_m);
    thisDiff=nomSum-(values_p[i]*xsec_p + values_m[i]*xsec_m)/(xsec_p + xsec_m);

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

  Double_t centralDiff=fabs( (nomVal_p*xsec_p + nomVal_m*xsec_m)/(xsec_p+xsec_m)-nomSum);

  return max((sqrt(posUnc)+centralDiff)/nomSum*100,
	     (sqrt(negUnc)+centralDiff)/nomSum*100);

}

Double_t getWzRatioUncert(vector<Double_t> &values_p, 
			  vector<Double_t> &values_m, 
			  vector<Double_t> &values_z, 
			  Double_t nomVal_p,
			  Double_t nomVal_m,
			  Double_t nomVal_z,
			  Double_t xsec_p, // check that order is right!!                                                             
			  Double_t xsec_m) {

  if ( values_p.size() != values_m.size() ||  values_p.size() != values_z.size() ) {
    cout <<  "Replica arrays are not the same size!" << endl;
    return 999;
  }

  Double_t nomRatio= (values_p[0]*xsec_p + values_m[0]*xsec_m)/(xsec_p + xsec_m)/values_z[0];
  Double_t prevDiff=0, thisDiff=0;
  Double_t posUnc=0, negUnc=0;

  for (UInt_t i=2; i<values_p.size()+1; i+=2) {
    prevDiff=nomRatio-(values_p[i-1]*xsec_p + values_m[i-1]*xsec_m)/(xsec_p + xsec_m)/values_z[i-1];
    thisDiff=nomRatio-(values_p[i]*xsec_p + values_m[i]*xsec_m)/(xsec_p + xsec_m)/values_z[i];
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

  Double_t centralDiff=fabs( (nomVal_p*xsec_p + nomVal_m*xsec_m)/(xsec_p+xsec_m)/nomVal_z-nomRatio);

  return max((sqrt(posUnc)+centralDiff)/nomRatio*100,
	     (sqrt(negUnc)+centralDiff)/nomRatio*100);

}

#endif
