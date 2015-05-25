//================================================================================================
//
// Compute cross sections and produce summary plots
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#endif

void signalyields_w(const TString infilename="input_w.txt", const TString outputDir=".")
{   

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  double nWe_exp,          nWep_exp,          nWem_exp,          nWm_exp,          nWmp_exp,          nWmm_exp;
  double nWe,              nWep,              nWem,              nWm,              nWmp,              nWmm;
  double nWe_err,          nWep_err,          nWem_err,          nWm_err,          nWmp_err,          nWmm_err;
  double rWe,              rWep,              rWem,              rWm,              rWmp,              rWmm;
  double nWe_comb,         nWep_comb,         nWem_comb,         nWm_comb,         nWmp_comb,         nWmm_comb;
  double nWe_diff,         nWep_diff,         nWem_diff,         nWm_diff,         nWmp_diff,         nWmm_diff;
  double rWe_syserr,       rWep_syserr,       rWem_syserr,       rWm_syserr,       rWmp_syserr,       rWmm_syserr;
  double nWe_comb_syserr,  nWep_comb_syserr,  nWem_comb_syserr,  nWm_comb_syserr,  nWmp_comb_syserr,  nWmm_comb_syserr;
  double rWe_staterr,      rWep_staterr,      rWem_staterr,      rWm_staterr,      rWmp_staterr,      rWmm_staterr;
  double nWe_comb_staterr, nWep_comb_staterr, nWem_comb_staterr, nWm_comb_staterr, nWmp_comb_staterr, nWmm_comb_staterr;

  ifstream ifs;
  ifs.open(infilename.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;

    if(state==0) {
      // Expected yields
      stringstream ss(line);
      string label;
      ss >> label >> nWe_exp >> nWep_exp >> nWem_exp >> nWm_exp >> nWmp_exp >> nWmm_exp;
      state++;
    } else if(state==1) {
      // Yields from MitEWK code
      stringstream ss(line);
      string label;
      ss >> label >> nWe >> nWep >> nWem >> nWm >> nWmp >> nWmm;
      state++;
    } else if(state==2) {
      // Uncertainties on yields from MitEWK code
      stringstream ss(line);
      string label;
      ss >> label >> nWe_err >> nWep_err >> nWem_err >> nWm_err >> nWmp_err >> nWmm_err;
      state++;
    } else if(state==3) {
      // r value from Combine
      stringstream ss(line);
      string label;
      ss >> label >> rWe >> rWep >> rWem >> rWm >> rWmp >> rWmm;
      state++;
      // Yields from Combine
      nWe_comb  = rWe*nWe_exp;
      nWep_comb = rWep*nWep_exp;
      nWem_comb = rWem*nWem_exp; 
      nWm_comb  = rWm*nWm_exp;
      nWmp_comb = rWmp*nWmp_exp;
      nWmm_comb = rWmm*nWmm_exp; 
      // Differences in yields from MitEWK code and Combine
      nWe_diff  = TMath::Abs(nWe - nWe_comb);
      nWep_diff = TMath::Abs(nWep - nWep_comb);
      nWem_diff = TMath::Abs(nWem - nWem_comb);
      nWm_diff  = TMath::Abs(nWm - nWm_comb);
      nWmp_diff = TMath::Abs(nWmp - nWmp_comb);
      nWmm_diff = TMath::Abs(nWmm - nWmm_comb); 
   } else if(state==4){
     // Systematic uncertainties on r value from Combine
     stringstream ss(line);
     string label;
     ss >> label >> rWe_syserr >> rWep_syserr >> rWem_syserr >> rWm_syserr >> rWmp_syserr >> rWmm_syserr;
     state++;
     // Uncertainties on yields from Combine
     nWe_comb_syserr  = rWe_syserr*nWe_exp;
     nWep_comb_syserr = rWep_syserr*nWep_exp;
     nWem_comb_syserr = rWem_syserr*nWem_exp;
     nWm_comb_syserr  = rWm_syserr*nWm_exp;
     nWmp_comb_syserr = rWmp_syserr*nWmp_exp;
     nWmm_comb_syserr = rWmm_syserr*nWmm_exp;
   } else{
     // Statistical uncertainties on r value from Combine
     stringstream ss(line);
     string label;
     ss >> label >> rWe_staterr >> rWep_staterr >> rWem_staterr >> rWm_staterr >> rWmp_staterr >> rWmm_staterr;
     state++;
     // Uncertainties on yields from Combine
     nWe_comb_staterr  = rWe_staterr*nWe_exp;
     nWep_comb_staterr = rWep_staterr*nWep_exp;
     nWem_comb_staterr = rWem_staterr*nWem_exp;
     nWm_comb_staterr  = rWm_staterr*nWm_exp;
     nWmp_comb_staterr = rWmp_staterr*nWmp_exp;
     nWmm_comb_staterr = rWmm_staterr*nWmm_exp;
   }
  }
  ifs.close();

  cout << endl;
  cout << setprecision(0) << fixed; 

  cout << "Channel               Expected Yield        Yield (MitEWK Code)   Yield (Combine)                      Difference" << endl << endl;
  cout << "Wenu                  " << nWe_exp << "                 " << nWe << " +/- " << nWe_err << "         " << nWe_comb << " +/- " << nWe_comb_staterr << "(stat) +/- " << nWe_comb_syserr << "(sys)     " << nWe_diff << endl;
  cout << "Wenu_p                " << nWep_exp << "                 " << nWep << " +/- " << nWep_err << "         " << nWep_comb << " +/- " << nWep_comb_staterr << "(stat) +/- " << nWep_comb_syserr << "(sys)     " << nWep_diff << endl;
  cout << "Wenu_m                " << nWem_exp << "                 " << nWem << " +/- " << nWem_err << "         " << nWem_comb << " +/- " << nWem_comb_staterr << "(stat) +/- " << nWem_comb_syserr << "(sys)     " << nWem_diff << endl;
  cout << "Wmunu                 " << nWm_exp << "                 " << nWm << " +/- " << nWm_err << "         " << nWm_comb << " +/- " << nWm_comb_staterr << "(stat) +/- " << nWm_comb_syserr << "(sys)     " << nWm_diff << endl;
  cout << "Wmunu_p               " << nWmp_exp << "                 " << nWmp << " +/- " << nWmp_err << "         " << nWmp_comb << " +/- " << nWmp_comb_staterr << "(stat) +/- " << nWmp_comb_syserr << "(sys)     " << nWmp_diff << endl;
  cout << "Wmunu_m               " << nWmm_exp << "                 " << nWmm << " +/- " << nWmm_err << "         " << nWmm_comb << " +/- " << nWmm_comb_staterr << "(stat) +/- " << nWmm_comb_syserr << "(sys)     " << nWmm_diff << endl;
  cout << endl;
}
