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

void signalyields_z(const TString infilename="input_z.txt", const TString outputDir=".")
{   

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  double rZee, Nz_init_Zee, effHLT_Zee, effSel_Zee, effGsf,         nZee;
  double rZmm, Nz_init_Zmm, effHLT_Zmm, effSel_Zmm, effSta, effTrk, nZmm;

  double Nz_Zee, effZ_Zee, nZee_comb, nZee_diff;
  double Nz_Zmm, effZ_Zmm, nZmm_comb, nZmm_diff;

  ifstream ifs;
  ifs.open(infilename.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;

    if(state==0) {
      // Information for Zee
      stringstream ss(line);
      string label;
      ss >> label >> rZee >> Nz_init_Zee >> effGsf >> effHLT_Zee >> effSel_Zee >> nZee;
      state++;
      Nz_Zee = Nz_init_Zee*rZee;
      effZ_Zee = (1-(1-effHLT_Zee)*(1-effHLT_Zee))*effGsf*effGsf*effSel_Zee*effSel_Zee;
      nZee_comb = effZ_Zee*Nz_Zee;
      nZee_diff = TMath::Abs(nZee-nZee_comb);
    } else if(state==1) {
      // Information for Zmm
      stringstream ss(line);
      string label;
      ss >> label >> rZmm >> Nz_init_Zmm >> effHLT_Zmm >> effSel_Zmm >> effSta >> effTrk >> nZmm;
      Nz_Zmm = Nz_init_Zmm*rZmm;
      effZ_Zmm = (1-(1-effHLT_Zmm)*(1-effHLT_Zmm))*effTrk*effTrk*effSta*effSta*effSel_Zmm*effSel_Zmm;
      nZmm_comb = effZ_Zmm*Nz_Zmm;
      nZmm_diff = TMath::Abs(nZmm-nZmm_comb);
      state++;
   }
  }
  ifs.close();

  cout<<endl;
  cout<<setprecision(3)<<fixed;

  cout<<setw(8)<<"Channel"<<endl;
  cout<<setw(16)<<"effHLT"<<setw(8)<<"effSel"<<setw(8)<<"effGsf"<<setw(8)<<"effSta"<<setw(8)<<"effTrk"<<setw(8)<<"effZ"<<setw(12)<<"Nz"<<endl;
  cout<<setw(8)<<"Zee"<<setw(8)<<effHLT_Zee<<setw(8)<<effSel_Zee<<setw(8)<<effGsf<<setw(8)<<"-"<<setw(8)<<"-"<<setw(8)<<effZ_Zee<<setw(12)<<Nz_Zee<<endl;
  cout<<setw(8)<<"Zmm"<<setw(8)<<effHLT_Zmm<<setw(8)<<effSel_Zmm<<setw(8)<<"-"<<setw(8)<<effSta<<setw(8)<<effTrk<<setw(8)<<effZ_Zmm<<setw(12)<<Nz_Zmm<<endl;

  cout<<endl;
  cout<<setprecision(0) << fixed;
  cout<<setw(25)<<"Yield (EWK code)"<<setw(20)<<"Yield (Combine)"<<setw(15)<<"Difference" << endl;
  cout<<setw(8)<<"Zee"<<setw(17)<<nZee<<setw(20)<<nZee_comb<<setw(15)<<nZee_diff<<endl;
  cout<<setw(8)<<"Zmm"<<setw(17)<<nZmm<<setw(20)<<nZmm_comb<<setw(15)<<nZmm_diff<<endl;
  cout<<endl;
}
