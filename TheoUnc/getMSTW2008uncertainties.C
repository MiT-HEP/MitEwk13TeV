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

#include "pdfUncertMSTW2008.hh"

using namespace std;

#endif

void getMSTW2008uncertainties(TString indir="/Users/jaylawhorn/Desktop/wz-8tev-acc") {

  Double_t W_p_mu = xsecUncert(indir+"/wpm_parse_mstw2008nlo68cl.txt",
			       indir+"/wpm_parse_mstw2008nlo68cl_asmz+68cl.txt",
			       indir+"/wpm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
			       indir+"/wpm_parse_mstw2008nlo68cl_asmz-68cl.txt",
			       indir+"/wpm_parse_mstw2008nlo68cl_asmz-68clhalf.txt");

  Double_t W_m_mu = xsecUncert(indir+"/wmm_parse_mstw2008nlo68cl.txt",
			       indir+"/wmm_parse_mstw2008nlo68cl_asmz+68cl.txt",
			       indir+"/wmm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
			       indir+"/wmm_parse_mstw2008nlo68cl_asmz-68cl.txt",
			       indir+"/wmm_parse_mstw2008nlo68cl_asmz-68clhalf.txt");
  
  Double_t W_pOverM_mu = 
    chRatioUncert(indir+"/wpm_parse_mstw2008nlo68cl.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz-68clhalf.txt");
  

  Double_t W_mu = 
    chSumUncert(indir+"/wpm_parse_mstw2008nlo68cl.txt",
		indir+"/wpm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		indir+"/wpm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		indir+"/wpm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		indir+"/wpm_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		indir+"/wmm_parse_mstw2008nlo68cl.txt",
		indir+"/wmm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		indir+"/wmm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		indir+"/wmm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		indir+"/wmm_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		5.82288,
		3.94813);

  Double_t Z_mumu = xsecUncert(indir+"/zmm_parse_mstw2008nlo68cl.txt",
			       indir+"/zmm_parse_mstw2008nlo68cl_asmz+68cl.txt",
			       indir+"/zmm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
			       indir+"/zmm_parse_mstw2008nlo68cl_asmz-68cl.txt",
			       indir+"/zmm_parse_mstw2008nlo68cl_asmz-68clhalf.txt");

  Double_t Z_W_mu = 
    wzRatioUncert(indir+"/wpm_parse_mstw2008nlo68cl.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wpm_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wmm_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  indir+"/zmm_parse_mstw2008nlo68cl.txt",
		  indir+"/zmm_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/zmm_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/zmm_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/zmm_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  5.82288,
		  3.94813);

  Double_t W_p_el = xsecUncert(indir+"/wpe_parse_mstw2008nlo68cl.txt",
			       indir+"/wpe_parse_mstw2008nlo68cl_asmz+68cl.txt",
			       indir+"/wpe_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
			       indir+"/wpe_parse_mstw2008nlo68cl_asmz-68cl.txt",
			       indir+"/wpe_parse_mstw2008nlo68cl_asmz-68clhalf.txt");

  Double_t W_m_el = xsecUncert(indir+"/wme_parse_mstw2008nlo68cl.txt",
			       indir+"/wme_parse_mstw2008nlo68cl_asmz+68cl.txt",
			       indir+"/wme_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
			       indir+"/wme_parse_mstw2008nlo68cl_asmz-68cl.txt",
			       indir+"/wme_parse_mstw2008nlo68cl_asmz-68clhalf.txt");

  Double_t W_pOverM_el = 
    chRatioUncert(indir+"/wpe_parse_mstw2008nlo68cl.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  indir+"/wme_parse_mstw2008nlo68cl.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz-68clhalf.txt");

  Double_t W_el = 
    chSumUncert(indir+"/wpe_parse_mstw2008nlo68cl.txt",
		indir+"/wpe_parse_mstw2008nlo68cl_asmz+68cl.txt",
		indir+"/wpe_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		indir+"/wpe_parse_mstw2008nlo68cl_asmz-68cl.txt",
		indir+"/wpe_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		indir+"/wme_parse_mstw2008nlo68cl.txt",
		indir+"/wme_parse_mstw2008nlo68cl_asmz+68cl.txt",
		indir+"/wme_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		indir+"/wme_parse_mstw2008nlo68cl_asmz-68cl.txt",
		indir+"/wme_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		5.82288,
		3.94813);

  Double_t Z_elel = xsecUncert(indir+"/zee_parse_mstw2008nlo68cl.txt",
			       indir+"/zee_parse_mstw2008nlo68cl_asmz+68cl.txt",
			       indir+"/zee_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
			       indir+"/zee_parse_mstw2008nlo68cl_asmz-68cl.txt",
			       indir+"/zee_parse_mstw2008nlo68cl_asmz-68clhalf.txt");

  Double_t Z_W_el = 
    wzRatioUncert(indir+"/wpe_parse_mstw2008nlo68cl.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wpe_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  indir+"/wme_parse_mstw2008nlo68cl.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/wme_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  indir+"/zee_parse_mstw2008nlo68cl.txt",
		  indir+"/zee_parse_mstw2008nlo68cl_asmz+68cl.txt",
		  indir+"/zee_parse_mstw2008nlo68cl_asmz+68clhalf.txt",
		  indir+"/zee_parse_mstw2008nlo68cl_asmz-68cl.txt",
		  indir+"/zee_parse_mstw2008nlo68cl_asmz-68clhalf.txt",
		  5.82288,
		  3.94813);

  cout << "W+    m " << W_p_mu << endl;
  cout << "W-    m " << W_m_mu << endl;
  cout << "W+/W- m " << W_pOverM_mu << endl;
  cout << "W     m " << W_mu << endl;
  cout << "Z     m " << Z_mumu << endl;
  cout << "Z/W   m " << Z_W_mu << endl;
  cout << endl;
  cout << "W+    e " << W_p_el << endl;
  cout << "W-    e " << W_m_el << endl;
  cout << "W+/W- e " << W_pOverM_el << endl;
  cout << "W     e " << W_el << endl;
  cout << "Z     e " << Z_elel << endl;
  cout << "Z/W   e " << Z_W_el << endl;
  cout << endl;

}
