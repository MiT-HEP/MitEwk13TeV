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

void calcAccUnc(TString input="test.txt") {

  ifstream ifs;
  ifs.open(input);
  assert(ifs.is_open());
  string line;

  Float_t nTotal, nPassPre, nPassPost, nAccPre, nAccPost;
  getline(ifs,line);
  stringstream ss(line);
  ss >> nTotal >> nPassPre >> nPassPost >> nAccPre >> nAccPost;

  vector<Float_t>  vTotal, vPassPre, vPassPost, vAccPre, vAccPost;
  while (getline(ifs,line)) {
    stringstream ss(line);
    Float_t  total,  passPre,  passPost,  accPre,  accPost;
    ss >> total >> passPre >> passPost >> accPre >> accPost;
    vTotal.push_back(total);
    vPassPre.push_back(passPre);
    vPassPost.push_back(passPost);
    vAccPre.push_back(accPre);
    vAccPost.push_back(accPost);
  }

  Float_t sum1=0, sum2=0, sum3=0, sum4=0;
  for (Int_t i=0; i<vAccPre.size()/2; i++) {
    sum1+=TMath::Max( TMath::Max(vAccPre[2*i]-nAccPre, vAccPre[2*i+1]-nAccPre), 0.0)*TMath::Max( TMath::Max(vAccPre[2*i]-nAccPre, vAccPre[2*i+1]-nAccPre), 0.0);
    sum2+=TMath::Max( TMath::Max(nAccPre-vAccPre[2*i], nAccPre-vAccPre[2*i+1]), 0.0)*TMath::Max( TMath::Max(nAccPre-vAccPre[2*i], nAccPre-vAccPre[2*i+1]), 0.0);
    sum3+=TMath::Max( TMath::Max(vAccPost[2*i]-nAccPost, vAccPost[2*i+1]-nAccPost), 0.0)*TMath::Max( TMath::Max(vAccPost[2*i]-nAccPost, vAccPost[2*i+1]-nAccPost), 0.0);
    sum4+=TMath::Max( TMath::Max(nAccPost-vAccPost[2*i], nAccPost-vAccPost[2*i+1]), 0.0)*TMath::Max( TMath::Max(nAccPost-vAccPost[2*i], nAccPost-vAccPost[2*i+1]), 0.0);
  }
  
  cout << "acc, pre: " << nAccPre << " + " << sqrt(sum1)/1.64485 << " - " << sqrt(sum2)/1.64485 << endl;
  cout << "rel unc: " << (sqrt(sum1)/1.64485)/nAccPre << ", " << (sqrt(sum2)/1.64485)/nAccPre << endl;
  cout << "acc, post: " << nAccPost << " + " << sqrt(sum3)/1.64485 << " - " << sqrt(sum4)/1.64485 << endl;
  cout << "rel unc: " << (sqrt(sum3)/1.64485)/nAccPost << ", " << (sqrt(sum4)/1.64485)/nAccPost << endl;
  
}
