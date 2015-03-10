#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void MergeNtuples(const TString input) 
{
  gBenchmark->Start("MergeNtuples");

  TString outfilename;          // output of merged files  
  vector<TString> infilenames;  // list input ntuple files to be stored
  
  // 
  // parse input file
  //  
  ifstream ifs;
  ifs.open(input.Data()); 
  assert(ifs.is_open());
  string line;
  getline(ifs,line); 
  outfilename = line;
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();

  TTree::SetMaxTreeSize(kMaxLong64);
        
  //
  // Combine TTrees from each file
  //
  TChain chain("Events");
  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Adding " << infilenames[ifile] << endl;
    chain.Add(infilenames[ifile]);
  }
  cout << "Merging..." << endl;
  chain.Merge(outfilename,"fast");
  std::cout << outfilename << " created!" << std::endl;
  
  gBenchmark->Show("MergeNtuples");
}
