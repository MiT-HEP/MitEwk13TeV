#ifndef CONFPARSE_HH
#define CONFPARSE_HH

#include <vector>                   // STL vector class
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TString.h>                // ROOT string class
#include "MitEwk13TeV/Utils/CSample.hh"  // helper class to handle samples

void confParse(const TString    conf,      // input conf file
               vector<TString>  &snamev,   // vector to store output of sample names
	       vector<CSample*> &samplev   // vector to store output of sample handles
) {
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    if(line[0]=='$') {
      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color, linecol;
      ss >> chr >> sname >> color >> linecol;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      samplev.back()->linecol = linecol;
      continue;
    }
    
    string fname;
    string json;
    Double_t xsec;
    stringstream ss(line);
    ss >> fname >> xsec >> json;
    samplev.back()->fnamev.push_back(fname);
    samplev.back()->typev.push_back(0);
    samplev.back()->xsecv.push_back(xsec);
    samplev.back()->jsonv.push_back(json);
  }
  ifs.close();
}
#endif
