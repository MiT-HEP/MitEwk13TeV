//================================================================================================
//
// Compute Z->ee acceptance at full selection level
//
//  * outputs results summary text file
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TH1D.h>                   // histogram class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "Math/LorentzVector.h"     // 4-vector class

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"

// helper functions for lepton ID selection
#include "EWKAna/Utils/LeptonIDCuts.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void computeAccSelZee(const TString conf,             // input file
                      const TString outputDir         // output directory
) {
  gBenchmark->Start("computeAccSelZee");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.5;
  const Double_t ELE_MASS   = 0.000511;  
  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;
  
  // Data/MC Z efficiency scale factor from simultaneous fit
  Double_t zEffScale = 0.988542, zEffScaleErr = 0.011674;
    

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   colorv;  // plot color per input file
  vector<Int_t>   linev;   // plot line style per input file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    ss >> fname >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
  }
  ifs.close();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);

  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo   *gen   = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  
  TFile *infile=0;
  TTree *eventTree=0;
  
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBBv, nSelBEv, nSelEEv;
  vector<Double_t> accv, accBBv, accBEv, accEEv;
  vector<Double_t> accErrv, accErrBBv, accErrBEv, accErrEEv;
    
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",      &gen);         TBranch *genBr      = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBBv.push_back(0);
    nSelBEv.push_back(0);
    nSelEEv.push_back(0);
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      if(gen->vmass<MASS_LOW || gen->vmass>MASS_HIGH) continue;

      infoBr->GetEntry(ientry);     
    
      Double_t weight=1;
      nEvtsv[ifile]+=weight;
    
      // trigger requirement               
      ULong_t trigger = kHLT_Ele22_CaloIdL_CaloIsoVL;
      ULong_t trigObj = kHLT_Ele22_CaloIdL_CaloIsoVL_EleObj;   
      if(!(info->triggerBits & trigger)) continue;
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
    
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      for(Int_t i1=0; i1<electronArr->GetEntriesFast(); i1++) {
  	const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[i1]);
	
	// check ECAL gap
	if(fabs(ele1->scEta)>=ETA_BARREL && fabs(ele1->scEta)<=ETA_ENDCAP) continue;
        
	if(ele1->scEt	     < PT_CUT)	     continue;  // lepton pT cut
        if(fabs(ele1->scEta) > ETA_CUT)	     continue;  // lepton |eta| cut
        if(!passEleID(ele1,info->rhoLowEta)) continue;  // lepton selection

        LorentzVector vEle1(ele1->pt, ele1->eta, ele1->phi, ELE_MASS);
	Bool_t isB1 = (fabs(ele1->scEta)<ETA_BARREL) ? kTRUE : kFALSE;

        for(Int_t i2=i1+1; i2<electronArr->GetEntriesFast(); i2++) {          
	  const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[i2]);
	  
	  // check ECAL gap
	  if(fabs(ele2->scEta)>=ETA_BARREL && fabs(ele2->scEta)<=ETA_ENDCAP) continue;
        
          if(ele2->scEt        < PT_CUT)       continue;  // lepton pT cut
          if(fabs(ele2->scEta) > ETA_CUT)      continue;  // lepton |eta| cut
	  if(!passEleID(ele2,info->rhoLowEta)) continue;  // lepton selection

          LorentzVector vEle2(ele2->pt, ele2->eta, ele2->phi, ELE_MASS);  
          Bool_t isB2 = (fabs(ele2->scEta)<ETA_BARREL) ? kTRUE : kFALSE;

          // trigger match
	  if(!(ele1->hltMatchBits & trigObj) && !(ele2->hltMatchBits & trigObj)) continue;
	  
	  // mass window
          LorentzVector vDilep = vEle1 + vEle2;
          if((vDilep.M()<MASS_LOW) || (vDilep.M()>MASS_HIGH)) continue;
          
          
          /******** We have a Z candidate! HURRAY! ********/
       
          nSelv[ifile]+=weight;
  	  if(isB1 && isB2)        nSelBBv[ifile]+=weight;
	  else if(!isB1 && !isB2) nSelEEv[ifile]+=weight;
	  else                    nSelBEv[ifile]+=weight;
        }
      }      
    }
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);     accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBBv.push_back(nSelBBv[ifile]/nEvtsv[ifile]); accErrBBv.push_back(sqrt(accBBv[ifile]*(1.-accBBv[ifile])/nEvtsv[ifile]));
    accBEv.push_back(nSelBEv[ifile]/nEvtsv[ifile]); accErrBEv.push_back(sqrt(accBEv[ifile]*(1.-accBEv[ifile])/nEvtsv[ifile]));
    accEEv.push_back(nSelEEv[ifile]/nEvtsv[ifile]); accErrEEv.push_back(sqrt(accEEv[ifile]*(1.-accEEv[ifile])/nEvtsv[ifile]));
    
    delete infile;
    infile=0, eventTree=0;  
  }  
  delete info;
  delete gen;
  delete electronArr;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> e e" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  cout << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "   ================================================" << endl;
    cout << "    Label: " << labelv[ifile] << endl;
    cout << "     File: " << fnamev[ifile] << endl;
    cout << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     barrel-barrel: " << setw(12) << nSelBBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBBv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    cout << "     barrel-endcap: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    cout << "     endcap-endcap: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    cout << "             total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    cout << "     efficiency corrected: " << accv[ifile]*zEffScale;
    cout << " +/- " << accv[ifile]*zEffScale*sqrt(accErrv[ifile]*accErrv[ifile]/accv[ifile]/accv[ifile] + zEffScaleErr*zEffScaleErr/zEffScale/zEffScale) << endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/sel.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> e e" << endl;
  txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  txtfile << "  pT > " << PT_CUT << endl;
  txtfile << "  |eta| < " << ETA_CUT << endl;
  txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
  txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
  txtfile << endl;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "   ================================================" << endl;
    txtfile << "    Label: " << labelv[ifile] << endl;
    txtfile << "     File: " << fnamev[ifile] << endl;
    txtfile << endl;
    txtfile << "    *** Acceptance ***" << endl;
    txtfile << "     barrel-barrel: " << setw(12) << nSelBBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBBv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    txtfile << "     barrel-endcap: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    txtfile << "     endcap-endcap: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    txtfile << "             total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    txtfile << "     efficiency corrected: " << accv[ifile]*zEffScale;
    txtfile << " +/- " << accv[ifile]*zEffScale*sqrt(accErrv[ifile]*accErrv[ifile]/accv[ifile]/accv[ifile] + zEffScaleErr*zEffScaleErr/zEffScale/zEffScale) << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccSelZee"); 
}
