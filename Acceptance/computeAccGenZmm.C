//================================================================================================
// Not used for 13 TeV measurement.
//
// Compute Z->mumu acceptance at generator level
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

#include "../Utils/MyTools.hh"

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#endif


//=== MAIN MACRO ================================================================================================= 

void computeAccGenZmm(const TString conf      , // input file
                      const TString outputDir , // output directory
                      const TString outputName, // output filename
                      const bool    doDressed=0, // do dressed
                      const bool    doborn=0 ) {
  gBenchmark->Start("computeAccGenZmm");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  // bool doPtWeights = true;
  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 13;
  
  const Int_t NPDF = 100;
  const Int_t NQCD = 6;

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
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  
  TFile *infile=0;
  TTree *eventTree=0;
   
  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv ;
  vector<Double_t> nSelv  , accv  , accErrv  ;
  vector<Double_t> nSelBBv, accBBv, accErrBBv;
  vector<Double_t> nSelBEv, accBEv, accErrBEv; 
  vector<Double_t> nSelEEv, accEEv, accErrEEv;
  
  vector<Double_t> nEvtsv_pT, nSelv_pT, accv_pT;
  
  vector<Double_t> nEvtsv_QCD, nSelv_QCD;
  vector<Double_t> nEvtsv_PDF, nSelv_PDF;
  
  TString sqrts = "13TeV";
  if(conf.Contains("5")) sqrts = "5TeV";
  TFile *rf = new TFile("/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction/Z_pT/zPt_Normal"+sqrts+".root");
  TH1D *hh_diff = (TH1D*)rf->Get("hZptRatio");
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
    
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
    eventTree->SetBranchAddress("GenEvtInfo",  &gen);        TBranch *genBr  = eventTree->GetBranch("GenEvtInfo"); 
    eventTree->SetBranchAddress("GenParticle", &genPartArr); TBranch *partBr = eventTree->GetBranch("GenParticle");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBBv.push_back(0);
    nSelBEv.push_back(0);
    nSelEEv.push_back(0);
    nEvtsv_pT.push_back(0);
    nSelv_pT.push_back(0);
    
    for(int i=0;i<NQCD;++i) {nSelv_QCD.push_back(0);nEvtsv_QCD.push_back(0);}
    for(int i=0;i<NPDF;++i) {nSelv_PDF.push_back(0);nEvtsv_PDF.push_back(0);}
    
    TLorentzVector *vec  =new TLorentzVector(0,0,0,0);
    TLorentzVector *lep1 =new TLorentzVector(0,0,0,0);
    TLorentzVector *lep2 =new TLorentzVector(0,0,0,0);
    TLorentzVector *lep3 =new TLorentzVector(0,0,0,0);
    TLorentzVector *lep4 =new TLorentzVector(0,0,0,0);
    TLorentzVector *gph  =new TLorentzVector(0,0,0,0);
    
    //
    // loop over events
    //      
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      
      genBr->GetEntry(ientry);
      genPartArr->Clear(); partBr->GetEntry(ientry);
          
      Int_t lepq1=-99, lepq2=-99;
      
      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue; 
      toolbox::fillGenBorn(genPartArr, BOSON_ID, vec, lep1, lep2, lep3, lep4);
      
      double ptWeight = 1;
      for(int i = 0; i <= hh_diff->GetNbinsX();++i){
        if(vec->Pt() > hh_diff->GetBinLowEdge(i) && vec->Pt() < hh_diff->GetBinLowEdge(i+1)){ ptWeight = hh_diff->GetBinContent(i); break;}
      }
      
      if(doDressed){
        for(Int_t i=0; i<genPartArr->GetEntriesFast(); i++) {
          const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
          if(fabs(genloop->pdgId)!=22) continue;
          gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
          if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep3->Eta(),lep3->Phi())<0.1)  lep3->operator+=(*gph);
          if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep4->Eta(),lep4->Phi())<0.1)  lep4->operator+=(*gph);
        }
      }
      
      if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;
    
      Double_t weight=gen->weight;
      nEvtsv[ifile]+=weight;
      // nEvtsv[ifile]+=weight*ptWeight;
      nEvtsv_pT[ifile]+=weight*ptWeight;
      
      // -------------------------------------------------
      // I'm not sure which indexing is correct for Z's
      // -------------------------------------------------
      nEvtsv_QCD[0]+=weight*gen->lheweight[0];
      nEvtsv_QCD[1]+=weight*gen->lheweight[1];
      nEvtsv_QCD[2]+=weight*gen->lheweight[2];
      nEvtsv_QCD[3]+=weight*gen->lheweight[3];
      nEvtsv_QCD[4]+=weight*gen->lheweight[5];
      nEvtsv_QCD[5]+=weight*gen->lheweight[7];
      
      for(int npdf=0; npdf<NPDF; npdf++)   nEvtsv_PDF[npdf]+=weight*gen->lheweight[8+npdf];

      if(lep1->Pt()        < PT_CUT)  continue;
      if(lep2->Pt()        < PT_CUT)  continue;
      if(fabs(lep1->Eta()) > ETA_CUT) continue;
      if(fabs(lep2->Eta()) > ETA_CUT) continue;

      TLorentzVector dilep;

      if(doborn)
	dilep=(*lep1)+(*lep2);
      else
	dilep=(*lep3)+(*lep4);

      if(dilep.M()<MASS_LOW || dilep.M()>MASS_HIGH) continue;
      //std::cout << dilep.M() << " " << vec->M() << std::endl;
      
      Bool_t isB1 = (fabs(lep1->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      Bool_t isB2 = (fabs(lep2->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      
      nSelv[ifile]+=weight;
      // nSelv[ifile]+=weight*ptWeight;
      nSelv_pT[ifile]+=weight*ptWeight;
      if(isB1 && isB2)        nSelBBv[ifile]+=weight;
      else if(!isB1 && !isB2) nSelEEv[ifile]+=weight;
      else		                nSelBEv[ifile]+=weight;       
      
      nSelv_QCD[0]+=weight*gen->lheweight[0];
      nSelv_QCD[1]+=weight*gen->lheweight[1];
      nSelv_QCD[2]+=weight*gen->lheweight[2];
      nSelv_QCD[3]+=weight*gen->lheweight[3];
      nSelv_QCD[4]+=weight*gen->lheweight[5];
      nSelv_QCD[5]+=weight*gen->lheweight[7];
      for(int npdf=0; npdf<NPDF; npdf++) nSelv_PDF[npdf]+=weight*gen->lheweight[8+npdf];      
    }

    delete vec;
    delete lep1;
    delete lep2;
    delete lep3;
    delete lep4;
    delete gph;
    
    // compute acceptances
    accv.push_back  (nSelv[ifile]  /nEvtsv[ifile]); accErrv.push_back  (sqrt(accv  [ifile]*(1.-accv  [ifile])/nEvtsv[ifile]));
    accBBv.push_back(nSelBBv[ifile]/nEvtsv[ifile]); accErrBBv.push_back(sqrt(accBBv[ifile]*(1.-accBBv[ifile])/nEvtsv[ifile]));
    accBEv.push_back(nSelBEv[ifile]/nEvtsv[ifile]); accErrBEv.push_back(sqrt(accBEv[ifile]*(1.-accBEv[ifile])/nEvtsv[ifile]));
    accEEv.push_back(nSelEEv[ifile]/nEvtsv[ifile]); accErrEEv.push_back(sqrt(accEEv[ifile]*(1.-accEEv[ifile])/nEvtsv[ifile]));
    
    accv_pT.push_back  (nSelv_pT[ifile]  /nEvtsv_pT[ifile]);
    

    std::cout << "nselv " << nSelv[ifile] << "  nevtsv " << nEvtsv[ifile] << std::endl;
    
    delete infile;
    infile=0, eventTree=0;  
  }
  
  delete gen;
  
  
  char masterOutput[600];
  for(uint ifile = 0; ifile < fnamev.size(); ++ifile){// go through info per file
    sprintf(masterOutput,"%s/%s.txt",outputDir.Data(),outputName.Data());
    ofstream txtfile;
    txtfile.open(masterOutput);
    txtfile << "acc " << nSelv[ifile]/nEvtsv[ifile] << endl;    
    for(int j = 0; j < NPDF; ++j) txtfile << "pdf" << j << " " << nSelv_PDF[j]/nEvtsv_PDF[j] << endl;
    for(int j = 0; j < NQCD; ++j) txtfile << "qcd" << j << " " << nSelv_QCD[j]/nEvtsv_QCD[j] << endl;
    txtfile.close();
  }
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================    
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
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
    cout << "     b-b: " << setw(12) << nSelBBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBBv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    cout << "     b-e: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    cout << "     e-e: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrEEv[ifile] << endl;
    cout << "   total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile]   << endl;
    cout << " with pt: " << setw(12) << accv_pT[ifile] << endl;
    cout << " pt diff: " << setw(12) << 100*fabs(accv[ifile]/accv_pT[ifile] - 1 ) << endl;
    cout << endl;
  }
  
  char txtfname[300];
  sprintf(txtfname,"%s/gen.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << " Z -> mu mu" << endl;
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
    txtfile << "     b-b: " << setw(12) << nSelBBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBBv[ifile] << " +/- " << accErrBBv[ifile] << endl;
    txtfile << "     b-e: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    txtfile << "     b-e: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrEEv[ifile] << endl;
    txtfile << "   total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile]   << endl;
    txtfile << " with pt: " << setw(12) << accv_pT[ifile] << endl;
    txtfile << " pt diff: " << setw(12) << 100*fabs(accv[ifile]/accv_pT[ifile] - 1 ) << endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccGenZmm_Sys"); 
}
