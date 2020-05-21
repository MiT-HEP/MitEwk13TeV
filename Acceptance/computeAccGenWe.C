//================================================================================================
// Not used for 13 TeV measurement.
//
// Compute W->enu acceptance at generator level
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
#include <TLorentzVector.h>

#include "../Utils/MyTools.hh"

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#endif

//=== MAIN MACRO ================================================================================================= 

void computeAccGenWe(const TString conf,       // input file
                     const TString outputDir, // output directory
                     const TString outputName, // output filename
                     const bool doDressed=0,
		     const Int_t   charge=0,      // 0 = inclusive, +1 = W+, -1 = W-
		     const bool doborn=0
		     ) {
  gBenchmark->Start("computeAccGenWe");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t ETA_BARREL = 1.4442;
  const Double_t ETA_ENDCAP = 1.566;
  // const Double_t ETA_BARREL = 10.;
  // const Double_t ETA_ENDCAP = 10.;

  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 11;
  
  const Int_t NPDF = 100;
  const Int_t NQCD = 6;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString> fnamev;  // file name per input file
  vector<TString> labelv;  // TLegend label per input file
  vector<Int_t>   xsecv;  // plot color per input file
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
    Int_t xsec, linesty;
    stringstream ss(line);
    ss >> fname >> xsec >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    xsecv.push_back(xsec);
    linev.push_back(linesty);
  }
  ifs.close();
  
  int NFILES=fnamev.size();

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  
  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr     = new TClonesArray("baconhep::TGenParticle");
  
  TFile *infile=0;
  TTree *eventTree=0;

  // Variables to store acceptances and uncertainties (per input file)
  vector<Double_t> nEvtsv, nSelv, nSelBv, nSelEv;
  vector<Double_t> accv, accBv, accEv;
  vector<Double_t> accErrv, accErrBv, accErrEv;
  
  
  vector<Double_t> nEvtsv_pT, nSelv_pT, accv_pT;
  
  vector<vector<Double_t>> nEvtsv_QCD, nSelv_QCD;
  vector<vector<Double_t>> nEvtsv_PDF, nSelv_PDF;
  
  TString sqrts = "13TeV";
  if(conf.Contains("5")) sqrts = "5TeV";
  TFile *rf = new TFile("/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction/Z_pT/zPt_Normal"+sqrts+".root");
  TH1D *hh_diff = (TH1D*)rf->Get("hZptRatio");
  //
  
  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    
    std::cout << "cross section is ... "  << xsecv[ifile] << std::endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events");
    assert(eventTree);
    eventTree->SetBranchAddress("GenEvtInfo", &gen); TBranch *genBr = eventTree->GetBranch("GenEvtInfo");
    eventTree->SetBranchAddress("GenParticle", &genPartArr); TBranch *partBr = eventTree->GetBranch("GenParticle");
    std::cout << "blah" << std::endl;
    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBv.push_back(0);
    nSelEv.push_back(0);
    nEvtsv_pT.push_back(0);
    nSelv_pT.push_back(0);
    
    vector<Double_t> tempQCD_Selv, tempQCD_Evtsv;
    vector<Double_t> tempPDF_Selv, tempPDF_Evtsv;
    
    for(int i=0;i<NQCD;++i) {tempQCD_Selv.push_back(0);tempQCD_Evtsv.push_back(0);}
    for(int i=0;i<NPDF;++i) {tempPDF_Selv.push_back(0);tempPDF_Evtsv.push_back(0);}
    
    TLorentzVector *vec=new TLorentzVector(0,0,0,0);
    TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
    TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
    TLorentzVector *lep3 =new TLorentzVector(0,0,0,0);
    TLorentzVector *lep4 =new TLorentzVector(0,0,0,0);
    TLorentzVector *gph  =new TLorentzVector(0,0,0,0);

    // loop over events
    //
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); partBr->GetEntry(ientry);
      
      Int_t lepq1=-99;
      Int_t lepq2=-99;
      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      if (charge== -1 &&     toolbox::flavor(genPartArr, BOSON_ID) != LEPTON_ID) continue; // check for a e- from W
      if (charge==  1 &&     toolbox::flavor(genPartArr, BOSON_ID) !=-LEPTON_ID) continue; // check for a e+ from W
      if (charge==  0 && fabs(toolbox::flavor(genPartArr, BOSON_ID))!= LEPTON_ID) continue; // check flavor

      // the function returns: lep1, lep3 are the paricles, lep2, lep4 are the anti-particles
      toolbox::fillGenBorn(genPartArr, BOSON_ID, vec, lep1, lep2, lep3, lep4);

      double ptWeight = 1;
      for(int i = 0; i <= hh_diff->GetNbinsX();++i){
        if(vec->Pt() > hh_diff->GetBinLowEdge(i) && vec->Pt() < hh_diff->GetBinLowEdge(i+1)){ ptWeight = hh_diff->GetBinContent(i); break;}
      }

      if(charge==1) {
	//For W+->e+vu decay, change things up so that lep1 and lep3 are the charged particles
        TLorentzVector *tmp = lep1;
        lep1=lep2;
        lep2=tmp;
	tmp=lep3;
	lep3=lep4;
	lep4=tmp;
      }
      
      
      if(doDressed){
	for(Int_t i=0; i<genPartArr->GetEntriesFast(); i++) {
	  const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
	  if(fabs(genloop->pdgId)!=22) continue;
	  gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	  if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep3->Eta(),lep3->Phi())<0.1)
	    {
	      lep3->operator+=(*gph);
	    }
	}
      }
      
      Double_t weight=gen->weight;
      nEvtsv[ifile]+=weight;
      // nEvtsv[ifile]+=weight*ptWeight;
      nEvtsv_pT[ifile]+=weight*ptWeight;
      
    
      tempQCD_Evtsv[0]+=weight*gen->lheweight[1];
      tempQCD_Evtsv[1]+=weight*gen->lheweight[2];
      tempQCD_Evtsv[2]+=weight*gen->lheweight[3];
      tempQCD_Evtsv[3]+=weight*gen->lheweight[4];
      tempQCD_Evtsv[4]+=weight*gen->lheweight[6];
      tempQCD_Evtsv[5]+=weight*gen->lheweight[8];
      for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Evtsv[npdf]+=weight*gen->lheweight[9+npdf];

      for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Evtsv[npdf]+=weight*gen->lheweight[9+npdf];
    
      Bool_t isBarrel=kTRUE;
      if (doborn) 
	{
	  if (lep1->Pt()        < PT_CUT ) continue;
	  if (fabs(lep1->Eta()) > ETA_CUT) continue;
	  isBarrel = (fabs(lep1->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
	}
      else 
	{
	  if (lep3->Pt()        < PT_CUT ) continue;
	  if (fabs(lep3->Eta()) > ETA_CUT) continue;
	  isBarrel = (fabs(lep3->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
	}

      nSelv[ifile]+=weight;
      // nSelv[ifile]+=weight*ptWeight;
      nSelv_pT[ifile]+=weight*ptWeight;
      if(isBarrel) nSelBv[ifile]+=weight;
      else	       nSelEv[ifile]+=weight;      
      
      tempQCD_Selv[0]+=weight*gen->lheweight[1];
      tempQCD_Selv[1]+=weight*gen->lheweight[2];
      tempQCD_Selv[2]+=weight*gen->lheweight[3];
      tempQCD_Selv[3]+=weight*gen->lheweight[4];
      tempQCD_Selv[4]+=weight*gen->lheweight[6];
      tempQCD_Selv[5]+=weight*gen->lheweight[8];
      for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Selv[npdf]+=weight*gen->lheweight[9+npdf];
    }
    
    delete vec;
    delete lep1;
    delete lep2;
    delete lep3;
    delete lep4;
    delete gph;


    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.-accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.-accEv[ifile])/nEvtsv[ifile]));
    nSelv_PDF.push_back(tempPDF_Selv);
    nEvtsv_PDF.push_back(tempPDF_Evtsv);
    nSelv_QCD.push_back(tempQCD_Selv);
    nEvtsv_QCD.push_back(tempQCD_Evtsv);
    
    accv_pT.push_back  (nSelv_pT[ifile]  /nEvtsv_pT[ifile]);
    
    std::cout << "nselv " << nSelv[ifile] << "  nevtsv " << nEvtsv[ifile] << std::endl;
    
  
    
    delete infile;
    infile=0, eventTree=0;  
  }
  // vector<Double_t> accv_PDF
  double accTot=0;
  double accNum=0, accDnm=0;
  std::cout << "here" << std::endl;
  
  // Print full set for efficiency calculations
  char masterOutput[600];
  // just start printing....
  for(uint ifile = 0; ifile < fnamev.size(); ++ifile){// go through info per file
    sprintf(masterOutput,"%s/%s_%d.txt",outputDir.Data(),outputName.Data(),ifile);
    ofstream txtfile;
    txtfile.open(masterOutput);
    txtfile << "acc " << nSelv[ifile]/nEvtsv[ifile] << endl;
    
    for(int j = 0; j < NPDF; ++j) txtfile << "pdf" << j << " " << nSelv_PDF[ifile][j]/nEvtsv_PDF[ifile][j] << endl;
    for(int j = 0; j < NQCD; ++j) txtfile << "qcd" << j << " " << nSelv_QCD[ifile][j]/nEvtsv_QCD[ifile][j] << endl;
    txtfile.close();
  }

  delete gen;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  if(charge== 0) cout << " W -> e nu"  << endl;
  if(charge==-1) cout << " W- -> e nu" << endl;
  if(charge== 1) cout << " W+ -> e nu" << endl;
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
    cout << "            b: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile] << endl;
    cout << "            e: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile] << endl;
    cout << "        total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile]  << endl;
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
  if(charge== 0) txtfile << " W -> e nu"  << endl;
  if(charge==-1) txtfile << " W- -> e nu" << endl;
  if(charge== 1) txtfile << " W+ -> e nu" << endl;
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
    txtfile << "            b: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile] << endl;
    txtfile << "            e: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile] << endl;
    txtfile << "        total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile]  << endl;
    txtfile << " with pt: " << setw(12) << accv_pT[ifile] << endl;
    txtfile << " pt diff: " << setw(12) << 100*fabs(accv[ifile]/accv_pT[ifile] - 1 ) << endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccGenWe_Sys"); 
}
