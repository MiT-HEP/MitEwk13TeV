//================================================================================================
// Not used for 13 TeV measurement.
//
// Compute W->munu acceptance at generator level
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

void computeAccGenWm_Sys(const TString conf,             // input file
                     const TString outputDir,        // output directory
                     const bool doDressed=0,
		     const Int_t   charge =0           // 0 = inclusive, +1 = W+, -1 = W-
) {
  gBenchmark->Start("computeAccGenWm_Sys");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4;
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;
 
  const Int_t BOSON_ID  = 24;
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
  vector<Double_t> nEvtsv, nSelv, nSelBv, nSelEv;
  vector<Double_t> accv, accBv, accEv;
  vector<Double_t> accErrv, accErrBv, accErrEv;
  
  vector<Double_t> nEvtsv_QCD, nSelv_QCD;
  vector<Double_t> accv_QCD;
  vector<Double_t> accErrv_QCD;
  
  vector<Double_t> nEvtsv_PDF, nSelv_PDF;
  vector<Double_t> accv_PDF;
  vector<Double_t> accErrv_PDF;
  
  // vector<vector<Double_t>> nEvtsv_QCD, nSelv_QCD;
  // vector<vector<Double_t>> accv_QCD;
  // vector<vector<Double_t>> accErrv_QCD;
  
  // vector<vector<Double_t>> nEvtsv_PDF, nSelv_PDF;
  // vector<vector<Double_t>> accv_PDF;
  // vector<vector<Double_t>> accErrv_PDF;

  double accv_uncPDF = 0, accv_uncPDF_num = 0, accv_uncPDF_dnm = 0;
  double accv_uncQCD = 0, accv_uncQCD_num = 0, accv_uncQCD_dnm = 0;
  double totalXsec = 0;

  //
  // loop through files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << " ..." << endl;
    infile = TFile::Open(fnamev[ifile]); 
    assert(infile);
  
    eventTree = (TTree*)infile->Get("Events"); 
    assert(eventTree);
    eventTree->SetBranchAddress("GenEvtInfo", &gen); TBranch *genBr = eventTree->GetBranch("GenEvtInfo"); 
    eventTree->SetBranchAddress("GenParticle", &genPartArr); TBranch *partBr = eventTree->GetBranch("GenParticle");

    nEvtsv.push_back(0);
    nSelv.push_back(0);
    nSelBv.push_back(0);
    nSelEv.push_back(0);
    
    for(int i=0;i<NQCD;++i) {nSelv_QCD.push_back(0);nEvtsv_QCD.push_back(0);}
    for(int i=0;i<NPDF;++i) {nSelv_PDF.push_back(0);nEvtsv_PDF.push_back(0);}
    
    //
    // loop over events
    //    
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; ientry<(uint)(0.1*eventTree->GetEntries()); ientry++) {
      if(ientry%1000000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); partBr->GetEntry(ientry);

      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
	  TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
	  TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      Int_t lepq1=-99;
	  Int_t lepq2=-99;
      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      if (charge==-1 && toolbox::flavor(genPartArr, BOSON_ID)!=LEPTON_ID) continue;
      if (charge==1 && toolbox::flavor(genPartArr, BOSON_ID)!=-LEPTON_ID) continue;
      if (charge==0 && fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
      // int mparam = fabs(1-fabs(charge));
      toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,&lepq1,&lepq2,1);
      
      TLorentzVector *gph=new TLorentzVector(0,0,0,0);
      if(doDressed){
          for(Int_t i=0; i<genPartArr->GetEntriesFast(); i++) {
            const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
            if(fabs(genloop->pdgId)!=22) continue;
            gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep1->Eta(),lep1->Phi())<0.1)
              {
            lep1->operator+=(*gph);
              }
            if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep2->Eta(),lep2->Phi())<0.1)
              {
            lep2->operator+=(*gph);
              }
          }
      }

      Double_t weight=gen->weight;
      nEvtsv[ifile]+=weight;
      // nEvtsv_QCD[ifile][0]+=weight*gen->lheweight[1];
      // nEvtsv_QCD[ifile][1]+=weight*gen->lheweight[2];
      // nEvtsv_QCD[ifile][2]+=weight*gen->lheweight[3];
      // nEvtsv_QCD[ifile][3]+=weight*gen->lheweight[4];
      // nEvtsv_QCD[ifile][4]+=weight*gen->lheweight[6];
      // nEvtsv_QCD[ifile][5]+=weight*gen->lheweight[8];
      // for(int npdf=0; npdf<NPDF; npdf++) nEvtsv_PDF[ifile][npdf]+=weight*gen->lheweight[9+npdf];
      nEvtsv_QCD[0]+=weight*gen->lheweight[1];
      nEvtsv_QCD[1]+=weight*gen->lheweight[2];
      nEvtsv_QCD[2]+=weight*gen->lheweight[3];
      nEvtsv_QCD[3]+=weight*gen->lheweight[4];
      nEvtsv_QCD[4]+=weight*gen->lheweight[6];
      nEvtsv_QCD[5]+=weight*gen->lheweight[8];
      for(int npdf=0; npdf<NPDF; npdf++) nEvtsv_PDF[npdf]+=weight*gen->lheweight[9+npdf];
    
      Bool_t isBarrel=kTRUE;
      if (lep1) {
        if (fabs(lep1->Eta())>ETA_BARREL && fabs(lep1->Eta())<ETA_ENDCAP) continue;
	
        if (lep1->Pt() < PT_CUT) continue;
        if (fabs(lep1->Eta()) > ETA_CUT) continue;
        isBarrel = (fabs(lep1->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      }
      else if (lep2) {
        if (fabs(lep2->Eta())>ETA_BARREL && fabs(lep2->Eta())<ETA_ENDCAP) continue;

        if (lep2->Pt() < PT_CUT) continue;
        if (fabs(lep2->Eta()) > ETA_CUT) continue;
        isBarrel = (fabs(lep2->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      }
      else continue;

      nSelv[ifile]+=weight;
      if(isBarrel) nSelBv[ifile]+=weight;
      else         nSelEv[ifile]+=weight;
            // Get the values with the QCD and PDF weights:
      // QCD first
      // nSelv_QCD[ifile][0]+=weight*gen->lheweight[1];
      // nSelv_QCD[ifile][1]+=weight*gen->lheweight[2];
      // nSelv_QCD[ifile][2]+=weight*gen->lheweight[3];
      // nSelv_QCD[ifile][3]+=weight*gen->lheweight[4];
      // nSelv_QCD[ifile][4]+=weight*gen->lheweight[6];
      // nSelv_QCD[ifile][5]+=weight*gen->lheweight[8];
      // for(int npdf=0; npdf<NPDF; npdf++) nSelv_PDF[ifile][npdf]+=weight*gen->lheweight[9+npdf];
      nSelv_QCD[0]+=weight*gen->lheweight[1];
      nSelv_QCD[1]+=weight*gen->lheweight[2];
      nSelv_QCD[2]+=weight*gen->lheweight[3];
      nSelv_QCD[3]+=weight*gen->lheweight[4];
      nSelv_QCD[4]+=weight*gen->lheweight[6];
      nSelv_QCD[5]+=weight*gen->lheweight[8];
      for(int npdf=0; npdf<NPDF; npdf++) nSelv_PDF[npdf]+=weight*gen->lheweight[9+npdf];
    }
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.-accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.-accEv[ifile])/nEvtsv[ifile]));
    
     std::cout << "nselv " << nSelv[ifile] << "  nevtsv " << nEvtsv[ifile] << std::endl;
    
    
    for(int npdf=0; npdf<NPDF; npdf++){
        accv_PDF.push_back(nSelv_PDF[npdf]/nEvtsv_PDF[npdf]);     
        accErrv_PDF.push_back(sqrt(accv_PDF[npdf]*(1.-accv_PDF[npdf])/nEvtsv_PDF[npdf]));
        std::cout << "nselvpdf " << nSelv_PDF[npdf] << "  nevtsvpdf " << nEvtsv_PDF[npdf] << " ratio " << nSelv_PDF[npdf]/nEvtsv_PDF[npdf] << std::endl;
        std::cout << "accv " << accv[ifile] << "  accvpdf " << accv_PDF[npdf] << std::endl;
        std::cout << "diff " << accv_PDF[npdf]-accv[ifile] << "  pct diff " << 100*(accv_PDF[npdf]-accv[ifile])/accv[ifile] << std::endl;
        // a = sum(((acc_i-accTot)/accTot)^2)
        // unc = sqrt(a/NPDF)
        accv_uncPDF+=(accv_PDF[npdf]-accv[ifile])*(accv_PDF[npdf]-accv[ifile])/(NPDF*accv[ifile]*accv[ifile]);
        accv_uncPDF_num+=(nSelv_PDF[npdf]-nSelv[ifile])*(nSelv_PDF[npdf]-nSelv[ifile])/(NPDF*nSelv[ifile]*nSelv[ifile]);
        accv_uncPDF_dnm+=(nEvtsv_PDF[npdf]-nEvtsv[ifile])*(nEvtsv_PDF[npdf]-nEvtsv[ifile])/(NPDF*nEvtsv[ifile]*nEvtsv[ifile]);
    }
    accv_uncPDF=sqrt(accv_uncPDF);
    accv_uncPDF_num=sqrt(accv_uncPDF_num);
    accv_uncPDF_dnm=sqrt(accv_uncPDF_dnm);
    for(int nqcd=0; nqcd<NQCD; nqcd++){
        accv_QCD.push_back(nSelv_QCD[nqcd]/nEvtsv_QCD[nqcd]);
        // std::cout << "nselvqcd " << nSelv_QCD[nqcd] << "  nevtsvqcd " << nEvtsv_QCD[nqcd] << " ratio " << nSelv_QCD[nqcd]/nEvtsv_QCD[nqcd] << std::endl;
        // std::cout << "accv " << accv[ifile] << "  accvqcd " << accv_QCD[nqcd] << std::endl;
        std::cout << "qcd a " << nqcd <<" " << accv_QCD[nqcd] << "  nom " << accv[ifile] << "  pct diff " << 100*fabs(accv_QCD[nqcd]-accv[ifile])/(accv[ifile]) << std::endl;
        std::cout << "qcd n " << nqcd <<" " << nSelv_QCD[nqcd] << "  nom " << nSelv[ifile] << "  pct diff n " << 100*fabs(nSelv_QCD[nqcd]-nSelv[ifile])/nSelv[ifile] << std::endl;
        std::cout << "qcd d " << nqcd <<" " << nEvtsv_QCD[nqcd] << "  nom " << nEvtsv[ifile] << "  pct diff d " << 100*fabs(nEvtsv_QCD[nqcd]-nEvtsv[ifile])/nEvtsv[ifile] << std::endl;
        accErrv_QCD.push_back(sqrt(accv_QCD[nqcd]*(1.-accv_QCD[nqcd])/nEvtsv_QCD[nqcd]));
        // a = abs((acc_i-accTot)/accTot)
        if(fabs(accv_QCD[nqcd]-accv[ifile])/(accv[ifile]) > accv_uncQCD) accv_uncQCD = fabs(accv_QCD[nqcd]-accv[ifile])/(accv[ifile]);
        if(fabs(nSelv_QCD[nqcd]-nSelv[ifile])/nSelv[ifile] > accv_uncQCD_num) accv_uncQCD_num = fabs(nSelv_QCD[nqcd]-nSelv[ifile])/nSelv[ifile];
        if(fabs(nEvtsv_QCD[nqcd]-nEvtsv[ifile])/nEvtsv[ifile] > accv_uncQCD_dnm) accv_uncQCD_dnm = fabs(nEvtsv_QCD[nqcd]-nEvtsv[ifile])/nEvtsv[ifile];
    }
    
    
    delete infile;
    infile=0, eventTree=0;  
  }

  delete gen;
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  if(charge== 0) cout << " W -> mu nu"  << endl;
  if(charge==-1) cout << " W- -> mu nu" << endl;
  if(charge== 1) cout << " W+ -> mu nu" << endl;
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
    cout << "            barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile] << endl;
    cout << "            endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile] << endl;
    cout << "             total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile] << endl;
    cout << "   qcd uncertainty: " << setw(12) << accv_uncQCD << std::endl;
    cout << " qcd uncertainty n: " << setw(12) << accv_uncQCD_num << std::endl;
    cout << " qcd uncertainty d: " << setw(12) << accv_uncQCD_dnm << std::endl;
    cout << "   pdf uncertainty: " << setw(12) << accv_uncPDF << std::endl;
    cout << " pdf uncertainty n: " << setw(12) << accv_uncPDF_num << std::endl;
    cout << " pdf uncertainty d: " << setw(12) << accv_uncPDF_dnm << std::endl;
    cout << endl;
  }
  
  char txtfname[100];
  sprintf(txtfname,"%s/gen.txt",outputDir.Data());
  ofstream txtfile;
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  if(charge== 0) txtfile << " W -> mu nu"  << endl;
  if(charge==-1) txtfile << " W- -> mu nu" << endl;
  if(charge== 1) txtfile << " W+ -> mu nu" << endl;
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
    txtfile << "            barrel: " << setw(12) << nSelBv[ifile] << " / " << nEvtsv[ifile] << " = " << accBv[ifile] << " +/- " << accErrBv[ifile] << endl;
    txtfile << "            endcap: " << setw(12) << nSelEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEv[ifile] << " +/- " << accErrEv[ifile] << endl;
    txtfile << "             total: " << setw(12) << nSelv[ifile]  << " / " << nEvtsv[ifile] << " = " << accv[ifile]  << " +/- " << accErrv[ifile] << endl;
    txtfile << "   qcd uncertainty: " << setw(12) << accv_uncQCD << std::endl;
    txtfile << " qcd uncertainty n: " << setw(12) << accv_uncQCD_num << std::endl;
    txtfile << " qcd uncertainty d: " << setw(12) << accv_uncQCD_dnm << std::endl;
    txtfile << "   pdf uncertainty: " << setw(12) << accv_uncPDF << std::endl;
    txtfile << " pdf uncertainty n: " << setw(12) << accv_uncPDF_num << std::endl;
    txtfile << " pdf uncertainty d: " << setw(12) << accv_uncPDF_dnm << std::endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccGenWm_Sys"); 
}
