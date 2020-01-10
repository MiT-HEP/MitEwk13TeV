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
  
  // vector<Double_t> nEvtsv_QCD, nSelv_QCD;
  vector<Double_t> accv_QCD;
  vector<Double_t> accErrv_QCD;
  
  // vector<Double_t> nEvtsv_PDF, nSelv_PDF;
  vector<Double_t> accv_PDF;
  vector<Double_t> accErrv_PDF;
  
  vector<vector<Double_t>> nEvtsv_QCD, nSelv_QCD;
  // vector<vector<Double_t>> accv_QCD;
  // vector<vector<Double_t>> accErrv_QCD;
  
  vector<vector<Double_t>> nEvtsv_PDF, nSelv_PDF;
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
    std::cout << "cross section is ... "  << xsecv[ifile] << std::endl;
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
    
    vector<Double_t> tempQCD_Selv, tempQCD_Evtsv;
    vector<Double_t> tempPDF_Selv, tempPDF_Evtsv;
    
    for(int i=0;i<NQCD;++i) {tempQCD_Selv.push_back(0);tempQCD_Evtsv.push_back(0);}
    for(int i=0;i<NPDF;++i) {tempPDF_Selv.push_back(0);tempPDF_Evtsv.push_back(0);}
    
    //
    // loop over events
    //    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    // for(UInt_t ientry=0; ientry<(uint)(0.10*eventTree->GetEntries()); ientry++) {
      if(ientry%100000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
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
      // std::cout << "---------------" << std::endl;
    //  std::cout << "evtnum " << ientry <<  " lep1 " << lep1->Pt() << "  lepq1 " << lepq1 << std::endl;
    //  std::cout << "evtnum " << ientry << " lep2 " << lep2->Pt() << "  lepq2 " << lepq2 << std::endl;
      
      if(charge!=0&&charge!=lepq1) lep1=lep2;
      // if(charge==0&&toolbox::flavor(genPartArr, BOSON_ID)*lepq1<0){
        // // std::cout << "lep1 ! " << BOSON_ID*glepq1 << std::endl;
        // genLep->SetPtEtaPhiM(glep1->Pt(),glep1->Eta(),glep1->Phi(),glep1->M());
      // }
      if(charge==0&&toolbox::flavor(genPartArr, BOSON_ID)*lepq2<0){
        // std::cout << "lep2 ! " << BOSON_ID*glepq2 << std::endl;
        //genLep->SetPtEtaPhiM(glep2->Pt(),glep2->Eta(),glep2->Phi(),glep2->M());
        lep1=lep2;
      }
      //if(charge==0&&charge==lepq1) lep1=lep2;
     
      
      TLorentzVector *gph=new TLorentzVector(0,0,0,0);
      if(doDressed){
          for(Int_t i=0; i<genPartArr->GetEntriesFast(); i++) {
            const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
            if(fabs(genloop->pdgId)!=22) continue;
            gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep1->Eta(),lep1->Phi())<0.1)
              {
                lep1->operator+=(*gph);
              } else if(toolbox::deltaR(gph->Eta(),gph->Phi(),lep2->Eta(),lep2->Phi())<0.1)
              {
                lep2->operator+=(*gph);
              }
          }
      }
      delete lep2;

      Double_t weight=gen->weight;
      nEvtsv[ifile]+=weight;
      // -------------------------------------------------
      // Clean this up? 
      // -------------------------------------------------
      // nEvtsv_QCD[ifile][0]+=weight*gen->lheweight[1];
      // nEvtsv_QCD[ifile][1]+=weight*gen->lheweight[2];
      // nEvtsv_QCD[ifile][2]+=weight*gen->lheweight[3];
      // nEvtsv_QCD[ifile][3]+=weight*gen->lheweight[4];
      // nEvtsv_QCD[ifile][4]+=weight*gen->lheweight[6];
      // nEvtsv_QCD[ifile][5]+=weight*gen->lheweight[8];
      // for(int npdf=0; npdf<NPDF; npdf++) nEvtsv_PDF[ifile][npdf]+=weight*gen->lheweight[9+npdf];
      
      tempQCD_Evtsv[0]+=weight*gen->lheweight[1];
      tempQCD_Evtsv[1]+=weight*gen->lheweight[2];
      tempQCD_Evtsv[2]+=weight*gen->lheweight[3];
      tempQCD_Evtsv[3]+=weight*gen->lheweight[4];
      tempQCD_Evtsv[4]+=weight*gen->lheweight[6];
      tempQCD_Evtsv[5]+=weight*gen->lheweight[8];
      
      for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Evtsv[npdf]+=weight*gen->lheweight[9+npdf];
    
    // // 2015
      // tempQCD_Evtsv[0]+=weight*gen->lheweight[1]/gen->lheweight[0];
      // tempQCD_Evtsv[1]+=weight*gen->lheweight[2]/gen->lheweight[0];
      // tempQCD_Evtsv[2]+=weight*gen->lheweight[3]/gen->lheweight[0];
      // tempQCD_Evtsv[3]+=weight*gen->lheweight[4]/gen->lheweight[0];
      // tempQCD_Evtsv[4]+=weight*gen->lheweight[6]/gen->lheweight[0];
      // tempQCD_Evtsv[5]+=weight*gen->lheweight[8]/gen->lheweight[0];
      
      // for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Evtsv[npdf]+=weight*gen->lheweight[9+npdf]/gen->lheweight[0];
    
      Bool_t isBarrel=kTRUE;
      if (lep1) {
        if (fabs(lep1->Eta())>ETA_BARREL && fabs(lep1->Eta())<ETA_ENDCAP) continue;
	
        if (lep1->Pt() < PT_CUT) continue;
        if (fabs(lep1->Eta()) > ETA_CUT) continue;
        isBarrel = (fabs(lep1->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      } else if (lep2) {
        if (fabs(lep2->Eta())>ETA_BARREL && fabs(lep2->Eta())<ETA_ENDCAP) continue;
std::cout << "eta 2 cut " << std::endl;
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
      
      tempQCD_Selv[0]+=weight*gen->lheweight[1];
      tempQCD_Selv[1]+=weight*gen->lheweight[2];
      tempQCD_Selv[2]+=weight*gen->lheweight[3];
      tempQCD_Selv[3]+=weight*gen->lheweight[4];
      tempQCD_Selv[4]+=weight*gen->lheweight[6];
      tempQCD_Selv[5]+=weight*gen->lheweight[8];
      for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Selv[npdf]+=weight*gen->lheweight[9+npdf];
      
      // tempQCD_Selv[0]+=weight*gen->lheweight[1]/gen->lheweight[0];
      // tempQCD_Selv[1]+=weight*gen->lheweight[2]/gen->lheweight[0];
      // tempQCD_Selv[2]+=weight*gen->lheweight[3]/gen->lheweight[0];
      // tempQCD_Selv[3]+=weight*gen->lheweight[4]/gen->lheweight[0];
      // tempQCD_Selv[4]+=weight*gen->lheweight[6]/gen->lheweight[0];
      // tempQCD_Selv[5]+=weight*gen->lheweight[8]/gen->lheweight[0];
      // for(int npdf=0; npdf<NPDF; npdf++) tempPDF_Selv[npdf]+=weight*gen->lheweight[9+npdf]/gen->lheweight[0];
    }
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);   accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBv.push_back(nSelBv[ifile]/nEvtsv[ifile]); accErrBv.push_back(sqrt(accBv[ifile]*(1.-accBv[ifile])/nEvtsv[ifile]));
    accEv.push_back(nSelEv[ifile]/nEvtsv[ifile]); accErrEv.push_back(sqrt(accEv[ifile]*(1.-accEv[ifile])/nEvtsv[ifile]));
    nSelv_PDF.push_back(tempPDF_Selv);
    nEvtsv_PDF.push_back(tempPDF_Evtsv);
    nSelv_QCD.push_back(tempQCD_Selv);
    nEvtsv_QCD.push_back(tempQCD_Evtsv);
    
    
   std::cout << "nselv " << nSelv[ifile] << "  nevtsv " << nEvtsv[ifile] << std::endl;
    
    

    
    
    delete infile;
    infile=0, eventTree=0;  
  }
  
  // vector<Double_t> accv_PDF
  double accTot=0;
  double accNum=0, accDnm=0;
  std::cout << "here" << std::endl;

  
  
  for(int i=0;i<NPDF;++i) {accv_PDF.push_back(0);}
  for(int i=0;i<NQCD;++i) {accv_QCD.push_back(0);}
  
  for(int ifile=0;ifile<NFILES;ifile++){
    std::cout << "in loop" << std::endl;
    accNum += accv[ifile]*xsecv[ifile];
    accDnm += xsecv[ifile];
    for(int ipdf=0; ipdf<NPDF; ipdf++){
      
      accv_PDF[ipdf]+=xsecv[ifile]*nSelv_PDF[ifile][ipdf]/nEvtsv_PDF[ifile][ipdf];
    }
  }
  accTot=accNum/accDnm;
  
    
  
  char txtfnamePDFs[100];
  sprintf(txtfnamePDFs,"%s/pdf_vars.txt",outputDir.Data());
  ofstream txtfile1;
  txtfile1.open(txtfnamePDFs);
  txtfile1 << accTot << std::endl;
  for(int ipdf=0; ipdf < NPDF; ipdf++){
    accv_PDF[ipdf]=accv_PDF[ipdf]/accDnm;
    std::cout << "accv " << accTot << "  accvpdf " << accv_PDF[ipdf] << std::endl;
    std::cout << "diff " << accv_PDF[ipdf]-accTot << "  pct diff " << 100*(accv_PDF[ipdf]-accTot)/accTot    << std::endl;
    accv_uncPDF+=(accv_PDF[ipdf]-accTot)*(accv_PDF[ipdf]-accTot)/(NPDF*accTot*accTot);
    txtfile1 << accv_PDF[ipdf] << std::endl;
  }
  txtfile1.close();
  
  accv_uncPDF=sqrt(accv_uncPDF);

  
  sprintf(txtfnamePDFs,"%s/qcd_vars.txt",outputDir.Data());
  ofstream txtfile2;
  txtfile2.open(txtfnamePDFs);
  txtfile2 << accTot << std::endl;
  for(int ifile=0;ifile<NFILES;ifile++){
    std::cout << "file -------------- " << ifile << std::endl;
    accv_QCD.push_back(0);
    for(int iqcd=0; iqcd<NQCD; iqcd++){
      std::cout << "iqcd " << iqcd << "  num " << nSelv_QCD[ifile][iqcd] << "  denom " << nEvtsv_QCD[ifile][iqcd] << " acc " << nSelv_QCD[ifile][iqcd]/nEvtsv_QCD[ifile][iqcd]  << std::endl;
      accv_QCD[iqcd]+=xsecv[ifile]*nSelv_QCD[ifile][iqcd]/(nEvtsv_QCD[ifile][iqcd]*accDnm);
      // accv_uncQCD+=(accv_QCD[iqcd]-accv[ifile])*(accv_QCD[iqcd]-accv[ifile])/(NQCD*accv[ifile]*accv[ifile]);
      
    }
  }
  for(int iqcd=0; iqcd<NQCD; iqcd++){
    std::cout << fabs(accv_QCD[iqcd]-accTot)/(accTot) << std::endl;
    if(fabs(accv_QCD[iqcd]-accTot)/(accTot) > accv_uncQCD) accv_uncQCD = fabs(accv_QCD[iqcd]-accTot)/(accTot);
    txtfile2 << accv_QCD[iqcd] << std::endl;
  }
  txtfile2.close();

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
    // cout << " qcd uncertainty n: " << setw(12) << accv_uncQCD_num << std::endl;
    // cout << " qcd uncertainty d: " << setw(12) << accv_uncQCD_dnm << std::endl;
    cout << "   pdf uncertainty: " << setw(12) << accv_uncPDF << std::endl;
    // cout << " pdf uncertainty n: " << setw(12) << accv_uncPDF_num << std::endl;
    // cout << " pdf uncertainty d: " << setw(12) << accv_uncPDF_dnm << std::endl;
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
    // txtfile << " qcd uncertainty n: " << setw(12) << accv_uncQCD_num << std::endl;
    // txtfile << " qcd uncertainty d: " << setw(12) << accv_uncQCD_dnm << std::endl;
    txtfile << "   pdf uncertainty: " << setw(12) << accv_uncPDF << std::endl;
    // txtfile << " pdf uncertainty n: " << setw(12) << accv_uncPDF_num << std::endl;
    // txtfile << " pdf uncertainty d: " << setw(12) << accv_uncPDF_dnm << std::endl;
    txtfile << endl;
  }
  txtfile.close();  
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccGenWm_Sys"); 
}
