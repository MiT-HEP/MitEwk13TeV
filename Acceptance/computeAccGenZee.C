//================================================================================================
// Not used for 13 TeV measurement.
//
// Compute Z->ee acceptance at generator level
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

void computeAccGenZee(const TString conf,             // input file
                      const TString outputDir,         // output directory
                       const TString outputName, // output filename
                      const bool doDressed=0
) {
  gBenchmark->Start("computeAccGenZee_Sys");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW   = 60;
  const Double_t MASS_HIGH  = 120;
  const Double_t PT_CUT     = 25;
  const Double_t ETA_CUT    = 2.4; 
  // const Double_t ETA_BARREL = 1.4442;
  // const Double_t ETA_ENDCAP = 1.566;
  const Double_t ETA_BARREL = 10.;
  const Double_t ETA_ENDCAP = 10.;

  const Int_t BOSON_ID  = 23;
  const Int_t LEPTON_ID = 11;
  
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
  vector<Double_t> nEvtsv, nSelv, nSelBBv, nSelBEv, nSelEEv;
  vector<Double_t> accv, accBBv, accBEv, accEEv;
  vector<Double_t> accErrv, accErrBBv, accErrBEv, accErrEEv;
  
  vector<Double_t> nEvtsv_QCD, nSelv_QCD;
  vector<Double_t> accv_QCD;
  vector<Double_t> accErrv_QCD;
  
  vector<Double_t> nEvtsv_PDF, nSelv_PDF;
  vector<Double_t> accv_PDF;
  vector<Double_t> accErrv_PDF;
    
  double accv_uncPDF = 0, accv_uncPDF_num = 0, accv_uncPDF_dnm = 0;
  double accv_uncQCD = 0, accv_uncQCD_num = 0, accv_uncQCD_dnm = 0;
    
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
    nSelBBv.push_back(0);
    nSelBEv.push_back(0);
    nSelEEv.push_back(0);
    
    for(int i=0;i<NQCD;++i) {nSelv_QCD.push_back(0);nEvtsv_QCD.push_back(0);}
    for(int i=0;i<NPDF;++i) {nSelv_PDF.push_back(0);nEvtsv_PDF.push_back(0);}
    
    //
    // loop over events
    //
    // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    for(UInt_t ientry=0; (uint)ientry<(eventTree->GetEntries())*0.1; ientry++) {
    // for(UInt_t ientry=0; ientry<100; ientry++) {
      if(ientry%1000000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
      genBr->GetEntry(ientry);
      genPartArr->Clear(); partBr->GetEntry(ientry);

      TLorentzVector *vec=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
      TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
      Int_t lepq1=-99;
      Int_t lepq2=-99;
      if (fabs(toolbox::flavor(genPartArr, BOSON_ID))!=LEPTON_ID) continue;
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

      if(vec->M()<MASS_LOW || vec->M()>MASS_HIGH) continue;
      Double_t weight=gen->weight;
      nEvtsv[ifile]+=weight; //std::cout << weight << std::endl;
      
      // -------------------------------------------------
      // I'm not sure which indexing is correct for Z's 
      // -------------------------------------------------
      
      // standard 13 tev ??
      nEvtsv_QCD[0]+=weight*gen->lheweight[0];
      nEvtsv_QCD[1]+=weight*gen->lheweight[1];
      nEvtsv_QCD[2]+=weight*gen->lheweight[2];
      nEvtsv_QCD[3]+=weight*gen->lheweight[3];
      nEvtsv_QCD[4]+=weight*gen->lheweight[5];
      nEvtsv_QCD[5]+=weight*gen->lheweight[7];
      for(int npdf=0; npdf<NPDF; npdf++) nEvtsv_PDF[npdf]+=weight*gen->lheweight[8+npdf];
            
      // // minlo
      // nEvtsv_QCD[0]+=weight*gen->lheweight[1];
      // nEvtsv_QCD[1]+=weight*gen->lheweight[2];
      // nEvtsv_QCD[2]+=weight*gen->lheweight[3];
      // nEvtsv_QCD[3]+=weight*gen->lheweight[4];
      // nEvtsv_QCD[4]+=weight*gen->lheweight[6];
      // nEvtsv_QCD[5]+=weight*gen->lheweight[8];
      // for(int npdf=0; npdf<NPDF; npdf++) nEvtsv_PDF[npdf]+=weight*gen->lheweight[9+npdf];

      if(lep1->Pt() < PT_CUT)         continue;
      if(lep2->Pt() < PT_CUT)         continue;
      if(fabs(lep1->Eta()) > ETA_CUT) continue;
      if(fabs(lep2->Eta()) > ETA_CUT) continue;
      if(fabs(lep1->Eta())>ETA_BARREL && fabs(lep1->Eta())<ETA_ENDCAP) continue;
      if(fabs(lep2->Eta())>ETA_BARREL && fabs(lep2->Eta())<ETA_ENDCAP) continue;

      TLorentzVector dilep=(*lep1)+(*lep2);

      if(dilep.M()<MASS_LOW || dilep.M()>MASS_HIGH) continue;
      Bool_t isB1 = (fabs(lep1->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      Bool_t isB2 = (fabs(lep2->Eta())<ETA_BARREL) ? kTRUE : kFALSE;
      
      nSelv[ifile]+=weight;
      if(isB1 && isB2)        { nSelBBv[ifile]+=weight; }
      else if(!isB1 && !isB2) { nSelEEv[ifile]+=weight; }
      else		      { nSelBEv[ifile]+=weight; } 
      // Get the values with the QCD and PDF weights:
      // QCD first
      nSelv_QCD[0]+=weight*gen->lheweight[0];
      nSelv_QCD[1]+=weight*gen->lheweight[1];
      nSelv_QCD[2]+=weight*gen->lheweight[2];
      nSelv_QCD[3]+=weight*gen->lheweight[3];
      nSelv_QCD[4]+=weight*gen->lheweight[5];
      nSelv_QCD[5]+=weight*gen->lheweight[7];
      for(int npdf=0; npdf<NPDF; npdf++) nSelv_PDF[npdf]+=weight*gen->lheweight[8+npdf];
      // cout << "weight = " << weight << endl;
      // cout << gen->lheweight[0]  << " " << gen->lheweight[1]  << " " << gen->lheweight[2]  << " " << gen->lheweight[3]  << " " << gen->lheweight[4]  << " " << gen->lheweight[5]  << " " << gen->lheweight[6]  << " " << gen->lheweight[7]  << " " << gen->lheweight[8]  << " " << endl;
      // nSelv_QCD[0]+=weight*gen->lheweight[1];
      // nSelv_QCD[1]+=weight*gen->lheweight[2];
      // nSelv_QCD[2]+=weight*gen->lheweight[3];
      // nSelv_QCD[3]+=weight*gen->lheweight[4];
      // nSelv_QCD[4]+=weight*gen->lheweight[6];
      // nSelv_QCD[5]+=weight*gen->lheweight[8];
      // for(int npdf=0; npdf<NPDF; npdf++)  nSelv_PDF[npdf]+=weight*gen->lheweight[9+npdf];
    }
    
    // compute acceptances
    accv.push_back(nSelv[ifile]/nEvtsv[ifile]);     accErrv.push_back(sqrt(accv[ifile]*(1.-accv[ifile])/nEvtsv[ifile]));
    accBBv.push_back(nSelBBv[ifile]/nEvtsv[ifile]); accErrBBv.push_back(sqrt(accBBv[ifile]*(1.-accBBv[ifile])/nEvtsv[ifile]));
    accBEv.push_back(nSelBEv[ifile]/nEvtsv[ifile]); accErrBEv.push_back(sqrt(accBEv[ifile]*(1.-accBEv[ifile])/nEvtsv[ifile]));
    accEEv.push_back(nSelEEv[ifile]/nEvtsv[ifile]); accErrEEv.push_back(sqrt(accEEv[ifile]*(1.-accEEv[ifile])/nEvtsv[ifile]));
    
    std::cout << "nselv " << nSelv[ifile] << "  nevtsv " << nEvtsv[ifile] << std::endl;
    
  // char txtfnamePDFs[500];
  // sprintf(txtfnamePDFs,"%s/pdf_vars.txt",outputDir.Data());
  // ofstream txtfile1;
  // txtfile1.open(txtfnamePDFs);
  // txtfile1 << accv[0] << std::endl;
    
    // double sumCheck=0;
    // double sqCheck = 0;
    // double maxPDF = 0;
    // for(int npdf=0; npdf<NPDF; npdf++){
        // accv_PDF.push_back(nSelv_PDF[npdf]/nEvtsv_PDF[npdf]);     
        // std::cout << "accv _pdf" << npdf << "  = " << accv_PDF[npdf] << std::endl;
        // accErrv_PDF.push_back(sqrt(accv_PDF[npdf]*(1.-accv_PDF[npdf])/nEvtsv_PDF[npdf]));
        // std::cout << "nselvpdf " << nSelv_PDF[npdf] << "  nevtsvpdf " << nEvtsv_PDF[npdf] << " ratio " << nSelv_PDF[npdf]/nEvtsv_PDF[npdf] << std::endl;
        // std::cout << "accv " << accv[ifile] << "  accvpdf " << accv_PDF[npdf] << std::endl;
        // std::cout << "diff " << accv_PDF[npdf]-accv[ifile] << "  pct diff " << 100*(accv_PDF[npdf]-accv[ifile])/accv[ifile] << std::endl;
        // // a = sum(((acc_i-accTot)/accTot)^2)
        // // unc = sqrt(a/NPDF)
        // double pctDiff = (accv_PDF[npdf]-accv[ifile])/accv[ifile];
        // accv_uncPDF+=pctDiff*pctDiff;
        // sumCheck+=fabs(pctDiff);
        // sqCheck+=accv_uncPDF;
        // (fabs(pctDiff)>maxPDF)?maxPDF=fabs(pctDiff):maxPDF=maxPDF;
        // accv_uncPDF_num+=(nSelv_PDF[npdf]-nSelv[ifile])*(nSelv_PDF[npdf]-nSelv[ifile])/(NPDF*nSelv[ifile]*nSelv[ifile]);
        // accv_uncPDF_dnm+=(nEvtsv_PDF[npdf]-nEvtsv[ifile])*(nEvtsv_PDF[npdf]-nEvtsv[ifile])/(NPDF*nEvtsv[ifile]*nEvtsv[ifile]);
        // txtfile1 << accv_PDF[npdf] << std::endl;
    // }
    // txtfile1.close();
  // sprintf(txtfnamePDFs,"%s/qcd_vars.txt",outputDir.Data());
  // ofstream txtfile2;
  // txtfile2.open(txtfnamePDFs);
  // txtfile2 << accv[0] << std::endl;
    
    // std::cout << "wum check " << sumCheck << std::endl;
    // std::cout << "sq check " << sqCheck << std::endl;
    // std::cout << "maxpdf  " << maxPDF << std::endl;
    // accv_uncPDF=sqrt(accv_uncPDF/NPDF);
    // accv_uncPDF_num=sqrt(accv_uncPDF_num);
    // accv_uncPDF_dnm=sqrt(accv_uncPDF_dnm);
    // for(int nqcd=0; nqcd<NQCD; nqcd++){
        // accv_QCD.push_back(nSelv_QCD[nqcd]/nEvtsv_QCD[nqcd]);
        // // std::cout << "nselvqcd " << nSelv_QCD[nqcd] << "  nevtsvqcd " << nEvtsv_QCD[nqcd] << " ratio " << nSelv_QCD[nqcd]/nEvtsv_QCD[nqcd] << std::endl;
        // // std::cout << "accv " << accv[ifile] << "  accvqcd " << accv_QCD[nqcd] << std::endl;
        // std::cout << "qcd a " << nqcd <<" " << accv_QCD[nqcd] << "  nom " << accv[ifile] << "  pct diff " << 100*fabs(accv_QCD[nqcd]-accv[ifile])/(accv[ifile]) << std::endl;
        // std::cout << "qcd n " << nqcd <<" " << nSelv_QCD[nqcd] << "  nom " << nSelv[ifile] << "  pct diff n " << 100*fabs(nSelv_QCD[nqcd]-nSelv[ifile])/nSelv[ifile] << std::endl;
        // std::cout << "qcd d " << nqcd <<" " << nEvtsv_QCD[nqcd] << "  nom " << nEvtsv[ifile] << "  pct diff d " << 100*fabs(nEvtsv_QCD[nqcd]-nEvtsv[ifile])/nEvtsv[ifile] << std::endl;
        // accErrv_QCD.push_back(sqrt(accv_QCD[nqcd]*(1.-accv_QCD[nqcd])/nEvtsv_QCD[nqcd]));
        // // a = abs((acc_i-accTot)/accTot)
        // if(fabs(accv_QCD[nqcd]-accv[ifile])/(accv[ifile]) > accv_uncQCD) accv_uncQCD = fabs(accv_QCD[nqcd]-accv[ifile])/(accv[ifile]);
        // if(fabs(nSelv_QCD[nqcd]-nSelv[ifile])/nSelv[ifile] > accv_uncQCD_num) accv_uncQCD_num = fabs(nSelv_QCD[nqcd]-nSelv[ifile])/nSelv[ifile];
        // if(fabs(nEvtsv_QCD[nqcd]-nEvtsv[ifile])/nEvtsv[ifile] > accv_uncQCD_dnm) accv_uncQCD_dnm = fabs(nEvtsv_QCD[nqcd]-nEvtsv[ifile])/nEvtsv[ifile];
        // txtfile2 << accv_QCD[nqcd] << std::endl;
    // }
    // txtfile2.close();
    
    
    delete infile;
    infile=0, eventTree=0;  
  } 
  
  delete gen;
  
  
  
  // Print full set for efficiency calculations
  char masterOutput[600];
  // just start printing....
  for(uint ifile = 0; ifile < fnamev.size(); ++ifile){// go through info per file
    sprintf(masterOutput,"%s/%s.txt",outputDir.Data(),outputName.Data());
    ofstream txtfile;
    txtfile.open(masterOutput);
    txtfile << "acc " << nSelv[ifile]/nEvtsv[ifile] << endl;
    
    for(int j = 0; j < NPDF; ++j){
      txtfile << "pdf" << j << " " << nSelv_PDF[j]/nEvtsv_PDF[j] << endl;
    }
    for(int j = 0; j < NQCD; ++j){
      txtfile << "qcd" << j << " " << nSelv_QCD[j]/nEvtsv_QCD[j] << endl;
    }
    txtfile.close();
  }
    
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
    cout << "     barrel-endcap: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    cout << "     endcap-endcap: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrEEv[ifile] << endl;
    cout << "             total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    // cout << "   qcd uncertainty: " << setw(12) << accv_uncQCD << std::endl;
    // cout << " qcd uncertainty n: " << setw(12) << accv_uncQCD_num << std::endl;
    // cout << " qcd uncertainty d: " << setw(12) << accv_uncQCD_dnm << std::endl;
    // cout << "   pdf uncertainty: " << setw(12) << accv_uncPDF << std::endl;
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
    txtfile << "     barrel-endcap: " << setw(12) << nSelBEv[ifile] << " / " << nEvtsv[ifile] << " = " << accBEv[ifile] << " +/- " << accErrBEv[ifile] << endl;
    txtfile << "     endcap-endcap: " << setw(12) << nSelEEv[ifile] << " / " << nEvtsv[ifile] << " = " << accEEv[ifile] << " +/- " << accErrEEv[ifile] << endl;
    txtfile << "             total: " << setw(12) << nSelv[ifile]   << " / " << nEvtsv[ifile] << " = " << accv[ifile]   << " +/- " << accErrv[ifile] << endl;
    // txtfile << "   qcd uncertainty: " << setw(12) << accv_uncQCD << std::endl;
    // txtfile << " qcd uncertainty n: " << setw(12) << accv_uncQCD_num << std::endl;
    // txtfile << " qcd uncertainty d: " << setw(12) << accv_uncQCD_dnm << std::endl;
    // txtfile << "   pdf uncertainty: " << setw(12) << accv_uncPDF << std::endl;
    // txtfile << " pdf uncertainty n: " << setw(12) << accv_uncPDF_num << std::endl;
    // txtfile << " pdf uncertainty d: " << setw(12) << accv_uncPDF_dnm << std::endl;
    txtfile << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("computeAccGenZee_Sys"); 
}
