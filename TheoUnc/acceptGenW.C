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
#include "LHAPDF/LHAPDF.h"

using namespace std;

#endif

enum PROC {wme=0, wpe, wmm, wpm};

Bool_t isProc(Int_t proc, Double_t genV_id, Double_t genV_m, Double_t genL1_id, Double_t genL2_id);

Bool_t acceptEndcap(Int_t proc, Double_t genV_id, Double_t genV_m, Double_t genL1_id, Double_t genL2_id, Double_t genL1_pt, Double_t genL1_eta, Double_t genL2_pt, Double_t genL2_eta);

Bool_t acceptBarrel(Int_t proc, Double_t genV_id, Double_t genV_m, Double_t genL1_id, Double_t genL2_id, Double_t genL1_pt, Double_t genL1_eta, Double_t genL2_pt, Double_t genL2_eta);

//void acceptGenW(TString input="wellnu_ct10_wp_v2.root",
void acceptGenW(TString input="/afs/cern.ch/work/j/jlawhorn/theo-unc-06-10/WmJetsToLNu.root",
		TString outputDir="idfk/",
		//TString pdfName="NNPDF30_nlo_nf_5_pdfas",
		//TString pdfName="CT14nlo",
		TString pdfName="MMHT2014nlo68cl",
		Int_t setMin=0, 
		Int_t setMax=0,
		Int_t proc=2) {

  TString procName[4]={"wme", "wpe", "wmm", "wpm"};
  char output[150];
  sprintf(output, "%s%s_%s.root", outputDir.Data(), procName[proc].Data(), pdfName.Data());
  cout << "proc: " << procName[proc] << endl;
  TChain chain("Events");
  chain.Add(input);
  
  //PDF info
  Double_t id_1,      id_2,       x_1,        x_2;
  Double_t xPDF_1,    xPDF_2,     scalePDF,   weight;
  //Generator level V+l info
  Double_t genV_id,   genL1_id,   genL2_id;
  Double_t genV_pt,   genV_eta,   genV_phi,   genV_m;
  Double_t genVf_pt,  genVf_eta,  genVf_phi,  genVf_m;
  Double_t genL1_pt,  genL1_eta,  genL1_phi,  genL1_m;
  Double_t genL2_pt,  genL2_eta,  genL2_phi,  genL2_m;
  Double_t genL1f_pt, genL1f_eta, genL1f_phi, genL1f_m;
  Double_t genL2f_pt, genL2f_eta, genL2f_phi, genL2f_m;

  chain.SetBranchAddress("id_1",       &id_1);
  chain.SetBranchAddress("id_2",       &id_2);
  chain.SetBranchAddress("x_1",        &x_1);
  chain.SetBranchAddress("x_2",        &x_2);
  chain.SetBranchAddress("xPDF_1",     &xPDF_1);
  chain.SetBranchAddress("xPDF_2",     &xPDF_2);
  chain.SetBranchAddress("scalePDF",   &scalePDF);
  chain.SetBranchAddress("weight",     &weight);
  chain.SetBranchAddress("genV_pt",    &genV_pt);
  chain.SetBranchAddress("genV_eta",   &genV_eta);
  chain.SetBranchAddress("genV_phi",   &genV_phi);
  chain.SetBranchAddress("genV_m",     &genV_m);
  chain.SetBranchAddress("genV_id",    &genV_id);
  chain.SetBranchAddress("genVf_pt",   &genVf_pt);
  chain.SetBranchAddress("genVf_eta",  &genVf_eta);
  chain.SetBranchAddress("genVf_phi",  &genVf_phi);
  chain.SetBranchAddress("genVf_m",    &genVf_m);
  chain.SetBranchAddress("genL1_pt",   &genL1_pt);
  chain.SetBranchAddress("genL1_eta",  &genL1_eta);
  chain.SetBranchAddress("genL1_phi",  &genL1_phi);
  chain.SetBranchAddress("genL1_m",    &genL1_m);
  chain.SetBranchAddress("genL1_id",   &genL1_id);
  chain.SetBranchAddress("genL2_pt",   &genL2_pt);
  chain.SetBranchAddress("genL2_eta",  &genL2_eta);
  chain.SetBranchAddress("genL2_phi",  &genL2_phi);
  chain.SetBranchAddress("genL2_m",    &genL2_m);
  chain.SetBranchAddress("genL2_id",   &genL2_id);
  chain.SetBranchAddress("genL1f_pt",  &genL1f_pt);
  chain.SetBranchAddress("genL1f_eta", &genL1f_eta);
  chain.SetBranchAddress("genL1f_phi", &genL1f_phi);
  chain.SetBranchAddress("genL1f_m",   &genL1f_m);
  chain.SetBranchAddress("genL2f_pt",  &genL2f_pt);
  chain.SetBranchAddress("genL2f_eta", &genL2f_eta);
  chain.SetBranchAddress("genL2f_phi", &genL2f_phi);
  chain.SetBranchAddress("genL2f_m",   &genL2f_m);

  TFile *outFile = new TFile(output, "recreate");
  LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(292200);
  //LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(25100);
  //LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(13100);
  //LHAPDF::PDF* nomPdf = LHAPDF::mkPDF(11000);

  for (Int_t iPdfSet=setMin; iPdfSet<setMax+1; iPdfSet++) {
    
    LHAPDF::PDF* testPdf = LHAPDF::mkPDF(pdfName.Data(),iPdfSet);
    
    char histname[100];
    sprintf(histname,"dEta_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dEta = new TH1D(histname, "", 20, -5, 5); dEta->Sumw2();
    sprintf(histname,"dPt_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dPt = new TH1D(histname, "", 25, 25, 100); dPt->Sumw2();
    
    sprintf(histname,"dPreB_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dPreB = new TH1D(histname, "", 1, 0, 2); dPreB->Sumw2();
    sprintf(histname,"dPreE_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dPreE = new TH1D(histname, "", 1, 0, 2); dPreE->Sumw2();
    
    sprintf(histname,"dPostB_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dPostB = new TH1D(histname, "", 1, 0, 2); dPostB->Sumw2();
    sprintf(histname,"dPostE_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dPostE = new TH1D(histname, "", 1, 0, 2); dPostE->Sumw2();
    
    sprintf(histname,"dTot_%s_%i",pdfName.Data(),iPdfSet);
    TH1D* dTot = new TH1D(histname, "", 1, 0, 2); dTot->Sumw2();
    
    for (Int_t i=0; i<chain.GetEntries(); i++) {
      chain.GetEntry(i);
      //cout  << genV_id << ", " << genL1_id << ", " << genL2_id << ", " << genL1f_pt << ", " << genL1f_eta << endl;
      if (!isProc(proc, genV_id, genV_m, genL1_id, genL2_id)) continue;
      if ((proc==0||proc==2) && genL2f_pt==0 && genL2f_eta==0) continue;
      if ((proc==1||proc==3) && genL1f_pt==0 && genL1f_eta==0) continue;

      Int_t fid_1 = ( id_1==0 ? 21 : id_1 );
      Int_t fid_2 = ( id_2==0 ? 21 : id_2 );
      //cout << genV_id << ", " << genL1f_pt << ", " << genL2f_pt << endl;
      Double_t newWeight = weight*testPdf->xfxQ(fid_1, x_1, scalePDF)*testPdf->xfxQ(fid_2, x_2, scalePDF)/(nomPdf->xfxQ(fid_1, x_1, scalePDF)*nomPdf->xfxQ(fid_2, x_2, scalePDF));
      if (newWeight < 3e-10) continue;

      dTot->Fill(1.0,newWeight);
      
      Bool_t passBar =acceptBarrel(proc, genV_id, genV_m, genL1_id, genL2_id, genL1_pt, genL1_eta, genL2_pt, genL2_eta);
      Bool_t passBarF=acceptBarrel(proc, genV_id, genV_m, genL1_id, genL2_id, genL1f_pt, genL1f_eta, genL2f_pt, genL2f_eta);
      Bool_t passEnd =acceptEndcap(proc, genV_id, genV_m, genL1_id, genL2_id, genL1_pt, genL1_eta, genL2_pt, genL2_eta);
      Bool_t passEndF=acceptEndcap(proc, genV_id, genV_m, genL1_id, genL2_id, genL1f_pt, genL1f_eta, genL2f_pt, genL2f_eta);
      
      if      (passBar) dPreB->Fill(1.0, newWeight);
      else if (passEnd) dPreE->Fill(1.0, newWeight);
      
      if      (passBarF) dPostB->Fill(1.0, newWeight);
      else if (passEndF) dPostE->Fill(1.0, newWeight);
      
      if (passBarF || passEndF) {
	if (proc==wme || proc==wmm) {
	  dEta->Fill(genL2f_eta, newWeight);  
	  dPt->Fill(genL2f_pt, newWeight);  
	}
	else {
	  dEta->Fill(genL1f_eta, newWeight);  
	  dPt->Fill(genL1f_pt, newWeight);        
	}
      }
    }

    Double_t acc=(dPreB->Integral()+dPreE->Integral())/(dTot->Integral());
    //cout << dPreB->Integral()/dTot->Integral() << endl;
    //cout << dPreE->Integral()/dTot->Integral() << endl;
    cout << "Pre-FSR acceptance: " << acc << " +/- " << sqrt(acc*(1-acc)/dTot->GetEntries()) << endl;
    acc=(dPostB->Integral()+dPostE->Integral())/(dTot->Integral());
    //cout << dPostB->Integral()/dTot->Integral() << endl;
    //cout << dPostE->Integral()/dTot->Integral() << endl;
    cout << "Post-FSR acceptance: " << acc << " +/- " << sqrt(acc*(1-acc)/dTot->GetEntries()) << endl;

    outFile->Write();

  }

  outFile->Close();
}

Bool_t isProc(Int_t proc, Double_t genV_id, Double_t genV_m, Double_t genL1_id, Double_t genL2_id) {
  if (proc==wme) {
    if (genV_id==-24 && genL2_id==11) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wpe) {
    if (genV_id==24 && genL1_id==-11) return kTRUE;
    else return kFALSE;
  }
  if (proc==wmm) {
    if (genV_id==-24 && genL2_id==13) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wpm) {
    if (genV_id==24 && genL1_id==-13) return kTRUE;
    else return kFALSE;
  }
  return kFALSE;
}

Bool_t acceptEndcap(Int_t proc, Double_t genV_id, Double_t genV_m, Double_t genL1_id, Double_t genL2_id, Double_t genL1_pt, Double_t genL1_eta, Double_t genL2_pt, Double_t genL2_eta) {
  
  if (proc==wme) {
    if (genL2_pt>25 && fabs(genL2_eta)>1.566 && fabs(genL2_eta)<2.5) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wpe) {
    if (genL1_pt>25 && fabs(genL1_eta)>1.566 && fabs(genL1_eta)<2.5) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wmm) {
    if (genL2_pt>25 && fabs(genL2_eta)>1.2 && fabs(genL2_eta)<2.4) return kTRUE;
    //if (genL2_pt>25 && fabs(genL2_eta)>1.2 && fabs(genL2_eta)<2.1) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wpm) {
    if (genL1_pt>25 && fabs(genL1_eta)>1.2 && fabs(genL1_eta)<2.4) return kTRUE;
    //if (genL1_pt>25 && fabs(genL1_eta)>1.2 && fabs(genL1_eta)<2.1) return kTRUE;
    else return kFALSE;
  }
  return kFALSE;
}

Bool_t acceptBarrel(Int_t proc, Double_t genV_id, Double_t genV_m, Double_t genL1_id, Double_t genL2_id, Double_t genL1_pt, Double_t genL1_eta, Double_t genL2_pt, Double_t genL2_eta) {
  
  if (proc==wme) {
    if (genL2_pt>25 && fabs(genL2_eta)<1.4442) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wpe) {
    if (genL1_pt>25 && fabs(genL1_eta)<1.4442) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wmm) {
    if (genL2_pt>25 && fabs(genL2_eta)<1.2) return kTRUE;
    else return kFALSE;
  }
  else if (proc==wpm) {
    if (genL1_pt>25 && fabs(genL1_eta)<1.2) return kTRUE;
    else return kFALSE;
  }
  return kFALSE;
}
