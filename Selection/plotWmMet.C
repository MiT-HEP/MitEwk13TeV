//================================================================================================
//
// Make plots of various distributions after Wmunu selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include "TLorentzVector.h"               // 4-vector class

#include "ConfParse.hh"                   // input conf file parser
#include "../Utils/CSample.hh"            // helper class to handle samples
#include "../Utils/MyTools.hh"            // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hMC, const TString name);

// make webpage
void makeHTML(const TString outDir);

// make LaTex table of bkg. to sig. ratios
void printRatiosLatex(vector<double> ratios,ostream& os);

//=== MAIN MACRO ================================================================================================= 

void plotWmMet(const TString  conf,      // input file
               const TString  inputDir,  // input directory
	       const TString  outputDir, // output directory
	       const Double_t lumi       // integrated luminosity (/fb)
) {
  gBenchmark->Start("plotWmMet");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  const TString format("png");
  
  const Double_t PT_CUT  = 25;
  const Double_t ETA_CUT = 2.4;
  
  const Double_t ETA_BARREL = 1.2;
  const Double_t ETA_ENDCAP = 1.2;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir + TString("/plots");

  
  //
  // Create histograms
  //
  enum { eInc, ePos, eNeg };
  vector<TH1D*> hMetv[3];
  vector<TH1D*> hMet2v[3];
   

  vector<TH1D*> hTrkMetv[3];
  vector<TH1D*> hTrkMet2v[3];

  char hname[100];
  for(UInt_t ich=0; ich<3; ich++) {
    for(UInt_t isam=0; isam<samplev.size(); isam++) {
      sprintf(hname,"hMet_ch%i_%i",ich,isam);      hMetv[ich].push_back(new TH1D(hname,"",50,0,100));        hMetv[ich][isam]->Sumw2();
      sprintf(hname,"hMet2_ch%i_%i",ich,isam);     hMet2v[ich].push_back(new TH1D(hname,"",50,0,100));       hMet2v[ich][isam]->Sumw2();

      sprintf(hname,"hTrkMet_ch%i_%i",ich,isam);      hTrkMetv[ich].push_back(new TH1D(hname,"",50,0,100));        hTrkMetv[ich][isam]->Sumw2();
      sprintf(hname,"hTrkMet2_ch%i_%i",ich,isam);     hTrkMet2v[ich].push_back(new TH1D(hname,"",50,0,100));       hTrkMet2v[ich][isam]->Sumw2();
    }
  }
  
  //
  // Declare output ntuple variables
  //  
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t genVPt, genVPhi, genVy, genVMass;
  Float_t genLepPt, genLepPhi;
  Float_t scale1fb;
  Float_t met, metPhi, tkMet, tkMetPhi, sumEt, mt, u1, u2;
  Int_t   q;
  TLorentzVector *lep=0;
  Double_t weightPDF;
  ///// muon specific /////
  Float_t trkIso, emIso, hadIso;
  Float_t pfChIso, pfGamIso, pfNeuIso, pfCombIso;
  Float_t d0, dz;
  Float_t muNchi2;
  UInt_t nPixHits, nTkLayers, nValidHits, nMatch, typeBits;
  
  TFile *infile=0;
  TTree *intree=0;
        
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(isam==0 && !hasData) continue;    
    
    // Read input file and get the TTrees
    TString infilename = inputDir + TString("/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << infilename << "..." << endl;
    infile = new TFile(infilename);	    assert(infile);
    intree = (TTree*)infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("runNum",     &runNum);       // event run number
    intree->SetBranchAddress("lumiSec",    &lumiSec);      // event lumi section
    intree->SetBranchAddress("evtNum",     &evtNum);       // event number
    intree->SetBranchAddress("npv",        &npv);          // number of primary vertices
    intree->SetBranchAddress("npu",        &npu);          // number of in-time PU events (MC)
    intree->SetBranchAddress("genVPt",     &genVPt);       // GEN boson pT
    intree->SetBranchAddress("genVPhi",    &genVPhi);      // GEN boson phi
    intree->SetBranchAddress("genVy",      &genVy);        // GEN boson rapidity
    intree->SetBranchAddress("genVMass",   &genVMass);     // GEN boson mass    
    intree->SetBranchAddress("genLepPt",   &genLepPt);     // GEN lepton pT
    intree->SetBranchAddress("genLepPhi",  &genLepPhi);    // GEN lepton phi
    intree->SetBranchAddress("scale1fb",   &scale1fb);     // event weight per 1/fb (MC)
    intree->SetBranchAddress("met",        &met);          // MET
    intree->SetBranchAddress("metPhi",     &metPhi);       // phi(MET)
    intree->SetBranchAddress("tkMet",      &tkMet);        // track MET
    intree->SetBranchAddress("tkMetPhi",   &tkMetPhi);     // phi(track MET)
    intree->SetBranchAddress("sumEt",      &sumEt);        // Sum ET
    intree->SetBranchAddress("mt",         &mt);           // transverse mass
    intree->SetBranchAddress("u1",         &u1);           // parallel component of recoil
    intree->SetBranchAddress("u2",         &u2);           // perpendicular component of recoil
    intree->SetBranchAddress("q",          &q);            // lepton charge
    intree->SetBranchAddress("lep",        &lep);          // lepton 4-vector
    intree->SetBranchAddress("weightPDF",  &weightPDF);    // PDF scale factor
    ///// muon specific /////
    intree->SetBranchAddress("trkIso",     &trkIso);       // track isolation of lepton
    intree->SetBranchAddress("emIso",      &emIso);        // ECAL isolation of lepton
    intree->SetBranchAddress("hadIso",     &hadIso);       // HCAL isolation of lepton
    intree->SetBranchAddress("pfChIso",    &pfChIso);      // PF charged hadron isolation of lepton
    intree->SetBranchAddress("pfGamIso",   &pfGamIso);     // PF photon isolation of lepton
    intree->SetBranchAddress("pfNeuIso",   &pfNeuIso);     // PF neutral hadron isolation of lepton
    intree->SetBranchAddress("pfCombIso",  &pfCombIso);    // PF combined isolation of lepton
    intree->SetBranchAddress("d0",         &d0);           // transverse impact parameter of lepton
    intree->SetBranchAddress("dz",         &dz);           // longitudinal impact parameter of lepton
    intree->SetBranchAddress("muNchi2",    &muNchi2);      // muon fit normalized chi^2 of lepton
    intree->SetBranchAddress("nPixHits",   &nPixHits);     // number of pixel hits of muon
    intree->SetBranchAddress("nTkLayers",  &nTkLayers);    // number of tracker layers of muon
    intree->SetBranchAddress("nMatch",     &nMatch);       // number of matched segments of muon  
    intree->SetBranchAddress("nValidHits", &nValidHits);   // number of valid muon hits of muon 
    intree->SetBranchAddress("typeBits",   &typeBits);     // number of valid muon hits of muon 
    
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      if(lep->Pt()        < PT_CUT)  continue;	
      if(fabs(lep->Eta()) > ETA_CUT) continue;

      Double_t weight = 1;
      if(isam!=0) {
        weight *= scale1fb*lumi;
      }
      
      for(UInt_t ich=0; ich<3; ich++) {
        if(ich==eInc || (ich==ePos && q>0) || (ich==eNeg && q<0)) {
	  hMetv[ich][isam]   ->Fill(met,       weight);
	  hMet2v[ich][isam]  ->Fill(met,       weight);

	  hTrkMetv[ich][isam]   ->Fill(tkMet,       weight);
	  hTrkMet2v[ich][isam]  ->Fill(tkMet,       weight);
	}
      }
    }
    delete infile;
    infile=0, intree=0;    
  }  

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  // string buffers
  char pname[100];       // plot name
  char xlabel[100];      // x-axis label
  char ylabel[100];      // y-axis label
  char lumitext[100];    // lumi label
  char sigtext[100];     // label for signal in legend

  // label for lumi
  if(lumi<0.1)    sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi*1000.);
  else if(lumi<1) sprintf(lumitext,"%.0f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi*1000.);
  else            sprintf(lumitext,"%.2f fb^{-1}  at  #sqrt{s} = 13 TeV",lumi);

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->SetPad(0,0.3,1.0,1.0);
  c->SetTopMargin(0.1);
  c->SetLeftMargin(0.18);  
  c->SetRightMargin(0.07);  
  c->SetTickx(1);
  c->SetTicky(1);

  for(UInt_t ich=0; ich<3; ich++) {
 
    sprintf(xlabel,"#slash{E}_{T} [GeV]");
    if(ich==eInc) sprintf(sigtext,"W#rightarrow#mu#nu");
    else if(ich==ePos) sprintf(sigtext,"W^{+}#rightarrow#mu^{+}#nu");
    else if(ich==eNeg) sprintf(sigtext,"W^{-}#rightarrow#mu^{-}#nu");

    //
    // PF MET
    //    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetv[ich][0]->GetBinWidth(1));
    sprintf(pname,"pfmet_ch%i",ich);
    CPlot plotMet_pf(pname,"",xlabel,ylabel);
//    plotMet.AddHist1D(hMetv[ich][0],samplev[0]->label,"E"); // No data to plot

    // Get yields
    double ndy_pf=0, nwtaunu_pf=0, ndiboson_pf=0, nttbar_pf=0, nsig_pf=0;

    Bool_t ewkDrawn = kFALSE;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      if(samplev[isam]->label=="Z#rightarrow#tau^{+}#tau^{-}") ndy_pf += hMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="W#rightarrow#tau#nu") nwtaunu_pf += hMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="WW" || samplev[isam]->label=="WZ" || samplev[isam]->label=="ZZ") ndiboson_pf += hMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="Top") nttbar_pf += hMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="W#rightarrow#mu#nu") nsig_pf += hMetv[ich][isam]->Integral();

      if(samplev[isam]->label=="W#rightarrow#mu#nu") {
        plotMet_pf.AddToStack(hMetv[ich][isam],sigtext,samplev[isam]->color,samplev[isam]->linecol);
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && !ewkDrawn) {
        plotMet_pf.AddToStack(hMetv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
        ewkDrawn=kTRUE;
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && ewkDrawn) {
        // only put "EWK" in legend once
        plotMet_pf.AddToStack(hMetv[ich][isam],samplev[isam]->color,samplev[isam]->linecol);
      } else {
        plotMet_pf.AddToStack(hMetv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
      }
    }

    // Save ratios
    vector<double> bkgyields_pf;
    vector<double> ratios_pf;

    bkgyields_pf.push_back(ndy_pf);
    bkgyields_pf.push_back(nwtaunu_pf);
    bkgyields_pf.push_back(ndiboson_pf);
    bkgyields_pf.push_back(nttbar_pf);

    cout << "W->munu pf met --- signal yield is " << nsig_pf << endl;
    for(UInt_t k=0;k<bkgyields_pf.size();k++) ratios_pf.push_back(bkgyields_pf[k]/nsig_pf*100);

    // To make a LaTex table
    if(ich==eInc) {
      ofstream latexfile_pf;
      char latexfname_pf[100];    
      sprintf(latexfname_pf,"%s/pfmetlatex.txt",outputDir.Data());
      latexfile_pf.open(latexfname_pf);
      assert(latexfile_pf.is_open());
      printRatiosLatex(ratios_pf,latexfile_pf);
    }

//    plotMet_pf.SetLegend(0.68,0.57,0.84,0.70);
    plotMet_pf.SetLegend(0.7,0.4,0.9,0.7); 
    plotMet_pf.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMet_pf.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotMet_pf.Draw(c,kTRUE,format,1);

    sprintf(ylabel,"Events / %.1f GeV/c",hMet2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"pfmetlog_ch%i",ich);
    CPlot plotMet2_pf(pname,"",xlabel,ylabel);
//    plotMet2.AddHist1D(hMet2v[ich][0],samplev[0]->label,"E"); // No data to plot
    ewkDrawn = kFALSE;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      if(samplev[isam]->label=="W#rightarrow#mu#nu") {
        plotMet2_pf.AddToStack(hMet2v[ich][isam],sigtext,samplev[isam]->color,samplev[isam]->linecol);
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && !ewkDrawn) {
        plotMet2_pf.AddToStack(hMet2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
        ewkDrawn=kTRUE;
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && ewkDrawn) {
        // only put "EWK" in legend once
        plotMet2_pf.AddToStack(hMet2v[ich][isam],samplev[isam]->color,samplev[isam]->linecol);
      } else {
        plotMet2_pf.AddToStack(hMet2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
      }
    }

//    plotMet2_pf.SetLegend(0.68,0.57,0.84,0.70);
    plotMet2_pf.SetLegend(0.22,0.59,0.42,0.87); 
    plotMet2_pf.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMet2_pf.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotMet2_pf.SetYRange(0.1,1e4);
    plotMet2_pf.SetLogy();
    plotMet2_pf.Draw(c,kTRUE,format,1);
    
    //
    // Track MET
    //    
    sprintf(ylabel,"Events / %.1f GeV/c",hMetv[ich][0]->GetBinWidth(1));
    sprintf(pname,"trkmet_ch%i",ich);
    CPlot plotMet_trk(pname,"",xlabel,ylabel);
//    plotMet.AddHist1D(hMetv[ich][0],samplev[0]->label,"E");

    // Get yields
    double ndy_trk=0, nwtaunu_trk=0, ndiboson_trk=0, nttbar_trk=0, nsig_trk=0;

    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      if(samplev[isam]->label=="Z#rightarrow#tau^{+}#tau^{-}") ndy_trk += hTrkMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="W#rightarrow#tau#nu") nwtaunu_trk += hTrkMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="WW" || samplev[isam]->label=="WZ" || samplev[isam]->label=="ZZ") ndiboson_trk += hTrkMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="Top") nttbar_trk += hTrkMetv[ich][isam]->Integral();
      else if(samplev[isam]->label=="W#rightarrow#mu#nu") nsig_trk += hTrkMetv[ich][isam]->Integral();

      if(samplev[isam]->label=="W#rightarrow#mu#nu") {
        plotMet_trk.AddToStack(hTrkMetv[ich][isam],sigtext,samplev[isam]->color,samplev[isam]->linecol);
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && !ewkDrawn) {
        plotMet_trk.AddToStack(hTrkMetv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
        ewkDrawn=kTRUE;
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && ewkDrawn) {
        // only put "EWK" in legend once
        plotMet_trk.AddToStack(hTrkMetv[ich][isam],samplev[isam]->color,samplev[isam]->linecol);
      } else {
        plotMet_trk.AddToStack(hTrkMetv[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
      }
    }

    // Save ratios
    vector<double> bkgyields_trk;
    vector<double> ratios_trk;

    bkgyields_trk.push_back(ndy_trk);
    bkgyields_trk.push_back(nwtaunu_trk);
    bkgyields_trk.push_back(ndiboson_trk);
    bkgyields_trk.push_back(nttbar_trk);

    cout << "W->munu trk met --- signal yield is " << nsig_trk << endl;
    for(UInt_t k=0;k<bkgyields_trk.size();k++) ratios_trk.push_back(bkgyields_trk[k]/nsig_trk*100);

    // To make a LaTex table
    if(ich==eInc) {
      ofstream latexfile_trk;
      char latexfname_trk[100];    
      sprintf(latexfname_trk,"%s/trkmetlatex.txt",outputDir.Data());
      latexfile_trk.open(latexfname_trk);
      assert(latexfile_trk.is_open());
      printRatiosLatex(ratios_trk,latexfile_trk);
    }

//    plotMet_trk.SetLegend(0.68,0.57,0.84,0.70);
    plotMet_trk.SetLegend(0.7,0.4,0.9,0.7); 
    plotMet_trk.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMet_trk.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotMet_trk.Draw(c,kTRUE,format,1);

    sprintf(ylabel,"Events / %.1f GeV/c",hTrkMet2v[ich][0]->GetBinWidth(1));
    sprintf(pname,"trkmetlog_ch%i",ich);
    CPlot plotMet2_trk(pname,"",xlabel,ylabel);
//    plotMet2.AddHist1D(hMet2v[ich][0],samplev[0]->label,"E");
    ewkDrawn = kFALSE;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      if(samplev[isam]->label=="W#rightarrowe#nu") {
        plotMet2_trk.AddToStack(hTrkMet2v[ich][isam],sigtext,samplev[isam]->color,samplev[isam]->linecol);
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && !ewkDrawn) {
        plotMet2_trk.AddToStack(hTrkMet2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
        ewkDrawn=kTRUE;
      } else if(samplev[isam]->label=="EWK+t#bar{t}" && ewkDrawn) {
        // only put "EWK" in legend once
        plotMet2_trk.AddToStack(hTrkMet2v[ich][isam],samplev[isam]->color,samplev[isam]->linecol);
      } else {
        plotMet2_trk.AddToStack(hTrkMet2v[ich][isam],samplev[isam]->label,samplev[isam]->color,samplev[isam]->linecol);
      }
    }

//    plotMet2_trk.SetLegend(0.68,0.57,0.84,0.70);
    plotMet2_trk.SetLegend(0.22,0.59,0.42,0.87);
    plotMet2_trk.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotMet2_trk.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotMet2_trk.SetYRange(0.1,1e4);
    plotMet2_trk.SetLogy();
    plotMet2_trk.Draw(c,kTRUE,format,1);

  }

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
//  makeHTML(outputDir);
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("plotWmMet"); 
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hMC, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hMC->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hMC->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.48);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  
  return hDiff;
}

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wmunu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
      
  for(UInt_t ich=0; ich<3; ich++) {
    if(ich==0) htmlfile << "<h3 style=\"text-align:left; color:DD6600;\"><a target=\"_blank\" href=\"w_all.html\">W<sup>&plusmn;</sup></a></h3>" << endl;
    if(ich==1) htmlfile << "<h3 style=\"text-align:left; color:DD6600;\"><a target=\"_blank\" href=\"w_pos.html\">W<sup>+</sup></a></h3>" << endl;
    if(ich==2) htmlfile << "<h3 style=\"text-align:left; color:DD6600;\"><a target=\"_blank\" href=\"w_neg.html\">W<sup>-</sup></a></h3>" << endl;
    
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_ch" << ich << ".png\"><img src=\"plots/lpt_ch" << ich << ".png\" alt=\"plots/lpt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_ch" << ich << ".png\"><img src=\"plots/met_ch" << ich << ".png\" alt=\"plots/met_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_ch" << ich << ".png\"><img src=\"plots/mt_ch" << ich << ".png\" alt=\"plots/mt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
  }
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
  
  for(UInt_t ich=0; ich<3; ich++) {
    if(ich==0) sprintf(htmlfname,"%s/w_all.html",outDir.Data());
    if(ich==1) sprintf(htmlfname,"%s/w_pos.html",outDir.Data());
    if(ich==2) sprintf(htmlfname,"%s/w_neg.html",outDir.Data());
    
    htmlfile.open(htmlfname);
    htmlfile << "<!DOCTYPE html" << endl;
    htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
    htmlfile << "<html>" << endl;
    htmlfile << "<head><title>Wmunu</title></head>" << endl;
    htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_ch" << ich << ".png\"><img src=\"plots/lpt_ch" << ich << ".png\" alt=\"plots/lpt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_ch" << ich << ".png\"><img src=\"plots/met_ch" << ich << ".png\" alt=\"plots/met_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_ch" << ich << ".png\"><img src=\"plots/mt_ch" << ich << ".png\" alt=\"plots/mt_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lptlog_ch" << ich << ".png\"><img src=\"plots/lptlog_ch" << ich << ".png\" alt=\"plots/lptlog_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metlog_ch" << ich << ".png\"><img src=\"plots/metlog_ch" << ich << ".png\" alt=\"plots/metlog_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtlog_ch" << ich << ".png\"><img src=\"plots/mtlog_ch" << ich << ".png\" alt=\"plots/mtlog_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/leta_ch" << ich << ".png\"><img src=\"plots/leta_ch" << ich << ".png\" alt=\"plots/leta_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lphi_ch" << ich << ".png\"><img src=\"plots/lphi_ch" << ich << ".png\" alt=\"plots/lphi_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metphi_ch" << ich << ".png\"><img src=\"plots/metphi_ch" << ich << ".png\" alt=\"plots/metphi_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/npv_ch" << ich << ".png\"><img src=\"plots/npv_ch" << ich << ".png\" alt=\"plots/npv_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">Barrel</h3>" << endl;
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_barrel_ch" << ich << ".png\"><img src=\"plots/lpt_barrel_ch" << ich << ".png\" alt=\"plots/lpt_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_barrel_ch" << ich << ".png\"><img src=\"plots/met_barrel_ch" << ich << ".png\" alt=\"plots/met_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_barrel_ch" << ich << ".png\"><img src=\"plots/mt_barrel_ch" << ich << ".png\" alt=\"plots/mt_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lptlog_barrel_ch" << ich << ".png\"><img src=\"plots/lptlog_barrel_ch" << ich << ".png\" alt=\"plots/lptlog_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metlog_barrel_ch" << ich << ".png\"><img src=\"plots/metlog_barrel_ch" << ich << ".png\" alt=\"plots/metlog_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtlog_barrel_ch" << ich << ".png\"><img src=\"plots/mtlog_barrel_ch" << ich << ".png\" alt=\"plots/mtlog_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lphi_barrel_ch" << ich << ".png\"><img src=\"plots/lphi_barrel_ch" << ich << ".png\" alt=\"plots/lphi_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metphi_barrel_ch" << ich << ".png\"><img src=\"plots/metphi_barrel_ch" << ich << ".png\" alt=\"plots/metphi_barrel_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">Endcap</h3>" << endl;
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lpt_endcap_ch" << ich << ".png\"><img src=\"plots/lpt_endcap_ch" << ich << ".png\" alt=\"plots/lpt_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met_endcap_ch" << ich << ".png\"><img src=\"plots/met_endcap_ch" << ich << ".png\" alt=\"plots/met_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt_endcap_ch" << ich << ".png\"><img src=\"plots/mt_endcap_ch" << ich << ".png\" alt=\"plots/mt_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lptlog_endcap_ch" << ich << ".png\"><img src=\"plots/lptlog_endcap_ch" << ich << ".png\" alt=\"plots/lptlog_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metlog_endcap_ch" << ich << ".png\"><img src=\"plots/metlog_endcap_ch" << ich << ".png\" alt=\"plots/metlog_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtlog_endcap_ch" << ich << ".png\"><img src=\"plots/mtlog_endcap_ch" << ich << ".png\" alt=\"plots/mtlog_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lphi_endcap_ch" << ich << ".png\"><img src=\"plots/lphi_endcap_ch" << ich << ".png\" alt=\"plots/lphi_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metphi_endcap_ch" << ich << ".png\"><img src=\"plots/metphi_endcap_ch" << ich << ".png\" alt=\"plots/metphi_endcap_ch" << ich << ".png\" width=\"100%\"></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;
    htmlfile << "</table>" << endl;
    htmlfile << "<hr />" << endl;
    
    htmlfile << "</body>" << endl;
    htmlfile << "</html>" << endl;
    htmlfile.close();  
  }  
}

// make LaTex table of bkg. to sig. ratios
void printRatiosLatex(vector<double> ratios, ostream& os) {

  os.precision(3);
  for(UInt_t j=0;j<ratios.size();j++)
    os << fixed << ratios[j] << endl;

}

