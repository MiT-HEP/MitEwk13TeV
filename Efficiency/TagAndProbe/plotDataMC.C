#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBenchmark.h>             // class to track macro running statistics
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>      // graph class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <fstream>
#include <sstream>

#include "CPlot.hh"	     // helper class for plots
#include "MitStyleRemix.hh"  // style settings for drawing
#endif

void plotDataMC(const TString outdir   = "Data/extra",
                const TString mcfname  = "MC/eff.root",
	        const TString datfname = "Data/eff.root",
                const TString fname    = "noNameGiven",
                const double efflow    = 0.70,
                const double effhigh   = 1.20,
                const double lumi      = 7.3,
                const TString conf     = "muhlt.bins",
                const TString xaxislabel="",   // 'Supercluster' or 'Muon'
                const TString yaxislabel=""    // use different labels for each efficiency
) {
  gBenchmark->Start("plotDataMC");
  
  // bin edges for kinematic variables
  vector<Double_t> ptBinEdgesv;
  vector<Double_t> etaBinEdgesv;
  
  //
  // parse binning file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  Int_t opts[6];
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    Double_t edge;
    stringstream ss(line);
    if(state==0) {
      ss >> opts[0] >> opts[1] >> opts[2] >> opts[3] >> opts[4] >> opts[5];
    } else {
      ss >> edge;
      if(state==1)      { ptBinEdgesv.push_back(edge);  }
      else if(state==2) { etaBinEdgesv.push_back(edge); }
   }
  }
  ifs.close();

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  vector<TString> etalabelv;
  vector<TString> ptlabelv;
  
  // eta bins
  etalabelv.push_back("0 < |#eta| < 1.2");
  etalabelv.push_back("1.2 < |#eta| < 2.4");
  
  CPlot::sOutDir = outdir;
  TString format = "png";

  // y-axis range for scale factors
  Double_t scalelow = 0.80, scalehigh = 1.05;
  
  char lumitext[100]; // lumi label
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);    

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TFile mcfile(mcfname);      
  TFile datafile(datfname);
  
  TH2F *hMCEff=0,   *hMCErrl=0,   *hMCErrh=0;
  TH2F *hDataEff=0, *hDataErrl=0, *hDataErrh=0;
  
  TGraphAsymmErrors *grMCEffEta=0, *grDataEffEta=0, *grScaleEta=0;
  TGraphAsymmErrors *grMCEffPt=0,  *grDataEffPt=0,  *grScalePt=0;
  
  vector<TGraphAsymmErrors*> mceff_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> eff_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> mceff_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> eff_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> scale_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> scale_vs_eta_per_ptv;
  
  grMCEffEta = (TGraphAsymmErrors*)mcfile.Get("grEffEta");
  grMCEffPt  = (TGraphAsymmErrors*)mcfile.Get("grEffPt");

  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");  

  grDataEffEta = (TGraphAsymmErrors*)datafile.Get("grEffEta");
  grDataEffPt  = (TGraphAsymmErrors*)datafile.Get("grEffPt");
  
  hDataEff  = (TH2F*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2F*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2F*)datafile.Get("hErrhEtaPt");
    
  if(grMCEffEta && grDataEffEta) {
    grScaleEta = new TGraphAsymmErrors(grMCEffEta->GetN());
    for(Int_t i=0; i<grMCEffEta->GetN(); i++) {
      Double_t mcval   = grMCEffEta->GetY()[i];
      Double_t dataval = grDataEffEta->GetY()[i];
      if(i==2 || i==8) { grMCEffEta->RemovePoint(i); grDataEffEta->RemovePoint(i); } // For removing the ECAL gap bins for electrons
      Double_t scale   = dataval/mcval;
      grScaleEta->SetPoint(i,grMCEffEta->GetX()[i],scale);
      
      Double_t mcerrl   = grMCEffEta->GetErrorYlow(i);
      Double_t mcerrh   = grMCEffEta->GetErrorYhigh(i);
      Double_t dataerrl = grDataEffEta->GetErrorYlow(i);
      Double_t dataerrh = grDataEffEta->GetErrorYhigh(i);
      grScaleEta->SetPointError(i, 0, 0,
	  		        scale*sqrt(mcerrl*mcerrl/mcval/mcval + dataerrl*dataerrl/dataval/dataval),
			        scale*sqrt(mcerrh*mcerrh/mcval/mcval + dataerrh*dataerrh/dataval/dataval));
    }
  }
  
  if(grMCEffPt && grDataEffPt) {
    grScalePt = new TGraphAsymmErrors(grMCEffPt->GetN());
    for(Int_t i=0; i<grMCEffPt->GetN(); i++) {
      Double_t mcval   = grMCEffPt->GetY()[i];
      Double_t dataval = grDataEffPt->GetY()[i];
      Double_t scale   = dataval/mcval;
      grScalePt->SetPoint(i, grMCEffPt->GetX()[i],scale);
      //if(i==grMCEffPt->GetN()-1)
        //grScalePt->SetPoint(i,165,scale);
    
      Double_t mcerrl   = grMCEffPt->GetErrorYlow(i);
      Double_t mcerrh   = grMCEffPt->GetErrorYhigh(i);
      Double_t dataerrl = grDataEffPt->GetErrorYlow(i);
      Double_t dataerrh = grDataEffPt->GetErrorYhigh(i);
      grScalePt->SetPointError(i, 0, 0,
			        scale*sqrt(mcerrl*mcerrl/mcval/mcval + dataerrl*dataerrl/dataval/dataval),
			        scale*sqrt(mcerrh*mcerrh/mcval/mcval + dataerrh*dataerrh/dataval/dataval));
    }
  }
  
  if(hMCEff->GetEntries()>0 && hDataEff->GetEntries()>0) {
    const Int_t nx = hMCEff->GetNbinsX();
    const Int_t ny = hMCEff->GetNbinsY(); 

    for(Int_t ix=1; ix<=nx; ix++) {
      Double_t xval[ny], xerr[ny];
      Double_t mceffval[ny], mcefferrl[ny], mcefferrh[ny];
      Double_t effval[ny],   efferrl[ny],   efferrh[ny];
      Double_t scaleval[ny], scaleerrl[ny], scaleerrh[ny];
      
      for(Int_t iy=1; iy<=ny; iy++) {
        xval[iy-1] = 0.5*(hMCEff->GetYaxis()->GetBinLowEdge(iy+1) + hMCEff->GetYaxis()->GetBinLowEdge(iy));
        xerr[iy-1] = 0.5*(hMCEff->GetYaxis()->GetBinLowEdge(iy+1) - hMCEff->GetYaxis()->GetBinLowEdge(iy));
	if(iy==ny) {
	  xval[iy-1] = 125;
          xerr[iy-1] = 25;
        }
	
        Double_t mceff  = hMCEff->GetBinContent(hMCEff->GetBin(ix,iy));
        Double_t mcerrl = hMCErrl->GetBinContent(hMCErrl->GetBin(ix,iy));
        Double_t mcerrh = hMCErrh->GetBinContent(hMCErrh->GetBin(ix,iy));
	mceffval[iy-1]  = mceff;
        mcefferrl[iy-1] = mcerrl;
        mcefferrh[iy-1] = mcerrh;
        
        Double_t dataeff  = hDataEff->GetBinContent(hDataEff->GetBin(ix,iy));
        Double_t dataerrl = hDataErrl->GetBinContent(hDataErrl->GetBin(ix,iy));
        Double_t dataerrh = hDataErrh->GetBinContent(hDataErrh->GetBin(ix,iy));
        effval[iy-1]  = dataeff;
        efferrl[iy-1] = dataerrl;
        efferrh[iy-1] = dataerrh;
        
        Double_t scale = dataeff/mceff;
        scaleval[iy-1]  = scale;
        scaleerrl[iy-1] = (scale>0) ? scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff) : 0;
        scaleerrh[iy-1] = (scale>0) ? scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff) : 0;
      }
      
      mceff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,mceffval,xerr,xerr,mcefferrl,mcefferrh));
      eff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,effval,xerr,xerr,efferrl,efferrh));
      scale_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,scaleval,0,0,scaleerrl,scaleerrh));
    }  
  }

    
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  if(grMCEffEta && grDataEffEta) {
    CPlot plotEffEta("effeta_"+fname,"",xaxislabel+" #eta",yaxislabel+" efficiency");
    plotEffEta.AddGraph(grMCEffEta,   "MC","",  kRed, kOpenSquare);
    plotEffEta.AddGraph(grDataEffEta,"Data","",kBlue,kFullDotLarge);
    plotEffEta.SetYRange(efflow,effhigh);
    //plotEffEta.AddTextBox(lumitext,0.58,0.70,0.93,0.76,0);
    plotEffEta.AddTextBox(lumitext,0.62,0.79,0.95,0.87,0);
    plotEffEta.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotEffEta.Draw(c,kTRUE,format);
    
    CPlot plotScaleEta("scaleeta_"+fname,"","|#eta|","scale factor");
    plotScaleEta.AddGraph(grScaleEta,"",kBlue,kFullDotLarge);
    plotScaleEta.AddLine(0,1.0,2.7,1.0,kBlack,7);
    plotScaleEta.SetYRange(scalelow,scalehigh);
    plotScaleEta.SetXRange(0,2.7);
    plotScaleEta.Draw(c,kTRUE,format);
  }
  
  if(grMCEffPt && grDataEffPt) {
    CPlot plotEffPt("effpt_"+fname,"",xaxislabel+" p_{T} [GeV]",yaxislabel+" efficiency");
    plotEffPt.AddGraph(grMCEffPt,   "MC","",  kRed, kOpenSquare);
    plotEffPt.AddGraph(grDataEffPt,"Data","",kBlue,kFullDotLarge);
    plotEffPt.SetYRange(efflow,effhigh);
    plotEffPt.SetXRange(20,85);
    //plotEffPt.AddTextBox(lumitext,0.58,0.70,0.93,0.76,0);
    plotEffPt.AddTextBox(lumitext,0.62,0.79,0.95,0.87,0);
    plotEffPt.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotEffPt.Draw(c,kTRUE,format);
  
    CPlot plotScalePt("scalept_"+fname,"","p_{T} [GeV/c]","scale factor");
    plotScalePt.AddGraph(grScalePt,"",kBlue,kFullDotLarge);
    plotScalePt.AddLine(0,1.0,150,1.0,kBlack,7);
    plotScalePt.SetYRange(scalelow,scalehigh);
    plotScalePt.Draw(c,kTRUE,format);
  }
  
  if(mceff_vs_pt_per_etav.size()>0 && eff_vs_pt_per_etav.size()>0) {

    // To make an eta-phi scale factor table using LaTex
    ofstream latexfile;
    char latexfname[100];    
    sprintf(latexfname,"%s/scalefactors.txt",outdir.Data());
    latexfile.open(latexfname);
    assert(latexfile.is_open());

    for(int a=0; a<etaBinEdgesv.size()-1; a++) {
      if(a==0) latexfile << "& ";
      latexfile << "& $" << etaBinEdgesv[a] << "< \\eta<" << etaBinEdgesv[a+1] << "$ ";
    }
    latexfile << "\\\\";
    latexfile << endl;

    vector<double> sfv;
    vector<double> sf_errv;
    for(UInt_t w=0; w<scale_vs_pt_per_etav.size(); w++) {
      TGraphAsymmErrors* g = (TGraphAsymmErrors*)scale_vs_pt_per_etav[w];
      for(UInt_t d=0; d<g->GetN(); d++) {
        double pt, sf;
        g->GetPoint(d,pt,sf);
        sfv.push_back(sf);
        double sf_err = g->GetErrorY(d);
        sf_errv.push_back(sf_err);
        //cout<<"eta from "<<etaBinEdgesv[w]<<" to "<<etaBinEdgesv[w+1]<<" and pt from "<<ptBinEdgesv[d]<<" to "<<ptBinEdgesv[d+1]<<" ---> scalefactor is "<<sf<<" +/- "<<sf_err<<endl;
      }
    }

    // For each bin in pT
    //cout<<"ptBinEdgesv.size() is "<<ptBinEdgesv.size()<<endl;
    for(UInt_t b=0; b<ptBinEdgesv.size()-1; b++) {
      //cout<<"pt bin "<<ptBinEdgesv[b]<<" to "<<ptBinEdgesv[b+1]<<endl;
      latexfile << "& $" << ptBinEdgesv[b] << "<p_{T}<" << ptBinEdgesv[b+1] << "$ ";
      // Get scale factors for each bin in eta
      for(UInt_t c=0; c<etaBinEdgesv.size()-1; c++) {
        //cout<<"---eta bin "<<etaBinEdgesv[c]<<" to "<<etaBinEdgesv[c+1]<<": b = "<<sfv[b+c*(ptBinEdgesv.size()-1)]<<" ± "<<sf_errv[b+c*(ptBinEdgesv.size()-1)]<<endl;
        ios_base::fmtflags flags = latexfile.flags();
	latexfile.precision(3);
        latexfile << "& $" << fixed << sfv[b+c*(ptBinEdgesv.size()-1)] << " \\pm " << sf_errv[b+c*(ptBinEdgesv.size()-1)] << "$ ";
        latexfile.flags(flags);
      }
      latexfile << " \\\\" << endl;
     }

  }

  gBenchmark->Show("plotDataMC");
}
