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
#include <TGaxis.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <fstream>
#include <sstream>

#include "CPlot.hh"	     // helper class for plots
#include "MitStyleRemix.hh"  // style settings for drawing
#endif

void plotDataMC_singlepTbins(const TString outdir   = "./DataMC",
                const TString mcfname  = "MC/eff.root",
	        const TString datfname = "Data/eff.root",
                const TString fname    = "noNameGiven",
                const double efflow    = 0.70,
                const double effhigh   = 1.20,
                const double lumi      = 7.3,
                const TString conf     = "mupteta.bins",
                const TString xaxislabel="",   // 'Supercluster' or 'Muon'
                const TString yaxislabel="",   // use different labels for each efficiency
                const bool doAbsEta=kFALSE
) {
  gBenchmark->Start("plotDataMC_singlepTbins");

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

  TString format = "png";

  char lumitext[300]; // lumi label
  sprintf(lumitext,"%.0f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);    

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TFile mcfile(mcfname);      
  TFile datafile(datfname);
  
  TH2F *hMCEff=0,   *hMCErrl=0,   *hMCErrh=0;
  TH2F *hDataEff=0, *hDataErrl=0, *hDataErrh=0;

  vector<TGraphAsymmErrors*> mceff_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> mceff_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> eff_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> eff_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> scale_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> scale_vs_pt_per_etav;
  
  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");  

  hDataEff  = (TH2F*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2F*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2F*)datafile.Get("hErrhEtaPt");

  if(hMCEff->GetEntries()>0 && hDataEff->GetEntries()>0) {

    const Int_t nx = hMCEff->GetNbinsX();
    const Int_t ny = hMCEff->GetNbinsY(); 
    
    // for plotting
    for(Int_t iy=1; iy<=ny; iy++) {
      Double_t xval[nx], xerr[nx];
      Double_t mceffval[nx], mcefferrl[nx], mcefferrh[nx];
      Double_t effval[nx],   efferrl[nx],   efferrh[nx];
      Double_t scaleval[nx], scaleerrl[nx], scaleerrh[nx];

      for(Int_t ix=1; ix<=nx; ix++) {
        xval[ix-1] = 0.5*(hMCEff->GetXaxis()->GetBinLowEdge(ix) + hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
        xerr[ix-1] = 0.5*(hMCEff->GetXaxis()->GetBinLowEdge(ix) - hMCEff->GetXaxis()->GetBinLowEdge(ix+1));

        Double_t mceff  = hMCEff->GetBinContent(ix,iy);
        Double_t mcerrl = hMCErrl->GetBinContent(ix,iy);
        Double_t mcerrh = hMCErrh->GetBinContent(ix,iy);
	mceffval[ix-1]  = mceff;
        mcefferrl[ix-1] = mcerrl;
        mcefferrh[ix-1] = mcerrh;

        Double_t dataeff  = hDataEff->GetBinContent(ix,iy);
        Double_t dataerrl = hDataErrl->GetBinContent(ix,iy);
        Double_t dataerrh = hDataErrh->GetBinContent(ix,iy);
        effval[ix-1]  = dataeff;
        efferrl[ix-1] = dataerrl;
        efferrh[ix-1] = dataerrh;
        
        Double_t scale = dataeff/mceff;
        scaleval[ix-1]  = scale;
        scaleerrl[ix-1] = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
        scaleerrh[ix-1] = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      }
      
      mceff_vs_eta_per_ptv.push_back(new TGraphAsymmErrors(nx,xval,mceffval,xerr,xerr,mcefferrl,mcefferrh));
      eff_vs_eta_per_ptv.push_back(new TGraphAsymmErrors(nx,xval,effval,xerr,xerr,efferrl,efferrh));
      scale_vs_eta_per_ptv.push_back(new TGraphAsymmErrors(nx,xval,scaleval,xerr,xerr,scaleerrl,scaleerrh));
    }

    // for the latex table
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

  TCanvas *c = MakeCanvas("c","c",1200,1200);

  c->SetLeftMargin(0.16);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.16);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.30556,1.0,1.0);
  c->cd(1)->SetTopMargin(0.08);
  c->cd(1)->SetBottomMargin(0.0);//0.05);
  c->cd(1)->SetLeftMargin(0.16);  
  c->cd(1)->SetRightMargin(0.05);

  c->cd(2);
  c->cd(2)->SetPad(0.0,0.0,1.0,0.30556);
  c->cd(2)->SetLeftMargin(0.16);
  c->cd(2)->SetRightMargin(0.05);
  c->cd(2)->SetTopMargin(0.0);//0.05);
  c->cd(2)->SetBottomMargin(0.35);

  TGaxis::SetMaxDigits(3);

  // Data/MC legend
  TLegend* legend = new TLegend(0.50,0.24,0.93,0.4);

  // Luminosity and 'CMS Preliminary' text boxes
  TPaveText *lumitb = new TPaveText(0.62,0.92,0.94,0.99,"NDC");
  lumitb->SetBorderSize(0);
  lumitb->SetFillColor(0);
  lumitb->AddText(lumitext);

  TPaveText *cmstb = new TPaveText(0.19,0.83,0.54,0.89,"NDC");
  cmstb->SetBorderSize(0);
  cmstb->SetFillColor(0);
  cmstb->AddText("CMS Preliminary");

  if(mceff_vs_eta_per_ptv.size()>0 && eff_vs_eta_per_ptv.size()>0) {
    for(UInt_t ig=0; ig<mceff_vs_eta_per_ptv.size(); ig++) {

      double etalow = -2.5, etahigh = 2.5;
      if(xaxislabel.CompareTo("Supercluster")==0) { 
        etalow = -2.5; etahigh = 2.5;
        if(doAbsEta) etalow = 0.0;
      } else if(xaxislabel.CompareTo("Muon")==0) {
        etalow = -2.4; etahigh = 2.4;
        if(doAbsEta) etalow = 0.0;
      }

      // Format axes
      TH1D* hDummyEff = new TH1D("hDummyEff","",10,etalow,etahigh);
      hDummyEff->GetXaxis()->SetLabelSize(0);
      hDummyEff->GetXaxis()->SetLabelOffset(0.005);
      hDummyEff->GetXaxis()->SetTitleSize(0.05);
      hDummyEff->GetXaxis()->SetTitleOffset(1.4);
      hDummyEff->GetXaxis()->SetTickLength(0.03);
      hDummyEff->GetXaxis()->SetNdivisions(5,5,0,kFALSE);

      hDummyEff->GetYaxis()->SetLabelSize(0.05);
      hDummyEff->GetYaxis()->SetLabelOffset(0.005);
      hDummyEff->GetYaxis()->SetTitle(yaxislabel+" efficiency");
      hDummyEff->GetYaxis()->SetTitleSize(0.05);
      hDummyEff->GetYaxis()->SetTitleOffset(1.4);
      hDummyEff->GetYaxis()->SetTickLength(0.03);
      hDummyEff->GetYaxis()->SetRangeUser(efflow,effhigh);

      TH1D* hDummyScale = new TH1D("hDummyScale","",10,etalow,etahigh);

      hDummyScale->GetXaxis()->SetLabelSize(0.05*(0.69444/0.30556));
      hDummyScale->GetXaxis()->SetLabelOffset(0.005/0.30556);
      hDummyScale->GetXaxis()->SetTitleSize(0.05*(0.69444/0.30556));
      hDummyScale->GetXaxis()->SetTitleOffset(1.4);
      hDummyScale->GetXaxis()->SetTitle(xaxislabel+" #eta");
      if(doAbsEta) hDummyScale->GetXaxis()->SetTitle(xaxislabel+" |#eta|");
      hDummyScale->GetXaxis()->SetTickLength(0.03);//*(0.69444/0.30556));
      hDummyScale->GetXaxis()->SetNdivisions(5,5,0,kFALSE);

      hDummyScale->GetYaxis()->SetLabelSize(0.05*(0.69444/0.30556));
      hDummyScale->GetYaxis()->SetLabelOffset(0.005*(0.69444/0.30556));
      hDummyScale->GetYaxis()->SetTitleSize(0.05*(0.69444/0.30556));
      hDummyScale->GetYaxis()->SetTitle("Data / MC");
      hDummyScale->GetYaxis()->CenterTitle();
      hDummyScale->GetYaxis()->SetTitleOffset(0.6);
      hDummyScale->GetYaxis()->SetTickLength(0.03);
      hDummyScale->GetYaxis()->SetRangeUser(0.8,1.2);
      hDummyScale->GetYaxis()->SetNdivisions(5);

      // Format efficiency graphs
      TGraphAsymmErrors *grMCEffEta   = mceff_vs_eta_per_ptv[ig];
      TGraphAsymmErrors *grDataEffEta = eff_vs_eta_per_ptv[ig];

      double w, d;
      for(int l=0;l<grMCEffEta->GetN();l++) {
        grMCEffEta->GetPoint(l,w,d);
        grMCEffEta->SetPointEXlow(l,grMCEffEta->GetErrorX(l));
        grMCEffEta->SetPointEXhigh(l,grMCEffEta->GetErrorX(l));
        // Don't draw error bars above 1
        if(d+grMCEffEta->GetErrorYhigh(l)>1) grMCEffEta->SetPointEYhigh(l,0);

        grDataEffEta->GetPoint(l,w,d);
        grDataEffEta->SetPointEXlow(l,grDataEffEta->GetErrorX(l));
        grDataEffEta->SetPointEXhigh(l,grDataEffEta->GetErrorX(l));
        // Don't draw error bars above 1
        if(d+grDataEffEta->GetErrorYhigh(l)>1) grDataEffEta->SetPointEYhigh(l,0);
      }

      grMCEffEta->SetMarkerColor(kRed);
      grDataEffEta->SetMarkerColor(kBlue);
      grMCEffEta->SetLineColor(kRed);
      grDataEffEta->SetLineColor(kBlue);
      grMCEffEta->SetLineWidth(2);
      grDataEffEta->SetLineWidth(2);
      grMCEffEta->SetMarkerStyle(kOpenSquare);
      grDataEffEta->SetMarkerStyle(kFullDotLarge);
      grMCEffEta->SetMarkerSize(1.2);
      grDataEffEta->SetMarkerSize(1.2);

      TGraphAsymmErrors *grScaleEta = scale_vs_eta_per_ptv[ig];

      // Remove ECAL gap
      if(xaxislabel.CompareTo("Supercluster")==0) {
        grMCEffEta->RemovePoint(2);   grMCEffEta->RemovePoint(8);
        grDataEffEta->RemovePoint(2); grDataEffEta->RemovePoint(8);
        grScaleEta->RemovePoint(2);   grScaleEta->RemovePoint(8);
        // Comment out the above 3 lines and uncomment the next 3 lines for electron supercluster efficiencies
        //grMCEffEta->RemovePoint(1);
        //grDataEffEta->RemovePoint(1);
        //grScaleEta->RemovePoint(1);
      }

      // Add entries to legend
      if(ig==0) {
        legend->AddEntry(grMCEffEta,"MC","L");
        legend->AddEntry(grDataEffEta,"Data","L");
      }

      // Add pT bin label
      char binlabely[100];
      sprintf(binlabely,"%i GeV/c < p_{T} < %i GeV/c",Int_t(ptBinEdgesv[ig]),Int_t(ptBinEdgesv[ig+1]));
      TPaveText *ptbintext = new TPaveText(0.21,0.75,0.51,0.80,"NDC");
      ptbintext->SetTextColor(kBlack);
      ptbintext->SetFillStyle(0);
      ptbintext->SetBorderSize(0);
      ptbintext->SetTextSize(0.042);
      ptbintext->AddText(binlabely);

      // Draw efficiency plot
      c->cd(1);
      hDummyEff->Draw("axis");
      grMCEffEta->Draw("p same");
      grDataEffEta->Draw("p same");
      lumitb->Draw("same");
      cmstb->Draw("same");
      legend->Draw("same");
      ptbintext->Draw("same");

      // Add dashed line at 1 in scale factor graph
      TLine *line = new TLine(-2.5,1,2.5,1);
      line->SetLineColor(kBlack);
      line->SetLineStyle(3);
      line->SetLineWidth(2);

      // Draw ratio plot
      c->cd(2);
      hDummyScale->Draw("axis");
      gPad->SetTickx(1);
      grScaleEta->Draw("e p same");
      line->Draw("same");

      char pname[300];
      sprintf(pname,outdir+"/"+fname+"_eta_pt%i.",ig);

      gSystem->mkdir(outdir,true);
      c->SaveAs(pname+format);

      delete hDummyEff;
      delete hDummyScale;
    }

    // To make an eta-phi scale factor table using LaTex
    ofstream latexfile;
    char latexfname[300];    
    sprintf(latexfname,"%s/%s_scalefactors.txt",outdir.Data(),fname.Data());
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

  gBenchmark->Show("plotDataMC_singlepTbins");
}
