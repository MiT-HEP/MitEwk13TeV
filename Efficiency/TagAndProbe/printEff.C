#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TGraphAsymmErrors.h>      // graphs
#include <TH2F.h>                   // 2D histograms
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#endif

//=== MAIN MACRO =================================================================================================

void printEff(const TString outdir   = "Data/extra",
              const TString mcfname  = "MC/eff.root",
	      const TString datfname = "Data/eff.root"
) {
  gBenchmark->Start("printEff");

    
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   

  TFile mcfile(mcfname);
  TFile datafile(datfname);
  
  gSystem->mkdir(outdir,true);
    
  TH2F *hMCEff=0,   *hMCErrl=0,   *hMCErrh=0;
  TH2F *hDataEff=0, *hDataErrl=0, *hDataErrh=0;
  
  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");  
  
  hDataEff  = (TH2F*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2F*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2F*)datafile.Get("hErrhEtaPt");
    
  //
  // Print out text file
  //
  ofstream txtfile;
  TString outfname = outdir + TString("/eff_table.txt");
  txtfile.open(outfname.Data());
  assert(txtfile.is_open());
    
  const Int_t nx = hMCEff->GetNbinsX();
  const Int_t ny = hMCEff->GetNbinsY();
  
  txtfile << " pT        ";
  txtfile << " eta           ";
  txtfile << "    MC T&P           ";
  txtfile << "    Data T&P         ";  
  txtfile << "    Scale factor     ";
  txtfile << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------------------------";
  txtfile << endl;
  
  for(Int_t ix=1; ix<=nx; ix++) {
    for(Int_t iy=1; iy<=ny; iy++) {    

      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      
      Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
            
      Double_t scale     = (hDataEff->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleerrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleerrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);

      txtfile << "[" << setw(4) << hMCEff->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << hMCEff->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      txtfile << "[" << setw(6) << hMCEff->GetXaxis()->GetBinLowEdge(ix) << "," << setw(6) << hMCEff->GetXaxis()->GetBinLowEdge(ix+1) << "]";       
      ios_base::fmtflags flags = txtfile.flags();
      txtfile.precision(4);
      txtfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleerrl,scaleerrh);
      txtfile << endl;
      txtfile.flags(flags);
    }
    txtfile << endl;
  }
  txtfile.close();
  
  //
  // Print out LaTeX tables
  //
  ofstream mc_latex;
  outfname = outdir + TString("/mc_eff_table.tex");
  mc_latex.open(outfname.Data());
  assert(mc_latex.is_open());
  
  ofstream dat_latex;
  outfname = outdir + TString("/dat_eff_table.tex");
  dat_latex.open(outfname.Data());
  assert(dat_latex.is_open());  
  
  ofstream sf_latex;
  outfname = outdir + TString("/sf_table.tex");
  sf_latex.open(outfname.Data());
  assert(sf_latex.is_open());
  
  mc_latex << "\\begin{table}[!ht]" << endl;
  mc_latex << "\\begin{center}" << endl;
  mc_latex << "\\begin{tabular}{c";
  for(int ix=0; ix<nx; ix++) mc_latex << "|c";
  mc_latex << "}" << endl;
  mc_latex << "\\hline";
  for(int ix=1; ix<=nx; ix++)
    mc_latex << " & $" << hDataEff->GetXaxis()->GetBinLowEdge(ix) << " < |\\eta| < " << hDataEff->GetXaxis()->GetBinLowEdge(ix+1) << "$";
  mc_latex << "  \\\\" << endl; 
  mc_latex << "\\hline" << endl;  

  dat_latex << "\\begin{table}[!ht]" << endl;
  dat_latex << "\\begin{center}" << endl;
  dat_latex << "\\begin{tabular}{c";
  for(int ix=0; ix<nx; ix++) dat_latex << "|c";
  dat_latex << "}" << endl;
  dat_latex << "\\hline";
  for(int ix=1; ix<=nx; ix++)
    dat_latex << " & $" << hDataEff->GetXaxis()->GetBinLowEdge(ix) << " < |\\eta| < " << hDataEff->GetXaxis()->GetBinLowEdge(ix+1) << "$";
  dat_latex << "  \\\\" << endl; 
  dat_latex << "\\hline" << endl;  
  
  sf_latex << "\\begin{table}[!ht]" << endl;
  sf_latex << "\\begin{center}" << endl;
  sf_latex << "\\begin{tabular}{c";
  for(int ix=0; ix<nx; ix++) sf_latex << "|c";
  sf_latex << "}" << endl;
  sf_latex << "\\hline";
  for(int ix=1; ix<=nx; ix++)
    sf_latex << " & $" << hDataEff->GetXaxis()->GetBinLowEdge(ix) << " < |\\eta| < " << hDataEff->GetXaxis()->GetBinLowEdge(ix+1) << "$";
  sf_latex << "  \\\\" << endl; 
  sf_latex << "\\hline" << endl;  
    
  for(Int_t iy=1; iy<=ny; iy++) {
    mc_latex  << "$" << setw(3) << hMCEff->GetYaxis()->GetBinLowEdge(iy) << " < p_T < " << setw(3) << hMCEff->GetYaxis()->GetBinLowEdge(iy+1) << "$";
    dat_latex << "$" << setw(3) << hMCEff->GetYaxis()->GetBinLowEdge(iy) << " < p_T < " << setw(3) << hMCEff->GetYaxis()->GetBinLowEdge(iy+1) << "$";
    sf_latex  << "$" << setw(3) << hMCEff->GetYaxis()->GetBinLowEdge(iy) << " < p_T < " << setw(3) << hMCEff->GetYaxis()->GetBinLowEdge(iy+1) << "$";
    for(Int_t ix=1; ix<=nx; ix++) {
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      
      Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
            
      Double_t scale     = (hDataEff->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleerrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleerrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
                     
      ios_base::fmtflags flags = mc_latex.flags();
      mc_latex.precision(4);
      mc_latex << " & $" << setw(6) << fixed << mceff << " \\pm " << TMath::Max(mcerrl,mcerrh) << "$";
      mc_latex.flags(flags);

      flags = dat_latex.flags();
      dat_latex.precision(4);
      dat_latex << " & $" << setw(6) << fixed << dataeff << " \\pm " << TMath::Max(dataerrl,dataerrh) << "$";
      dat_latex.flags(flags);
      
      flags = sf_latex.flags();
      sf_latex.precision(4);
      sf_latex << " & $" << setw(6) << fixed << scale << " \\pm " << TMath::Max(scaleerrl,scaleerrh) << "$";
      sf_latex.flags(flags);
    }
    mc_latex  << "  \\\\" << endl;
    dat_latex << "  \\\\" << endl;
    sf_latex  << "  \\\\" << endl;
  }    
  mc_latex << "\\hline" << endl;
  mc_latex << "\\end{tabular}" << endl;
  mc_latex << "\\caption{Write some stuff}" << endl;
  mc_latex << "\\label{tab:mytable}" << endl;
  mc_latex << "\\end{center}" << endl;
  mc_latex << "\\end{table}" << endl;
  mc_latex.close();

  dat_latex << "\\hline" << endl;
  dat_latex << "\\end{tabular}" << endl;
  dat_latex << "\\caption{Write some stuff}" << endl;
  dat_latex << "\\label{tab:mytable}" << endl;
  dat_latex << "\\end{center}" << endl;
  dat_latex << "\\end{table}" << endl;
  dat_latex.close();
  
  sf_latex << "\\hline" << endl;
  sf_latex << "\\end{tabular}" << endl;
  sf_latex << "\\caption{Write some stuff}" << endl;
  sf_latex << "\\label{tab:mytable}" << endl;
  sf_latex << "\\end{center}" << endl;
  sf_latex << "\\end{table}" << endl;
  sf_latex.close();
    
  cout << " <> Output saved to " << outdir << endl;
  
  gBenchmark->Show("printEff"); 
}
