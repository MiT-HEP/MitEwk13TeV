//================================================================================================
//
// Perform fit to extract W->munu signal
//
//  * outputs plots and fit results summary
//
//________________________________________________________________________________________________

  ... (This code can be found by following the instructions here: http://www.cmsaf.mit.edu/twiki/bin/view/CmsHep/WZBosonCrossSection)

//=== MAIN MACRO ================================================================================================= 

void fitWm(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {
	
  ...
  
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  // Remove this constraint by commenting out the "cewk(p,m).setConstant(kTRUE)" lines
  //

  ...
  
  //
  // Construct PDFs for fitting
  //

  ...
  
  // For Combine

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet), hDataMet);
  RooDataHist antiMet("antiMet", "antiMet", RooArgSet(pfmet), hAntiDataMet);
  RooDataHist dataMetp("dataMetp", "dataMetp", RooArgSet(pfmet), hDataMetp);
  RooDataHist antiMetp("antiMetp", "antiMetp", RooArgSet(pfmet), hAntiDataMetp);
  RooDataHist dataMetm("dataMetm", "dataMetm", RooArgSet(pfmet), hDataMetm);
  RooDataHist antiMetm("antiMetm", "antiMetm", RooArgSet(pfmet), hAntiDataMetm);

  RooWorkspace combine_workspace("combine_workspace");

  combine_workspace.import(dataMet);       combine_workspace.import(dataMetp);       combine_workspace.import(dataMetm);
  combine_workspace.import(pdfWm);         combine_workspace.import(pdfWmp);         combine_workspace.import(pdfWmm);
  combine_workspace.import(pdfEWK);        combine_workspace.import(pdfEWKp);        combine_workspace.import(pdfEWKm);
  combine_workspace.import(*(qcd.model));  combine_workspace.import(*(qcdp.model));  combine_workspace.import(*(qcdm.model));
//  combine_workspace.import(qcd);           combine_workspace.import(qcdp);           combine_workspace.import(qcdm);

  combine_workspace.import(pdfWm_5GeV);    combine_workspace.import(pdfWmp_5GeV);    combine_workspace.import(pdfWmm_5GeV);
  combine_workspace.import(pdfWm_10GeV);   combine_workspace.import(pdfWmp_10GeV);   combine_workspace.import(pdfWmm_10GeV);
  combine_workspace.import(pdfWm_15GeV);   combine_workspace.import(pdfWmp_15GeV);   combine_workspace.import(pdfWmm_15GeV);

  combine_workspace.import(antiMet);       combine_workspace.import(antiMetp);       combine_workspace.import(antiMetm);
  combine_workspace.import(apdfWm);        combine_workspace.import(apdfWmp);        combine_workspace.import(apdfWmm);
  combine_workspace.import(apdfEWK);       combine_workspace.import(apdfEWKp);       combine_workspace.import(apdfEWKm);
  combine_workspace.import(*(aqcd.model)); combine_workspace.import(*(aqcdp.model)); combine_workspace.import(*(aqcdm.model));
//  combine_workspace.import(aqcd);          combine_workspace.import(aqcdp);          combine_workspace.import(aqcdm);

  combine_workspace.import(apdfWm_5GeV);  combine_workspace.import(apdfWmp_5GeV);  combine_workspace.import(apdfWmm_5GeV);
  combine_workspace.import(apdfWm_10GeV); combine_workspace.import(apdfWmp_10GeV); combine_workspace.import(apdfWmm_10GeV);
  combine_workspace.import(apdfWm_15GeV); combine_workspace.import(apdfWmp_15GeV); combine_workspace.import(apdfWmm_15GeV);

  combine_workspace.writeToFile("Wmunu_pdfTemplates.root");
  
  RooAbsPdf *wm       = combine_workspace.pdf("wm");
  RooAbsPdf *wm_5GeV  = combine_workspace.pdf("wm_5GeV");
  RooAbsPdf *wm_10GeV = combine_workspace.pdf("wm_10GeV");
  RooAbsPdf *wm_15GeV = combine_workspace.pdf("wm_15GeV");
  RooPlot *plot = combine_workspace.var("pfmet")->frame();  plot->SetTitle("Wmunu - MC Shape Template (Selection)");
  wm->plotOn(plot,RooFit::LineColor(1));
  wm_5GeV->plotOn(plot,RooFit::LineColor(2));
  wm_10GeV->plotOn(plot,RooFit::LineColor(3));
  wm_15GeV->plotOn(plot,RooFit::LineColor(4));
  TCanvas *wmcanvas = new TCanvas();  wmcanvas->cd();  plot->Draw();
  wmcanvas->Print("/home/cmedlock/public_html/WmShapePlot.png");

  RooAbsPdf *wmp       = combine_workspace.pdf("wmp");
  RooAbsPdf *wmp_5GeV  = combine_workspace.pdf("wmp_5GeV");
  RooAbsPdf *wmp_10GeV = combine_workspace.pdf("wmp_10GeV");
  RooAbsPdf *wmp_15GeV = combine_workspace.pdf("wmp_15GeV");
  RooPlot *plotp = combine_workspace.var("pfmet")->frame();  plotp->SetTitle("Wmunu_p - MC Shape Template (Selection)");
  wmp->plotOn(plotp,RooFit::LineColor(1));
  wmp_5GeV->plotOn(plotp,RooFit::LineColor(2));
  wmp_10GeV->plotOn(plotp,RooFit::LineColor(3));
  wmp_15GeV->plotOn(plotp,RooFit::LineColor(4));
  TCanvas *wmpcanvas = new TCanvas();  wmpcanvas->cd();  plotp->Draw();
  wmpcanvas->Print("/home/cmedlock/public_html/WmpShapePlot.png");

  RooAbsPdf *wmm       = combine_workspace.pdf("wmm");
  RooAbsPdf *wmm_5GeV  = combine_workspace.pdf("wmm_5GeV");
  RooAbsPdf *wmm_10GeV = combine_workspace.pdf("wmm_10GeV");
  RooAbsPdf *wmm_15GeV = combine_workspace.pdf("wmm_15GeV");
  RooPlot *plotm = combine_workspace.var("pfmet")->frame();  plotm->SetTitle("Wmunu_m - MC Shape Template (Selection)");
  wmm->plotOn(plotm,RooFit::LineColor(1));
  wmm_5GeV->plotOn(plotm,RooFit::LineColor(2));
  wmm_10GeV->plotOn(plotm,RooFit::LineColor(3));
  wmm_15GeV->plotOn(plotm,RooFit::LineColor(4));
  TCanvas *wmmcanvas = new TCanvas();  wmmcanvas->cd();  plotm->Draw();
  wmmcanvas->Print("/home/cmedlock/public_html/WmmShapePlot.png");

  RooAbsPdf *awm       = combine_workspace.pdf("awm");
  RooAbsPdf *awm_5GeV  = combine_workspace.pdf("awm_5GeV");
  RooAbsPdf *awm_10GeV = combine_workspace.pdf("awm_10GeV");
  RooAbsPdf *awm_15GeV = combine_workspace.pdf("awm_15GeV");
  RooPlot *aplot = combine_workspace.var("pfmet")->frame();  aplot->SetTitle("Wmunu - MC Shape Template (Anti-Selection)");
  awm->plotOn(aplot,RooFit::LineColor(1));
  awm_5GeV->plotOn(aplot,RooFit::LineColor(2));
  awm_10GeV->plotOn(aplot,RooFit::LineColor(3));
  awm_15GeV->plotOn(aplot,RooFit::LineColor(4));
  TCanvas *awmcanvas = new TCanvas();  awmcanvas->cd();  aplot->Draw();
  awmcanvas->Print("/home/cmedlock/public_html/aWmShapePlot.png");

  RooAbsPdf *awmp       = combine_workspace.pdf("awmp");
  RooAbsPdf *awmp_5GeV  = combine_workspace.pdf("awmp_5GeV");
  RooAbsPdf *awmp_10GeV = combine_workspace.pdf("awmp_10GeV");
  RooAbsPdf *awmp_15GeV = combine_workspace.pdf("awmp_15GeV");
  RooPlot *aplotp = combine_workspace.var("pfmet")->frame();  aplotp->SetTitle("Wmunu_m - MC Shape Template (Anti-Selection)");
  awmp->plotOn(aplotp,RooFit::LineColor(1));
  awmp_5GeV->plotOn(aplotp,RooFit::LineColor(2));
  awmp_10GeV->plotOn(aplotp,RooFit::LineColor(3));
  awmp_15GeV->plotOn(aplotp,RooFit::LineColor(4));
  TCanvas *awmpcanvas = new TCanvas();  awmpcanvas->cd();  aplotp->Draw();
  awmpcanvas->Print("/home/cmedlock/public_html/aWmpShapePlot.png");

  RooAbsPdf *awmm       = combine_workspace.pdf("awmm");
  RooAbsPdf *awmm_5GeV  = combine_workspace.pdf("awmm_5GeV");
  RooAbsPdf *awmm_10GeV = combine_workspace.pdf("awmm_10GeV");
  RooAbsPdf *awmm_15GeV = combine_workspace.pdf("awmm_15GeV");
  RooPlot *aplotm = combine_workspace.var("pfmet")->frame();  aplotm->SetTitle("Wmunu_m - MC Shape Template (Anti-Selection)");
  awmm->plotOn(aplotm,RooFit::LineColor(1));
  awmm_5GeV->plotOn(aplotm,RooFit::LineColor(2));
  awmm_10GeV->plotOn(aplotm,RooFit::LineColor(3));
  awmm_15GeV->plotOn(aplotm,RooFit::LineColor(4));
  TCanvas *awmmcanvas = new TCanvas();  awmmcanvas->cd();  aplotm->Draw();
  awmmcanvas->Print("/home/cmedlock/public_html/aWmmShapePlot.png");

  ...
  
}
