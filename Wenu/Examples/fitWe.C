//================================================================================================
//
// Perform fit to extract W->enu signal
//
//  * outputs plots and fit results summary
//
//________________________________________________________________________________________________

  ... (This code can be found by following the instructions here: http://www.cmsaf.mit.edu/twiki/bin/view/CmsHep/WZBosonCrossSection)

//=== MAIN MACRO ================================================================================================= 

void fitWe(const TString  outputDir,   // output directory
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
) {

  ...
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
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

  RooDataHist dataMet("dataMet", "dataMet", RooArgSet(pfmet),hDataMet);
  RooDataHist dataMetp("dataMetp","dataMetp",RooArgSet(pfmet),hDataMetp);
  RooDataHist dataMetm("dataMetm","dataMetm",RooArgSet(pfmet),hDataMetm);

  RooWorkspace combine_workspace("combine_workspace");

  combine_workspace.import(dataMet);
  combine_workspace.import(dataMetp);
  combine_workspace.import(dataMetm);

  combine_workspace.import(pdfWe);
  combine_workspace.import(pdfWep);
  combine_workspace.import(pdfWem);
  combine_workspace.import(pdfEWK);
  combine_workspace.import(pdfEWKp);
  combine_workspace.import(pdfEWKm);
  combine_workspace.import(*(qcd.model));
  combine_workspace.import(*(qcdp.model));
  combine_workspace.import(*(qcdm.model));

  combine_workspace.writeToFile("Wenu_pdfTemplates.root");

  ...
