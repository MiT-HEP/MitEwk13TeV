//================================================================================================
//
// Compute cross sections and produce summary plots
//
//________________________________________________________________________________________________

#include "MitStyleRemix.hh"         // style settings for drawing

// root -q -L xsecSummary.C

void roundXsec(Double_t x, Double_t y, Double_t z, Double_t aa, Double_t tmp[4])
{
  int iPres = 10;;
  tmp[0] = x;//TMath::Nint(x/Double_t(iPres))*iPres;
  tmp[1] =y;// TMath::Nint(y/Double_t(iPres))*iPres;
  tmp[2] = z;//TMath::Nint(z/Double_t(iPres))*iPres;
  tmp[3]=  aa;//TMath::Nint(aa/Double_t(iPres))*iPres;
}

void roundXsec(const Double_t &x, const Double_t &y,Double_t tmp[2])
{
  int iPres = 10;
  tmp[0] =x;// TMath::Nint(x/Double_t(iPres))*iPres;
  tmp[1] = y;//TMath::Nint(y/Double_t(iPres))*iPres;
}

void xsecSummary_ele()
{   
  //TRandom* gRandom = new TRandom3();
  double rand = 0.901493; //gRandom->Gaus(1.0,0.05); 
  //--------------------------------------------------------------------------------------------------------------

  const Double_t theo_wp = 11571;   const Double_t theo_wp_unc = 28;
  const Double_t theo_wm = 8550;   const Double_t theo_wm_unc = 20;
  const Double_t theo_w  = 20121;   const Double_t theo_w_unc  = 47;
  const Double_t theo_z  = 1944;   const Double_t theo_z_unc  = 14;
  const Double_t theo_wr = 1.35322;   const Double_t theo_wr_unc = 0.001;
  const Double_t theo_wpr= 5.95535;   const Double_t theo_wpr_unc= 0.005;
  const Double_t theo_wmr= 4.4006;   const Double_t theo_wmr_unc= 0.0033;
  const Double_t theo_wz = 10.35541;  const Double_t theo_wz_unc = 0.006;

  const Double_t th_wp = 11571*0.0070;  
  const Double_t th_wm = 8550*0.0063;  
  const Double_t th_w  = 20121*0.0062;  
  const Double_t th_z  = 1944*0.0048;  
  const Double_t th_wr = 1.35322*0.0037;  
  const Double_t th_wpr= 5.95535*0.0071;  
  const Double_t th_wmr= 4.4006*0.0057; 
  const Double_t th_wz = 10.35541*0.0064; 

  // input parameter
  double blind=1;
  Double_t xsec_wp = 12167; blind=theo_wp/xsec_wp; xsec_wp=xsec_wp*blind;  const Double_t xsec_wp_stat = 0;   const Double_t xsec_wp_sys =sqrt(pow(91*blind,2)+th_wp*th_wp);   const Double_t xsec_wp_lumi = 207*blind;
  Double_t xsec_wm = 9006;   blind=theo_wm/xsec_wm; xsec_wm=xsec_wm*blind;   const Double_t xsec_wm_stat = 0;   const Double_t xsec_wm_sys = sqrt(pow(58*blind,2)+th_wm*th_wm);   const Double_t xsec_wm_lumi = 153*blind;
  Double_t xsec_w  = 21166;  blind=theo_w/xsec_w; xsec_w=xsec_w*blind;  const Double_t xsec_w_stat  = 0;   const Double_t xsec_w_sys  = sqrt(pow(133*blind,2)+th_w*th_w);   const Double_t xsec_w_lumi  = 360*blind;
  Double_t xsec_z  = 1986;   blind=theo_z/xsec_z; xsec_z=xsec_z*blind;   const Double_t xsec_z_stat  = 0;    const Double_t xsec_z_sys  = sqrt(pow(15*blind,2)+th_z*th_z);    const Double_t xsec_z_lumi  = 34*blind;
  Double_t xsec_wr = 1.35123; blind=theo_wr/xsec_wr; xsec_wr=xsec_wr*blind;   const Double_t xsec_wr_stat = 0.0;   const Double_t xsec_wr_sys = sqrt(pow(0.0101492*blind,2)+th_wr*th_wr); 
  Double_t xsec_wpr= 6.12483;  blind=theo_wpr/xsec_wpr; xsec_wpr=xsec_wpr*blind;  const Double_t xsec_wpr_stat= 0.0;   const Double_t xsec_wpr_sys= sqrt(pow(0.05215*blind,2)+th_wmr*th_wmr); 
  Double_t xsec_wmr= 4.53371;  blind=theo_wmr/xsec_wmr; xsec_wmr=xsec_wmr*blind;  const Double_t xsec_wmr_stat= 0.0;   const Double_t xsec_wmr_sys= sqrt(pow(0.0387496*blind,2)+th_wpr*th_wpr); 
  Double_t xsec_wz = 10.6556;  blind=theo_wz/xsec_wz; xsec_wz=xsec_wz*blind;  const Double_t xsec_wz_stat = 0.0;   const Double_t xsec_wz_sys = sqrt(pow(0.080893*blind,2)+th_wz*th_wz); ; 

  const Double_t ratio_wp     = xsec_wp/theo_wp;   const Double_t ratio_wp_the = theo_wp_unc*xsec_wp/(theo_wp*theo_wp); 
  const Double_t ratio_wp_exp = sqrt(xsec_wp_stat*xsec_wp_stat+xsec_wp_sys*xsec_wp_sys)/(theo_wp);

  const Double_t ratio_wm     = xsec_wm/theo_wm;   const Double_t ratio_wm_the = theo_wm_unc*xsec_wm/(theo_wm*theo_wm); 
  const Double_t ratio_wm_exp = sqrt(xsec_wm_stat*xsec_wm_stat+xsec_wm_sys*xsec_wm_sys)/(theo_wm);

  const Double_t ratio_w      = xsec_w/theo_w;     const Double_t ratio_w_the  = theo_w_unc*xsec_w/(theo_w*theo_w); 
  const Double_t ratio_w_exp  = sqrt(xsec_w_stat*xsec_w_stat+xsec_w_sys*xsec_w_sys)/(theo_w);

  const Double_t ratio_z      = xsec_z/theo_z;     const Double_t ratio_z_the  = theo_z_unc*xsec_z/(theo_z*theo_z); 
  const Double_t ratio_z_exp  = sqrt(xsec_z_stat*xsec_z_stat+xsec_z_sys*xsec_z_sys)/(theo_z);

  const Double_t ratio_wr     = xsec_wr/theo_wr;   const Double_t ratio_wr_the = theo_wr_unc*xsec_wr/(theo_wr*theo_wr); 
  const Double_t ratio_wr_exp = sqrt(xsec_wr_stat*xsec_wr_stat+xsec_wr_sys*xsec_wr_sys)/(theo_wr);

  const Double_t ratio_wpr    = xsec_wpr/theo_wpr; const Double_t ratio_wpr_the= theo_wpr_unc*xsec_wpr/(theo_wpr*theo_wpr); 
  const Double_t ratio_wpr_exp= sqrt(xsec_wpr_stat*xsec_wpr_stat+xsec_wpr_sys*xsec_wpr_sys)/(theo_wpr);

  const Double_t ratio_wmr    = xsec_wmr/theo_wmr;   const Double_t ratio_wmr_the = theo_wmr_unc*xsec_wmr/(theo_wmr*theo_wmr); 
  const Double_t ratio_wmr_exp= sqrt(xsec_wmr_stat*xsec_wmr_stat+xsec_wmr_sys*xsec_wmr_sys)/(theo_wmr);

  const Double_t ratio_wz     = xsec_wz/theo_wz;   const Double_t ratio_wz_the = theo_wz_unc*xsec_wz/(theo_wz*theo_wz); 
  const Double_t ratio_wz_exp = sqrt(xsec_wz_stat*xsec_wz_stat+xsec_wz_sys*xsec_wz_sys)/(theo_wz);

  Double_t test = sqrt(xsec_wp_stat*xsec_wp_stat+xsec_wp_sys*xsec_wp_sys)/(xsec_wp);

  //==============================================================================================================  
  
  //--------------------------------------------------------------------------------------------------------------
  // plotting parameter
  Double_t yshift, ydrift; // used for point and error bar
  yshift = 3.75;
  ydrift = -0.245;

  Double_t td, tp, tw; // used for text 
  td = 0.025;
  tp = 0.708;
  tw = 0.035;

  Double_t range_max, range_min;
  //range_min = 0.82;
  //range_max = 1.27;
  range_min = 0.82;
  range_max = 1.25;

  Int_t expColor, theColor;  // colors used
  Int_t prodColor, decayColor, cmbColor;  // colors used
  expColor  = kBlue;
  theColor  = kGreen+2;
  prodColor = kBlue+2;
  decayColor= kRed+2;
  cmbColor  = kBlack;
  //==============================================================================================================  

  //--------------------------------------------------------------------------------------------------------------
  // make plots
  TCanvas *c = MakeCanvas("c","c",800,600);
  c->SetTickx(1);
  c->SetTicky(0);
  c->SetFrameFillStyle(0);
  c->SetFrameLineWidth(2);
  c->SetFrameBorderMode(0);  
  c->SetLeftMargin(0.07);
  
  gStyle->SetEndErrorSize(8);

  Double_t xval, errl, errh, yval;
  xval = ratio_wp;
  yval = yshift+0*ydrift;
  errl = ratio_wp_exp;
  errh = sqrt(ratio_wp_exp*ratio_wp_exp+ratio_wp_the*ratio_wp_the);
  TGraphAsymmErrors grWP(1,&xval,&yval,&errl,&errl,0,0);
 
  grWP.SetTitle("");
  grWP.GetXaxis()->SetTitle("");
  grWP.GetXaxis()->SetTitleSize(0.05);
  grWP.GetYaxis()->SetTitle("");
  grWP.GetYaxis()->SetRangeUser(0,5);
  grWP.GetXaxis()->SetLimits(range_min,range_max);
  grWP.GetXaxis()->SetNdivisions(506);
  grWP.GetYaxis()->SetNdivisions(0);
  grWP.SetMarkerStyle(kFullCircle);
  grWP.SetMarkerSize(1);
  grWP.SetLineWidth(2);
  grWP.SetMarkerColor(kBlack);
  grWP.SetLineColor(expColor);
  grWP.Draw("AP");

  // lumi uncertainty band
  TBox lumi_box(0.983,yshift+7.*ydrift,1.017,3.98);
  //TBox lumi_box(0.974,yshift+7.*ydrift,1.026,3.98);
  lumi_box.SetLineColor(796);
  lumi_box.SetFillColor(796);
  lumi_box.Draw();
  c->RedrawAxis();

  TLine theory_line(1,0.02,1,3.98);
  theory_line.SetLineColor(kRed);
  theory_line.SetLineStyle(1);
  theory_line.SetLineWidth(3);
  theory_line.Draw();

  TGraphAsymmErrors grWP2(1,&xval,&yval,&errh,&errh,0,0);
  grWP2.SetMarkerStyle(kFullCircle);
  grWP2.SetMarkerSize(1);
  grWP2.SetLineWidth(2);
  grWP2.SetMarkerColor(kBlack);
  grWP2.SetLineColor(theColor);
  grWP2.Draw("EPSAME");
  grWP.Draw("EPSAME");
 
  xval = ratio_wm;
  yval = yshift+2*ydrift;
  errl = ratio_wm_exp;
  errh = sqrt(ratio_wm_exp*ratio_wm_exp+ratio_wm_the*ratio_wm_the);
  TGraphAsymmErrors grWM(1,&xval,&yval,&errl,&errl,0,0);
  grWM.SetMarkerStyle(kFullCircle);
  grWM.SetMarkerSize(1);
  grWM.SetLineWidth(2);
  grWM.SetMarkerColor(kBlack);
  grWM.SetLineColor(expColor);
  TGraphAsymmErrors grWM2(1,&xval,&yval,&errh,&errh,0,0);
  grWM2.SetMarkerStyle(kFullCircle);
  grWM2.SetMarkerSize(1);
  grWM2.SetLineWidth(2);
  grWM2.SetMarkerColor(kBlack);
  grWM2.SetLineColor(theColor);
  grWM2.Draw("EPSAME");
  grWM.Draw("EPSAME");

  xval = ratio_w;
  yval = yshift+4*ydrift;
  errl = ratio_w_exp;
  errh = sqrt(ratio_w_exp*ratio_w_exp+ratio_w_the*ratio_w_the);
  TGraphAsymmErrors grW(1,&xval,&yval,&errl,&errl,0,0);
  grW.SetMarkerStyle(kFullCircle);
  grW.SetMarkerSize(1);
  grW.SetLineWidth(2);
  grW.SetMarkerColor(kBlack);
  grW.SetLineColor(expColor);
  TGraphAsymmErrors grW2(1,&xval,&yval,&errh,&errh,0,0);
  grW2.SetMarkerStyle(kFullCircle);
  grW2.SetMarkerSize(1);
  grW2.SetLineWidth(2);
  grW2.SetMarkerColor(kBlack);
  grW2.SetLineColor(theColor);
  grW2.Draw("EPSAME");
  grW.Draw("EPSAME");

  xval = ratio_z;
  yval = yshift+6*ydrift;
  errl = ratio_z_exp;
  errh = sqrt(ratio_z_exp*ratio_z_exp+ratio_z_the*ratio_z_the);
  TGraphAsymmErrors grZ(1,&xval,&yval,&errl,&errl,0,0);
  grZ.SetMarkerSize(1);
  grZ.SetLineWidth(2);
  grZ.SetMarkerColor(kBlack);
  grZ.SetLineColor(expColor);
  TGraphAsymmErrors grZ2(1,&xval,&yval,&errh,&errh,0,0);
  grZ2.SetMarkerStyle(kFullCircle);
  grZ2.SetMarkerSize(1);
  grZ2.SetLineWidth(2);
  grZ2.SetMarkerColor(kBlack);
  grZ2.SetLineColor(theColor);
  grZ2.Draw("EPSAME");
  grZ.Draw("EPSAME");

  xval = ratio_wr;
  yval = yshift+8*ydrift;
  errl = ratio_wr_exp;
  errh = sqrt(ratio_wr_exp*ratio_wr_exp+ratio_wr_the*ratio_wr_the);
  TGraphAsymmErrors grWR(1,&xval,&yval,&errl,&errl,0,0);
  grWR.SetMarkerStyle(kFullCircle);
  grWR.SetMarkerSize(1);
  grWR.SetLineWidth(2);
  grWR.SetMarkerColor(kBlack);
  grWR.SetLineColor(expColor);
  TGraphAsymmErrors grWR2(1,&xval,&yval,&errh,&errh,0,0);
  grWR2.SetMarkerStyle(kFullCircle);
  grWR2.SetMarkerSize(1);
  grWR2.SetLineWidth(2);
  grWR2.SetMarkerColor(kBlack);
  grWR2.SetLineColor(theColor);
  grWR2.Draw("EPSAME");
  grWR.Draw("EPSAME");

  xval = ratio_wpr;
  yval = yshift+10*ydrift;
  errl = ratio_wpr_exp;
  errh = sqrt(ratio_wpr_exp*ratio_wpr_exp+ratio_wpr_the*ratio_wpr_the);
  TGraphAsymmErrors grWPR(1,&xval,&yval,&errl,&errl,0,0);
  grWPR.SetMarkerStyle(kFullCircle);
  grWPR.SetMarkerSize(1);
  grWPR.SetLineWidth(2);
  grWPR.SetMarkerColor(kBlack);
  grWPR.SetLineColor(expColor);
  TGraphAsymmErrors grWPR2(1,&xval,&yval,&errh,&errh,0,0);
  grWPR2.SetMarkerStyle(kFullCircle);
  grWPR2.SetMarkerSize(1);
  grWPR2.SetLineWidth(2);
  grWPR2.SetMarkerColor(kBlack);
  grWPR2.SetLineColor(theColor);
  grWPR2.Draw("EPSAME");
  grWPR.Draw("EPSAME");

  xval = ratio_wmr;
  yval = yshift+12*ydrift;
  errl = ratio_wmr_exp;
  errh = sqrt(ratio_wmr_exp*ratio_wmr_exp+ratio_wmr_the*ratio_wmr_the);
  TGraphAsymmErrors grWMR(1,&xval,&yval,&errl,&errl,0,0);
  grWMR.SetMarkerStyle(kFullCircle);
  grWMR.SetMarkerSize(1);
  grWMR.SetLineWidth(2);
  grWMR.SetMarkerColor(kBlack);
  grWMR.SetLineColor(expColor);
  TGraphAsymmErrors grWMR2(1,&xval,&yval,&errh,&errh,0,0);
  grWMR2.SetMarkerStyle(kFullCircle);
  grWMR2.SetMarkerSize(1);
  grWMR2.SetLineWidth(2);
  grWMR2.SetMarkerColor(kBlack);
  grWMR2.SetLineColor(theColor);
  grWMR2.Draw("EPSAME");
  grWMR.Draw("EPSAME");

  xval = ratio_wz;
  yval = yshift+14*ydrift;
  errl = ratio_wz_exp;
  errh = sqrt(ratio_wz_exp*ratio_wz_exp+ratio_wz_the*ratio_wz_the);
  TGraphAsymmErrors grWZ(1,&xval,&yval,&errl,&errl,0,0);
  grWZ.SetMarkerStyle(kFullCircle);
  grWZ.SetMarkerSize(1);
  grWZ.SetLineWidth(2);
  grWZ.SetMarkerColor(kBlack);
  grWZ.SetLineColor(expColor);
  TGraphAsymmErrors grWZ2(1,&xval,&yval,&errh,&errh,0,0);
  grWZ2.SetMarkerStyle(kFullCircle);
  grWZ2.SetMarkerSize(1);
  grWZ2.SetLineWidth(2);
  grWZ2.SetMarkerColor(kBlack);
  grWZ2.SetLineColor(theColor);
  grWZ2.Draw("EPSAME");
  grWZ.Draw("EPSAME");


  // legend
  xval = range_min+(0.06*(range_max-range_min));
  yval = 4.6;
  errl = 0.015*(range_max-range_min);
  errh = 0.03*(range_max-range_min);
  TGraphAsymmErrors grLeg(1,&xval,&yval,&errl,&errl,0,0);
  grLeg.SetMarkerStyle(kFullCircle);
  grLeg.SetMarkerSize(1);
  grLeg.SetLineWidth(2);
  grLeg.SetMarkerColor(kBlack);
  grLeg.SetLineColor(expColor);
  TGraphAsymmErrors grLeg2(1,&xval,&yval,&errh,&errh,0,0);
  grLeg2.SetMarkerStyle(kFullCircle);
  grLeg2.SetMarkerSize(1);
  grLeg2.SetLineWidth(2);
  grLeg2.SetMarkerColor(kBlack);
  grLeg2.SetLineColor(theColor);
  grLeg2.Draw("EPSAME");
  grLeg.Draw("EPSAME");

  TPaveText tb4(0.14,0.84,0.6,0.875,"NDC");
  tb4.SetFillStyle(0);
  tb4.SetBorderSize(0);
  tb4.SetTextAlign(12);
  tb4.AddText("Observation, uncertainty (exp., exp. #oplus theory)");
  tb4.Draw(); 

  TPaveText tb6(0.14,0.795,0.6,0.83,"NDC");
  tb6.SetFillStyle(0);
  tb6.SetBorderSize(0);
  tb6.SetTextAlign(12);
  tb6.AddText("Uncertainty (lumi)");
  tb6.Draw(); 

  TBox exp_box_leg2(range_min+(0.04*(range_max-range_min)),4.18,range_min+(0.08*(range_max-range_min)),4.38);
  exp_box_leg2.SetLineColor(796);
  exp_box_leg2.SetFillColor(796);
  exp_box_leg2.Draw();
  
  // split lines

  TLine split_line0(range_min,yshift+1*ydrift,range_max,yshift+1*ydrift);
  split_line0.SetLineColor(kBlack);
  split_line0.SetLineStyle(2);
  split_line0.SetLineWidth(1);
  split_line0.Draw();

  TLine split_line1(range_min,yshift+3*ydrift,range_max,yshift+3*ydrift);
  split_line1.SetLineColor(kBlack);
  split_line1.SetLineStyle(2);
  split_line1.SetLineWidth(1);
  split_line1.Draw();

  TLine split_line2(range_min,yshift+5*ydrift,range_max,yshift+5*ydrift);
  split_line2.SetLineColor(kBlack);
  split_line2.SetLineStyle(2);
  split_line2.SetLineWidth(1);
  split_line2.Draw();

  TLine split_line3(range_min,yshift+7*ydrift,range_max,yshift+7*ydrift);
  split_line3.SetLineColor(kBlack);
  split_line3.SetLineStyle(2);
  split_line3.SetLineWidth(1);
  split_line3.Draw();

  TLine split_line4(range_min,yshift+9*ydrift,range_max,yshift+9*ydrift);
  split_line4.SetLineColor(kBlack);
  split_line4.SetLineStyle(2);
  split_line4.SetLineWidth(1);
  split_line4.Draw();
  
  TLine split_line5(range_min,yshift+11*ydrift,range_max,yshift+11*ydrift);
  split_line5.SetLineColor(kBlack);
  split_line5.SetLineStyle(2);
  split_line5.SetLineWidth(1);
  split_line5.Draw();
  
  TLine split_line6(range_min,yshift+13*ydrift,range_max,yshift+13*ydrift);
  split_line6.SetLineColor(kBlack);
  split_line6.SetLineStyle(2);
  split_line6.SetLineWidth(1);
  split_line6.Draw();
  
  TLine split_line7(range_min,yshift+(-1)*ydrift,range_max,yshift+(-1)*ydrift);
  split_line7.SetLineColor(kBlack);
  split_line7.SetLineStyle(2);
  split_line7.SetLineWidth(1);
  split_line7.Draw();



  TPaveText tb1(0.08,0.93,0.34,0.99,"NDC");
  tb1.SetFillStyle(0);
  tb1.SetBorderSize(0);
  tb1.SetTextAlign(12);
  tb1.AddText("#bf{CMS}");
  tb1.Draw();

  TPaveText tb2(0.75,0.93,0.95,0.99,"NDC");
  tb2.SetFillStyle(0);
  tb2.SetBorderSize(0);
  tb2.SetTextAlign(12);
  tb2.AddText("200 pb^{-1} (13 TeV)");
  tb2.Draw(); 

  TPaveText tb3(0.6,0.83,0.95,0.88,"NDC");
  tb3.SetFillStyle(0);
  tb3.SetBorderSize(0);
  tb3.SetTextAlign(12);
  tb3.AddText("Theory: FEWZ (NNLO), NNPDF3.1");
  tb3.Draw(); 

  TPaveText tb7(0.6,0.795,0.95,0.83,"NDC");
  tb7.SetFillStyle(0);
  tb7.SetBorderSize(0);
  tb7.SetTextAlign(12);
  tb7.AddText("Observation: NNPDF3.1");
  tb7.Draw(); 

  TPaveText textwp(0.10,tp-(0*td),0.50,tp+tw-(0*td),"NDC");
  textwp.SetFillStyle(0);
  textwp.SetBorderSize(0);
  textwp.SetTextAlign(12);
  textwp.AddText("#bf{W^{+}#rightarrowl^{+}#nu}");
  textwp.Draw(); 

  char buffer[200]; 
  Double_t tmp[4];
  Double_t tmp2[4];
  roundXsec(xsec_wp, xsec_wp_stat, xsec_wp_sys, xsec_wp_lumi, tmp);
  sprintf(buffer,"%.0f #pm %.0f_{syst}", tmp[0], tmp[2]);
  sprintf(buffer,"%s #pm %.0f_{lum} pb",buffer,tmp[3]);
  TPaveText resultwp(0.60,tp-(-0.5*td),1.00,tp+tw-(-0.5*td),"NDC");
  resultwp.SetFillStyle(0);
  resultwp.SetBorderSize(0);
  resultwp.SetTextAlign(12);
  resultwp.AddText(buffer);
  resultwp.Draw(); 

  roundXsec(theo_wp, theo_wp_unc, tmp2);
  sprintf(buffer,"%.0f #pm %.0f pb", tmp2[0], tmp2[1]);
  TPaveText theorywp(0.60,tp-(0.6*td),1.00,tp+tw-(0.6*td),"NDC");
  theorywp.SetFillStyle(0);
  theorywp.SetBorderSize(0);
  theorywp.SetTextAlign(12);
  theorywp.AddText(buffer);
  theorywp.Draw(); 

  TPaveText textwm(0.10,tp-(3*td),0.50,tp+tw-(3*td),"NDC");
  textwm.SetFillStyle(0);
  textwm.SetBorderSize(0);
  textwm.SetTextAlign(12);
  textwm.AddText("#bf{W^{-}#rightarrowl^{-}#nu}");
  textwm.Draw(); 

  roundXsec(xsec_wm, xsec_wm_stat, xsec_wm_sys,xsec_wm_lumi, tmp);
  sprintf(buffer,"%.0f #pm %.0f_{syst}", tmp[0], tmp[2]);
  sprintf(buffer,"%s #pm %.0f_{lum} pb",buffer,tmp[3]);
  TPaveText resultwm(0.60,tp-(2.5*td),1.00,tp+tw-(2.5*td),"NDC");
  resultwm.SetFillStyle(0);
  resultwm.SetBorderSize(0);
  resultwm.SetTextAlign(12);
  resultwm.AddText(buffer);
  resultwm.Draw(); 

  roundXsec(theo_wm, theo_wm_unc, tmp2);
  sprintf(buffer,"%.0f #pm %.0f pb", tmp2[0], tmp2[1]);
  TPaveText theorywm(0.60,tp-(3.6*td),1.00,tp+tw-(3.6*td),"NDC");
  theorywm.SetFillStyle(0);
  theorywm.SetBorderSize(0);
  theorywm.SetTextAlign(12);
  theorywm.AddText(buffer);
  theorywm.Draw(); 

  TPaveText textw(0.10,tp-(6*td),0.50,tp+tw-(6*td),"NDC");
  textw.SetFillStyle(0);
  textw.SetBorderSize(0);
  textw.SetTextAlign(12);
  textw.AddText("#bf{W#rightarrowl#nu}");
  textw.Draw(); 

  roundXsec(xsec_w, xsec_w_stat, xsec_w_sys, xsec_w_lumi, tmp);
  sprintf(buffer,"%.0f #pm %.0f_{syst}", tmp[0], tmp[2]);
  sprintf(buffer,"%s #pm %.0f_{lum} pb",buffer,tmp[3]);
  TPaveText resultw(0.60,tp-(5.5*td),1.00,tp+tw-(5.5*td),"NDC");
  resultw.SetFillStyle(0);
  resultw.SetBorderSize(0);
  resultw.SetTextAlign(12);
  resultw.AddText(buffer);
  resultw.Draw(); 

  roundXsec(theo_w, theo_w_unc, tmp2);
  sprintf(buffer,"%.0f #pm %.0f pb", tmp2[0], tmp2[1]);
  TPaveText theoryw(0.60,tp-(6.6*td),1.00,tp+tw-(6.6*td),"NDC");
  theoryw.SetFillStyle(0);
  theoryw.SetBorderSize(0);
  theoryw.SetTextAlign(12);
  theoryw.AddText(buffer);
  theoryw.Draw(); 

  TPaveText textz(0.10,tp-(9*td),0.50,tp+tw-(9*td),"NDC");
  textz.SetFillStyle(0);
  textz.SetBorderSize(0);
  textz.SetTextAlign(12);
  textz.AddText("#bf{Z#rightarrowl^{+}l^{-}}");
  textz.Draw(); 

  roundXsec(xsec_z, xsec_z_stat, xsec_z_sys,xsec_z_lumi, tmp);
  sprintf(buffer,"%.0f #pm %.0f_{syst}", tmp[0], tmp[2]);
  sprintf(buffer,"%s #pm %.0f_{lum} pb",buffer,tmp[3]);
  TPaveText resultz(0.60,tp-(8.5*td),1.00,tp+tw-(8.5*td),"NDC");
  resultz.SetFillStyle(0);
  resultz.SetBorderSize(0);
  resultz.SetTextAlign(12);
  resultz.AddText(buffer);
  resultz.Draw(); 

  roundXsec(theo_z, theo_z_unc, tmp2);
  sprintf(buffer,"%.0f #pm %.0f pb", tmp2[0], tmp2[1]);
  TPaveText theoryz(0.60,tp-(9.6*td),1.00,tp+tw-(9.6*td),"NDC");
  theoryz.SetFillStyle(0);
  theoryz.SetBorderSize(0);
  theoryz.SetTextAlign(12);
  theoryz.AddText(buffer);
  theoryz.Draw(); 

  TPaveText textwr(0.10,tp-(12*td),0.50,tp+tw-(12*td),"NDC");
  textwr.SetFillStyle(0);
  textwr.SetBorderSize(0);
  textwr.SetTextAlign(12);
  textwr.AddText("#bf{W^{+}#rightarrowl^{+}#nu / W^{-}#rightarrowl^{-}#nu}");
  textwr.Draw(); 

  sprintf(buffer,"%.3f #pm %.3f_{syst}", xsec_wr, xsec_wr_sys);
  TPaveText resultwr(0.60,tp-(11.5*td),1.00,tp+tw-(11.5*td),"NDC");
  resultwr.SetFillStyle(0);
  resultwr.SetBorderSize(0);
  resultwr.SetTextAlign(12);
  resultwr.AddText(buffer);
  resultwr.Draw(); 

  sprintf(buffer,"%.3f #pm %.3f", theo_wr, theo_wr_unc);
  TPaveText theorywr(0.60,tp-(12.6*td),1.00,tp+tw-(12.6*td),"NDC");
  theorywr.SetFillStyle(0);
  theorywr.SetBorderSize(0);
  theorywr.SetTextAlign(12);
  theorywr.AddText(buffer);
  theorywr.Draw(); 

  TPaveText textwpr(0.10,tp-(15*td),0.50,tp+tw-(15*td),"NDC");
  textwpr.SetFillStyle(0);
  textwpr.SetBorderSize(0);
  textwpr.SetTextAlign(12);
  textwpr.AddText("#bf{W^{+}#rightarrowl^{+}#nu / Z#rightarrowl^{+}l^{-}}");
  textwpr.Draw(); 

  sprintf(buffer,"%.2f #pm %.2f_{syst}", xsec_wpr, xsec_wpr_sys);
  TPaveText resultwpr(0.60,tp-(14.5*td),1.00,tp+tw-(14.5*td),"NDC");
  resultwpr.SetFillStyle(0);
  resultwpr.SetBorderSize(0);
  resultwpr.SetTextAlign(12);
  resultwpr.AddText(buffer);
  resultwpr.Draw(); 

  sprintf(buffer,"%.2f #pm %.2f", theo_wpr, theo_wpr_unc);
  TPaveText theorywpr(0.60,tp-(15.6*td),1.00,tp+tw-(15.6*td),"NDC");
  theorywpr.SetFillStyle(0);
  theorywpr.SetBorderSize(0);
  theorywpr.SetTextAlign(12);
  theorywpr.AddText(buffer);
  theorywpr.Draw(); 

  TPaveText textwmr(0.10,tp-(18*td),0.50,tp+tw-(18*td),"NDC");
  textwmr.SetFillStyle(0);
  textwmr.SetBorderSize(0);
  textwmr.SetTextAlign(12);
  textwmr.AddText("#bf{W^{-}#rightarrowl^{-}#nu / Z#rightarrowl^{+}l^{-}}");
  textwmr.Draw(); 

  sprintf(buffer,"%.2f #pm %.2f_{syst}", xsec_wmr, xsec_wmr_sys);
  TPaveText resultwmr(0.60,tp-(17.5*td),1.00,tp+tw-(17.5*td),"NDC");
  resultwmr.SetFillStyle(0);
  resultwmr.SetBorderSize(0);
  resultwmr.SetTextAlign(12);
  resultwmr.AddText(buffer);
  resultwmr.Draw(); 

  sprintf(buffer,"%.2f #pm %.2f", theo_wmr, theo_wmr_unc);
  TPaveText theorywmr(0.60,tp-(18.6*td),1.00,tp+tw-(18.6*td),"NDC");
  theorywmr.SetFillStyle(0);
  theorywmr.SetBorderSize(0);
  theorywmr.SetTextAlign(12);
  theorywmr.AddText(buffer);
  theorywmr.Draw(); 

  TPaveText textwz(0.10,tp-(21*td),0.50,tp+tw-(21*td),"NDC");
  textwz.SetFillStyle(0);
  textwz.SetBorderSize(0);
  textwz.SetTextAlign(12);
  textwz.AddText("#bf{W#rightarrowl#nu / Z#rightarrowl^{+}l^{-}}");
  textwz.Draw(); 

  sprintf(buffer,"%.2f #pm  %.2f_{syst}", xsec_wz,  xsec_wz_sys);
  TPaveText resultwz(0.60,tp-(20.5*td),1.00,tp+tw-(20.5*td),"NDC");
  resultwz.SetFillStyle(0);
  resultwz.SetBorderSize(0);
  resultwz.SetTextAlign(12);
  resultwz.AddText(buffer);
  resultwz.Draw(); 

  sprintf(buffer,"%.2f #pm %.2f", theo_wz, theo_wz_unc);
  TPaveText theorywz(0.60,tp-(21.6*td),1.00,tp+tw-(21.6*td),"NDC");
  theorywz.SetFillStyle(0);
  theorywz.SetBorderSize(0);
  theorywz.SetTextAlign(12);
  theorywz.AddText(buffer);
  theorywz.Draw(); 

  TPaveText tb5(0.28,0.02,0.98,0.12,"NDC");
  tb5.SetFillStyle(0);
  tb5.SetBorderSize(0);
  tb5.SetTextAlign(12);
  tb5.AddText("ratio (exp./th.) of total cross sections and ratios");
  tb5.Draw();

  c->SaveAs("xsecSummary13TeV_ele.png");
  c->SaveAs("xsecSummary13TeV_ele.pdf");
}
