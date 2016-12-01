//================================================================================================
//
// Compute cross sections and produce summary plots
//
//________________________________________________________________________________________________

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "MitStyleRemix.hh"         // style settings for drawing
#include "TMath.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
// root -q -L xsecPDFSummary.C
void roundXsec(Double_t &x, Double_t &y, Double_t &z)
{
  Double_t xmin = TMath::Min(TMath::Min(x,y),z);
  int iPres = 1;
  if (xmin > 10) iPres = 10;
  //if (xmin > 100) iPres = 100;

  cout << xmin << endl;

  x = TMath::Nint(x/Double_t(iPres))*iPres;
  y = TMath::Nint(y/Double_t(iPres))*iPres;
  z = TMath::Nint(z/Double_t(iPres))*iPres;
}

void xsecPDFSummary(int itype = 1)
{   
  //-------------------------------------------------------------------------------------------------------------
  //const int itype = 3;  // defined which measurement to plot
  bool lumi     = false;
  bool electron = false;
  bool muon     = false;

  // input parameter
  Double_t xsec[24];
  Double_t xsec_stat[24];
  Double_t xsec_sys[24];
  Double_t xsec_lumi[24];

  Double_t plot_xsec[4];

/*  // total cross sections
  xsec[0] = 11297;   xsec_stat[0] = 8;     xsec_sys[0] = 191;     xsec_lumi[0] = 305;  // wp - 0
  xsec[1] = 8279;    xsec_stat[1] = 8;     xsec_sys[1] = 122;     xsec_lumi[1] = 224;  // wm - 1
  xsec[2] = 19638;   xsec_stat[2] = 11;     xsec_sys[2] = 267;     xsec_lumi[2] = 530;  // w  - 2
  xsec[3] = 1897;    xsec_stat[3] = 1;     xsec_sys[3] = 30;      xsec_lumi[3] = 51;   // z  - 3
  xsec[4] = 1.365;   xsec_stat[4] = 0.002;  xsec_sys[4] = 0.022;   xsec_lumi[4] = 0;     // wr - 4
  xsec[5] = 5.960;   xsec_stat[5] = 0.006;  xsec_sys[5] = 0.098;   xsec_lumi[5] = 0;     // wpr- 5
  xsec[6] = 4.366;   xsec_stat[6] = 0.005;  xsec_sys[6] = 0.068;   xsec_lumi[6] = 0;     // wmr- 6
  xsec[7] = 10.335;  xsec_stat[7] = 0.009;  xsec_sys[7] = 0.147;   xsec_lumi[7] = 0;     // wz - 7
*/


xsec[0]   = 11356.357;  xsec_stat[0] =   9.186;  xsec_sys[0] = 197.279;   xsec_lumi[0] = 306.622;
xsec[1]   =  8284.005;  xsec_stat[1] =   7.887;  xsec_sys[1] = 123.894;   xsec_lumi[1] = 223.668; 
xsec[2]   = 19708.254;  xsec_stat[2] =  11.055;  xsec_sys[2] = 265.308;   xsec_lumi[2] = 532.123; 
xsec[3]   =  1904.184;  xsec_stat[3] =   1.301;  xsec_sys[3] =  30.401;   xsec_lumi[3] =  51.413; 
xsec[4]   =     1.372;  xsec_stat[4] =   0.002;  xsec_sys[4] =   0.022;   xsec_lumi[4] =   0.000; 
xsec[5]   =     5.966;  xsec_stat[5] =   0.006;  xsec_sys[5] =   0.100;   xsec_lumi[5] =   0.000; 
xsec[6]   =     4.351;  xsec_stat[6] =   0.005;  xsec_sys[6] =   0.069;   xsec_lumi[6] =   0.000; 
xsec[7]   =    10.330;  xsec_stat[7] =   0.009;  xsec_sys[7] =   0.147;   xsec_lumi[7] =   0.000; 

// preapproval
// // fiducial muon
xsec[8]   =  5097.780;   xsec_stat[8]  =   3.451;   xsec_sys[8]  =  47.816;   xsec_lumi[8]  =  137.640;
xsec[9]   =  3865.767;   xsec_stat[9]  =   3.936;   xsec_sys[9]  =  37.902;   xsec_lumi[9]  =  104.376;
xsec[10]  =  8963.546;   xsec_stat[10] =   5.243;   xsec_sys[10] =  80.152;   xsec_lumi[10] =  242.016; 
xsec[11]  =   690.960;   xsec_stat[11] =   0.599;   xsec_sys[11] =   8.748;   xsec_lumi[11] =   18.656;
xsec[12]  =     1.319;   xsec_stat[12] =   0.002;   xsec_sys[12] =   0.009;   xsec_lumi[12] =      0.0;
xsec[13]  =     7.378;   xsec_stat[13] =   0.008;   xsec_sys[13] =   0.052;   xsec_lumi[13] =      0.0;
xsec[14]  =     5.595;   xsec_stat[14] =   0.007;   xsec_sys[14] =   0.043;   xsec_lumi[14] =      0.0;
xsec[15]  =    12.973;   xsec_stat[15] =   0.014;   xsec_sys[15] =   0.085;   xsec_lumi[15] =      0.0;
//fiducial electron
xsec[16]  =  4803.922;   xsec_stat[16] =   8.733;   xsec_sys[16] = 113.877;   xsec_lumi[16] = 129.706;
xsec[17]  =  3563.283;   xsec_stat[17] =   6.477;   xsec_sys[17] =  68.823;   xsec_lumi[17] =  96.209;
xsec[18]  =  8367.205;   xsec_stat[18] =  10.858;   xsec_sys[18] = 176.191;   xsec_lumi[18] = 225.915;
xsec[19]  =   634.777;   xsec_stat[19] =   0.702;   xsec_sys[19] =  13.744;   xsec_lumi[19] =  17.139;
xsec[20]  =     1.348;   xsec_stat[20] =   0.003;   xsec_sys[20] =   0.016;   xsec_lumi[20] =     0.0;
xsec[21]  =     7.568;   xsec_stat[21] =   0.016;   xsec_sys[21] =   0.151;   xsec_lumi[21] =     0.0;
xsec[22]  =     5.613;   xsec_stat[22] =   0.012;   xsec_sys[22] =   0.092;   xsec_lumi[22] =     0.0;
xsec[23]  =    13.181;   xsec_stat[23] =   0.022;   xsec_sys[23] =   0.232;   xsec_lumi[23] =     0.0;


// from xinmei's update
/*  // fiducial muon cross sections
  xsec[8]  = 5086;   xsec_stat[8] = 4;     xsec_sys[8] = 53;      xsec_lumi[8] = 137; // wp - 8
  xsec[9]  = 3861;   xsec_stat[9] = 4;    xsec_sys[9] = 39 ;      xsec_lumi[9] = 104; // wm - 9
  xsec[10] = 8948;  xsec_stat[10] = 5;   xsec_sys[10] = 92;      xsec_lumi[10] = 242; // w  - 10
  xsec[11] = 689;  xsec_stat[11] = 1;   xsec_sys[11] = 9 ;     xsec_lumi[11] = 19; // z  - 11
  xsec[12] = 1.317;  xsec_stat[12] = 0.002; xsec_sys[12] = 0.005;   xsec_lumi[12] = 0;     // wr - 12
  xsec[13] = 7.381;  xsec_stat[13] = 0.008; xsec_sys[13] = 0.061;   xsec_lumi[13] = 0;     // wpr- 13
  xsec[14] = 5.604;  xsec_stat[14] = 0.007;   xsec_sys[14] = 0.046;   xsec_lumi[14] = 0;     // wmr- 14
  xsec[15] = 12.985;  xsec_stat[15] = 0.014;  xsec_sys[15] = 0.107;    xsec_lumi[15] = 0;     // wz - 15

  // fiducial electron cross sections
  xsec[16] = 4780;  xsec_stat[16] = 6;   xsec_sys[16] = 98;     xsec_lumi[16] = 129; // wp - 8
  xsec[17] = 3575;  xsec_stat[17] = 6;   xsec_sys[17] = 63;      xsec_lumi[17] = 97; // wm - 9
  xsec[18] = 8355;  xsec_stat[18] = 9;    xsec_sys[18] = 158;      xsec_lumi[18] = 226; // w  - 10
  xsec[19] = 631;  xsec_stat[19] = 1;   xsec_sys[19] = 14;     xsec_lumi[19] = 17; // z  - 11
  xsec[20] = 1.337;  xsec_stat[20] = 0.003; xsec_sys[20] = 0.016;   xsec_lumi[20] = 0;     // wr - 12
  xsec[21] = 7.570;  xsec_stat[21] = 0.013;  xsec_sys[21] = 0.127;   xsec_lumi[21] = 0;     // wpr- 13
  xsec[22] = 5.662;  xsec_stat[22] = 0.012; xsec_sys[22] = 0.083;   xsec_lumi[22] = 0;     // wmr- 14
  xsec[23] = 13.232;  xsec_stat[23] = 0.02; xsec_sys[23] = 0.205;   xsec_lumi[23] = 0;     // wz - 15
*/
  plot_xsec[0] = xsec[itype];
  plot_xsec[1] = xsec_stat[itype];
  plot_xsec[2] = xsec_sys[itype];
  plot_xsec[3] = xsec_lumi[itype];

  Double_t nnpdf_xsec[24];
  Double_t nnpdf_pdf_p[24];
  Double_t nnpdf_pdf_m[24];
  Double_t nnpdf_scale_p[24];
  Double_t nnpdf_scale_m[24];
  Double_t plot_nnpdf[5];

  nnpdf_xsec[0]= 11328.8;  nnpdf_pdf_p[0] = 305.653; nnpdf_pdf_m[0] = 261.311; nnpdf_scale_p[0] = 323.940; nnpdf_scale_m[0] = 268.752;  
  nnpdf_xsec[1]= 8369.09;  nnpdf_pdf_p[1] = 227.035; nnpdf_pdf_m[1] = 207.479; nnpdf_scale_p[1] = 244.286; nnpdf_scale_m[1] = 213.113;   
  nnpdf_xsec[2]= 19697.6;  nnpdf_pdf_p[2] = 527.621; nnpdf_pdf_m[2] = 460.716; nnpdf_scale_p[2] = 563.479; nnpdf_scale_m[2] = 473.946;   
  nnpdf_xsec[3]= 1867.66;  nnpdf_pdf_p[3] = 45.6562; nnpdf_pdf_m[3] = 42.377;  nnpdf_scale_p[3] = 47.3550; nnpdf_scale_m[3] = 43.4859; 
  nnpdf_xsec[4]= 1.35367;  nnpdf_pdf_p[4] = 0.01039; nnpdf_pdf_m[4] = 0.01229; nnpdf_scale_p[4] = 0.01054; nnpdf_scale_m[4] = 0.01230; 
  nnpdf_xsec[5]= 6.06494;  nnpdf_pdf_p[5] = 0.04001; nnpdf_pdf_m[5] = 0.04714; nnpdf_scale_p[5] = 0.04361; nnpdf_scale_m[5] = 0.04716; 
  nnpdf_xsec[6]= 4.48044;  nnpdf_pdf_p[6] = 0.01934; nnpdf_pdf_m[6] = 0.01767; nnpdf_scale_p[6] = 0.02685; nnpdf_scale_m[6] = 0.01779; 
  nnpdf_xsec[7]= 10.54538; nnpdf_pdf_p[7] = 0.05935; nnpdf_pdf_m[7] = 0.064818;  nnpdf_scale_p[7] = 0.07045706; nnpdf_scale_m[7] = 0.06494966;   

  nnpdf_xsec[8 ]= 5032.41;  nnpdf_pdf_p[8 ] = 137.775;  nnpdf_pdf_m[8 ] = 116.078;  nnpdf_scale_p[8 ] = 175.219;  nnpdf_scale_m[8 ] = 155.715;  
  nnpdf_xsec[9 ]= 3838.06;  nnpdf_pdf_p[9 ] = 104.118;  nnpdf_pdf_m[9 ] = 95.1500;  nnpdf_scale_p[9 ] = 161.613;  nnpdf_scale_m[9 ] = 116.483;   
  nnpdf_xsec[10]= 8870.35;  nnpdf_pdf_p[10] = 237.602;  nnpdf_pdf_m[10] = 207.473;  nnpdf_scale_p[10] = 351.531;  nnpdf_scale_m[10] = 243.280;   
  nnpdf_xsec[11]= 677.112;  nnpdf_pdf_p[11] = 16.5525;  nnpdf_pdf_m[11] = 15.3634;  nnpdf_scale_p[11] = 25.4884;  nnpdf_scale_m[11] = 25.4884; 
  nnpdf_xsec[12]= 1.31121;  nnpdf_pdf_p[12] = 0.01007;  nnpdf_pdf_m[12] = 0.01191;  nnpdf_scale_p[12] = 0.03362;  nnpdf_scale_m[12] = 0.03362; 
  nnpdf_xsec[13]= 7.43114;  nnpdf_pdf_p[13] = 0.04902;  nnpdf_pdf_m[13] = 0.05776;  nnpdf_scale_p[13] = 0.16733;  nnpdf_scale_m[13] = 0.16732; 
  nnpdf_xsec[14]= 5.66750;  nnpdf_pdf_p[14] = 0.02446;  nnpdf_pdf_m[14] = 0.02236;  nnpdf_scale_p[14] = 0.11440;  nnpdf_scale_m[14] = 0.11440; 
  nnpdf_xsec[15]= 13.0986;  nnpdf_pdf_p[15] = 0.07372;  nnpdf_pdf_m[15] = 0.08051;  nnpdf_scale_p[15] = 0.24345;  nnpdf_scale_m[15] = 0.24344;   

  nnpdf_xsec[16]= 4873.83;  nnpdf_pdf_p[16] = 131.497;  nnpdf_pdf_m[16] = 112.420;  nnpdf_scale_p[16] = 158.733;  nnpdf_scale_m[16] = 138.355;  
  nnpdf_xsec[17]= 3691.73;  nnpdf_pdf_p[17] = 100.149;  nnpdf_pdf_m[17] = 91.5219;  nnpdf_scale_p[17] = 152.004;  nnpdf_scale_m[17] = 107.208;   
  nnpdf_xsec[18]= 8565.44;  nnpdf_pdf_p[18] = 229.434;  nnpdf_pdf_m[18] = 200.341;  nnpdf_scale_p[18] = 341.621;  nnpdf_scale_m[18] = 238.048;   
  nnpdf_xsec[19]= 623.513;  nnpdf_pdf_p[19] = 15.2422;  nnpdf_pdf_m[19] = 14.1473;  nnpdf_scale_p[19] = 23.6945;  nnpdf_scale_m[19] = 17.6491; 
  nnpdf_xsec[20]= 1.32022;  nnpdf_pdf_p[20] = 0.01014;  nnpdf_pdf_m[20] = 0.01199;  nnpdf_scale_p[20] = 0.02912;  nnpdf_scale_m[20] = 0.02724; 
  nnpdf_xsec[21]= 7.81566;  nnpdf_pdf_p[21] = 0.05156;  nnpdf_pdf_m[21] = 0.06075;  nnpdf_scale_p[21] = 0.17256;  nnpdf_scale_m[21] = 0.16316; 
  nnpdf_xsec[22]= 5.92005;  nnpdf_pdf_p[22] = 0.02555;  nnpdf_pdf_m[22] = 0.02335;  nnpdf_scale_p[22] = 0.12022;  nnpdf_scale_m[22] = 0.11486; 
  nnpdf_xsec[23]= 13.7357;  nnpdf_pdf_p[23] = 0.07731;  nnpdf_pdf_m[23] = 0.08443;  nnpdf_scale_p[23] = 0.26121;  nnpdf_scale_m[23] = 0.24455;   

  plot_nnpdf[0] = nnpdf_xsec[itype];  // value
  plot_nnpdf[1] = nnpdf_pdf_p[itype];   // pdf uncertainty up
  plot_nnpdf[2] = nnpdf_pdf_m[itype];   // pdf uncertainty down
  plot_nnpdf[3] = nnpdf_scale_p[itype];  // total uncertainty up 
  plot_nnpdf[4] = nnpdf_scale_m[itype];  // total uncertainty down

  Double_t ct14_xsec[24];
  Double_t ct14_pdf_p[24];
  Double_t ct14_pdf_m[24];
  Double_t ct14_scale_p[24];
  Double_t ct14_scale_m[24];
  Double_t plot_ct14[5];

  ct14_xsec[0]= 11501.7;   ct14_pdf_p[0] = 316.582;  ct14_pdf_m[0] = 307.738;  ct14_scale_p[0] = 334.801;  ct14_scale_m[0] = 314.273;  	
  ct14_xsec[1]= 8520.06;   ct14_pdf_p[1] = 213.365;  ct14_pdf_m[1] = 236.613;  ct14_scale_p[1] = 232.274;  ct14_scale_m[1] = 241.747;  	
  ct14_xsec[2]= 20021.8;   ct14_pdf_p[2] = 523.149;  ct14_pdf_m[2] = 536.751;  ct14_scale_p[2] = 560.454;  ct14_scale_m[2] = 548.523;  	
  ct14_xsec[3]= 1897.28;   ct14_pdf_p[3] = 45.2047;  ct14_pdf_m[3] = 50.3253;  ct14_scale_p[3] = 46.9736;  ct14_scale_m[3] = 51.2927;	
  ct14_xsec[4]= 1.34996;   ct14_pdf_p[4] = 0.01340;  ct14_pdf_m[4] = 0.01365;  ct14_scale_p[4] = 0.01352;  ct14_scale_m[4] = 0.01366;	
  ct14_xsec[5]= 6.0622;    ct14_pdf_p[5] = 0.05542;  ct14_pdf_m[5] = 0.05876;  ct14_scale_p[5] = 0.05806;  ct14_scale_m[5] = 0.05878;	
  ct14_xsec[6]= 4.49067;   ct14_pdf_p[6] = 0.02316;  ct14_pdf_m[6] = 0.02929;  ct14_scale_p[6] = 0.02974;  ct14_scale_m[6] = 0.02936;	
  ct14_xsec[7]= 10.55287;  ct14_pdf_p[7] = 0.07858;  ct14_pdf_m[7] = 0.08805;  ct14_scale_p[7] = 0.08781;  ct14_scale_m[7] = 0.08813;     

  ct14_xsec[8 ]= 5089.50;   ct14_pdf_p[8 ] = 140.087;  ct14_pdf_m[8 ] = 136.174;  ct14_scale_p[8 ] = 179.363;  ct14_scale_m[8 ] = 171.936;  	
  ct14_xsec[9 ]= 3911.85;   ct14_pdf_p[9 ] = 97.9632;  ct14_pdf_m[9 ] = 108.636;  ct14_scale_p[9 ] = 166.928;  ct14_scale_m[9 ] = 128.421;  	
  ct14_xsec[10]= 9001.37;   ct14_pdf_p[10] = 235.196;  ct14_pdf_m[10] = 241.311;  ct14_scale_p[10] = 371.940;  ct14_scale_m[10] = 273.590;  	
  ct14_xsec[11]= 687.684;   ct14_pdf_p[11] = 16.3847;  ct14_pdf_m[11] = 18.2408;  ct14_scale_p[11] = 27.2987;  ct14_scale_m[11] = 21.3386;	
  ct14_xsec[12]= 1.30105;   ct14_pdf_p[12] = 0.01291;  ct14_pdf_m[12] = 0.01315;  ct14_scale_p[12] = 0.03483;  ct14_scale_m[12] = 0.03230;	
  ct14_xsec[13]= 7.40092;   ct14_pdf_p[13] = 0.06765;  ct14_pdf_m[13] = 0.07174;  ct14_scale_p[13] = 0.17833;  ct14_scale_m[13] = 0.16363;	
  ct14_xsec[14]= 5.68844;   ct14_pdf_p[14] = 0.02933;  ct14_pdf_m[14] = 0.03709;  ct14_scale_p[14] = 0.11964;  ct14_scale_m[14] = 0.11355;	
  ct14_xsec[15]= 13.0894;   ct14_pdf_p[15] = 0.09746;  ct14_pdf_m[15] = 0.10921;  ct14_scale_p[15] = 0.26238;  ct14_scale_m[15] = 0.23871;     

  ct14_xsec[16]= 4936.44;   ct14_pdf_p[16] = 135.874;  ct14_pdf_m[16] = 132.078;  ct14_scale_p[16] = 163.005;  ct14_scale_m[16] = 155.295;  	
  ct14_xsec[17]= 3763.55;   ct14_pdf_p[17] = 94.2494;  ct14_pdf_m[17] = 104.518;  ct14_scale_p[17] = 157.133;  ct14_scale_m[17] = 119.011;  	
  ct14_xsec[18]= 8700.00;   ct14_pdf_p[18] = 227.322;  ct14_pdf_m[18] = 233.232;  ct14_scale_p[18] = 361.606;  ct14_scale_m[18] = 267.303;  	
  ct14_xsec[19]= 633.413;   ct14_pdf_p[19] = 15.0917;  ct14_pdf_m[19] = 16.8012;  ct14_scale_p[19] = 25.3599;  ct14_scale_m[19] = 19.9297;	
  ct14_xsec[20]= 1.31165;   ct14_pdf_p[20] = 0.01302;  ct14_pdf_m[20] = 0.01326;  ct14_scale_p[20] = 0.03064;  ct14_scale_m[20] = 0.02768;	
  ct14_xsec[21]= 7.79339;   ct14_pdf_p[21] = 0.07124;  ct14_pdf_m[21] = 0.07124;  ct14_scale_p[21] = 0.18460;  ct14_scale_m[21] = 0.16883;	
  ct14_xsec[22]= 5.94170;   ct14_pdf_p[22] = 0.03064;  ct14_pdf_m[22] = 0.03875;  ct14_scale_p[22] = 0.12566;  ct14_scale_m[22] = 0.11934;	
  ct14_xsec[23]= 13.7351;   ct14_pdf_p[23] = 0.10227;  ct14_pdf_m[23] = 0.11460;  ct14_scale_p[23] = 0.28083;  ct14_scale_m[23] = 0.25652;     

  plot_ct14[0] = ct14_xsec[itype];  // value
  plot_ct14[1] = ct14_pdf_p[itype];   // pdf uncertainty up
  plot_ct14[2] = ct14_pdf_m[itype];   // pdf uncertainty down
  plot_ct14[3] = ct14_scale_p[itype];  // total uncertainty up 
  plot_ct14[4] = ct14_scale_m[itype];  // total uncertainty down

  Double_t mmht_xsec[24];
  Double_t mmht_pdf_p[24];
  Double_t mmht_pdf_m[24];
  Double_t mmht_scale_p[24];
  Double_t mmht_scale_m[24];
  Double_t plot_mmht[5];

  mmht_xsec[0]= 11578.1;   mmht_pdf_p[0] = 239.980;  mmht_pdf_m[0] = 200.904;  mmht_scale_p[0] = 263.848; mmht_scale_m[0] = 210.907;  
  mmht_xsec[1]= 8588.01;   mmht_pdf_p[1] = 161.253;  mmht_pdf_m[1] = 163.612;  mmht_scale_p[1] = 185.914; mmht_scale_m[1] = 171.068;  
  mmht_xsec[2]= 20166.8;   mmht_pdf_p[2] = 382.784;  mmht_pdf_m[2] = 371.484;  mmht_scale_p[2] = 433.052; mmht_scale_m[2] = 388.538;  
  mmht_xsec[3]= 1915.66;   mmht_pdf_p[3] = 36.5145;  mmht_pdf_m[3] = 35.2492;  mmht_scale_p[3] = 38.7239; mmht_scale_m[3] = 36.6432;
  mmht_xsec[4]= 1.34825;   mmht_pdf_p[4] = 0.01100;  mmht_pdf_m[4] = 0.00795;  mmht_scale_p[4] = 0.01114; mmht_scale_m[4] = 0.00796;
  mmht_xsec[5]= 6.04473;   mmht_pdf_p[5] = 0.04467;  mmht_pdf_m[5] = 0.05213;  mmht_scale_p[5] = 0.04790; mmht_scale_m[5] = 0.05214;
  mmht_xsec[6]= 4.48338;   mmht_pdf_p[6] = 0.02454;  mmht_pdf_m[6] = 0.03728;  mmht_scale_p[6] = 0.03082; mmht_scale_m[6] = 0.03733;
  mmht_xsec[7]= 10.52811;  mmht_pdf_p[7] = 0.06921;  mmht_pdf_m[7] = 0.08940;  mmht_scale_p[7] = 0.07871; mmht_scale_m[7] = 0.08948;  

  mmht_xsec[8 ]= 5117.34;   mmht_pdf_p[ 8] = 106.067;   mmht_pdf_m[ 8] = 88.7966;  mmht_scale_p[ 8] = 154.708;  mmht_scale_m[ 8] = 136.033;  
  mmht_xsec[9 ]= 3932.64;   mmht_pdf_p[ 9] = 73.8413;   mmht_pdf_m[ 9] = 74.9216;  mmht_scale_p[ 9] = 130.936;  mmht_scale_m[ 9] = 99.4816;  
  mmht_xsec[10]= 9050.30;   mmht_pdf_p[10] = 171.782;   mmht_pdf_m[10] = 166.711;  mmht_scale_p[10] = 283.464;  mmht_scale_m[10] = 206.356;  
  mmht_xsec[11]= 693.777;   mmht_pdf_p[11] = 13.2241;   mmht_pdf_m[11] = 12.7659;  mmht_scale_p[11] = 21.3940;  mmht_scale_m[11] = 16.1562;
  mmht_xsec[12]= 1.30132;   mmht_pdf_p[12] = 0.01061;   mmht_pdf_m[12] = 0.00767;  mmht_scale_p[12] = 0.03235;  mmht_scale_m[12] = 0.03051;
  mmht_xsec[13]= 7.37705;   mmht_pdf_p[13] = 0.05451;   mmht_pdf_m[13] = 0.06361;  mmht_scale_p[13] = 0.17144;  mmht_scale_m[13] = 0.16117;
  mmht_xsec[14]= 5.66886;   mmht_pdf_p[14] = 0.03103;   mmht_pdf_m[14] = 0.04713;  mmht_scale_p[14] = 0.12402;  mmht_scale_m[14] = 0.11774;
  mmht_xsec[15]= 13.0459;   mmht_pdf_p[15] = 0.08576;   mmht_pdf_m[15] = 0.11078;  mmht_scale_p[15] = 0.26035;  mmht_scale_m[15] = 0.24139;  

  mmht_xsec[16]= 4961.98;   mmht_pdf_p[16] = 102.847;   mmht_pdf_m[16] = 86.1007;  mmht_scale_p[16] = 137.006;  mmht_scale_m[16] = 118.972;  
  mmht_xsec[17]= 3783.49;   mmht_pdf_p[17] = 71.0408;   mmht_pdf_m[17] = 72.0801;  mmht_scale_p[17] = 123.199;  mmht_scale_m[17] = 92.0309;  
  mmht_xsec[18]= 8745.77;   mmht_pdf_p[18] = 166.002;   mmht_pdf_m[18] = 161.102;  mmht_scale_p[18] = 280.102;  mmht_scale_m[18] = 207.816;  
  mmht_xsec[19]= 638.943;   mmht_pdf_p[19] = 12.1789;   mmht_pdf_m[19] = 11.7569;  mmht_scale_p[19] = 20.5419;  mmht_scale_m[19] = 15.9734;
  mmht_xsec[20]= 1.31155;   mmht_pdf_p[20] = 0.01070;   mmht_pdf_m[20] = 0.00773;  mmht_scale_p[20] = 0.02770;  mmht_scale_m[20] = 0.02550;
  mmht_xsec[21]= 7.76695;   mmht_pdf_p[21] = 0.05739;   mmht_pdf_m[21] = 0.06698;  mmht_scale_p[21] = 0.17583;  mmht_scale_m[21] = 0.16471;
  mmht_xsec[22]= 5.92192;   mmht_pdf_p[22] = 0.03241;   mmht_pdf_m[22] = 0.04923;  mmht_scale_p[22] = 0.12937;  mmht_scale_m[22] = 0.12280;
  mmht_xsec[23]= 13.6888;   mmht_pdf_p[23] = 0.08998;   mmht_pdf_m[23] = 0.11624;  mmht_scale_p[23] = 0.27622;  mmht_scale_m[23] = 0.25656;  

  plot_mmht[0] = mmht_xsec[itype];  // value
  plot_mmht[1] = mmht_pdf_p[itype];   // pdf uncertainty up
  plot_mmht[2] = mmht_pdf_m[itype];   // pdf uncertainty down
  plot_mmht[3] = mmht_scale_p[itype];  // total uncertainty up 
  plot_mmht[4] = mmht_scale_m[itype];  // total uncertainty down

  Double_t abm_xsec[24];
  Double_t abm_pdf_p[24];
  Double_t abm_pdf_m[24];
  Double_t abm_scale_p[24];
  Double_t abm_scale_m[24];
  Double_t plot_abm[5];

  abm_xsec[0]= 11725.7;   abm_pdf_p[0] = 94.2999;   abm_pdf_m[0] = 112.752;   abm_scale_p[0] = 145.694;  abm_scale_m[0] = 130.146;   
  abm_xsec[1]= 8554.61;   abm_pdf_p[1] = 65.3948;   abm_pdf_m[1] = 75.098;    abm_scale_p[1] = 113.011;  abm_scale_m[1] = 90.0871;   
  abm_xsec[2]= 20280.3;   abm_pdf_p[2] = 158.667;   abm_pdf_m[2] = 185.508;   abm_scale_p[2] = 258.165;  abm_scale_m[2] = 217.993;   
  abm_xsec[3]= 1920.11;   abm_pdf_p[3] = 15.0242;   abm_pdf_m[3] = 16.7455;   abm_scale_p[3] = 19.8174;  abm_scale_m[3] = 19.5216; 
  abm_xsec[4]= 1.37069;   abm_pdf_p[4] = 0.002657;  abm_pdf_m[4] = 0.004092;  abm_scale_p[4] = 0.00321;  abm_scale_m[4] = 0.00411; 
  abm_xsec[5]= 6.10678;   abm_pdf_p[5] = 0.012549;  abm_pdf_m[5] = 0.013113;  abm_scale_p[5] = 0.02151;  abm_scale_m[5] = 0.01316; 
  abm_xsec[6]= 4.45527;   abm_pdf_p[6] = 0.007909;  abm_pdf_m[6] = 0.000587;  abm_scale_p[6] = 0.02013;  abm_scale_m[6] = 0.00621; 
  abm_xsec[7]= 10.5621;   abm_pdf_p[7] = 0.020458;  abm_pdf_m[7] = 0.018983;  abm_scale_p[7] = 0.04164;  abm_scale_m[7] = 0.01937;   

  abm_xsec[8 ]= 5190.05;   abm_pdf_p[8 ] = 41.7392;  abm_pdf_m[8 ] = 49.9065;  abm_scale_p[8 ] = 121.611;  abm_scale_m[8 ] = 118.106;  
  abm_xsec[9 ]= 4002.31;   abm_pdf_p[9 ] = 30.5952;  abm_pdf_m[9 ] = 35.1349;  abm_scale_p[9 ] = 94.5485;  abm_scale_m[9 ] = 78.3829;  
  abm_xsec[10]= 9192.36;   abm_pdf_p[10] = 71.9183;  abm_pdf_m[10] = 84.0844;  abm_scale_p[10] = 195.183;  abm_scale_m[10] = 156.216;  
  abm_xsec[11]= 703.505;   abm_pdf_p[11] = 5.50468;  abm_pdf_m[11] = 6.13535;  abm_scale_p[11] = 14.7878;  abm_scale_m[11] = 12.8825;  
  abm_xsec[12]= 1.29676;   abm_pdf_p[12] = 0.00251;  abm_pdf_m[12] = 0.00387;  abm_scale_p[12] = 0.02981;  abm_scale_m[12] = 0.02966;  
  abm_xsec[13]= 7.37742;   abm_pdf_p[13] = 0.01516;  abm_pdf_m[13] = 0.01584;  abm_scale_p[13] = 0.14972;  abm_scale_m[13] = 0.14745;  
  abm_xsec[14]= 5.68910;   abm_pdf_p[14] = 0.01009;  abm_pdf_m[14] = 0.00749;  abm_scale_p[14] = 0.11062;  abm_scale_m[14] = 0.10759;  
  abm_xsec[15]= 13.0665;   abm_pdf_p[15] = 0.02530;  abm_pdf_m[15] = 0.02348;  abm_scale_p[15] = 0.21933;  abm_scale_m[15] = 0.21319;  
				      		         		      			   
  abm_xsec[16]= 5029.57;   abm_pdf_p[16] = 40.4486;  abm_pdf_m[16] = 48.3633;  abm_scale_p[16] = 100.452;  abm_scale_m[16] = 96.4453;  
  abm_xsec[17]= 3847.03;   abm_pdf_p[17] = 29.4082;  abm_pdf_m[17] = 33.7717;  abm_scale_p[17] = 83.6617;  abm_scale_m[17] = 66.4563;  
  abm_xsec[18]= 8904.45;   abm_pdf_p[18] = 69.6658;  abm_pdf_m[18] = 81.4508;  abm_scale_p[18] = 192.494;  abm_scale_m[18] = 155.580;  
  abm_xsec[19]= 647.812;   abm_pdf_p[19] = 5.06890;  abm_pdf_m[19] = 5.64964;  abm_scale_p[19] = 13.9480;  abm_scale_m[19] = 12.2411;  
  abm_xsec[20]= 1.30739;   abm_pdf_p[20] = 0.00253;  abm_pdf_m[20] = 0.00390;  abm_scale_p[20] = 0.02528;  abm_scale_m[20] = 0.02507;  
  abm_xsec[21]= 7.76394;   abm_pdf_p[21] = 0.01595;  abm_pdf_m[21] = 0.01667;  abm_scale_p[21] = 0.15579;  abm_scale_m[21] = 0.15337;  
  abm_xsec[22]= 5.93849;   abm_pdf_p[22] = 0.01054;  abm_pdf_m[22] = 0.00782;  abm_scale_p[22] = 0.11527;  abm_scale_m[22] = 0.11210;  
  abm_xsec[23]= 13.7454;   abm_pdf_p[23] = 0.02662;  abm_pdf_m[23] = 0.02470;  abm_scale_p[23] = 0.23746;  abm_scale_m[23] = 0.23119;  

  plot_abm[0] = abm_xsec[itype];  // value
  plot_abm[1] = abm_pdf_p[itype];   // pdf uncertainty up
  plot_abm[2] = abm_pdf_m[itype];   // pdf uncertainty down
  plot_abm[3] = abm_scale_p[itype];  // total uncertainty up 
  plot_abm[4] = abm_scale_m[itype];  // total uncertainty down

  Double_t hera_xsec[24];
  Double_t hera_pdf_p[24];
  Double_t hera_pdf_m[24];
  Double_t hera_scale_p[24];
  Double_t hera_scale_m[24];
  Double_t plot_hera[5];

  hera_xsec[0]= 11775.9;   hera_pdf_p[0] = 561.572;  hera_pdf_m[0] = 237.466;  hera_scale_p[0] = 572.541;  hera_scale_m[0] = 246.275;  
  hera_xsec[1]= 8703.65;   hera_pdf_p[1] = 384.758;  hera_pdf_m[1] = 162.646;  hera_scale_p[1] = 396.021;  hera_scale_m[1] = 170.343;  
  hera_xsec[2]= 20479.4;   hera_pdf_p[2] = 942.095;  hera_pdf_m[2] = 389.003;  hera_scale_p[2] = 964.280;  hera_scale_m[2] = 405.820;  
  hera_xsec[3]= 1929.64;   hera_pdf_p[3] = 90.7990;  hera_pdf_m[3] = 40.6933;  hera_scale_p[3] = 91.7231;  hera_scale_m[3] = 41.9241;
  hera_xsec[4]= 1.35297;   hera_pdf_p[4] = 0.01345;  hera_pdf_m[4] = 0.01259;  hera_scale_p[4] = 0.01357;  hera_scale_m[4] = 0.01259;
  hera_xsec[5]= 6.10249;   hera_pdf_p[5] = 0.06203;  hera_pdf_m[5] = 0.05512;  hera_scale_p[5] = 0.06444;  hera_scale_m[5] = 0.05513;
  hera_xsec[6]= 4.51043;   hera_pdf_p[6] = 0.03613;  hera_pdf_m[6] = 0.03299;  hera_scale_p[6] = 0.04070;  hera_scale_m[6] = 0.03306;
  hera_xsec[7]= 10.61292;  hera_pdf_p[7] = 0.09816;  hera_pdf_m[7] = 0.08811;  hera_scale_p[7] = 0.10514;  hera_scale_m[7] = 0.08819;  

  hera_xsec[8 ]= 5160.55;   hera_pdf_p[8 ] = 246.097;  hera_pdf_m[8 ] = 104.064;  hera_scale_p[8 ] = 271.041;  hera_scale_m[8 ] = 148.856; 
  hera_xsec[9 ]= 4017.21;   hera_pdf_p[9 ] = 177.587;  hera_pdf_m[9 ] = 75.0700;  hera_scale_p[9 ] = 209.742;  hera_scale_m[9 ] = 102.866; 
  hera_xsec[10]= 9177.69;   hera_pdf_p[10] = 422.193;  hera_pdf_m[10] = 174.329;  hera_scale_p[10] = 484.158;  hera_scale_m[10] = 218.331; 
  hera_xsec[11]= 694.662;   hera_pdf_p[11] = 32.6872;  hera_pdf_m[11] = 14.6494;  hera_scale_p[11] = 37.8158;  hera_scale_m[11] = 18.4314; 
  hera_xsec[12]= 1.28459;   hera_pdf_p[12] = 0.01277;  hera_pdf_m[12] = 0.01195;  hera_scale_p[12] = 0.03402;  hera_scale_m[12] = 0.03149; 
  hera_xsec[13]= 7.42867;   hera_pdf_p[13] = 0.07550;  hera_pdf_m[13] = 0.06710;  hera_scale_p[13] = 0.18013;  hera_scale_m[13] = 0.16215; 
  hera_xsec[14]= 5.78287;   hera_pdf_p[14] = 0.04632;  hera_pdf_m[14] = 0.04229;  hera_scale_p[14] = 0.12812;  hera_scale_m[14] = 0.11701; 
  hera_xsec[15]= 13.2115;   hera_pdf_p[15] = 0.12219;  hera_pdf_m[15] = 0.10968;  hera_scale_p[15] = 0.27397;  hera_scale_m[15] = 0.24068; 
		   		        		    			  		        
  hera_xsec[16]= 5026.59;   hera_pdf_p[16] = 239.709;  hera_pdf_m[16] = 101.363;  hera_scale_p[16] = 256.648;  hera_scale_m[16] = 131.119; 
  hera_xsec[17]= 3865.39;   hera_pdf_p[17] = 170.875;  hera_pdf_m[17] = 72.2330;  hera_scale_p[17] = 198.917;  hera_scale_m[17] = 92.9248; 
  hera_xsec[18]= 8891.92;   hera_pdf_p[18] = 409.046;  hera_pdf_m[18] = 168.900;  hera_scale_p[18] = 470.780;  hera_scale_m[18] = 215.271; 
  hera_xsec[19]= 614.202;   hera_pdf_p[19] = 30.1717;  hera_pdf_m[19] = 13.5220;  hera_scale_p[19] = 35.0650;  hera_scale_m[19] = 17.3377; 
  hera_xsec[20]= 1.30039;   hera_pdf_p[20] = 0.01293;  hera_pdf_m[20] = 0.01210;  hera_scale_p[20] = 0.02995;  hera_scale_m[20] = 0.02696; 
  hera_xsec[21]= 7.83912;   hera_pdf_p[21] = 0.07967;  hera_pdf_m[21] = 0.07080;  hera_scale_p[21] = 0.18690;  hera_scale_m[21] = 0.16757; 
  hera_xsec[22]= 6.02825;   hera_pdf_p[22] = 0.04829;  hera_pdf_m[22] = 0.04409;  hera_scale_p[22] = 0.13423;  hera_scale_m[22] = 0.12271; 
  hera_xsec[23]= 13.8673;   hera_pdf_p[23] = 0.12825;  hera_pdf_m[23] = 0.11513;  hera_scale_p[23] = 0.29294;  hera_scale_m[23] = 0.25873;  

  plot_hera[0] = hera_xsec[itype];  // value
  plot_hera[1] = hera_pdf_p[itype];   // pdf uncertainty up
  plot_hera[2] = hera_pdf_m[itype];   // pdf uncertainty down
  plot_hera[3] = hera_scale_p[itype];  // total uncertainty up 
  plot_hera[4] = hera_scale_m[itype];  // total uncertainty down

  // set axis title, 
  TPaveText tb5(0.82,0.02,0.96,0.12,"NDC");
  char file_name[200]; 
  switch (itype)
    {
    case 0: 
      tb5.AddText("#sigma^{tot}_{W^{+}} [pb]"); lumi = true; sprintf(file_name,"pdf-wp-tot"); break;
    case 1: 
      tb5.AddText("#sigma^{tot}_{W^{-}} [pb]"); lumi = true; sprintf(file_name,"pdf-wm-tot");break;
    case 2: 
      tb5.AddText("#sigma^{tot}_{W} [pb]"); lumi = true; sprintf(file_name,"pdf-w-tot");break;
    case 3: 
      tb5.AddText("#sigma^{tot}_{Z} [pb]"); lumi = true; sprintf(file_name,"pdf-z-tot");break;
    case 4: 
      tb5.AddText("#sigma^{tot}_{W^{+}}/#sigma^{tot}_{W^{-}}"); lumi = false;  sprintf(file_name,"pdf-rpm-tot");break;  
    case 5: 
      tb5.AddText("#sigma^{tot}_{W^{+}}/#sigma^{tot}_{Z}"); lumi = false; sprintf(file_name,"pdf-wpr-tot");break;  
    case 6: 
      tb5.AddText("#sigma^{tot}_{W^{-}}/#sigma^{tot}_{Z}"); lumi = false; sprintf(file_name,"pdf-wmr-tot");break;  
    case 7: 
      tb5.AddText("#sigma^{tot}_{W}/#sigma^{tot}_{Z}"); lumi = false; sprintf(file_name,"pdf-wz-tot");break;  
    case 8: 
      tb5.AddText("#sigma^{fid}_{W^{+}} [pb]"); lumi = true; muon = true; sprintf(file_name,"pdf-wp-mfid"); break;
    case 9: 
      tb5.AddText("#sigma^{fid}_{W^{-}} [pb]"); lumi = true; muon = true; sprintf(file_name,"pdf-wm-mfid"); break;
    case 10: 
      tb5.AddText("#sigma^{fid}_{W} [pb]"); lumi = true; muon = true; sprintf(file_name,"pdf-w-mfid"); break;
    case 11: 
      tb5.AddText("#sigma^{fid}_{Z} [pb]"); lumi = true; muon = true; sprintf(file_name,"pdf-z-mfid"); break;
    case 12: 
      tb5.AddText("#sigma^{fid}_{W^{+}}/#sigma^{fid}_{W^{-}}"); lumi = false; muon = true; sprintf(file_name,"pdf-rpm-mfid"); break;  
    case 13: 
      tb5.AddText("#sigma^{fid}_{W^{+}}/#sigma^{fid}_{Z}"); lumi = false; muon = true; sprintf(file_name,"pdf-wpr-mfid"); break;  
    case 14: 
      tb5.AddText("#sigma^{fid}_{W^{-}}/#sigma^{fid}_{Z}"); lumi = false; muon = true; sprintf(file_name,"pdf-wmr-mfid"); break;  
    case 15: 
      tb5.AddText("#sigma^{fid}_{W}/#sigma^{fid}_{Z}"); lumi = false; muon = true; sprintf(file_name,"pdf-wz-mfid"); break;  
    case 16: 
      tb5.AddText("#sigma^{fid}_{W^{+}} [pb]"); lumi = true; electron = true; sprintf(file_name,"pdf-wp-efid"); break;
    case 17: 
      tb5.AddText("#sigma^{fid}_{W^{-}} [pb]"); lumi = true; electron = true; sprintf(file_name,"pdf-wm-efid"); break;
    case 18: 
      tb5.AddText("#sigma^{fid}_{W} [pb]"); lumi = true; electron = true; sprintf(file_name,"pdf-w-efid"); break;
    case 19: 
      tb5.AddText("#sigma^{fid}_{Z} [pb]"); lumi = true; electron = true; sprintf(file_name,"pdf-z-efid"); break;
    case 20: 
      tb5.AddText("#sigma^{fid}_{W^{+}}/#sigma^{fid}_{W^{-}}"); lumi = false; electron = true; sprintf(file_name,"pdf-rpm-efid"); break;  
    case 21: 
      tb5.AddText("#sigma^{fid}_{W^{+}}/#sigma^{fid}_{Z}"); lumi = false; electron = true; sprintf(file_name,"pdf-wpr-efid"); break;  
    case 22: 
      tb5.AddText("#sigma^{fid}_{W^{-}}/#sigma^{fid}_{Z}"); lumi = false; electron = true; sprintf(file_name,"pdf-wmr-efid"); break;  
    case 23: 
      tb5.AddText("#sigma^{fid}_{W}/#sigma^{fid}_{Z}"); lumi = false; electron = true; sprintf(file_name,"pdf-wz-efid"); break;  
    default:
      tb5.AddText("What is this?"); sprintf(file_name,"pdf-dummy"); break;
    }

  //==============================================================================================================  
  
  //--------------------------------------------------------------------------------------------------------------
  // plotting parameter
  Double_t yshift, ydrift; // used for point and error bar
  yshift = 3.50;
  ydrift = -0.72;
  

  Double_t td, tp, tw; // used for text 
  td = -0.037;
  tp = 0.70;
  tw = 0.04;
      
  

  Double_t range_max, range_min;
  //range_min = 0.82;
  //range_max = 1.27;
  if (lumi) 
    {
      range_min = plot_xsec[0]*0.85;
      range_max = plot_xsec[0]*1.12;
    }
  else
    {
      range_min = plot_xsec[0]*0.93;
      range_max = plot_xsec[0]*1.06;
    }
  Int_t expColor, theColor;  // colors used
  Int_t prodColor, decayColor, cmbColor;  // colors used
  expColor  = 33;
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

  Double_t xval, errlp, errhp, errlm, errhm, yval;
  xval = plot_nnpdf[0];
  yval = yshift + (0*ydrift);
  errlp = plot_nnpdf[1] ;
  errlm = plot_nnpdf[2] ;
  errhp = plot_nnpdf[3] ;
  errhm = plot_nnpdf[4] ;
  TGraphAsymmErrors grWP(1,&xval,&yval,&errlm,&errlp,0,0);
 
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
  grWP.SetLineColor(kBlack);
  grWP.Draw("AP");

  // lumi uncertainty band
  TBox lumi_box(0.973*plot_xsec[0],0.01,1.027*plot_xsec[0],4.);
  lumi_box.SetLineColor(796);
  lumi_box.SetFillColor(796);
  if (lumi) lumi_box.Draw();
  c->RedrawAxis("SAME");

  Double_t unc = sqrt(plot_xsec[1]*plot_xsec[1]+plot_xsec[2]*plot_xsec[2]);

  TBox exp_box(plot_xsec[0]-unc,0.01,plot_xsec[0]+unc,4.);
  exp_box.SetLineColor(expColor);
  exp_box.SetFillColor(expColor);
  exp_box.Draw();
  c->RedrawAxis("SAME");

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
  tb2.AddText("2.3 fb^{-1} (13 TeV)");
  tb2.Draw(); 

  TPaveText tb3(0.70,0.83,0.95,0.88,"NDC");
  tb3.SetFillStyle(0);
  tb3.SetBorderSize(0);
  tb3.SetTextAlign(12);
  tb3.AddText("Theory: FEWZ (NNLO)");
  tb3.Draw(); 

  TPaveText tb7(0.70,0.785,0.95,0.835,"NDC");
  tb7.SetFillStyle(0);
  tb7.SetBorderSize(0);
  tb7.SetTextAlign(12);
  if (electron)   
    {
      TPaveText tb7(0.705,0.785,0.89,0.835,"NDC");
      tb7.SetFillStyle(0);
      tb7.SetBorderSize(0);
      tb7.SetTextAlign(12);
      tb7.AddText("Electron channel");
      tb7.Draw(); 
    }
  else if (muon) 
    {
      TPaveText tb7(0.705,0.785,0.87,0.835,"NDC");
      tb7.SetFillStyle(0);
      tb7.SetBorderSize(0);
      tb7.SetTextAlign(12);
      tb7.AddText("Muon channel");
      tb7.Draw(); 
    }
  else 
    {
      TPaveText tb7(0.70,0.785,0.95,0.835,"NDC");
      tb7.SetFillStyle(0);
      tb7.SetBorderSize(0);
      tb7.SetTextAlign(12);
      tb7.AddText("Observation: NNPDF3.0");
      tb7.Draw(); 
    }

  TLine theory_line(plot_xsec[0],0.01,plot_xsec[0],4.);
  theory_line.SetLineColor(kRed);
  theory_line.SetLineStyle(1);
  theory_line.SetLineWidth(3);
  theory_line.Draw();

  // legend
  TLine theory_line_leg(range_min+(0.04*(range_max-range_min)),4.49,range_min+(0.04*(range_max-range_min)),4.69);
  theory_line_leg.SetLineColor(kRed);
  theory_line_leg.SetLineStyle(1);
  theory_line_leg.SetLineWidth(3);
  theory_line_leg.Draw();

  if (lumi)
    {
      TBox exp_box_leg(range_min+(0.02*(range_max-range_min)),4.195,range_min+(0.06*(range_max-range_min)),4.395);
      exp_box_leg.SetLineColor(796);
      exp_box_leg.SetFillColor(796);
      exp_box_leg.Draw();

      TBox exp_box_leg2(range_min+(0.03*(range_max-range_min)),4.195,range_min+(0.05*(range_max-range_min)),4.395);
      exp_box_leg2.SetLineColor(expColor);
      exp_box_leg2.SetFillColor(expColor);
      exp_box_leg2.Draw();
    }
  else
    {
      TBox exp_box_leg2(range_min+(0.02*(range_max-range_min)),4.195,range_min+(0.06*(range_max-range_min)),4.395);
      exp_box_leg2.SetLineColor(expColor);
      exp_box_leg2.SetFillColor(expColor);
      exp_box_leg2.Draw();
    }

  TPaveText tb4(0.12,0.84,0.3,0.88,"NDC");
  tb4.SetFillStyle(0);
  tb4.SetBorderSize(0);
  tb4.SetTextAlign(12);
  tb4.AddText("Observation");
  tb4.Draw(); 

  TPaveText tb6(0.11,0.79,0.5,0.83,"NDC");
  tb6.SetFillStyle(0);
  tb6.SetBorderSize(0);
  tb6.SetTextAlign(12);
  if (lumi)
    tb6.AddText("Uncertainty (exp., exp. #oplus lumi)");
  else
    tb6.AddText("Uncertainty");
  tb6.Draw();

  TPaveText tb8(0.08,0.17,0.5,0.20,"NDC");
  tb8.SetFillStyle(0);
  tb8.SetBorderSize(0);
  tb8.SetTextAlign(12);
  tb8.AddText("(inner uncertainty: PDF only)");
  tb8.Draw();

  tb5.SetFillStyle(0);
  tb5.SetBorderSize(0);
  tb5.SetTextAlign(12);
  tb5.Draw();

  TGraphAsymmErrors grWP2(1,&xval,&yval,&errhm,&errhp,0,0);
  grWP2.SetMarkerStyle(kFullCircle);
  grWP2.SetMarkerSize(1);
  grWP2.SetLineWidth(2);
  grWP2.SetMarkerColor(kBlack);
  grWP2.SetLineColor(kBlack);
  grWP.Draw("EPSAME");
  grWP2.Draw("EPSAME");

  // add CT14 result
  xval = plot_ct14[0];
  yval = yshift + (1*ydrift);
  errlp = plot_ct14[1];
  errlm = plot_ct14[2];
  errhp = plot_ct14[3];
  errhm = plot_ct14[4];

  TGraphAsymmErrors grCT(1,&xval,&yval,&errlm,&errlp,0,0);
  grCT.SetMarkerStyle(kFullCircle);
  grCT.SetMarkerSize(1);
  grCT.SetLineWidth(2);
  grCT.SetMarkerColor(kBlack);
  grCT.SetLineColor(kBlack);
  TGraphAsymmErrors grCT2(1,&xval,&yval,&errhm,&errhp,0,0);
  grCT2.SetMarkerStyle(kFullCircle);
  grCT2.SetMarkerSize(1);
  grCT2.SetLineWidth(2);
  grCT2.SetMarkerColor(kBlack);
  grCT2.SetLineColor(kBlack);
  grCT2.Draw("EPSAME");
  grCT.Draw("EPSAME");
  // add MMHT result
  xval = plot_mmht[0];
  yval = yshift + (2*ydrift);
  errlp = plot_mmht[1];
  errlm = plot_mmht[2];
  errhp = plot_mmht[3];
  errhm = plot_mmht[4];
  TGraphAsymmErrors grMM(1,&xval,&yval,&errlm,&errlp,0,0);
  grMM.SetMarkerStyle(kFullCircle);
  grMM.SetMarkerSize(1);
  grMM.SetLineWidth(2);
  grMM.SetMarkerColor(kBlack);
  grMM.SetLineColor(kBlack);
  TGraphAsymmErrors grMM2(1,&xval,&yval,&errhm,&errhp,0,0);
  grMM2.SetMarkerStyle(kFullCircle);
  grMM2.SetMarkerSize(1);
  grMM2.SetLineWidth(2);
  grMM2.SetMarkerColor(kBlack);
  grMM2.SetLineColor(kBlack);
  grMM2.Draw("EPSAME");
  grMM.Draw("EPSAME");

  xval = plot_abm[0];
  yval = yshift + (3*ydrift);
  errlp = plot_abm[1];
  errlm = plot_abm[2];
  errhp = plot_abm[3];
  errhm = plot_abm[4];
  TGraphAsymmErrors grAB(1,&xval,&yval,&errlm,&errlp,0,0);
  grAB.SetMarkerStyle(kFullCircle);
  grAB.SetMarkerSize(1);
  grAB.SetLineWidth(2);
  grAB.SetMarkerColor(kBlack);
  grAB.SetLineColor(kBlack);
  TGraphAsymmErrors grAB2(1,&xval,&yval,&errhm,&errhp,0,0);
  grAB2.SetMarkerStyle(kFullCircle);
  grAB2.SetMarkerSize(1);
  grAB2.SetLineWidth(2);
  grAB2.SetMarkerColor(kBlack);
  grAB2.SetLineColor(kBlack);
  grAB2.Draw("EPSAME");
  grAB.Draw("EPSAME");

  xval = plot_hera[0];
  yval = yshift + (4*ydrift);
  errlp = plot_hera[1];
  errlm = plot_hera[2];
  errhp = plot_hera[3];
  errhm = plot_hera[4];
  TGraphAsymmErrors grHE(1,&xval,&yval,&errlm,&errlp,0,0);
  grHE.SetMarkerStyle(kFullCircle);
  grHE.SetMarkerSize(1);
  grHE.SetLineWidth(2);
  grHE.SetMarkerColor(kBlack);
  grHE.SetLineColor(kBlack);
  TGraphAsymmErrors grHE2(1,&xval,&yval,&errhm,&errhp,0,0);
  grHE2.SetMarkerStyle(kFullCircle);
  grHE2.SetMarkerSize(1);
  grHE2.SetLineWidth(2);
  grHE2.SetMarkerColor(kBlack);
  grHE2.SetLineColor(kBlack);
  grHE2.Draw("EPSAME");
  grHE.Draw("EPSAME");

  // add text
  char buffer[200]; 

  //sprintf(buffer,"%.2f #pm %.2f_{stat} #pm %.2f_{syst}", xsec_wp, xsec_wp_stat, xsec_wp_sys);
  //sprintf(buffer,"%s #pm %.2f_{lum} nb",buffer,xsec_wp_lumi);
  //TPaveText resultwp(0.08,0.65,0.40,0.70,"NDC");
  //resultwp.SetFillStyle(0);
  //resultwp.SetBorderSize(0);
  //resultwp.SetTextAlign(12);
  //resultwp.AddText(buffer);
  //resultwp.Draw(); 

  sprintf(buffer,"#bf{NNPDF3.0}");
  TPaveText tNNPDF(0.08,tp+(0*td),0.40,tp+tw+(0*td),"NDC");
  tNNPDF.SetFillStyle(0);
  tNNPDF.SetBorderSize(0);
  tNNPDF.SetTextAlign(12);
  tNNPDF.AddText(buffer);
  tNNPDF.Draw(); 
  if (itype == 4)
    sprintf(buffer,"%.3f^{+%.3f}_{-%.3f}", plot_nnpdf[0], plot_nnpdf[3], plot_nnpdf[4]);
  else 
    sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_nnpdf[0], plot_nnpdf[3], plot_nnpdf[4]);
  if (lumi){ 
    roundXsec(plot_nnpdf[0], plot_nnpdf[3], plot_nnpdf[4]);
    sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_nnpdf[0], plot_nnpdf[3], plot_nnpdf[4]);
    sprintf(buffer,"%s pb",buffer);
  }

  TPaveText theoryNNPDF(0.10,tp+(1*td),0.42,tp+tw+(1*td),"NDC");
  theoryNNPDF.SetFillStyle(0);
  theoryNNPDF.SetBorderSize(0);
  theoryNNPDF.SetTextAlign(12);
  theoryNNPDF.AddText(buffer);
  theoryNNPDF.Draw(); 

  sprintf(buffer,"#bf{CT14}");
  TPaveText tCT14(0.08,tp+(3*td),0.40,tp+tw+(3*td),"NDC");
  tCT14.SetFillStyle(0);
  tCT14.SetBorderSize(0);
  tCT14.SetTextAlign(12);
  tCT14.AddText(buffer);
  tCT14.Draw(); 
  if (itype == 4)
    sprintf(buffer,"%.3f^{+%.3f}_{-%.3f}", plot_ct14[0], plot_ct14[3], plot_ct14[4]);
  else
    sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_ct14[0], plot_ct14[3], plot_ct14[4]);
  if (lumi) {
    roundXsec(plot_ct14[0], plot_ct14[3], plot_ct14[4]);
    sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_ct14[0], plot_ct14[3], plot_ct14[4]);
    sprintf(buffer,"%s pb",buffer);
  }
  TPaveText theoryCT14(0.10,tp+(4*td),0.42,tp+tw+(4*td),"NDC");
  theoryCT14.SetFillStyle(0);
  theoryCT14.SetBorderSize(0);
  theoryCT14.SetTextAlign(12);
  theoryCT14.AddText(buffer);
  theoryCT14.Draw(); 

  sprintf(buffer,"#bf{MMHT2014}");
  TPaveText tMMHT(0.08,tp+(6*td),0.40,tp+tw+(6*td),"NDC");
  tMMHT.SetFillStyle(0);
  tMMHT.SetBorderSize(0);
  tMMHT.SetTextAlign(12);
  tMMHT.AddText(buffer);
  tMMHT.Draw(); 
  if (itype == 4)
    sprintf(buffer,"%.3f^{+%.3f}_{-%.3f}", plot_mmht[0], plot_mmht[3], plot_mmht[4]);
  else
    sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_mmht[0], plot_mmht[3], plot_mmht[4]);
  
  if (lumi) {
    roundXsec(plot_mmht[0], plot_mmht[3], plot_mmht[4]);
    sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_mmht[0], plot_mmht[3], plot_mmht[4]);
    sprintf(buffer,"%s pb",buffer);
  }
  TPaveText theoryMMHT(0.10,tp+(7*td),0.42,tp+tw+(7*td),"NDC");
  theoryMMHT.SetFillStyle(0);
  theoryMMHT.SetBorderSize(0);
  theoryMMHT.SetTextAlign(12);
  theoryMMHT.AddText(buffer);
  theoryMMHT.Draw(); 

  sprintf(buffer,"#bf{ABM12LHC}");
  TPaveText tABM(0.08,tp+(9*td),0.40,tp+tw+(9*td),"NDC");
  tABM.SetFillStyle(0);
  tABM.SetBorderSize(0);
  tABM.SetTextAlign(12);
  tABM.AddText(buffer);
  tABM.Draw(); 

  if (itype == 4)
    sprintf(buffer,"%.3f^{+%.3f}_{-%.3f}", plot_abm[0], plot_abm[3], plot_abm[4]);
  else
    sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_abm[0], plot_abm[3], plot_abm[4]);
  if (lumi) {
    roundXsec(plot_abm[0], plot_abm[3], plot_abm[4]);
    sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_abm[0], plot_abm[3], plot_abm[4]);
    sprintf(buffer,"%s pb",buffer);
  }
  TPaveText theoryABM(0.10,tp+(10*td),0.42,tp+tw+(10*td),"NDC");
  theoryABM.SetFillStyle(0);
  theoryABM.SetBorderSize(0);
  theoryABM.SetTextAlign(12);
  theoryABM.AddText(buffer);
  theoryABM.Draw(); 

  sprintf(buffer,"#bf{HERAPDF15}");
  TPaveText tHERA(0.08,tp+(12*td),0.40,tp+tw+(12*td),"NDC");
  tHERA.SetFillStyle(0);
  tHERA.SetBorderSize(0);
  tHERA.SetTextAlign(12);
  tHERA.AddText(buffer);
  tHERA.Draw(); 
  if (itype == 4)
    sprintf(buffer,"%.3f^{+%.3f}_{-%.3f}", plot_hera[0], plot_hera[3], plot_hera[4]);
  else
    sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_hera[0], plot_hera[3], plot_hera[4]);
  if (lumi) {
    roundXsec(plot_hera[0], plot_hera[3], plot_hera[4]);
    sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_hera[0], plot_hera[3], plot_hera[4]);
    sprintf(buffer,"%s pb",buffer);
  }

  TPaveText theoryHERA(0.10,tp+(13*td),0.42,tp+tw+(13*td),"NDC");
  theoryHERA.SetFillStyle(0);
  theoryHERA.SetBorderSize(0);
  theoryHERA.SetTextAlign(12);
  theoryHERA.AddText(buffer);
  theoryHERA.Draw(); 

  sprintf(buffer,"%s.png",file_name);
  c->SaveAs(buffer);
  sprintf(buffer,"%s.pdf",file_name);
  c->SaveAs(buffer);

  // make latex table for note
  if (itype <= 3)
    {
  cout << file_name << "& $";
  sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_nnpdf[0], plot_nnpdf[3], plot_nnpdf[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_ct14[0], plot_ct14[3], plot_ct14[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_mmht[0], plot_mmht[3], plot_mmht[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_abm[0], plot_abm[3], plot_abm[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.0f^{+%.0f}_{-%.0f}", plot_hera[0], plot_hera[3], plot_hera[4]);
  cout << buffer << "$" << endl;
    }
  else
    {
  cout << file_name << "& $";
  sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_nnpdf[0], plot_nnpdf[3], plot_nnpdf[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_ct14[0], plot_ct14[3], plot_ct14[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_mmht[0], plot_mmht[3], plot_mmht[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_abm[0], plot_abm[3], plot_abm[4]);
  cout << buffer << "$ & $";
  sprintf(buffer,"%.2f^{+%.2f}_{-%.2f}", plot_hera[0], plot_hera[3], plot_hera[4]);
  cout << buffer << "$" << endl;
    }
    
}
