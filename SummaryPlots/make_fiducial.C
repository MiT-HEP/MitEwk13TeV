//================================================================================================
//
// Compute fiducial cross-sections
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#endif

void make_fiducial(const TString infilename="make_fiducial.txt")
{   

  Double_t wpe, wpm, wpl, wpe_stat, wpm_stat, wpl_stat, wpe_sys, wpm_sys, wpl_sys, wpe_lumi, wpm_lumi, wpl_lumi;
  Double_t wme, wmm, wml, wme_stat, wmm_stat, wml_stat, wme_sys, wmm_sys, wml_sys, wme_lumi, wmm_lumi, wml_lumi;
  Double_t we, wm, wl, we_stat, wm_stat, wl_stat, we_sys, wm_sys, wl_sys, we_lumi, wm_lumi, wl_lumi;
  Double_t ze, zm, zl, ze_stat, zm_stat, zl_stat, ze_sys, zm_sys, zl_sys, ze_lumi, zm_lumi, zl_lumi;
  Double_t wre, wrm, wrl, wre_stat, wrm_stat, wrl_stat, wre_sys, wrm_sys, wrl_sys, wre_lumi, wrm_lumi, wrl_lumi;
  Double_t wzre, wzrm, wzrl, wzre_stat, wzrm_stat, wzrl_stat, wzre_sys, wzrm_sys, wzrl_sys, wzre_lumi, wzrm_lumi, wzrl_lumi;
  Double_t theo_wp, theo_wm, theo_w, theo_z, theo_wr, theo_wzr, theo_wp_unc, theo_wm_unc, theo_w_unc, theo_z_unc, theo_wr_unc, theo_wzr_unc;
  Double_t acc_ze, acc_we, acc_wpe, acc_wme, acc_zm, acc_wm, acc_wpm, acc_wmm, acc_wre, acc_wzre, acc_wrm, acc_wzrm;
  Double_t acc_ze_unc, acc_we_unc, acc_wpe_unc, acc_wme_unc, acc_zm_unc, acc_wm_unc, acc_wpm_unc, acc_wmm_unc;
  Double_t acc_wre_unc, acc_wzre_unc, acc_wrm_unc, acc_wzrm_unc;
  
  //--------------------------------------------------------------------------------------------------------------
  // Parse text file
  //==============================================================================================================  
  
  ifstream ifs;
  ifs.open(infilename.Data());
  assert(ifs.is_open());
  string line;
  string label;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='*') continue;

    // w+(e)
    if (state==0) {
      stringstream ss(line);
      ss >> label >> wpe >> label >> wpe_stat >> label >> label >> wpe_sys >> label >> label >> wpe_lumi;
      state++; continue;
    }
    // w+(m)
    if (state==1) {
      stringstream ss(line);
      ss >> label >> wpm >> label >> wpm_stat >> label >> label >> wpm_sys >> label >> label >> wpm_lumi;
      state++; continue;
    }
    // w+(l)
    if (state==2) {
      stringstream ss(line);
      ss >> label >> wpl >> label >> wpl_stat >> label >> label >> wpl_sys >> label >> label >> wpl_lumi;
      state++; continue;
    }
    // w-(e)
    if (state==3) {
      stringstream ss(line);
      ss >> label >> wme >> label >> wme_stat >> label >> label >> wme_sys >> label >> label >> wme_lumi;
      state++; continue;
    }
    // w-(m)
    if (state==4) {
      stringstream ss(line);
      ss >> label >> wmm >> label >> wmm_stat >> label >> label >> wmm_sys >> label >> label >> wmm_lumi;
      state++; continue;
    }
    // w-(l)
    if (state==5) {
      stringstream ss(line);
      ss >> label >> wml >> label >> wml_stat >> label >> label >> wml_sys >> label >> label >> wml_lumi;
      state++; continue;
    }
    // w(e)
    if (state==6) {
      stringstream ss(line);
      ss >> label >> we >> label >> we_stat >> label >> label >> we_sys >> label >> label >> we_lumi;
      state++; continue;
    }
    // w(m)
    if (state==7) {
      stringstream ss(line);
      ss >> label >> wm >> label >> wm_stat >> label >> label >> wm_sys >> label >> label >> wm_lumi;
      state++; continue;
    }
    // w(l)
    if (state==8) {
      stringstream ss(line);
      ss >> label >> wl >> label >> wl_stat >> label >> label >> wl_sys >> label >> label >> wl_lumi;
      state++; continue;
    }
    // z(e)
    if (state==9) {
      stringstream ss(line);
      ss >> label >> ze >> label >> ze_stat >> label >> label >> ze_sys >> label >> label >> ze_lumi;
      state++; continue;
    }
    // z(m)
    if (state==10) {
      stringstream ss(line);
      ss >> label >> zm >> label >> zm_stat >> label >> label >> zm_sys >> label >> label >> zm_lumi;
      state++; continue;
    }
    // z(l)
    if (state==11) {
      stringstream ss(line);
      ss >> label >> zl >> label >> zl_stat >> label >> label >> zl_sys >> label >> label >> zl_lumi;
      state++; continue;
    }
    // w+/w-(e)
    if (state==12) {
      stringstream ss(line);
      ss >> label >> wre >> label >> wre_stat >> label >> label >> wre_sys >> label >> label >> wre_lumi;
      state++; continue;
    }
    // w+/w-(m)
    if (state==13) {
      stringstream ss(line);
      ss >> label >> wrm >> label >> wrm_stat >> label >> label >> wrm_sys >> label >> label >> wrm_lumi;
      state++; continue;
    }
    // w+/w-(l)
    if (state==14) {
      stringstream ss(line);
      ss >> label >> wrl >> label >> wrl_stat >> label >> label >> wrl_sys >> label >> label >> wrl_lumi;
      state++; continue;
    }
    // w/z(e)
    if (state==15) {
      stringstream ss(line);
      ss >> label >> wzre >> label >> wzre_stat >> label >> label >> wzre_sys >> label >> label >> wzre_lumi;
      state++; continue;
    }
    // w/z(m)
    if (state==16) {
      stringstream ss(line);
      ss >> label >> wzrm >> label >> wzrm_stat >> label >> label >> wzrm_sys >> label >> label >> wzrm_lumi;
      state++; continue;
    }
    // w/z(l)
    if (state==17) {
      stringstream ss(line);
      ss >> label >> wzrl >> label >> wzrl_stat >> label >> label >> wzrl_sys >> label >> label >> wzrl_lumi;
      state++; continue;
    }
    // theo w+
    if (state==18) {
      stringstream ss(line);
      ss >> label >> theo_wp >> label >> theo_wp_unc;
      state++; continue;
    }
    // theo w-
    if (state==19) {
      stringstream ss(line);
      ss >> label >> theo_wm >> label >> theo_wm_unc;
      state++; continue;
    }
    // theo w
    if (state==20) {
      stringstream ss(line);
      ss >> label >> theo_w >> label >> theo_w_unc;
      state++; continue;
    }
    // theo z
    if (state==21) {
      stringstream ss(line);
      ss >> label >> theo_z >> label >> theo_z_unc;
      state++; continue;
    }
    // theo w+/w-
    if (state==22) {
      stringstream ss(line);
      ss >> label >> theo_wr >> label >> theo_wr_unc;
      state++; continue;
    }
    // theo w/z
    if (state==23) {
      stringstream ss(line);
      ss >> label >> theo_wzr >> label >> theo_wzr_unc;
      state++; continue;
    }
    // acc ze
    if (state==24) {
      stringstream ss(line);
      ss >> label >> acc_ze >> label >> acc_ze_unc;
      state++; continue;
    }
    // acc we
    if (state==25) {
      stringstream ss(line);
      ss >> label >> acc_we >> label >> acc_we_unc;
      state++; continue;
    }
    // acc wpe
    if (state==26) {
      stringstream ss(line);
      ss >> label >> acc_wpe >> label >> acc_wpe_unc;
      state++; continue;
    }
    // acc wme
    if (state==27) {
      stringstream ss(line);
      ss >> label >> acc_wme >> label >> acc_wme_unc;
      state++; continue;
    }
    // acc zm
    if (state==28) {
      stringstream ss(line);
      ss >> label >> acc_zm >> label >> acc_zm_unc;
      state++; continue;
    }
    // acc wm
    if (state==29) {
      stringstream ss(line);
      ss >> label >> acc_wm >> label >> acc_wm_unc;
      state++; continue;
    }
    // acc wpm
    if (state==30) {
      stringstream ss(line);
      ss >> label >> acc_wpm >> label >> acc_wpm_unc;
      state++; continue;
    }
    // acc wmm
    if (state==31) {
      stringstream ss(line);
      ss >> label >> acc_wmm >> label >> acc_wmm_unc;
      state++; continue;
    }
    // sys wpe
    if (state==32) {
      stringstream ss(line);
      ss >> label >> wpe_sys;
      state++; continue;
    }
    // sys wme
    if (state==33) {
      stringstream ss(line);
      ss >> label >> wme_sys;
      state++; continue;
    }
    // sys we
    if (state==34) {
      stringstream ss(line);
      ss >> label >> we_sys;
      state++; continue;
    }
    // sys wre
    if (state==35) {
      stringstream ss(line);
      ss >> label >> wre_sys;
      state++; continue;
    }
    // sys ze
    if (state==36) {
      stringstream ss(line);
      ss >> label >> ze_sys;
      state++; continue;
    }
    // sys wzre
    if (state==37) {
      stringstream ss(line);
      ss >> label >> wzre_sys;
      state++; continue;
    }
    // sys wpm
    if (state==38) {
      stringstream ss(line);
      ss >> label >> wpm_sys;
      state++; continue;
    }
    // sys wmm
    if (state==39) {
      stringstream ss(line);
      ss >> label >> wmm_sys;
      state++; continue;
    }
    // sys wm
    if (state==40) {
      stringstream ss(line);
      ss >> label >> wm_sys;
      state++; continue;
    }
    // sys wm
    if (state==41) {
      stringstream ss(line);
      ss >> label >> wrm_sys;
      state++; continue;
    }
    // sys wm
    if (state==42) {
      stringstream ss(line);
      ss >> label >> zm_sys;
      state++; continue;
    }
    // sys zm
    if (state==43) {
      stringstream ss(line);
      ss >> label >> wzrm_sys;
      state++; continue;
    }

    // theo wpe
    if (state==44) {
      stringstream ss(line);
      ss >> label >> acc_wpe_unc;
      state++; continue;
    }
    // theo wme
    if (state==45) {
      stringstream ss(line);
      ss >> label >> acc_wme_unc;
      state++; continue;
    }
    // theo we
    if (state==46) {
      stringstream ss(line);
      ss >> label >> acc_we_unc;
      state++; continue;
    }
    // theo wre
    if (state==47) {
      stringstream ss(line);
      ss >> label >> acc_wre_unc;
      state++; continue;
    }
    // theo ze
    if (state==48) {
      stringstream ss(line);
      ss >> label >> acc_ze_unc;
      state++; continue;
    }
    // theo wzre
    if (state==49) {
      stringstream ss(line);
      ss >> label >> acc_wzre_unc;
      state++; continue;
    }
    // theo wpm
    if (state==50) {
      stringstream ss(line);
      ss >> label >> acc_wpm_unc;
      state++; continue;
    }
    // theo wmm
    if (state==51) {
      stringstream ss(line);
      ss >> label >> acc_wmm_unc;
      state++; continue;
    }
    // theo wm
    if (state==52) {
      stringstream ss(line);
      ss >> label >> acc_wm_unc;
      state++; continue;
    }
    // theo wm
    if (state==53) {
      stringstream ss(line);
      ss >> label >> acc_wrm_unc;
      state++; continue;
    }
    // theo wm
    if (state==54) {
      stringstream ss(line);
      ss >> label >> acc_zm_unc;
      state++; continue;
    }
    // theo zm
    if (state==55) {
      stringstream ss(line);
      ss >> label >> acc_wzrm_unc;
      state++; continue;
    }

  }
  ifs.close();
  //cout << setprecision(3) << fixed; 

  //
  // FIDUCIAL CROSS SECTIONS
  //

  wpe*=acc_wpe; wpe_stat*=acc_wpe; wpe_lumi*=acc_wpe; Double_t wpe_new = wpe_sys*wpe/100;
  wme*=acc_wme; wme_stat*=acc_wme; wme_lumi*=acc_wme; Double_t wme_new = wme_sys*wme/100;

  we*=acc_we; we_stat*=acc_we; we_lumi*=acc_we; Double_t we_new = we_sys*we/100;
  wpm*=acc_wpm; wpm_stat*=acc_wpm; wpm_lumi*=acc_wpm; Double_t wpm_new = wpm_sys*wpm/100;

  wmm*=acc_wmm; wmm_stat*=acc_wmm; wmm_lumi*=acc_wmm; Double_t wmm_new = wmm_sys*wmm/100;
  wm*=acc_wm; wm_stat*=acc_wm; wm_lumi*=acc_wm; Double_t wm_new = wm_sys*wm/100;

  ze*=acc_ze; ze_stat*=acc_ze; ze_lumi*=acc_ze; Double_t ze_new = ze_sys*ze/100;
  zm*=acc_zm; zm_stat*=acc_zm; zm_lumi*=acc_zm; Double_t zm_new = zm_sys*zm/100;

  wre*=acc_wpe/acc_wme; wre_stat*=acc_wpe/acc_wme; wre_lumi*=acc_wpe/acc_wme;
  wre_sys*=wre/100;
  wre_lumi=0;

  wrm*=acc_wpm/acc_wmm; wrm_stat*=acc_wpm/acc_wmm; wrm_lumi*=acc_wpm/acc_wmm;
  wrm_sys*=wrm/100;
  wrm_lumi=0;

  wzre*=acc_we/acc_ze; wzre_stat*=acc_we/acc_ze; wzre_lumi*=acc_we/acc_ze;
  wzre_sys*=wzre/100;
  wzre_lumi=0;

  wzrm*=acc_wm/acc_zm; wzrm_stat*=acc_wm/acc_zm; wzrm_lumi*=acc_wm/acc_zm;
  wzrm_sys*=wzrm/100;
  wzrm_lumi=0;

  // 
  // COMBINE CHANNELS
  //

  Double_t wpe_weight = (1/(wpe_new*wpe_new+wpe_stat*wpe_stat))/(1/(wpe_new*wpe_new+wpe_stat*wpe_stat)+1/(wpm_new*wpm_new+wpm_stat*wpm_stat));
  Double_t wpm_weight = (1/(wpm_new*wpm_new+wpm_stat*wpm_stat))/(1/(wpe_new*wpe_new+wpe_stat*wpe_stat)+1/(wpm_new*wpm_new+wpm_stat*wpm_stat));

  cout << "wpe " << wpe_weight << " " << wpm_weight << endl;

  wpl_stat=TMath::Sqrt(wpe_weight*wpe_weight*wpe_stat*wpe_stat+wpm_weight*wpm_weight*wpm_stat*wpm_stat);
  wpl_sys=TMath::Sqrt(wpe_weight*wpe_weight*wpe_new*wpe_new+wpm_weight*wpm_weight*wpm_new*wpm_new);
  wpl_lumi=wpe_weight*wpe_lumi+wpm_weight*wpm_lumi;
  wpl = wpe*wpe_weight+wpm*wpm_weight;

  Double_t wme_weight = (1/(wme_new*wme_new+wme_stat*wme_stat))/(1/(wme_new*wme_new+wme_stat*wme_stat)+1/(wmm_new*wmm_new+wmm_stat*wmm_stat));
  Double_t wmm_weight = (1/(wmm_new*wmm_new+wmm_stat*wmm_stat))/(1/(wme_new*wme_new+wme_stat*wme_stat)+1/(wmm_new*wmm_new+wmm_stat*wmm_stat));

  cout << "wme " << wme_weight << " " << wmm_weight << endl;

  wml_stat=TMath::Sqrt(wme_weight*wme_weight*wme_stat*wme_stat+wmm_weight*wmm_weight*wmm_stat*wmm_stat);
  wml_sys=TMath::Sqrt(wme_weight*wme_weight*wme_new*wme_new+wmm_weight*wmm_weight*wmm_new*wmm_new);
  wml_lumi=wme_weight*wme_lumi+wmm_weight*wmm_lumi;
  wml = wme*wme_weight+wmm*wmm_weight;

  Double_t we_weight = (1/(we_new*we_new+we_stat*we_stat))/(1/(we_new*we_new+we_stat*we_stat)+1/(wm_new*wm_new+wm_stat*wm_stat));
  Double_t wm_weight = (1/(wm_new*wm_new+wm_stat*wm_stat))/(1/(we_new*we_new+we_stat*we_stat)+1/(wm_new*wm_new+wm_stat*wm_stat));

  cout << "we " << we_weight << " " << wm_weight << endl;

  wl_stat=TMath::Sqrt(we_weight*we_weight*we_stat*we_stat+wm_weight*wm_weight*wm_stat*wm_stat);
  wl_sys=TMath::Sqrt(we_weight*we_weight*we_new*we_new+wm_weight*wm_weight*wm_new*wm_new);
  wl_lumi=we_weight*we_lumi+wm_weight*wm_lumi;
  wl = we*we_weight+wm*wm_weight;

  Double_t ze_weight = (1/(ze_new*ze_new+ze_stat*ze_stat))/(1/(ze_new*ze_new+ze_stat*ze_stat)+1/(zm_new*zm_new+zm_stat*zm_stat));
  Double_t zm_weight = (1/(zm_new*zm_new+zm_stat*zm_stat))/(1/(ze_new*ze_new+ze_stat*ze_stat)+1/(zm_new*zm_new+zm_stat*zm_stat));

  cout << "ze " << ze_weight << " " << zm_weight << endl;

  zl_stat=TMath::Sqrt(ze_weight*ze_weight*ze_stat*ze_stat+zm_weight*zm_weight*zm_stat*zm_stat);
  zl_sys=TMath::Sqrt(ze_weight*ze_weight*ze_new*ze_new+zm_weight*zm_weight*zm_new*zm_new);
  zl_lumi=ze_weight*ze_lumi+zm_weight*zm_lumi;
  zl = ze*ze_weight+zm*zm_weight;

  Double_t wre_weight = (1/(wre_sys*wre_sys+wre_stat*wre_stat))/(1/(wre_sys*wre_sys+wre_stat*wre_stat)+1/(wrm_sys*wrm_sys+wrm_stat*wrm_stat));
  Double_t wrm_weight = (1/(wrm_sys*wrm_sys+wrm_stat*wrm_stat))/(1/(wre_sys*wre_sys+wre_stat*wre_stat)+1/(wrm_sys*wrm_sys+wrm_stat*wrm_stat));

  cout << "wre " << wre_weight << " " << wrm_weight << endl;

  wrl_stat=TMath::Sqrt(wre_weight*wre_weight*wre_stat*wre_stat+wrm_weight*wrm_weight*wrm_stat*wrm_stat);
  wrl_sys=TMath::Sqrt(wre_weight*wre_weight*wre_sys*wre_sys+wrm_weight*wrm_weight*wrm_sys*wrm_sys);
  wrl_lumi=wre_weight*wre_lumi+wrm_weight*wrm_lumi;
  wrl = wre*wre_weight+wrm*wrm_weight;

  Double_t wzre_weight = (1/(wzre_sys*wzre_sys+wzre_stat*wzre_stat))/(1/(wzre_sys*wzre_sys+wzre_stat*wzre_stat)+1/(wzrm_sys*wzrm_sys+wzrm_stat*wzrm_stat));
  Double_t wzrm_weight = (1/(wzrm_sys*wzrm_sys+wzrm_stat*wzrm_stat))/(1/(wzre_sys*wzre_sys+wzre_stat*wzre_stat)+1/(wzrm_sys*wzrm_sys+wzrm_stat*wzrm_stat));

  cout << "wzre " << wzre_weight << " " << wzrm_weight << endl;

  wzrl_stat=TMath::Sqrt(wzre_weight*wzre_weight*wzre_stat*wzre_stat+wzrm_weight*wzrm_weight*wzrm_stat*wzrm_stat);
  wzrl_sys=TMath::Sqrt(wzre_weight*wzre_weight*wzre_sys*wzre_sys+wzrm_weight*wzrm_weight*wzrm_sys*wzrm_sys);
  wzrl_lumi=wzre_weight*wzre_lumi+wzrm_weight*wzrm_lumi;
  wzrl = wzre*wzre_weight+wzrm*wzrm_weight;

  cout << endl;
  cout << "    MEASURED FIDUCIAL XSECS " << endl;
  cout << endl;
  cout << "W+(e)    " << wpe << " +/- " << wpe_stat << "_{stat} +/- " << wpe_new << "_{sys} +/- " << wpe_lumi << "_{lumi}" << endl;
  cout << "W+(m)    " << wpm << " +/- " << wpm_stat << "_{stat} +/- " << wpm_new << "_{sys} +/- " << wpm_lumi << "_{lumi}" << endl;
  cout << "W+       " << wpl << " +/- " << wpl_stat << "_{stat} +/- " << wpl_sys << "_{sys} +/- " << wpl_lumi << "_{lumi}" << endl;
  cout << endl;
  cout << "W-(e)    " << wme << " +/- " << wme_stat << "_{stat} +/- " << wme_new << "_{sys} +/- " << wme_lumi << "_{lumi}" << endl;
  cout << "W-(m)    " << wmm << " +/- " << wmm_stat << "_{stat} +/- " << wmm_new << "_{sys} +/- " << wmm_lumi << "_{lumi}" << endl;
  cout << "W-       " << wml << " +/- " << wml_stat << "_{stat} +/- " << wml_sys << "_{sys} +/- " << wml_lumi << "_{lumi}" << endl;
  cout << endl;
  cout << "W(e)     " << we << " +/- " << we_stat << "_{stat} +/- " << we_new << "_{sys} +/- " << we_lumi << "_{lumi}" << endl; 
  cout << "W(m)     " << wm << " +/- " << wm_stat << "_{stat} +/- " << wm_new << "_{sys} +/- " << wm_lumi << "_{lumi}" << endl; 
  cout << "W        " << wl << " +/- " << wl_stat << "_{stat} +/- " << wl_sys << "_{sys} +/- " << wl_lumi << "_{lumi}" << endl;
  cout << endl;
  cout << "Z(e)     " << ze << " +/- " << ze_stat << "_{stat} +/- " << ze_new << "_{sys} +/- " << ze_lumi << "_{lumi}" << endl;
  cout << "Z(m)     " << zm << " +/- " << zm_stat << "_{stat} +/- " << zm_new << "_{sys} +/- " << zm_lumi << "_{lumi}" << endl;
  cout << "Z        " << zl << " +/- " << zl_stat << "_{stat} +/- " << zl_sys << "_{sys} +/- " << zl_lumi << "_{lumi}" << endl;
  cout << endl;
  cout << "W+/W-(e) " << wre << " +/- " << wre_stat << "_{stat} +/- " << wre_sys << "_{sys} +/- " << wre_lumi << "_{lumi}" << endl;
  cout << "W+/W-(m) " << wrm << " +/- " << wrm_stat << "_{stat} +/- " << wrm_sys << "_{sys} +/- " << wrm_lumi << "_{lumi}" << endl;
  cout << "W+/W-    " << wrl << " +/- " << wrl_stat << "_{stat} +/- " << wrl_sys << "_{sys} +/- " << wrl_lumi << "_{lumi}" << endl;
  cout << endl;
  cout << "W/Z(e)   " << wzre << " +/- " << wzre_stat << "_{stat} +/- " << wzre_sys << "_{sys} +/- " << wzre_lumi << "_{lumi}" << endl;
  cout << "W/Z(m)   " << wzrm << " +/- " << wzrm_stat << "_{stat} +/- " << wzrm_sys << "_{sys} +/- " << wzrm_lumi << "_{lumi}" << endl;
  cout << "W/Z      " << wzrl << " +/- " << wzrl_stat << "_{stat} +/- " << wzrl_sys << "_{sys} +/- " << wzrl_lumi << "_{lumi}" << endl;
  cout << endl;
  cout << "    THEORETICAL FIDUCIAL XSECS" << endl;

  acc_wpe_unc*=acc_wpe/100;
  acc_wpm_unc*=acc_wpm/100;
  acc_wme_unc*=acc_wme/100;
  acc_wmm_unc*=acc_wmm/100;
  acc_we_unc*=acc_we/100;
  acc_wm_unc*=acc_wm/100;
  acc_ze_unc*=acc_ze/100;
  acc_zm_unc*=acc_zm/100;

  acc_wre=acc_wpe/acc_wme;
  acc_wrm=acc_wpm/acc_wmm;
  acc_wzre=acc_we/acc_ze;
  acc_wzrm=acc_wm/acc_zm;

  acc_wre_unc*=acc_wre/100;
  acc_wrm_unc*=acc_wrm/100;
  acc_wzre_unc*=acc_wzre/100;
  acc_wzrm_unc*=acc_wzrm/100;

  Double_t acc_wpl = wpe_weight*acc_wpe+wpm_weight*acc_wpm;
  Double_t acc_wml = wme_weight*acc_wme+wmm_weight*acc_wmm;
  Double_t acc_wl = we_weight*acc_we+wm_weight*acc_wm;
  Double_t acc_zl = ze_weight*acc_ze+zm_weight*acc_zm;

  cout << endl;
  Double_t tf_wpe_unc = theo_wp*acc_wpe*TMath::Sqrt(theo_wp_unc*theo_wp_unc/(theo_wp*theo_wp)+acc_wpe_unc*acc_wpe_unc/(acc_wpe*acc_wpe));
  cout << "W+(e)    " << theo_wp*acc_wpe << " +/- " << tf_wpe_unc << endl;
  Double_t tf_wme_unc = theo_wm*acc_wme*TMath::Sqrt(theo_wm_unc*theo_wm_unc/(theo_wm*theo_wm)+acc_wme_unc*acc_wme_unc/(acc_wme*acc_wme));
  cout << "W-(e)    " << theo_wm*acc_wme << " +/- " << tf_wme_unc << endl;
  Double_t tf_we_unc = theo_w*acc_we*TMath::Sqrt(theo_w_unc*theo_w_unc/(theo_w*theo_w)+acc_we_unc*acc_we_unc/(acc_we*acc_we));
  cout << "W(e)     " << theo_w*acc_we << " +/- " << tf_we_unc << endl;
  cout << endl;
  Double_t tf_wpm_unc = theo_wp*acc_wpm*TMath::Sqrt(theo_wp_unc*theo_wp_unc/(theo_wp*theo_wp)+acc_wpm_unc*acc_wpm_unc/(acc_wpm*acc_wpm));
  cout << "W+(m)    " << theo_wp*acc_wpm << " +/- " << tf_wpm_unc << endl;
  Double_t tf_wmm_unc = theo_wm*acc_wmm*TMath::Sqrt(theo_wm_unc*theo_wm_unc/(theo_wm*theo_wm)+acc_wmm_unc*acc_wmm_unc/(acc_wmm*acc_wmm));
  cout << "W-(m)    " << theo_wm*acc_wmm << " +/- " << tf_wmm_unc << endl;
  Double_t tf_wm_unc = theo_w*acc_wm*TMath::Sqrt(theo_w_unc*theo_w_unc/(theo_w*theo_w)+acc_wm_unc*acc_wm_unc/(acc_wm*acc_wm));
  cout << "W(m)     " << theo_w*acc_wm << " +/- " << tf_wm_unc << endl;
  cout << endl;
  cout << "W+(l)    " << theo_wp*acc_wpl << " +/- " << TMath::Sqrt(tf_wpe_unc*tf_wpe_unc*wpe_weight*wpe_weight+tf_wpm_unc*tf_wpm_unc*wpm_weight*wpm_weight) << endl;
  cout << "W-(l)    " << theo_wm*acc_wml << " +/- " << TMath::Sqrt(tf_wme_unc*tf_wme_unc*wme_weight*wme_weight+tf_wmm_unc*tf_wmm_unc*wmm_weight*wmm_weight) << endl;
  cout << "W(l)     " << theo_w*acc_wl << " +/- " << TMath::Sqrt(tf_we_unc*tf_we_unc*we_weight*we_weight+tf_wm_unc*tf_wm_unc*wm_weight*wm_weight) << endl;
  cout << endl;
  Double_t tf_ze_unc = theo_z*acc_ze*TMath::Sqrt(theo_z_unc*theo_z_unc/(theo_z*theo_z)+acc_ze_unc*acc_ze_unc/(acc_ze*acc_ze));
  cout << "Z(e)     " << theo_z*acc_ze << " +/- " << tf_ze_unc << endl;
  Double_t tf_zm_unc = theo_z*acc_zm*TMath::Sqrt(theo_z_unc*theo_z_unc/(theo_z*theo_z)+acc_zm_unc*acc_zm_unc/(acc_zm*acc_zm));
  cout << "Z(m)     " << theo_z*acc_zm << " +/- " << tf_zm_unc << endl;
  cout << "Z(l)     " << theo_z*acc_zl << " +/- " << TMath::Sqrt(tf_ze_unc*tf_ze_unc*ze_weight*ze_weight+tf_zm_unc*tf_zm_unc*zm_weight*zm_weight) << endl;
  cout << endl;
  Double_t tf_wre = theo_wr*acc_wpe/acc_wme;
  Double_t tf_wre_unc = theo_wr*acc_wpe/acc_wme*TMath::Sqrt(theo_wr_unc*theo_wr_unc/(theo_wr*theo_wr)+acc_wre_unc*acc_wre_unc/(acc_wre*acc_wre));
  cout << "W+/W-(e) " << theo_wr*acc_wpe/acc_wme << " +/- " << tf_wre_unc << endl;
  Double_t tf_wrm = theo_wr*acc_wpm/acc_wmm;
  Double_t tf_wrm_unc = theo_wr*acc_wpm/acc_wmm*TMath::Sqrt(theo_wr_unc*theo_wr_unc/(theo_wr*theo_wr)+acc_wrm_unc*acc_wrm_unc/(acc_wrm*acc_wrm));
  cout << "W+/W-(m) " << theo_wr*acc_wpm/acc_wmm << " +/- " << tf_wrm_unc << endl;
  cout << "W+/W-(l) " << tf_wre*wre_weight+tf_wrm*wrm_weight << " +/- " << TMath::Sqrt(tf_wre_unc*tf_wre_unc*wre_weight*wre_weight+tf_wrm_unc*tf_wrm_unc*wrm_weight*wrm_weight) << endl;
  cout << endl;
  Double_t tf_wzre = theo_wzr*acc_we/acc_ze;
  Double_t tf_wzre_unc = theo_wzr*acc_we/acc_ze*TMath::Sqrt(theo_wzr_unc*theo_wzr_unc/(theo_wzr*theo_wzr)+acc_wzre_unc*acc_wzre_unc/(acc_wzre*acc_wzre));
  cout << "W/Z(e)   " << theo_wzr*acc_we/acc_ze << " +/- " << tf_wzre_unc << endl;
  Double_t tf_wzrm = theo_wzr*acc_wm/acc_zm;
  Double_t tf_wzrm_unc = theo_wzr*acc_wm/acc_zm*TMath::Sqrt(theo_wzr_unc*theo_wzr_unc/(theo_wzr*theo_wzr)+acc_wzrm_unc*acc_wzrm_unc/(acc_wzrm*acc_wzrm));
  cout << "W/Z(m)   " << theo_wzr*acc_wm/acc_zm << " +/- " << tf_wzrm_unc << endl;
  cout << "W/Z(l)   " << tf_wzre*wzre_weight+tf_wzrm*wzrm_weight << " +/- " << TMath::Sqrt(tf_wzre_unc*tf_wzre_unc*wzre_weight*wzre_weight+tf_wzrm_unc*tf_wzrm_unc*wzrm_weight*wzrm_weight) << endl;
  cout << endl;

}
