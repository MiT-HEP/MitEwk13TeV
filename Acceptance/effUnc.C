#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMinuit.h"

// partial directory names for each of the options
TString inDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/";
TString outDir = "/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Acceptance/";
// function to do some division

void writeTable(TString filename, int s, string lep, vector<vector<double>> vals);

void writeUncTable(TString filename, vector<double> e, vector<double> m);
void writeMegaUncTable(TString filename, vector<double> e13, vector<double> m13, vector<double> e5, vector<double> m5, vector<double> e135, vector<double> m135);

void writeAccTable(TString filename, vector<double>& vAcc13, vector<double>& vAcc13D);

vector<double> combineNjets(vector<double>& v0j, vector<double>& v1j, vector<double>& v2j);
vector<double> makeRatios(vector<double>& num, vector<double>& dnm);
vector<vector<double>> makeRatios(vector<vector<double>>& num, vector<vector<double>>& dnm);

void readFileM(TString filename, vector<double> &all, vector<double> &sel, vector<double> &sta);
void readFileE(TString filename, vector<double> &sel);

vector<double> computeRow(int r, vector<vector<double>> vals);

void effUnc(){
  // This is horrible, but too bad.
  vector<vector<double>> vWmpJets_all, vWmmJets_all;
  vector<vector<double>> vWmpJets_sel, vWmmJets_sel;
  vector<vector<double>> vWmpJets_sta, vWmmJets_sta;
  vector<vector<double>> vWepJets_gsf, vWemJets_gsf;
  for(int i = 0; i < 3; ++i){
    vector<double> empty;
    vWmpJets_all.push_back(empty);
    vWmmJets_all.push_back(empty);
    
    vWmpJets_sel.push_back(empty);
    vWmmJets_sel.push_back(empty);
    
    vWmpJets_sta.push_back(empty);
    vWmmJets_sta.push_back(empty);
    
    vWepJets_gsf.push_back(empty);
    vWemJets_gsf.push_back(empty);
  }
  
  
  vector<vector<double>> vMu13_sel, vMu13_sta, vMu13_all, vEle13_gsf;
  vector<vector<double>> vMu5_sel, vMu5_sta, vMu5_all, vEle5_gsf;
  // w+, w-, w, w+/w-, z, w+/z, w-/z, w/z
  for(int i = 0; i < 8; ++i){
    vector<double> empty;
    vMu13_sel.push_back(empty);
    vMu13_sta.push_back(empty);
    vMu13_all.push_back(empty);
    vEle13_gsf.push_back(empty);
    
    vMu5_sel.push_back(empty);
    vMu5_sta.push_back(empty);
    vMu5_all.push_back(empty);
    vEle5_gsf.push_back(empty);
  }
  
  
  vector<double> vWmp13_sel, vWmp13_sta;
  vector<double> vWmm13_sel, vWmm13_sta;

  vector<double> vWep13_sel, vWep13_sta;
  vector<double> vWem13_sel, vWem13_sta;
  
  vector<double> vZmm13_sel, vZmm13_sta;
  vector<double> vZee13_gsf;
  
  vector<double> vWmp5_sel, vWmp5_sta;
  vector<double> vWmm5_sel, vWmm5_sta;

  vector<double> vWep5_gsf;
  vector<double> vWem5_gsf;
  
  vector<double> vZmm5_sel, vZmm5_sta;
  vector<double> vZee5_gsf;
  
  std::cout << "load 13 TeV W+ " << std::endl;
  readFileM("Wmp_0j_13TeV_wEff/sel_nums_only.txt", vWmpJets_all[0], vWmpJets_sel[0], vWmpJets_sta[0]);
  readFileM("Wmp_1j_13TeV_wEff/sel_nums_only.txt", vWmpJets_all[1], vWmpJets_sel[1], vWmpJets_sta[1]);
  readFileM("Wmp_2j_13TeV_wEff/sel_nums_only.txt", vWmpJets_all[2], vWmpJets_sel[2], vWmpJets_sta[2]);
  std::cout << "load 13 TeV W- " << std::endl;
  readFileM("Wmm_0j_13TeV_wEff/sel_nums_only.txt", vWmmJets_all[0], vWmmJets_sel[0], vWmmJets_sta[0]);
  readFileM("Wmm_1j_13TeV_wEff/sel_nums_only.txt", vWmmJets_all[1], vWmmJets_sel[1], vWmmJets_sta[1]);
  readFileM("Wmm_2j_13TeV_wEff/sel_nums_only.txt", vWmmJets_all[2], vWmmJets_sel[2], vWmmJets_sta[2]);
  std::cout << "load 13 TeV Z " << std::endl;
  readFileM("Zmm_13TeV_wEff/sel_nums_only.txt",  vMu13_all[4], vMu13_sel[4], vMu13_sta[4]);
  
  readFileE("Wep_0j_13TeV_wEff/sel_nums_only.txt", vWepJets_gsf[0]);
  readFileE("Wep_1j_13TeV_wEff/sel_nums_only.txt", vWepJets_gsf[1]);
  readFileE("Wep_2j_13TeV_wEff/sel_nums_only.txt", vWepJets_gsf[2]);
  std::cout << "load 13 TeV W- " << std::endl;
  readFileE("Wem_0j_13TeV_wEff/sel_nums_only.txt", vWemJets_gsf[0]);
  readFileE("Wem_1j_13TeV_wEff/sel_nums_only.txt", vWemJets_gsf[1]);
  readFileE("Wem_2j_13TeV_wEff/sel_nums_only.txt", vWemJets_gsf[2]);
  readFileE("Zee_13TeV_wEff/sel_nums_only.txt", vEle13_gsf[4]);
  
  std::cout << "load 5 TeV W+ " << std::endl;
  readFileM("Wmp_0j_13TeV_wEff/sel_nums_only.txt", vMu5_all[0], vMu5_sta[0], vMu5_sel[0]);
  std::cout << "load 5 TeV W- " << std::endl;
  readFileM("Wmm_0j_13TeV_wEff/sel_nums_only.txt", vMu5_all[1], vMu5_sta[1], vMu5_sel[1]);
  std::cout << "load 5 TeV Z " << std::endl;
  readFileM("Zmm_13TeV_wEff/sel_nums_only.txt", vMu5_all[4], vMu5_sel[4], vMu5_sta[4]);
  
  readFileE("Wep_0j_13TeV_wEff/sel_nums_only.txt",vEle5_gsf[0]);
  readFileE("Wem_0j_13TeV_wEff/sel_nums_only.txt", vEle5_gsf[1]);
  readFileE("Zee_13TeV_wEff/sel_nums_only.txt", vEle5_gsf[4]);
  
  cout << "Done reading files " << endl;

  // new contain inclusive acceptance: main / fsr / mc / bkg / tag
  vMu13_all[0] = combineNjets(vWmpJets_all[0], vWmpJets_all[1], vWmpJets_all[2]);
  vMu13_all[1] = combineNjets(vWmmJets_all[0], vWmmJets_all[1], vWmmJets_all[2]);
  
  vMu13_sel[0] = combineNjets(vWmpJets_sel[0], vWmpJets_sel[1], vWmpJets_sel[2]);
  vMu13_sel[1] = combineNjets(vWmmJets_sel[0], vWmmJets_sel[1], vWmmJets_sel[2]);
  
  vMu13_sta[0] = combineNjets(vWmpJets_sta[0], vWmpJets_sta[1], vWmpJets_sta[2]);
  vMu13_sta[1] = combineNjets(vWmmJets_sta[0], vWmmJets_sta[1], vWmmJets_sta[2]);
  
  
  vEle13_gsf[0] = combineNjets(vWepJets_gsf[0], vWepJets_gsf[1], vWepJets_gsf[2]);
  vEle13_gsf[1] = combineNjets(vWemJets_gsf[0], vWemJets_gsf[1], vWemJets_gsf[2]);
  
  cout << "Done Combine Jets " << endl;
 
  for(int i = 0; i < vMu13_all[0].size(); ++i){
    vMu13_all[2].push_back((vMu13_all[0][i]+vMu13_all[1][i])/2.0);
    vMu13_sel[2].push_back((vMu13_sel[0][i]+vMu13_sel[1][i])/2.0);
    vMu13_sta[2].push_back((vMu13_sta[0][i]+vMu13_sta[1][i])/2.0);
    
    vEle13_gsf[2].push_back((vEle13_gsf[0][i]+vEle13_gsf[1][i])/2.0);
    
    vMu5_all[2].push_back((vMu5_all[0][i]+vMu5_all[1][i])/2.0);
    vMu5_sel[2].push_back((vMu5_sel[0][i]+vMu5_sel[1][i])/2.0);
    vMu5_sta[2].push_back((vMu5_sta[0][i]+vMu5_sta[1][i])/2.0);
    
    vEle5_gsf[2].push_back((vEle5_gsf[0][i]+vEle5_gsf[1][i])/2.0);
    
  }

  cout << "Start Ratios...  " << endl;
  cout << "W 13 all  " << endl;
  vMu13_all[3]   = makeRatios(vMu13_all[0],vMu13_all[1]);
  vMu13_all[5]   = makeRatios(vMu13_all[0],vMu13_all[4]);
  vMu13_all[6]   = makeRatios(vMu13_all[1],vMu13_all[4]);
  vMu13_all[7]   = makeRatios(vMu13_all[2],vMu13_all[4]);
  cout << "W 13 sel  " << endl;
  vMu13_sel[3]   = makeRatios(vMu13_sel[0],vMu13_sel[1]);
  vMu13_sel[5]   = makeRatios(vMu13_sel[0],vMu13_sel[4]);
  vMu13_sel[6]   = makeRatios(vMu13_sel[1],vMu13_sel[4]);
  vMu13_sel[7]   = makeRatios(vMu13_sel[2],vMu13_sel[4]);
  cout << "W 13 sta  " << endl;
  vMu13_sta[3]   = makeRatios(vMu13_sta[0],vMu13_sta[1]);
  vMu13_sta[5]   = makeRatios(vMu13_sta[0],vMu13_sta[4]);
  vMu13_sta[6]   = makeRatios(vMu13_sta[1],vMu13_sta[4]);
  vMu13_sta[7]   = makeRatios(vMu13_sta[2],vMu13_sta[4]);
  cout << "W 13 gsf  " << endl;
  vEle13_gsf[3]   = makeRatios(vEle13_gsf[0],vEle13_gsf[1]);
  vEle13_gsf[5]   = makeRatios(vEle13_gsf[0],vEle13_gsf[4]);
  vEle13_gsf[6]   = makeRatios(vEle13_gsf[1],vEle13_gsf[4]);
  vEle13_gsf[7]   = makeRatios(vEle13_gsf[2],vEle13_gsf[4]);
  
  cout << "W 5 all  " << endl;
  vMu5_all[3]   = makeRatios(vMu5_all[0],vMu5_all[1]);
  vMu5_all[5]   = makeRatios(vMu5_all[0],vMu5_all[4]);
  vMu5_all[6]   = makeRatios(vMu5_all[1],vMu5_all[4]);
  vMu5_all[7]   = makeRatios(vMu5_all[2],vMu5_all[4]);
  cout << "W 5 sel  " << endl;
  vMu5_sel[3]   = makeRatios(vMu5_sel[0],vMu5_sel[1]);
  vMu5_sel[5]   = makeRatios(vMu5_sel[0],vMu5_sel[4]);
  vMu5_sel[6]   = makeRatios(vMu5_sel[1],vMu5_sel[4]);
  vMu5_sel[7]   = makeRatios(vMu5_sel[2],vMu5_sel[4]);
  cout << "W 5 sta  " << endl;
  vMu5_sta[3]   = makeRatios(vMu5_sta[0],vMu5_sta[1]);
  vMu5_sta[5]   = makeRatios(vMu5_sta[0],vMu5_sta[4]);
  vMu5_sta[6]   = makeRatios(vMu5_sta[1],vMu5_sta[4]);
  vMu5_sta[7]   = makeRatios(vMu5_sta[2],vMu5_sta[4]);
  cout << "W 5 gsf  " << endl;
  vEle5_gsf[3]   = makeRatios(vEle5_gsf[0],vEle5_gsf[1]);
  vEle5_gsf[5]   = makeRatios(vEle5_gsf[0],vEle5_gsf[4]);
  vEle5_gsf[6]   = makeRatios(vEle5_gsf[1],vEle5_gsf[4]);
  vEle5_gsf[7]   = makeRatios(vEle5_gsf[2],vEle5_gsf[4]);
  
  // do double ratios...
  
  // compute rows/columns?
  vector<vector<double>> vMu13_AllUnc, vEle13_AllUnc;
  vector<vector<double>> vMu13_StaUnc;
  vector<vector<double>> vMu5_AllUnc, vEle5_AllUnc;
  for(int i = 0; i < 4; i++){
    vector<double> empty;
    vMu13_StaUnc.push_back(empty);
    vMu13_AllUnc.push_back(empty);
    vEle13_AllUnc.push_back(empty);
    vMu5_AllUnc.push_back(empty);
    vEle5_AllUnc.push_back(empty);
  }
  
  cout << "Compute Rows Mu13 " << endl;
  vMu13_StaUnc[0] = computeRow(2,vMu13_sta); // row number, vector
  vMu13_StaUnc[1] = computeRow(3,vMu13_sta); // row number, vector
  vMu13_StaUnc[2] = computeRow(4,vMu13_sta); // row number, vector
  vMu13_StaUnc[3] = computeRow(5,vMu13_sta); // row number, vector
  
  vMu13_AllUnc[0] = computeRow(2,vMu13_all); // row number, vector
  vMu13_AllUnc[1] = computeRow(3,vMu13_all); // row number, vector
  vMu13_AllUnc[2] = computeRow(4,vMu13_all); // row number, vector
  vMu13_AllUnc[3] = computeRow(5,vMu13_all); // row number, vector
  
  cout << "Compute Rows Ele 13  " << endl;
  vEle13_AllUnc[0] = computeRow(2,vEle13_gsf); // row number, vector
  vEle13_AllUnc[1] = computeRow(3,vEle13_gsf); // row number, vector
  vEle13_AllUnc[2] = computeRow(4,vEle13_gsf); // row number, vector
  vEle13_AllUnc[3] = computeRow(5,vEle13_gsf); // row number, vector
  
  cout << "Compute Rows Mu5 " << endl;
  vMu5_AllUnc[0] = computeRow(2,vMu5_all); // row number, vector
  vMu5_AllUnc[1] = computeRow(3,vMu5_all); // row number, vector
  vMu5_AllUnc[2] = computeRow(4,vMu5_all); // row number, vector
  vMu5_AllUnc[3] = computeRow(5,vMu5_all); // row number, vector
  
  cout << "Compute Rows Ele 5  " << endl;
  vEle5_AllUnc[0] = computeRow(2,vEle5_gsf); // row number, vector
  vEle5_AllUnc[1] = computeRow(3,vEle5_gsf); // row number, vector
  vEle5_AllUnc[2] = computeRow(4,vEle5_gsf); // row number, vector
  vEle5_AllUnc[3] = computeRow(5,vEle5_gsf); // row number, vector
    
  // Write Tables
  cout << "Write Tables " << endl;
  writeTable("TTTTTTTT_13TeV_MuEffUnc_Summary", 13, "muon", vMu13_AllUnc);
  writeTable("TTTTTTTT_13TeV_MuEffUnc_Summary_Standalone", 13, "muon", vMu13_StaUnc);
  writeTable("TTTTTTTT_13TeV_EleEffUnc_Summary", 13, "electron", vEle13_AllUnc);
 
 
  
  // writeSummaryTable("TESTTEST_13TeV_MuEff_SysUnc_Summary", 13, "muon", vWmp13, vWmm13, vZmm13)
  // writeSummaryTable("TESTTEST_5TeV_EleEff_SysUnc_Summary",  5, "muon", vWmp5 , vWmm5  , vZmm5)
  
  // writeAccTable("TESTTEST_13TeV_Acc.txt",vAccEle13,vAccMu13);
  // writeAccTable("TESTTEST_5TeV_Acc.txt",vAccEle5,vAccMu5);
  
  // // write each of the uncertainty tables
  // writeUncTable("TESTTEST_5TeV.txt",vUncEle5,vUncMu5);
  // writeUncTable("TESTTEST_13TeV.txt",vUncEle13,vUncMu13);
  // writeUncTable("TESTTEST_13to5TeV.txt",vUncEle135,vUncMu135);
  // writeMegaUncTable("TESTTEST_Mega.txt",vUncEle13,vUncMu13,vUncEle5,vUncMu5,vUncEle135,vUncMu135);

}


void readFileM(TString filename, vector<double> &all, vector<double> &sel, vector<double> &sta){
  double acc1,acc2,acc3;
  char infilename[250];
  sprintf(infilename,"%s/%s",inDir.Data(),filename.Data());
  ifstream myfile(infilename);
  assert(myfile);
  myfile>>acc1>>acc2;
  // cout << "acc1 " << acc1 << "  acc2 " << acc2 << endl; 
  sel.push_back(acc1);
  sel.push_back(acc1+acc2);
  sta.push_back(acc1);
  sta.push_back(acc1+acc2);
  all.push_back(acc1);
  all.push_back(acc1+acc2);
  while(myfile>>acc3>>acc2>>acc1){
    all.push_back(acc3);
    sel.push_back(acc2);
    sta.push_back(acc1);
  }
  // cout << "Done reading File!" << endl;
}

void readFileE(TString filename, vector<double> &sel){
  double acc1,acc2;
  char infilename[250];
  sprintf(infilename,"%s/%s",inDir.Data(),filename.Data());
  ifstream myfile(infilename);
  assert(myfile);
  myfile>>acc1>>acc2;
  sel.push_back(acc1);
  sel.push_back(acc1+acc2);
  while(myfile>>acc2){
    sel.push_back(acc2);
  }
  // cout << "Done reading File!" << endl;
}

vector<double> combineNjets(vector<double>& v0j, vector<double>& v1j, vector<double>& v2j){
  vector<double> result;
  cout << "Combine N Jets files " << endl;
  // get x-sections for 0,1,2 jet
  double avg = 0;
  for(uint i = 0; i < v0j.size(); ++i){
    avg +=51410*v0j[i] + 8417*v1j[i] + 3306*v2j[i];
    avg = avg/(51410+8417+3306);
    result.push_back(avg);
    cout << "avg = " << avg << endl;
  }
  cout << "Done combining N Jets files " << endl;
  return result;
}

vector<double> makeRatios(vector<double>& num, vector<double>& dnm){
  vector<double> result;
  cout << "Make Ratios" << endl;
  // get x-sections for 0,1,2 jet
  double rat = 0;
  for(uint i = 0; i < num.size(); ++i){
    rat = num[i]/dnm[i];
    cout << rat << endl;
    result.push_back(rat);
  }
  cout << "Done ratios " << endl;
  return result;
}

vector<double> sumMu(vector<double>& sel, vector<double>& sta){
  vector<double> result;
  // get x-sections for 0,1,2 jet
  for(uint i = 0; i < sel.size(); ++i){
    result.push_back(sel[i] + sta[i]);
  }
  return result;
}

vector<double> computeRow(int r, vector<vector<double>> vals){
  cout << "Compute Row Begin" << endl;
  vector<double> result;
  cout << "Size of megavec " << vals.size() << endl;
  for(int i = 0; i < vals.size(); ++i){
    cout << "Yes Hello" << endl;
    cout << vals[i][0] << endl;
    double pctDiff = 100*(vals[i][0]-vals[i][r])/vals[i][0];
    cout << vals[i][0] << " " << vals[i][r] << " " << pctDiff << endl;
    result.push_back(pctDiff);
  }
  cout << "Compute Row End" << endl;
  return result;
}

void writeTable(TString filename, int s, string lep, vector<vector<double>> vals){
  
  char txtfname[250];
  sprintf(txtfname,"%s/%s.tex",outDir.Data(),filename.Data());
  // ofstream txtfile;
  // txtfile.open(txtfname);
  FILE *txtfile;
  txtfile = fopen(txtfname,"w");
  fprintf(txtfile,"\\begin{table}%[htbp]\n");
  fprintf(txtfile,"\\begin{center}\n");
  fprintf(txtfile,"\\scalebox{0.7}{\n");
  fprintf(txtfile,"\\begin{tabular}{ccccccccc}\n");
  fprintf(txtfile,"\\hline\n");
  fprintf(txtfile,"Source & $W^+$& $W^-$ & $W$ & $W^+/W^-$ & $Z$ & $W^+/Z$&$W^-/Z$ &$W/Z$  \\\\\n");
  fprintf(txtfile,"\\hline \\hline\n");
  fprintf(txtfile,"FSR [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f\\\\\n", 
         vals[0][0],vals[0][1],vals[0][2],vals[0][3],vals[0][4],vals[0][5],vals[0][6],vals[0][7]);
  fprintf(txtfile,"MC [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f\\\\\n", 
         vals[1][0],vals[1][1],vals[1][2],vals[1][3],vals[1][4],vals[1][5],vals[1][6],vals[1][7]);
  fprintf(txtfile,"Background [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f\\\\\n", 
         vals[2][0],vals[2][1],vals[2][2],vals[2][3],vals[2][4],vals[2][5],vals[2][6],vals[2][7]);
  fprintf(txtfile,"Tag p_T [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f\\\\\n", 
         vals[3][0],vals[3][1],vals[3][2],vals[3][3],vals[3][4],vals[3][5],vals[3][6],vals[3][7]);
  fprintf(txtfile,"\\hline\n");
  fprintf(txtfile,"\\end{tabular}}\n");
  fprintf(txtfile,"\\end{center}\n");
  fprintf(txtfile,"\\caption{Uncertainties on the lepton efficiency scale factors for the %s channel in %d TeV}\n",lep.c_str(),s);
  fprintf(txtfile,"\\label{tab:Eff:Unc:%s:summary:%dTeV}\n",lep.c_str(), s);
  fprintf(txtfile,"\\end{table}\n");

  fclose(txtfile);
  
}

void writeAccTable(TString filename, vector<double>& e, vector<double>& m){
  char txtfname[250];
  sprintf(txtfname,"%s/%s.tex",outDir.Data(),filename.Data());
  // ofstream txtfile;
  // txtfile.open(txtfname);
  
  FILE *txtfile;
  txtfile = fopen(txtfname,"w");
  fprintf(txtfile,"\\begin{center}\n");
  fprintf(txtfile,"\\scalebox{0.8}{\n");
  fprintf(txtfile,"\\begin{tabular}{|c|c|c|c|c|}\n");
  fprintf(txtfile,"\\\\hline\n");
  fprintf(txtfile,"Process & $A_{Gen}(\\mathrm{Post-FSR})$ & $A_{Gen}(\\mathrm{Dressed}$) \\\\");
  fprintf(txtfile,"\\hline \\hline\n");
  fprintf(txtfile,"$W\\rightarrow e^+\\nu$     & %.3f & %.3f \\\\\n",e[0],  e[4]);
  fprintf(txtfile,"$W\\rightarrow e^-\\nu$     & %.3f & %.3f \\\\\n",e[1],  e[5]);
  fprintf(txtfile,"$W\\rightarrow e\\nu$       & %.3f & %.3f \\\\\n",e[2],  e[6]);
  fprintf(txtfile,"$Z\\rightarrow ee$          & %.3f & %.3f \\\\\n",e[3],  e[7]);
  fprintf(txtfile,"\\hline\n" );
  fprintf(txtfile,"$W\\rightarrow \\mu^+\\nu$   & %.3f & %.3f \\\\\n",m[0],  m[4]);
  fprintf(txtfile,"$W\\rightarrow \\mu^-\\nu$   & %.3f & %.3f \\\\\n",m[1],  m[5]);
  fprintf(txtfile,"$W\\rightarrow \\mu\\nu$     & %.3f & %.3f \\\\\n",m[2],  m[6]);
  fprintf(txtfile,"$Z\\rightarrow \\mu\\mu$     & %.3f & %.3f \\\\\n",m[3],  m[7]);
  fprintf(txtfile,"\\hline\n" );
  fprintf(txtfile,"\\end{tabular} }\n" );
  fprintf(txtfile,"\\end{center}\n" );
  // fprintf(txtfile,"*" );
  fclose(txtfile);
  // txtfile.close();
}


void writeSummaryTable(TString filename, int s, string lep, vector<double> wp, vector<double> wm, vector<double> z){
  char txtfname[250];
  sprintf(txtfname,"%s/%s.tex",outDir.Data(),filename.Data());
  // ofstream txtfile;
  // txtfile.open(txtfname);
  FILE *txtfile;
  txtfile = fopen(txtfname,"w");
  fprintf(txtfile,"\\begin{table}%[htbp]\n");
  fprintf(txtfile,"\\begin{center}\n");
  fprintf(txtfile,"\\scalebox{0.7}{\n");
  fprintf(txtfile,"\\begin{tabular}{ccccccccc}\n");
  fprintf(txtfile,"\\hline\n");
  fprintf(txtfile,"Source & $W^+$& $W^-$ & $W$ & $W^+/W^-$ & $Z$ & $W^+/Z$&$W^+/Z$ &$W/Z$  \\\n");
  fprintf(txtfile,"\\hline \\hline\n");
  fprintf(txtfile,"Binning [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f\\\n", 
         wp[0], wm[0], (wm[0]+wp[0])/2, abs(wm[0]-wp[0]), z[0], abs(wp[0]-z[0]), abs(wm[0]-z[0]), abs((wp[0]+wm[0])/2-z[0]));
  fprintf(txtfile,"Signal Shape [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f \\\n");
  fprintf(txtfile,"Background Shape [\\%] & %.3f  & %.3f & %.3f & %.3f & %.3f& %.3f& %.3f& %.3f  \\\n");
  fprintf(txtfile,"\\hline\n");
  fprintf(txtfile,"\\end{tabular}}\n");
  fprintf(txtfile,"\\end{center}\n");
  fprintf(txtfile,"\\caption{Summary of the propagated %s efficiency systematic uncertainties at %d TeV.}\n",lep.c_str(),s);
  fprintf(txtfile,"\\label{tab:Eff:Unc:%s:summary:%dTeV}\n",lep.c_str(), s);
  fprintf(txtfile,"\\end{table}\n");

  fclose(txtfile);
}

void writeUncTable(TString filename, vector<double> e, vector<double> m){
  char txtfname[250];
  sprintf(txtfname,"%s/%s.tex",outDir.Data(),filename.Data());
  // ofstream txtfile;
  // txtfile.open(txtfname);
  FILE *txtfile;
  txtfile = fopen(txtfname,"w");
  fprintf(txtfile,"\\begin{center}\n");
  fprintf(txtfile,"\\scalebox{0.8}{\n");
  fprintf(txtfile,"\\begin{tabular}{|c|c|c|c|c|}\n");
  fprintf(txtfile,"\\hline\n");
  fprintf(txtfile,"Process &\\multicolumn{2}{|c|}{QCD[\\%]}&\\multicolumn{2}{|c|}{PDF[\\%]}\\\\\\cline{2-5}\n");
  fprintf(txtfile,"( TeV) & Ele & Mu & Ele & Mu \\\\\n");
  fprintf(txtfile,"\\hline \\hline\n");
  fprintf(txtfile,"$W^+$     & %.3f & %.3f & %.3f & %.3f\\\\\n",e[8],  m[8],  e[0],m[0] );
  fprintf(txtfile,"$W^-$     & %.3f & %.3f & %.3f & %.3f\\\\\n",e[9],  m[9],  e[1],m[1] );
  fprintf(txtfile,"$W$       & %.3f & %.3f & %.3f & %.3f\\\\\n",e[10], m[10], e[2],m[2] );
  fprintf(txtfile,"$Z$       & %.3f & %.3f & %.3f & %.3f\\\\\n",e[11], m[11], e[3],m[3] );
  fprintf(txtfile,"$W^+/W^-$ & %.3f & %.3f & %.3f & %.3f\\\\\n",e[12], m[12], e[4],m[4] );
  fprintf(txtfile,"$W^+/Z$   & %.3f & %.3f & %.3f & %.3f\\\\\n",e[13], m[13], e[5],m[5] );
  fprintf(txtfile,"$W^-/Z$   & %.3f & %.3f & %.3f & %.3f\\\\\n",e[14], m[14], e[6],m[6] );
  fprintf(txtfile,"$W/Z$     & %.3f & %.3f & %.3f & %.3f\\\\\n",e[15], m[15], e[7],m[7] );
  fprintf(txtfile,"\\hline\n" );
  fprintf(txtfile,"\\end{tabular} }\n" );
  fprintf(txtfile,"\\end{center}\n" );
  // fprintf(txtfile,"*" );
  fclose(txtfile);
  // txtfile.close();
}

void writeMegaUncTable(TString filename, vector<double> e13, vector<double> m13, vector<double> e5, vector<double> m5, vector<double> e135, vector<double> m135){
  char txtfname[250];
  sprintf(txtfname,"%s/%s.tex",outDir.Data(),filename.Data());
  // ofstream txtfile;
  // txtfile.open(txtfname);
  FILE *txtfile;
  txtfile = fopen(txtfname,"w");
  fprintf(txtfile,"\\begin{center}\n");
  fprintf(txtfile,"\\scalebox{0.8}{\n");
  fprintf(txtfile,"\\begin{tabular}{|c||c|c|c|c||c|c|c|c|}\n");
  fprintf(txtfile,"\\hline\n");
  fprintf(txtfile," &\\multicolumn{4}{|c|}{13 TeV}&\\multicolumn{4}{|c|}{5 TeV}\\\\\\cline{2-9}\n");
  fprintf(txtfile,"Process &\\multicolumn{2}{|c|}{QCD[\\%]}&\\multicolumn{2}{|c||}{PDF[\\%]}&\\multicolumn{2}{|c|}{QCD[\\%]}&\\multicolumn{2}{|c|}{PDF[\\%]}\\\\\\cline{2-9}\n");
  fprintf(txtfile,"  & Ele & Mu & Ele & Mu & Ele & Mu & Ele & Mu \\\\\n");
  fprintf(txtfile,"\\hline \\hline\n");
  fprintf(txtfile,"$W^+$     & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[8],  m13[8],  e13[0],m13[0],e5[8],  m5[8],  e5[0],m5[0] );
  fprintf(txtfile,"$W^-$     & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[9],  m13[9],  e13[1],m13[1],e5[9],  m5[9],  e5[1],m5[1] );
  fprintf(txtfile,"$W$       & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[10], m13[10], e13[2],m13[2],e5[10], m5[10], e5[2],m5[2] );
  fprintf(txtfile,"$Z$       & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[11], m13[11], e13[3],m13[3],e5[11], m5[11], e5[3],m5[3] );
  fprintf(txtfile,"$W^+/W^-$ & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[12], m13[12], e13[4],m13[4],e5[12], m5[12], e5[4],m5[4] );
  fprintf(txtfile,"$W^+/Z$   & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[13], m13[13], e13[5],m13[5],e5[13], m5[13], e5[5],m5[5] );
  fprintf(txtfile,"$W^-/Z$   & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[14], m13[14], e13[6],m13[6],e5[14], m5[14], e5[6],m5[6] );
  fprintf(txtfile,"$W/Z$     & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n",e13[15], m13[15], e13[7],m13[7],e5[15], m5[15], e5[7],m5[7] );
  fprintf(txtfile,"\\hline\n" );
  fprintf(txtfile,"\\end{tabular} }\n" );
  fprintf(txtfile,"\\end{center}\n" );
  // fprintf(txtfile,"*" );
  fclose(txtfile);
}
