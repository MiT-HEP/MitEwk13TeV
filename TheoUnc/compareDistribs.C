#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <TH1.h>

//#include "MitStyleRemix.hh"
#include "compare.hh"

using namespace std;

#endif

void compareWmm();
void compareWpm();
void compareZmm();

void compareDistribs() {

  compareWmm();
  compareWpm();
  compareZmm();

}

void compareWmm() {

  vector<Config> conf_list;

  conf_list.push_back(Config("baseComp-05-06/CT10-wmm-pythia6-powheg1-v2-bacon.root", "ct10", "CT10 + POWHEG-1 + PYTHIA 6"));
  conf_list.push_back(Config("baseComp-05-06/wmm-pythia6-powheg1.root", "pow1py6", "POWHEG-1 + PYTHIA 6"));
  conf_list.push_back(Config("baseComp-05-06/wmm-pythia8-powheg1.root", "pow1py8", "POWHEG-1 + PYTHIA 8"));
  conf_list.push_back(Config("baseComp-05-06/wmmLOEWK-pythia8-powheg2.root", "loewk", "POWHEG-2 LO + PYTHIA 8"));
  conf_list.push_back(Config("baseComp-05-06/wmmNLOEWK-pythia8-powheg2.root", "nloewk", "POWHEG-2 NLO + PYTHIA 8"));
  
  Config tune4C("baseComp-05-06/wmm-tune4C-powheg1-bacon.root","tune4C","POWHEG-1 + PYTHIA 8 + 4C tune");
  Config monash("baseComp-05-06/wmm-tuneMonash-bacon.root","monash","POWHEG-1 + PYTHIA 8 + Monash tune");

  Config amc("baseComp-05-06/WJets-bacon.root","amc","aMC@NLO + PYTHIA 8");
   
  for (UInt_t i=0; i<conf_list.size(); i++) {
    calcAcc(wmm, conf_list[i]);
  }
  calcAcc(wmm, amc);

  /*
  Int_t nP=50, lP=0, hP=100;
  Int_t nE=20, lE=-10, hE=10;

  for (UInt_t i=0; i<conf_list.size()-1; i++) {
    drawTwoPt(wmm, "genVf_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(wmm, "genVf_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
    drawTwoPt(wmm, "genL2f_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(wmm, "genL2f_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
  }

  drawThreePt(wmm, "genVf_pt", conf_list[2], tune4C, monash, nP, lP, hP);
  drawThreeEta(wmm, "genVf_eta", conf_list[2], tune4C, monash, nE, lE, hE);
  drawThreePt(wmm, "genL2f_pt", conf_list[2], tune4C, monash, nP, lP, hP);
  drawThreeEta(wmm, "genL2f_eta", conf_list[2], tune4C, monash, nE, lE, hE);

  drawThreePt(wmm, "genVf_pt", conf_list[0], conf_list[4], amc, nP, lP, hP);
  drawThreeEta(wmm, "genVf_eta", conf_list[0], conf_list[4], amc, nE, lE, hE);
  drawThreePt(wmm, "genL2f_pt", conf_list[0], conf_list[4], amc, nP, lP, hP);
  drawThreeEta(wmm, "genL2f_eta", conf_list[0], conf_list[4], amc, nE, lE, hE);
  */
}

void compareWpm() {

  vector<Config> conf_list;

  conf_list.push_back(Config("baseComp-05-06/CT10-wpm-pythia6-powheg1-v2-bacon.root", "ct10", "CT10 + POWHEG-1 + PYTHIA 6"));
  conf_list.push_back(Config("baseComp-05-06/wpm-pythia6-powheg1.root", "pow1py6", "POWHEG-1 + PYTHIA 6"));
  conf_list.push_back(Config("baseComp-05-06/wpm-pythia8-powheg1.root", "pow1py8", "POWHEG-1 + PYTHIA 8"));
  conf_list.push_back(Config("baseComp-05-06/wpmLOEWK-pythia8-powheg2.root", "loewk", "POWHEG-2 LO + PYTHIA 8"));
  conf_list.push_back(Config("baseComp-05-06/wpmNLOEWK-pythia8-powheg2.root", "nloewk", "POWHEG-2 NLO + PYTHIA 8"));
  
  //Config tune4C("baseComp-05-06/wpm-tune4C-powheg1-bacon.root","tune4C","POWHEG-1 + PYTHIA 8 + 4C tune");
  Config monash("baseComp-05-06/wpm-tuneMonash-bacon.root","monash","POWHEG-1 + PYTHIA 8 + Monash tune");  

  Config amc("baseComp-05-06/WJets-bacon.root","amc","aMC@NLO + PYTHIA 8");

  for (UInt_t i=0; i<conf_list.size(); i++) {
    calcAcc(wpm, conf_list[i]);
  }
  calcAcc(wpm, amc);
  /* 
  Int_t nP=50, lP=0, hP=100;
  Int_t nE=20, lE=-10, hE=10;

  for (UInt_t i=0; i<conf_list.size()-1; i++) {
    drawTwoPt(wpm, "genVf_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(wpm, "genVf_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
    drawTwoPt(wpm, "genL1f_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(wpm, "genL1f_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
  }

  drawTwoPt(wpm, "genVf_pt", conf_list[2], monash, nP, lP, hP);
  drawTwoEta(wpm, "genVf_eta", conf_list[2], monash, nE, lE, hE);
  drawTwoPt(wpm, "genL1f_pt", conf_list[2], monash, nP, lP, hP);
  drawTwoEta(wpm, "genL1f_eta", conf_list[2], monash, nE, lE, hE);

  drawThreePt(wpm, "genVf_pt", conf_list[0], conf_list[4], amc, nP, lP, hP);
  drawThreeEta(wpm, "genVf_eta", conf_list[0], conf_list[4], amc, nE, lE, hE);
  drawThreePt(wpm, "genL1f_pt", conf_list[0], conf_list[4], amc, nP, lP, hP);
  drawThreeEta(wpm, "genL1f_eta", conf_list[0], conf_list[4], amc, nE, lE, hE);
  */
}

void compareZmm() {

  vector<Config> conf_list;

  conf_list.push_back(Config("baseComp-05-06/CT10-zmm-pythia6-powheg1-v2-bacon.root", "ct10", "CT10 + POWHEG-1 + PYTHIA 6"));
  conf_list.push_back(Config("baseComp-05-06/zmm-pythia6-powheg1.root", "pow1py6", "POWHEG-1 + PYTHIA 6"));
  conf_list.push_back(Config("baseComp-05-06/zmm-pythia8-powheg1.root", "pow1py8", "POWHEG-1 + PYTHIA 8"));
  conf_list.push_back(Config("baseComp-05-06/zeeLOEWK-pythia8-powheg2.root", "loewk", "POWHEG-2 LO + PYTHIA 8"));
  conf_list.push_back(Config("baseComp-05-06/zeeNLOEWK-pythia8-powheg2.root", "nloewk", "POWHEG-2 NLO + PYTHIA 8"));
  
  Config tune4C("baseComp-05-06/zmm-tune4C-bacon.root","tune4C","POWHEG-1 + PYTHIA 8 + 4C tune");
  Config monash("baseComp-05-06/zmm-tuneMonash-bacon.root","monash","POWHEG-1 + PYTHIA 8 + Monash tune");

  Config amc("baseComp-05-06/DYJets-bacon.root","amc","aMC@NLO + PYTHIA 8");
  Config ct10("baseComp-05-06/CT10-zee-pythia6-powheg1-v2-bacon.root","ct10ee","CT10 + POWHEG-1 + PYTHIA 6");

  for (UInt_t i=0; i<conf_list.size(); i++) {
    if (i>2) calcAcc(zee, conf_list[i]);
    else calcAcc(zmm, conf_list[i]);
  }

  calcAcc(zmm, amc);
  calcAcc(zee, amc);
  calcAcc(zee, ct10);

  /*
  Int_t nP=50, lP=0, hP=100;
  Int_t nE=20, lE=-10, hE=10;

  for (UInt_t i=0; i<conf_list.size()-3; i++) {
    drawTwoPt(zmm, "genVf_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(zmm, "genVf_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
    drawTwoPt(zmm, "genL1f_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(zmm, "genL1f_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
    drawTwoPt(zmm, "genL2f_pt", conf_list[i], conf_list[i+1], nP, lP, hP);
    drawTwoEta(zmm, "genL2f_eta", conf_list[i], conf_list[i+1], nE, lE, hE);
  }

  drawTwoPt(zee, "genVf_pt", conf_list[3], conf_list[4], nP, lP, hP);
  drawTwoEta(zee, "genVf_eta", conf_list[3], conf_list[4], nE, lE, hE);
  drawTwoPt(zee, "genL1f_pt", conf_list[3], conf_list[4], nP, lP, hP);
  drawTwoEta(zee, "genL1f_eta", conf_list[3], conf_list[4], nE, lE, hE);
  drawTwoPt(zee, "genL2f_pt", conf_list[3], conf_list[4], nP, lP, hP);
  drawTwoEta(zee, "genL2f_eta", conf_list[3], conf_list[4], nE, lE, hE);

  drawThreePt(zmm, "genVf_pt", conf_list[2], tune4C, monash, nP, lP, hP);
  drawThreeEta(zmm, "genVf_eta", conf_list[2], tune4C, monash, nE, lE, hE);
  drawThreePt(zmm, "genL1f_pt", conf_list[2], tune4C, monash, nP, lP, hP);
  drawThreeEta(zmm, "genL1f_eta", conf_list[2], tune4C, monash, nE, lE, hE);
  drawThreePt(zmm, "genL2f_pt", conf_list[2], tune4C, monash, nP, lP, hP);
  drawThreeEta(zmm, "genL2f_eta", conf_list[2], tune4C, monash, nE, lE, hE);

  drawThreePt(zee, "genVf_pt", ct10, conf_list[4], amc, nP, lP, hP);
  drawThreeEta(zee, "genVf_eta", ct10, conf_list[4], amc, nE, lE, hE);
  drawThreePt(zee, "genL1f_pt", ct10, conf_list[4], amc, nP, lP, hP);
  drawThreeEta(zee, "genL1f_eta", ct10, conf_list[4], amc, nE, lE, hE);
  drawThreePt(zee, "genL2f_pt", ct10, conf_list[4], amc, nP, lP, hP);
  drawThreeEta(zee, "genL2f_eta", ct10, conf_list[4], amc, nE, lE, hE);
  */
}
