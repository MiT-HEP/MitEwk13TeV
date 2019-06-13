#include <stdio.h>
#include <iostream>
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

const int NBpt  = 2;
const float ptrange[NBpt+1]   = {25., 40., 8000.};
const int NBeta = 12;
const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};


// const int NBeta = 12;
// const float etarange[NBeta+1] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4};
// const int NBpt = 8;
// const float ptrange[NBpt+1] = {25,30,35,40,45,50,60,80,8000};

// TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/testReweights_v2_2/results/";
// TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV_v0/results/TOYS/";
TString masterDir="/afs/cern.ch/user/s/sabrandt/lowPU/CMSSW_9_4_12/src/MitEwk13TeV/Efficiency/LowPU2017ID_13TeV_v1/results/TOYS/";
TString sigDirFSR="_POWxPythia_v2_POWxPhotos_v2/";
TString sigDirMC="_aMCxPythia_v2_minloxPythia_v2/";
// TString bkgDir="_aMCxPythia_v1_POWBKG_v1/"; // should be exp vs Powerlaw
TString bkgDir="_POWBKG_v2_aMCxPythia_v2/"; // should be exp vs Powerlaw

TString subfolder="";

// TString sigDirMC="_aMCxPythia_v0_POWBKG_v0/";
// TString sigDirFSR="_aMCxPythia_v0_POWBKG_v0/";
// TString bkgDir="_aMCxPythia_v0_POWBKG_v0/"; // should be exp vs Powerlaw

// For the closure tests
// TString sigDirMC="_aMCxPythia_v0_aMCxPythia_v0/";
// TString sigDirFSR="_aMCxPythia_v0_aMCxPythia_v0/";
// TString bkgDir="_aMCxPythia_v0_aMCxPythia_v0/"; 

// const string charges[2] = {"Negative","Positive"};
const string charges[2] = {"Combined","Combined"};

/* old
const float StaSigSys[NBpt*NBeta] = {1.00, 1.28, 0.52, 0.69, 1.48, 1.11, 1.00, 0.97,
				     0.35, 0.22, 0.27, 0.11, 0.11, 0.19, 0.23, 0.34};
const float StaBkgSys[NBpt*NBeta] = {0.10, 0.33, 0.39, 1.30, 0.05, 0.14, 0.10, 0.22,
				     0.13, 0.05, 0.06, 0.02, 0.03, 0.02, 0.01, 0.26};
const float SITSigSys[NBpt*NBeta] = {0.03, 0.16, 1.52, 0.02, 1.10, 0.77, 0.58, 0.01,
				     0.02, 0.02, 0.81, 0.12, 0.01, 0.45, 0.05, 0.04};
const float SITBkgSys[NBpt*NBeta] = {0.13, 0.01, 0.07, 0.00, 0.01, 0.01, 0.15, 0.00,
				     0.11, 0.00, 0.07, 0.01, 0.00, 0.00, 0.00, 0.15};
*/
const float StaSigSys[NBpt*NBeta] = {0.50, 0.93, 0.74, 0.64, 0.67, 0.73, 0.55, 0.41,
                                     0.15, 0.26, 0.17, 0.11, 0.08, 0.21, 0.17, 0.13};
const float StaBkgSys[NBpt*NBeta] = {0.25, 0.76, 0.31, 0.22, 0.30, 0.16, 0.56, 0.46,
                                     0.02, 0.03, 0.01, 0.01, 0.01, 0.01, 0.03, 0.01};
const float SITSigSys[NBpt*NBeta] = {0.38, 0.57, 1.25, 0.78, 0.79, 1.00, 0.63, 0.51,
                                     0.46, 0.28, 0.32, 0.29, 0.35, 0.61, 0.21, 0.39};
const float SITBkgSys[NBpt*NBeta] = {0.03, 0.08, 0.14, 0.11, 0.11, 0.12, 0.07, 0.06,
                                     0.01, 0.00, 0.02, 0.00, 0.00, 0.00, 0.00, 0.01};
void makeTH2DMu(TString effTypeSig = "MuStaEff"){
// void makeTH2DMu(TString effTypeSig = "MuSITEff"){
    // alternate input option is MuStaEff, HLT should not have fits?MuSITEff
    bool doAbs=false; 
    vector<double> vSigMCNeg;
    vector<double> vSigMCPos;
    vector<double> vSigFSRNeg;
    vector<double> vSigFSRPos;
    vector<double> vBkgNeg; 
    vector<double> vBkgPos;
    
  for(int i = 0; i < NBeta*NBpt; ++i){
    double value=0;
    char infilename[200];
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirFSR.Data(),subfolder.Data(),charges[0].c_str(),i);
    ifstream infile1(infilename); assert (infile1); infile1>>value; doAbs?value=fabs(value):value=value;vSigFSRNeg.push_back(value); value=0;
    std::cout << infilename << vSigFSRNeg.back() << std::endl;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirFSR.Data(),subfolder.Data(),charges[1].c_str(),i);
    ifstream infile2(infilename); assert (infile2);infile2>>value; doAbs?value=fabs(value):value=value;vSigFSRPos.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirMC.Data(),subfolder.Data(),charges[0].c_str(),i);
    ifstream infile3(infilename); assert (infile3);infile3>>value; doAbs?value=fabs(value):value=value;vSigMCNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),sigDirMC.Data(),subfolder.Data(),charges[1].c_str(),i);
    ifstream infile4(infilename); assert (infile4);infile4>>value;doAbs?value=fabs(value):value=value;vSigMCPos.push_back(value); value=0;    
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),bkgDir.Data(),subfolder.Data(),charges[0].c_str(),i);
    ifstream infile5(infilename); assert (infile5);infile5>>value; doAbs?value=fabs(value):value=value;vBkgNeg.push_back(value); value=0;
    sprintf(infilename,"%s%s%s%s%s/Sig_pull_%d.txt",masterDir.Data(),effTypeSig.Data(),bkgDir.Data(),subfolder.Data(),charges[1].c_str(),i);
    ifstream infile6(infilename); infile6>>value; doAbs?value=fabs(value):value=value;vBkgPos.push_back(value); value=0; 
  }
    
    char histname[100];
    sprintf(histname,"h%sSigMCNeg",effTypeSig.Data());
	TH2D *hSigMCNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sSigMCPos",effTypeSig.Data());
	TH2D *hSigMCPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sSigFSRNeg",effTypeSig.Data());
    TH2D *hSigFSRNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sSigFSRPos",effTypeSig.Data());
    TH2D *hSigFSRPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sBkgNeg",effTypeSig.Data());
    TH2D *hBkgNeg = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname,"h%sBkgPos",effTypeSig.Data());
    TH2D *hBkgPos = new TH2D(histname,histname, NBeta, etarange, NBpt, ptrange);

	for(int ipt=0; ipt<NBpt; ipt++){

		for(int ieta=0; ieta<NBeta; ieta++){
			hSigFSRNeg->SetBinContent(ieta+1, ipt+1, 1.+vSigFSRNeg[ipt*NBeta+ieta]);
			hSigFSRPos->SetBinContent(ieta+1, ipt+1, 1.+vSigFSRPos[ipt*NBeta+ieta]);
			hSigMCNeg->SetBinContent(ieta+1, ipt+1, 1.+vSigMCNeg[ipt*NBeta+ieta]);
			hSigMCPos->SetBinContent(ieta+1, ipt+1, 1.+vSigMCPos[ipt*NBeta+ieta]);
			hBkgNeg->SetBinContent(ieta+1, ipt+1, 1.+vBkgNeg[ipt*NBeta+ieta]);
			hBkgPos->SetBinContent(ieta+1, ipt+1, 1.+vBkgPos[ipt*NBeta+ieta]);
		}
	}
    char outfilename[100];
    sprintf(outfilename,"SysUnc_%s.root",effTypeSig.Data());
	TFile *f = new TFile(outfilename, "RECREATE");
	hSigFSRNeg->Write();
	hSigFSRPos->Write();
	hSigMCNeg->Write();
	hSigMCPos->Write();
	hBkgNeg->Write();
	hBkgPos->Write();
	f->Close();
}
