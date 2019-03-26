#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

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
const int NBeta = 12;
const float ptrange[NBpt+1]   = {25., 40., 8000.};
const float etarange[NBeta+1] = {-2.4, -2.1, -1.2, -0.9, -0.3, -0.2, 0., 0.2, 0.3, 0.9, 1.2, 2.1, 2.4};
/* 50ns
const float StaSigSys[NBpt*NBeta] = {1.00, 1.28, 0.52, 0.69, 1.48, 1.11, 1.00, 0.97,
				     0.35, 0.22, 0.27, 0.11, 0.11, 0.19, 0.23, 0.34};
const float StaBkgSys[NBpt*NBeta] = {0.10, 0.33, 0.39, 1.30, 0.05, 0.14, 0.10, 0.22,
				     0.13, 0.05, 0.06, 0.02, 0.03, 0.02, 0.01, 0.26};
const float SITSigSys[NBpt*NBeta] = {0.03, 0.16, 1.52, 0.02, 1.10, 0.77, 0.58, 0.01,
				     0.02, 0.02, 0.81, 0.12, 0.01, 0.45, 0.05, 0.04};
const float SITBkgSys[NBpt*NBeta] = {0.13, 0.01, 0.07, 0.00, 0.01, 0.01, 0.15, 0.00,
				     0.11, 0.00, 0.07, 0.01, 0.00, 0.00, 0.00, 0.15};
*/
/* 25ns 74X 
const float StaSigSys[NBpt*NBeta] = {0.51, 0.77, 0.66, 0.25, 0.05, 0.40, 0.25, 1.44,
				     0.94, 0.83, 0.28, 0.43, 0.21, 0.38, 0.04, 0.16,
				     0.33, 0.07, 0.11, 0.29, 0.04, 0.21, 0.02, 0.13}; 
const float StaBkgSys[NBpt*NBeta] = {0.78, 0.31, 0.40, 0.27, 0.95, 0.03, 0.11, 1.05,
				     0.23, 0.51, 0.44, 0.52, 0.01, 0.22, 0.01, 0.34,
				     0.33, 0.12, 0.00, 0.17, 0.00, 0.12, 0.11, 0.01};
const float SITSigSys[NBpt*NBeta] = {0.36, 0.57, 1.25, 0.92, 0.65, 0.70, 0.79, 0.69,
				     0.89, 0.99, 0.52, 0.50, 0.45, 0.27, 0.32, 0.26,
				     0.38, 0.33, 0.35, 0.29, 0.34, 0.49, 0.20, 0.39};
const float SITBkgSys[NBpt*NBeta] = {0.03, 0.07, 0.14, 0.12, 0.08, 0.08, 0.09, 0.06,
				     0.12, 0.13, 0.06, 0.05, 0.01, 0.00, 0.00, 0.00,
				     0.02, 0.00, 0.00, 0.02, 0.00, 0.00, 0.00, 0.00};
*/
const float StaSigSys[NBpt*NBeta] = {0.17, 0.26, 0.62, 0.19, 0.73, 0.15, 0.03, 1.57,
                                     0.24, 0.28, 0.13, 0.21, 0.19, 0.04, 0.09, 0.04,
                                     0.25, 0.03, 0.03, 0.46, 0.05, 0.03, 0.06, 0.05}; 
const float StaBkgSys[NBpt*NBeta] = {0.10, 0.75, 0.65, 0.63, 0.91, 0.27, 0.31, 0.96,
                                     0.70, 0.63, 0.33, 0.40, 0.14, 0.11, 0.18, 0.17,
                                     0.50, 0.10, 0.39, 0.30, 0.13, 0.13, 0.22, 0.09};
const float SITSigSys[NBpt*NBeta] = {0.00, 0.43, 0.91, 0.71, 0.33, 0.40, 0.50, 0.45,
                                     0.63, 0.68, 0.37, 0.21, 0.46, 0.21, 0.27, 0.03,
                                     0.33, 0.26, 0.27, 0.26, 0.28, 0.38, 0.18, 0.35};
const float SITBkgSys[NBpt*NBeta] = {0.00, 0.05, 0.11, 0.09, 0.03, 0.04, 0.06, 0.03,
                                     0.08, 0.09, 0.03, 0.00, 0.00, 0.00, 0.00, 0.00,
                                     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
void makeTH2DMu_2016(){

	TH2D *h = new TH2D("h","", NBeta, etarange, NBpt, ptrange);

	for(int ieta=0; ieta<NBeta; ieta++){

		for(int ipt=0; ipt<NBpt; ipt++){
			h->SetBinContent(ieta+1, ipt+1, 1.+SITSigSys[ipt*NBeta+ieta]*0.01);
		}
	}

	TFile *f = new TFile("MuSITSigSys.root", "RECREATE");
	h->Write();
	f->Close();
}
