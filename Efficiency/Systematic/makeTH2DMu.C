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
const int NBeta = 8;
const float ptrange[NBpt+1]   = {25., 40., 8000.};
const float etarange[NBeta+1] = {-2.4, -2.1, -1.2, -0.9, 0., 0.9, 1.2, 2.1, 2.4};


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
void makeTH2DMu(){

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
