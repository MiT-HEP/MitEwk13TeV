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

const int NBeta = 6;
const float etarange[NBeta+1] = {-2.5, -1.5, -0.8, 0, 0.8, 1.5, 2.5};
const float sysdt[NBeta] = {0.020, 0.0017, 0.00055, 0.00055, 0.0017, 0.020}; 
const float sysmc[NBeta] = {0.018, 0.0022, 0.00058, 0.00058, 0.0022, 0.018};

void makeTH1DEle(){

	TH1D *hdt = new TH1D("hdt","", NBeta, etarange);
	TH1D *hmc = new TH1D("hmc","", NBeta, etarange);


	for(int ieta=0; ieta<NBeta; ieta++){
		hdt->SetBinContent(ieta+1, sysdt[ieta]);
		hmc->SetBinContent(ieta+1, sysmc[ieta]);
	}

	TFile *f = new TFile("EleChargeMisIDSys.root", "RECREATE");
	hdt->Write();
	hmc->Write();
	f->Close();
}
