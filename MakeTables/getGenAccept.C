void calcAccW(TString fname);
void calcAccZ(TString fname);

void getGenAccept() {

  calcAccW("/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc/wme_NNPDF30_nlo_as_0118.root");
  calcAccW("/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc/wpe_NNPDF30_nlo_as_0118.root");
  calcAccW("/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc/wmm_NNPDF30_nlo_as_0118.root");
  calcAccW("/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc/wpm_NNPDF30_nlo_as_0118.root");
  calcAccZ("/afs/cern.ch/work/j/jlawhorn/public/wz-pdf/NNPDF30_amc/zee_NNPDF30_nlo_as_0118.root");

}

void calcAccW(TString fname) {

  TFile *f = new TFile(fname);

  TH1D *tot   = (TH1D*) f->Get("dTot_NNPDF30_nlo_as_0118_0");
  TH1D *preB  = (TH1D*) f->Get("dPreB_NNPDF30_nlo_as_0118_0");
  TH1D *preE  = (TH1D*) f->Get("dPreE_NNPDF30_nlo_as_0118_0");
  TH1D *postB = (TH1D*) f->Get("dPostB_NNPDF30_nlo_as_0118_0");
  TH1D *postE = (TH1D*) f->Get("dPostE_NNPDF30_nlo_as_0118_0");

  cout << setprecision(4) << (preB->Integral()+preE->Integral())/(tot->Integral()) << "\\\% & ";
  cout << setprecision(4) << (postB->Integral()+postE->Integral())/(tot->Integral()) << "\\\% & " << endl;

}

void calcAccZ(TString fname) {

  TFile *f = new TFile(fname);

  TH1D *tot    = (TH1D*) f->Get("dTot_NNPDF30_nlo_as_0118_0");
  TH1D *preBB  = (TH1D*) f->Get("dPreBB_NNPDF30_nlo_as_0118_0");
  TH1D *preBE  = (TH1D*) f->Get("dPreBB_NNPDF30_nlo_as_0118_0");
  TH1D *preEE  = (TH1D*) f->Get("dPreEE_NNPDF30_nlo_as_0118_0");
  TH1D *postBB = (TH1D*) f->Get("dPostBB_NNPDF30_nlo_as_0118_0");
  TH1D *postBE = (TH1D*) f->Get("dPostBB_NNPDF30_nlo_as_0118_0");
  TH1D *postEE = (TH1D*) f->Get("dPostEE_NNPDF30_nlo_as_0118_0");

  cout << setprecision(4) << (preBB->Integral()+preBE->Integral()+preEE->Integral())/(tot->Integral()) << "\\\% & ";
  cout << setprecision(4) << (postBB->Integral()+postBE->Integral()+postEE->Integral())/(tot->Integral()) << "\\\% & " << endl;

}
