#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getEleScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.5)    { return 1.0/0.998; }
    else if (fabs(eta) < 1.0)    { return 1.0/0.998; }
    else if (fabs(eta) < 1.4442) { return 1.0/1.003; }
    else if (fabs(eta) < 2.0)    { return 1.0/0.996; }
    else                         { return 1.0/1.017; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.5)    { return 1.0/(0.998+0.001); }
    else if (fabs(eta) < 1.0)    { return 1.0/(0.998+0.002); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.003+0.001); }
    else if (fabs(eta) < 2.0)    { return 1.0/(0.996+0.001); }
    else                         { return 1.0/(1.017+0.001); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.5)    { return 1.0/(0.998-0.001); }
    else if (fabs(eta) < 1.0)    { return 1.0/(0.998-0.002); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.003-0.001); }
    else if (fabs(eta) < 2.0)    { return 1.0/(0.996-0.001); }
    else                         { return 1.0/(1.017-0.001); }
  }
  else return -1.0;
}

Double_t getEleResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.5)    { return 0.8; }
    else if (fabs(eta) < 1.0)    { return 0.9; }
    else if (fabs(eta) < 1.4442) { return 1.7; }
    else if (fabs(eta) < 2.0)    { return 0.1; }
    else                         { return 0.2; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.5)    { return 0.8+0.2; }
    else if (fabs(eta) < 1.0)    { return 0.9+0.2; }
    else if (fabs(eta) < 1.4442) { return 1.7+0.2; }
    else if (fabs(eta) < 2.0)    { return 0.1+0.3; }
    else                         { return 0.2+0.1; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.5)    { return 0.8-0.2; }
    else if (fabs(eta) < 1.0)    { return 0.9-0.2; }
    else if (fabs(eta) < 1.4442) { return 1.7-0.2; }
    else if (fabs(eta) < 2.0)    { return 0.1-0.3; }
    else                         { return 0.2-0.1; }
  }
  else return -1.0;
}

Double_t getMuScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.8) { return 1.0/1.002; }
    else if (fabs(eta) < 1.6) { return 1.0/1.000; }
    else                      { return 1.0/1.002; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.8) { return 1.0/(1.002+0.001); }
    else if (fabs(eta) < 1.6) { return 1.0/(1.000+0.001); }
    else                      { return 1.0/(1.002+0.001); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.8) { return 1.0/(1.002-0.001); }
    else if (fabs(eta) < 1.6) { return 1.0/(1.000-0.001); }
    else                      { return 1.0/(1.002-0.001); }
  }
  else return -1.0;
}

/* Data->MC scale correction
$0 < |\eta| < 0.8$ & $1.00162$ \pm $0.00113235$ \\
$0.8 < |\eta| < 1.6$ & $1.00033$ \pm $0.000966038$ \\
$1.6 < |\eta| < 2.4$ & $1.00182$ \pm $0.000943297$ \\

  MC->Data resolution correction [GeV]
0 < |\eta| < 0.8 & $0.452368$ \pm $0.237725$ \\
0.8 < |\eta| < 1.6 & $0.529896$ \pm $0.222717$ \\
1.6 < |\eta| < 2.4 & $0.655471$ \pm $0.276792$ \\*/

Double_t getMuResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.8) { return 0.5; }
    else if (fabs(eta) < 1.6) { return 0.5; }
    else                      { return 0.7; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.8) { return 0.5+0.2; }
    else if (fabs(eta) < 1.6) { return 0.5+0.2; }
    else                      { return 0.7+0.3; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.8) { return 0.5-0.2; }
    else if (fabs(eta) < 1.6) { return 0.5-0.2; }
    else                      { return 0.7-0.3; }
  }
  else return -1.0;
}

#endif
