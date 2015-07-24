#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.5)    { return 0.998; }
    else if (fabs(eta) < 1.0)    { return 0.998; }
    else if (fabs(eta) < 1.4442) { return 1.003; }
    else if (fabs(eta) < 2.0)    { return 0.996; }
    else                         { return 1.017; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.5)    { return 0.998+0.001; }
    else if (fabs(eta) < 1.0)    { return 0.998+0.002; }
    else if (fabs(eta) < 1.4442) { return 1.003+0.001; }
    else if (fabs(eta) < 2.0)    { return 0.996+0.001; }
    else                         { return 1.017+0.001; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.5)    { return 0.998-0.001; }
    else if (fabs(eta) < 1.0)    { return 0.998-0.002; }
    else if (fabs(eta) < 1.4442) { return 1.003-0.001; }
    else if (fabs(eta) < 2.0)    { return 0.996-0.001; }
    else                         { return 1.017-0.001; }
  }
  else return -1.0;
}

Double_t getResCorr(const Double_t eta, const Int_t sigma)
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

#endif
