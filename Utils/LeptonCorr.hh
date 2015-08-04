#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getEleScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 1.0/1.004; }
    else if (fabs(eta) < 0.8)    { return 1.0/1.002; }
    else if (fabs(eta) < 1.4442) { return 1.0/1.001; }
    else                         { return 1.0/1.014; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(1.004+0.005 ); }
    else if (fabs(eta) < 0.8)    { return 1.0/(1.002+0.006 ); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.001+0.001 ); }
    else                         { return 1.0/(1.014+0.017 ); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(1.004-0.005 ); }
    else if (fabs(eta) < 0.8)    { return 1.0/(1.002-0.006 ); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.001-0.001 ); }
    else                         { return 1.0/(1.014-0.017 ); }
  }
  else return -1.0;
}

Double_t getEleResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 0.43; }
    else if (fabs(eta) < 0.8)    { return 0.61; }
    else if (fabs(eta) < 1.4442) { return 1.00; }
    else                         { return 1.40; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 0.43+0.07; }
    else if (fabs(eta) < 0.8)    { return 0.61+0.06; }
    else if (fabs(eta) < 1.4442) { return 1.00+0.05; }
    else                         { return 1.40+0.07; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 0.43-0.07; }
    else if (fabs(eta) < 0.8)    { return 0.61-0.06; }
    else if (fabs(eta) < 1.4442) { return 1.00-0.05; }
    else                         { return 1.40-0.07; }
  }
  else return -1.0;
}

Double_t getMuScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 1.2) { return 1.0/1.0014 ; }
    else if (fabs(eta) < 2.1) { return 1.0/1.00002; }
    else                      { return 1.0/1.002 ; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 1.2) { return 1.0/1.0014 +0.0003; }
    else if (fabs(eta) < 2.1) { return 1.0/1.00002+0.00007; }
    else                      { return 1.0/1.002  +0.002; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 1.2) { return 1.0/1.0014 -0.0003; }
    else if (fabs(eta) < 2.1) { return 1.0/1.00002-0.00007; }
    else                      { return 1.0/1.002  -0.002; }
  }
  else return -1.0;
}
Double_t getMuResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 1.2) { return 0.43; }
    else if (fabs(eta) < 2.1) { return 0.65; }
    else                      { return 2.0; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 1.2) { return 0.43+0.05; }
    else if (fabs(eta) < 2.1) { return 0.65+0.07; }
    else                      { return 2.0 +0.2 ; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 1.2) { return 0.43-0.05; }
    else if (fabs(eta) < 2.1) { return 0.65-0.07; }
    else                      { return 2.0 -0.2 ; }
  }
  else return -1.0;
}

#endif
