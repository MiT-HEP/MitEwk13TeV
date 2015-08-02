#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getEleScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 1.0/1.003; }
    else if (fabs(eta) < 0.8)    { return 1.0/1.002; }
    else if (fabs(eta) < 1.4442) { return 1.0/1.001; }
    else                         { return 1.0/1.012; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(1.003+0.005 ); }
    else if (fabs(eta) < 0.8)    { return 1.0/(1.002+0.009 ); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.001+0.001 ); }
    else                         { return 1.0/(1.012+0.004 ); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(1.003-0.005 ); }
    else if (fabs(eta) < 0.8)    { return 1.0/(1.002-0.009 ); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.001-0.001 ); }
    else                         { return 1.0/(1.012-0.004 ); }
  }
  else return -1.0;
}

Double_t getEleResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 0.58; }
    else if (fabs(eta) < 0.8)    { return 0.6 ; }
    else if (fabs(eta) < 1.4442) { return 1.10; }
    else                         { return 1.37; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 0.58+0.08; }
    else if (fabs(eta) < 0.8)    { return 0.6 +0.1 ; }
    else if (fabs(eta) < 1.4442) { return 1.10+0.08; }
    else                         { return 1.37+0.08; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 0.58-0.08; }
    else if (fabs(eta) < 0.8)    { return 0.6 -0.1 ; }
    else if (fabs(eta) < 1.4442) { return 1.10-0.08; }
    else                         { return 1.37-0.08; }
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
