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
    if      (fabs(eta) < 1.2) { return 1.0/1.00146 ; }
    else if (fabs(eta) < 2.1) { return 1.0/1.00039; }
    else                      { return 1.0/1.00118 ; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 1.2) { return 1.0/1.00146 +0.00040583; }
    else if (fabs(eta) < 2.1) { return 1.0/1.00039+0.000582759; }
    else                      { return 1.0/1.00118  +0.00166661; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 1.2) { return 1.0/1.00146 -0.00040583; }
    else if (fabs(eta) < 2.1) { return 1.0/1.00039-0.000582759; }
    else                      { return 1.0/1.00118  -0.00166661; }
  }
  else return -1.0;
}
Double_t getMuResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 1.2) { return 0.379459; }
    else if (fabs(eta) < 2.1) { return 0.637196; }
    else                      { return 1.92794; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 1.2) { return 0.379459+0.0547912; }
    else if (fabs(eta) < 2.1) { return 0.637196+0.0703625; }
    else                      { return 1.92794 +0.119745 ; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 1.2) { return 0.379459-0.0547912; }
    else if (fabs(eta) < 2.1) { return 0.637196-0.0703625; }
    else                      { return 1.92794 -0.119745 ; }
  }
  else return -1.0;
}

#endif
