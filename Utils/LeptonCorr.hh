#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getEleScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 1.0/1.00363; }
    else if (fabs(eta) < 0.8)    { return 1.0/1.00181; }
    else if (fabs(eta) < 1.4442) { return 1.0/1.00128; }
    else                         { return 1.0/1.01365; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(1.00363+0.00043495); }
    else if (fabs(eta) < 0.8)    { return 1.0/(1.00181+0.000579655); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.00128+0.000923351); }
    else                         { return 1.0/(1.01365+0.0006267); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(1.00363-0.00043495); }
    else if (fabs(eta) < 0.8)    { return 1.0/(1.00181-0.000579655); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.00128-0.000923351); }
    else                         { return 1.0/(1.01365-0.0006267); }
  }
  else return -1.0;
}

Double_t getEleResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 0.430059; }
    else if (fabs(eta) < 0.8)    { return 0.614233 ; }
    else if (fabs(eta) < 1.4442) { return 1.00031; }
    else                         { return 1.39975; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 0.430059+0.128463; }
    else if (fabs(eta) < 0.8)    { return 0.614233+0.107999 ; }
    else if (fabs(eta) < 1.4442) { return 1.00031+0.0879275; }
    else                         { return 1.39975+0.0847074; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 0.430059-0.128463; }
    else if (fabs(eta) < 0.8)    { return 0.614233-0.107999 ; }
    else if (fabs(eta) < 1.4442) { return 1.00031-0.0879275; }
    else                         { return 1.39975-0.0847074; }
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
    if      (fabs(eta) < 1.2) { return 1.0/(1.00146+0.00040583); }
    else if (fabs(eta) < 2.1) { return 1.0/(1.00039+0.000582759); }
    else                      { return 1.0/(1.00118+0.00166661); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 1.2) { return 1.0/(1.00146-0.00040583); }
    else if (fabs(eta) < 2.1) { return 1.0/(1.00039-0.000582759); }
    else                      { return 1.0/(1.00118-0.00166661); }
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
    else                      { return 1.92794+0.119745; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 1.2) { return 0.379459-0.0547912; }
    else if (fabs(eta) < 2.1) { return 0.637196-0.0703625; }
    else                      { return 1.92794-0.119745; }
  }
  else return -1.0;
}

#endif
