#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getEleScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 1.0/0.998094; }
    else if (fabs(eta) < 0.8)    { return 1.0/0.998722; }
    else if (fabs(eta) < 1.4442) { return 1.0/1.00281; }
    else                         { return 1.0/1.00499; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(0.998094+0.000270935; }
    else if (fabs(eta) < 0.8)    { return 1.0/(0.998722+0.000284622; }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.00281 +4.52316e-05; }
    else                         { return 1.0/(1.00499 +0.000110672; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(0.998094-0.000270935; }
    else if (fabs(eta) < 0.8)    { return 1.0/(0.998722-0.000284622; }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.00281 -4.52316e-05; }
    else                         { return 1.0/(1.00499 -0.000110672; }
  }
  else return -1.0;
}

Double_t getEleResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 0.342735; }
    else if (fabs(eta) < 0.8)    { return 0.37259; }
    else if (fabs(eta) < 1.4442) { return 0.647278; }
    else                         { return 1.08814; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 0.342735+0.0220807; }
    else if (fabs(eta) < 0.8)    { return 0.37259 +0.0224161; }
    else if (fabs(eta) < 1.4442) { return 0.647278+0.0151588; }
    else                         { return 1.08814 +0.0122041; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 0.342735-0.0220807; }
    else if (fabs(eta) < 0.8)    { return 0.37259 -0.0224161; }
    else if (fabs(eta) < 1.4442) { return 0.647278-0.0151588; }
    else                         { return 1.08814 -0.0122041; }
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
