#ifndef EWKANA_UTILS_LEPTONCORR_HH
#define EWKANA_UTILS_LEPTONCORR_HH

Double_t getEleScaleCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 1.0/0.998062; }
    else if (fabs(eta) < 0.8)    { return 1.0/0.998656; }
    else if (fabs(eta) < 1.4442) { return 1.0/1.00281; }
    else                         { return 1.0/1.00496; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(0.998062+8.98619e-05); }
    else if (fabs(eta) < 0.8)    { return 1.0/(0.998656+8.62067e-05); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.00281 +7.5253e-05); }
    else                         { return 1.0/(1.00496 +0.000114103); }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 1.0/(0.998062-8.98619e-05); }
    else if (fabs(eta) < 0.8)    { return 1.0/(0.998656-8.62067e-05); }
    else if (fabs(eta) < 1.4442) { return 1.0/(1.00281 -7.5253e-05); }
    else                         { return 1.0/(1.00496 -0.000114103); }
  }
  else return -1.0;
}

Double_t getEleResCorr(const Double_t eta, const Int_t sigma)
{
  if (sigma==0) {
    if      (fabs(eta) < 0.4)    { return 0.361086; }
    else if (fabs(eta) < 0.8)    { return 0.371181; }
    else if (fabs(eta) < 1.4442) { return 0.658034; }
    else                         { return 1.09921; }
  }
  else if (sigma==1) {
    if      (fabs(eta) < 0.4)    { return 0.361086+0.0196288; }
    else if (fabs(eta) < 0.8)    { return 0.371181+0.020926; }
    else if (fabs(eta) < 1.4442) { return 0.658034+0.0149622; }
    else                         { return 1.09921 +0.0121135; }
  }
  else if (sigma==-1) {
    if      (fabs(eta) < 0.4)    { return 0.361086-0.0196288; }
    else if (fabs(eta) < 0.8)    { return 0.371181-0.020926; }
    else if (fabs(eta) < 1.4442) { return 0.658034-0.0149622; }
    else                         { return 1.09921 -0.0121135; }
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
