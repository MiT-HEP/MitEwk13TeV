#ifndef CEFFUSER1D_HH
#define CEFFUSER1D_HH

#include <TGraphAsymmErrors.h>
#include <iostream>

class CEffUser1D
{
public:
  CEffUser1D();
  ~CEffUser1D();
  
  void    loadEff(TGraphAsymmErrors* gr);
  Float_t getEff(const Double_t x);
  Float_t getErrLow(const Double_t x);
  Float_t getErrHigh(const Double_t x);    
  void    printEff(std::ostream& os);
  void    printErrLow(std::ostream& os);
  void    printErrHigh(std::ostream& os);

protected:
  Int_t   getBin(const Double_t x);
  void    print(const Double_t *yval, std::ostream& os);  
  
  TGraphAsymmErrors *graph;
};

#endif
