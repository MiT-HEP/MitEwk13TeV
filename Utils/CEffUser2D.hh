#ifndef CEFFUSER2D_HH
#define CEFFUSER2D_HH

#include <TH2D.h>
#include <iostream>

class CEffUser2D
{
public:
  CEffUser2D();
  ~CEffUser2D();
  
  void    loadEff(TH2D* eff, TH2D* errl, TH2D* errh);
  Float_t getEff(const Double_t x, const Double_t y);
  Float_t getErrLow(const Double_t x, const Double_t y);
  Float_t getErrHigh(const Double_t x, const Double_t y);    
  void    printEff(std::ostream& os);
  void    printErrLow(std::ostream& os);
  void    printErrHigh(std::ostream& os);

protected:
  void    printHist2D(const TH2D* h,std::ostream& os);
  Float_t getValue(const TH2D* h, const Double_t x, const Double_t y);

  TH2D *hEff;
  TH2D *hErrl;
  TH2D *hErrh;
};

#endif
