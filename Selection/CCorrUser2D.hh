#ifndef CCORRUSER2D_HH
#define CCORRUSER2D_HH

#include <TH2D.h>
#include <iostream>

class CCorrUser2D
{
public:
  CCorrUser2D();
  ~CCorrUser2D();
  
  // void    loadCorr(TH2D* corr, TH2D* errl, TH2D* errh);
  void    loadCorr(TH2D* corr);
  Float_t getCorr(const Double_t x, const Double_t y);
  Float_t getErrLow(const Double_t x, const Double_t y);
  Float_t getErrHigh(const Double_t x, const Double_t y);    
  void    printCorr(std::ostream& os);
  void    printErrLow(std::ostream& os);
  void    printErrHigh(std::ostream& os);

protected:
  void    printHist2D(const TH2D* h,std::ostream& os);
  Float_t getValue(const TH2D* h, const Double_t x, const Double_t y);

  TH2D *hCorr;
  TH2D *hErrl;
  TH2D *hErrh;
};

#endif
