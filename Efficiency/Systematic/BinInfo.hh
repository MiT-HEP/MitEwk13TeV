#ifndef BIN_INFO_HH
#define BIN_INFO_HH

struct BinInfo
{

  UInt_t nEvents, nBkgFail, nBkgPass, iBin, absEta;
  Float_t ptLo, ptHi;
  Float_t etaLo, etaHi;
  Float_t phiLo, phiHi;
  Float_t npvLo, npvHi;

};

//"nEvents/i:nBkgFail:nBkgPassibin:absEta:ptLo/F:ptHi:etaLo:etaHi:phiLo:phiHi:npvLo:npvHi"

#endif
