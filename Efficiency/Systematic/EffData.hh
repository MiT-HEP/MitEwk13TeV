#ifndef EFF_DATA_HH
#define EFF_DATA_HH

struct EffData
{
  Float_t mass, pt, eta, phi, weight;
  Int_t q;
  UInt_t npv, npu, pass;
  UInt_t runNum, lumiSec, evtNum;
};

// "mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum"

#endif
