#include "CCorrUser2D.hh"
#include <iomanip>

using namespace std;

CCorrUser2D::CCorrUser2D():hCorr(0),hErrl(0),hErrh(0){}
CCorrUser2D::~CCorrUser2D(){}

//--------------------------------------------------------------------------------------------------
// void CCorrUser2D::loadCorr(TH2D* corr, TH2D* errl, TH2D* errh)
void CCorrUser2D::loadCorr(TH2D* corr)
{
  hCorr = corr;
  // hErrl = errl;
  // hErrh = errh;
}

//--------------------------------------------------------------------------------------------------
Float_t CCorrUser2D::getCorr(const Double_t x, const Double_t y)
{
  if(!hCorr) {
    cout << "Corriciency table not loaded! Aborting..." << endl;
    assert(0);
  }
  return getValue(hCorr,x,y);
}

//--------------------------------------------------------------------------------------------------
Float_t CCorrUser2D::getErrLow(const Double_t x, const Double_t y)
{
  if(!hErrl) {
    cout << "Low errors table not loaded! Aborting..." << endl;
    assert(0);
  }
  return getValue(hErrl,x,y);
}

//--------------------------------------------------------------------------------------------------
Float_t CCorrUser2D::getErrHigh(const Double_t x, const Double_t y)
{
  if(!hErrh) {
    cout << "High errors table not loaded! Aborting..." << endl;
    assert(0);
  }
  return getValue(hErrh,x,y);
}

//--------------------------------------------------------------------------------------------------  
void CCorrUser2D::printCorr(ostream& os)
{
  if(!hCorr) {
    cout << "Corriciency table not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "Corriciency Table:" << endl;
  os << "-----------------" << endl;
  printHist2D(hCorr,os);   
}

//--------------------------------------------------------------------------------------------------
void CCorrUser2D::printErrLow(ostream& os)
{
  if(!hErrl) {
    cout << "Error (low) table not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "Low Errors Table:" << endl;
  os << "-----------------" << endl;
  printHist2D(hErrl,os); 
}

//--------------------------------------------------------------------------------------------------
void CCorrUser2D::printErrHigh(ostream& os)
{
  if(!hErrh) {
    cout << "Error (high) table not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "High Errors Table:" << endl;
  os << "------------------" << endl;
  printHist2D(hErrh,os);
}

//--------------------------------------------------------------------------------------------------
void CCorrUser2D::printHist2D(const TH2D* h, ostream& os) {
  const Int_t nx = h->GetNbinsX();
  const Int_t ny = h->GetNbinsY();
  
  for(Int_t iy=0; iy<=ny; iy++) {
    for(Int_t ix=0; ix<=nx; ix++) {
      if(ix==0 && iy==0) {
        os << setw(11) << "";
      } else if(ix==0) {
        os << "[" << setw(4) << h->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << h->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      } else if(iy==0) { 
        os << "[" << setw(4) << h->GetXaxis()->GetBinLowEdge(ix) << "," << setw(4) << h->GetXaxis()->GetBinLowEdge(ix+1) << "]";
      } else {
        ios_base::fmtflags flags = os.flags();
	os.precision(7);
	os << " " << setw(9) << fixed << h->GetBinContent(h->GetBin(ix,iy)) << " ";
	os.flags(flags);
      }
    }
    os << endl;
  }
}

//--------------------------------------------------------------------------------------------------
Float_t CCorrUser2D::getValue(const TH2D* h, const Double_t x, const Double_t y)
{
  Int_t ix=0;
  Int_t iy=0;
  const Int_t nx = h->GetNbinsX();
  const Int_t ny = h->GetNbinsY();
  
  for(Int_t i=1; i<=nx; i++) {
    if((x >= h->GetXaxis()->GetBinLowEdge(i)) && (x < h->GetXaxis()->GetBinLowEdge(i+1))) {
      ix=i;
      break;
    }
  }
  
  for(Int_t i=1; i<=ny; i++) {
    if((y >= h->GetYaxis()->GetBinLowEdge(i)) && (y < h->GetYaxis()->GetBinLowEdge(i+1))) {
      iy=i;
      break;
    }
  }
  
  if(ix>0 && iy>0)
    return h->GetBinContent(h->GetBin(ix,iy));
  else 
    return -1;
}
