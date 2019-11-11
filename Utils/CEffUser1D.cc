#include "CEffUser1D.hh"
#include <iomanip>

using namespace std;

CEffUser1D::CEffUser1D():graph(0){}
CEffUser1D::~CEffUser1D(){}

//--------------------------------------------------------------------------------------------------
void CEffUser1D::loadEff(TGraphAsymmErrors* gr)
{
  graph = gr;
}

//--------------------------------------------------------------------------------------------------
Float_t CEffUser1D::getEff(const Double_t x)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  return graph->GetY()[getBin(x)];
}

//--------------------------------------------------------------------------------------------------
Float_t CEffUser1D::getErrLow(const Double_t x)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  return graph->GetEYlow()[getBin(x)];
}

//--------------------------------------------------------------------------------------------------
Float_t CEffUser1D::getErrHigh(const Double_t x)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  return graph->GetEYhigh()[getBin(x)];
}

//--------------------------------------------------------------------------------------------------  
void CEffUser1D::printEff(ostream& os)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "Efficiency:" << endl;
  os << "-----------" << endl;
  print(graph->GetY(),os);
}

//--------------------------------------------------------------------------------------------------
void CEffUser1D::printErrLow(ostream& os)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "Low Errors:" << endl;
  os << "-----------" << endl;
  print(graph->GetEYlow(),os);
}

//--------------------------------------------------------------------------------------------------
void CEffUser1D::printErrHigh(ostream& os)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  os << "High Errors:" << endl;
  os << "------------" << endl;
  print(graph->GetEYhigh(),os);
}

//--------------------------------------------------------------------------------------------------
void CEffUser1D::print(const Double_t *yval, std::ostream& os)
{
  const Double_t *xval  = graph->GetX();
  const Double_t *xerrl = graph->GetEXlow();
  const Double_t *xerrh = graph->GetEXhigh();
  
  for(Int_t i=0; i<graph->GetN(); i++) {
    os << "[" << setw(4) << xval[i]-xerrl[i] << "," << setw(4) << xval[i]+xerrh[i] << "]";
  }
  os << endl;
  for(Int_t i=0; i<graph->GetN(); i++) {
    ios_base::fmtflags flags = os.flags();
    os.precision(7);
    os << " " << setw(9) << fixed << yval[i] << " ";
    os.flags(flags);
  }
  os << endl; 
}

//--------------------------------------------------------------------------------------------------
Int_t CEffUser1D::getBin(const Double_t x)
{
  if(!graph) {
    cout << "Efficiency graph not loaded! Aborting..." << endl;
    assert(0);
  }
  const Double_t *xval  = graph->GetX();
  const Double_t *xerrl = graph->GetEXlow();
  const Double_t *xerrh = graph->GetEXhigh();
  for(Int_t i=0; i<graph->GetN(); i++) {
    if((x >= (xval[i]-xerrl[i])) && (x < (xval[i]+xerrh[i]))) return i; 
  }
  return -1;
}
