#include "TROOT.h"
#include "TH1D.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooVoigtianShape.h"
#include "RooKeysPdf.h"

class CSignalModel
{
public:
  CSignalModel():model(0){}
  virtual ~CSignalModel(){ delete model; }
  RooAbsPdf *model;
};

class CBreitWignerConvCrystalBall : public CSignalModel
{
public:
  CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass);
  ~CBreitWignerConvCrystalBall();
  RooRealVar     *mass, *width;
  RooBreitWigner *bw;
  RooRealVar     *mean, *sigma, *alpha, *n;
  RooCBShape     *cb;
};

class CMCTemplateConvGaussian : public CSignalModel
{
public:
  CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *sigma0=0, int intOrder=1);
  ~CMCTemplateConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
};

class CVoigtianCBShape : public CSignalModel
{
public:
  CVoigtianCBShape(RooRealVar &m, const Bool_t pass);
  ~CVoigtianCBShape();
  RooRealVar *mass, *width;
  RooRealVar *sigma, *alpha, *n;
};

class CMCDatasetConvGaussian : public CSignalModel
{
public:
  CMCDatasetConvGaussian(RooRealVar &m, TTree* tree, const Bool_t pass, RooRealVar *sigma0=0);
  ~CMCDatasetConvGaussian();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TTree       *inTree;
  RooDataSet  *dataSet;
  RooKeysPdf  *keysPdf;
};

//--------------------------------------------------------------------------------------------------
CBreitWignerConvCrystalBall::CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  sprintf(vname,"mass%s",name);
  mass = new RooRealVar(vname,vname,91,80,100);    
  mass->setVal(91.1876);
  mass->setConstant(kTRUE);
  
  sprintf(vname,"width%s",name);
  width = new RooRealVar(vname,vname,2.5,0.1,10);    
  width->setVal(2.4952);
  width->setConstant(kTRUE);
  
  sprintf(vname,"bw%s",name);
  bw = new RooBreitWigner(vname,vname,m,*mass,*width);

  if(pass) {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,5);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"n%s",name);     n     = new RooRealVar(vname,vname,1,0,10);
  } else {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,5);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"n%s",name);     n     = new RooRealVar(vname,vname,1,0,10);
  }  
//  n->setVal(1.0);
//  n->setConstant(kTRUE);
  
  sprintf(vname,"cb%s",name);
  cb = new RooCBShape(vname,vname,m,*mean,*sigma,*alpha,*n);
        
  sprintf(vname,"signal%s",name);
  model = new RooFFTConvPdf(vname,vname,m,*bw,*cb);
}

CBreitWignerConvCrystalBall::~CBreitWignerConvCrystalBall()
{
  delete mass;
  delete width;
  delete bw;
  delete mean;
  delete sigma;
  delete alpha;
  delete n;
  delete cb;
}

//--------------------------------------------------------------------------------------------------
CMCTemplateConvGaussian::CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *sigma0, int intOrder)
{  
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];  
  
  if(pass) {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  } else {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  }


  sprintf(vname,"inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);
  
  sprintf(vname,"dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
  sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);
  sprintf(vname,"signal%s",name);   model    = new RooFFTConvPdf(vname,vname,m,*histPdf,*gaus);
}

CMCTemplateConvGaussian::~CMCTemplateConvGaussian()
{
  delete mean;
  //delete sigma;
  delete gaus;
  delete inHist;
  delete dataHist;
  delete histPdf;
}

//--------------------------------------------------------------------------------------------------
CVoigtianCBShape::CVoigtianCBShape(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  sprintf(vname,"mass%s",name);
  mass = new RooRealVar(vname,vname,91,80,100);
    
  sprintf(vname,"width%s",name);
  width = new RooRealVar(vname,vname,2.5,0.1,10);    
  width->setVal(2.4952);
  width->setConstant(kTRUE);
  
  if(pass) {
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,10);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
  } else {
    sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,10);
    sprintf(vname,"alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
  }
  
  sprintf(vname,"n%s",name);     
  n = new RooRealVar(vname,vname,1,0,10);
  n->setVal(1.0);
  
  sprintf(vname,"signal%s",name);
  model = new RooVoigtianShape(vname,vname,m,*mass,*sigma,*alpha,*n,*width,0);  
}

CVoigtianCBShape::~CVoigtianCBShape()
{
  delete mass;
  delete width;
  delete sigma;
  delete alpha;
  delete n;
}

//--------------------------------------------------------------------------------------------------
CMCDatasetConvGaussian::CMCDatasetConvGaussian(RooRealVar &m, TTree* tree, const Bool_t pass, RooRealVar *sigma0)
{  
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];  

  if(pass) {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  } else {
    sprintf(vname,"mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { sigma = sigma0; }
    else       { sprintf(vname,"sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5); }
    sprintf(vname,"gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);
  }
  
  sprintf(vname,"inTree_%s",tree->GetName());
  inTree = (TTree*)tree->Clone(vname);
  
  sprintf(vname,"dataSet%s",name); dataSet = new RooDataSet(vname,vname,inTree,RooArgSet(m));
  sprintf(vname,"keysPdf%s",name); keysPdf = new RooKeysPdf(vname,vname,m,*dataSet,RooKeysPdf::NoMirror,1);
  sprintf(vname,"signal%s",name);  model   = new RooFFTConvPdf(vname,vname,m,*keysPdf,*gaus);
}

CMCDatasetConvGaussian::~CMCDatasetConvGaussian()
{
  delete mean;
  delete sigma;
  delete gaus;
  delete inTree;
  delete dataSet;
  delete keysPdf;
}
