#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooCMSShape.h"
#include "RooGenericPdf.h"

class CBackgroundModel
{
public:
  CBackgroundModel():model(0){}
  virtual ~CBackgroundModel() { delete model; }
  RooAbsPdf *model;
};

class CExponential : public CBackgroundModel
{
public:
  CExponential(RooRealVar &m, const Bool_t pass);
  CExponential(RooRealVar &m, const Bool_t pass, Double_t tVal);
  ~CExponential();
  RooRealVar *t;
};

class CErfExpo : public CBackgroundModel
{
public:
  CErfExpo(RooRealVar &m, const Bool_t pass);
  ~CErfExpo();
  RooRealVar *alfa, *beta, *gamma, *peak; 
};

class CDoubleExp : public CBackgroundModel
{
public:
  CDoubleExp(RooRealVar &m, const Bool_t pass);
  ~CDoubleExp();
  RooExponential *exp1, *exp2;
  RooRealVar *t1, *t2, *frac;
};

class CLinearExp : public CBackgroundModel
{
public:
  CLinearExp(RooRealVar &m, const Bool_t pass);
  ~CLinearExp();
  RooRealVar *a, *t;
};

class CQuadraticExp : public CBackgroundModel
{
public:
  CQuadraticExp(RooRealVar &m, const Bool_t pass);
  ~CQuadraticExp();
  RooRealVar *a1, *a2, *t;
};

class CMCTemplateConvGaussian2 : public CBackgroundModel
{
public:
  CMCTemplateConvGaussian2(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *sigma0=0, int intOrder=1);
  ~CMCTemplateConvGaussian2();
  RooRealVar  *bkgmean, *bkgsigma;
  RooGaussian *bkggaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
  };

//--------------------------------------------------------------------------------------------------
CExponential::CExponential(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  sprintf(vname,"t%s",name);
  if(pass)
    t = new RooRealVar(vname,vname,-0.1,-1.,1.);
  else
    t = new RooRealVar(vname,vname,-0.1,-1.,1.);
      
  sprintf(vname,"background%s",name);
  model = new RooExponential(vname,vname,m,*t);
}

CExponential::CExponential(RooRealVar &m, const Bool_t pass, Double_t tVal)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  sprintf(vname,"t%s",name);
  if(pass)
    t = new RooRealVar(vname,vname,-0.1,-1.,0.);
  else
    t = new RooRealVar(vname,vname,-0.1,-1.,0.);

  t->setVal(tVal);
  t->setConstant(kTRUE);
      
  sprintf(vname,"background%s",name);
  model = new RooExponential(vname,vname,m,*t);
}

CExponential::~CExponential()
{
  delete t;
}

//--------------------------------------------------------------------------------------------------
CErfExpo::CErfExpo(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
  
  if(pass) {
    sprintf(vname,"alfa%s",name);  alfa  = new RooRealVar(vname,vname,50,5,200);
    sprintf(vname,"beta%s",name);  beta  = new RooRealVar(vname,vname,0.01,0,10);
    sprintf(vname,"gamma%s",name); gamma = new RooRealVar(vname,vname,0.1,0,1);
  } else {
    sprintf(vname,"alfa%s",name);  alfa  = new RooRealVar(vname,vname,50,5,200);
    sprintf(vname,"beta%s",name);  beta  = new RooRealVar(vname,vname,0.01,0,10);
    sprintf(vname,"gamma%s",name); gamma = new RooRealVar(vname,vname,0.1,0,1);
  }  
  
  sprintf(vname,"peak%s",name);  
  peak = new RooRealVar(vname,vname,91.1876,85,97); 
  peak->setVal(91.1876);
  peak->setConstant(kTRUE);  
  
  sprintf(vname,"background%s",name);
  model = new RooCMSShape(vname,vname,m,*alfa,*beta,*gamma,*peak);
}

CErfExpo::~CErfExpo()
{
  delete alfa;
  delete beta;
  delete gamma;
  delete peak;
}

//--------------------------------------------------------------------------------------------------
CDoubleExp::CDoubleExp(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];
 
  if(pass) {
    sprintf(vname,"t1%s",name);   t1   = new RooRealVar(vname,vname,-0.20,-1.,0.);
    sprintf(vname,"t2%s",name);   t2   = new RooRealVar(vname,vname,-0.05,-1.,0.);
    sprintf(vname,"frac%s",name); frac = new RooRealVar(vname,vname, 0.50, 0.,1.);
  } else {
    sprintf(vname,"t1%s",name);   t1   = new RooRealVar(vname,vname,-0.20,-1.,0.);
    sprintf(vname,"t2%s",name);   t2   = new RooRealVar(vname,vname,-0.05,-1.,0.);
    sprintf(vname,"frac%s",name); frac = new RooRealVar(vname,vname, 0.50, 0.,1.);
  }
    
  sprintf(vname,"exp1%s",name);
  exp1 = new RooExponential(vname,vname,m,*t1);
  sprintf(vname,"exp2%s",name);
  exp2 = new RooExponential(vname,vname,m,*t2);
  sprintf(vname,"background%s",name);
  model = new RooAddPdf(vname,vname,RooArgList(*exp1,*exp2),RooArgList(*frac));
}

CDoubleExp::~CDoubleExp()
{
  delete exp1;
  delete exp2;
  delete t1;
  delete t2;
  delete frac;
}

//--------------------------------------------------------------------------------------------------
CLinearExp::CLinearExp(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char aname[50];
  sprintf(aname,"a%s",name);
  a = new RooRealVar(aname,aname,-0,-100.,100.);
  //a->setConstant(kTRUE);
  
  char tname[50];
  sprintf(tname,"t%s",name);
  t = new RooRealVar(tname,tname,-1e-6,-10.,10.);
  //t->setConstant(kTRUE); 
  
  char formula[200];
  sprintf(formula,"(1+%s*m)*exp(%s*m)",aname,tname);
 
  char vname[50]; sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,formula,RooArgList(m,*a,*t));
}

CLinearExp::~CLinearExp()
{
  delete a;
  delete t;
}

//--------------------------------------------------------------------------------------------------
CQuadraticExp::CQuadraticExp(RooRealVar &m, const Bool_t pass)
{
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");

  char a1name[50]; 
  sprintf(a1name,"a1%s",name);
  a1 = new RooRealVar(a1name,a1name,0,-10,10.);
  //a1->setConstant(kTRUE);
  
  char a2name[50]; 
  sprintf(a2name,"a2%s",name);
  a2 = new RooRealVar(a2name,a2name,0.0,-10,10);
  //a2->setConstant(kTRUE);
  
  char tname[50];
  sprintf(tname,"t%s",name);
  t = new RooRealVar(tname,tname,-1e-6,-10.,1.); 
  //t->setConstant(kTRUE); 
  
  char formula[200];
  sprintf(formula,"(1+%s*m+%s*m*m)*exp(%s*m)",a1name,a2name,tname);
 
  char vname[50]; sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,formula,RooArgList(m,*a1,*a2,*t));
}

CQuadraticExp::~CQuadraticExp()
{
  delete a1;
  delete a2;
  delete t;
}

CMCTemplateConvGaussian2::CMCTemplateConvGaussian2(RooRealVar &m, TH1D* hist, const Bool_t pass, RooRealVar *sigma0, int intOrder)
{  
  char name[10];
  if(pass) sprintf(name,"%s","Pass");
  else     sprintf(name,"%s","Fail");
  
  char vname[50];  
  
  if(pass) {
    sprintf(vname,"bkgmean%s",name);  bkgmean  = new RooRealVar(vname,vname,0,-10,10);
    if(sigma0) { bkgsigma = sigma0; }
    else { sprintf(vname,"bkgsigma%s",name); bkgsigma = new RooRealVar(vname,vname,2,0,15); }
    sprintf(vname,"bkggaus%s",name);  bkggaus  = new RooGaussian(vname,vname,m,*bkgmean,*bkgsigma);
  } else {
    sprintf(vname,"bgkmean%s",name); bkgmean  = new RooRealVar(vname,vname,0,-30,30);
    if(sigma0) { bkgsigma = sigma0; }
    else { sprintf(vname,"bkgsigma%s",name); bkgsigma = new RooRealVar(vname,vname,2,0,35); }
    sprintf(vname,"bkggaus%s",name);  bkggaus  = new RooGaussian(vname,vname,m,*bkgmean,*bkgsigma);
  }


  sprintf(vname,"inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);
  
  sprintf(vname,"dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
  sprintf(vname,"histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);
  sprintf(vname,"background%s",name);   model    = new RooFFTConvPdf(vname,vname,m,*histPdf,*bkggaus);
}

CMCTemplateConvGaussian2::~CMCTemplateConvGaussian2()
{
  delete bkgmean;
  delete bkgsigma;
  delete bkggaus;
  delete inHist;
  delete dataHist;
  delete histPdf;
}
