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
  CExponential(RooRealVar &m, const Bool_t pass, const int ibin);
  ~CExponential();
  RooRealVar *t;
};

class CErfExpo : public CBackgroundModel
{
public:
  CErfExpo(RooRealVar &m, const Bool_t pass, const int ibin);
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
  CQuadraticExp(RooRealVar &m, const Bool_t pass, const int ibin);
  ~CQuadraticExp();
  RooRealVar *a1, *a2, *t;
};

class CQuadratic : public CBackgroundModel
{
public:
  CQuadratic(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2);
  ~CQuadratic();
  RooRealVar *a0, *a1, *a2;
};

class CPower : public CBackgroundModel
{
public:
  CPower(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0);
  ~CPower();
  RooRealVar *t;
};

//--------------------------------------------------------------------------------------------------
CExponential::CExponential(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  
  char vname[50];
  
  sprintf(vname,"t%s",name);
  if(pass)
    t = new RooRealVar(vname,vname,-0.1,-1.,0.);
  else
    t = new RooRealVar(vname,vname,-0.1,-1.,0.);
      
  sprintf(vname,"background%s",name);
  model = new RooExponential(vname,vname,m,*t);
}

CExponential::~CExponential()
{
  delete t;
}

//--------------------------------------------------------------------------------------------------
CErfExpo::CErfExpo(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  
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
  a = new RooRealVar(aname,aname,-0,-10.,10.);
  //a->setConstant(kTRUE);
  
  char tname[50];
  sprintf(tname,"t%s",name);
  t = new RooRealVar(tname,tname,-1e-6,-10.,0.);
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
CQuadraticExp::CQuadraticExp(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);

  char a1name[50]; 
  sprintf(a1name,"a1%s",name);
  a1 = new RooRealVar(a1name,a1name,0.,-10,10.);
  
  char a2name[50]; 
  sprintf(a2name,"a2%s",name);
  a2 = new RooRealVar(a2name,a2name,0.,-10,10);
  
  char tname[50];
  sprintf(tname,"t%s",name);
  t = new RooRealVar(tname,tname,-1e-6,-10.,0.); 
  
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

//--------------------------------------------------------------------------------------------------
CQuadratic::CQuadratic(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2)
{ 
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);

  char a0name[50];
  sprintf(a0name,"a0%s",name); 
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
    a0 = new RooRealVar(a0name,a0name,p0,p0-e0,p0+e0);
  }else{
    a0 = new RooRealVar(a0name,a0name,0.,-10.,10.);
  }

  char a1name[50]; 
  sprintf(a1name,"a1%s",name);
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
    a1 = new RooRealVar(a1name,a1name,p1,p1-e1,p1+e1);
  }else{
    a1 = new RooRealVar(a1name,a1name,0.,-10.,10.);
  }

  char a2name[50]; 
  sprintf(a2name,"a2%s",name);
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
      a2 = new RooRealVar(a2name,a2name,p2,p2-e2,p2+e2);
  }else{
    a2 = new RooRealVar(a2name,a2name,0.,-10.,10.);
  }
  
  char formula[200];
  sprintf(formula,"(%s+%s*m+%s*m*m)",a0name, a1name,a2name);
  
  char vname[50]; sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,formula,RooArgList(m,*a0,*a1,*a2));
}

CQuadratic::~CQuadratic()
{
  delete a0;
  delete a1;
  delete a2;
}

//--------------------------------------------------------------------------------------------------
CPower::CPower(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);

  char vname[50];

  sprintf(vname,"t%s",name);
  if((p0!=0.)||(e0!=0.)){
    t = new RooRealVar(vname,vname,p0,p0-e0,p0+e0);
  }else{
    t = new RooRealVar(vname,vname,1.0,0.05,2.0);
  }

  sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,"pow(@0,-@1)",RooArgSet(m,*t));
}

CPower::~CPower()
{
  delete t;
}

