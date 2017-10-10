#include <TH1D.h>
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooBifurGauss.h"

class CPepeModel0
{
public:
  CPepeModel0():model(0){}
  CPepeModel0(const char *name, RooRealVar &x);
  ~CPepeModel0() {
    delete sigma;
    delete model;
  }
  RooRealVar *sigma;
  RooGenericPdf *model;
};

class CPepeModel1
{
public:
  CPepeModel1():model(0){}
  CPepeModel1(const char *name, RooRealVar &x, RooRealVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModel1() {
    //delete sigma;
    //delete a1;
    delete model;
  }
  RooRealVar *sigma, *a1;
  RooGenericPdf *model;
};

class CPepeModel2
{
public:
  CPepeModel2():model(0){}
  CPepeModel2(const char *name, RooRealVar &x,RooRealVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModel2() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete model;
  }
  RooRealVar *a1, *a2, *a3;
  RooGenericPdf *model;
};

// Define the class for special Pepe2, with isolation dependence
class CPepeModel2isobins
{
public:
  CPepeModel2isobins():model(0){}
  CPepeModel2isobins(const char *name, RooRealVar &x, RooRealVar &iso, RooFormulaVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModel2isobins() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete model;
  }
  RooFormulaVar *a1;//, *a2;
  RooRealVar *a3, *a2, *c1, *c2, *isoVal;
  // RooRealVar c1("c1", "c1",0,  -5.0, 5.0);
  // RooRealVar c2("c2", "c2",0,  -5.0, 5.0);
  RooGenericPdf *model;
};

class CPepeModel2dim
{
public:
  CPepeModel2dim():model(0){}
  CPepeModel2dim(const char *name, RooRealVar &x, RooRealVar &iso, RooFormulaVar *sigma1=0, RooFormulaVar *sigma0=0);
  ~CPepeModel2dim() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete model;
  }
  // The iso-dependent variable (formula var)
  RooFormulaVar *a1, *a2;//, *a3;
  // The parameters to define the dependence on isolation for each of the above:
  RooRealVar *c1, *c2;
  RooRealVar *d1, *d2;
  // RooRealVar *e1, *e2;
  RooRealVar *a3;//, *a2;
  // RooRealVar c1("c1", "c1",0,  -5.0, 5.0);
  // RooRealVar c2("c2", "c2",0,  -5.0, 5.0);
  RooGenericPdf *model;
};

class CPepeModel3
{
public:
  CPepeModel3():model(0){}
  CPepeModel3(const char *name, RooRealVar &x,RooRealVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModel3() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete a4;
    delete model;
  }
  RooRealVar *a1, *a2, *a3, *a4;
  RooGenericPdf *model;
};

class CHistModel
{
public:
  CHistModel():model(0){}
  CHistModel(const char *name, RooRealVar &x, TH1D* hist, int intOrder=1);
  virtual ~CHistModel() { 
    delete inHist;
    delete dataHist;
    delete model;
  }
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *model;
};


//--------------------------------------------------------------------------------------------------
CPepeModel0::CPepeModel0(const char *name, RooRealVar &x)
{
  char sigmaName[50]; sprintf(sigmaName,"sigma_%s",name); sigma = new RooRealVar(sigmaName,sigmaName,10,5,20);
  
  // f(x) = x*exp[-x^2 / s^2]
  char formula[200];
  sprintf(formula,
          "(%s/%s/%s)*exp(-%s*%s/%s/%s)",
	  x.GetName(),sigmaName,sigmaName,
	  x.GetName(),x.GetName(),
	  sigmaName,sigmaName);
  
  char vname[50];
  sprintf(vname,"pepe0Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*sigma));
}

//--------------------------------------------------------------------------------------------------
CPepeModel1::CPepeModel1(const char *name, RooRealVar &x, RooRealVar *sigma1,RooRealVar *sigma0)
{
  char sigmaName[50];
  if(sigma0) {
    sprintf(sigmaName,"%s",sigma0->GetName());
    sigma = sigma0;
  } else {
    sprintf(sigmaName,"sigma_%s",name);
    sigma = new RooRealVar(sigmaName,sigmaName,10,1,30);
  }

  //sprintf(sigmaName,"sigma_%s",name);
  //sigma = new RooRealVar(sigmaName,sigmaName,10,1,20);
  
  char a1Name[50];
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    a1 = sigma1;
  } else {
    sprintf(a1Name,"a1_%s",name);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
    a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  
  // f(x) = x*exp[-x^2 / (s + ax)^2] = x*exp[-x*x/(s*s + 2*a*x + a*a*x*x)]
  char formula[200];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s + 2*%s*%s*%s + %s*%s*%s*%s))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  sigmaName,sigmaName,
	  sigmaName,a1Name,x.GetName(),
	  a1Name,a1Name,x.GetName(),x.GetName());
  
  char vname[50];
  sprintf(vname,"pepe1Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*sigma,*a1));
}

//--------------------------------------------------------------------------------------------------
/*CPepeModel2::CPepeModel2(const char *name, RooRealVar &x)
{
  char a1Name[50]; sprintf(a1Name, "a1_%s", name); a1 = new RooRealVar(a1Name,a1Name,0.1,-1,1);
  //char a2Name[50]; sprintf(a2Name, "a2_%s", name); a2 = new RooRealVar(a2Name,a2Name,1,0,10);
  //char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,100,70,400);
  char a2Name[50]; sprintf(a2Name, "a2_%s", name); a2 = new RooRealVar(a2Name,a2Name,1,0,20);
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,100,70,600);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c]
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s + %s*%s + %s))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  a1Name,x.GetName(),x.GetName(),
	  a2Name,x.GetName(),
	  a3Name);
  
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3));
}
*/
//--------------------------------------------------------------------------------------------------
CPepeModel2::CPepeModel2(const char *name, RooRealVar &x,  RooRealVar *sigma1, RooRealVar *sigma0)
{
  char a1Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    a1 = sigma1;
  } else {
    sprintf(a1Name,"a1_%s",name);
    a1 = new RooRealVar(a1Name,a1Name,4.0,-10.0,20.0);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma0;
  } else {
    sprintf(a2Name,"a2_%s",name);
    a2 = new RooRealVar(a2Name,a2Name,6.0,0.0,20);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c]
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*0.01 + %s*%s + %s*100))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  a1Name,x.GetName(),x.GetName(),
	  a2Name,x.GetName(),
	  a3Name);
  
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3));
}

//--------------------------------------------------------------------------------------------------
// Pepe2, with the a1 and a2 parameters as a function of isolation
// test with a1 first, then parametrize a2 as well
CPepeModel2isobins::CPepeModel2isobins(const char *name, RooRealVar &x, RooRealVar &iso,  RooFormulaVar *sigma1, RooRealVar *sigma0)
{
  char a1Name[50]; 
  char c1Name[50]; 
  char c2Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    sprintf(c1Name,"c1_%s",name);
    sprintf(c2Name,"c2_%s",name);
    c1 = (RooRealVar*)sigma1->getParameter(0);
    c2 = (RooRealVar*)sigma1->getParameter(1);
    // c2 = new RooRealVar("c2", "c2", 0, -5.0, 5.0);
    // isoVal=new RooRealVar("iso","iso",0);
    // isoVal->setVal(iso);
    a1 = new RooFormulaVar(a1Name,a1Name,"@0*@2+@1",RooArgSet(*c1,*c2,iso));
    a1 = sigma1; 
  } else {
    sprintf(a1Name,"a1_%s",name);
    sprintf(c1Name,"c1_%s",name);
    sprintf(c2Name,"c2_%s",name);
    // isoVal = new RooRealVar("iso","iso",0);
    // isoVal->setVal(iso);
    c1 = new RooRealVar("c1", "c1", 0, -5.0, 5.0);
    c2 = new RooRealVar("c2", "c2", 0, -5.0, 5.0);
    a1 = new RooFormulaVar(a1Name,a1Name,"@0*@2+@1",RooArgSet(*c1,*c2,iso));
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma0;
  } else {
    sprintf(a2Name,"a2_%s",name);
    // Keep a2 as is for now, update to RooFormulaVar once a1 works
    // a2 = new RooFormulaVar(a2Name,a2Name,6.0,0.0,20);
    a2 = new RooRealVar(a2Name,a2Name,6.0,0.0,20);
  }
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c]
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*0.01 + %s*%s + %s*100))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  a1Name,x.GetName(),x.GetName(),
	  a2Name,x.GetName(),
	  a3Name);
  
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3));
}


//--------------------------------------------------------------------------------------------------
// Pepe2, with the a1 and a2 parameters as a function of isolation
// test with a1 first, then parametrize a2 as well
CPepeModel2dim::CPepeModel2dim(const char *name, RooRealVar &x, RooRealVar &iso,  RooFormulaVar *sigma1, RooFormulaVar *sigma0)
{
  char a1Name[50]; 
  char c1Name[50]; 
  char c2Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    //sprintf(c1Name,"c1_%s",name);
   // sprintf(c2Name,"c2_%s",name);
    //c1 = (RooRealVar*)sigma1->getParameter(0);
    //c2 = (RooRealVar*)sigma1->getParameter(1);
    // c2 = new RooRealVar("c2", "c2", 0, -5.0, 5.0);
    // isoVal=new RooRealVar("iso","iso",0);
    // isoVal->setVal(iso);
    //a1 = new RooFormulaVar(a1Name,a1Name,"@0*@2+@1",RooArgSet(*c1,*c2,iso));
    a1 = sigma1; 
  } else {
    sprintf(a1Name,"a1_%s",name);
    sprintf(c1Name,"c1_%s",name);
    sprintf(c2Name,"c2_%s",name);
    // isoVal = new RooRealVar("iso","iso",0);
    // isoVal->setVal(iso);
    c1 = new RooRealVar(c1Name, c1Name, 2.0, 0.0, 10.0);
    c2 = new RooRealVar(c2Name, c2Name, 0.5, -20.0, 20.0);
    a1 = new RooFormulaVar(a1Name,a1Name,"@0*@2+@1",RooArgSet(*c1,*c2,iso));
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];  
  char d1Name[50]; 
  char d2Name[50]; 
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma0;
  } else {
    sprintf(a2Name,"a2_%s",name);
    // Keep a2 as is for now, update to RooFormulaVar once a1 works
    // a2 = new RooFormulaVar(a2Name,a2Name,6.0,0.0,20);
    // a2 = new RooRealVar(a2Name,a2Name,6.0,0.0,20);
	
    sprintf(d1Name,"d1_%s",name);
    sprintf(d2Name,"d2_%s",name);
    // isoVal = new RooRealVar("iso","iso",0);
    // isoVal->setVal(iso);
    d1 = new RooRealVar(d1Name, d1Name, 10.0, -10.0, 30.0);
    d2 = new RooRealVar(d2Name, d2Name, -1.0, -20.0, 20.0);
    a2 = new RooFormulaVar(a2Name,a2Name,"@0*@2+@1",RooArgSet(*d1,*d2,iso));
  }
  char a3Name[50]; 
  // char e1Name[50]; 
  // char e2Name[50]; 
  sprintf(a3Name, "a3_%s", name); 
  a3 = new RooRealVar(a3Name,a3Name,5.0,-10,6.0);
  // sprintf(e1Name,"e1_%s",name);
  // sprintf(e2Name,"e2_%s",name);
  // e1 = new RooRealVar(e1Name, e1Name, 0, -50.0, 20.0);
  // e2 = new RooRealVar(e2Name, e2Name, 0, -50.0, 20.0);
  // a3 = new RooFormulaVar(a3Name,a3Name,"@0*@2+@1",RooArgSet(*e1,*e2,iso));
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c]
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*0.01 + %s*%s + %s*100))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  a1Name,x.GetName(),x.GetName(),
	  a2Name,x.GetName(),
	  a3Name);
  
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3));
}


CPepeModel3::CPepeModel3(const char *name, RooRealVar &x,  RooRealVar *sigma1, RooRealVar *sigma0)
{
  char a1Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    a1 = sigma1;
  } else {
    sprintf(a1Name,"a1_%s",name);
    a1 = new RooRealVar(a1Name,a1Name,0.5,0.0,6.0);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma1;
  } else {
    sprintf(a2Name,"a2_%s",name);
    a2 = new RooRealVar(a2Name,a2Name,-7.0,-10.0,10.0);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,10.0,0.0,20);
  char a4Name[50]; sprintf(a4Name, "a4_%s", name); a4 = new RooRealVar(a4Name,a4Name,2.9,0.3,6.0);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c]
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*%s*0.001 + %s*%s*%s*0.01 + %s*%s + %s*100))",
      x.GetName(),
      x.GetName(),x.GetName(),
      a1Name,x.GetName(),x.GetName(),x.GetName(),
      a2Name,x.GetName(),x.GetName(),
      a3Name,x.GetName(),
      a4Name
      );
  
  char vname[50];
  sprintf(vname,"pepe3Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3,*a4));
}

//--------------------------------------------------------------------------------------------------
CHistModel::CHistModel(const char *name, RooRealVar &x, TH1D* hist, int intOrder)
{
  char vname[100];
  
  sprintf(vname,"inHist_%s",name);   inHist   = (TH1D*)hist->Clone(vname);  
  sprintf(vname,"dataHist_%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(x),inHist);
  sprintf(vname,"histPdf_%s",name);  model    = new RooHistPdf(vname,vname,x,*dataHist,intOrder); 
}
