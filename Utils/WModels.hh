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
  CPepeModel2isobins(const char *name, RooRealVar &x, double iso, RooRealVar *lin1=0, RooRealVar *lin2=0, RooRealVar *lin3=0, RooRealVar *off1=0, RooRealVar *off2=0, RooRealVar *off3=0);
  ~CPepeModel2isobins() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    // delete c1, c2, d1, d2;
    // delete a3;
    delete model;
  }
  RooRealVar *c1, *c2, *c3, *d1, *d2, *d3;
  // RooRealVar *a3;
  RooGenericPdf *model;
};

// Define the class for special Pepe2, with isolation dependence
class CPepeModel2isobinsConst
{
public:
  CPepeModel2isobinsConst():model(0){}
  CPepeModel2isobinsConst(const char *name, RooRealVar &x, double iso, RooRealVar *off1=0, RooRealVar *off2=0, RooRealVar *off3=0);
  ~CPepeModel2isobinsConst() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    // delete c1, c2, d1, d2;
    // delete a3;
    delete model;
  }
  RooRealVar *d1, *d2, *d3;
  // RooRealVar *a3;
  RooGenericPdf *model;
};

// Define the class for special Pepe2, with isolation dependence
class CPepeModel2isobinsQuad
{
public:
  CPepeModel2isobinsQuad():model(0){}
  CPepeModel2isobinsQuad(const char *name, RooRealVar &x, double iso, RooRealVar *qu1=0, RooRealVar *qu2=0, RooRealVar *qu3=0, RooRealVar *lin1=0, RooRealVar *lin2=0, RooRealVar *lin3=0, RooRealVar *off1=0, RooRealVar *off2=0, RooRealVar *off3=0);
  ~CPepeModel2isobinsQuad() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    // delete c1, c2, d1, d2;
    // delete a3;
    delete model;
  }
  RooRealVar *b1, *b2, *b3, *c1, *c2, *c3, *d1, *d2, *d3;
  // RooRealVar *a3;
  RooGenericPdf *model;
};

class CPepeModel2isobinsMuons
{
public:
  CPepeModel2isobinsMuons():model(0){}
  CPepeModel2isobinsMuons(const char *name, RooRealVar &x, double iso, RooRealVar *lin1=0, RooRealVar *lin2=0, RooRealVar *lin3=0, RooRealVar *off1=0, RooRealVar *off2=0, RooRealVar *off3=0);
  ~CPepeModel2isobinsMuons() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    // delete c1, c2, d1, d2;
    // delete a3;
    delete model;
  }
  RooRealVar *c1, *c2, *c3, *d1, *d2, *d3;
  // RooRealVar *a3;
  RooGenericPdf *model;
};

class CPepeModel2isobins2
{
public:
  CPepeModel2isobins2():model(0){}
  CPepeModel2isobins2(const char *name, RooRealVar &x, double iso, RooRealVar *lin2=0, RooRealVar *lin3=0, RooRealVar *off1=0, RooRealVar *off2=0, RooRealVar *off3=0);
  ~CPepeModel2isobins2() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    // delete c1, c2, d1, d2;
    // delete a3;
    delete model;
  }
  RooRealVar *c2, *c3, *d1, *d2, *d3;
  // RooRealVar *a3;
  RooGenericPdf *model;
};


// delete this?
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

CPepeModel2isobins::CPepeModel2isobins(const char *name, RooRealVar &x, double iso,  RooRealVar *lin1, RooRealVar *lin2, RooRealVar *lin3, RooRealVar *off1, RooRealVar *off2, RooRealVar *off3)
{
    
   std::cout << "New Pepe model with isolation = " << iso << std::endl;
  char c1Name[50]; 
  char d1Name[50]; 
  if(lin1 && off1) {
    sprintf(c1Name,"%s",lin1->GetName());
    sprintf(d1Name,"%s",off1->GetName());
    c1 = lin1;
    d1 = off1;
  } else {
    sprintf(c1Name,"c1_%s",name);
    sprintf(d1Name,"d1_%s",name);
    // c1 = new RooRealVar(c1Name, c1Name, -2.5, -10.0, 10.0);
    // d1 = new RooRealVar(d1Name, d1Name, 2.5, -2.0, 4.0);    
    c1 = new RooRealVar(c1Name, c1Name, 1.0, -100.0, 100.0);
    d1 = new RooRealVar(d1Name, d1Name, 0.5, -200.0, 400.0);
  }
  char c2Name[50]; 
  char d2Name[50]; 
  char a2Name[50];
  if(lin2 && off2) {
    sprintf(c2Name,"%s",lin2->GetName());
    sprintf(d2Name,"%s",off2->GetName());
    c2 = lin2;
    d2 = off2;
  } else {
    sprintf(d2Name,"d2_%s",name);
    sprintf(c2Name,"c2_%s",name);
    // c2 = new RooRealVar(c2Name, c2Name, 2.0, -10.0, 10.0);
    // d2 = new RooRealVar(d2Name, d2Name, 6.0, 3.0, 12.0);    
    c2 = new RooRealVar(c2Name, c2Name, -3.0, -100.0, 200.0);
    d2 = new RooRealVar(d2Name, d2Name, 8.0, -500.0, 200.0);
  }
  // char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  char c3Name[50]; 
  char d3Name[50]; 
  char a3Name[50];
  if(lin3 && off3) {
    sprintf(c3Name,"%s",lin3->GetName());
    sprintf(d3Name,"%s",off3->GetName());
    c3 = lin3;
    d3 = off3;
  } else {
    sprintf(d3Name,"d3_%s",name);
    sprintf(c3Name,"c3_%s",name);
    // c3 = new RooRealVar(c3Name, c3Name, 3.0, 0.0, 10.0);
    // d3 = new RooRealVar(d3Name, d3Name, 1.5, 0.5, 3.0);    
    c3 = new RooRealVar(c3Name, c3Name, 3.0, -30, 500.0);
    d3 = new RooRealVar(d3Name, d3Name, 1.0, -20, 60);
  }
  
  // char formula[300];
  // sprintf(formula,
          // "%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + %s*100))",
	  // x.GetName(),
	  // x.GetName(),x.GetName(),
	  // c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  // c2Name,iso,d2Name,x.GetName(),
	  // a3Name);
      
  char formula[300];
  sprintf(formula,
          "(%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + (%s*%f+%s)*100)))*(((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + (%s*%f+%s)*100) > 0)",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  c2Name,iso,d2Name,x.GetName(),
	  c3Name,iso,d3Name,
	  c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  c2Name,iso,d2Name,x.GetName(),
	  c3Name,iso,d3Name);

  std::cout << "the formula is  " << formula << std::endl;
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*c1,*c2,*c3,*d1,*d2,*d3));
}


//--------------------------------------------------------------------------------------------------
// Pepe2, with the a1 and a2 parameters as a function of isolation
// test with a1 first, then parametrize a2 as well

CPepeModel2isobinsConst::CPepeModel2isobinsConst(const char *name, RooRealVar &x, double iso,  RooRealVar *off1, RooRealVar *off2, RooRealVar *off3)
{
    
   std::cout << "New Pepe model with isolation = " << iso << std::endl;
  char d1Name[50]; 
  if(off1) {
    sprintf(d1Name,"%s",off1->GetName());
    d1 = off1;
  } else {
    sprintf(d1Name,"d1_%s",name);
    d1 = new RooRealVar(d1Name, d1Name, 0.5, -2.0, 4.0);
  }
  char d2Name[50]; 
  if(off2) {
    sprintf(d2Name,"%s",off2->GetName());
    d2 = off2;
  } else {
    sprintf(d2Name,"d2_%s",name);
    d2 = new RooRealVar(d2Name, d2Name, 8.0, 6.0, 12.0);
  }
  // char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  char d3Name[50]; 
  if(off3) {
    sprintf(d3Name,"%s",off3->GetName());
    d3 = off3;
  } else {
    sprintf(d3Name,"d3_%s",name);
    d3 = new RooRealVar(d3Name, d3Name, 1.0, 0.5, 3.0);
  }
  
  // char formula[300];
  // sprintf(formula,
          // "%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + %s*100))",
	  // x.GetName(),
	  // x.GetName(),x.GetName(),
	  // c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  // c2Name,iso,d2Name,x.GetName(),
	  // a3Name);
      
  char formula[300];
  sprintf(formula,
          "(%s*exp(-%s*%s/((%s)*%s*%s*0.01 + (%s)*%s + (%s)*100)))*(((%s)*%s*%s*0.01 + (%s)*%s + (%s)*100) > 0)",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  d1Name,x.GetName(),x.GetName(),
	  d2Name,x.GetName(),
	  d3Name,
	  d1Name,x.GetName(),x.GetName(),
	  d2Name,x.GetName(),
	  d3Name);

  std::cout << "the formula is  " << formula << std::endl;
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*d1,*d2,*d3));
}


CPepeModel2isobinsQuad::CPepeModel2isobinsQuad(const char *name, RooRealVar &x, double iso,  RooRealVar *qu1, RooRealVar *qu2, RooRealVar *qu3, RooRealVar *lin1, RooRealVar *lin2, RooRealVar *lin3, RooRealVar *off1, RooRealVar *off2, RooRealVar *off3)
{
    
   std::cout << "New Pepe model with isolation = " << iso << std::endl;
  char b1Name[50]; 
  char c1Name[50]; 
  char d1Name[50]; 
  if(lin1 && off1 && qu1) {
	std::cout << "initializing 1 params" << std::endl;
    sprintf(b1Name,"%s",qu1->GetName());
    sprintf(c1Name,"%s",lin1->GetName());
    sprintf(d1Name,"%s",off1->GetName());
	std::cout << "b1 name =  "<< b1Name << std::endl;
	std::cout << "c1 name =  "<< c1Name << std::endl;
    b1 = qu1;
    c1 = lin1;
    d1 = off1;
	std::cout << "done" << std::endl;
  } else {
    sprintf(b1Name,"b1_%s",name);
    sprintf(c1Name,"c1_%s",name);
    sprintf(d1Name,"d1_%s",name);
    b1 = new RooRealVar(b1Name, b1Name, 0.0, -15.0, 15.0);
    c1 = new RooRealVar(c1Name, c1Name, 1.0, -10.0, 20.0);
    d1 = new RooRealVar(d1Name, d1Name, 0.5, -4.0, 4.0);
  }
  char b2Name[50]; 
  char c2Name[50]; 
  char d2Name[50]; 
  if(lin2 && off2 && qu2) {
	std::cout << "initializing 2 params" << std::endl;
    sprintf(b2Name,"%s",qu2->GetName());
    sprintf(c2Name,"%s",lin2->GetName());
    sprintf(d2Name,"%s",off2->GetName());
	std::cout << "b2 name =  "<< b2Name << std::endl;
    b2 = qu2;
    c2 = lin2;
    d2 = off2;
	std::cout << "done" << std::endl;
  } else {
    sprintf(b2Name,"b2_%s",name);
    sprintf(d2Name,"d2_%s",name);
    sprintf(c2Name,"c2_%s",name);
    b2 = new RooRealVar(b2Name, b2Name, 0.0, -10.0, 15.0);
    c2 = new RooRealVar(c2Name, c2Name, -3.0, -10.0, 5.0);
    d2 = new RooRealVar(d2Name, d2Name, 8.0, 6.0, 12.0);
  }
  // char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  char b3Name[50]; 
  char c3Name[50]; 
  char d3Name[50]; 
  // char a3Name[50];
  if(lin3 && off3 && qu3) {
	std::cout << "initializing 3 params" << std::endl;
    sprintf(b3Name,"%s",qu3->GetName());
    sprintf(c3Name,"%s",lin3->GetName());
    sprintf(d3Name,"%s",off3->GetName());
	std::cout << "b3 name =  "<< b3Name << std::endl;
    b3 = qu3;
    c3 = lin3;
    d3 = off3;
	std::cout << "done" << std::endl;
  } else {
    sprintf(b3Name,"b3_%s",name);
    sprintf(d3Name,"d3_%s",name);
    sprintf(c3Name,"c3_%s",name);
    b3 = new RooRealVar(b3Name, b3Name, 0.0, -10.0, 10.0);
    c3 = new RooRealVar(c3Name, c3Name, 3.0, 0.0, 7.0);
    d3 = new RooRealVar(d3Name, d3Name, 1.0, 0.1, 3.0);
  }
  
  // char formula[300];
  // sprintf(formula,
          // "%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + %s*100))",
	  // x.GetName(),
	  // x.GetName(),x.GetName(),
	  // c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  // c2Name,iso,d2Name,x.GetName(),
	  // a3Name);
      
  char formula[500];
  sprintf(formula,
          "(%s*exp(-%s*%s/((%s*%f*%f+%s*%f+%s)*%s*%s*0.01 + (%s*%f*%f+%s*%f+%s)*%s + (%s*%f*%f+%s*%f+%s)*100)))*(((%s*%f*%f+%s*%f+%s)*%s*%s*0.01 + (%s*%f*%f+%s*%f+%s)*%s + (%s*%f*%f+%s*%f+%s)*100) > 0)",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  b1Name,iso,iso,c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  b2Name,iso,iso,c2Name,iso,d2Name,x.GetName(),
	  b3Name,iso,iso,c3Name,iso,d3Name,
	  b1Name,iso,iso,c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  b2Name,iso,iso,c2Name,iso,d2Name,x.GetName(),
	  b3Name,iso,iso,c3Name,iso,d3Name);

  std::cout << "the formula is  " << formula << std::endl;
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  // split RooArgSet into 2 lines because you can only put up to 9 arguments at once
  RooArgSet qcdVars(x,*b1,*b2,*b3,*c1,*c2,*c3);
  std::cout << "test 1 " << std::endl;
  qcdVars.add(RooArgSet(*d1,*d2,*d3));
  std::cout << "test 2 " << std::endl;
  model = new RooGenericPdf(vname,vname,formula,qcdVars);
  std::cout << "test 3 " << std::endl;
}


CPepeModel2isobinsMuons::CPepeModel2isobinsMuons(const char *name, RooRealVar &x, double iso,  RooRealVar *lin1, RooRealVar *lin2, RooRealVar *lin3, RooRealVar *off1, RooRealVar *off2, RooRealVar *off3)
{
    
   std::cout << "New Pepe model with isolation = " << iso << std::endl;
  char c1Name[50]; 
  char d1Name[50]; 
  if(lin1 && off1) {
    sprintf(c1Name,"%s",lin1->GetName());
    sprintf(d1Name,"%s",off1->GetName());
    c1 = lin1;
    d1 = off1;
  } else {
    sprintf(c1Name,"c1_%s",name);
    sprintf(d1Name,"d1_%s",name);
    // c1 = new RooRealVar(c1Name, c1Name, 1.0, -5.0, 15.0);
    c1 = new RooRealVar(c1Name, c1Name, 0.5, -10.0, 10.0);
    // d1 = new RooRealVar(d1Name, d1Name, 0.5, -5.0, 10.0);
    d1 = new RooRealVar(d1Name, d1Name, 4.0, -30.0, 30.0);
  }
  char c2Name[50]; 
  char d2Name[50]; 
  char a2Name[50];
  if(lin2 && off2) {
    sprintf(c2Name,"%s",lin2->GetName());
    sprintf(d2Name,"%s",off2->GetName());
    c2 = lin2;
    d2 = off2;
  } else {
    sprintf(d2Name,"d2_%s",name);
    sprintf(c2Name,"c2_%s",name);
    // c2 = new RooRealVar(c2Name, c2Name, -3.0, -10.0, 10.0);
    c2 = new RooRealVar(c2Name, c2Name, 0.1, -15.0, 10.0);
    // d2 = new RooRealVar(d2Name, d2Name, 8.0, 0.0, 15.0);
    d2 = new RooRealVar(d2Name, d2Name, 6.0, -40.0, 40.0);
  }
  // char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  char c3Name[50]; 
  char d3Name[50]; 
  char a3Name[50];
  if(lin3 && off3) {
    sprintf(c3Name,"%s",lin3->GetName());
    sprintf(d3Name,"%s",off3->GetName());
    c3 = lin3;
    d3 = off3;
  } else {
    sprintf(d3Name,"d3_%s",name);
    sprintf(c3Name,"c3_%s",name);
    // c3 = new RooRealVar(c3Name, c3Name, 3.0, -5.0, 10.0);
    c3 = new RooRealVar(c3Name, c3Name, 0.1, -10.0, 10.0);
    // d3 = new RooRealVar(d3Name, d3Name, 1.0, -5.0, 10.0);
    d3 = new RooRealVar(d3Name, d3Name, 6.0, -40.0, 40.0);
  }
  
  // char formula[300];
  // sprintf(formula,
          // "%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + %s*100))",
	  // x.GetName(),
	  // x.GetName(),x.GetName(),
	  // c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  // c2Name,iso,d2Name,x.GetName(),
	  // a3Name);
	  
	  //(%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + (%s*%f+%s)*100)))",//*(((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + (%s*%f+%s)*100) > 0)",
      
  char formula[300];
  sprintf(formula,
          "(%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + (%s*%f+%s)*100)))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  c2Name,iso,d2Name,x.GetName(),
	  c3Name,iso,d3Name);
	  //c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  //c2Name,iso,d2Name,x.GetName(),
	  //c3Name,iso,d3Name);

  std::cout << "the formula is  " << formula << std::endl;
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*c1,*c2,*c3,*d1,*d2,*d3));
}

CPepeModel2isobins2::CPepeModel2isobins2(const char *name, RooRealVar &x, double iso, RooRealVar *lin2, RooRealVar *lin3, RooRealVar *off1, RooRealVar *off2, RooRealVar *off3)
{
    
   std::cout << "New Pepe model with isolation = " << iso << std::endl;
  char c1Name[50]; 
  char d1Name[50]; 
  if(off1) {
    // sprintf(c1Name,"%s",lin1->GetName());
    sprintf(d1Name,"%s",off1->GetName());
    // c1 = lin1;
    d1 = off1;
  } else {
    // sprintf(c1Name,"c1_%s",name);
    sprintf(d1Name,"d1_%s",name);
    // c1 = new RooRealVar(c1Name, c1Name, 0.5, -10.0, 10.0);
    d1 = new RooRealVar(d1Name, d1Name, 4.0, -30.0, 30.0);
  }
  char c2Name[50]; 
  char d2Name[50]; 
  char a2Name[50];
  if(lin2 && off2) {
    sprintf(c2Name,"%s",lin2->GetName());
    sprintf(d2Name,"%s",off2->GetName());
    c2 = lin2;
    d2 = off2;
  } else {
    sprintf(d2Name,"d2_%s",name);
    sprintf(c2Name,"c2_%s",name);
    c2 = new RooRealVar(c2Name, c2Name, 0.1, -10.0, 10.0);
    d2 = new RooRealVar(d2Name, d2Name, 6.0, -40.0, 40.0);
  }
  // char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  char c3Name[50]; 
  char d3Name[50]; 
  char a3Name[50];
  if(lin3 && off3) {
    sprintf(c3Name,"%s",lin3->GetName());
    sprintf(d3Name,"%s",off3->GetName());
    c3 = lin3;
    d3 = off3;
  } else {
    sprintf(d3Name,"d3_%s",name);
    sprintf(c3Name,"c3_%s",name);
    c3 = new RooRealVar(c3Name, c3Name, 0.1, -10.0, 10.0);
    d3 = new RooRealVar(d3Name, d3Name, 6.0, -40.0, 40.0);
  }
  
  // char formula[300];
  // sprintf(formula,
          // "%s*exp(-%s*%s/((%s*%f+%s)*%s*%s*0.01 + (%s*%f+%s)*%s + %s*100))",
	  // x.GetName(),
	  // x.GetName(),x.GetName(),
	  // c1Name,iso,d1Name,x.GetName(),x.GetName(),
	  // c2Name,iso,d2Name,x.GetName(),
	  // a3Name);
      
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/((%s)*%s*%s*0.01 + (%s*%f+%s)*%s + (%s*%f+%s)*100))",
	  x.GetName(),
	  x.GetName(),x.GetName(),
	  d1Name,x.GetName(),x.GetName(),
	  c2Name,iso,d2Name,x.GetName(),
	  c3Name,iso,d3Name);

  std::cout << "the formula is  " << formula << std::endl;
  char vname[50];
  sprintf(vname,"pepe2Pdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*c2,*c3,*d1,*d2,*d3));
}


//--------------------------------------------------------------------------------------------------
// Pepe2, with the a1 and a2 parameters as a function of isolation
// test with a1 first, then parametrize a2 as well
// Delete this? 
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
