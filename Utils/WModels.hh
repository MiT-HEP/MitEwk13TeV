#include <TH1D.h>
#include "RooRealVar.h"
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

class CPepeModelSq
{
public:
  CPepeModelSq():model(0){}
  CPepeModelSq(const char *name, RooRealVar &x,RooRealVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModelSq() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete a4;
    delete a5;
    delete model;
  }
  RooRealVar *a1, *a2, *a3, *a4, *a5;
  RooGenericPdf *model;
};

class CPepeModelGaussian
{
public:
  CPepeModelGaussian():model(0){}
  CPepeModelGaussian(const char *name, RooRealVar &x,RooRealVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModelGaussian() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete a4;
    delete model;
  }
  RooRealVar *a1, *a2, *a3, *a4;
  RooGenericPdf *model;
};

class CPepeModelBifGaus
{
public:
  CPepeModelBifGaus():model(0){}
  CPepeModelBifGaus(const char *name, RooRealVar &x,RooRealVar *sigma1=0, RooRealVar *sigma0=0);
  ~CPepeModelBifGaus() {
//     delete a1; // temporary fix to prevent segfault when using simultaneous fit
    delete a2;
    delete a3;
    delete m;
    delete sl;
    delete sr;
    delete model;
  }
  RooRealVar *a1, *a2, *a3, *m, *sl, *sr;
  RooAddPdf *model;
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
    a1 = new RooRealVar(a1Name,a1Name,4.0,-10.0,10.0);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma1;
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

CPepeModelSq::CPepeModelSq(const char *name, RooRealVar &x,  RooRealVar *sigma1, RooRealVar *sigma0)
{
  char a1Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    a1 = sigma1;
  } else {
    sprintf(a1Name,"a1_%s",name);
    a1 = new RooRealVar(a1Name,a1Name,4.0,-10.0,10.0);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma1;
  } else {
    sprintf(a2Name,"a2_%s",name);
    a2 = new RooRealVar(a2Name,a2Name,6.0,0.0,20);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.3,6.0);
  char a4Name[50]; sprintf(a4Name, "a4_%s", name); a4 = new RooRealVar(a4Name,a4Name,2.9,0.3,6.0);
  char a5Name[50]; sprintf(a5Name, "a5_%s", name); a5 = new RooRealVar(a5Name,a5Name,1.0,1500.0);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c]
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*0.01 + %s*%s + %s*100 + %s*%s*%s*%s*0.001-0.01*%s*sqrt(%s)))",
      x.GetName(),
      x.GetName(),x.GetName(),
      a1Name,x.GetName(),x.GetName(),
      a2Name,x.GetName(),
      a3Name,
      a4Name,x.GetName(),x.GetName(),x.GetName(),
      a5Name,x.GetName()
      );
  
  char vname[50];
  sprintf(vname,"pepeSqPdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3,*a4,*a5));
}

//--------------------------------------------------------------------------------------------------
CPepeModelGaussian::CPepeModelGaussian(const char *name, RooRealVar &x,  RooRealVar *sigma1, RooRealVar *sigma0)
{
  char a1Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    a1 = sigma1;
  } else {
    sprintf(a1Name,"a1_%s",name);
    a1 = new RooRealVar(a1Name,a1Name,4.0,-10.0,10.0);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma1;
  } else {
    sprintf(a2Name,"a2_%s",name);
    a2 = new RooRealVar(a2Name,a2Name,6.0,0.0,200);
    //a1 = new RooRealVar(a1Name,a1Name,0.01,-1,1);
  }
  //char a2Name[50]; sprintf(a2Name, "a2_%s", name); a2 = new RooRealVar(a2Name,a2Name,1,0,10);
  //char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,100,70,400);
  //char a2Name[50]; sprintf(a2Name, "a2_%s", name); a2 = new RooRealVar(a2Name,a2Name,1,0,20);
  //char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,200,50,600);
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.0,60.0);
  char a4Name[50]; sprintf(a4Name, "a4_%s", name); a4 = new RooRealVar(a4Name,a4Name,40.0,0.0,1000.0);
//   char a4Name[50]; sprintf(a4Name, "a4_%s", name); a4 = new RooRealVar(a4Name,a4Name,20,0.0,50.0);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c] * exp(-x^2/d)
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*0.01 + %s*%s + %s*100))*exp(-%s*%s/(%s*10))",
      x.GetName(),
      x.GetName(),x.GetName(),
      a1Name,x.GetName(),x.GetName(),
      a2Name,x.GetName(),
      a3Name,
      x.GetName(),x.GetName(),
      a4Name);
  std::cout << "Model: " << formula << std::endl;
  
  char vname[50];
  sprintf(vname,"pepeGPdf_%s",name);  
  model = new RooGenericPdf(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3,*a4));
}

CPepeModelBifGaus::CPepeModelBifGaus(const char *name, RooRealVar &x,  RooRealVar *sigma1, RooRealVar *sigma0)
{
  char a1Name[50]; 
  if(sigma1) {
    sprintf(a1Name,"%s",sigma1->GetName());
    a1 = sigma1;
  } else {
    sprintf(a1Name,"a1_%s",name);
    a1 = new RooRealVar(a1Name,a1Name,4.0,-10.0,10.0);
  }
  char a2Name[50];
  if(sigma0) {
    sprintf(a2Name,"%s",sigma0->GetName());
    a2 = sigma1;
  } else {
    sprintf(a2Name,"a2_%s",name);
    a2 = new RooRealVar(a2Name,a2Name,6.0,0.0,200);
  }
  char a3Name[50]; sprintf(a3Name, "a3_%s", name); a3 = new RooRealVar(a3Name,a3Name,2.9,0.0,60.0);
  
  // f(x) = x*exp[-x^2 / a*x*x + b*x + c] * exp(-x^2/d)
  char formula[300];
  sprintf(formula,
          "%s*exp(-%s*%s/(%s*%s*%s*0.01 + %s*%s + %s*100))",
      x.GetName(),
      x.GetName(),x.GetName(),
      a1Name,x.GetName(),x.GetName(),
      a2Name,x.GetName(),
      a3Name);//
//       x.GetName(),x.GetName(),
//       a4Name);
  std::cout << "Model: " << formula << std::endl;
  
  char vname[50];
  sprintf(vname,"pepeBiGPdf_%s",name);  
  RooGenericPdf model_pepe(vname,vname,formula,RooArgSet(x,*a1,*a2,*a3));
  
  char mName[50]; sprintf(mName, "m_%s", name); m = new RooRealVar(mName,mName,10.0,0.0,150.0);
  char srName[50]; sprintf(srName, "sr_%s", name); sr = new RooRealVar(srName,srName,10.0,0.0,300.0);
  char slName[50]; sprintf(slName, "sl_%s", name); sl = new RooRealVar(slName,slName,40.0,0.0,1000.0);
  
  RooBifurGauss model_bifur("model_bifur","model_bifur",x,*m,*sr,*sl);
  model = new RooAddPdf(vname, vname, RooArgList(model_bifur,model_pepe));
  
}

//--------------------------------------------------------------------------------------------------
CHistModel::CHistModel(const char *name, RooRealVar &x, TH1D* hist, int intOrder)
{
  char vname[100];
  
  sprintf(vname,"inHist_%s",name);   inHist   = (TH1D*)hist->Clone(vname);  
  sprintf(vname,"dataHist_%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(x),inHist);
  sprintf(vname,"histPdf_%s",name);  model    = new RooHistPdf(vname,vname,x,*dataHist,intOrder); 
}
