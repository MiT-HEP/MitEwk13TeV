#include "PdfDiagonalizer.h"


PdfDiagonalizer::PdfDiagonalizer(const char *name, RooWorkspace *w, RooFitResult &result) :
    // eigen_gl(new TVectorD()),
    name_(name),
    parameters_(result.floatParsFinal())
{
    int n = parameters_.getSize();
    
    TMatrixDSym cov(result.covarianceMatrix());
    std::cout << "covariance " << std::endl;
    TMatrixDSymEigen eigen(cov);
    std::cout << "eigenvecs " << std::endl;

    const TMatrixD& vectors = eigen.GetEigenVectors();
    const TVectorD& values  = eigen.GetEigenValues();

    char buff[10240];

    // create unit gaussians per eigen-vector
    for (int i = 0; i < n; ++i) {
        snprintf(buff,sizeof(buff),"%s_eig%d[-5,5]", name, i); // this is not actually a Gaussian??! -> Just creates a variable with given range, value=0
        eigenVars_.add(*w->factory(buff));
    }
    // put them in a list, with a one at the end to set the mean
    RooArgList eigvVarsPlusOne(eigenVars_);
    if (w->var("_one_") == 0) w->factory("_one_[1]");
    eigvVarsPlusOne.add(*w->var("_one_"));

	
    // then go create the linear combinations
    // each is equal to the transpose matrx times the square root of the eigenvalue (so that we get unit gaussians)
    for (int i = 0; i < n; ++i) {   
        RooArgList coeffs;
        for (int j = 0; j < n; ++j) {
            snprintf(buff,sizeof(buff),"%s_eigCoeff_%d_%d[%g]", name, i, j, vectors(i,j)*sqrt(values(j)));
            coeffs.add(*w->factory(buff)); 
        }
        snprintf(buff,sizeof(buff),"%s_eigBase_%d[%g]", name, i, (dynamic_cast<RooAbsReal*>(parameters_.at(i)))->getVal());
        coeffs.add(*w->factory(buff)); 
        snprintf(buff,sizeof(buff),"%s_eigLin_%d", name, i);
        RooAddition *add = new RooAddition(buff,buff,coeffs,eigvVarsPlusOne);
        w->import(*add,RooFit::Silence());
        replacements_.add(*add);
    }
}

RooAbsPdf *PdfDiagonalizer::diagonalize(RooAbsPdf &pdf)
{
  std::cout << "start diagonalize function " << std::endl;
  if (!pdf.dependsOn(parameters_)) return 0;

  // now do the customization
  RooCustomizer custom(pdf, name_.c_str());
  for (int i = 0, n = parameters_.getSize(); i < n; ++i) { 
    if (pdf.dependsOn(*parameters_.at(i))) {
      custom.replaceArg(*parameters_.at(i), *replacements_.at(i));
    }
  }

  RooAbsPdf *ret = dynamic_cast<RooAbsPdf *>(custom.build());
  ret->SetName((std::string(pdf.GetName()) + "_" + name_).c_str());
  std::cout << "print doDiago pdf value" << std::endl;
  ret->Print();
  return ret;
}

RooAbsPdf *PdfDiagonalizer::diagonalizeWithEigenVariations(RooAbsPdf &pdf, RooFitResult &result, int index_eig, int nSigma)
{
    
    // Check that the fit status is reasonable and that the covariance matrix is OK
    if(result.status() != 0 || result.covQual() !=3){return &pdf;}
    TMatrixDSym cov(result.covarianceMatrix()); 
    TMatrixDSymEigen eigen(cov);
    
    const TVectorD& values  = eigen.GetEigenValues();
    if (!pdf.dependsOn(eigenVars_)) return 0;

    RooArgList  eigenShiftedVars_; // these are the diagonal parameters    

    // now do the customization
    RooCustomizer custom(pdf, name_.c_str());
    for (int i = 0, n = eigenVars_.getSize(); i < n; ++i) { 
        if (pdf.dependsOn(*eigenVars_.at(i))) {
            
            RooRealVar *myVar = new RooRealVar("name","name",0);
            myVar->SetName(Form("%s_%i",eigenVars_.at(i)->GetName(),nSigma));
            double content=(dynamic_cast<RooAbsReal*>(eigenVars_.at(i)))->getVal();
            if(i==index_eig) content =+ nSigma*sqrt(values(i));
            myVar->setVal(content);
            eigenShiftedVars_.add(*myVar);

            custom.replaceArg(*eigenVars_.at(i), *eigenShiftedVars_.at(i));
        }
    }

    RooAbsPdf *ret = dynamic_cast<RooAbsPdf *>(custom.build());
    ret->SetName(std::string(pdf.GetName()).c_str());
    
    if(ret->getVal() <= 0 || ret->getVal() > 1|| std::isnan(ret->getVal()) || !std::isfinite(ret->getVal())) {
        return &pdf;
    }
    return ret;
}