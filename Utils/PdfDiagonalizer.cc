#include "PdfDiagonalizer.h"


PdfDiagonalizer::PdfDiagonalizer(const char *name, RooWorkspace *w, RooFitResult &result) :
    // eigen_gl(new TVectorD()),
    name_(name),
    parameters_(result.floatParsFinal())
{
    // std::cout << "INSIDE DIAGONALIZER" << std::endl;
    int n = parameters_.getSize();

    // std::cout << "parameters_.getSize()= " << parameters_.getSize() << std::endl;
    // parameters_.Print();
    // std::cout << "parameters_.at(0).GetName()= " << parameters_.at(0)->GetName() << std::endl;
    
    
    TMatrixDSym cov(result.covarianceMatrix());
    std::cout << "covariance " << std::endl;
	// cov.Print();
    TMatrixDSymEigen eigen(cov);
    std::cout << "eigenvecs " << std::endl;
    // eigen.GetEigenVectors().Print();

    const TMatrixD& vectors = eigen.GetEigenVectors();
    const TVectorD& values  = eigen.GetEigenValues();
    // values.Print();

    char buff[10240];

    // create unit gaussians per eigen-vector
    for (int i = 0; i < n; ++i) {
        // snprintf(buff,sizeof(buff),"%s_eig%s[-5,5]", name, parameters_.at(i)->GetName());
        snprintf(buff,sizeof(buff),"%s_eig%d[-5,5]", name, i); // this is not actually a Gaussian??! -> Just creates a variable with given range, value=0
        // snprintf(buff,sizeof(buff),"RooGaussian::%s_eig%d(x[-5,5],mean[0],sigma[1])", name, i);
        eigenVars_.add(*w->factory(buff));
		// eigenVars_.at(i)->Print();
		// std::cout << "norm = "<< ((RooGaussian*)eigenVars_.at(i))->getNorm() << std::endl;
    }
    // put them in a list, with a one at the end to set the mean
    RooArgList eigvVarsPlusOne(eigenVars_);
    if (w->var("_one_") == 0) w->factory("_one_[1]");
    eigvVarsPlusOne.add(*w->var("_one_"));
	// eigvVarsPlusOne.Print();

	
    // then go create the linear combinations
    // each is equal to the transpose matrx times the square root of the eigenvalue (so that we get unit gaussians)
    for (int i = 0; i < n; ++i) {   
        RooArgList coeffs;
        for (int j = 0; j < n; ++j) {
            // snprintf(buff,sizeof(buff),"%s_eigCoeff_%s_%d[%g]", name, parameters_.at(i)->GetName(), j, vectors(i,j)*sqrt(values(j)));
            snprintf(buff,sizeof(buff),"%s_eigCoeff_%d_%d[%g]", name, i, j, vectors(i,j)*sqrt(values(j)));
            coeffs.add(*w->factory(buff)); 
			// coeffs.at(j)->Print();
        }
        // snprintf(buff,sizeof(buff),"%s_eigBase_%s[%g]", name, parameters_.at(i)->GetName(), (dynamic_cast<RooAbsReal*>(parameters_.at(i)))->getVal());
        snprintf(buff,sizeof(buff),"%s_eigBase_%d[%g]", name, i, (dynamic_cast<RooAbsReal*>(parameters_.at(i)))->getVal());
        coeffs.add(*w->factory(buff)); 
		// coeffs.at(i)->Print();
        // snprintf(buff,sizeof(buff),"%s_eigLin_%s", name, parameters_.at(i)->GetName());
        snprintf(buff,sizeof(buff),"%s_eigLin_%d", name, i);
        RooAddition *add = new RooAddition(buff,buff,coeffs,eigvVarsPlusOne);
		// add->Print();
        w->import(*add,RooFit::Silence());
        replacements_.add(*add);
        // replacements_.Print();
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
            // myVar->SetName(Form("%s",eigenVars_.at(i)->GetName()));
            // myVar->setVal(nSigma*sqrt(values(i)));
            //    myVar.setVal(*eigenVars_->at(i) + nSigma*values(i));     

            double content=(dynamic_cast<RooAbsReal*>(eigenVars_.at(i)))->getVal();
            // std::cout << " name= " << myVar->GetName() << "  myVar= " <<  myVar->getVal() << "  central value= " << content;
            // content =+ nSigma*sqrt(values(i));
            if(i==index_eig) content =+ nSigma*sqrt(values(i));
            // std::cout << "  value inserted= " << content << std::endl;
            myVar->setVal(content);

            
            // std::cout << " name=" << myVar->GetName() << " myVar= " <<  myVar->getVal() << "central value=  " << *eigenVars_.at(i) << std::endl;
            eigenShiftedVars_.add(*myVar);

            custom.replaceArg(*eigenVars_.at(i), *eigenShiftedVars_.at(i));
            // std::cout << "eigen var new val = " << (*eigenVars_.at(i)).getVal(); << std::endl;
        }
    }
    // std::cout << "test 2" << std::endl;
    // eigenVars_.Print();

    RooAbsPdf *ret = dynamic_cast<RooAbsPdf *>(custom.build());
    ret->SetName(std::string(pdf.GetName()).c_str());
    // std::cout << "value = " << ret->getVal() << std::endl;
    // std::cout << "log = " << ret->getLogVal() << std::endl;
    // std::cout << "eval = " << ret->evaluate() << std::endl;
	// for(int i = 0; i < )
    
    if(ret->getVal() <= 0 || ret->getVal() > 1|| std::isnan(ret->getVal()) || !std::isfinite(ret->getVal())) {
        // std::cout << "blah returning PDF " << std::endl;
        return &pdf;
    }
    // std::cout << "print var" << std::endl;
	// ret->Print();
    // std::cout << "norm of eigenvar " << ret->getNorm() << std::endl;;
    return ret;
}