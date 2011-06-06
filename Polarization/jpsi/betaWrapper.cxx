#include "Math/DistFunc.h"
#include "RooRealVar.h"
#ifndef __CINT__
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h"


RooAbsPdf* makeBetaPdf(const char* name,
		       RooRealVar& x, 
		       RooRealVar& alpha, 
		       RooRealVar& beta) {
  using namespace RooFit;
  return bindPdf(name,ROOT::Math::beta_pdf,x,alpha,beta);
}
