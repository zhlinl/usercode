#ifndef __RooPolarizationConstraint_h__
#define __RooPolarizationConstraint_h__

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "TH1.h"

class RooPolarizationConstraint: public RooAbsPdf {
 public:
  RooPolarizationConstraint();
  RooPolarizationConstraint(const char* name,
			    const char* title,			   
			    RooAbsReal& lambda_theta,
			    RooAbsReal& lambda_phi,
			    RooAbsReal& lambda_thetaphi);
  
  RooPolarizationConstraint(const RooPolarizationConstraint& other,
			    const char* name);
  virtual ~RooPolarizationConstraint();

  virtual TObject* clone(const char* name) const {
    return new RooPolarizationConstraint(*this,name);
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars,
			      RooArgSet& analVars,
			      const char*) const;

  Double_t analyticalIntegral(Int_t code,
			      const char* rN) const;

 protected:  
  Double_t evaluate() const;
  RooRealProxy lTheta, lPhi, lThetaPhi; // parameters
    
  ClassDef(RooPolarizationConstraint,1) // RooPolarizationConstraint
};

#endif
