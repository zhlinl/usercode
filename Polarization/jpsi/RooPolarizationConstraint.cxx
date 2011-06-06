#include "RooPolarizationConstraint.h"
#include <string>
#include <cmath>

ClassImp(RooPolarizationConstraint)

RooPolarizationConstraint::RooPolarizationConstraint() {
}


RooPolarizationConstraint::RooPolarizationConstraint(const char* name,
						     const char* title,
						     RooAbsReal& lambda_theta,
						     RooAbsReal& lambda_phi,
						     RooAbsReal& lambda_thetaphi):
  RooAbsPdf(name,title),
  lTheta("lTheta","Lambda Theta",this,lambda_theta),
  lPhi("lPhi","Lambda Phi",this,lambda_phi),
  lThetaPhi("lThetaPhi","Lambda Theta-Phi",this,lambda_thetaphi)
{}

RooPolarizationConstraint::RooPolarizationConstraint(const RooPolarizationConstraint& other,
						     const char* name):
  RooAbsPdf(other,name),
  lTheta("lTheta",this,other.lTheta),
  lPhi("lPhi",this,other.lPhi),
  lThetaPhi("lThetaPhi",this,other.lThetaPhi)
{}

RooPolarizationConstraint::~RooPolarizationConstraint() {
}

Int_t RooPolarizationConstraint::getAnalyticalIntegral(RooArgSet& allVars,
						       RooArgSet& analVars,
						       const char*) const {
  return 1;
}

Double_t RooPolarizationConstraint::analyticalIntegral(Int_t code,
						       const char* rN) const {
  assert(code == 1);
  return 1.0;
}

Double_t RooPolarizationConstraint::evaluate() const {

  Double_t e1 = std::min(0.,1+lTheta-2*fabs(lPhi));
  Double_t e2 = std::min(0.,1-lPhi-2*fabs(lThetaPhi));
  Double_t e3 = lPhi <= -1./3. ? std::min(0.,1-((1.+2*lPhi)*(1.+2*lPhi) + 2.*lThetaPhi*lThetaPhi )) : 0;

  Double_t r = std::exp(-100.*(e1*e1+e2*e2+e3*e3)); 

  return r;
}


