#include "RooPolarizationPdf.h"
#include <string>
#include <math.h>
#include "RooFormulaVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "TString.h"

ClassImp(RooPolarizationPdf)

RooPolarizationPdf::RooPolarizationPdf():
  hasmaps(false),
  pi(atan2(0.,-1.)) {
  initialize();
}


RooPolarizationPdf::RooPolarizationPdf(const char* name,
				       const char* title,
				       RooAbsReal& _costh,
				       RooAbsReal& _phi,
				       RooAbsReal& lambda_theta,
				       RooAbsReal& lambda_phi,
				       RooAbsReal& lambda_thetaphi):
  RooAbsPdf(name,title),
  costh("costh","Dependent",this,_costh),
  phi("phi","Dependent",this,_phi),
  lTheta("lTheta","Lambda Theta",this,lambda_theta),
  lPhi("lPhi","Lambda Phi",this,lambda_phi),
  lThetaPhi("lThetaPhi","Lambda Theta-Phi",this,lambda_thetaphi),   
  hasmaps(false),
  pi(atan2(0.,-1.)) {
  initialize();
}

RooPolarizationPdf::RooPolarizationPdf(const char* name,
				       const char* title,
				       RooAbsReal& _costh,
				       RooAbsReal& _phi,
				       RooAbsReal& lambda_theta,
				       RooAbsReal& lambda_phi,
				       RooAbsReal& lambda_thetaphi,
				       RooAbsReal& _map_u,
				       RooAbsReal& _map_th,
				       RooAbsReal& _map_ph,
				       RooAbsReal& _map_thph):
  RooAbsPdf(name,title),
  costh("costh","Dependent",this,_costh),
  phi("phi","Dependent",this,_phi),
  lTheta("lTheta","Lambda Theta",this,lambda_theta),
  lPhi("lPhi","Lambda Phi",this,lambda_phi),  
  lThetaPhi("lThetaPhi","Lambda Theta-Phi",this,lambda_thetaphi),  
  map_uniform("map_uniform","Acceptance Map",this,_map_u),  
  map_theta("map_theta","Acceptance Map",this,_map_th),  
  map_phi("map_phi","Acceptance Map",this,_map_ph),  
  map_thetaphi("map_thetaphi","Acceptance Map",this,_map_thph),    
  hasmaps(true),
  pi(atan2(0.,-1.)) {
  initialize();
}

RooPolarizationPdf::RooPolarizationPdf(const RooPolarizationPdf& other,
				       const char* name):
  RooAbsPdf(other,name),
  costh("costh",this,other.costh),
  phi("phi",this,other.phi),
  lTheta("lTheta",this,other.lTheta),
  lPhi("lPhi",this,other.lPhi),
  lThetaPhi("lThetaPhi",this,other.lThetaPhi),  
  map_uniform("map_uniform",this,other.map_uniform),
  map_theta("map_theta",this,other.map_theta),
  map_phi("map_phi",this,other.map_phi),
  map_thetaphi("map_thetaphi",this,other.map_thetaphi),  
  hasmaps(other.hasmaps),
  pi(other.pi) { 
  initialize();
}

RooPolarizationPdf::~RooPolarizationPdf() {    
}

void RooPolarizationPdf::initialize() {
  //std::cout << "called initialize!" << std::endl;
  if(hasmaps) {    
    
    //std::cout << map_uniform << std::endl;
    //std::cout << map_theta << std::endl;
    //std::cout << map_phi << std::endl;
    //std::cout << map_thetaphi << std::endl;    

    //cache the normalizations of each map component!
    //uniform
    //std::cout << "make uniform integral" << std::endl;
    RooAbsReal* i = map_uniform.arg().createIntegral(RooArgSet(costh.arg(),phi.arg()));
    //std::cout << "made uniform integral" << std::endl;
    norm_uniform = i->getVal();
    delete i;

    //theta
    i = map_theta.arg().createIntegral(RooArgSet(costh.arg(),phi.arg()));
    norm_theta = i->getVal();
    delete i;

    //phi
    i = map_phi.arg().createIntegral(RooArgSet(costh.arg(),phi.arg()));
    norm_phi = i->getVal();
    delete i;

    //thetaphi
    i = map_thetaphi.arg().createIntegral(RooArgSet(costh.arg(),phi.arg()));
    norm_thetaphi = i->getVal();
    delete i;
    
    //std::cout << "Normalizations: " << norm_uniform << ' ' << norm_theta 
    //          << ' ' << norm_phi << ' ' << norm_thetaphi << std::endl;

  } else {    
    norm_uniform=norm_theta=norm_phi=norm_thetaphi=1;    
  }
}

Double_t RooPolarizationPdf::evaluate() const {
  Double_t result = 0.0;  
    
  result = 1. + lTheta*costh*costh 
    + lPhi*(1.-costh)*(1.+costh)*cos(2.*phi*pi/180.)
    + lThetaPhi*2.*sqrt((1.-costh)*(1.+costh))*costh*cos(phi*pi/180.);
  
  if( hasmaps ) {
    //std::cout << "prev result " << result << std::endl;
    result *= map_uniform;
    //std::cout << "after map modulation" << result << std::endl;
  }
  
  if(result < 0) return 1e-120;

  return result;
}

Int_t RooPolarizationPdf::getAnalyticalIntegral(RooArgSet& allVars,
						RooArgSet& analVars,
						const char*) const {
  if (matchArgs(allVars,analVars,costh,phi)) return 3;
  //if (matchArgs(allVars,analVars,phi))       return 2;
  //if (matchArgs(allVars,analVars,costh))     return 1;
  return 0;
}

Double_t RooPolarizationPdf::analyticalIntegral(Int_t code,
						const char* rN) const {  
  assert(code == 1 || code == 2 || code == 3);
  Double_t result = 0.0;

  Double_t b = costh.max(rN);
  Double_t a = costh.min(rN);
  Double_t p2 = phi.max(rN);
  Double_t p1 = phi.min(rN);
    
  if(hasmaps) {
    result = norm_uniform+norm_theta*lTheta+norm_phi*lPhi+norm_thetaphi*lThetaPhi;
  } else {
    result = coreIntegral(code,b,a,p2,p1);
  }
  
  return result;
}

Double_t RooPolarizationPdf::coreIntegral(Int_t code,
					  const Double_t& b, 
					  const Double_t& a,
					  const Double_t& p2, 
					  const Double_t& p1) const {    
  Double_t result = 0.;

  switch(code){
  case 3:
    {      
      Double_t r1mb2 = sqrt((1.-b)*(1.+b));
      Double_t r1ma2 = sqrt((1.-a)*(1.+a));
      Double_t r1mb23 = r1mb2*r1mb2*r1mb2;
      Double_t r1ma23 = r1ma2*r1ma2*r1ma2;
      
      Double_t dsin2th = sin(2.*p2*pi/180.)-sin(2.*p1*pi/180.);
      Double_t dsinth = sin(p2*pi/180.)-sin(p1*pi/180.);
      
      result = (b - a + lTheta*(b*b*b-a*a*a)/3)*(p2-p1)
	+ lPhi*(a*a*a - b*b*b + 3.*(b-a))*dsin2th*30./pi
	+ lThetaPhi*( -r1mb23 + r1ma23 )*dsinth*120/pi;
    }
    break;
  case 2:
    {
      Double_t dsin2th = sin(2.*p2*pi/180.)-sin(2.*p1*pi/180.);
      Double_t dsinth = sin(p2*pi/180.)-sin(p1*pi/180.);
      
      result = (1. + lTheta*costh*costh)*(p2 - p1)
	+ lPhi*(1.-costh)*(1.+costh)*dsin2th*90./pi
	+ lThetaPhi*sqrt((1.-costh)*(1.+costh))*costh*dsinth*180./pi;
    }
    break;
  case 1:
    {      
      Double_t r1mb2 = sqrt((1.-b)*(1.+b));
      Double_t r1ma2 = sqrt((1.-a)*(1.+a));
      Double_t r1mb23 = r1mb2*r1mb2*r1mb2;
      Double_t r1ma23 = r1ma2*r1ma2*r1ma2;
      
      result = b - a + lTheta*(b*b*b-a*a*a)/3.
	+ lPhi*(a*a*a - b*b*b + 3.*(b-a))*cos(2.*phi*pi/180.)/3.
	+ lThetaPhi*2.*( -r1mb23 + r1ma23 )/3.*cos(phi*pi/180.);
    }
    break;
  default:
    break;
  }

  return result;
}
