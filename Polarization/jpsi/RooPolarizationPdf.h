#ifndef __RooPolarizationPdf_h__
#define __RooPolarizationPdf_h__

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooHistFunc.h"
#include "RooFormulaVar.h"
#include "RooRealProxy.h"

class RooPolarizationPdf: public RooAbsPdf {
 public:
  RooPolarizationPdf();
  RooPolarizationPdf(const char* name,
		     const char* title,
		     RooAbsReal& _costh,  RooAbsReal& _phi, //observables
		     RooAbsReal& lambda_theta,
		     RooAbsReal& lambda_phi,
		     RooAbsReal& lambda_thetaphi);
  RooPolarizationPdf(const char* name,
		     const char* title,
		     RooAbsReal& _costh, RooAbsReal& _phi, //observables
		     RooAbsReal& lambda_theta,
		     RooAbsReal& lambda_phi,
		     RooAbsReal& lambda_thetaphi,
		     RooAbsReal& map_uniform,
		     RooAbsReal& map_theta,
		     RooAbsReal& map_phi,
		     RooAbsReal& map_thetaphi); //Acceptance Map

  RooPolarizationPdf(const RooPolarizationPdf& other,
		     const char* name);
  virtual ~RooPolarizationPdf();

  virtual TObject* clone(const char* name) const {
    return new RooPolarizationPdf(*this,name);
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, 
			      RooArgSet& analVars,
			      const char* n="") const;
  Double_t analyticalIntegral(Int_t code,const char* rN="") const;
  
  

 protected:  
  void initialize();
  Double_t coreIntegral(Int_t code,
			const Double_t& b, 
			const Double_t& a,
			const Double_t& p2, 
			const Double_t& p1) const;
  
  Double_t evaluate() const;
  RooRealProxy costh, phi; // observables
  RooRealProxy lTheta, lPhi, lThetaPhi; // parameters 
  RooRealProxy map_uniform, map_theta, map_phi, map_thetaphi; // acceptance * efficiency map
  Double_t norm_uniform, norm_theta, norm_phi, norm_thetaphi;  
  bool hasmaps;

  const Double_t pi;  
  ClassDef(RooPolarizationPdf,1) // RooPolarizationPdf
};

#endif
