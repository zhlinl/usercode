#include "PolarizationModel.h"

#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "TDirectory.h" 

#include <sstream>

namespace JPsiPolarization {
  PolarizationModel::PolarizationModel(const std::string n,
				       const std::string& f) {
    compName_ = n;
    fName_ = f;
    varsLoaded = false;
    model_ = NULL;     
  }  

  PolarizationModel::~PolarizationModel() {
    if(!varsLoaded) {
      for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	  i != vars_.end();
	  ++i)
	delete i->second;
    }
    
    if(model_) delete model_;
    
  }

  void PolarizationModel::loadParameters(TDirectory& dir) {    
    std::stringstream theta_name,phi_name,thetaphi_name;
    std::stringstream theta_desc,phi_desc,thetaphi_desc;
    
    theta_name << compName_ << "lambda_theta_" << fName_;
    phi_name << compName_ << "lambda_phi_" << fName_;
    thetaphi_name << compName_ << "lambda_thetaphi_" << fName_;

    theta_desc << "#lamdba_{#theta} in the " << fName_ << " Frame";
    phi_desc << "#lamdba_{#phi} in the " << fName_ << " Frame";
    thetaphi_desc << "#lamdba_{#theta #phi} " << fName_ << " Frame" ;  
    
    dir.GetObject(theta_name.str().c_str(),vars_[theta_name.str().c_str()]);
    if(!vars_[theta_name.str().c_str()]) 
      vars_[theta_name.str().c_str()] = new RooRealVar(theta_name.str().c_str(),
						       theta_desc.str().c_str(),
						       0,-2,2);
    
    dir.GetObject(phi_name.str().c_str(),vars_[phi_name.str().c_str()]);
    if(!vars_[phi_name.str().c_str()]) 
      vars_[phi_name.str().c_str()] = new RooRealVar(phi_name.str().c_str(),
						     phi_desc.str().c_str(),
						     0,-2,2);

    dir.GetObject(thetaphi_name.str().c_str(),vars_[thetaphi_name.str().c_str()]);
    if(!vars_[thetaphi_name.str().c_str()]) 
      vars_[thetaphi_name.str().c_str()] = new RooRealVar(thetaphi_name.str().c_str(),
							  thetaphi_desc.str().c_str(),
							  0);

    varsLoaded = true;
  }

  void PolarizationModel::saveParameters(TDirectory& dir) {
    for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	i != vars_.end();
	++i) {
      dir.Add(i->second->Clone(),true);
    }    
    dir.Write();
  }

  void PolarizationModel::fix(const std::string& s) {
    for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	i != vars_.end();
	++i) {
      if(!s.size() || i->first.find(s) != std::string::npos) {
	std::cout << "Fixing: ";
	i->second->Print();
	i->second->setConstant(true);
      }
    }
  }

  void PolarizationModel::unfix(const std::string& s) {
    for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	i != vars_.end();
	++i) {
      if(!s.size() || i->first.find(s) != std::string::npos) {
	std::cout << "Unfixing: ";
	i->second->Print();
	i->second->setConstant(false);
      }
    }
  }

  void PolarizationModel::initModel(RooRealVar& costh,
				    RooRealVar& phi) {
    std::stringstream theta_name,phi_name,thetaphi_name;
    std::stringstream theta_desc,phi_desc,thetaphi_desc;

    std::stringstream modelName,modelDesc;
    std::stringstream polFunc;
    
    theta_name << compName_ << "lambda_theta_" << fName_;
    phi_name << compName_ << "lambda_phi_" << fName_;
    thetaphi_name << compName_ << "lambda_thetaphi_" << fName_;
    
    theta_desc << "#lamdba_{#theta} in the " << fName_ << " Frame";
    phi_desc << "#lamdba_{#phi} in the " << fName_ << " Frame";
    thetaphi_desc << "#lamdba_{#theta #phi} " << fName_ << " Frame" ;    
    
    RooRealVar *lambda_theta; 
    RooRealVar *lambda_phi; 
    RooRealVar *lambda_thetaphi;

    if(varsLoaded) {
      lambda_theta = vars_[theta_name.str().c_str()];
      lambda_phi = vars_[phi_name.str().c_str()];
      lambda_thetaphi = vars_[thetaphi_name.str().c_str()];
    } else {    
      lambda_theta = new RooRealVar(theta_name.str().c_str(),theta_desc.str().c_str(),0,-2,2);
      vars_[theta_name.str().c_str()] = lambda_theta;
      
      lambda_phi = new RooRealVar(phi_name.str().c_str(),phi_desc.str().c_str(),0,-2,2);
      vars_[phi_name.str().c_str()] = lambda_phi;
      
      lambda_thetaphi = new RooRealVar(thetaphi_name.str().c_str(),thetaphi_desc.str().c_str(),0);
      vars_[thetaphi_name.str().c_str()] = lambda_thetaphi;
    }


    modelName << compName_ <<"polFunc_" << fName_;
    modelDesc << compName_ << " Polarization PDF for " << fName_ << " Frame";
    polFunc << "3/(4*pi*(3+" << compName_ << "lambda_theta_" << fName_ 
	    << "))*(1+" << compName_ << "lambda_theta_" << fName_
	    << "*pow(costh_" << fName_
	    << ",2)+" << compName_ << "lambda_phi_" << fName_ 
	    << "*(1-pow(costh_" << fName_
	    << ",2))*cos(2*(phi_" << fName_ 
	    << "/180*pi))+" << compName_ << "lambda_thetaphi_" << fName_ 
	    << "*sin(2*acos(costh_" << fName_ 
	    << "))*cos(phi_" << fName_ << "/180*pi))";

    model_ = new RooGenericPdf(modelName.str().c_str(),
			      modelDesc.str().c_str(),
			      polFunc.str().c_str(),
			      RooArgSet(*lambda_theta,*lambda_phi,*lambda_thetaphi,costh,phi)) ;
  }
}

