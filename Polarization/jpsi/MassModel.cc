#include "MassModel.h"

#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "TDirectory.h"

namespace JPsiPolarization {
  MassModel::MassModel() {
    varsLoaded = false;
    promptModel_ = NULL;
    nonpromptModel_ = NULL;
    bkgModel_ = NULL;
  }  

  MassModel::~MassModel() {
    if(!varsLoaded) {
      for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	  i != vars_.end();
	  ++i)
	delete i->second;
    }
    
    if(promptModel_) delete promptModel_;
    if(nonpromptModel_) delete nonpromptModel_;
    if(bkgModel_) delete bkgModel_;
  }

  void MassModel::loadParameters(TDirectory& dir) {
    if(!varsLoaded) {
      for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	  i != vars_.end();
	  ++i) {	
	delete i->second;
      }    
    }
    
    dir.GetObject("CBm",vars_["CBm"]);
    if(!vars_["CBm"]) vars_["CBm"] = new RooRealVar("CBm","Prompt Model Mass",3.1,3.05,3.15);
    
    dir.GetObject("CBs",vars_["CBs"]);
    if(!vars_["CBs"]) vars_["CBs"] = new RooRealVar("CBs","Prompt Model Sigma",0.02,0.008,0.1);
    
    dir.GetObject("CBa",vars_["CBa"]);	
    if(!vars_["CBa"]) vars_["CBa"] = new RooRealVar("CBa","Prompt Model CB Tail Parameter",0.5,0.,3);
    
    dir.GetObject("CBn",vars_["CBn"]);
    if(!vars_["CBn"]) vars_["CBn"] = new RooRealVar("CBn","Prompt Model CB Power Law Parameter",10,0,30);
      
    dir.GetObject("bkgLambda",vars_["bkgLambda"]);
    if(!vars_["bkgLambda"]) vars_["bkgLambda"] =  new RooRealVar("bkgLambda","Background Decay Parameter",-1,-3,0.1);
    
    /*
    dir.GetObject("FracGauss",vars_["FracGauss"]);
    if(!vars_["FracGauss"]) vars_["FracGauss"] =  new RooRealVar("FracGauss","Prompt Model CB/Gaussian Fraction",0.99,0,1);

    dir.GetObject("gSig",vars_["gSig"]);
    if(!vars_["gSig"]) vars_["gSig"] = new RooRealVar("gSig","Prompt Model Sigma",.1,0.0,1.0);
    */
    varsLoaded = true;
  }

  void MassModel::saveParameters(TDirectory& dir) {
    
    for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	i != vars_.end();
	++i) {
      dir.Add(i->second->Clone(),true);
    }    
    dir.Write();
  }

  void MassModel::fix(const std::string& s) {
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

  void MassModel::unfix(const std::string& s) {
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

  void MassModel::initModels(RooRealVar& var,
			     const bool& p,
			     const bool& np,
			     const bool& bkg) {
    if(p) initPromptModel(var);
    if(np) initNonPromptModel(var);
    if(bkg) initBackgroundModel(var);
  }

  void MassModel::initPromptModel(RooRealVar& var) {
    RooRealVar *CBm;
    RooRealVar *CBs;
    RooRealVar *CBa;
    RooRealVar *CBn;
    RooRealVar *frac;
    RooRealVar *gSig;

    if(varsLoaded) {
      CBm  = vars_["CBm"];
      CBs  = vars_["CBs"];
      CBa  = vars_["CBa"];
      CBn  = vars_["CBn"];
      //frac = vars_["FracGauss"];
      //gSig = vars_["gSig"];
    } else {
      CBm = new RooRealVar("CBm","Prompt Model Mass",3.1,3.05,3.15);
      vars_["CBm"] = CBm;
      CBs = new RooRealVar("CBs","Prompt Model Sigma",0.02,0.008,0.1);
      vars_["CBs"] = CBs;
      CBa = new RooRealVar("CBa","Prompt Model CB Tail Parameter",0.5,0.,3);
      vars_["CBa"] = CBa;
      CBn = new RooRealVar("CBn","Prompt Model CB Power Law Parameter",10,0,30);
      vars_["CBn"] = CBn;
      /*
      frac = new RooRealVar("FracGauss","Prompt Model CB/Gaussian Fraction",0.99,0,1);
      vars_["FracGauss"] = frac;
      gSig = new RooRealVar("gSig","Prompt Model Sigma",.1,0.0,1.0);
      vars_["gSig"] = gSig;
      */
    }

    RooCBShape *cb = new RooCBShape("promptMassCBShape","Prompt Mass CB Shape",
				    var,
				    *CBm,*CBs,*CBa,*CBn);
    pdfs_.push_back(cb);    
    /*
    RooGaussian *gaus = new RooGaussian("promptMassGaus","Prompt Mass Gaussian",
					var,
					*CBm,*gSig);
    pdfs_.push_back(gaus);
    */
    promptModel_ = new RooCBShape("promptMassCBShape","Prompt Mass CB Shape",
				  var,
				  *CBm,*CBs,*CBa,*CBn);//new RooAddPdf("promptMassModel","Prompt Mass PDF", RooArgList(*cb,*gaus),*frac);
  }

  void MassModel::initNonPromptModel(RooRealVar& var) {
        
    RooRealVar *CBm = vars_["CBm"];
    RooRealVar *CBs = vars_["CBs"];
    RooRealVar *CBa = vars_["CBa"];
    RooRealVar *CBn = vars_["CBn"];
    //RooRealVar *frac = vars_["FracGauss"];
    

    RooCBShape *cb = new RooCBShape("nonPromptMassCBShape","Non-Prompt Mass CB Shape",
				    var,
				    *CBm,*CBs,*CBa,*CBn);
    pdfs_.push_back(cb);
    /*
    RooRealVar *gSig = vars_["gSig"];

    RooGaussian *gaus = new RooGaussian("nonPromptMassGaus","Non-Prompt Mass Gaussian",
					var,
					*CBm,*gSig);
    pdfs_.push_back(gaus);
    */
    nonpromptModel_ = new RooCBShape("nonPromptMassCBShape","Non-Prompt Mass CB Shape",
				     var,
				     *CBm,*CBs,*CBa,*CBn);//new RooAddPdf("nonPromptMassModel","Non-Prompt Mass PDF", RooArgList(*cb,*gaus),*frac);
  }

  void MassModel::initBackgroundModel(RooRealVar& var) {
    RooRealVar *lambda,*p0,*p1,*p2;

    if(varsLoaded) {
      lambda = vars_["bkgLambda"];
    } else {
      
      lambda = new RooRealVar("bkgLambda","Background Decay Parameter",
			      0,-3,3);
      vars_["bkgLambda"] = lambda;
      /*
      vars_["p0"] = new RooRealVar("p0","First Order Chebychev Coefficient",0,-20,20);      
      vars_["p1"] = new RooRealVar("p1","Second Order Chebychev Coefficient",0,-20,20);      
      vars_["p2"] = new RooRealVar("p2","Third Order Chebychev Coefficient",0,-20,20);      
      */
    }

    bkgModel_ = new RooExponential("bkgMassModel","Background Mass PDF",var,*lambda);
    //bkgModel_ = new RooChebychev("bkgMassModel","Background Mass PDF",var,
    //			 RooArgList(*vars_["p0"],*vars_["p1"],*vars_["p2"]));
  }
}
