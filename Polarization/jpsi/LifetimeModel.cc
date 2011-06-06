#include "LifetimeModel.h"

#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooGaussModel.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooHistPdf.h"
#include "RooDecay.h"
#include "RooDataHist.h"
#include "RooRealConstant.h"
#include "RooConstVar.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TDirectory.h"

namespace JPsiPolarization {

  LifetimeModel::LifetimeModel() {
    isTemplate_ = false;
    varsLoaded = false;
    promptModel_ = NULL;
    nonPromptModel_ = NULL;
    bkgModel_ = NULL;
  }

  LifetimeModel::~LifetimeModel() {
    
    for(std::map<std::string,RooDataHist*>::iterator i = templates_.begin();
	i != templates_.end();
	++i) 
      delete i->second;

    for(std::map<std::string,RooAddModel*>::iterator i = addmodels_.begin();
	i != addmodels_.end();
	++i) 
      delete i->second;
    
    for(std::vector<RooDecay*>::iterator i = dmodels_.begin();
	i != dmodels_.end();
	++i) 
      delete *i;

    for(std::vector<RooGaussModel*>::iterator i = gmodels_.begin();
	i != gmodels_.end();
	++i) 
      delete *i;

    for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	i != vars_.end();
	++i) {
      delete i->second;
    }

    if(promptModel_) delete promptModel_;
    if(nonPromptModel_) delete nonPromptModel_;
    if(bkgModel_) delete bkgModel_;
  }
 
  void LifetimeModel::loadParameters(TDirectory& dir) {
    if(isTemplate_) {
      RooHistPdf *prMod, *nprMod, *bkgMod;

      dir.GetObject("promptCtHist",templates_["promptCtHist"]);
      dir.GetObject("nonPromptCtHist",templates_["nonPromptCtHist"]);
      dir.GetObject("bkgCtHist",templates_["bkgCtHist"]);      

      dir.GetObject("promptCtPdf",prMod);
      if(prMod) promptModel_ = prMod;
      dir.GetObject("nonPromptCtPdf",nprMod);
      if(nprMod) nonPromptModel_ = nprMod;
      dir.GetObject("bkgCtPdf",bkgMod);  
      if(bkgMod) bkgModel_ = bkgMod;

    } else {
      //Gauss 1
      dir.GetObject("ctMeanPrompt1",vars_["ctMeanPrompt1"]);
      if(!vars_["ctMeanPrompt1"]) vars_["ctMeanPrompt1"] = new RooRealVar("ctMeanPrompt1",
									  "Prompt C-Tau Mean for Gaussian 1",
									  0,-0.15,0.15);
      dir.GetObject("ctMeanPrompt2",vars_["ctMeanPrompt2"]);
      if(!vars_["ctMeanPrompt2"]) vars_["ctMeanPrompt2"] = new RooRealVar("ctMeanPrompt2",
									  "Prompt C-Tau Mean for Gaussian 2",
									  -1.22534e-3,-0.1,0.1);
      dir.GetObject("ctMeanPrompt3",vars_["ctMeanPrompt3"]);
      if(!vars_["ctMeanPrompt3"]) vars_["ctMeanPrompt3"] = new RooRealVar("ctMeanPrompt3",
									  "Prompt C-Tau Mean for Gaussian 3",
									  1.308879e-2,-0.1,0.1);
      
      dir.GetObject("ctSigmaPrompt1",vars_["ctSigmaPrompt1"]);
      if(!vars_["ctSigmaPrompt1"]) vars_["ctSigmaPrompt1"] = new RooRealVar("ctSigmaPrompt1",
									    "Prompt C-Tau Sigma for Gaussian 1",
									    .876883e-1,0.02,1);
      
      dir.GetObject("ctCoefPrompt1",vars_["ctCoefPrompt1"]);
      if(!vars_["ctCoefPrompt1"]) vars_["ctCoedPrompt1"] = new RooRealVar("ctCoefPrompt1",
									  "Prompt C-Tau Normalization for Gaussian 1",
									  .7992e-2,0,1);
      
      //Gauss 2
      dir.GetObject("ctSigmaPrompt2",vars_["ctSigmaPrompt2"]);
      if(!vars_["ctSigmaPrompt2"]) vars_["ctSigmaPrompt2"] = new RooRealVar("ctSigmaPrompt2",
									    "Prompt C-Tau Sigma for Gaussian 2",
									    1.42654e-1,0.02,1);
      
      dir.GetObject("ctCoefPrompt2",vars_["ctCoefPrompt2"]);
      if(!vars_["ctCoefPrompt2"]) vars_["ctCoedPrompt2"] = new RooRealVar("ctCoefPrompt2",
									  "Prompt C-Tau Normalization for Gaussian 2",
									  9.999e-1,0,1);
      
      //Gauss 3

      
      dir.GetObject("ctSigmaPrompt3",vars_["ctSigmaPrompt3"]);
      if(!vars_["ctSigmaPrompt3"]) vars_["ctSigmaPrompt3"] = new RooRealVar("ctSigmaPrompt3",
									    "Prompt C-Tau Sigma for Gaussian 3",
									    2.8883e-1,0.02,1);
      
      dir.GetObject("ctCoefPrompt3",vars_["ctCoefPrompt3"]);
      if(!vars_["ctCoefPrompt3"]) vars_["ctCoefPrompt3"] = new RooRealVar("ctCoefPrompt3",
									  "Prompt C-Tau Normalization for Gaussian 3",
									  5.91949e-1,0,1);
      
      dir.GetObject("ctCoefPromptTot",vars_["ctCoefPromptTot"]);
      if(!vars_["ctCoefPromptTot"]) vars_["ctCoefPromptTot"] = new RooRealVar("ctCoefPromptTot",
									      "Coefficient for Prompt Contribution",
									      .1,0,1);
      
      //non prompt
      //SS
      dir.GetObject("ssTauNonPrompt",vars_["ssTauNonPrompt"]);
      if(!vars_["ssTauNonPrompt"]) vars_["ssTauNonPrompt"] = new RooRealVar("ssTauNonPrompt",
									    "Single Sided Non-Prompt Decay parameter",
									    1,0.01,5);
      
      dir.GetObject("ssCoefNonPrompt",vars_["ssCoefNonPrompt"]);
      if(!vars_["ssCoefNonPrompt"]) vars_["ssCoefNonPrompt"] = new RooRealVar("ssCoefNonPrompt",
									      "Single Sided Non-Prompt Coefficient",
									      .5,0,1);
      
      //Flipped
      dir.GetObject("fTauNonPrompt",vars_["fTauNonPrompt"]);
      if(!vars_["fTauNonPrompt"]) vars_["fTauNonPrompt"] = new RooRealVar("fTauNonPrompt",
									  "Flipped Non-Prompt Decay parameter",
									  1,0.01,5);

      dir.GetObject("fCoefNonPrompt",vars_["fCoefNonPrompt"]);
            if(!vars_["fCoefNonPrompt"]) vars_["fCoefNonPrompt"] = new RooRealVar("fCoefNonPrompt",
      									      "Flipped Non-Prompt Coefficient",
      									      .5,0,1);

      //DS
      dir.GetObject("dsTauNonPrompt",vars_["dsTauNonPrompt"]);
      if(!vars_["dsTauNonPrompt"]) vars_["dsTauNonPrompt"] = new RooRealVar("dsTauNonPrompt",
									    "Double Sided Non-Prompt Decay parameter",
									    1,0.01,5);
      
      //bkg
      //SS
      dir.GetObject("ssTauBkg",vars_["ssTauBkg"]);
      if(!vars_["ssTauBkg"]) vars_["ssTauBkg"] = new RooRealVar("ssTauBkg",
								"Single Sided Non-Prompt Decay parameter",
								1,0.01,5);
      
      dir.GetObject("ssCoefBkg",vars_["ssCoefBkg"]);
      if(!vars_["ssCoefBkg"]) vars_["ssCoefBkg"] = new RooRealVar("ssCoefBkg",
								  "Single Sided Non-Prompt Coefficient",
								  .1,0,1);
      
      //Flipped
      dir.GetObject("fTauBkg",vars_["fTauBkg"]);
      if(!vars_["fTauBkg"]) vars_["fTauBkg"] = new RooRealVar("fTauBkg",
							      "Flipped Non-Prompt Decay parameter",
							      1,0.01,5);
      
      dir.GetObject("fCoefBkg",vars_["fCoefBkg"]);
      if(!vars_["fCoefBkg"]) vars_["fCoefBkg"] = new RooRealVar("fCoefBkg",
								"Flipped Background Coefficient",
								0.1,0,1);
      
      //DS
      dir.GetObject("dsTauBkg",vars_["dsTauBkg"]);
      if(!vars_["dsTauBkg"]) vars_["dsTauBkg"] = new RooRealVar("dsTauBkg",
								"Double Sided Non-Prompt Decay parameter",
								1,0.01,5);

      dir.GetObject("dsCoefBkg",vars_["dsCoefBkg"]);
      if(!vars_["dsCoefBkg"]) vars_["dsCoefBkg"] = new RooRealVar("dsCoefBkg",
     								"Doublesided Background Coefficient",
     								0.1,0,1);


    }
    varsLoaded = true;
  }

  void LifetimeModel::saveParameters(TDirectory& dir) {
    if(isTemplate_) {
      if(promptModel_) dir.Add(promptModel_->Clone(),true);
      if(nonPromptModel_) dir.Add(nonPromptModel_->Clone(),true);
      if(bkgModel_) dir.Add(bkgModel_->Clone(),true);      
    } else {
      for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	  i != vars_.end();
	  ++i) {
	if(i->second) 
	  dir.Add(i->second->Clone(),true);
      }
    }
  }

  void LifetimeModel::fix(const std::string& s) {
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

  void LifetimeModel::unfix(const std::string& s) {
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

  void LifetimeModel::initModels(RooRealVar& var,
				 const TH1F* tempPr,
				 const TH1F* tempNPr,
				 const TH1F* tempBkg,
				 const bool& doSig, 
				 const bool& doNP, 
				 const bool& doBKG) {
    if(doSig) initPromptModel(var,tempPr);
    if(doNP)  initNonPromptModel(var,tempNPr);
    if(doBKG) initBackgroundModel(var,tempBkg);
  }

  void LifetimeModel::initModels(RooRealVar& var,
				 RooRealVar& varErr,
				 const TH1F* tempPr,
				 const TH1F* tempNPr,
				 const TH1F* tempBkg,
				 const bool& doSig, 
				 const bool& doNP, 
				 const bool& doBKG) {
    if(doSig) initPromptModel(var,varErr,tempPr);
    if(doNP)  initNonPromptModel(var,tempNPr);
    if(doBKG) initBackgroundModel(var,tempBkg);
  }

  // initialize the lifetime PDF for the Prompt component
  void LifetimeModel::initPromptModel(RooRealVar& var, const TH1F* temp) {
    if(isTemplate_) {
      if(!varsLoaded) {
	templates_["promptCtHist"] = new RooDataHist("promptCtHist",
						    "Prompt C-Tau Template Histogram",
						    RooArgSet(var),
						    temp);
	
	promptModel_ = new RooHistPdf("promptCtPdf",
				      "Template PDF for Prompt C-Tau",
				      RooArgSet(var),
				      *templates_["promptCtHist"]);
      }
    } else {
      // Gauss 1
      RooRealVar *ctMeanPrompt1;
      RooRealVar *ctMeanPrompt2;
      RooRealVar *ctMeanPrompt3;
      RooRealVar *ctSigmaPrompt1;
      RooRealVar *ctCoefPrompt1;
      // Gauss 2 (same mean as Gauss 1)
      RooRealVar *ctSigmaPrompt2;
      RooRealVar *ctCoefPrompt2;
      //Gauss 3
      RooRealVar *ctSigmaPrompt3;
      RooRealVar *ctCoefPrompt3;
      RooRealVar *ctCoefPromptTot;
      if(varsLoaded) {
	//Gauss 1
	ctMeanPrompt1 = vars_["ctMeanPrompt1"] ;
	ctSigmaPrompt1 = vars_["ctSigmaPrompt1"];
	ctCoefPrompt1 = vars_["ctCoefPrompt1"];      
	//Gauss 2
	ctSigmaPrompt2 = vars_["ctSigmaPrompt2"];
	ctCoefPrompt2 = vars_["ctCoefPrompt2"];
	ctMeanPrompt2 = vars_["ctMeanPrompt2"] ;
	//Gauss 3
	ctSigmaPrompt3 = vars_["ctSigmaPrompt3"];
	ctCoefPrompt3 = vars_["ctCoefPrompt3"];   
	ctMeanPrompt3 = vars_["ctMeanPrompt3"] ;

	ctCoefPromptTot = vars_["ctCoefPromptTot"];
      } else {
	//Gauss 1
	ctMeanPrompt1 = new RooRealVar("ctMeanPrompt1",
				       "Prompt C-Tau Mean for Gaussian 1",
				       0,-.15,.15);
	vars_["ctMeanPrompt1"] = ctMeanPrompt1;
	ctMeanPrompt2 = new RooRealVar("ctMeanPrompt2",
				       "Prompt C-Tau Mean for Gaussian 2",
				       -1.25534e-3,-.1,.1);
	vars_["ctMeanPrompt2"] = ctMeanPrompt2;
	ctMeanPrompt3 = new RooRealVar("ctMeanPrompt3",
				       "Prompt C-Tau Mean for Gaussian 3",
				       1.30879e-2,-.1,.1);
	vars_["ctMeanPrompt3"] = ctMeanPrompt3;
	ctSigmaPrompt1 = new RooRealVar("ctSigmaPrompt1",
					"Prompt C-Tau Sigma for Gaussian 1",
					1.05,0.01,1);
	vars_["ctSigmaPrompt1"] = ctSigmaPrompt1;
	ctCoefPrompt1 = new RooRealVar("ctCoefPrompt1",
				       "Prompt C-Tau Normalization for Gaussian 1",
				       7.92205e-2,0,1);
	vars_["ctCoefPrompt1"] = ctCoefPrompt1;
	//Gauss 2
	ctSigmaPrompt2 = new RooRealVar("ctSigmaPrompt2",
					"Prompt C-Tau Sigma for Gaussian 2",
					1.42654e-1,0.01,1);
	vars_["ctSigmaPrompt2"] = ctSigmaPrompt2;
	ctCoefPrompt2 = new RooRealVar("ctCoefPrompt2",
				       "Prompt C-Tau Normalization for Gaussian 2",
				       9.9999e-1,0,1);
	vars_["ctCoefPrompt2"] = ctCoefPrompt2;	
	//Gauss 3
	ctSigmaPrompt3 = new RooRealVar("ctSigmaPrompt3",
					"Prompt C-Tau Sigma for Gaussian 3",
					2.8883e-1,0.01,1);
	vars_["ctSigmaPrompt3"] = ctSigmaPrompt3;
	ctCoefPrompt3 = new RooRealVar("ctCoefPrompt3",
				       "Prompt C-Tau Normalization for Gaussian 3",
				       5.91949e-1,0,1);
	vars_["ctCoefPrompt3"] = ctCoefPrompt3;
	
	ctCoefPromptTot = new RooRealVar("ctCoefPromptTot",
					 "Coefficient for Prompt Contribution",
					 .1,0,1);
	vars_["ctCoefPromptTot"] = ctCoefPromptTot;      
      }
      
      ctGaussPrompt1 = new RooGaussModel("ctGaussPrompt1",
					 "Prompt C-Tau Gaussian Model 1",
					 var,
					 *ctMeanPrompt1,*ctSigmaPrompt1);
      gmodels_.push_back(ctGaussPrompt1);
      //Gauss 2
      ctGaussPrompt2 = new RooGaussModel("ctGaussPrompt2",
					 "Prompt C-Tau Gaussian Model 2",
					 var,
					 *ctMeanPrompt2,*ctSigmaPrompt2);
      gmodels_.push_back(ctGaussPrompt2);    
      // Gauss 3
      ctGaussPrompt3 = new RooGaussModel("ctGaussPrompt3",
					 "Prompt C-Tau Gaussian Model 3",
					 var,
					 *ctMeanPrompt3,*ctSigmaPrompt3);
      gmodels_.push_back(ctGaussPrompt3);
      // Create composite model of the prompt lifetime
      RooAddModel *ctPromptModel = new RooAddModel("ctPromptModel",
						   "C-Tau Composite Model for Prompt",
						   RooArgList(*ctGaussPrompt1,*ctGaussPrompt2,*ctGaussPrompt3),
						   RooArgList(*ctCoefPrompt1,*ctCoefPrompt2,*ctCoefPrompt3));
      addmodels_["ResModel"] = ctPromptModel;
      
      
      promptModel_ = new RooAddPdf("promptCtPdf",
				   "PDF for Prompt C-Tau",
				   RooArgList(*ctGaussPrompt1,*ctGaussPrompt2,*ctGaussPrompt3),
				   RooArgList(*ctCoefPrompt1,*ctCoefPrompt2,*ctCoefPrompt3));
    }
  }

  void LifetimeModel::initPromptModel(RooRealVar& var, RooRealVar& varErr, const TH1F* temp) {
    RooRealVar *ctMeanPrompt1;
    RooRealVar *ctSigmaPrompt1;

    if(varsLoaded) {
	//Gauss 1
	ctMeanPrompt1 = vars_["ctMeanPrompt1"] ;
	ctSigmaPrompt1 = vars_["ctSigmaPrompt1"];	
      } else {
	//Gauss 1
	ctMeanPrompt1 = new RooRealVar("ctMeanPrompt1",
				       "Prompt C-Tau Mean for Gaussian 1",
				      0,-.15,.15);
	vars_["ctMeanPrompt1"] = ctMeanPrompt1;	
	ctSigmaPrompt1 = new RooRealVar("ctSigmaPrompt1",
					"Prompt C-Tau Sigma for Gaussian 1",
					1.05,0.001,5);
	vars_["ctSigmaPrompt1"] = ctSigmaPrompt1;	
      }
    ctGaussPrompt1 = new RooGaussModel("ctGaussPrompt1",
				       "Prompt C-Tau Gaussian Model 1",
				       var,
				       *ctMeanPrompt1,*ctSigmaPrompt1, 
				       RooRealConstant::value(1), varErr);
    gmodels_.push_back(ctGaussPrompt1);

    RooAddModel *ctPromptModel = new RooAddModel("ctPromptModel",
						 "C-Tau Composite Model for Prompt",
						 RooArgList(*ctGaussPrompt1),
						 RooArgList());

    addmodels_["ResModel"] = ctPromptModel;
        
    promptModel_ = new RooAddPdf("promptCtPdf",
				 "PDF for Prompt C-Tau",
				 RooArgList(*ctGaussPrompt1),
				 RooArgList());
  }

  void LifetimeModel::initNonPromptModel(RooRealVar& var, const TH1F* temp) {
    
    if(!promptModel_) {
      std::cout << "Cannot call initNonPromptModel() before calling initPromptModel()" << std::endl;
    }

    if(isTemplate_) {
      if(!varsLoaded) {
	templates_["nonPromptCtHist"] = new RooDataHist("nonPromptCtHist",
						       "Non-Prompt C-Tau Template Histogram",
						       RooArgSet(var),
						       temp);
	
	nonPromptModel_ = new RooHistPdf("nonPromptCtPdf",
					 "Template PDF for Non-Prompt C-Tau",
					 RooArgSet(var),
					 *templates_["nonPromptCtHist"]);
      }
    } else {
      
      // Single Sided decay model
      RooRealVar *ssTau;
      RooRealVar *ssCoef;
      //Flipped decay model
      RooRealVar *fTau;
      RooRealVar *fCoef;
      //Double Sided decay model
      RooRealVar *dsTau;
      
      if(varsLoaded) {
	//SS
	ssTau  = vars_["ssTauNonPrompt"] ;
	ssCoef = vars_["ssCoefNonPrompt"];       
	//Flipped
	fTau = vars_["fTauNonPrompt"];
	fCoef = vars_["fCoefNonPrompt"];
	//DS
	dsTau = vars_["dsTauNonPrompt"];
      } else {
	//SS
	ssTau = new RooRealVar("ssTauNonPrompt",
			       "Single Sided Non-Prompt Decay parameter",
			       1,0.001,10);
	vars_["ssTauNonPrompt"] = ssTau;
	ssCoef = new RooRealVar("ssCoefNonPrompt",
				"Single Sided Non-Prompt Coefficient",
				.5,0,1);
	vars_["ssCoefNonPrompt"] = ssCoef;
	//Flipped
	fTau = new RooRealVar("fTauNonPrompt",
			      "Flipped Non-Prompt Decay parameter",
			      1,0.001,10);
	vars_["fTauNonPrompt"] = fTau;
	fCoef = new RooRealVar("fCoefNonPrompt",
					"Flipped Non-Prompt Coefficient",
					.5,0,1);
		vars_["fCoefNonPrompt"] = fCoef;

	//DS
	dsTau = new RooRealVar("dsTauNonPrompt",
			       "Double Sided Non-Prompt Decay parameter",
			       1,0.01,5);
	vars_["dsTauNonPrompt"] = dsTau;
      }
      
      ssDecayNP = new RooDecay("ssDecayNonPrompt",
				       "Single Sided Non-Prompt Decay Model",
				       var,*ssTau,*addmodels_["ResModel"],
				       RooDecay::SingleSided);
      dmodels_.push_back(ssDecayNP);
      //Flipped decay model
      fDecayNP = new RooDecay("fDecayNonPrompt",
				      "Flipped Non-Prompt Decay Model",
				      var,*fTau,*addmodels_["ResModel"],
				      RooDecay::Flipped);
      dmodels_.push_back(fDecayNP);
      //Double Sided decay model
      dsDecayNP = new RooDecay("dsDecayNonPrompt",
				       "Double Sided Non-Prompt Decay Model",
				       var,*dsTau,*addmodels_["ResModel"],
				       RooDecay::DoubleSided);
      dmodels_.push_back(dsDecayNP);
      
      nonPromptModel_ = new RooAddPdf("nonPromptCtPdf",
				      "Non-Prompt Decay PDF",
				      RooArgSet(*ssDecayNP),//,*fDecayNP),
				      RooArgSet());//*ssCoef));//,*fCoef));
    }
  }
  
  void LifetimeModel::initBackgroundModel(RooRealVar& var, const TH1F* temp) {

    if(!promptModel_) {
      std::cout << "Cannot call initBackgroundModel() before calling initPromptModel()" << std::endl;
    }

    if(isTemplate_) {
      if(!varsLoaded) {
	templates_["bkgCtHist"] = new RooDataHist("bkgCtHist",
						 "Background C-Tau Template Histogram",
						 RooArgSet(var),
						 temp);
	
	bkgModel_ = new RooHistPdf("bkgCtPdf",
				   "Template PDF for Non-Prompt C-Tau",
				   RooArgSet(var),
				   *templates_["bkgCtHist"]);
      }
    } else {
      
      // Single Sided decay model
      RooRealVar *ssTau;
      RooRealVar *ssCoef;
      //Flipped decay model
      RooRealVar *fTau;
      RooRealVar *fCoef;
      //Double Sided decay model
      RooRealVar *dsTau;
      RooRealVar *dsCoef;

      
      if(varsLoaded) {
	//SS
	ssTau  = vars_["ssTauBkg"] ;
	ssCoef = vars_["ssCoefBkg"];       
	//Flipped
	fTau = vars_["fTauBkg"];   
	fCoef = vars_["fCoefBkg"];
	//DS
	dsTau = vars_["dsTauBkg"];
	dsCoef = vars_["dsCoefBkg"];

      } else {
	//SS
	ssTau = new RooRealVar("ssTauBkg",
			       "Single Sided Non-Prompt Decay parameter",
			       1,0.001,10);
	vars_["ssTauBkg"] = ssTau;
	ssCoef = new RooRealVar("ssCoefBkg",
				"Single Sided Non-Prompt Coefficient",
				.35,0,1);
	vars_["ssCoefBkg"] = ssCoef;
	//Flipped
	fTau = new RooRealVar("fTauBkg",
			      "Flipped Non-Prompt Decay parameter",
			      1,0.001,10);
	vars_["fTauBkg"] = fTau;
	fCoef = new RooRealVar("fCoefBkg",
			       "Flipped Background Coefficient",
			       0.35,0,1);
	vars_["fCoefBkg"] = fCoef;
	//DS
	dsTau = new RooRealVar("dsTauBkg",
			       "Double Sided Non-Prompt Decay parameter",
			       1,0.001,10);
	vars_["dsTauBkg"] = dsTau;
	dsCoef = new RooRealVar("dsCoefBkg",
			       "Doublesided Background Coefficient",
			       0.35,0,1);
	vars_["dsCoefBkg"] = dsCoef;
      }
      
      //Single Sided decay model
      ssDecayBK = new RooDecay("ssDecayBkg",
				       "Single Sided Background Decay Model",
				       var,*ssTau,*addmodels_["ResModel"],
				       RooDecay::SingleSided);
      dmodels_.push_back(ssDecayBK);
      //Flipped decay model
      fDecayBK = new RooDecay("fDecayBkg",
				      "Flipped Non-Prompt Decay Model",
				      var,*fTau,*addmodels_["ResModel"],
				      RooDecay::Flipped);
      dmodels_.push_back(fDecayBK);
      //Double Sided decay model    
      dsDecayBK = new RooDecay("dsDecayBkg",
				       "Double Sided Background Decay Model",
				       var,*dsTau,*addmodels_["ResModel"],
				       RooDecay::DoubleSided);
      dmodels_.push_back(dsDecayBK);
      
      bkgModel_ = new RooAddPdf("bkgCtPdf",
				"Non-Prompt Decay PDF",
				RooArgSet(*ssDecayBK,*fDecayBK,*dsDecayBK),
				RooArgSet(*ssCoef,*fCoef));
    }
  }
}
