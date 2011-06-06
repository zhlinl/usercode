#include "CompositeModelBuilder.h"

#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooAbsReal.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TDirectory.h"

#include <cmath>

namespace JPsiPolarization {
  CompositeModelBuilder::CompositeModelBuilder(const std::string& frame,
					       const std::string& range) {
    fName_ = frame;
    range_ = range;

    useP = true;
    useNP = true;
    useBkg = true;
    useMass = true;
    useLifetime = true;
    usePol = true;
    useAcc = true;
    
    mass_ = new MassModel();
    lifet_ = new LifetimeModel();
    promptPol_ = new PolarizationModel("prompt",fName_);
    nonPromptPol_ = new PolarizationModel("nonPrompt",fName_);
    acc_ = new AcceptanceMaps(fName_);

    promptAccHist_ = NULL;
    nonPromptAccHist_ = NULL;
    bkgShapeHist_ = NULL;

    promptLifetimeTemplate_ = NULL;
    nonPromptLifetimeTemplate_ = NULL;
    bkgLifetimeTemplate_ = NULL;

    promptModel_ = NULL;
    nonPromptModel_ = NULL;
    bkgModel_ = NULL;
    compositeModel_ = NULL;
  }

  CompositeModelBuilder::~CompositeModelBuilder() {
    
    for(std::map<std::string,RooRealVar*>::iterator i = vars_.begin();
	i != vars_.end();
	++i)
      delete i->second;

    for(std::vector<RooAbsPdf*>::iterator i = pdfs_.begin();
	i != pdfs_.end();
	++i)
      delete *i;

    delete mass_;
    delete lifet_;
    delete promptPol_;
    delete nonPromptPol_;
    delete acc_;
    delete promptModel_;
    delete nonPromptModel_;
    delete bkgModel_;
    delete compositeModel_;
  }
  
  void CompositeModelBuilder::loadParameters(TDirectory& dir) {
    TDirectory *mass, *lifet, *promptPol, *nonPromptPol, *acc;

    if(useMass && (mass = dir.GetDirectory("MassModel")))
      mass_->loadParameters(*mass);
    if(useLifetime && (lifet = dir.GetDirectory("LifetimeModel")))
      lifet_->loadParameters(*lifet);
    if(usePol && (promptPol = dir.GetDirectory("PromptPolarizationModel")))
      promptPol_->loadParameters(*promptPol);
    if(usePol && (nonPromptPol = dir.GetDirectory("NonPromptPolarizationModel")))
      nonPromptPol_->loadParameters(*nonPromptPol);
    if(acc_ && (acc = dir.GetDirectory("AcceptanceMaps")))
      acc_->loadParameters(*acc);

    dir.GetObject("nPrompt",vars_["nPrompt"]);
    if(!vars_["nPrompt"]) vars_["nPrompt"] = new RooRealVar("nPrompt","Fitted Number of Prompt J/Psis",5000,0,1000000);
    
    dir.GetObject("nNonPrompt",vars_["nNonPrompt"]);
    if(!vars_["nNonPrompt"]) vars_["nNonPrompt"] = new RooRealVar("nNonPrompt","Fitted Number of Non-Prompt J/Psis",5000,0,1000000);
    
    dir.GetObject("nBackground",vars_["nBackground"]);
    if(!vars_["nBackground"]) vars_["nBackground"] = new RooRealVar("nBackground","Fitted Number of Background Events",500,0,1000000);

  }

  void CompositeModelBuilder::saveParameters(TDirectory& dir) {
    TDirectory *mass, *lifet, *promptPol, *nonPromptPol, *acc;
    if(useMass) {
      if(!(mass = dir.GetDirectory("MassModel")))
	mass = dir.mkdir("MassModel");
      mass_->saveParameters(*mass);
    }
    if(useLifetime) {
      if(!(lifet = dir.GetDirectory("LifetimeModel")))
	lifet = dir.mkdir("LifetimeModel");
      lifet_->saveParameters(*lifet);
    }
    if(usePol && promptPol_->model()) {
      if(!(promptPol = dir.GetDirectory("PromptPolarizationModel")))
	promptPol = dir.mkdir("PromptPolarizationModel");
      promptPol_->saveParameters(*promptPol);
    }
    if(usePol && nonPromptPol_->model()) {
      if(!(nonPromptPol = dir.GetDirectory("NonPromptPolarizationModel")))
	nonPromptPol = dir.mkdir("NonPromptPolarizationModel");
      nonPromptPol_->saveParameters(*nonPromptPol);
    }
    if(acc_) {
      if(!(acc = dir.GetDirectory("AcceptanceMaps")))
	acc = dir.mkdir("AcceptanceMaps");
      acc_->saveParameters(*acc);
    }    
  }

  void CompositeModelBuilder::fix(const std::string& pdfName) {
    if(pdfName.find("mass") != std::string::npos) mass_->fix();
    if(pdfName.find("lifetime") != std::string::npos) lifet_->fix();
    if(pdfName.find("nBackground") != std::string::npos) vars_["nBackground"]->setConstant(true);
    if(pdfName.find("nPrompt") != std::string::npos) vars_["nPrompt"]->setConstant(true);
    if(pdfName.find("nNonPrompt") != std::string::npos) vars_["nNonPrompt"]->setConstant(true);
    
  }

  void CompositeModelBuilder::unfix(const std::string& pdfName) {
    if(pdfName.find("mass") != std::string::npos) mass_->unfix();
    if(pdfName.find("lifetime") != std::string::npos) lifet_->unfix();
    if(pdfName.find("nBackground") != std::string::npos) vars_["nBackground"]->setConstant(false);
    if(pdfName.find("nPrompt") != std::string::npos) vars_["nPrompt"]->setConstant(false);
    if(pdfName.find("nNonPrompt") != std::string::npos) vars_["nNonPrompt"]->setConstant(false);
    
  }

  void CompositeModelBuilder::saveParameter(const std::string& name,TDirectory& dir) {       
    if(name.find("nBackground") != std::string::npos) dir.Add(vars_["nBackground"]->Clone(),true);
    if(name.find("nPrompt") != std::string::npos)     dir.Add(vars_["nPrompt"]->Clone(),true);
    if(name.find("nNonPrompt") != std::string::npos)  dir.Add(vars_["nNonPrompt"]->Clone(),true);

    dir.Write();
  }

  void CompositeModelBuilder::initModel(RooRealVar& mass,
					RooRealVar& ct,
					RooRealVar& costh,
					RooRealVar& phi) {
    if(useMass && !(mass_->prompt() || mass_->nonPrompt() || mass_->background())) {
      std::cout << "Initializing Mass Model" << std::endl;
      mass_->initModels(mass,useP,useNP,useBkg);
    }
    if(useLifetime && !(lifet_->prompt() || lifet_->nonPrompt() || lifet_->background())) {
      std::cout << "Initializing Lifetime Model" << std::endl;
      lifet_->initModels(ct,promptLifetimeTemplate_,
			 promptLifetimeTemplate_,
			 promptLifetimeTemplate_,
			 useP,useNP,useBkg);
    }
    if(usePol)      {
      if(useP && !promptPol_->model()) {
	std::cout << "Initializing Prompt Polarization Model" << std::endl;
	promptPol_->initModel(costh,phi);
      }
      if(useNP && !nonPromptPol_->model()) {
	std::cout << "Initializing Non-Prompt Polarization Model" << std::endl;
	nonPromptPol_->initModel(costh,phi);
      }
    }
    if(useAcc && !(acc_->prompt() || acc_->nonPrompt() || acc_->background())) {
      std::cout << "Initializing Acceptance Histograms" << std::endl;
      if(useP)   acc_->setPromptHist(promptAccHist_);
      if(useNP)  acc_->setNonPromptHist(nonPromptAccHist_);
      if(useBkg) acc_->setBackgroundHist(bkgShapeHist_);
      acc_->initAcceptanceMaps(costh,phi);
    }
    buildModel();
  }

  void CompositeModelBuilder::initModel(RooRealVar& mass,
					RooRealVar& ct,
					RooRealVar& cterr,
					RooRealVar& costh,
					RooRealVar& phi) { 
    if(useMass && !(mass_->prompt() || mass_->nonPrompt() || mass_->background())) {
      std::cout << "Initializing Mass Model" << std::endl;
      mass_->initModels(mass,useP,useNP,useBkg);
    }
    if(useLifetime && !(lifet_->prompt() || lifet_->nonPrompt() || lifet_->background())) {
      std::cout << "Initializing Lifetime Model" << std::endl;
      lifet_->initModels(ct,cterr,
			 promptLifetimeTemplate_,
			 promptLifetimeTemplate_,
			 promptLifetimeTemplate_,
			 useP,useNP,useBkg);
    }
    if(usePol)      {
      if(useP && !promptPol_->model()) {
	std::cout << "Initializing Prompt Polarization Model" << std::endl;
	promptPol_->initModel(costh,phi);
      }
      if(useNP && !nonPromptPol_->model()) {
	std::cout << "Initializing Non-Prompt Polarization Model" << std::endl;
	nonPromptPol_->initModel(costh,phi);
      }
    }
    if(useAcc && !(acc_->prompt() || acc_->nonPrompt() || acc_->background())) {
      std::cout << "Initializing Acceptance Histograms" << std::endl;
      if(useP)   acc_->setPromptHist(promptAccHist_);
      if(useNP)  acc_->setNonPromptHist(nonPromptAccHist_);
      if(useBkg) acc_->setBackgroundHist(bkgShapeHist_);
      acc_->initAcceptanceMaps(costh,phi);
    }
    buildModel(ct);
  }

  void CompositeModelBuilder::buildModel() {
    
    RooRealVar* prompt;
    RooRealVar* nonprompt;
    RooRealVar* bkg;
    if(!vars_.size()) {
	prompt = new RooRealVar("nPrompt","Fitted Number of Prompt J/Psis",5000,0,1000000);
	nonprompt = new RooRealVar("nNonPrompt","Fitted Number of Non-Prompt J/Psis",5000,0,1000000);
	bkg = new RooRealVar("nBackground","Fitted Number of Background Events",500,0,1000000);
	vars_["nPrompt"] = prompt;
	vars_["nNonPrompt"] = nonprompt;
	vars_["nBackground"] = bkg;
    } else {
      prompt = vars_["nPrompt"];
      nonprompt = vars_["nNonPrompt"];
      bkg = vars_["nBackground"];
      
      for(std::vector<RooAbsPdf*>::iterator i = pdfs_.begin();
	  i != pdfs_.end();
	  ++i) {
	delete *i;
      }
      
      pdfs_.clear();

      delete promptModel_;
      delete nonPromptModel_;
      delete bkgModel_;
      delete compositeModel_;
    }

    RooArgList promptpdfpieces, nonpromptpdfpieces, bkgpdfpieces;
    RooArgList promptPolArgs, nonPromptPolArgs;
    RooProdPdf *promptPolProd, *nonPromptPolProd;
    RooProdPdf *promptbasepdf, *nonpromptbasepdf, *bkgbasepdf;

    if(useP) {
      if(mass_->prompt()) promptpdfpieces.add(*mass_->prompt());
      if(lifet_->prompt()) promptpdfpieces.add(*lifet_->prompt());
      if(promptPol_->model()) promptPolArgs.add(*promptPol_->model());
      if(acc_->prompt()) promptPolArgs.add(*acc_->prompt());
      promptPolProd = new RooProdPdf("promptPolProd","Product of Prompt Polarization PDF and Acceptance Map",promptPolArgs);
      pdfs_.push_back(promptPolProd);
      promptpdfpieces.add(*promptPolProd);
    } 

    if(useNP) {
      if(mass_->nonPrompt()) nonpromptpdfpieces.add(*mass_->nonPrompt());
      if(lifet_->nonPrompt()) nonpromptpdfpieces.add(*lifet_->nonPrompt());
      if(nonPromptPol_->model()) nonPromptPolArgs.add(*nonPromptPol_->model());
      if(acc_->nonPrompt()) nonPromptPolArgs.add(*acc_->nonPrompt());
      nonPromptPolProd = new RooProdPdf("nonPromptPolProd","Product of NonPrompt Polarization PDF and Acceptance Map",nonPromptPolArgs);
      pdfs_.push_back(nonPromptPolProd);
      nonpromptpdfpieces.add(*nonPromptPolProd);
    }

    if(useBkg) {
      if(mass_->background()) bkgpdfpieces.add(*mass_->background());
      if(lifet_->background()) bkgpdfpieces.add(*lifet_->background());
      if(acc_->background()) bkgpdfpieces.add(*acc_->background());
    }

    promptbasepdf = new RooProdPdf("prodPromptPdf","Product of all Prompt PDFs",
				   promptpdfpieces);
    pdfs_.push_back(promptbasepdf);
    nonpromptbasepdf = new RooProdPdf("prodNonPromptPdf","Product of all Non-Prompt PDFs",
				      nonpromptpdfpieces);
    pdfs_.push_back(nonpromptbasepdf);
    bkgbasepdf = new RooProdPdf("prodBkgPdf","Product of all Background PDFs",
				bkgpdfpieces);
    pdfs_.push_back(bkgbasepdf);
    
    //create the extended PDFs for prompt, nonprompt and bkg
    if(!range_.size()) {
      promptModel_ = new RooExtendPdf("promptModel","The Prompt Contribution to the PDF",
				      *promptbasepdf,*prompt);
      nonPromptModel_ = new RooExtendPdf("nonPromptModel","The Non-Prompt Contribution to the PDF",
				       *nonpromptbasepdf,*nonprompt);
      bkgModel_ = new RooExtendPdf("bkgModel","The Background Contribution to the PDF",
				   *bkgbasepdf,*bkg);
    } else {
      promptModel_ = new RooExtendPdf("promptModel","The Prompt Contribution to the PDF",
				      *promptbasepdf,*prompt,range_.c_str());
      nonPromptModel_ = new RooExtendPdf("nonPromptModel","The Non-Prompt Contribution to the PDF",
					 *nonpromptbasepdf,*nonprompt,range_.c_str());
      bkgModel_ = new RooExtendPdf("bkgModel","The Background Contribution to the PDF",
				   *bkgbasepdf,*bkg,range_.c_str());
    }
    
    RooArgList pdfPieces;
    if(useP) pdfPieces.add(*promptModel_);
    if(useNP) pdfPieces.add(*nonPromptModel_);
    if(useBkg) pdfPieces.add(*bkgModel_);

    //create the full PDF for fitting
    compositeModel_ = new RooAddPdf("compositePDF","The Composite Model for Fitting J/Psi Mass/Lifetime/Polarization",
				    pdfPieces);
    
  }

  void CompositeModelBuilder::buildModel(RooRealVar& ct) {
    
    RooRealVar* prompt;
    RooRealVar* nonprompt;
    RooRealVar* bkg;
    if(!vars_.size()) {
	prompt = new RooRealVar("nPrompt","Fitted Number of Prompt J/Psis",5000,0,1000000);
	nonprompt = new RooRealVar("nNonPrompt","Fitted Number of Non-Prompt J/Psis",5000,0,1000000);
	bkg = new RooRealVar("nBackground","Fitted Number of Background Events",500,0,1000000);
	vars_["nPrompt"] = prompt;
	vars_["nNonPrompt"] = nonprompt;
	vars_["nBackground"] = bkg;
    } else {
      prompt = vars_["nPrompt"];
      nonprompt = vars_["nNonPrompt"];
      bkg = vars_["nBackground"];
      
      for(std::vector<RooAbsPdf*>::iterator i = pdfs_.begin();
	  i != pdfs_.end();
	  ++i) {
	delete *i;
      }
      
      pdfs_.clear();

      delete promptModel_;
      delete nonPromptModel_;
      delete bkgModel_;
      delete compositeModel_;
    }

    RooArgList promptpdfpieces, nonpromptpdfpieces, bkgpdfpieces;
    RooArgList promptPolArgs, nonPromptPolArgs;
    RooArgSet pCondPdfs, pNormObs, npCondPdfs, npNormObs, bCondPdfs,bNormObs;
    RooProdPdf *promptPolProd, *nonPromptPolProd;
    RooProdPdf *promptbasepdf, *nonpromptbasepdf, *bkgbasepdf;

    if(useP) {
      if(mass_->prompt()) promptpdfpieces.add(*mass_->prompt());
      if(lifet_->prompt()) {
	//promptpdfpieces.add(*lifet_->prompt());
	pCondPdfs.add(*lifet_->prompt());
	pNormObs.add(ct);
      }
      if(promptPol_->model()) promptPolArgs.add(*promptPol_->model());
      if(acc_->prompt()) promptPolArgs.add(*acc_->prompt());
      promptPolProd = new RooProdPdf("promptPolProd","Product of Prompt Polarization PDF and Acceptance Map",promptPolArgs);
      pdfs_.push_back(promptPolProd);
      promptpdfpieces.add(*promptPolProd);
    } 

    if(useNP) {
      if(mass_->nonPrompt()) nonpromptpdfpieces.add(*mass_->nonPrompt());
      if(lifet_->nonPrompt()) {
	//nonpromptpdfpieces.add(*lifet_->nonPrompt());
	npCondPdfs.add(*lifet_->nonPrompt());
	npNormObs.add(ct);
      }
      if(nonPromptPol_->model()) nonPromptPolArgs.add(*nonPromptPol_->model());
      if(acc_->nonPrompt()) nonPromptPolArgs.add(*acc_->nonPrompt());
      nonPromptPolProd = new RooProdPdf("nonPromptPolProd","Product of NonPrompt Polarization PDF and Acceptance Map",nonPromptPolArgs);
      pdfs_.push_back(nonPromptPolProd);
      nonpromptpdfpieces.add(*nonPromptPolProd);
    }

    if(useBkg) {
      if(mass_->background()) bkgpdfpieces.add(*mass_->background());
      if(lifet_->background()) {
	//bkgpdfpieces.add(*lifet_->background());
	bCondPdfs.add(*lifet_->background());
	bNormObs.add(ct);
      }
      if(acc_->background()) bkgpdfpieces.add(*acc_->background());
    }

    promptbasepdf = new RooProdPdf("prodPromptPdf","Product of all Prompt PDFs",
				   promptpdfpieces,
				   RooFit::Conditional(pCondPdfs,pNormObs));
    pdfs_.push_back(promptbasepdf);
    nonpromptbasepdf = new RooProdPdf("prodNonPromptPdf","Product of all Non-Prompt PDFs",
				      nonpromptpdfpieces,
				      RooFit::Conditional(npCondPdfs,npNormObs));
    pdfs_.push_back(nonpromptbasepdf);
    bkgbasepdf = new RooProdPdf("prodBkgPdf","Product of all Background PDFs",
				bkgpdfpieces,
				RooFit::Conditional(bCondPdfs,bNormObs));
    pdfs_.push_back(bkgbasepdf);
      
    //create the extended PDFs for prompt, nonprompt and bkg
    if(!range_.size()) {
      promptModel_ = new RooExtendPdf("promptModel","The Prompt Contribution to the PDF",
				      *promptbasepdf,*prompt);
      nonPromptModel_ = new RooExtendPdf("nonPromptModel","The Non-Prompt Contribution to the PDF",
				       *nonpromptbasepdf,*nonprompt);
      bkgModel_ = new RooExtendPdf("bkgModel","The Background Contribution to the PDF",
				   *bkgbasepdf,*bkg);
    } else {
      promptModel_ = new RooExtendPdf("promptModel","The Prompt Contribution to the PDF",
				      *promptbasepdf,*prompt,range_.c_str());
      nonPromptModel_ = new RooExtendPdf("nonPromptModel","The Non-Prompt Contribution to the PDF",
					 *nonpromptbasepdf,*nonprompt,range_.c_str());
      bkgModel_ = new RooExtendPdf("bkgModel","The Background Contribution to the PDF",
				   *bkgbasepdf,*bkg,range_.c_str());
    }
    
    RooArgList pdfPieces;
    if(useP) pdfPieces.add(*promptModel_);
    if(useNP) pdfPieces.add(*nonPromptModel_);
    if(useBkg) pdfPieces.add(*bkgModel_);

    //create the full PDF for fitting
    compositeModel_ = new RooAddPdf("compositePDF","The Composite Model for Fitting J/Psi Mass/Lifetime/Polarization",
				    pdfPieces);
    
  }

  
  
  void CompositeModelBuilder::Print() {
    if(compositeModel_) compositeModel_->Print();
    if(promptModel_) promptModel_->Print();
    if(nonPromptModel_) nonPromptModel_->Print();
    if(bkgModel_) bkgModel_->Print();    
    for(std::vector<RooAbsPdf*>::const_iterator i = pdfs_.begin();
	i != pdfs_.end();
	++i)
      (*i)->Print();
    for(std::map<std::string,RooRealVar*>::const_iterator i = vars_.begin();
	i != vars_.end();
	++i) {
      i->second->printValue(std::cout);
      std::cout << std::endl;
    }
  }

}
