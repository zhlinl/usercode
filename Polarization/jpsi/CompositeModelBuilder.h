#ifndef _JPsiPolarization_CompositeModelBuilder_h_
#define _JPsiPolarization_CompositeModelBuilder_h_

#include <vector>
#include <string>
#include <map>

#include "MassModel.h"
#include "LifetimeModel.h"
#include "PolarizationModel.h"
#include "AcceptanceMaps.h"

class RooRealVar;
class RooExtendPdf;
class RooAbsPdf;
class RooAbsReal;
class RooAddPdf;
class TH2F;
class TH1F;
class TDirectory;

namespace JPsiPolarization {
  class CompositeModelBuilder {
  public:
    // provide NULL histogram pointers to skip using that acceptance map
    CompositeModelBuilder(const std::string& frame = "",
			  const std::string& range = "");

    ~CompositeModelBuilder();

    void loadParameters(TDirectory&);
    void saveParameters(TDirectory&);

    void fix(const std::string&);
    void unfix(const std::string&);

    void saveParameter(const std::string&,TDirectory&);

    void setPromptAccHist(const TH2F * p) { promptAccHist_ = p; }
    void setNonPromptAccHist(const TH2F* p) { nonPromptAccHist_ = p; }
    void setBackgroundShapeHist(const TH2F* p) { bkgShapeHist_ = p; }

    void setPromptLifetimeTemplate(const TH1F * p) { promptLifetimeTemplate_ = p; }
    void setNonPromptLifetimeTemplate(const TH1F* p) { nonPromptLifetimeTemplate_ = p; }
    void setBackgroundLifetimeTemplate(const TH1F* p) { bkgLifetimeTemplate_ = p; }

    void setUsePrompt(const bool& use=true)    { useP = use; }   
    void setUseNonPrompt(const bool& use=true) { useNP = use; }
    void setUseBkg(const bool& use=true)       { useBkg = use; } 
    
    void setUseMass(const bool& use=true) { useMass = use; }
    void setUseLifetime(const bool& use=true) { useLifetime = use; }
    void setUsePol(const bool& use=true) { usePol = use; }
    void setUseAcceptanceMaps(const bool& use=true) { useAcc = use; }
    
    void initModel(RooRealVar& mass,
		   RooRealVar& ct,
		   RooRealVar& costh,
		   RooRealVar& phi);

    void initModel(RooRealVar& mass,
		   RooRealVar& ct,
		   RooRealVar& cterr,
		   RooRealVar& costh,
		   RooRealVar& phi);

    //for special manipulation of the consituent models
    MassModel* getMassModel() { return mass_; }
    LifetimeModel* getLifetimeModel() { return lifet_; }
    PolarizationModel* getPromptPolarizationModel() { return promptPol_; }
    PolarizationModel* getNonPromptPolarizationModel() { return nonPromptPol_; }
    AcceptanceMaps* getAcceptanceMaps() { return acc_; }

    RooAddPdf* model() { return compositeModel_; }
    RooExtendPdf* getPromptModel() { return promptModel_; }
    RooExtendPdf* getNonPromptModel() { return nonPromptModel_; }
    RooExtendPdf* getBkgModel() { return bkgModel_; }
    
    RooRealVar* promptNorm() { return vars_["nPrompt"]; }
    RooRealVar* nonPromptNorm() { return vars_["nNonPrompt"]; }
    RooRealVar* backgroundNorm() { return vars_["nBackground"]; }

    void Print();

  private:    
    CompositeModelBuilder(const CompositeModelBuilder&) {}
    CompositeModelBuilder& operator=(const CompositeModelBuilder&) {}

    void buildModel();
    void buildModel(RooRealVar&);

    std::map<std::string,RooRealVar*> vars_;
    std::vector<RooAbsPdf*> pdfs_;
    
    std::string fName_;
    std::string range_;

    bool useP, useNP, useBkg;
    bool useMass,useLifetime,usePol,useAcc;

    MassModel* mass_;
    LifetimeModel* lifet_;
    PolarizationModel* promptPol_;
    PolarizationModel* nonPromptPol_;
    AcceptanceMaps* acc_;    

    const TH2F *promptAccHist_;
    const TH2F *nonPromptAccHist_;
    const TH2F *bkgShapeHist_;

    const TH1F *promptLifetimeTemplate_;
    const TH1F *nonPromptLifetimeTemplate_;
    const TH1F *bkgLifetimeTemplate_;

    // themse are all really prodpdfs wrapped in extendpdfs
    RooExtendPdf* promptModel_;
    RooExtendPdf* nonPromptModel_;
    RooExtendPdf* bkgModel_;
    // the full model for fitting
    RooAddPdf*  compositeModel_;
			  
  };
}

#endif
