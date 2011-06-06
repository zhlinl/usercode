#ifndef _JPsiPolarization_LifetimeModel_h_
#define _JPsiPolarization_LifetimeModel_h_

#include <vector>
#include <map>
#include <string>

class RooRealVar;
class RooAddModel;
class RooAddPdf;
class RooAbsPdf;
class RooHistPdf;
class RooDataHist;
class TH1F;
class RooDecay;
class RooGaussModel;
class TDirectory;

namespace JPsiPolarization {
  class LifetimeModel {
  public:
    LifetimeModel();
    ~LifetimeModel();

    void setTemplateMethod(const bool& f=true) { isTemplate_ = f; }
    
    void loadParameters(TDirectory&);
    void saveParameters(TDirectory&);
    
    void fix(const std::string& s = "");
    void unfix(const std::string& s = "");

    void initModels(RooRealVar&,
		    const TH1F* tempPr = NULL,
		    const TH1F* tempNPr = NULL,
		    const TH1F* tempBkg = NULL,
		    const bool& p=true, 
		    const bool& np=true, 
		    const bool& bkg=true);

    void initModels(RooRealVar& ct,
		    RooRealVar& ct_err,
		    const TH1F* tempPr = NULL,
		    const TH1F* tempNPr = NULL,
		    const TH1F* tempBkg = NULL,
		    const bool& p=true, 
		    const bool& np=true, 
		    const bool& bkg=true);

    RooAbsPdf* prompt() { return promptModel_; }
    RooAbsPdf* nonPrompt() { return nonPromptModel_; }
    RooAbsPdf* background() { return bkgModel_; }

    RooGaussModel* promptGauss1() { return ctGaussPrompt1; }
    RooGaussModel* promptGauss2() { return ctGaussPrompt2; }
    RooGaussModel* promptGauss3() { return ctGaussPrompt3; }

    RooDecay* nonPromptDecayss() { return ssDecayNP; }
    RooDecay* nonPromptDecayf() { return fDecayNP; }
    RooDecay* nonPromptDecayds() { return dsDecayNP; }

    RooDecay* backgroundDecayss() { return ssDecayBK; }
    RooDecay* backgroundDecayf() { return fDecayBK; }
    RooDecay* backgroundDecayds() { return dsDecayBK; }

  private:
    LifetimeModel(const LifetimeModel& p) {}
    LifetimeModel& operator=(const LifetimeModel& p) {}
    
    void initPromptModel(RooRealVar&, const TH1F* temp = NULL);
    void initPromptModel(RooRealVar&,RooRealVar&, const TH1F* temp = NULL);
    void initNonPromptModel(RooRealVar&, const TH1F* temp = NULL);
    void initBackgroundModel(RooRealVar&, const TH1F* temp = NULL);

    std::map<std::string,RooRealVar*> vars_;
    std::map<std::string,RooAddModel*> addmodels_;
    std::map<std::string,RooDataHist*> templates_;
    std::vector<RooAddPdf*> addpdfs_;
    std::vector<RooGaussModel*> gmodels_;
    std::vector<RooDecay*> dmodels_;
    
    bool isTemplate_,varsLoaded;

    RooAbsPdf *promptModel_;
    RooAbsPdf *nonPromptModel_;
    RooAbsPdf *bkgModel_;
    
    RooGaussModel *ctGaussPrompt1;
    RooGaussModel *ctGaussPrompt2;
    RooGaussModel *ctGaussPrompt3;

    RooDecay *ssDecayBK;
    RooDecay *fDecayBK;
    RooDecay *dsDecayBK;

    RooDecay *ssDecayNP;
    RooDecay *fDecayNP;
    RooDecay *dsDecayNP;

  };
}

#endif
