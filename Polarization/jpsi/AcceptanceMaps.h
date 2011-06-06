#ifndef _JPsiPolarization_AcceptanceMaps_h_
#define _JPsiPolarization_AcceptanceMaps_h_

#include <vector>
#include <map>
#include <string>

class RooRealVar;
class RooDataHist;
class RooHistPdf;
class TH2F;
class TDirectory;

namespace JPsiPolarization {
  class AcceptanceMaps {
  public:
    AcceptanceMaps(const std::string& framename);
    ~AcceptanceMaps();
    
    void loadParameters(TDirectory&);
    void saveParameters(TDirectory&);

    void setPromptHist(const TH2F * p) { promptAccHist_ = p; }
    void setNonPromptHist(const TH2F* p) { nonPromptAccHist_ = p; }
    void setBackgroundHist(const TH2F* p) { bkgAccHist_ = p; }
    
    void initAcceptanceMaps(RooRealVar& costh,
			    RooRealVar& phi,
			    const bool& p=true, 
			    const bool& np=true, 
			    const bool& bkg=true);

    RooHistPdf* prompt() { return promptAcc_; }
    RooHistPdf* nonPrompt() { return nonPromptAcc_; }
    RooHistPdf* background() {return bkgAcc_; }
    
  private:
    AcceptanceMaps(const AcceptanceMaps& p) {}
    AcceptanceMaps& operator=(const AcceptanceMaps& p) {}
    
    void initPromptAccMap(RooRealVar& ,RooRealVar& );
    void initNonPromptAccMap(RooRealVar& , RooRealVar& );
    void initBackgroundAccMap(RooRealVar& , RooRealVar& );
    
    std::string fName_;

    //these are not owned by this class
    const TH2F *promptAccHist_, *nonPromptAccHist_, *bkgAccHist_;

    std::map<std::string,RooDataHist*> data_;

    bool varsLoaded;

    RooHistPdf *promptAcc_;
    RooHistPdf *nonPromptAcc_;
    RooHistPdf *bkgAcc_;
    
  };
}

#endif
