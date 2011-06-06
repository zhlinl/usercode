#ifndef _JPsiPolarization_MassModel_h_
#define _JPsiPolarization_MassModel_h_

#include <vector>
#include <string>
#include <map>

class RooRealVar;
class RooAddPdf;
class RooExponential;
class RooChebychev;
class RooCBShape;
class RooAbsPdf;
class TDirectory;

namespace JPsiPolarization {
  class MassModel {
  public:
    MassModel();
    ~MassModel();
    
    void loadParameters(TDirectory&);
    void saveParameters(TDirectory&);

    void fix(const std::string& s = "");
    void unfix(const std::string& s = "");

    void initModels(RooRealVar&,
		    const bool& p = true, 
		    const bool& np = true, 
		    const bool& bkg = true);
    
    RooCBShape* prompt() { return promptModel_; }
    RooCBShape* nonPrompt() { return nonpromptModel_; }
    RooExponential* background() {return bkgModel_; }
    
  private:
    MassModel(const MassModel& p) {}
    MassModel& operator=(const MassModel& p) {}
    
    void initPromptModel(RooRealVar&);
    void initNonPromptModel(RooRealVar&);
    void initBackgroundModel(RooRealVar&);

    std::map<std::string,RooRealVar*> vars_;
    std::vector<RooAbsPdf*> pdfs_;
    
    bool varsLoaded;

    RooCBShape* promptModel_;
    RooCBShape* nonpromptModel_;
    RooExponential* bkgModel_;
    
  };
}

#endif
