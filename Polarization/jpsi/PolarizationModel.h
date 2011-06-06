#ifndef _JPsiPolarization_PolarizationModel_h_
#define _JPsiPolarization_PolarizationModel_h_

#include <vector>
#include <string>
#include <map>

class RooRealVar;
class RooGenericPdf;
class TDirectory;

namespace JPsiPolarization {
  class PolarizationModel {
  public:
    PolarizationModel(const std::string componentname,
		      const std::string& framename);
    ~PolarizationModel();
    
    void loadParameters(TDirectory&);
    void saveParameters(TDirectory&);

    void fix(const std::string& s = "");
    void unfix(const std::string& s = "");

    void initModel(RooRealVar& costh,
		   RooRealVar& phi);
		   
    RooGenericPdf* model() { return model_; }
           
  private:
    PolarizationModel(const PolarizationModel& p) {}
    PolarizationModel& operator=(const PolarizationModel& p) {}
        
    std::string compName_,fName_;

    std::map<std::string,RooRealVar*> vars_;    
    
    bool varsLoaded;

    RooGenericPdf* model_;
  };
}

#endif
