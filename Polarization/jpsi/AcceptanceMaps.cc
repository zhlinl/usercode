#include "AcceptanceMaps.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "TH2F.h"
#include "TDirectory.h"

#include <sstream>

namespace JPsiPolarization {
  AcceptanceMaps::AcceptanceMaps(const std::string& framename) {
    fName_ = framename;
    promptAccHist_ = NULL;
    nonPromptAccHist_ = NULL;
    bkgAccHist_ = NULL;
    varsLoaded = false;
    promptAcc_ = NULL;    
    nonPromptAcc_ = NULL;
    bkgAcc_ = NULL;
  }  

  AcceptanceMaps::~AcceptanceMaps() {   
    if(promptAcc_) delete promptAcc_;
    if(nonPromptAcc_) delete nonPromptAcc_;
    if(bkgAcc_) delete bkgAcc_;    

    for(std::map<std::string,RooDataHist*>::iterator i = data_.begin();
	i != data_.end();
	++i) 
      delete i->second;
    }

  void AcceptanceMaps::loadParameters(TDirectory& dir) {

    std::stringstream hnameprompt,hnamenonprompt,hnamebkg;

    hnameprompt << "promptAccHist_" << fName_;
    hnamenonprompt << "nonPromptAccHist_" << fName_;
    hnamebkg << "bkgAccHist_" << fName_;
    
    if(!promptAccHist_) {
      dir.GetObject(hnameprompt.str().c_str(),data_[hnameprompt.str()]);
      
    }
    if(!nonPromptAccHist_)
      dir.GetObject(hnamenonprompt.str().c_str(),data_[hnamenonprompt.str()]);
    if(!bkgAccHist_)
      dir.GetObject(hnamebkg.str().c_str(),data_[hnamebkg.str()]);

    varsLoaded = true;
  }

  void AcceptanceMaps::saveParameters(TDirectory& dir) {
    for(std::map<std::string,RooDataHist*>::iterator i = data_.begin();
	i != data_.end();
	++i) {
      if(i->second)
	dir.Add(i->second->Clone(),true);
    }
    dir.Write();
  }

  void AcceptanceMaps::initAcceptanceMaps(RooRealVar& costh, RooRealVar& phi,
					  const bool& p,
					  const bool& np,
					  const bool& b) {
    std::stringstream hnameprompt,hnamenonprompt,hnamebkg;

    hnameprompt << "promptAccHist_" << fName_;
    hnamenonprompt << "nonPromptAccHist_" << fName_;
    hnamebkg << "bkgAccHist_" << fName_;

    if(p) {
      if(!promptAccHist_ && !data_[hnameprompt.str()]) {
	std::cout << "No histogram for prompt pdf defined!" << std::endl;
      } else {
      initPromptAccMap(costh,phi);
      }
    }
    if(np) {
      if(!nonPromptAccHist_ && !data_[hnamenonprompt.str()]) {
	std::cout << "No histogram for non-prompt pdf defined!" << std::endl;
      } else {
      initNonPromptAccMap(costh,phi);
      }
    }
    if(b) {
      if(!bkgAccHist_ && !data_[hnamebkg.str()]) {
	std::cout << "No histogram for background pdf defined!" << std::endl;
      } else {
      initBackgroundAccMap(costh,phi);
      }
    }
  }

  void AcceptanceMaps::initPromptAccMap(RooRealVar& costh,
					RooRealVar& phi) {
    std::stringstream name,desc;
    std::stringstream hname,hdesc;
    name << "promptAccHistPdf_" << fName_;
    desc << "Prompt Acceptance PDF for " << fName_ << " Frame";
    
    hname << "promptAccHist_" << fName_;
    hdesc << "Prompt Acceptance Histogram for " << fName_ << " Frame";
    
    RooDataHist* dataHist_;

    if(varsLoaded) {
      dataHist_ = data_[hname.str()];
    } else {
      dataHist_ = new RooDataHist(hname.str().c_str(),
				  hdesc.str().c_str(),
				  RooArgSet(costh,phi),
				  promptAccHist_);
      
      data_[hname.str()] = dataHist_;
    }
    
    promptAcc_ = new RooHistPdf(name.str().c_str(),desc.str().c_str(),
				RooArgSet(costh,phi),
				*dataHist_,
				1);
    promptAcc_->setUnitNorm(true);

    //promptAcc_->Print();
  }
  
  void AcceptanceMaps::initNonPromptAccMap(RooRealVar& costh,
					   RooRealVar& phi) {  
    std::stringstream name,desc;
    std::stringstream hname,hdesc;
    name << "nonPromptAccHistPdf_" << fName_;
    desc << "Non-Prompt Acceptance PDF for " << fName_ << " Frame";
    
    hname << "nonPromptAccHist_" << fName_;
    hdesc << "Non-Prompt Acceptance Histogram for " << fName_ << " Frame";
    
    RooDataHist* dataHist_;
    
    if(varsLoaded) {
      dataHist_ = data_[hname.str()];
    } else {
      dataHist_ = new RooDataHist(hname.str().c_str(),
				  hdesc.str().c_str(),
				  RooArgSet(costh,phi),
				  nonPromptAccHist_);
      
      data_[hname.str()] = dataHist_;
    }

    nonPromptAcc_ = new RooHistPdf(name.str().c_str(),desc.str().c_str(),
				   RooArgSet(costh,phi),
				   *dataHist_,
				   1);
    nonPromptAcc_->setUnitNorm(true);

    //nonPromptAcc_->Print();
  }

  void AcceptanceMaps::initBackgroundAccMap(RooRealVar& costh,
					    RooRealVar& phi) {  
        
    std::stringstream name,desc;
    std::stringstream hname,hdesc;
    name << "bkgAccHistPdf_" << fName_;
    desc << "Background Acceptance PDF for " << fName_ << " Frame";
    
    hname << "bkgAccHist_" << fName_;
    hdesc << "Background Acceptance Histogram for " << fName_ << " Frame";
    
    RooDataHist* dataHist_;

    if(varsLoaded) {
      dataHist_ = data_[hname.str()];
    } else {
      dataHist_ = new RooDataHist(hname.str().c_str(),
				  hdesc.str().c_str(),
				  RooArgSet(costh,phi),
				  bkgAccHist_);
      
      data_[hname.str()] = dataHist_;
    }

    bkgAcc_ = new RooHistPdf(name.str().c_str(),desc.str().c_str(),
			     RooArgSet(costh,phi),
			     *dataHist_,
			     1);
    bkgAcc_->setUnitNorm(true);
    //bkgAcc_->Print();
  }
}
