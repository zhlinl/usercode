#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

using namespace edm;
using namespace std;
using namespace reco;

class YieldTree : public edm::EDAnalyzer {
  public:
    explicit YieldTree(const edm::ParameterSet&);
    ~YieldTree();
  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    InputTag upsilonsTag;
    TTree* t;
    Float_t invariantMass;
    Float_t upsPt;
    Float_t upsRapidity;
    Float_t muPlusPt;
    Float_t muPlusEta;
    Float_t muMinusPt;
    Float_t muMinusEta;
};

YieldTree::YieldTree(const edm::ParameterSet& pset):
  upsilonsTag( pset.getUntrackedParameter<InputTag>("UpsilonCandidates", InputTag()) )
{
  edm::Service<TFileService> fs;
  t = fs->make<TTree>("upsilonYield","upsilonYield");
  t->Branch("invariantMass",&invariantMass,"invariantMass/F");
  t->Branch("upsPt",&upsPt,"upsPt/F");
  t->Branch("upsRapidity",&upsRapidity,"upsRapidity/F");
  t->Branch("muPlusPt",&muPlusPt,"muPlusPt/F");
  t->Branch("muPlusEta",&muPlusEta,"muPlusEta/F");
  t->Branch("muMinusPt",&muMinusPt,"muMinusPt/F");
  t->Branch("muMinusEta",&muMinusEta,"muMinusEta/F");
}

void YieldTree::beginJob(){
}

void YieldTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<CompositeCandidateCollection> upsilonsHandle;
  iEvent.getByLabel(upsilonsTag, upsilonsHandle);
  if(upsilonsHandle.isValid()){
    for(CompositeCandidateCollection::const_iterator p=upsilonsHandle->begin(); p!= upsilonsHandle->end(); ++p){
      if( p->numberOfDaughters()==2 ){
        invariantMass = p->mass();
        upsPt = p->pt();
        upsRapidity = p->rapidity();
        for(size_t i=0; i<2; i++){
          if(p->daughter(i)->charge()>0){
            muPlusPt = p->daughter(i)->pt();
            muPlusEta = p->daughter(i)->eta();
          }else{
            muMinusPt = p->daughter(i)->pt();
            muMinusEta = p->daughter(i)->eta();
          }
        }
        t->Fill();
      }
    }
  }else{
    LogError("YieldTree")<<"Could not access: "<<upsilonsTag<<" from the event!"<<endl;
  }
}

void YieldTree::endJob(){
}

YieldTree::~YieldTree(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(YieldTree);

