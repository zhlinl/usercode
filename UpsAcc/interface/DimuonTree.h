#ifndef DimuonTree_H
#define DimuonTree_H

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#define NMAX 1000

using namespace std;

class DimuonTree : public edm::EDAnalyzer {
  public:
    explicit DimuonTree(const edm::ParameterSet&);
    ~DimuonTree();

  private:
    virtual void beginRun(edm::Run &r, const edm::EventSetup & c); 
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob() ;
    int getMotherId(const reco::Candidate*);
    TLorentzVector lrz(float pt, float eta, float phi, float UpsMass);

    string theRootFileName;
    TFile* hFile;
    TTree* tree_Upsilon;

    Int_t genUpsSize;
    Float_t genUpsPt[NMAX];
    Float_t genUpsEta[NMAX];
    Float_t genUpsPhi[NMAX];
    Float_t genUpsRapidity[NMAX];

    Int_t genMuSize;
    Float_t genMuPt[NMAX];
    Float_t genMuEta[NMAX];
    Float_t genMuPhi[NMAX];
    Int_t genMuCharge[NMAX];

    Int_t genPhotonSize;
    Float_t genPhotonEn[NMAX];

    Int_t recoMuSize;
    Float_t recoMuPt[NMAX];
    Float_t recoMuEta[NMAX];
    Float_t recoMuPhi[NMAX];
    Float_t recoMuD0[NMAX];
    Float_t recoMuDZ[NMAX];
    Int_t recoMuCharge[NMAX];
};
#endif
