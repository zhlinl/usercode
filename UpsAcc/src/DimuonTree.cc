#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "MRudolph/UpsAcc/interface/DimuonTree.h"

using namespace edm;
using namespace reco;
using namespace std;

DimuonTree::DimuonTree(const edm::ParameterSet& iConfig){
  theRootFileName = iConfig.getParameter<string>("OutputFile");
  hFile =  new TFile(theRootFileName.c_str(),"RECREATE");
}

void DimuonTree::beginRun(edm::Run &r, const edm::EventSetup & c) {
}

void DimuonTree::beginJob() {
  //cout << "doing begin job" << endl;
  hFile->cd("");
  tree_Upsilon = new TTree("UpsTree","UpsTree");
  tree_Upsilon->Branch("genMuSize",&genMuSize,"genMuSize/I");
  tree_Upsilon->Branch("genUpsSize",&genUpsSize,"genUpsSize/I");
  tree_Upsilon->Branch("genPhotonSize",&genPhotonSize,"genPhotonSize/I");
  tree_Upsilon->Branch("recoMuSize",&recoMuSize,"recoMuSize/I");
  tree_Upsilon->Branch("genMuPt",genMuPt,"genMuPt[genMuSize]/F");
  tree_Upsilon->Branch("genMuEta",genMuEta,"genMuEta[genMuSize]/F");
  tree_Upsilon->Branch("genMuPhi",genMuPhi,"genMuPhi[genMuSize]/F");
  tree_Upsilon->Branch("genMuCharge",genMuCharge,"genMuCharge[genMuSize]/I");
  tree_Upsilon->Branch("genUpsPt",genUpsPt,"genUpsPt[genUpsSize]/F");
  tree_Upsilon->Branch("genUpsEta",genUpsEta,"genUpsEta[genUpsSize]/F");
  tree_Upsilon->Branch("genUpsPhi",genUpsPhi,"genUpsPhi[genUpsSize]/F");
  tree_Upsilon->Branch("genUpsRapidity",genUpsRapidity,"genUpsRapidity[genUpsSize]/F");
  tree_Upsilon->Branch("genPhotonEn",genPhotonEn,"genPhotonEn[genPhotonSize]/F");
  tree_Upsilon->Branch("recoMuPt",recoMuPt,"recoMuPt[recoMuSize]/F");
  tree_Upsilon->Branch("recoMuEta",recoMuEta,"recoMuEta[recoMuSize]/F");
  tree_Upsilon->Branch("recoMuPhi",recoMuPhi,"recoMuPhi[recoMuSize]/F");
  tree_Upsilon->Branch("recoMuD0",recoMuD0,"recoMuD0[recoMuSize]/F");
  tree_Upsilon->Branch("recoMuDZ",recoMuDZ,"recoMuDZ[recoMuSize]/F");
  tree_Upsilon->Branch("recoMuCharge",recoMuCharge,"recoMuCharge[recoMuSize]/I");
}

void DimuonTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<BeamSpot> hBSProd;
  iEvent.getByLabel(InputTag("offlineBeamSpot"), hBSProd);
  math::XYZPoint bspnt;
  if(hBSProd.isValid())
    bspnt = hBSProd->position();
  else
    bspnt.SetXYZ(0.,0.,0.);

  genMuSize = 0;
  genUpsSize = 0;
  genPhotonSize = 0;
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel(InputTag("genParticles"), genParticles);
  if(genParticles.isValid()){
    //cout << "gen particles is valid" << endl;
    for(GenParticleCollection::const_iterator p=genParticles->begin(); p!= genParticles->end(); ++p){
      if( p->status() == 1 && abs(p->pdgId()) == 13 && (getMotherId(&(*p)) == 553 || getMotherId(&(*p)) == 100553 || getMotherId(&(*p)) == 200553 || getMotherId(&(*p)) == 888553 || getMotherId(&(*p)) == 999553) ){
        genMuPt[genMuSize] = p->pt();
        genMuEta[genMuSize] = p->eta();
        genMuPhi[genMuSize] = p->phi();
        genMuCharge[genMuSize] = p->charge();
        genMuSize++;
      }else if( p->status() == 2 && (p->pdgId() == 553 || p->pdgId() == 100553 || p->pdgId() == 200553 || p->pdgId() == 888553  || p->pdgId() == 999553 ) ){
        genUpsPt[genUpsSize] = p->pt();
        genUpsEta[genUpsSize] = p->eta();
        genUpsPhi[genUpsSize] = p->phi();
        TLorentzVector genUps; genUps.SetPtEtaPhiM (p->pt(), p->eta(), p->phi(), p->mass());
        genUpsRapidity[genUpsSize] = genUps.Rapidity();
        genUpsSize++;
      } else if( p->pdgId() ==22) {
	genPhotonEn[genPhotonSize] = p->et();
	genPhotonSize++;
      }
    }
  }else{
    LogError("DimuonTree")<<"Could not access GenParticleCollection"<<endl;
//    exit(1);
  }

  recoMuSize = 0;
  Handle<TrackCollection> recoTracks;
  iEvent.getByLabel(InputTag("generalTracks"), recoTracks);
  if(recoTracks.isValid()){
    //cout << "reco tracks is valid " << recoTracks->size() << " tracks" << endl;
    for(TrackCollection::const_iterator tk = recoTracks->begin(); tk != recoTracks->end(); ++tk){
      //cout << "going to save mu " << recoMuSize << " " << &*tk << endl;
      recoMuPt[recoMuSize] = tk->pt();
      //cout << "saved mu " << recoMuSize << " pt" << endl;
      recoMuEta[recoMuSize] = tk->eta();
      //cout << "saved mu " << recoMuSize << " eta" << endl;
      recoMuPhi[recoMuSize] = tk->phi();
      recoMuD0[recoMuSize] = -tk->dxy(bspnt);
      recoMuDZ[recoMuSize] = tk->dz(bspnt);
      //cout << "saved mu " << recoMuSize << " phi" << endl;
      recoMuCharge[recoMuSize] = tk->charge();
      //cout << "saved mu " << recoMuSize << " q" << endl;
      recoMuSize++;
    }
  }else{
    LogDebug("DimuonTree")<<"Could not access recoTrackCollection"<<endl;
//    exit(2);
  }
  //cout << "starting fill" << endl;
  tree_Upsilon->Fill();
}

void DimuonTree::endJob(){
  hFile->Write(); 
  hFile->Close();
}

DimuonTree::~DimuonTree(){
}

int DimuonTree::getMotherId( const Candidate * p ){
  const Candidate* mother = p->mother();
  if( mother ){
    if( mother->pdgId() == p->pdgId() ){
      return getMotherId(mother);
    }else{
      return mother->pdgId();
    }
  }else{
    return 0;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(DimuonTree);
