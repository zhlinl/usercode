// -*- C++ -*-
//
// Package:    UpsilonAnalyzerPAT
// Class:      UpsilonAnalyzerPAT
// 
//

// system include files
#include <memory>

// ROOT/Roofit include files
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include "TTree.h"
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h>
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/TrackReco/interface/TrackBase.h" 
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedTrackerVertex.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//
class UpsilonAnalyzerPAT : public edm::EDAnalyzer {
   public:
      explicit UpsilonAnalyzerPAT(const edm::ParameterSet&);
      ~UpsilonAnalyzerPAT();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void makeCuts() ;
      pair< unsigned int, const pat::CompositeCandidate* > theBestQQ();
      void fillHistosAndDS(unsigned int theCat, const pat::CompositeCandidate* aCand);
      bool selGlobalMuon(const pat::Muon* aMuon);
      bool selTrackerMuon(const pat::Muon* aMuon);
      bool selCaloMuon(const pat::Muon* aMuon);
      TLorentzVector lorentzMomentum(const GenParticleRef& genp) const;

      TFile* fOut;
      TTree* t;
      TTree* t1;
      TTree* t3;

      Float_t invariantMass;
      Float_t upsPt;
      Float_t upsRapidity;
      Float_t upsEta;
      Float_t upsPhi;
      Float_t genUpsP;
      Float_t genUpsPt;
      Float_t genUpsEta;
      Float_t genUpsPhi;
      Float_t genUpsRap;
      Float_t genMinMuPt;
      Float_t genMinMuEta;
      Float_t genMinMuPhi;
      Float_t genPosMuPt;
      Float_t genPosMuEta;
      Float_t genPosMuPhi;
      Float_t MCMinMuPt;
      Float_t MCPosMuPt;
      Float_t MCMinMuEta;
      Float_t MCPosMuEta;
      Float_t MCMinMuPhi;
      Float_t MCPosMuPhi;
      Float_t MCUpsPt;
      Float_t MCUpsEta;
      Float_t MCUpsPhi;
      Float_t MCUpsRap;
      Float_t MCUpsP;

      Int_t patMuCut;
      Float_t Mu3MuPt;
      Float_t Mu3MuEta;
      Float_t Mu3MuPhi;
      Float_t idMuPt;
      Float_t idMuEta;
      Float_t idMuPhi;
      Float_t L1OpenMuPt;
      Float_t L1OpenMuEta;
      Float_t L1OpenMuPhi;
      Int_t shareChamber;
 
      Float_t muPlusPt;
      Float_t muPlusP;
      Float_t muPlusEta;
      Float_t muPlusPhi;
      Float_t muMinusPt;
      Float_t muMinusP;
      Float_t muMinusEta;
      Float_t muMinusPhi;
      Float_t McInvMass;
      Int_t GenRef;
      Int_t SampleFlag;
      Int_t eventId;
      Int_t Category;
      Int_t hlt_L1MuOpen;
      Int_t hlt_Mu3;
      Int_t hlt_Mu5;
      Int_t hlt_DoubleMu0;
      Int_t hlt_DoubleMu3;
      Int_t hlt_L1DoubleMuOpen;

      Int_t L1MuOpen_Minus;
      Int_t Mu3_Minus;
      Int_t Mu5_Minus;
      Int_t DoubleMu0_Minus;
      Int_t DoubleMu3_Minus;
      Int_t L1DoubleMuOpen_Minus;

      Int_t L1MuOpen_Plus;
      Int_t Mu3_Plus;
      Int_t Mu5_Plus;
      Int_t DoubleMu0_Plus;
      Int_t DoubleMu3_Plus;
      Int_t L1DoubleMuOpen_Plus;

      Float_t vtxchi2;
      Float_t chamSegR;
      Float_t muMinusD0;
      Float_t muMinusDz;
      Int_t muMinusNhits;
      Int_t muMinusNSeg;
      Int_t muPlusNSeg;
      Int_t muMinusNPixelhits;
      Float_t muMinusTkNorChi2;
      Float_t muMinusGlbNorChi2;
      Float_t muPlusD0;
      Float_t muPlusDz;
      Int_t muPlusNhits;
      Int_t muPlusNPixelhits;
      Int_t muPlusIso03NTracks;
      Int_t muMinusIso03NTracks;
      Float_t muMinusIso03sumPt;
      Float_t muMinusIso03emEt;
      Float_t muMinusIso03hadEt;
      Float_t muMinusIso03hoEt;
      Int_t muMinusIso03NJets;
      Float_t muMinusIso03trackerVetoPt;
      Float_t muMinusIso03emVetoEt;
      Float_t muMinusIso03hadVetoEt;
      Float_t muMinusIso03hoVetoEt;
      Float_t muPlusIso03sumPt;
      Float_t muPlusIso03emEt;
      Float_t muPlusIso03hadEt;
      Float_t muPlusIso03hoEt;
      Int_t muPlusIso03NJets;
      Float_t muPlusIso03trackerVetoPt;
      Float_t muPlusIso03emVetoEt;
      Float_t muPlusIso03hadVetoEt;
      Float_t muPlusIso03hoVetoEt;

      Float_t muPlusGlbNorChi2;
      Float_t muPlusTkNorChi2;
      Int_t muMinusGlbNMuHits;
      Int_t muPlusGlbNMuHits;
      Float_t muMinusGlbTkChi2;
      Float_t muPlusGlbTkChi2;
      Int_t muMinusTrkArbirtated;
      Int_t muMinusTrkLSLoose;
      Int_t muMinusTrkLSTight;
      Int_t muMinusTrk2DLoose;
      Int_t muMinusTrk2DTight;
      Int_t muMinusTrkOSLoose;
      Int_t muMinusTrkOSTight;
      Int_t muMinusTrkLSLowPtLoose;  
      Int_t muMinusTrkLSLowPtTight;
      Int_t muMinusTrkLSAngLoose;
      Int_t muMinusTrkLSAngTight;
      Int_t muMinusTrkOSAngLoose;
      Int_t muMinusTrkOSAngTight;
      Int_t muPlusTrkArbirtated;
      Int_t muPlusTrkLSLoose;
      Int_t muPlusTrkLSTight;
      Int_t muPlusTrk2DLoose;
      Int_t muPlusTrk2DTight;
      Int_t muPlusTrkOSLoose;
      Int_t muPlusTrkOSTight;
      Int_t muPlusTrkLSLowPtLoose;
      Int_t muPlusTrkLSLowPtTight;  
      Int_t muPlusTrkLSAngLoose;
      Int_t muPlusTrkLSAngTight;
      Int_t muPlusTrkOSAngLoose;
      Int_t muPlusTrkOSAngTight;


      Int_t LS;
      Int_t FromKLoose;
      
      Handle<pat::MuonCollection> collMu;
      Handle<pat::CompositeCandidateCollection > collAll;
      Handle<pat::CompositeCandidateCollection > collCalo;
      Handle<GenParticleCollection> genParticles;
      Handle<TriggerResults> trigger;
      Handle<VertexCollection> priVtxs;
      Handle<PFDisplacedVertexCollection> PFDVtxs;
      // data members
      InputTag       _patMu;
      InputTag       _patUpsilon;
      InputTag       _patUpsilonWithCalo;
      InputTag       _genParticle;
      string         _theRootFile;
      int            _whichsample;
      bool           _onlythebest;
      bool           _applycuts;
      bool           _storeefficiency;
      bool           _useBS;
      bool           _useCalo;
      bool           _useMC; 
     bool           _removeSignal;
      bool           _RemovePairsSharingTrackerTrack;
      bool           _RemovePairsSharingChamber;
      bool           OS;
      bool           trackeronly;
      InputTag       _triggerresults;
      vector<unsigned int>                     _thePassedCats;
      vector<const pat::CompositeCandidate*>   _thePassedCands;

      // number of events
      unsigned int nEvents;
      unsigned int passedCandidates;


};    
      
//
UpsilonAnalyzerPAT::UpsilonAnalyzerPAT(const edm::ParameterSet& iConfig):
  _patMu(iConfig.getParameter<InputTag>("srcMu")),
  _patUpsilon(iConfig.getParameter<InputTag>("src")),
  _patUpsilonWithCalo(iConfig.getParameter<InputTag>("srcWithCaloMuons")),
  _genParticle(iConfig.getParameter<InputTag>("srcWithGenParticle")),
  _theRootFile(iConfig.getParameter<string>("theRootFile")),
  _whichsample(iConfig.getParameter<int>("whichsample")),
  _onlythebest(iConfig.getParameter<bool>("onlyTheBest")),		
  _applycuts(iConfig.getParameter<bool>("applyCuts")),			
  _storeefficiency(iConfig.getParameter<bool>("storeEfficiency")),	
  _useBS(iConfig.getParameter<bool>("useBeamSpot")),
  _useCalo(iConfig.getUntrackedParameter<bool>("useCaloMuons",false)),
  _useMC(iConfig.getUntrackedParameter<bool>("useMC",false)),
  _removeSignal(iConfig.getUntrackedParameter<bool>("removeSignalEvents",false)),
  _RemovePairsSharingTrackerTrack(iConfig.getUntrackedParameter<bool>("RemovePairsSharingTrackerTrack",false)),
  _RemovePairsSharingChamber(iConfig.getUntrackedParameter<bool>("RemovePairsSharingChamber",false)),
  OS(iConfig.getUntrackedParameter<bool>("OS",false)),
  trackeronly(iConfig.getUntrackedParameter<bool>("trackeronly",false)),
  _triggerresults(iConfig.getParameter<InputTag>("TriggerResultsLabel"))
{
   //now do what ever initialization is needed
  nEvents = 0;
  passedCandidates = 0;
  fOut = new TFile(_theRootFile.c_str(), "RECREATE");
  t3 = new TTree("PatMuTree","PatMuTree");
  t3->Branch("idMuPt",&idMuPt,"idMuPt/F");
  t3->Branch("idMuEta",&idMuEta,"idMuEta/F");
  t3->Branch("idMuPhi",&idMuPhi,"idMuPhi/F");
  t3->Branch("Mu3MuPt",&Mu3MuPt,"Mu3MuPt/F");
  t3->Branch("Mu3MuEta",&Mu3MuEta,"Mu3MuEta/F");
  t3->Branch("Mu3MuPhi",&Mu3MuPhi,"Mu3MuPhi/F");
  t3->Branch("L1OpenMuPt",&L1OpenMuPt,"L1OpenMuPt/F");
  t3->Branch("L1OpenMuEta",&L1OpenMuEta,"L1OpenMuEta/F");
  t3->Branch("L1OpenMuPhi",&L1OpenMuPhi,"L1OpenMuPhi/F");
  t3->Branch("patMuCut",&patMuCut,"patMuCut/I");

  t1 = new TTree("genTree","genTree");
  t1->Branch("MCPosMuPhi",&MCPosMuPhi,"MCPosMuPhi/F");
  t1->Branch("MCPosMuEta",&MCPosMuEta,"MCPosMuEta/F");
  t1->Branch("MCPosMuPt",&MCPosMuPt,"MCPosMuPt/F");
  t1->Branch("MCMinMuPt",&MCMinMuPt,"MCMinMuPt/F");
  t1->Branch("MCMinMuEta",&MCMinMuEta,"MCMinMuEta/F");
  t1->Branch("MCMinMuPhi",&MCMinMuPhi,"MCMinMuPhi/F");
  t1->Branch("MCUpsPhi",&MCUpsPhi,"MCUpsPhi/F");
  t1->Branch("MCUpsPt",&MCUpsPt,"MCUpsPt/F");
  t1->Branch("MCUpsP",&MCUpsP,"MCUpsP/F");
  t1->Branch("MCUpsEta",&MCUpsEta,"MCUpsEta/F");
  t1->Branch("MCUpsRap",&MCUpsRap,"MCUpsRap/F");
  
  t = new TTree("UpsilonTree","UpsilonTree");
  t->Branch("shareChamber",&shareChamber,"shareChamber/I");
  t->Branch("invariantMass",&invariantMass,"invariantMass/F");
  t->Branch("upsPt",&upsPt,"upsPt/F");
  t->Branch("upsRapidity",&upsRapidity,"upsRapidity/F");
  t->Branch("upsEta",&upsEta,"upsEta/F");
  t->Branch("upsPhi",&upsPhi,"upsPhi/F");
  t->Branch("genUpsP",&genUpsP,"genUpsP/F");
  t->Branch("genUpsPt",&genUpsPt,"genUpsPt/F");
  t->Branch("genUpsRap",&genUpsRap,"genUpsRap/F");
  t->Branch("genUpsEta",&genUpsEta,"genUpsEta/F");
  t->Branch("genUpsPhi",&genUpsPhi,"genUpsPhi/F"); 
  t->Branch("genMinMuPt",&genMinMuPt,"genMinMuPt/F");
  t->Branch("genMinMuEta",&genMinMuEta,"genMinMuEta/F");
  t->Branch("genMinMuPhi",&genMinMuPhi,"genMinMuPhi/F"); 
  t->Branch("genPosMuPt",&genPosMuPt,"genPosMuPt/F");
  t->Branch("genPosMuEta",&genPosMuEta,"genPosMuEta/F");
  t->Branch("genPosMuPhi",&genPosMuPhi,"genPosMuPhi/F");
  t->Branch("upsRapidity",&upsRapidity,"upsRapidity/F");
  t->Branch("muPlusPt",&muPlusPt,"muPlusPt/F");
  t->Branch("muPlusP",&muPlusP,"muPlusP/F");
  t->Branch("muPlusEta",&muPlusEta,"muPlusEta/F");
  t->Branch("muPlusPhi",&muPlusPhi,"muPlusPhi/F");
  t->Branch("muMinusP",&muMinusP,"muMinusP/F");
  t->Branch("muMinusPt",&muMinusPt,"muMinusPt/F");
  t->Branch("muMinusEta",&muMinusEta,"muMinusEta/F");
  t->Branch("muMinusPhi",&muMinusPhi,"muMinusPhi/F");
  t->Branch("muPlusTrkArbirtated",&muPlusTrkArbirtated,"muPlusTrkArbirtated/I");
  t->Branch("muPlusTrkLSLoose",&muPlusTrkLSLoose,"muPlusTrkLSLoose/I");
  t->Branch("muPlusTrkLSTight",&muPlusTrkLSTight,"muPlusTrkLSTight/I");
  t->Branch("muPlusTrk2DLoose",&muPlusTrk2DLoose,"muPlusTrk2DLoose/I");
  t->Branch("muPlusTrk2DTight",&muPlusTrk2DTight,"muPlusTrk2DTight/I");
  t->Branch("muPlusTrkOSLoose",&muPlusTrkOSLoose,"muPlusTrkOSLoose/I");
  t->Branch("muPlusTrkOSTight",&muPlusTrkOSTight,"muPlusTrkOSTight/I");
  t->Branch("muPlusTrkLSLowPtLoose",&muPlusTrkLSLowPtLoose,"muPlusTrkLSLowPtLoose/I");
  t->Branch("muPlusTrkLSLowPtTight",&muPlusTrkLSLowPtTight,"muPlusTrkLSLowPtTight/I");
  t->Branch("muPlusTrkLSAngLoose",&muPlusTrkLSAngLoose,"muPlusTrkLSAngLoose/I");
  t->Branch("muPlusTrkLSAngTight",&muPlusTrkLSAngTight,"muPlusTrkLSAngTight/I");
  t->Branch("muPlusTrkOSAngLoose",&muPlusTrkOSAngLoose,"muPlusTrkOSAngLoose/I");
  t->Branch("muPlusTrkOSAngTight",&muPlusTrkOSAngTight,"muPlusTrkOSAngTight/I");

  t->Branch("McInvMass",&McInvMass,"McInvMass/F");
  t->Branch("GenRef",&GenRef,"GenRef/I");
  t->Branch("LS",&LS,"LS/I");
  t->Branch("FromKLoose", &FromKLoose, "FromKLoose/I");
  t->Branch("SampleFlag",&SampleFlag,"SampleFlag/I");
  t->Branch("Category",&Category,"Category/I");
  t->Branch("eventId",&eventId,"eventId/I");
  t->Branch("hlt_L1MuOpen",&hlt_L1MuOpen,"hlt_L1MuOpen/I");
  t->Branch("hlt_Mu3",&hlt_Mu3,"hlt_Mu3/I");
  t->Branch("hlt_Mu5",&hlt_Mu5,"hlt_Mu5/I");
  t->Branch("hlt_DoubleMu0",&hlt_DoubleMu0,"hlt_DoubleMu0/I");
  t->Branch("hlt_DoubleMu3",&hlt_DoubleMu3,"hlt_DoubleMu3/I");
  t->Branch("hlt_L1DoubleMuOpen",&hlt_L1DoubleMuOpen,"hlt_L1DoubleMuOpen/I");
  t->Branch("L1MuOpen_Plus",&L1MuOpen_Plus,"L1MuOpen_Plus/I");
  t->Branch("Mu3_Plus",&Mu3_Plus,"Mu3_Plus/I");
  t->Branch("Mu5_Plus",&Mu5_Plus,"Mu5_Plus/I");
  t->Branch("DoubleMu0_Plus",&DoubleMu0_Plus,"DoubleMu0_Plus/I");
  t->Branch("DoubleMu3_Plus",&DoubleMu3_Plus,"DoubleMu3_Plus/I");
  t->Branch("L1DoubleMuOpen_Plus",&L1DoubleMuOpen_Plus,"L1DoubleMuOpen_Plus/I");
  t->Branch("L1MuOpen_Minus",&L1MuOpen_Minus,"L1MuOpen_Minus/I");
  t->Branch("Mu3_Minus",&Mu3_Minus,"Mu3_Minus/I");
  t->Branch("Mu5_Minus",&Mu5_Minus,"Mu5_Minus/I");
  t->Branch("DoubleMu0_Minus",&DoubleMu0_Minus,"DoubleMu0_Minus/I");
  t->Branch("DoubleMu3_Minus",&DoubleMu3_Minus,"DoubleMu3_Minus/I");
  t->Branch("L1DoubleMuOpen_Minus",&L1DoubleMuOpen_Minus,"L1DoubleMuOpen_Minus/I");

  t->Branch("vtxchi2",&vtxchi2,"vtxchi2/F");
  t->Branch("chamSegR",&chamSegR,"chamSegR/F");
  t->Branch("muMinusD0",&muMinusD0,"muMinusD0/F");
  t->Branch("muMinusDz",&muMinusDz,"muMinusDz/F");
  t->Branch("muMinusNhits",&muMinusNhits,"muMinusNhits/I");
  t->Branch("muMinusNSeg",&muMinusNSeg,"muMinusNSeg/I");
   t->Branch("muPlusNSeg",&muPlusNSeg,"muPlusNSeg/I");
 t->Branch("muMinusNPixelhits",&muMinusNPixelhits,"muMinusNPixelhits/I");
  t->Branch("muPlusGlbNMuHits",&muPlusGlbNMuHits,"muPlusGlbNMuHits/I");
  t->Branch("muMinusGlbNMuHits",&muMinusGlbNMuHits,"muMinusGlbNMuHits/I");
  t->Branch("muMinusGlbNorChi2",&muMinusGlbNorChi2,"muMinusGlbNorChi2/F");
  t->Branch("muMinusTkNorChi2",&muMinusTkNorChi2,"muMinusTkNorChi2/F");
  t->Branch("muMinusGlbTkChi2",&muMinusGlbTkChi2,"muMinusGlbTkChi2/F");
  t->Branch("muMinusTrkArbirtated",&muMinusTrkArbirtated,"muMinusTrkArbirtated/I");
  t->Branch("muMinusTrkLSLoose",&muMinusTrkLSLoose,"muMinusTrkLSLoose/I");
  t->Branch("muMinusTrkLSTight",&muMinusTrkLSTight,"muMinusTrkLSTight/I");
  t->Branch("muMinusTrk2DLoose",&muMinusTrk2DLoose,"muMinusTrk2DLoose/I");
  t->Branch("muMinusTrk2DTight",&muMinusTrk2DTight,"muMinusTrk2DTight/I");
  t->Branch("muMinusTrkOSLoose",&muMinusTrkOSLoose,"muMinusTrkOSLoose/I");
  t->Branch("muMinusTrkOSTight",&muMinusTrkOSTight,"muMinusTrkOSTight/I");
  t->Branch("muMinusTrkLSLowPtLoose",&muMinusTrkLSLowPtLoose,"muMinusTrkLSLowPtLoose/I");
  t->Branch("muMinusTrkLSLowPtTight",&muMinusTrkLSLowPtTight,"muMinusTrkLSLowPtTight/I");
  t->Branch("muMinusTrkLSAngLoose",&muMinusTrkLSAngLoose,"muMinusTrkLSAngLoose/I");
  t->Branch("muMinusTrkLSAngTight",&muMinusTrkLSAngTight,"muMinusTrkLSAngTight/I");
  t->Branch("muMinusTrkOSAngLoose",&muMinusTrkOSAngLoose,"muMinusTrkOSAngLoose/I");
  t->Branch("muMinusTrkOSAngTight",&muMinusTrkOSAngTight,"muMinusTrkOSAngTight/I");
  t->Branch("muPlusD0",&muPlusD0,"muPlusD0/F");
  t->Branch("muPlusDz",&muPlusDz,"muPlusDz/F");
  t->Branch("muPlusNhits",&muPlusNhits,"muPlusNhits/I");
  t->Branch("muPlusNPixelhits",&muPlusNPixelhits,"muPlusNPixelhits/I");
  t->Branch("muPlusGlbNorChi2",&muPlusGlbNorChi2,"muPlusGlbNorChi2/F");
  t->Branch("muPlusGlbTkChi2",&muPlusGlbTkChi2,"muPlusGlbTkChi2/F");
  t->Branch("muPlusTkNorChi2",&muPlusTkNorChi2,"muPlusTkNorChi2/F");
  t->Branch("muPlusIso03NTracks",&muPlusIso03NTracks,"muPlusIso03NTracks/I");
  t->Branch("muPlusIso03sumPt",&muPlusIso03sumPt,"muPlusIso03sumPt/F");
  t->Branch("muPlusIso03emEt",&muPlusIso03emEt,"muPlusIso03emEt/F");
  t->Branch("muPlusIso03hadEt",&muPlusIso03hadEt,"muPlusIso03hadEt/F");
  t->Branch("muPlusIso03hoEt",&muPlusIso03hoEt,"muPlusIso03hoEt/F");
  t->Branch("muPlusIso03NJets",&muPlusIso03NJets,"muPlusIso03NJets/I");
  t->Branch("muPlusIso03trackerVetoPt",&muPlusIso03trackerVetoPt,"muPlusIso03trackerVetoPt/F");
  t->Branch("muPlusIso03emVetoEt",&muPlusIso03emVetoEt,"muPlusIso03emVetoEt/F");
  t->Branch("muPlusIso03hadVetoEt",&muPlusIso03hadVetoEt,"muPlusIso03hadVetoEt/F");
  t->Branch("muPlusIso03hoVetoEt",&muPlusIso03hoVetoEt,"muPlusIso03hoVetoEt/F");

  t->Branch("muMinusIso03NTracks",&muMinusIso03NTracks,"muMinusIso03NTracks/I");
  t->Branch("muMinusIso03sumPt",&muMinusIso03sumPt,"muMinusIso03sumPt/F");
  t->Branch("muMinusIso03emEt",&muMinusIso03emEt,"muMinusIso03emEt/F");
  t->Branch("muMinusIso03hadEt",&muMinusIso03hadEt,"muMinusIso03hadEt/F");
  t->Branch("muMinusIso03hoEt",&muMinusIso03hoEt,"muMinusIso03hoEt/F");
  t->Branch("muMinusIso03NJets",&muMinusIso03NJets,"muMinusIso03NJets/I");
  t->Branch("muMinusIso03trackerVetoPt",&muMinusIso03trackerVetoPt,"muMinusIso03trackerVetoPt/F");
  t->Branch("muMinusIso03emVetoEt",&muMinusIso03emVetoEt,"muMinusIso03emVetoEt/F");
  t->Branch("muMinusIso03hadVetoEt",&muMinusIso03hadVetoEt,"muMinusIso03hadVetoEt/F");
  t->Branch("muMinusIso03hoVetoEt",&muMinusIso03hoVetoEt,"muMinusIso03hoVetoEt/F");

cout<<"done"<<endl;

}


UpsilonAnalyzerPAT::~UpsilonAnalyzerPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
UpsilonAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   nEvents++;

   iEvent.getByLabel(_triggerresults,trigger);

   try {iEvent.getByLabel(_patMu,collMu);}
   catch (...) {cout << "Mu not present in event!" << endl;}

   try {iEvent.getByLabel(_patUpsilon,collAll);} 
   catch (...) {cout << "upsilon not present in event!" << endl;}
   try {iEvent.getByLabel("offlinePrimaryVertices",priVtxs);}
   catch (...) {cout << "vtx not present in event!" << endl;}

   try {iEvent.getByLabel("particleFlowDisplacedVertex",PFDVtxs);}
   catch (...) {cout << "Displaced vtx not present in event!" << endl;}

   if (_useCalo) {
     try {iEvent.getByLabel(_patUpsilonWithCalo,collCalo);} 
     catch (...) {cout << "upsilon to calomuons not present in event!" << endl;}
   }
 
   if (_useMC){
      try {iEvent.getByLabel(_genParticle, genParticles);}
      catch (...) {cout << "genParticle not pressent in event!" << endl;}
   }
 
   _thePassedCats.clear();
   _thePassedCands.clear();

   SampleFlag = _whichsample;
   // APPLY CUTS
   cout<<"make cuts"<<endl;
   cout<<"event id"<<iEvent.id().event()<<endl;
   eventId = iEvent.id().event();
   this->makeCuts();

   // BEST Upsilon?
  if(collMu.isValid()){
    for(vector<pat::Muon>::const_iterator it=collMu->begin();
        it!=collMu->end();++it) {
        reco::GenParticleRef genMu = it->genParticleRef();
       const pat::Muon* it2 = dynamic_cast<const pat::Muon*>(&(*it));
       cout<<"before selection"<<endl;
       if(it2->muonID("TMLastStationAngTight")){ 
         cout<<"pass sele"<<endl;
        if(genMu.isAvailable() && abs(genMu->pdgId()) == 13){
           cout<<"genMu available"<<endl;
           if(genMu->numberOfMothers() > 0 &&  (genMu->motherRef()->pdgId() == 553 || genMu->motherRef()->pdgId() == 100553)){
             idMuPt = genMu->pt();
             idMuEta = genMu->eta();
             idMuPhi = genMu->phi();
             cout<<"fill idMu"<<endl;
             if(selTrackerMuon(it2)) patMuCut = 1;
             else patMuCut = 0;  
           if(!it->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()){
                cout<<"trigger"<<endl;
                L1OpenMuPt = genMu->pt();
                L1OpenMuEta = genMu->eta();
                L1OpenMuPhi = genMu->phi();
             }else{
                cout<<"trigger"<<endl;
                L1OpenMuPt = 0;
                L1OpenMuEta = 0;
                L1OpenMuPhi = 0;
             }
             if(!it->triggerObjectMatchesByPath("HLT_Mu3").empty()){
                Mu3MuPt = genMu->pt();
                Mu3MuEta = genMu->eta();
                Mu3MuPhi = genMu->phi();
             }else{
                Mu3MuPt = 0;
                Mu3MuEta = 0;
                Mu3MuPhi = 0;
             }    
           t3->Fill();
           }
        }
       }
    }
  }
  if(genParticles.isValid()){
  cout<<"generation valid"<<endl;
    for(GenParticleCollection::const_iterator p=genParticles->begin(); p!= genParticles->end(); ++p){
        if( p->status() == 1 && p->pdgId()==13 && (p->mother()->pdgId()== 553 || p->mother()->pdgId()== 100553 || p->mother()->pdgId()== 443)){
            MCMinMuPt = p->pt();
            MCMinMuEta = p->eta();
            MCMinMuPhi = p->phi();
        }else if( p->status() == 1 && p->pdgId()== -13 && (p->mother()->pdgId()== 553 || p->mother()->pdgId()== 100553 || p->mother()->pdgId()== 443)){
            MCPosMuPt = p->pt();
            MCPosMuEta = p->eta();
            MCPosMuPhi = p->phi();
        }else if( p->status() == 2 && (abs(p->pdgId()) == 553  || abs(p->pdgId()) == 100553)){
            MCUpsPt = p->pt();
            MCUpsEta = p->eta();
            MCUpsPhi = p->phi();
            MCUpsRap = p->rapidity();
        }
//    t1->Fill();
   }
  } 
    t1->Fill();
   if (_onlythebest) {  // yes, fill simply the best

     pair< unsigned int, const pat::CompositeCandidate* > theBest = theBestQQ();
     if (theBest.first < 10) fillHistosAndDS(theBest.first, theBest.second);

   } else {   // no, fill all passing cuts
    cout<<"get in"<<endl;
     for( unsigned int count = 0; count < _thePassedCands.size(); count++) { 
         cout<<"fill"<<_thePassedCands.size()<<"cat"<<_thePassedCats.at(count)<<"cand"<<_thePassedCands.at(count)<<endl;
       fillHistosAndDS(_thePassedCats.at(count), _thePassedCands.at(count)); 
     }

   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
UpsilonAnalyzerPAT::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
UpsilonAnalyzerPAT::endJob() {
 fOut->Write();
 fOut->Close(); 
}

void 
UpsilonAnalyzerPAT::fillHistosAndDS(unsigned int theCat, const pat::CompositeCandidate* aCand){
 
//  if(theCat == 2){ 

//  cout<<"event id"<<

  if(theCat == 0) Category = 0;
  else if (theCat == 1) Category = 1;
  else if (theCat == 2) Category = 2;
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));

        
  float theMass = aCand->mass();
   vtxchi2 = aCand->userFloat("vProb");
   chamSegR = aCand->userFloat("chamSegdeltaR");
//  cout<<"aCand->mass()"<<aCand->mass()<<"muon1"<<muon1->pt()<<endl;
  float theCtau; 
  if (_useBS) {theCtau = 10.*aCand->userFloat("ppdlBS");}
  else {theCtau = 10.*aCand->userFloat("ppdlPV");}

  shareChamber = aCand->userInt("shareChamber");
  reco::GenParticleRef genUps = aCand->genParticleRef();
//  cout<<"genUps is available"<<aCand->genParticleRef()->pdgId()<<endl;
  if(genUps.isAvailable()) cout<<"genUps is available"<<endl;
  if(genUps.isAvailable() && (genUps->pdgId() ==100443 || genUps->pdgId() == 443)) cout<<"gen jpsi"<<endl;
  bool isMatched = (genUps.isAvailable() && (genUps->pdgId() == 553 || genUps->pdgId() ==100553 || genUps->pdgId() ==100443 || genUps->pdgId() == 443));
  if (isMatched && _removeSignal) return;
  // Trigger Results inspection
  
  static const unsigned int NTRIGGERS = 6;
  string HLTbitNames[NTRIGGERS] = {"HLT_L1MuOpen", "HLT_Mu3", "HLT_Mu5", "HLT_DoubleMu0", "HLT_DoubleMu3","HLT_L1DoubleMuOpen"};
  unsigned int hltBits[NTRIGGERS];
/*
  HLTConfigProvider hltConfig;
  cout<<"trigresult"<<hltConfig.init(_triggerresults.process())<<endl;
  if (hltConfig.init(_triggerresults.process())) {
    // check if trigger name in config
    const unsigned int n(hltConfig.size());
    for (unsigned int ihlt = 0; ihlt < NTRIGGERS; ihlt++) {
      hltBits[ihlt] = 0;
      unsigned int triggerIndex( hltConfig.triggerIndex(HLTbitNames[ihlt]) );
      if (triggerIndex>=n) {
	cout << "TriggerName " << HLTbitNames[ihlt] << " not available in config!" << endl;
      } else {
	hltBits[ihlt] = triggerIndex;
      }
    }
  }
*/
//  cout<<"HLT_L1DoubleMuOpen"<<trigger->accept(hltBits[5])<<endl;

//  cout<<"hlt_Mu3"<<trigger->accept(hltBits[0])<<endl;
  cout<<"HLT_Mu3"<<muon2->triggerObjectMatchesByPath("HLT_Mu3").empty()<<endl;
  cout<<"HLT_Mu5"<<muon2->triggerObjectMatchesByPath("HLT_Mu5").empty()<<endl;
  cout<<"HLT_DoubleMu0"<<muon2->triggerObjectMatchesByPath("HLT_DoubleMu0").empty()<<endl;
  cout<<"HLT_DoubleMu3"<<muon2->triggerObjectMatchesByPath("HLT_DoubleMu3").empty()<<endl;
  cout<<"HLT_L1DoubleMuOpen"<<muon2->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()<<endl;

// use HLT object matching instead of trigger results
    if ( (!muon2->triggerObjectMatchesByPath(HLTbitNames[0]).empty() || !muon1->triggerObjectMatchesByPath(HLTbitNames[0]).empty())) hlt_L1MuOpen = 1;
    else hlt_L1MuOpen = 0;
    if ((!muon2->triggerObjectMatchesByPath(HLTbitNames[1]).empty() || !muon1->triggerObjectMatchesByPath(HLTbitNames[1]).empty())) hlt_Mu3 = 1;
    else hlt_Mu3 = 0;
    if ((!muon2->triggerObjectMatchesByPath(HLTbitNames[2]).empty() || !muon1->triggerObjectMatchesByPath(HLTbitNames[2]).empty())) hlt_Mu5 = 1;
    else hlt_Mu5 = 0;
    if ((!muon2->triggerObjectMatchesByPath(HLTbitNames[3]).empty() && !muon1->triggerObjectMatchesByPath(HLTbitNames[3]).empty())) hlt_DoubleMu0 = 1;
    else hlt_DoubleMu0 = 0;
    if ((!muon2->triggerObjectMatchesByPath(HLTbitNames[4]).empty() && !muon1->triggerObjectMatchesByPath(HLTbitNames[4]).empty())) hlt_DoubleMu3 = 1;
    else hlt_DoubleMu3 = 0;
//    if(trigger->accept(hltBits[5])) hlt_L1DoubleMuOpen = 1;
//    else hlt_L1DoubleMuOpen = 0;
    if ((!muon2->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty() && !muon1->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty())) hlt_L1DoubleMuOpen = 1;
   else hlt_L1DoubleMuOpen = 0;

  reco::GenParticleRef genMu1 = muon1->genParticleRef();
  reco::GenParticleRef genMu2 = muon2->genParticleRef();
	
  // Signal / background Upsilon 	
  cout<<"muon1->charge()"<<muon1->charge()<<"muon2->charge()"<<muon2->charge()<<endl;
  if (PFDVtxs.isValid()) {
      cout<<"PFDVtxs.isValid()"<<endl;
 
       FromKLoose = 0;
       for(PFDisplacedVertexCollection::const_iterator vtx = PFDVtxs->begin(); vtx != PFDVtxs->end(); ++vtx){
        if(vtx->isKplus_Loose()||vtx->isKminus_Loose()){
             cout<<"KLoose_True"<<endl;
             std::vector<reco::Track> refittedTracks = (*vtx).refittedTracks();
            for(unsigned it = 0; it < refittedTracks.size(); it++){
              cout<<"refittedTracks.size()"<<refittedTracks.size()<<endl;
              cout << "refitted track pt = " <<endl;
              cout<<"is"<<refittedTracks[it].pt()<<endl;
              cout<<"eta"<<refittedTracks[it].eta()<<endl;
              cout<<"phi"<<refittedTracks[it].phi()<<endl;

              TrackRef tk1 =  muon1->get<TrackRef>();
              TrackRef tk2 =  muon2->get<TrackRef>();
//              cout<<"vtx info px"<<(**it).px()<<endl;
//               cout<<","<<(*it)->eta()<<","<<(*it)->phi()<<endl;
              cout<<"mu1 info"<<tk1->pt()<<","<<tk1->eta()<<","<<tk1->phi()<<endl;
              cout<<"mu2 info"<<tk2->pt()<<","<<tk2->eta()<<","<<tk2->phi()<<endl;
                 if((fabs(refittedTracks[it].pt() - tk1->pt()  ) < 1e-4 &&
                    fabs(refittedTracks[it].eta() - tk1->eta() ) < 1e-4 &&
                    fabs(refittedTracks[it].phi() - tk1->phi() ) < 1e-4) ||
                    (fabs(refittedTracks[it].pt()  - tk2->pt()  ) < 1e-4 &&
                    fabs(refittedTracks[it].eta() - tk2->eta() ) < 1e-4 &&
                    fabs(refittedTracks[it].phi() - tk2->phi() ) < 1e-4)){
                     cout<<"KLoose_True match mu"<<endl;
                     FromKLoose = 1;
                 }
          }
//          }
        }//if
      }
  } 
  if( (OS && muon1->charge()*muon2->charge() < 0) ||
      (!OS &&  muon1->charge()*muon2->charge() >0) ){
   if(theMass < 20.0 && theMass > 0.0){
//   if(theMass < 12.0 && theMass > 7.0){
    if(muon1->charge() == -1){
       cout<<"muplus pt"<<muon1->pt()<<endl;
       muMinusPt = muon1->pt();
       muMinusP = muon1->p();
       muMinusEta = muon1->eta();
       muMinusPhi = muon1->phi();
       muPlusPt = muon2->pt();
       muPlusP = muon2->p();
       muPlusEta = muon2->eta();
       muPlusPhi = muon2->phi();
       muMinusD0 = muon1->innerTrack()->d0();
       muMinusDz = muon1->innerTrack()->dz();
       muMinusNhits = muon1->innerTrack()->found();
       muMinusNPixelhits = muon1->innerTrack()->hitPattern().numberOfValidPixelHits();
       muMinusNSeg = muon1->numberOfMatches(Muon::SegmentAndTrackArbitration);
       muPlusD0 = muon2->innerTrack()->d0();
       muPlusDz = muon2->innerTrack()->dz();
       muPlusNhits = muon2->innerTrack()->found();
       muPlusNSeg = muon2->numberOfMatches(Muon::SegmentAndTrackArbitration);
       muPlusNPixelhits = muon2->innerTrack()->hitPattern().numberOfValidPixelHits();
       cout<<"*************************"<<HLTbitNames[0]<<"="<<muon2->triggerObjectMatchesByPath(HLTbitNames[0]).empty()<<endl;
       cout<<"*************************"<<HLTbitNames[0]<<"="<<muon1->triggerObjectMatchesByPath(HLTbitNames[0]).empty()<<endl;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[0]).empty()){
          L1MuOpen_Plus = 0;
       }else   L1MuOpen_Plus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[1]).empty()){
          Mu3_Plus = 0;
       }else  Mu3_Plus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[2]).empty()){
          Mu5_Plus = 0;
       }else  Mu5_Plus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[3]).empty()){
          DoubleMu0_Plus = 0;
       }else  DoubleMu0_Plus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[4]).empty()){
         DoubleMu3_Plus = 0; 
       }else  DoubleMu3_Plus = 1;
       cout<<"*************************"<<"L1DoubleMuOpen_Plus="<<muon1->triggerObjectMatchesByPath("hltDoubleMuLevel1PathL1OpenFiltered").empty()<<endl;
       cout<<"*************************"<<"L1DoubleMuOpen_Plus="<<muon2->triggerObjectMatchesByPath("hltDoubleMuLevel1PathL1OpenFiltered").empty()<<endl;
       if(muon2->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()){
          L1DoubleMuOpen_Plus = 0;
       }else  L1DoubleMuOpen_Plus = 1;

       if(muon1->triggerObjectMatchesByPath(HLTbitNames[0]).empty()){
          L1MuOpen_Minus = 0;
       }else   L1MuOpen_Minus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[1]).empty()){
          Mu3_Minus = 0;
       }else  Mu3_Minus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[2]).empty()){
          Mu5_Minus = 0;
       }else  Mu5_Minus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[3]).empty()){
          DoubleMu0_Minus = 0;
       }else  DoubleMu0_Minus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[4]).empty()){
         DoubleMu3_Minus = 0;
       }else  DoubleMu3_Minus = 1;
       if(muon1->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()){
          L1DoubleMuOpen_Minus = 0;
       }else  L1DoubleMuOpen_Minus = 1;
//       muon->isIsolationValid() && muon->isolationR03().sumPt <= 5 && muon->isolationR03().emEt <= 4 && muon->isolationR03().hadEt <= 4 && muon->isolationR03().hoEt<= 1.5 && muon->isolationR03().nTracks <= 5 && muon->isolationR03().nJets <= 2
       if(muon1->isIsolationValid()){
         muMinusIso03NTracks = muon1->isolationR03().nTracks;
         muMinusIso03sumPt = muon1->isolationR03().sumPt;
         muMinusIso03emEt = muon1->isolationR03().emEt;
          muMinusIso03hadEt = muon1->isolationR03().hadEt;
         muMinusIso03hoEt = muon1->isolationR03().hoEt;
          muMinusIso03NJets = muon1->isolationR03().nJets;
          muMinusIso03trackerVetoPt = muon1->isolationR03().trackerVetoPt;
         muMinusIso03emVetoEt = muon1->isolationR03().emVetoEt;
         muMinusIso03hadVetoEt = muon1->isolationR03().hadVetoEt;
         muMinusIso03hoVetoEt = muon1->isolationR03().hoVetoEt;    
       }     
       if(muon2->isIsolationValid()){
         muPlusIso03NTracks = muon2->isolationR03().nTracks;
         muPlusIso03sumPt = muon2->isolationR03().sumPt;
         muPlusIso03emEt = muon2->isolationR03().emEt;
          muPlusIso03hadEt = muon2->isolationR03().hadEt;
         muPlusIso03hoEt = muon2->isolationR03().hoEt;
          muPlusIso03NJets = muon2->isolationR03().nJets;
          muPlusIso03trackerVetoPt = muon2->isolationR03().trackerVetoPt;
         muPlusIso03emVetoEt = muon2->isolationR03().emVetoEt;
         muPlusIso03hadVetoEt = muon2->isolationR03().hadVetoEt;
         muPlusIso03hoVetoEt = muon2->isolationR03().hoVetoEt;

       }

        for(int i = 0; i < muon2->innerTrack()->hitPattern().numberOfValidMuonCSCHits(); i++){
           muon2->outerTrack()->hitPattern().printHitPattern(i,std::cout);
        }
       if(muon1->isGlobalMuon()){
           cout<<"chi2"<<muon1->globalTrack()->chi2()<<"ndof"<<muon1->globalTrack()->ndof()<<endl;
          muMinusGlbNorChi2 = muon1->globalTrack()->chi2()/muon1->globalTrack()->ndof();
          muMinusGlbTkChi2 = muon1->innerTrack()->chi2()/muon1->innerTrack()->ndof();
          muMinusGlbNMuHits = muon1->globalTrack()->hitPattern().numberOfValidMuonHits();
       }else {
           cout<<"tkchi"<<muon1->innerTrack()->chi2()<<"ndof"<<muon1->innerTrack()->ndof();
          muMinusTkNorChi2 = muon1->innerTrack()->chi2()/muon1->innerTrack()->ndof();
       }
       if(muon2->isGlobalMuon()){
           cout<<"chi2+"<<muon2->globalTrack()->chi2()<<"ndof"<<muon2->globalTrack()->ndof()<<endl;
          muPlusGlbNorChi2 = muon2->globalTrack()->chi2()/muon2->globalTrack()->ndof();
          muPlusGlbTkChi2 = muon2->innerTrack()->chi2()/muon2->innerTrack()->ndof();
          muPlusGlbNMuHits = muon2->globalTrack()->hitPattern().numberOfValidMuonHits();
       }else {
           cout<<"tkchi2+"<<muon2->innerTrack()->chi2()<<"tkndof"<<muon2->innerTrack()->ndof()<<endl;
          muPlusTkNorChi2 = muon2->innerTrack()->chi2()/muon2->innerTrack()->ndof();
       }
//       if(genMu1.isAvailable()&&genMu1->pdgId() == 13 && genMu2.isAvailable()&&genMu2->pdgId() == -13){
//          cout<<"genref == 1"<<endl;
//          GenRef = 1;
//       }else GenRef = 0;
       if(genMu1.isAvailable()&&genMu1->pdgId() == 13){
         cout<<"mu1 pt"<<genMu1->pt()<<endl;
         genMinMuPt = genMu1->pt();
         genMinMuEta = genMu1->eta();
         genMinMuPhi = genMu1->phi();
       }
       if(genMu2.isAvailable()&&genMu2->pdgId() == -13){
       genPosMuPt = genMu2->pt();
       genPosMuEta = genMu2->eta();
       genPosMuPhi = genMu2->phi();
       }
          muMinusTrkArbirtated = ((muon1->isTrackerMuon() &&muon1->muonID("TrackerMuonArbitrated")) ? 1:0);
          muMinusTrkLSLoose = ((muon1->isTrackerMuon() &&muon1->muonID("TMLastStationLoose"))  ? 1:0);
          muMinusTrkLSTight = ((muon1->isTrackerMuon() &&muon1->muonID("TMLastStationTight")) ? 1:0);
          muMinusTrk2DLoose = ((muon1->isTrackerMuon() &&muon1->muonID("TM2DCompatibilityLoose"))? 1:0);
          muMinusTrk2DTight = ((muon1->isTrackerMuon() &&muon1->muonID("TM2DCompatibilityTight"))? 1:0);
          muMinusTrkOSLoose = ((muon1->isTrackerMuon() &&muon1->muonID("TMOneStationLoose"))? 1:0);
          muMinusTrkOSTight = ((muon1->isTrackerMuon() &&muon1->muonID("TMOneStationTight"))? 1:0);
          muMinusTrkLSLowPtLoose = ((muon1->isTrackerMuon() &&muon1->muonID("TMLastStationOptimizedLowPtLoose"))? 1:0);
          muMinusTrkLSLowPtTight = ((muon1->isTrackerMuon() &&muon1->muonID("TMLastStationOptimizedLowPtTight"))? 1:0);
          muMinusTrkLSAngLoose =((muon1->isTrackerMuon() && muon1->muonID("TMLastStationAngLoose"))? 1:0);
          muMinusTrkLSAngTight = ((muon1->isTrackerMuon() &&muon1->muonID("TMLastStationAngTight"))? 1:0);
          muMinusTrkOSAngLoose = ((muon1->isTrackerMuon() &&muon1->muonID("TMOneStationAngLoose"))? 1:0);
          muMinusTrkOSAngTight = ((muon1->isTrackerMuon() &&muon1->muonID("TMOneStationAngTight"))? 1:0);
          muPlusTrkArbirtated = ((muon2->isTrackerMuon() &&muon2->muonID("TrackerMuonArbitrated"))? 1:0);
          muPlusTrkLSLoose = ((muon2->isTrackerMuon() && muon2->muonID("TMLastStationLoose"))? 1:0);
          muPlusTrkLSTight = ((muon2->isTrackerMuon()&&muon2->muonID("TMLastStationTight"))? 1:0);
          muPlusTrk2DLoose = ((muon2->isTrackerMuon()&&muon2->muonID("TM2DCompatibilityLoose"))? 1:0);
          muPlusTrk2DTight = ((muon2->isTrackerMuon()&&muon2->muonID("TM2DCompatibilityTight"))? 1:0);
          muPlusTrkOSLoose = ((muon2->isTrackerMuon()&&muon2->muonID("TMOneStationLoose"))? 1:0);
          muPlusTrkOSTight = ((muon2->isTrackerMuon()&&muon2->muonID("TMOneStationTight"))? 1:0);
          muPlusTrkLSLowPtLoose = ((muon2->isTrackerMuon()&&muon2->muonID("TMLastStationOptimizedLowPtLoose"))? 1:0);
          muPlusTrkLSLowPtTight = ((muon2->isTrackerMuon()&&muon2->muonID("TMLastStationOptimizedLowPtTight"))? 1:0);
          muPlusTrkLSAngLoose = ((muon2->isTrackerMuon()&&muon2->muonID("TMLastStationAngLoose"))? 1:0);
          muPlusTrkLSAngTight = ((muon2->isTrackerMuon()&&muon2->muonID("TMLastStationAngTight"))? 1:0);
          muPlusTrkOSAngLoose = ((muon2->isTrackerMuon()&&muon2->muonID("TMOneStationAngLoose"))? 1:0);
          muPlusTrkOSAngTight = ((muon2->isTrackerMuon()&&muon2->muonID("TMOneStationAngTight"))? 1:0);

    }else{
       cout<<"fill"<<endl;
       muMinusPt = muon2->pt();
       muMinusP = muon2->p();
       cout<<"wrong??"<<endl;
       muMinusEta = muon2->eta();
       muMinusPhi = muon2->phi();
       muPlusPt = muon1->pt();
       muPlusP = muon1->p();
       muPlusEta = muon1->eta();
       muPlusPhi = muon1->phi();
       muMinusD0 = muon2->innerTrack()->d0();
       muMinusDz = muon2->innerTrack()->dz();
       muMinusNhits = muon2->innerTrack()->found();
       muMinusNSeg = muon2->numberOfMatches(Muon::SegmentAndTrackArbitration);
       muMinusNPixelhits = muon2->innerTrack()->hitPattern().numberOfValidPixelHits();
          muPlusTrkArbirtated = ((muon1->isTrackerMuon() &&muon1->muonID("TrackerMuonArbitrated"))? 1:0);
          muPlusTrkLSLoose = ((muon1->isTrackerMuon() && muon1->muonID("TMLastStationLoose"))? 1:0);
          muPlusTrkLSTight = ((muon1->isTrackerMuon()&&muon1->muonID("TMLastStationTight"))? 1:0);
          muPlusTrk2DLoose = ((muon1->isTrackerMuon()&&muon1->muonID("TM2DCompatibilityLoose"))? 1:0);
          muPlusTrk2DTight = ((muon1->isTrackerMuon()&&muon1->muonID("TM2DCompatibilityTight"))? 1:0);
          muPlusTrkOSLoose = ((muon1->isTrackerMuon()&&muon1->muonID("TMOneStationLoose"))? 1:0);
          muPlusTrkOSTight = ((muon1->isTrackerMuon()&&muon1->muonID("TMOneStationTight"))? 1:0);
          muPlusTrkLSLowPtLoose = ((muon1->isTrackerMuon()&&muon1->muonID("TMLastStationOptimizedLowPtLoose"))? 1:0);
          muPlusTrkLSLowPtTight = ((muon1->isTrackerMuon()&&muon1->muonID("TMLastStationOptimizedLowPtTight"))? 1:0);
          muPlusTrkLSAngLoose = ((muon1->isTrackerMuon()&&muon1->muonID("TMLastStationAngLoose"))? 1:0);
          muPlusTrkLSAngTight = ((muon1->isTrackerMuon()&&muon1->muonID("TMLastStationAngTight"))? 1:0);
          muPlusTrkOSAngLoose = ((muon1->isTrackerMuon()&&muon1->muonID("TMOneStationAngLoose"))? 1:0);
          muPlusTrkOSAngTight = ((muon1->isTrackerMuon()&&muon1->muonID("TMOneStationAngTight"))? 1:0);
          muMinusTrkArbirtated = ((muon2->isTrackerMuon() && muon2->muonID("TrackerMuonArbitrated"))? 1:0);
          muMinusTrkLSLoose = ((muon2->isTrackerMuon() &&muon2->muonID("TMLastStationLoose"))? 1:0);
          muMinusTrkLSTight = ((muon2->isTrackerMuon() &&muon2->muonID("TMLastStationTight"))? 1:0);
          muMinusTrk2DLoose = ((muon2->isTrackerMuon() &&muon2->muonID("TM2DCompatibilityLoose"))? 1:0);
          muMinusTrk2DTight = ((muon2->isTrackerMuon() &&muon2->muonID("TM2DCompatibilityTight"))? 1:0);
          muMinusTrkOSLoose = ((muon2->isTrackerMuon() &&muon2->muonID("TMOneStationLoose"))? 1:0);
          muMinusTrkOSTight = ((muon2->isTrackerMuon() &&muon2->muonID("TMOneStationTight"))? 1:0);
          muMinusTrkLSLowPtLoose = ((muon2->isTrackerMuon() &&muon2->muonID("TMLastStationOptimizedLowPtLoose"))? 1:0);
          muMinusTrkLSLowPtTight = ((muon2->isTrackerMuon() &&muon2->muonID("TMLastStationOptimizedLowPtTight"))? 1:0);
          muMinusTrkLSAngLoose = ((muon2->isTrackerMuon() &&muon2->muonID("TMLastStationAngLoose"))? 1:0);
          muMinusTrkLSAngTight = ((muon2->isTrackerMuon() &&muon2->muonID("TMLastStationAngTight"))? 1:0);
          muMinusTrkOSAngLoose = ((muon2->isTrackerMuon() &&muon2->muonID("TMOneStationAngLoose"))? 1:0);
          muMinusTrkOSAngTight = ((muon2->isTrackerMuon() &&muon2->muonID("TMOneStationAngTight"))? 1:0);
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[0]).empty()){
          L1MuOpen_Plus = 0;
       }else   L1MuOpen_Plus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[1]).empty()){
          Mu3_Plus = 0;
       }else  Mu3_Plus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[2]).empty()){
          Mu5_Plus = 0;
       }else  Mu5_Plus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[3]).empty()){
          DoubleMu0_Plus = 0;
       }else  DoubleMu0_Plus = 1;
       if(muon1->triggerObjectMatchesByPath(HLTbitNames[4]).empty()){
         DoubleMu3_Plus = 0;
       }else  DoubleMu3_Plus = 1;
       if(muon1->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()){
          L1DoubleMuOpen_Plus = 0;
       }else  L1DoubleMuOpen_Plus = 1;

       if(muon2->triggerObjectMatchesByPath(HLTbitNames[0]).empty()){
          L1MuOpen_Minus = 0;
       }else   L1MuOpen_Minus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[1]).empty()){
          Mu3_Minus = 0;
       }else  Mu3_Minus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[2]).empty()){
          Mu5_Minus = 0;
       }else  Mu5_Minus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[3]).empty()){
          DoubleMu0_Minus = 0;
       }else  DoubleMu0_Minus = 1;
       if(muon2->triggerObjectMatchesByPath(HLTbitNames[4]).empty()){
         DoubleMu3_Minus = 0;
       }else  DoubleMu3_Minus = 1;
       if(muon2->triggerObjectMatchesByFilter("hltDoubleMuLevel1PathL1OpenFiltered").empty()){
          L1DoubleMuOpen_Minus = 0;
       }else  L1DoubleMuOpen_Minus = 1;

       if(muon1->isIsolationValid()){
         muPlusIso03NTracks = muon1->isolationR03().nTracks;
         muPlusIso03sumPt = muon1->isolationR03().sumPt;
         muPlusIso03emEt = muon1->isolationR03().emEt;
          muPlusIso03hadEt = muon1->isolationR03().hadEt;
         muPlusIso03hoEt = muon1->isolationR03().hoEt;
          muPlusIso03NJets = muon1->isolationR03().nJets;
          muPlusIso03trackerVetoPt = muon1->isolationR03().trackerVetoPt;
         muPlusIso03emVetoEt = muon1->isolationR03().emVetoEt;
         muPlusIso03hadVetoEt = muon1->isolationR03().hadVetoEt;
         muPlusIso03hoVetoEt = muon1->isolationR03().hoVetoEt;
       }
       if(muon2->isIsolationValid()){
         muMinusIso03NTracks = muon2->isolationR03().nTracks;
         muMinusIso03sumPt = muon2->isolationR03().sumPt;
         muMinusIso03emEt = muon2->isolationR03().emEt;
          muMinusIso03hadEt = muon2->isolationR03().hadEt ;
         muMinusIso03hoEt = muon2->isolationR03().hoEt;
          muMinusIso03NJets = muon2->isolationR03().nJets;
          muMinusIso03trackerVetoPt = muon2->isolationR03().trackerVetoPt;
         muMinusIso03emVetoEt = muon2->isolationR03().emVetoEt;
         muMinusIso03hadVetoEt = muon2->isolationR03().hadVetoEt;
         muMinusIso03hoVetoEt = muon2->isolationR03().hoVetoEt;
       }
      if(muon2->isGlobalMuon()){
          muMinusGlbNorChi2  = muon2->globalTrack()->chi2()/muon2->globalTrack()->ndof();
          muMinusGlbNMuHits = muon2->globalTrack()->hitPattern().numberOfValidMuonHits();
          muMinusGlbTkChi2  = muon2->innerTrack()->chi2()/muon2->innerTrack()->ndof();
       }else {
          muMinusTkNorChi2  = muon2->innerTrack()->chi2()/muon2->innerTrack()->ndof();
        }
       muPlusD0 = muon1->innerTrack()->d0();
       muPlusDz = muon1->innerTrack()->dz();
       muPlusNhits = muon1->innerTrack()->found();
       muPlusNSeg = muon1->numberOfMatches(Muon::SegmentAndTrackArbitration);
       muPlusNPixelhits = muon1->innerTrack()->hitPattern().numberOfValidPixelHits();
        for(int i = 0; i < muon1->innerTrack()->hitPattern().numberOfValidMuonCSCHits(); i++){
           muon1->outerTrack()->hitPattern().printHitPattern(i,std::cout);
        }
       if(muon1->isGlobalMuon()){
          muPlusGlbNorChi2 = muon1->globalTrack()->chi2()/muon1->globalTrack()->ndof();
          muPlusGlbNMuHits = muon1->globalTrack()->hitPattern().numberOfValidMuonHits();
          muPlusGlbTkChi2 = muon1->innerTrack()->chi2()/muon1->innerTrack()->ndof();
       }else {
          muPlusTkNorChi2 = muon1->innerTrack()->chi2()/muon1->innerTrack()->ndof();
        }
//       cout<<"avalible"<<genMu2.isAvailable()<<"id"<<genMu2->pdgId()<<endl;
       if(genMu2.isAvailable()&&genMu2->pdgId() == 13){
       genMinMuPt = genMu2->pt();
       genMinMuEta = genMu2->eta();
       genMinMuPhi = genMu2->phi();
       cout<<"genMu2->pt()"<<genMu2->pt()<<endl;
       }
       if(genMu1.isAvailable()&&genMu1->pdgId() == -13){
       genPosMuPt = genMu1->pt();;
       genPosMuEta = genMu1->eta();
       genPosMuPhi = genMu1->phi();
       }
    }
       invariantMass = theMass;
       upsPt = aCand->pt();
       upsRapidity = aCand->rapidity();
       upsEta = aCand->eta();
       upsPhi = aCand->phi();
    if (isMatched) {
       GenRef = 1;
       cout<<genMu2.isAvailable()<<genMu2->pdgId()<<genMu1.isAvailable()<<genMu1->pdgId()<<endl;
       if(genMu2.isAvailable()&& abs(genMu2->pdgId()) == 13 && genMu1.isAvailable()&& abs(genMu1->pdgId()) == 13){
          McInvMass = (lorentzMomentum(genMu1)+lorentzMomentum(genMu2)).Mag();
          cout<<"lorentzMomentum(*genMu1)"<<McInvMass<<endl;
       }
          genUpsP = genUps->p();
       genUpsPt = genUps->pt();
       genUpsRap = genUps->rapidity();
          genUpsEta = genUps->eta();
          genUpsPhi = genUps->phi();
    }else GenRef = 0;
cout<<"before filling"<<endl;
    t->Fill();
   }//themass
  }//opposite charge
// }//theCat = 2
}
        
void UpsilonAnalyzerPAT::makeCuts() {
  if (collAll.isValid()) {
    cout<<"collAll.isValid"<<endl;
    for(vector<pat::CompositeCandidate>::const_iterator it=collAll->begin();
	it!=collAll->end();++it) {
//Remove sharing track
       const pat::CompositeCandidate* cand = &(*it);	
       // cout << "Now checking candidate of type " << theUpsilonCat << " with pt = " << cand->pt() << endl;
       const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
       const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));
       bool sameTrack = true;
       sameTrack &= fabs(muon1->innerTrack()->pt()  - muon2->innerTrack()->pt()  ) < 1e-4;
       sameTrack &= fabs(muon1->innerTrack()->eta() - muon2->innerTrack()->eta() ) < 1e-4;
       sameTrack &= fabs(muon1->innerTrack()->phi() - muon1->innerTrack()->phi() ) < 1e-4;
       sameTrack &= (muon1->innerTrack()->numberOfValidHits() - muon1->innerTrack()->numberOfValidHits()) == 0;
      if (_RemovePairsSharingTrackerTrack && sameTrack){
         cout<<"remove"<<endl;
         continue;
      }
      cout<<"sharechamber="<< cand->userInt("shareChamber")<<endl;
     if (_RemovePairsSharingChamber && cand->userInt("shareChamber") == 1){
         cout<<"remove sharing chamber"<<endl;
//         throw;
         continue;
      } 
      if ( (OS && muon1->charge()*muon2->charge() < 0)||
           (!OS &&  muon1->charge()*muon2->charge()>0) ) {	  
// For glb-glb,glb-trk:	  
        // global + global?
      if(!trackeronly){
	if (muon1->isGlobalMuon() && muon2->isGlobalMuon()) {
            cout<<"vertex probaility"<<cand->userFloat("vProb")<<endl;
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selGlobalMuon(muon2) &&
//                              fabs(cand->rapidity()) < 2.3 && //FIXME 
			      cand->userFloat("vProb") > 0.001 )){
	    _thePassedCats.push_back(0);  _thePassedCands.push_back(cand);
            continue;
	  }
	}
	
        // global + tracker? (x2)    
	if (muon1->isGlobalMuon() && muon2->isTrackerMuon() ) {
            cout<<"vertex probaility"<<cand->userFloat("vProb")<<endl;
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selTrackerMuon(muon2) 
  //                            fabs(cand->rapidity())<2.3 && //FIXME
			       && cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats.push_back(1);  _thePassedCands.push_back(cand);
	    continue;
	  }
	}

        if (muon2->isGlobalMuon() && muon1->isTrackerMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon2) &&
			      selTrackerMuon(muon1) 
    //                          fabs(cand->rapidity())<2.3 && //FIXME
			      && cand->userFloat("vProb") > 0.001
)){	    _thePassedCats.push_back(1);  _thePassedCands.push_back(cand);
	    continue;
	  }
	}
    }//if tracker only
        // tracker + tracker?  
        if (muon1->isTrackerMuon() && muon2->isTrackerMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon1) &&
			      selTrackerMuon(muon2) 
      //                        fabs(cand->rapidity())<2.3 &&//FIXME
			      && cand->userFloat("vProb") > 0.001
)){	    _thePassedCats.push_back(2);  _thePassedCands.push_back(cand);
	    continue;
	  }
	}
      }
    }
  }
  
  if (_useCalo && collCalo.isValid()) {

    for(vector<pat::CompositeCandidate>::const_iterator it=collCalo->begin();
	it!=collCalo->end();++it) {
      
      const pat::CompositeCandidate* cand = &(*it);
      
      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));
      
      if (muon1->charge()*muon2->charge() < 0 && !(muon1->isTrackerMuon()) && !(muon2->isTrackerMuon()) ) {

	// global + calo? (x2)
	if (muon1->isGlobalMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats.push_back(3);  _thePassedCands.push_back(cand);
            continue;
	  }
	}

	if (muon2->isGlobalMuon() && muon1->isCaloMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon2) &&
			      selCaloMuon(muon1) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats.push_back(3);  _thePassedCands.push_back(cand);
            continue;
	  }
	}
	
        // tracker + calo? (x2)    
	if (muon1->isTrackerMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats.push_back(4);  _thePassedCands.push_back(cand);
	    continue;
	  }
	}

        if (muon2->isTrackerMuon() && muon1->isCaloMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon2) &&
			      selCaloMuon(muon1) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats.push_back(4);  _thePassedCands.push_back(cand);
	    continue;
	  }
	}

        // calo + calo? 
        if (muon1->isCaloMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selCaloMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats.push_back(5);  _thePassedCands.push_back(cand);
	    continue;
	  }
	}
      }
    }
  }

  return;
}

pair< unsigned int, const pat::CompositeCandidate* > 
UpsilonAnalyzerPAT::theBestQQ() {

  unsigned int theBestCat = 99;
  const pat::CompositeCandidate* theBestCand = new pat::CompositeCandidate();

  for( unsigned int i = 0; i < _thePassedCands.size(); i++) { 
    if (_thePassedCats.at(i) < theBestCat) {
      theBestCat = _thePassedCats.at(i);
      theBestCand = _thePassedCands.at(i);
    }
  }

  pair< unsigned int, const pat::CompositeCandidate* > result = make_pair(theBestCat, theBestCand );
  return result;

}

bool
UpsilonAnalyzerPAT::selGlobalMuon(const pat::Muon* aMuon) {

  TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::Vertex *thePrimaryV = & (*priVtxs)[0];

  return (
//   aMuon->muonID("GlobalMuonPromptTight")&&//);
          iTrack->found() > 13 &&//11
          aMuon->muonID("TMLastStationAngTight")&&
          aMuon->muonID("TrackerMuonArbitrated")&&
	  aMuon->globalTrack()->chi2()/aMuon->globalTrack()->ndof() < 20.0 &&
	  (p.numberOfValidPixelHits() >1) && 
	  fabs(iTrack->dxy(thePrimaryV->position())) < 3.0 && //5.0 &&
          fabs(iTrack->dz(thePrimaryV->position())) < 20.0 &&
          iTrack->p() >= 2.5 && iTrack->pt() >= 1 );
//          ((fabs(aMuon->eta()) < 1.2 && fabs(aMuon->pt()) > 3.5) ||
//          (fabs(aMuon->eta()) > 1.2 && fabs(aMuon->eta()) < 2.4 && fabs(aMuon->pt()) > 3.0))); //FIXME
//)

}

bool 
UpsilonAnalyzerPAT::selTrackerMuon(const pat::Muon* aMuon) {
  
  TrackRef iTrack = aMuon->innerTrack();
  cout<<"iTrack is valid"<<iTrack.isAvailable()<<endl;
  const reco::Vertex *thePrimaryV = & (*priVtxs)[0];
  if(iTrack.isAvailable()){ const reco::HitPattern& p = iTrack->hitPattern();
  return (
          iTrack.isAvailable() &&
          iTrack->found() > 11 && //11
	  iTrack->chi2()/iTrack->ndof() < 5.0 &&
	  aMuon->muonID("TMLastStationAngTight")&&
          aMuon->muonID("TrackerMuonArbitrated")&&
	  (p.numberOfValidPixelHits() > 1) && 
	  fabs(iTrack->dxy(thePrimaryV->position())) < 3.0 &&//5.0
          fabs(iTrack->dz(thePrimaryV->position())) < 20.0 );
   }else return 0;
//          iTrack->p() >= 2.5 && iTrack->pt() >= 1 );
//          ((fabs(aMuon->eta()) < 1.2 && fabs(aMuon->pt()) > 3.5) || 
//          (fabs(aMuon->eta()) > 1.2 && fabs(aMuon->eta()) < 2.4 && fabs(aMuon->pt()) > 3.0))
//);
}

bool 
UpsilonAnalyzerPAT::selCaloMuon(const pat::Muon* aMuon) {
//useless now  
  TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::Vertex *thePrimaryV = & (*priVtxs)[0];

  return (aMuon->caloCompatibility() > 0.89 &&
	  iTrack->found() > 13 &&
	  iTrack->chi2()/iTrack->ndof() < 5.0 &&
	  (p.numberOfValidPixelHits() > 1) && 
	  fabs(iTrack->d0()) < 5.0 &&
          fabs(iTrack->dz()) < 20.0 &&
          ((fabs(aMuon->eta()) < 1.2 && fabs(aMuon->pt()) > 3.5) || 
          (fabs(aMuon->eta()) > 1.2 && fabs(aMuon->eta()) < 2.4 && fabs(aMuon->pt()) > 3.0)));
}

TLorentzVector UpsilonAnalyzerPAT::lorentzMomentum(const GenParticleRef& genp) const {

    double ereco = genp->energy();
    double pxreco = genp->px();
    double pyreco = genp->py();
    double pzreco = genp->pz();

    // cout << "pfcl = " << pxreco << " "  << pyreco << " "  << pzreco << " "  << ereco << endl;
    TLorentzVector lrzpreco(pxreco, pyreco, pzreco, ereco);
    return lrzpreco;

}

//define this as a plug-in
DEFINE_FWK_MODULE(UpsilonAnalyzerPAT);
