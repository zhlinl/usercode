// -*- C++ -*-
//
// Package:    JPsiAnalyzerPAT
// Class:      JPsiAnalyzerPAT
// 
/**\class JPsiAnalyzerPAT JPsiAnalyzerPAT.cc OctoberXTracking/JPsiAnalyzerPAT/src/JPsiAnalyzerPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Roberto Covarelli 
//         Created:  Fri Oct  9 04:59:40 PDT 2009
// $Id: JPsiAnalyzerPAT.cc,v 1.41 2010/11/26 11:00:06 covarell Exp $
//
// based on: Onia2MuMu package V00-11-00
// changes done by: FT

// system include files
#include <memory>
#include <fstream>
#include <ostream>
#include <math.h>

// ROOT/Roofit include files
#include <TStyle.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/BeamSpot/interface/BeamSpot.h>

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace RooFit;

//
// class declaration
//
class JPsiAnalyzerPAT : public edm::EDAnalyzer {
   public:
      explicit JPsiAnalyzerPAT(const edm::ParameterSet&);
      ~JPsiAnalyzerPAT();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void makeCuts(int sign) ;
      pair< unsigned int, const pat::CompositeCandidate* > theBestQQ(int sign);
      void fillTreeAndDS(unsigned int theCat, const pat::CompositeCandidate* aCand, const edm::Event&);
      bool isMuonInAccept(const pat::Muon* aMuon);
      bool selGlobalMuon(const pat::Muon* aMuon);
      bool selTrackerMuon(const pat::Muon* aMuon);
      bool selCaloMuon(const pat::Muon* aMuon);
      int getJpsiVarType(const double jpsivar, vector<double> vectbin);
      double CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode);

      // additional functions by f
      void resetDSVariables();
      void analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles);
      // void calcPol(TLorentzVector&, TLorentzVector&, std::vector< float >&, std::vector< float >& );
      void beginRun(const edm::Run &, const edm::EventSetup &);
      void hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup);
      void matchMuonToHlt(const pat::Muon*, const pat::Muon*);

      // ROOT tree 
      TTree* tree_;//data; //*recoData;
      TFile* fOut_;

      // SMALL dataset and RooRealVars
      TFile* fOut2_;
      RooDataSet* data;
      RooRealVar* Jpsi_Mass;      
      RooRealVar* Jpsi_Pt;
      RooRealVar* Jpsi_Rap; 
      RooRealVar* Jpsi_ct;
      RooRealVar* Jpsi_ctErr;
      RooRealVar* Jpsi_ctTrue;      			
      RooCategory* Jpsi_Type;
      RooCategory* Jpsi_PtType;
      RooCategory* Jpsi_RapType;
      RooCategory* Jpsi_MatchType;
      RooCategory* Jpsi_Sign;   

      //1.) J/psi variables RECO
       double JpsiMass, JpsiPt, JpsiRap;
       double JpsiPx, JpsiPy, JpsiPz;
      TLorentzVector* JpsiP;
      double Jpsict, JpsictErr, JpsiVprob;
      int JpsiType,  JpsiCharge, MCType; //GG, GT and TT

      //2.) muon variables RECO
       double muPosPx, muPosPy, muPosPz;
      TLorentzVector* muPosP;
       double muNegPx, muNegPy, muNegPz;
      TLorentzVector* muNegP;
			//-----------------------------------------------------
			//-----------additional Reco Muon Variables------------
			//-----------------------------------------------------
			//(1). Positive Muon
			double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;
			int  muPos_arbitrated, muPos_oneStationTight, muPos_lastStationAngTight, 
					 muPos_lastStationTight, muPos_lastStationOptimizedLowPtTight,
					 muPos_lastStationOptimizedBarrelLowPtTight,muPos_oneStationAngTight;
			int muPos_found, muPos_pixeLayers, muPos_nValidMuHits;
			//(2).Negative Muon
			double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;
			int muNeg_arbitrated, muNeg_oneStationTight, muNeg_lastStationAngTight,
					muNeg_lastStationTight, muNeg_lastStationOptimizedLowPtTight,
					muNeg_lastStationOptimizedBarrelLowPtTight, muNeg_oneStationAngTight;
			int muNeg_found, muNeg_pixeLayers, muNeg_nValidMuHits;


      //3.) J/psi variables GEN
      // double JpsiMass_Gen, JpsiPt_Gen, JpsiRap_Gen;
      // double JpsiPx_Gen,   JpsiPy_Gen, JpsiPz_Gen;
      TLorentzVector* JpsiP_Gen;
      double Jpsict_Gen;

      //4.)muon variables GEN
      // double muPosPx_Gen, muPosPy_Gen, muPosPz_Gen;
      TLorentzVector* muPosP_Gen;
      // double muNegPx_Gen, muNegPy_Gen, muNegPz_Gen;
      TLorentzVector* muNegP_Gen;

      //5.) Event related variables
      unsigned int eventNb, runNb, lumiBlock, nPriVtx;

      //6.) POL variables
      // std::vector<std::string> polVarNames_;
      // std::vector<std::string> polVarNamesGen_;
      // std::map<std::string, double> mapPolVarsToValue_;
      // std::map<std::string, double> mapPolVarsToValueGen_;

      //7.) TriggerNames Map
      std::map<std::string, int> mapTriggerNameToIntFired_;

      Handle<pat::CompositeCandidateCollection > collAll;
      Handle<pat::CompositeCandidateCollection > collCalo;
      // Handle<TriggerResults> trigger;

      // data members
      InputTag       _patJpsi;
      InputTag       _patJpsiWithCalo;
      bool           _writeTree;
      string         _treefilename; 
      bool           _writeDataSet; 
      string         _datasetname;
      string         _triggerForDataset;
      double         _massMin;
      double         _massMax;
      vector<double> _ptbinranges;
      vector<double> _etabinranges;
      bool           _onlythebest;
      bool           _applycuts;
      bool           _useBS;
      bool           _useCalo;
      bool           _removeSignal;
      bool           _removeMuons;
      bool           _storeWs;
      bool           _writeOutCands;
      int            _MassCorr;
      // bool           _JSON;
      int            _oniaPDG;
      InputTag       _genParticles;
      bool           _isMC;
      bool           _isPromptMC;

      InputTag      _triggerresults;
      vector<unsigned int>                     _thePassedCats[3];
      vector<const pat::CompositeCandidate*>   _thePassedCands[3];

      // number of events
      unsigned int nEvents;
      unsigned int passedTriggerResults_;
      unsigned int passedMuonSelectionCuts_;
      unsigned int passedTriggerMatch_;
      unsigned int passedTriggerResultsAnalyzer_;


      // limits 
      float JpsiMassMin;
      float JpsiMassMax;
      float JpsiCtMin;
      float JpsiCtMax;
      float JpsiPtMin;           // SET BY 
      float JpsiPtMax;           // DEFINITION
      float JpsiRapMin;          // OF BIN
      float JpsiRapMax;          // LIMITS

      math::XYZPoint RefVtx;
      ofstream* theTextFile;
      ofstream* JSON;
  
      int runtmp,lumitmp,count;
      int runmax,runmin;

      // Trigger Filter Studies
      edm::Handle< edm::TriggerResults> handleTriggerResults_;
      edm::InputTag tagTriggerResults_;
      HLTConfigProvider hltConfig_;
      bool hltConfigInit_;
      std::vector<std::string> HLTbitNames_;
      std::map< std::string, std::string> mapTriggerToLastFilter_;
};    
      
// constants, enums and typedefs
enum {CS, HX, PHX, sGJ, GJ1, GJ2};
//

//
// constructors and destructor
//
JPsiAnalyzerPAT::JPsiAnalyzerPAT(const edm::ParameterSet& iConfig):
  _patJpsi(iConfig.getParameter<InputTag>("src")),
  _patJpsiWithCalo(iConfig.getParameter<InputTag>("srcWithCaloMuons")),
  _writeTree(iConfig.getParameter<bool>("writeTree")),
  _treefilename(iConfig.getParameter<string>("treeFileName")),	
  _writeDataSet(iConfig.getParameter<bool>("writeDataSet")),
  _datasetname(iConfig.getParameter<string>("dataSetName")),
  _triggerForDataset(iConfig.getParameter<string>("triggerForDataset")),
  _massMin(iConfig.getParameter<double>("massMin")),
  _massMax(iConfig.getParameter<double>("massMax")),
  _ptbinranges(iConfig.getParameter< vector<double> >("pTBinRanges")),	
  _etabinranges(iConfig.getParameter< vector<double> >("etaBinRanges")),
  _onlythebest(iConfig.getParameter<bool>("onlyTheBest")),		
  _applycuts(iConfig.getParameter<bool>("applyCuts")),			
  _useBS(iConfig.getParameter<bool>("useBeamSpot")),
  _useCalo(iConfig.getUntrackedParameter<bool>("useCaloMuons",false)),
  _removeSignal(iConfig.getUntrackedParameter<bool>("removeSignalEvents",false)),
  _removeMuons(iConfig.getUntrackedParameter<bool>("removeTrueMuons",false)),
  _storeWs(iConfig.getUntrackedParameter<bool>("storeWrongSign",false)),
  _writeOutCands(iConfig.getUntrackedParameter<bool>("writeOutCandidates",false)),
  _MassCorr(iConfig.getParameter<int>("massCorrectionMode")),
  _oniaPDG(iConfig.getParameter<int>("oniaPDG")),
  _genParticles(iConfig.getParameter<InputTag>("genParticles")),
  _isMC(iConfig.getUntrackedParameter<bool>("isMC",false)),
  _isPromptMC(iConfig.getUntrackedParameter<bool>("isPromptMC",false) ),
  tagTriggerResults_(iConfig.getParameter<InputTag>("TriggerResultsLabel"))
{
   //now do what ever initialization is needed
  nEvents = 0; 
  // passedTriggerResults_=0;
  passedMuonSelectionCuts_=0;
  passedTriggerMatch_=0;
  // passedTriggerResultsAnalyzer_=0;

  /* count=0;
  runtmp=0;
  lumitmp=11111111;
  runmax=0;
  runmin=100000000;*/

  JpsiMassMin = _massMin;
  JpsiMassMax = _massMax;
  JpsiCtMin = -1.0;
  JpsiCtMax = 3.5;

  if (_writeOutCands) theTextFile = new ofstream("passedCandidates.txt");
  
  /* if (_JSON){
    JSON = new ofstream("PseudoJSON.txt");
    *JSON << "{";
    }*/

  // Add the Trigger you want to choose into HLTbitNames_
  // Mu + Track Trigger
  HLTbitNames_.push_back( "HLT_Mu0_Track0_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu3_Track0_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu5_Track0_Jpsi" );
  //
  HLTbitNames_.push_back( "HLT_Mu0_TkMu0_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu3_TkMu0_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu5_TkMu0_Jpsi" );
  //
  HLTbitNames_.push_back( "HLT_Mu0_TkMu0_OST_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu3_TkMu0_OST_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu5_TkMu0_OST_Jpsi" );
  //
  HLTbitNames_.push_back( "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1" );
  HLTbitNames_.push_back( "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1" );
  HLTbitNames_.push_back( "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1" );
  //
  HLTbitNames_.push_back( "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2" );
  HLTbitNames_.push_back( "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2" );
  HLTbitNames_.push_back( "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2" );
  //
  HLTbitNames_.push_back( "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3" );
  HLTbitNames_.push_back( "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3" );
  HLTbitNames_.push_back( "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3" );
  //Mixed MuX + Track Trigger
  HLTbitNames_.push_back( "HLT_Mu3_Track3_Jpsi" );
  HLTbitNames_.push_back( "HLT_Mu3_Track3_Jpsi_v2" );
  HLTbitNames_.push_back( "HLT_Mu3_Track5_Jpsi_v1" );
  HLTbitNames_.push_back( "HLT_Mu3_Track5_Jpsi_v2" );
  // Double Muon Trigger
  HLTbitNames_.push_back( "HLT_DoubleMu0" );
  HLTbitNames_.push_back( "HLT_DoubleMu0_Quarkonium_v1" );
  HLTbitNames_.push_back( "HLT_DoubleMu0_Quarkonium_LS_v1" );
  HLTbitNames_.push_back( "HLT_L1DoubleMuOpen" );
  HLTbitNames_.push_back( "HLT_L1DoubleMuOpen_Tight" );
  HLTbitNames_.push_back( "HLT_DoubleMu3" );
  // Single Muon Trigger
  HLTbitNames_.push_back( "HLT_Mu3" );
  HLTbitNames_.push_back( "HLT_Mu5" );
  HLTbitNames_.push_back( "HLT_Mu7" );
  HLTbitNames_.push_back( "HLT_Mu9" );
  HLTbitNames_.push_back( "HLT_Mu11" );
  for(std::vector<std::string>::iterator it = HLTbitNames_.begin(); it != HLTbitNames_.end(); ++it){
      mapTriggerNameToIntFired_[*it] = -9999;
  }

// MAP your Trigger TO the last Filter used in the path
// MuX + Track Trigger
  mapTriggerToLastFilter_["HLT_Mu0_Track0_Jpsi"]="hltMu0TrackJpsiTrackMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu3_Track0_Jpsi"]="hltMu3TrackJpsiTrackMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu5_Track0_Jpsi"]="hltMu5TrackJpsiTrackMassFiltered";
//
  mapTriggerToLastFilter_["HLT_Mu0_TkMu0_Jpsi"] = "hltMu0TkMuJpsiTkMuMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu3_TkMu0_Jpsi"] = "hltMu3TkMuJpsiTkMuMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu5_TkMu0_Jpsi"] = "hltMu5TkMuJpsiTkMuMassFiltered";
//
  mapTriggerToLastFilter_["HLT_Mu0_TkMu0_OST_Jpsi"] = "hltMu0TkMuJpsiTkMuMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu3_TkMu0_OST_Jpsi"] = "hltMu3TkMuJpsiTkMuMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu5_TkMu0_OST_Jpsi"] = "hltMu5TkMuJpsiTkMuMassFiltered";
//
  mapTriggerToLastFilter_["HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1"] = "hltMu0TkMuJpsiTkMuMassFilteredTight";
  mapTriggerToLastFilter_["HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1"] = "hltMu3TkMuJpsiTkMuMassFilteredTight";
  mapTriggerToLastFilter_["HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1"] = "hltMu5TkMuJpsiTkMuMassFilteredTight";
//
  mapTriggerToLastFilter_["HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2"] = "hltMu0TkMuJpsiTkMuMassFilteredTight";
  mapTriggerToLastFilter_["HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2"] = "hltMu3TkMuJpsiTkMuMassFilteredTight";
  mapTriggerToLastFilter_["HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2"] = "hltMu5TkMuJpsiTkMuMassFilteredTight";
//
  mapTriggerToLastFilter_["HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3"] = "hltMu0TkMuJpsiTkMuMassFilteredTight";
  mapTriggerToLastFilter_["HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3"] = "hltMu3TkMuJpsiTkMuMassFilteredTight";
  mapTriggerToLastFilter_["HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3"] = "hltMu5TkMuJpsiTkMuMassFilteredTight";
//
  mapTriggerToLastFilter_["HLT_Mu3_Track3_Jpsi"]    = "hltMu3Track3JpsiTrackMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu3_Track3_Jpsi_v2"] = "hltMu3Track3JpsiTrackMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu3_Track5_Jpsi_v1"] = "hltMu3Track5JpsiTrackMassFiltered";
  mapTriggerToLastFilter_["HLT_Mu3_Track5_Jpsi_v2"] = "hltMu3Track5JpsiTrackMassFiltered";
// Double Muon Trigger
  mapTriggerToLastFilter_["HLT_DoubleMu0"]                  = "hltDiMuonL3PreFiltered0";
  mapTriggerToLastFilter_["HLT_DoubleMu0_Quarkonium_v1"]    = "hltDoubleMu0QuarkoniumL3PreFiltered";
  mapTriggerToLastFilter_["HLT_DoubleMu0_Quarkonium_LS_v1"] = "hltDoubleMu0QuarkoniumLSL3PreFiltered";
  mapTriggerToLastFilter_["HLT_L1DoubleMuOpen"]             = "hltDoubleMuLevel1PathL1OpenFiltered";
  mapTriggerToLastFilter_["HLT_L1DoubleMuOpen_Tight"]       = "hltL1DoubleMuOpenTightL1Filtered";
  mapTriggerToLastFilter_["HLT_DoubleMu3"]                  = "hltDiMuonL3PreFiltered";
// Single Muon Trigger
  mapTriggerToLastFilter_["HLT_Mu3"]            = "hltSingleMu3L3Filtered3";
  mapTriggerToLastFilter_["HLT_Mu5"]            = "hltSingleMu5L3Filtered5";
  mapTriggerToLastFilter_["HLT_Mu7"]            = "hltSingleMu7L3Filtered7";
  mapTriggerToLastFilter_["HLT_Mu9"]            = "hltSingleMu9L3Filtered9";
  mapTriggerToLastFilter_["HLT_Mu11"]           = "hltSingleMu11L3Filtered11";
}


JPsiAnalyzerPAT::~JPsiAnalyzerPAT()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (_writeOutCands) theTextFile->close();
}

// ------------ method called once each job just before starting event loop  ------------
void
JPsiAnalyzerPAT::beginJob()
{
    //std::cout << "[JPsiAnalyzerPAT] --- beginJob " << std::endl;
  if (_writeTree) {
    fOut_ = new TFile(_treefilename.c_str(), "RECREATE");
    fOut_->cd();

    // TTree
    //load Branches
    tree_ = new TTree ("data", "CMSSW Quarkonia J/psi Polarization+Trigger Tree");

    JpsiP = new TLorentzVector();
    muPosP = new TLorentzVector();
    muNegP = new TLorentzVector();
    JpsiP_Gen = new TLorentzVector();
    muPosP_Gen = new TLorentzVector();
    muNegP_Gen = new TLorentzVector();

    // Event variables
    tree_->Branch("eventNb",             &eventNb,             "eventNb/I");
    tree_->Branch("runNb",               &runNb,               "runNb/I");
    tree_->Branch("lumiBlock",           &lumiBlock,           "lumiBlock/I");
    tree_->Branch("nPriVtx",             &nPriVtx,             "nPriVtx/I");

    // Jpsi Variables
    tree_->Branch("JpsiType",   &JpsiType,  "JpsiType/I");
    tree_->Branch("JpsiP",  "TLorentzVector", &JpsiP);
     tree_->Branch("JpsiMass",   &JpsiMass,  "JpsiMass/D");
     tree_->Branch("JpsiPt",     &JpsiPt,    "JpsiPt/D");
     tree_->Branch("JpsiRap",    &JpsiRap,   "JpsiRap/D");
    tree_->Branch("JpsiCharge", &JpsiCharge,"JpsiCharge/I");
     tree_->Branch("JpsiPx",     &JpsiPx,    "JpsiPx/D");
     tree_->Branch("JpsiPy",     &JpsiPy,    "JpsiPy/D");
     tree_->Branch("JpsiPz",     &JpsiPz,    "JpsiPz/D");
    tree_->Branch("Jpsict",     &Jpsict,    "Jpsict/D");
    tree_->Branch("JpsictErr",  &JpsictErr, "JpsictErr/D");
    tree_->Branch("JpsiVprob",  &JpsiVprob, "JpsiVprob/D");
    tree_->Branch("muPosP", "TLorentzVector", &muPosP);
    tree_->Branch("muNegP", "TLorentzVector", &muNegP);
     tree_->Branch("muPosPx",    &muPosPx,   "muPosPx/D");
     tree_->Branch("muPosPy",    &muPosPy,   "muPosPy/D");
     tree_->Branch("muPosPz",    &muPosPz,   "muPosPz/D");
     tree_->Branch("muNegPx",    &muNegPx,   "muNegPx/D");
     tree_->Branch("muNegPy",    &muNegPy,   "muNegPy/D");
     tree_->Branch("muNegPz",    &muNegPz,   "muNegPz/D");

		//-----------------------------------------------------
		//-----------additional Reco Muon Variables------------
		//-----------------------------------------------------
		 //1). Positive Muon
		 tree_->Branch("muPos_nchi2In", &muPos_nchi2In, "muPos_nchi2In/D");
		 tree_->Branch("muPos_dxy", &muPos_dxy, "muPos_dxy/D");
		 tree_->Branch("muPos_dz", &muPos_dz, "muPos_dz/D");
		 tree_->Branch("muPos_nchi2Gl", &muPos_nchi2Gl, "muPos_nchi2Gl/D");
		 tree_->Branch("muPos_arbitrated", &muPos_arbitrated, "muPos_arbitrated/I");
		 tree_->Branch("muPos_found", &muPos_found, "muPos_found/I");
		 tree_->Branch("muPos_pixeLayers", &muPos_pixeLayers, "muPos_pixeLayers/I");
		 tree_->Branch("muPos_nValidMuHits", &muPos_nValidMuHits, "muPos_nValidMuHits/I");
		 tree_->Branch("muPos_oneStationTight", &muPos_oneStationTight, "muPos_oneStationTight/I");
		 tree_->Branch("muPos_lastStationAngTight",&muPos_lastStationAngTight, "muPos_lastStationAngTight/I");
		 tree_->Branch("muPos_lastStationTight",&muPos_lastStationTight, "muPos_lastStationTight/I");
		 tree_->Branch("muPos_lastStationOptimizedLowPtTight",&muPos_lastStationOptimizedLowPtTight, "muPos_lastStationOptimizedLowPtTight/I");
		 tree_->Branch("muPos_lastStationOptimizedBarrelLowPtTight",&muPos_lastStationOptimizedBarrelLowPtTight, "muPos_lastStationOptimizedBarrelLowPtTight/I");
		 tree_->Branch("muPos_oneStationAngTight",&muPos_oneStationAngTight,"muPos_oneStationAngTight/I");
		 //2). Negative Muon
		 tree_->Branch("muNeg_nchi2In", &muNeg_nchi2In, "muNeg_nchi2In/D");
		 tree_->Branch("muNeg_dxy", &muNeg_dxy, "muNeg_dxy/D");
		 tree_->Branch("muNeg_dz", &muNeg_dz, "muNeg_dz/D");
		 tree_->Branch("muNeg_nchi2Gl", &muNeg_nchi2Gl, "muNeg_nchi2Gl/D");
		 tree_->Branch("muNeg_arbitrated", &muNeg_arbitrated, "muNeg_arbitrated/I");
		 tree_->Branch("muNeg_found", &muNeg_found, "muNeg_found/I");
		 tree_->Branch("muNeg_pixeLayers", &muNeg_pixeLayers, "muNeg_pixeLayers/I");
		 tree_->Branch("muNeg_nValidMuHits", &muNeg_nValidMuHits, "muNeg_nValidMuHits/I");
		 tree_->Branch("muNeg_oneStationTight", &muNeg_oneStationTight, "muNeg_oneStationTight/I");
		 tree_->Branch("muNeg_lastStationAngTight",&muNeg_lastStationAngTight, "muNeg_lastStationAngTight/I");
		 tree_->Branch("muNeg_lastStationTight",&muNeg_lastStationTight, "muNeg_lastStationTight/I");
		 tree_->Branch("muNeg_lastStationOptimizedLowPtTight",&muNeg_lastStationOptimizedLowPtTight, "muNeg_lastStationOptimizedLowPtTight/I");
		 tree_->Branch("muNeg_lastStationOptimizedBarrelLowPtTight",&muNeg_lastStationOptimizedBarrelLowPtTight, "muNeg_lastStationOptimizedBarrelLowPtTight/I");
		 tree_->Branch("muNeg_oneStationAngTight",&muNeg_oneStationAngTight,"muNeg_oneStationAngTight/I");


    //add HLT Variables to TTree
    for(std::vector< std::string >:: iterator it = HLTbitNames_.begin(); it != HLTbitNames_.end(); ++it){
        std::string hlt_name= *it;
        tree_->Branch(hlt_name.c_str(), &(mapTriggerNameToIntFired_[*it]), (hlt_name + "/I").c_str());
    }

    //add Generator Information
    if(_isMC){
        tree_->Branch("MCType",         &MCType,        "MCType/I");
	tree_->Branch("JpsiP_Gen",  "TLorentzVector", &JpsiP_Gen);
        // tree_->Branch("JpsiMass_Gen",   &JpsiMass_Gen,  "JpsiMass_Gen/D");
        // tree_->Branch("JpsiPt_Gen",     &JpsiPt_Gen,    "JpsiPt_Gen/D");
        // tree_->Branch("JpsiRap_Gen",    &JpsiRap_Gen,   "JpsiRap_Gen/D");
        // tree_->Branch("JpsiPx_Gen",     &JpsiPx_Gen,    "JpsiPx_Gen/D");
        // tree_->Branch("JpsiPy_Gen",     &JpsiPy_Gen,    "JpsiPy_Gen/D");
        // tree_->Branch("JpsiPz_Gen",     &JpsiPz_Gen,    "JpsiPz_Gen/D");
        tree_->Branch("Jpsict_Gen",     &Jpsict_Gen,    "Jpsict_Gen/D");
        tree_->Branch("muPosP_Gen",  "TLorentzVector", &muPosP_Gen);
        tree_->Branch("muNegP_Gen",  "TLorentzVector", &muNegP_Gen);
        // tree_->Branch("muPosPx_Gen",    &muPosPx_Gen,   "muPosPx_Gen/D");
        // tree_->Branch("muPosPy_Gen",    &muPosPy_Gen,   "muPosPy_Gen/D");
        // tree_->Branch("muPosPz_Gen",    &muPosPz_Gen,   "muPosPz_Gen/D");
        // tree_->Branch("muNegPx_Gen",    &muNegPx_Gen,   "muNegPx_Gen/D");
        // tree_->Branch("muNegPy_Gen",    &muNegPy_Gen,   "muNegPy_Gen/D");
        // tree_->Branch("muNegPz_Gen",    &muNegPz_Gen,   "muNegPz_Gen/D");
    }
  }
  if (_writeDataSet) {
    
     fOut2_ = new TFile(_datasetname.c_str(), "RECREATE");
     fOut2_->cd();

     Jpsi_PtType = new RooCategory("Jpsi_PtType","Category of Pt");
     Jpsi_RapType = new RooCategory("Jpsi_RapType","Category of Rap");

     JpsiPtMin = _ptbinranges[0];  cout << "Pt min = " << JpsiPtMin << endl;
     JpsiPtMax = _ptbinranges[_ptbinranges.size()-1];  cout << "Pt max = " << JpsiPtMax << endl;
     
     for(unsigned int i=0;i<_ptbinranges.size()-1;i++){
       char catname[100];
       sprintf(catname,"P%d",i+1);
       Jpsi_PtType->defineType(catname,i+1); 
       cout << "Pt bin " << i+1 << ": Min = " << _ptbinranges[i] << " Max = " << _ptbinranges[i+1] << endl;   
     }
     
     JpsiRapMin = _etabinranges[0];  cout << "Rap min = " << JpsiRapMin << endl;
     JpsiRapMax = _etabinranges[_etabinranges.size()-1];  cout << "Rap max = " << JpsiRapMax << endl;
     
     for(unsigned int i=0;i<_etabinranges.size()-1;i++){
       char catname[100];
       sprintf(catname,"E%d",i+1);
       Jpsi_RapType->defineType(catname,i+1); 
       cout << "Rap bin " << i+1 << ": Min = " << _etabinranges[i] << " Max = " << _etabinranges[i+1] << endl;   
     }
     
     Jpsi_Type = new RooCategory("Jpsi_Type","Category of Jpsi");
     Jpsi_MatchType = new RooCategory("Jpsi_MatchType","Category of matching");
     
     Jpsi_Type->defineType("GG",0);
     Jpsi_Type->defineType("GT",1);
     Jpsi_Type->defineType("TT",2);
     if (_useCalo) {
       Jpsi_Type->defineType("GC",3);  
       Jpsi_Type->defineType("TC",4);  
       Jpsi_Type->defineType("CC",5);  
     }
     
     Jpsi_MatchType->defineType("unmatched",0);
     Jpsi_MatchType->defineType("matched",1);
     
     Jpsi_Sign = new RooCategory("Jpsi_Sign","Sign of Jpsi dimuons");
     
     Jpsi_Sign->defineType("OS",0);
     Jpsi_Sign->defineType("SSP",1);
     Jpsi_Sign->defineType("SSM",2);

     Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
     Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/psi pt",JpsiPtMin,JpsiPtMax,"GeV/c");
     Jpsi_Rap = new RooRealVar("Jpsi_Rap","J/psi eta",-JpsiRapMax,JpsiRapMax);
     Jpsi_ct = new RooRealVar("Jpsi_ct","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
     Jpsi_ctErr = new RooRealVar("Jpsi_ctErr","J/psi ctau error",-1.,1.,"mm");
     Jpsi_ctTrue = new RooRealVar("Jpsi_ctTrue","J/psi ctau true",-100.,JpsiCtMax,"mm"); 		

     RooArgList varlist(*Jpsi_Mass,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap,*Jpsi_Type,*Jpsi_MatchType);
     varlist.add(*Jpsi_ctTrue);   varlist.add(*Jpsi_PtType);
     varlist.add(*Jpsi_RapType);  varlist.add(*Jpsi_ctErr);
     varlist.add(*Jpsi_Sign);

     data = new RooDataSet("data","A sample",varlist);
  }
}


double JPsiAnalyzerPAT::CorrectMass(const reco::Muon& mu1,const reco::Muon& mu2, int mode){  
  double CMass=0;
  const double mumass=0.105658;
  double k1,k2;
  double pt1=mu1.innerTrack()->pt();
  double pt2=mu2.innerTrack()->pt();
  double eta1=mu1.innerTrack()->eta();
  double eta2=mu2.innerTrack()->eta();
  if (mode==1){
    k1=1.0009;//constant scale correction
    k2=1.0009;
  }
  if (mode==2){
    k1=1.0019-0.0004*pt1;
    k2=1.0019-0.0004*pt2; // pt dependent correction
  }
  if (mode==3){
      double a0=0.00038; //3.8 * pow(10,-4);
      double a1=0.0;
      double a2=0.0003; //3.0 * pow(10,-4);
      double a3=0.0;

      k1=1+a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=1+a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  if (mode == 4){
      double a0=1.002;
      double a1=-0.002;
      double a2=0.001;
      double a3=-0.0001;

      k1=a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  math::XYZVector mom1=mu1.innerTrack()->momentum();
  math::XYZVector mom2=mu2.innerTrack()->momentum();
  mom1=k1*mom1; 
  mom2=k2*mom2;
  double E1=sqrt(mom1.mag2()+(mumass*mumass));
  double E2=sqrt(mom2.mag2()+(mumass*mumass));
  math::XYZVector momtot=mom1+mom2;
  CMass=sqrt((E1+E2)*(E1+E2)-momtot.mag2());
  return CMass;
}

// ------------ method called to for each event  ------------
void
JPsiAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   nEvents++;

   // reset TTree Variables
   if (_writeTree) this->resetDSVariables();

   // check HLT TriggerReuslts
   this->hltReport(iEvent, iSetup);

   // Event related infos
   eventNb= iEvent.id().event() ;
   runNb=iEvent.id().run() ;
   lumiBlock= iEvent.luminosityBlock() ;

   /* int Nrun=iEvent.id().run();
   int lumi=iEvent.luminosityBlock();
   if (Nrun>runmax) runmax=Nrun;     // this is only for printout at the end of job
   if (Nrun<runmin) runmin=Nrun;
   if (_JSON){                      // these lines write out a JSON file of runs analyzed
     if (Nrun!=runtmp){
       runtmp=Nrun; 
       if (count == 0)  *JSON << " \"" << Nrun << "\" :[[" << lumi <<", ";
       if (count!=0) *JSON << lumitmp <<"]],"<< " \"" << Nrun << "\" :[[" << lumi <<", ";
       lumitmp=lumi;
       count++;
     }
     if (lumi!=lumitmp) {
       if (Nrun==runtmp && lumi!=lumitmp+1) {
	 *JSON << lumitmp << "],["<< lumi << ", "; 
       }
       lumitmp=lumi;
     }
     }*/

   Handle<reco::VertexCollection> privtxs;
   iEvent.getByLabel("offlinePrimaryVertices", privtxs);
   nPriVtx = privtxs->size();
   VertexCollection::const_iterator privtx;

   if ( privtxs->begin() != privtxs->end() ) {
     privtx=privtxs->begin();
     RefVtx = privtx->position();
   } else {
     RefVtx.SetXYZ(0.,0.,0.);
   }
   // }

   try {iEvent.getByLabel(_patJpsi,collAll);} 
   catch (...) {cout << "J/psi not present in event!" << endl;}

   if (_useCalo) {
     try {iEvent.getByLabel(_patJpsiWithCalo,collCalo);} 
     catch (...) {cout << "J/psi to calomuons not present in event!" << endl;}
   }

   _thePassedCats[0].clear();      _thePassedCands[0].clear();
   _thePassedCats[1].clear();      _thePassedCands[1].clear();
   _thePassedCats[2].clear();      _thePassedCands[2].clear();

   // APPLY CUTS
   int lastSign = 0;
   this->makeCuts(0);
   if (_storeWs) {
     this->makeCuts(1);
     this->makeCuts(2);
     lastSign = 2;
   }

   // BEST J/PSI? 
if (_onlythebest)
{  // yes, fill simply the best (possibly wrong-sign)

	for (int iSign = 0; iSign <= lastSign; iSign++) {
		pair< unsigned int, const pat::CompositeCandidate* > theBest = theBestQQ(iSign);
		if (theBest.first < 10)
		{
			fillTreeAndDS(theBest.first, theBest.second, iEvent);
			passedMuonSelectionCuts_++;

			//! FILL GENERATOR COLLECTION
			Handle<reco::GenParticleCollection> genParticles;
			iEvent.getByLabel( _genParticles, genParticles );
			if ( genParticles.isValid() )
			{
				//std::cout << "------ analyze GENERATED JPsis:" << std::endl;
				this->analyzeGenerator( genParticles );
			}

			// Write all Branches to the Tree ONLY 
			// - for the best candidate
			// - for the opposite sign
			if (iSign == 0 && _writeTree) tree_->Fill();
		}
	}

} 
else
{   // no, fill all candidates passing cuts (possibly wrong-sign)

	for (int iSign = 0; iSign <= lastSign; iSign++)
	{
		for( unsigned int count = 0; count < _thePassedCands[iSign].size(); count++)
		{ 
			fillTreeAndDS(_thePassedCats[iSign].at(count), _thePassedCands[iSign].at(count),iEvent); 
		}
	}

}
}

void
JPsiAnalyzerPAT::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

    //init HLTConfigProvider
    const std::string pro = tagTriggerResults_.process();
    bool changed = true;

    //bool init(const edm::Run& iRun, const edm::EventSetup& iSetup, const std::string& processName, bool& changed);
    hltConfigInit_ = false;
    if( hltConfig_.init(iRun, iSetup, pro, changed) ) hltConfigInit_ = true;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalyzerPAT::endJob() {
  
  cout << "Total number of events = " << nEvents << endl;
  // cout << "Analyzed runs from  " << runmin << "  to  " << runmax << endl; 
  cout << "============================================================" << endl;
  // cout << "Total number of passed candidates TRIGGER RESULTS ANALYZER   = " << passedTriggerResultsAnalyzer_ << endl;
  cout << "Total number of passed candidates MUON SELECTION CUTS        = " << passedMuonSelectionCuts_ << endl;
  cout << "Total number of passed candidates TRIGGER MATCH              = " << passedTriggerMatch_ << endl;

  /* if (_JSON){
    cout << "JSON file produced" << endl;
    *JSON << lumitmp <<"]]}";
    JSON->close();
    }*/

  // Write TTree to File
  if (_writeTree) {
    fOut_->cd();
    tree_->Write();
    fOut_->Close();
  }
  if (_writeDataSet) {
    fOut2_->cd();
    data->Write();
    fOut2_->Close();
  }
}

//! Fill the TTree with all RECO variables
void 
JPsiAnalyzerPAT::fillTreeAndDS(unsigned int theCat, const pat::CompositeCandidate* aCand, const edm::Event& iEvent){
  
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aCand->daughter("muon2"));
  
  //1.) continue only if we have opposite sign muons <- REMOVED
  //
  if (muon1->charge()*muon2->charge() >= 0) {
    if(muon1->charge() == 0) {
      printf("pat::Muon with zero charge!\n");   return;
    }
    if(muon2->charge() == 0) {
      printf("pat::Muon with zero charge!\n");   return;
    }
    // return;
  }

  const pat::Muon *muonPos = 0, *muonNeg = 0;
  if(muon1->charge() > 0){ muonPos = muon1; muonNeg = muon2;}
  else if(muon1->charge() < 0){ muonPos = muon2; muonNeg = muon1;}

  //
  float theMass = aCand->mass();
  if (_MassCorr!=0){
    double CMass = CorrectMass(*muon1,*muon2,_MassCorr);
    if (CMass!=0.0) theMass = CMass;
  }

  float theRapidity = aCand->rapidity();
  // if (!_useRapidity) theRapidity = theRapidity;

  float theCtau; 
  if (_useBS) {theCtau = 10.*aCand->userFloat("ppdlBS");}
  else {theCtau = 10.*aCand->userFloat("ppdlPV");}

  float theCtauErr; 
  if (_useBS) {theCtauErr = 10.*aCand->userFloat("ppdlErrBS");}
  else {theCtauErr = 10.*aCand->userFloat("ppdlErrPV");}

  float theCharge = aCand->charge();

  // MC matching
  reco::GenParticleRef genJpsi = aCand->genParticleRef();
  bool isMatched = (genJpsi.isAvailable() && genJpsi->pdgId() == _oniaPDG);

  // Input DataSet Type: P, NP, BG J/psi
  if (isMatched && _isPromptMC) MCType= 0;
  if (isMatched && _isPromptMC == false) MCType=1;
  if (!isMatched && _isMC) MCType=2;

  if (isMatched && _removeSignal) return;

  reco::GenParticleRef genMu1 = muon1->genParticleRef();
  reco::GenParticleRef genMu2 = muon2->genParticleRef();
  bool isMuMatched = (genMu1.isAvailable() && genMu2.isAvailable() && 
		      genMu1->pdgId()*genMu2->pdgId() == -169 && 
		      genMu1->momentum().rho() > 2.5 && genMu2->momentum().rho() > 2.5);
  if (isMuMatched && _removeMuons) return;

  // PAT trigger match, 2 muons to the last filter used in the HLT path (new way)
  this->matchMuonToHlt(muon1, muon2);

  // some counter
  passedTriggerMatch_++;
    
   JpsiMass=theMass;

  if (_writeOutCands) *theTextFile << iEvent.id().run() << "\t" << iEvent.luminosityBlock() << "\t" << iEvent.id().event() << "\t" << theMass << "\n";

  // write out JPsi RECO information
   JpsiPt=aCand->pt();
   JpsiRap=theRapidity;
  JpsiCharge=theCharge;
   //std::cout << "[JPsiAnalyzerPAT::fillTreeAndDS] ----- JpsiCharge: " << theCharge << std::endl;
   JpsiPx=aCand->px();
   JpsiPy=aCand->py();
   JpsiPz=aCand->pz();
  JpsiP->SetPxPyPzE(aCand->px(),aCand->py(),aCand->pz(),aCand->energy());
  Jpsict=theCtau;
  JpsictErr=theCtauErr;
  Jpsict_Gen=10.*aCand->userFloat("ppdlTrue");
  JpsiType=theCat;
  JpsiVprob=aCand->userFloat("vProb");
  
  // write out Muon RECO information
  float f_muPosPx, f_muPosPy, f_muPosPz;
  float f_muNegPx, f_muNegPy, f_muNegPz;
  f_muPosPx = muonPos->px();
  f_muPosPy = muonPos->py();
  f_muPosPz = muonPos->pz();
  f_muNegPx = muonNeg->px();
  f_muNegPy = muonNeg->py();
  f_muNegPz = muonNeg->pz();
   muPosPx= f_muPosPx ;
   muPosPy= f_muPosPy ;
   muPosPz= f_muPosPz ;
   muNegPx= f_muNegPx ;
   muNegPy= f_muNegPy ;
   muNegPz= f_muNegPz ;
  
  Double_t muMass = 0.105658;
  
  Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
   //TLorentzVector *muPosP = new TLorentzVector();
  muPosP->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);
  
  Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
   //TLorentzVector *muNegP = new TLorentzVector();
  muNegP->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);

	//-----------------------------------------------------
	//-----------additional Reco Muon Variables------------
	//-----------------------------------------------------
	//1.Positive Muon
	if(muonPos->isTrackerMuon())
	{
		TrackRef iTrack =muonPos->innerTrack();
		const reco::HitPattern& p1=iTrack->hitPattern();
		muPos_found=iTrack->found();
		muPos_nchi2In=iTrack->chi2()/iTrack->ndof();
		muPos_pixeLayers=p1.pixelLayersWithMeasurement();
		muPos_arbitrated=muonPos->muonID("TrackerMuonArbitrated");
		muPos_oneStationTight=muonPos->muonID("TMOneStationTight");
		muPos_lastStationAngTight=muonPos->muonID("TMLastStationAngTight");
		muPos_lastStationTight=muonPos->muonID("TMLastStationTight");
		muPos_lastStationOptimizedLowPtTight=muonPos->muonID("TMLastStationOptimizedLowPtTight");
		muPos_lastStationOptimizedBarrelLowPtTight=muonPos->muonID("TMLastStationOptimizedBarrelLowPtTight");
		muPos_oneStationAngTight=muonPos->muonID("TMOneStationAngTight");
		muPos_dxy=iTrack->dxy(RefVtx);
		muPos_dz=iTrack->dz(RefVtx);
		if(muonPos->isGlobalMuon())
		{
			TrackRef gTrack =muonPos->globalTrack();
			const reco::HitPattern& q1=gTrack->hitPattern();
			muPos_nValidMuHits=q1.numberOfValidMuonHits();
			muPos_nchi2Gl=gTrack->chi2()/gTrack->ndof();
		}
	}
  //2.Negative Muobn
	if(muonNeg->isTrackerMuon())
	{
		TrackRef iTrack =muonNeg->innerTrack();
		const reco::HitPattern& p2=iTrack->hitPattern();
		muNeg_found=iTrack->found();
		muNeg_nchi2In=iTrack->chi2()/iTrack->ndof();
		muNeg_pixeLayers=p2.pixelLayersWithMeasurement();
		muNeg_arbitrated=muonNeg->muonID("TrackerMuonArbitrated");
		muNeg_oneStationTight=muonNeg->muonID("TMOneStationTight");
		muNeg_lastStationAngTight=muonNeg->muonID("TMLastStationAngTight");
		muNeg_lastStationTight=muonNeg->muonID("TMLastStationTight");
		muNeg_lastStationOptimizedLowPtTight=muonNeg->muonID("TMLastStationOptimizedLowPtTight");
		muNeg_lastStationOptimizedBarrelLowPtTight=muonNeg->muonID("TMLastStationOptimizedBarrelLowPtTight");
		muNeg_oneStationAngTight=muonNeg->muonID("TMOneStationAngTight");
		muNeg_dxy=iTrack->dxy(RefVtx);
		muNeg_dz=iTrack->dz(RefVtx);
		if(muonNeg->isGlobalMuon())
		{
			TrackRef gTrack =muonNeg->globalTrack();
			const reco::HitPattern& q2=gTrack->hitPattern();
			muNeg_nValidMuHits=q2.numberOfValidMuonHits();
			muNeg_nchi2Gl=gTrack->chi2()/gTrack->ndof();
		//cout<<"======================q2.numberOfValidMuonHits(): "<<q2.numberOfValidMuonHits()<<endl;
		//cout<<"======================gTrack->chi2()/gTrack->ndof(): "<<gTrack->chi2()/gTrack->ndof()<<endl;
		}
		/*
		cout<<"======================iTrack->found(): "<<iTrack->found()<<endl;
		cout<<"======================p2.pixelLayersWithMeasurement(): "<<p2.pixelLayersWithMeasurement()<<endl;
		cout<<"======================iTrack->chi2()/iTrack->ndof(): "<<iTrack->chi2()/iTrack->ndof()<<endl;
		cout<<"======================iTrack->dxy(RefVtx): "<<iTrack->dxy(RefVtx)<<endl;
		cout<<"======================iTrack->ddz(RefVtx): "<<iTrack->dz(RefVtx)<<endl;
		cout<<"======================muNeg_arbitrated: "<<muNeg_arbitrated<<endl;
		cout<<"======================muNeg_oneStaTight: "<<muNeg_oneStaTight<<endl;
		cout<<endl;
		*/
	}

  //write out Calculated Polarization variables

  //! Fill Polarization Variables;
  // std::vector< float > thisCosTh, thisPhi;
  // thisCosTh.resize(6); thisPhi.resize(6);
  // this->calcPol(*muPosP, *muNegP, thisCosTh, thisPhi);

  if (_writeDataSet) {

    if (theMass > JpsiMassMin && theMass < JpsiMassMax && 
	theCtau > JpsiCtMin && theCtau < JpsiCtMax && 
	aCand->pt() > JpsiPtMin && aCand->pt() < JpsiPtMax && 
	fabs(theRapidity) > JpsiRapMin && fabs(theRapidity) < JpsiRapMax &&
	isMuonInAccept(muon1) && isMuonInAccept(muon2) &&
	mapTriggerNameToIntFired_[_triggerForDataset] == 1) {
      

      int ss=999;
      if (muon1->charge() + muon2->charge() == 0) ss=0;
      if (muon1->charge() + muon2->charge() == 2) ss=1;
      if (muon1->charge() + muon2->charge() == -2) ss=2;

      Jpsi_Sign->setIndex(ss,kTRUE);
      
      Jpsi_Pt->setVal(aCand->pt()); 
      Jpsi_Rap->setVal(theRapidity); 
      Jpsi_Mass->setVal(theMass);
      Jpsi_ct->setVal(theCtau);
      Jpsi_ctErr->setVal(theCtauErr);
      // cout << "Type = " << theCat << " pt = " << aCand->pt() << " eta = " << theRapidity << endl;
      // cout << " PPDL = " << theCtau << " Mother = " << aCand->userInt("momPDGId") << " PPDL true = " << 10.*aCand->userFloat("ppdlTrue") << endl;
      Jpsi_Type->setIndex(theCat,kTRUE);
      Jpsi_MatchType->setIndex((int)isMatched,kTRUE);
      Jpsi_ctTrue->setVal(10.*aCand->userFloat("ppdlTrue"));
    
      Jpsi_PtType->setIndex(getJpsiVarType(aCand->pt(),_ptbinranges),kTRUE);
      Jpsi_RapType->setIndex(getJpsiVarType(fabs(theRapidity),_etabinranges),kTRUE);
      // Fill RooDataSet
      RooArgSet varlist_tmp(*Jpsi_Mass,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap,*Jpsi_Type,*Jpsi_MatchType);   // temporarily remove tag-and-probe weights
      varlist_tmp.add(*Jpsi_ctTrue);   varlist_tmp.add(*Jpsi_PtType);
      varlist_tmp.add(*Jpsi_RapType);  varlist_tmp.add(*Jpsi_ctErr);
      varlist_tmp.add(*Jpsi_Sign);
      data->add(varlist_tmp);
    }
  }
}
        
void JPsiAnalyzerPAT::makeCuts(int sign) {

  if (collAll.isValid()) {

    for(vector<pat::CompositeCandidate>::const_iterator it=collAll->begin();
	it!=collAll->end();++it) {
      
      const pat::CompositeCandidate* cand = &(*it);	
      // cout << "Now checking candidate of type " << theJpsiCat << " with pt = " << cand->pt() << endl;
      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));
 
      bool thisSign = (sign == 0 && muon1->charge() + muon2->charge() == 0) || 
	(sign == 1 && muon1->charge() + muon2->charge() == 2) || 
	(sign == 2 && muon1->charge() + muon2->charge() == -2);
			double vprob=0.001;

      if (thisSign) {	  
	  
        // global + global?
	if (muon1->isGlobalMuon() && muon2->isGlobalMuon() &&
	    muon1->isTrackerMuon() && muon2->isTrackerMuon()   ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selGlobalMuon(muon2) &&
			      cand->userFloat("vProb") > vprob )) {
	    _thePassedCats[sign].push_back(0);  _thePassedCands[sign].push_back(cand);
            continue;
	  }
	}
	
        // global + tracker? (x2)    
	if (muon1->isGlobalMuon() && muon2->isTrackerMuon() &&
	    muon1->isTrackerMuon()    ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selTrackerMuon(muon2) &&
			      cand->userFloat("vProb") > vprob )) {
	    _thePassedCats[sign].push_back(1);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        if (muon2->isGlobalMuon() && muon1->isTrackerMuon() &&
            muon2->isTrackerMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon2) &&
			      selTrackerMuon(muon1) &&
			      cand->userFloat("vProb") > vprob )) {
	    _thePassedCats[sign].push_back(1);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        // tracker + tracker?  
        if (muon1->isTrackerMuon() && muon2->isTrackerMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon1) &&
			      selTrackerMuon(muon2) &&
			      cand->userFloat("vProb") > vprob )) {
	    _thePassedCats[sign].push_back(2);  _thePassedCands[sign].push_back(cand);
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

      bool thisSign = (sign == 0 && muon1->charge() + muon2->charge() == 0) || 
	(sign == 1 && muon1->charge() + muon2->charge() == 2) || 
	(sign == 2 && muon1->charge() + muon2->charge() == -2);

      if (thisSign && !(muon1->isTrackerMuon() && muon2->isTrackerMuon()) ) {

	// global + calo? (x2)
	if (muon1->isGlobalMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats[sign].push_back(3);  _thePassedCands[sign].push_back(cand);
            continue;
	  }
	}

	if (muon2->isGlobalMuon() && muon1->isCaloMuon() ) {
	  if (!_applycuts || (selGlobalMuon(muon2) &&
			      selCaloMuon(muon1) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats[sign].push_back(3);  _thePassedCands[sign].push_back(cand);
            continue;
	  }
	}
	
        // tracker + calo? (x2)    
	if (muon1->isTrackerMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats[sign].push_back(4);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        if (muon2->isTrackerMuon() && muon1->isCaloMuon() ) {
	  if (!_applycuts || (selTrackerMuon(muon2) &&
			      selCaloMuon(muon1) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats[sign].push_back(4);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}

        // calo + calo? 
        if (muon1->isCaloMuon() && muon2->isCaloMuon() ) {
	  if (!_applycuts || (selCaloMuon(muon1) &&
			      selCaloMuon(muon2) &&
			      cand->userFloat("vProb") > 0.001 )) {
	    _thePassedCats[sign].push_back(5);  _thePassedCands[sign].push_back(cand);
	    continue;
	  }
	}
      }
    }
  }

  return;
}

pair< unsigned int, const pat::CompositeCandidate* > 
JPsiAnalyzerPAT::theBestQQ(int sign) {

  unsigned int theBestCat = 99;
  const pat::CompositeCandidate* theBestCand = new pat::CompositeCandidate();

  for( unsigned int i = 0; i < _thePassedCands[sign].size(); i++) { 
    if (_thePassedCats[sign].at(i) < theBestCat) {
      theBestCat = _thePassedCats[sign].at(i);
      theBestCand = _thePassedCands[sign].at(i);
    }
  }

  pair< unsigned int, const pat::CompositeCandidate* > result = make_pair(theBestCat, theBestCand );
  return result;

}

bool
JPsiAnalyzerPAT::isMuonInAccept(const pat::Muon* aMuon) {
   // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
   return (fabs(aMuon->eta()) < 2.4 &&
           ((fabs(aMuon->eta()) < 1.3 && aMuon->pt() > 3.3) ||
           (fabs(aMuon->eta()) > 1.3 && fabs(aMuon->eta()) < 2.2 && aMuon->p() > 2.9) ||
           (fabs(aMuon->eta()) > 2.2 && aMuon->pt() > 0.8)));

   // *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
   // by just returning TRUE
   //  return true;
}

	 bool
	 JPsiAnalyzerPAT::selGlobalMuon(const pat::Muon* aMuon) {

	 TrackRef iTrack = aMuon->innerTrack();
	 const reco::HitPattern& p = iTrack->hitPattern();

	 TrackRef gTrack = aMuon->globalTrack();
	 const reco::HitPattern& q = gTrack->hitPattern();

	 return (// isMuonInAccept(aMuon) &&
	 iTrack->found() > 11 &&
	 gTrack->chi2()/gTrack->ndof() < 20.0 &&
	 q.numberOfValidMuonHits() > 0 &&
	 iTrack->chi2()/iTrack->ndof() < 4.0 &&
	 aMuon->muonID("TrackerMuonArbitrated") &&
	 aMuon->muonID("TMOneStationTight") &&
	 p.pixelLayersWithMeasurement() > 1 &&
	 fabs(iTrack->dxy(RefVtx)) < 3.0 &&
	 fabs(iTrack->dz(RefVtx)) < 15.0 );
	 }

	 bool
	 JPsiAnalyzerPAT::selTrackerMuon(const pat::Muon* aMuon) {

	 TrackRef iTrack = aMuon->innerTrack();
	 const reco::HitPattern& p = iTrack->hitPattern();

	 return (// isMuonInAccept(aMuon) &&
	 iTrack->found() > 11 &&
	 iTrack->chi2()/iTrack->ndof() < 4.0 &&
	 aMuon->muonID("TrackerMuonArbitrated") &&
	 aMuon->muonID("TMOneStationTight") &&
	 p.pixelLayersWithMeasurement() > 1 &&
	 fabs(iTrack->dxy(RefVtx)) < 3.0 &&
	 fabs(iTrack->dz(RefVtx)) < 15.0 );
	 }

bool 
JPsiAnalyzerPAT::selCaloMuon(const pat::Muon* aMuon) {

	TrackRef iTrack = aMuon->innerTrack();
	const reco::HitPattern& p = iTrack->hitPattern();

	return (// isMuonInAccept(aMuon) &&
			aMuon->caloCompatibility() > 0.89 &&
			iTrack->found() > 11 &&
			iTrack->chi2()/iTrack->ndof() < 4.0 &&
			p.pixelLayersWithMeasurement() > 1 &&
			fabs(iTrack->dxy(RefVtx)) < 3.0 &&
			fabs(iTrack->dz(RefVtx)) < 15.0 );
}

int 
JPsiAnalyzerPAT::getJpsiVarType(const double jpsivar, vector<double> vectbin) {

  for(unsigned int i=0;i<vectbin.size()-1;i++) {
    if(jpsivar > vectbin[i] && jpsivar < vectbin[i+1]) return i+1;
  }

  return -999;
}

// reset the global DataSet variables
void
JPsiAnalyzerPAT::resetDSVariables(){

    //reset J/psi RECO variables
     JpsiMass=-9999.;
     JpsiPt=-9999.;
     JpsiRap=-9999.;
    JpsiCharge=-9999;
     JpsiPx=-9999.;
     JpsiPy=-9999.;
     JpsiPz=-9999.;
    Jpsict=-9999.;
    JpsictErr=-9999.;
    Jpsict_Gen=-9999.;
    JpsiVprob=-9999.;

    JpsiType=-1;

    //reset MUON RECO variables
     muPosPx=-9999.;
    muPosPy=-9999.;
    muPosPz=-9999.;
    muNegPx=-9999.;
    muNegPy=-9999.;
    muNegPz=-9999.;
		//-----------------------------------------------------
		//-----------additional Reco Muon Variables------------
		//-----------------------------------------------------
		//1).Positive Muon
		muPos_nchi2In=-9999.;
		muPos_dxy=-9999.;
		muPos_dz=-9999.;
		muPos_nchi2Gl=-9999.;
		muPos_arbitrated=0;
		muPos_found=-9999;
		muPos_pixeLayers=-9999;
		muPos_nValidMuHits=-9999;
		muPos_oneStationTight=0;
		muPos_lastStationAngTight=0;
		muPos_lastStationTight=0;
		muPos_lastStationOptimizedLowPtTight=0;
		muPos_lastStationOptimizedBarrelLowPtTight=0;
		muPos_oneStationAngTight=0;
		//2).Negtive Muon
		muNeg_nchi2In=-9999.;
		muNeg_dxy=-9999.;
		muNeg_dz=-9999.;
		muNeg_nchi2Gl=-9999.;
		muNeg_arbitrated=0;
		muNeg_found=-9999;
		muNeg_pixeLayers=-9999;
		muNeg_nValidMuHits=-9999;
		muNeg_oneStationTight=0;
		muNeg_lastStationAngTight=0;
		muNeg_lastStationTight=0;
		muNeg_lastStationOptimizedLowPtTight=0;
		muNeg_lastStationOptimizedBarrelLowPtTight=0;
		muNeg_oneStationAngTight=0;


    if(_isMC){
        MCType=-1;

        //reset J/psi GEN variables
        /* JpsiMass_Gen=-9999.;
        JpsiPt_Gen=-9999.;
        JpsiRap_Gen=-9999.;
        JpsiPx_Gen=-9999.;
        JpsiPy_Gen=-9999.;
        JpsiPz_Gen=-9999.; */

        //reset MUON GEN variables
        /* muPosPx_Gen=-9999.;
        muPosPy_Gen=-9999.;
        muPosPz_Gen=-9999.;
        muNegPx_Gen=-9999.;
        muNegPy_Gen=-9999.;
        muNegPz_Gen=-9999.; */
    }

    //reset EVENT information
    eventNb= 0 ;
    runNb= 0 ;
    nPriVtx= 0 ;
    lumiBlock= 0 ;

    //reset Trigger Variables
    for(std::map< std::string, int >::iterator clearIt= mapTriggerNameToIntFired_.begin(); clearIt != mapTriggerNameToIntFired_.end(); clearIt++){
        clearIt->second=0;
    }
}

//! fill Generator Information
void
JPsiAnalyzerPAT::analyzeGenerator(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
    using namespace trigger;

    std::vector < const reco::Candidate* > genMuons;
    //bool genjpsi= false;
    reco::Candidate::size_type nrD;

    //int count= 0;
    for( size_t i = 0; i < genParticles->size(); ++ i )
    {
        // std::cout << "analyzeGenerator: " << i << std::endl;
        const reco::Candidate & cand = (*genParticles)[ i ];
        int Mc_particleID = cand.pdgId();
        if (abs(Mc_particleID) == _oniaPDG && cand.status()==2 )//&& cand.pt() >= 1)
        {
//          std::cout << "------::analyzeGenerator:: gen JPsi's: pt=" << cand.pt() << "; eta=" << cand.eta() << "; phi=" << cand.phi() << std::endl;

            //Fill global TTree variables
            // JpsiMass_Gen=cand.mass();
            // JpsiPt_Gen=cand.pt();
            // JpsiRap_Gen=cand.rapidity();
            // JpsiPx_Gen=cand.px();
            // JpsiPy_Gen=cand.py();
            // JpsiPz_Gen=cand.pz();
          Double_t enGen = sqrt(cand.px()*cand.px() + cand.py()*cand.py() + cand.pz()*cand.pz() + cand.mass()*cand.mass());
	  JpsiP_Gen->SetPxPyPzE(cand.px(),cand.py(),cand.pz(),enGen);

            //Jpsict_Gen=10.*cand.userFloat("ppdlTrue"));

            nrD= cand.numberOfDaughters();
            int count_muon=0;
            for(reco::Candidate::size_type t=0; t < nrD; t++){
                const reco::Candidate* muon= cand.daughter(t);
                int pID = muon->pdgId();
//              std::cout << "------::analyzeGenerator:: gen JPsi's daughter pdgId: " << pID << std::endl;

                if (abs(pID) == 13 && cand.daughter(t)->status()==1)
                {
                    genMuons.push_back(muon);
//                  std::cout << "------::analyzeGenerator:: gen JPsi's daughter #: " << count_muon << std::endl;
//                  std::cout << " muon" << count_muon << " pt=     " << muon->pt() << std::endl;
//                  std::cout << " muon" << count_muon << " eta=     " << muon->eta() << std::endl;
//                  std::cout << " moun" << count_muon << " phi=     " << muon->phi() << std::endl;
                    count_muon++;
                }
            }


            if ( genMuons.empty() ) break;

            const reco::Candidate* muon1= genMuons.front();
            const reco::Candidate* muon2= genMuons.back();

            // look for opposite charge gen muon pair
            if (muon1->charge()*muon2->charge() <= 0){
                const reco::Candidate *muonPos = 0, *muonNeg = 0;

                if(muon1->charge() > 0){ muonPos = muon1; muonNeg = muon2;}
                else if(muon1->charge() < 0){ muonPos = muon2; muonNeg = muon1;}

                float f_muPosPx, f_muPosPy, f_muPosPz;
                float f_muNegPx, f_muNegPy, f_muNegPz;

                f_muPosPx = muonPos->px();
                f_muPosPy = muonPos->py();
                f_muPosPz = muonPos->pz();

                f_muNegPx = muonNeg->px();
                f_muNegPy = muonNeg->py();
                f_muNegPz = muonNeg->pz();

                // fill global TTree variables - gen muon
                // muPosPx_Gen=muonPos->px();
                // muPosPy_Gen=muonPos->py();
                // muPosPz_Gen=muonPos->pz();

                // muNegPx_Gen=muonNeg->px();
                // muNegPy_Gen=muonNeg->py();
                // muNegPz_Gen=muonNeg->pz();

                // fill Polarization variables - gen muons
                Double_t muMass = 0.105658;

                Double_t enMuPos = sqrt(f_muPosPx*f_muPosPx + f_muPosPy*f_muPosPy + f_muPosPz*f_muPosPz + muMass*muMass);
                // TLorentzVector *muPos = new TLorentzVector();
                muPosP_Gen->SetPxPyPzE(f_muPosPx, f_muPosPy, f_muPosPz, enMuPos);

                Double_t enMuNeg = sqrt(f_muNegPx*f_muNegPx + f_muNegPy*f_muNegPy + f_muNegPz*f_muNegPz + muMass*muMass);
                // TLorentzVector *muNeg = new TLorentzVector();
                muNegP_Gen->SetPxPyPzE(f_muNegPx, f_muNegPy, f_muNegPz, enMuNeg);

                //! Fill Polarization Variables;
                // std::vector< float > thisCosTh, thisPhi;
                // thisCosTh.resize(6); thisPhi.resize(6);
                // this->calcPol(*muPosP_Gen, *muNegP_Gen, thisCosTh, thisPhi);
     
            }
        } // end loop over genParticles
    }
}

void
JPsiAnalyzerPAT::hltReport(const edm::Event &iEvent ,const edm::EventSetup& iSetup)
{

    std::map<std::string, bool> mapTriggernameToTriggerFired;
    std::map<std::string, unsigned int> mapTriggernameToHLTbit;
    std::map<std::string, unsigned int> mapTriggerNameToPrescaleFac;

    for(std::vector<std::string>::const_iterator it= HLTbitNames_.begin(); it !=HLTbitNames_.end(); ++it){
        mapTriggernameToTriggerFired[*it]=false;
        mapTriggernameToHLTbit[*it]=1000;
        mapTriggerNameToPrescaleFac[*it]=0;
    }

    // HLTConfigProvider
    if ( hltConfigInit_ ) {
        
        //! Use HLTConfigProvider
      const unsigned int n= hltConfig_.size();
      for (std::map<std::string, unsigned int>::iterator it = mapTriggernameToHLTbit.begin(); it != mapTriggernameToHLTbit.end(); it++) {
	unsigned int triggerIndex= hltConfig_.triggerIndex( it->first );
	if (triggerIndex >= n) {
	  //std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " not available in config!" << std::endl;
            }
	else {
	  it->second= triggerIndex;
                //std::cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerName " << it->first << " available in config!" << std::endl;
	}
      }
    }

    // Get Trigger Results
    try {
    iEvent.getByLabel( tagTriggerResults_, handleTriggerResults_ );
    //cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResult is present in current event" << endl;
    }
    catch(...) {
    //cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults NOT present in current event" << endl;
    }
    if ( handleTriggerResults_.isValid() ){
    //cout << "[JPsiAnalyzerPAT::hltReport] --- J/psi TriggerResults IS valid in current event" << endl;

    // loop over Trigger Results to check if paths was fired
    for(std::vector< std::string >::iterator itHLTNames= HLTbitNames_.begin(); itHLTNames != HLTbitNames_.end(); itHLTNames++){
      const std::string triggerPathName =  *itHLTNames;
      //std::cout << "[FloJPsiAnalyzer::hltReport] --- TriggerName --- TriggerName LOOP" << std::endl;

      if ( mapTriggernameToHLTbit[triggerPathName] < 1000 && handleTriggerResults_->accept( mapTriggernameToHLTbit[triggerPathName] ) ){
          //std::cout << "[FloJPsiAnalyzer::hltReport] --- TriggerName " << triggerPathName << " fired!" << std::endl;
          mapTriggerNameToIntFired_[triggerPathName] = 2;
      }
    }
    }
    else cout << "[JPsiAnalyzerPAT::hltReport] --- TriggerResults NOT valid in current event" << endl;
}

void
JPsiAnalyzerPAT::matchMuonToHlt(const pat::Muon* muon1, const pat::Muon* muon2)
{
    //! Loop over Trigger Paths and match muons to last Filter/collection
    for ( std::map<std::string, int>::iterator it = mapTriggerNameToIntFired_.begin(); it != mapTriggerNameToIntFired_.end(); it ++ ) {

        std::string triggerName = it->first;

        //! just use Triggers which are in TriggerResults; value == 2
        if ( it->second != 2 ) continue;

        std::string hltLastFilterName = mapTriggerToLastFilter_[triggerName];

        const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter( hltLastFilterName );
        const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter( hltLastFilterName );
        bool pass1 = mu1HLTMatches.size() > 0;
        bool pass2 = mu2HLTMatches.size() > 0;

        // treat "MuX_TrackX" Trigger separately: Match by Tracker collection: hltMuTrackJpsiCtfTrackCands
        if (    triggerName == "HLT_Mu0_Track0_Jpsi" ||
                triggerName == "HLT_Mu3_Track0_Jpsi" ||
                triggerName == "HLT_Mu5_Track0_Jpsi" ||

                triggerName == "HLT_Mu3_Track3_Jpsi"    ||
                triggerName == "HLT_Mu3_Track3_Jpsi_v2" ||
                triggerName == "HLT_Mu3_Track5_Jpsi_v1" ||
                triggerName == "HLT_Mu3_Track5_Jpsi_v2"     )
        {
                bool matchedMu3[2] = {false, false}, matchedTrack[2] = {false, false};
                for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
                    if (mu1HLTMatches[k].collection() == "hltL3MuonCandidates::HLT") matchedMu3[0] = true;
                    if (mu1HLTMatches[k].collection() == "hltMuTrackJpsiCtfTrackCands::HLT" ) matchedTrack[0] = true;
                }
                for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
                    if (mu2HLTMatches[k].collection() == "hltL3MuonCandidates::HLT") matchedMu3[1] = true;
                    if (mu2HLTMatches[k].collection() == "hltMuTrackJpsiCtfTrackCands::HLT") matchedTrack[1] = true;
                }
                if( (matchedMu3[0] && matchedTrack[1]) || (matchedMu3[1] && matchedTrack[0]) ) {
                    mapTriggerNameToIntFired_[triggerName] = 1;
                    //std::cout << "[JPsiAnalyzerPAT::matchMuonToHlt] ---- ---- \"MuX + Track\" Trigger: " << triggerName << " FIRED and MATCHED" << std::endl;
                }
        }

        // treat "MuX_TkMuX" Trigger separately: Match by Tracker collection:hltMuTkMuJpsiTrackerMuonCands
        if (    triggerName == "HLT_Mu0_TkMu0_Jpsi" ||
                triggerName == "HLT_Mu3_TkMu0_Jpsi" ||
                triggerName == "HLT_Mu5_TkMu0_Jpsi" ||
             //
                triggerName == "HLT_Mu0_TkMu0_OST_Jpsi" ||
                triggerName == "HLT_Mu3_TkMu0_OST_Jpsi" ||
                triggerName == "HLT_Mu5_TkMu0_OST_Jpsi" ||
             //
                triggerName == "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1" ||
                triggerName == "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1" ||
                triggerName == "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1" ||
//
                triggerName == "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2" ||
                triggerName == "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2" ||
                triggerName == "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2" ||
//
                triggerName == "HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3" ||
                triggerName == "HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3" ||
                triggerName == "HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3" )
        {
                bool matchedMu3[2] = {false, false}, matchedTrack[2] = {false, false};
                for (unsigned k = 0; k < mu1HLTMatches.size(); ++k) {
                    if (mu1HLTMatches[k].collection() == "hltL3MuonCandidates::HLT") matchedMu3[0] = true;
                    if (mu1HLTMatches[k].collection() == "hltMuTkMuJpsiTrackerMuonCands::HLT") matchedTrack[0] = true;
                }
                for (unsigned k = 0; k < mu2HLTMatches.size(); ++k) {
                    if (mu2HLTMatches[k].collection() == "hltL3MuonCandidates::HLT") matchedMu3[1] = true;
                    if (mu2HLTMatches[k].collection() == "hltMuTkMuJpsiTrackerMuonCands::HLT") matchedTrack[1] = true;
                }
                if( (matchedMu3[0] && matchedTrack[1]) || (matchedMu3[1] && matchedTrack[0]) ) {
                    mapTriggerNameToIntFired_[triggerName] = 1;
                    //std::cout << "[JPsiAnalyzerPAT::matchMuonToHlt] ---- ---- \"MuX + Track\" Trigger: " << triggerName << " FIRED and MATCHED" << std::endl;
                }
        }

        // All the other Paths match by last filter:
        // double muon trigger:
        if ( triggerName == "HLT_DoubleMu0"                     && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_DoubleMu0_Quarkonium_v1"       && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_DoubleMu0_Quarkonium_LS_v1"    && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_L1DoubleMuOpen"                && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_L1DoubleMuOpen_Tight"          && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_DoubleMu3"                     && pass1 == true && pass2 == true ) mapTriggerNameToIntFired_[triggerName] = 1;

        // single muon trigger:
        if ( triggerName == "HLT_Mu3" && (pass1  == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_Mu5" && (pass1  == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_Mu7" && (pass1  == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_Mu9" && (pass1  == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
        if ( triggerName == "HLT_Mu11"&& (pass1  == true || pass2 == true) ) mapTriggerNameToIntFired_[triggerName] = 1;
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalyzerPAT);
