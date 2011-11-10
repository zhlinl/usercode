// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc B2MuMuX/Analyzer/src/Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Linlin, Zhang 
//         Created:  Tue Nov. 7 2011
// $Id: Analyzer.h,v 1.4 2011/11/07 20:18:52 zhlinl Exp $
//
//

#ifndef _Analyzer_h
#define _Analyzer_h


// system include files
#include <memory>
#include <fstream>
#include <ostream>
#include <iostream>
#include <math.h>

// ROOT/Roofit include files
#include <TStyle.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TClonesArray.h>
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include "FWCore/Common/interface/TriggerNames.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/BeamSpot/interface/BeamSpot.h>
#include <MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"


#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "TMath.h"
#include "Math/VectorUtil.h"

//#include "OniaMu.cc"

//
// class declaration
//
using namespace std;
using namespace edm;
using namespace reco;
using namespace RooFit;

const double muonMass = 0.105658367;
const double muonMass2  = muonMass*muonMass;

const double jpsiMass = 3.096916;
const double jpsiMass2  = jpsiMass*jpsiMass;

class Analyzer : public edm::EDAnalyzer {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      void MCAna(const edm::Event& );
      void hltReport(const edm::Event & ,const edm::EventSetup& );
      void Onia2MuMu(const edm::Event & ,const edm::EventSetup&);

      void OniaDitrack(const pat::CompositeCandidate* aCand,const edm::Event&, const edm::EventSetup&);
      void OniaTrack(const pat::CompositeCandidate* aCand,const edm::Event&, const edm::EventSetup&); 
      void OniaMu(const pat::CompositeCandidate* aCand,const edm::Event&, const edm::EventSetup&);

      bool selGlobalMuon(const pat::Muon* aMuon);
      bool selTrackerMuon(const pat::Muon* aMuon);
      bool selTrack(TrackRef atrkRef);
      bool isMuonInAccept(const pat::Muon* aMuon);

      void resetVariables();

      // ----------member data ---------------------------

      // ROOT tree 
      TTree* tree_;//data; //*recoData;
      TFile* fOut_;
      string   _treefilename; 

      //-----

      Handle< vector<pat::Muon> > allmuons;
      //----------------------InputTag------------------- 
      bool AnalyzeMC_;
      InputTag  muons_;
      StringCutObjectSelector<pat::Muon> higherPuritySelection_;
      StringCutObjectSelector<pat::Muon> lowerPuritySelection_; 
      StringCutObjectSelector<reco::Candidate, true> dimuonSelection_;
      bool           _applycuts;
      bool           _applyExpHitcuts;
      bool           _applyDiMuoncuts;

      bool TriMuReco_; 
      bool OniaTrackReco_;
      bool OniaDiTrackReco_;

      InputTag  _trackLabelPi;
      double PiPt_c;
      double OniaPiDR_c;

      //-hlt-
      InputTag tagTriggerResults_;
      vector<string> HLTBitNames_DoubleMu;
      vector<string> HLTLastFilterNames_DoubleMu;

      vector<unsigned int> *itriggerflag;
      vector<string> *itriggerNames;     
      vector<unsigned int> *Dimutriggerflag;
      vector<string> *DimutriggerNames;
      int  DimuTriggerResult[150];
      vector<bool> *JPsiMuonTrigMatch;
      //---------------------------------------------

      unsigned int eventNb, runNb, lumiBlock, nPriVtx, nEvents;
      Vertex thePrimaryV;
      Vertex theBeamSpotV;
      math::XYZPoint RefVtx;

      vector<pat::CompositeCandidate>   _thePassedCands;

      //--varialbles---
      int nGenKstar;
      TClonesArray* Gen_Jpsi_P4, *Gen_B_P4, *Gen_B_V3, *Gen_Kstar_P4, *Gen_Kaon_P4, *Gen_Pion_P4,
				*Gen_muonPos_P4, *Gen_muonNeg_P4;
      vector<int> *Gen_B_pdgid,*Gen_Kstar_pdgid,*Gen_Kaon_pdgid,*Gen_Pion_pdgid;
      
      int nonia;
      TClonesArray *OniaP4, *MuPP4, *MuMP4;
      vector<int> *OniaCats;
      vector<float> *OniaVtxCL;

      int nXcand;  
      TClonesArray *Pi1P4, *Pi2P4, *KstarP4, *xP4;      
      vector<int>   *xOindex, *Pi1NHits, *Pi1PixelHits, *Pi2NHits, *Pi2PixelHits;
      vector<float> *Pi1D0, *Pi1Dz,  *Pi1NormChi2;
      vector<float> *Pi2D0, *Pi2Dz,  *Pi2NormChi2;
      vector<float> *xVtxCL, *xVtxC2;   
      vector<double>  *xcosAlpha, *xlxySig, *xctauPV, *xctauErrPV, *xlxyPV;  

      int nOniaMu;
      TClonesArray *Mu3P4, *TrMP4;
      vector<int>  *TrMOindex,*TrMCats, *Mu3NHits, *Mu3PixelHits;
      vector<float> *Mu3D0, *Mu3Dz,  *Mu3NormChi2,*TrMVtxCL, *TrMVtxC2;
      vector<double>  *TrMcosAlpha, *TrMlxySig, *TrMctauPV, *TrMctauErrPV, *TrMlxyPV;  

      int nOniaTrack;
      TClonesArray *PiP4, *BP4;
      vector<int> *BOindex, *PiNHits, *PiPixelHits;
      vector<float> *PiD0, *PiDz, *PiNormChi2, *BVtxCL, *BVtxC2; 
      vector<double>  *BcosAlpha, *BlxySig, *BctauPV, *BctauErrPV, *BlxyPV;  

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig):
   _treefilename(iConfig.getParameter<string>("treeFileName")),
   AnalyzeMC_(iConfig.getUntrackedParameter<bool>("AnalyzeMC",false)),
   muons_(iConfig.getParameter<InputTag>("src")),
   higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),//Muon quaulity cut
   lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
   dimuonSelection_(iConfig.existsAs<std::string>("dimuonSelection") ? iConfig.getParameter<std::string>("dimuonSelection") : ""),
   _applycuts(iConfig.getParameter<bool>("applyCuts")),
   _applyExpHitcuts(iConfig.getUntrackedParameter<bool>("applyExpHitCuts",false)), //Muon quaulity cut
   _applyDiMuoncuts(iConfig.getUntrackedParameter<bool>("applyDiMuonCuts",false)), //Muon quaulity cut 

   TriMuReco_(iConfig.getUntrackedParameter<bool>("TriMuRECO",false)), //This is for three Muon selection don't need for Analyzer selection
   OniaTrackReco_(iConfig.getUntrackedParameter<bool>("OniaTrackRECO",true)), //Onia plus one track
   OniaDiTrackReco_(iConfig.getUntrackedParameter<bool>("OniaDiTrackRECO",true)),//Onia plus two tracks

   _trackLabelPi(iConfig.getParameter<InputTag>("TrackLabel")),
   PiPt_c(iConfig.getUntrackedParameter<double>("MinTrPt", 0.4)), //Pion pt cut
   OniaPiDR_c(iConfig.getUntrackedParameter<double>("OniaPiPiMaxDR", 1)), //deltaR cut 

   tagTriggerResults_(iConfig.getParameter<InputTag>("triggerResultsLabel")),
   HLTBitNames_DoubleMu(iConfig.getParameter< vector<string> >("HLTBitNames_DoubleMu")), //HLT trigger Menu
   HLTLastFilterNames_DoubleMu(iConfig.getParameter< vector<string> >("HLTLastFilterNames_DoubleMu")), //HLT trigger filter Menu
   itriggerflag(0), itriggerNames(0), Dimutriggerflag(0), DimutriggerNames(0), JPsiMuonTrigMatch(0),

   Gen_B_pdgid(0), Gen_Kstar_pdgid(0), Gen_Kaon_pdgid(0), Gen_Pion_pdgid(0),

   OniaCats(0), OniaVtxCL(0),

   xOindex(0), Pi1NHits(0), Pi1PixelHits(0), Pi2NHits(0), Pi2PixelHits(0),
   Pi1D0(0), Pi1Dz(0), Pi1NormChi2(0), Pi2D0(0), Pi2Dz(0), Pi2NormChi2(0),
   xVtxCL(0), xVtxC2(0), xcosAlpha(0), xlxySig(0), xctauPV(0), xctauErrPV(0), xlxyPV(0),

   TrMOindex(0), TrMCats(0), Mu3NHits(0), Mu3PixelHits(0), 
   Mu3D0(0), Mu3Dz(0), Mu3NormChi2(0), TrMVtxCL(0), TrMVtxC2(0), 
   TrMcosAlpha(0), TrMlxySig(0), TrMctauPV(0), TrMctauErrPV(0), TrMlxyPV(0),

    BOindex(0), PiNHits(0), PiPixelHits(0),
   PiD0(0), PiDz(0),  PiNormChi2(0), BVtxCL(0), BVtxC2(0),
   BcosAlpha(0), BlxySig(0), BctauPV(0), BctauErrPV(0), BlxyPV(0)
{
   nEvents=0;
   //now do what ever initialization is needed


}


Analyzer::~Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer::beginJob()
{
   // edm::Service<TFileService> fs;
   // tree_= fs->make<TTree>("data", "J/psi plus ntracks(Pi) Tree");
   fOut_ = new TFile(_treefilename.c_str(), "RECREATE");
   fOut_->cd();

   // TTree
   //load Branches
   tree_ = new TTree ("data", " J/psi plus ntracks(Pi) Tree");


   // Event variables
   tree_->Branch("eventNb",             &eventNb,             "eventNb/I");
   tree_->Branch("runNb",               &runNb,               "runNb/I");
   tree_->Branch("lumiBlock",           &lumiBlock,           "lumiBlock/I");
   tree_->Branch("nPriVtx",             &nPriVtx,             "nPriVtx/I");

   //add HLT Variables to TTree
   tree_->Branch("itriggerflag",        &itriggerflag);
   tree_->Branch("itriggerNames",       &itriggerNames);
   tree_->Branch("Dimutriggerflag",        &Dimutriggerflag);
   tree_->Branch("DimutriggerNames",       &DimutriggerNames);
   tree_->Branch("JPsiMuonTrigMatch",   &JPsiMuonTrigMatch);  

   if(AnalyzeMC_){
      Gen_Jpsi_P4=new TClonesArray("TLorentzVector", 10000);
      Gen_B_P4=new TClonesArray("TLorentzVector", 10000);
			Gen_B_V3=new TClonesArray("TVector3", 10000);
      Gen_Kstar_P4=new TClonesArray("TLorentzVector", 10000);
      Gen_Kaon_P4=new TClonesArray("TLorentzVector", 10000);
      Gen_Pion_P4=new TClonesArray("TLorentzVector", 10000);
      Gen_muonPos_P4=new TClonesArray("TLorentzVector", 10000);
      Gen_muonNeg_P4=new TClonesArray("TLorentzVector", 10000);

      tree_->Branch("nGenKstar",    &nGenKstar,    "nGenKstar/I");

      tree_->Branch("Gen_Jpsi_P4",  "TClonesArray",  &Gen_Jpsi_P4,   32000, 0);        
      tree_->Branch("Gen_B_P4",     "TClonesArray",  &Gen_B_P4,      32000, 0);        
			tree_->Branch("Gen_B_V3",     "TClonesArray",  &Gen_B_V3,      32000, 0);
      tree_->Branch("Gen_Kstar_P4",     "TClonesArray",  &Gen_Kstar_P4,      32000, 0);        
      tree_->Branch("Gen_Kaon_P4",     "TClonesArray",  &Gen_Kaon_P4,      32000, 0);        
      tree_->Branch("Gen_Pion_P4",     "TClonesArray",  &Gen_Pion_P4,      32000, 0);        
      tree_->Branch("Gen_muonPos_P4",     "TClonesArray",  &Gen_muonPos_P4,      32000, 0);        
      tree_->Branch("Gen_muonNeg_P4",     "TClonesArray",  &Gen_muonNeg_P4,      32000, 0);        

      tree_->Branch("Gen_B_pdgid",         &Gen_B_pdgid);  
      tree_->Branch("Gen_Kstar_pdgid",         &Gen_Kstar_pdgid);  
      tree_->Branch("Gen_Kaon_pdgid",         &Gen_Kaon_pdgid);  
      tree_->Branch("Gen_Pion_pdgid",         &Gen_Pion_pdgid);  

         }

   OniaP4=new TClonesArray("TLorentzVector", 10000);
   MuPP4=new TClonesArray("TLorentzVector", 10000);
   MuMP4=new TClonesArray("TLorentzVector", 10000);

   tree_->Branch("nonia",    &nonia,    "nonia/I");
   tree_->Branch("OniaP4",   "TClonesArray",   &OniaP4, 32000, 0);
   tree_->Branch("MuPP4",   "TClonesArray",    &MuPP4,  32000, 0);
   tree_->Branch("MuMP4",   "TClonesArray",    &MuMP4,  32000, 0);
   tree_->Branch("OniaCats",    &OniaCats);
   tree_->Branch("OniaVtxCL",   &OniaVtxCL);

   if (TriMuReco_){
      Mu3P4=new TClonesArray("TLorentzVector", 10000);
      TrMP4=new TClonesArray("TLorentzVector", 10000);

      tree_->Branch("nOniaMu",    &nOniaMu,    "nOniaMu/I");
      tree_->Branch("Mu3P4",  "TClonesArray",   &Mu3P4, 32000, 0);
      tree_->Branch("TrMP4",  "TClonesArray",   &TrMP4, 32000, 0);

      tree_->Branch("Mu3D0",        &Mu3D0);
      tree_->Branch("Mu3Dz",        &Mu3Dz);
      tree_->Branch("Mu3NHits",     &Mu3NHits);
      tree_->Branch("Mu3PixelHits", &Mu3PixelHits);
      tree_->Branch("Mu3NormChi2",  &Mu3NormChi2);

      tree_->Branch("TrMOindex",             &TrMOindex);
      tree_->Branch("TrMCats",               &TrMCats);
      tree_->Branch("TrMVtxCL",              &TrMVtxCL);
      tree_->Branch("TrMVtxC2",              &TrMVtxC2);
      

      tree_->Branch("TrMcosAlpha",           &TrMcosAlpha);
      tree_->Branch("TrMlxySig",             &TrMlxySig);
      tree_->Branch("TrMctauPV",             &TrMctauPV);
      tree_->Branch("TrMctauErrPV",          &TrMctauErrPV);
      tree_->Branch("TrMlxyPV",              &TrMlxyPV);
   }

   if(OniaTrackReco_){
      PiP4=new TClonesArray("TLorentzVector", 10000);
      BP4=new TClonesArray("TLorentzVector", 10000);

     
      tree_->Branch("nOniaTrack",    &nOniaTrack,    "nOniaTrack/I");
      tree_->Branch("PiP4",  "TClonesArray",   &PiP4, 32000, 0);
      tree_->Branch("BP4",   "TClonesArray",   &BP4,  32000, 0);

      tree_->Branch("PiD0",        &PiD0);
      tree_->Branch("PiDz",        &PiDz);
      tree_->Branch("PiNHits",     &PiNHits);
      tree_->Branch("PiPixelHits", &PiPixelHits);
      tree_->Branch("PiNormChi2",  &PiNormChi2);

      tree_->Branch("BOindex",             &BOindex);     
      tree_->Branch("BVtxCL",              &BVtxCL);
      tree_->Branch("BVtxC2",              &BVtxC2);
     
      tree_->Branch("BcosAlpha",           &BcosAlpha);
      tree_->Branch("BlxySig",             &BlxySig);
      tree_->Branch("BctauPV",             &BctauPV);
      tree_->Branch("BctauErrPV",          &BctauErrPV);
      tree_->Branch("BlxyPV",              &BlxyPV);
   }
   if(OniaDiTrackReco_){
      Pi1P4=new TClonesArray("TLorentzVector", 10000);
      Pi2P4=new TClonesArray("TLorentzVector", 10000);
      KstarP4=new TClonesArray("TLorentzVector", 10000);
      xP4=new TClonesArray("TLorentzVector", 10000);
    

      tree_->Branch("nXcand",    &nXcand,    "nXcand/I");
      tree_->Branch("Pi1P4",  "TClonesArray",  &Pi1P4, 32000, 0);
      tree_->Branch("Pi2P4",  "TClonesArray",  &Pi2P4, 32000, 0);
      tree_->Branch("KstarP4",    "TClonesArray",  &KstarP4,   32000, 0);
      tree_->Branch("xP4",    "TClonesArray",  &xP4,   32000, 0);

      tree_->Branch("Pi1D0",        &Pi1D0);
      tree_->Branch("Pi1Dz",        &Pi1Dz);
      tree_->Branch("Pi1NHits",     &Pi1NHits);
      tree_->Branch("Pi1PixelHits", &Pi1PixelHits);
      tree_->Branch("Pi1NormChi2",  &Pi1NormChi2);

      tree_->Branch("Pi2D0",        &Pi2D0);
      tree_->Branch("Pi2Dz",        &Pi2Dz);
      tree_->Branch("Pi2NHits",     &Pi2NHits);
      tree_->Branch("Pi2PixelHits", &Pi2PixelHits);
      tree_->Branch("Pi2NormChi2",  &Pi2NormChi2);


      tree_->Branch("xOindex",             &xOindex);    
      tree_->Branch("xVtxCL",              &xVtxCL);
      tree_->Branch("xVtxC2",              &xVtxC2);
     
      tree_->Branch("xcosAlpha",           &xcosAlpha);
      tree_->Branch("xlxySig",             &xlxySig);
      tree_->Branch("xctauPV",             &xctauPV);
      tree_->Branch("xctauErrPV",          &xctauErrPV);
      tree_->Branch("xlxyPV",              &xlxyPV);
   }


}

void 
Analyzer::resetVariables(){
   runNb = 0; eventNb=0; lumiBlock=0; nPriVtx=0;
  itriggerflag->clear(); itriggerNames->clear(); Dimutriggerflag->clear(); DimutriggerNames->clear(); JPsiMuonTrigMatch->clear();
   if(AnalyzeMC_){
      Gen_Jpsi_P4->Clear(); Gen_B_P4->Clear(); Gen_B_V3->Clear(); Gen_Kstar_P4->Clear(); Gen_Kaon_P4->Clear(); 
			Gen_Pion_P4->Clear(); Gen_muonPos_P4->Clear(); Gen_muonNeg_P4->Clear();
      Gen_B_pdgid->clear(); Gen_Kstar_pdgid->clear(); Gen_Kaon_pdgid->clear(); Gen_Pion_pdgid->clear();
   }
   OniaP4->Clear(); MuPP4->Clear(); MuMP4->Clear(); OniaCats->clear();  OniaVtxCL->clear(); 

   if (TriMuReco_){
      Mu3P4->Clear(); TrMP4->Clear(); 

      Mu3D0->clear(); Mu3Dz->clear(); Mu3NHits->clear(); Mu3PixelHits->clear(); Mu3NormChi2->clear();
      TrMOindex->clear();  TrMCats->clear(); 
      TrMVtxCL->clear(); TrMVtxC2->clear(); 
      TrMcosAlpha->clear(); TrMlxySig->clear(); TrMctauPV->clear(); TrMctauErrPV->clear(); TrMlxyPV->clear();  
   }

   if(OniaTrackReco_){
      PiP4->Clear(); BP4->Clear(); 

      PiD0->clear(); PiDz->clear(); PiNHits->clear(); PiPixelHits->clear(); PiNormChi2->clear();
      BOindex->clear(); 
      BVtxCL->clear(); BVtxC2->clear(); 
      BcosAlpha->clear(); BlxySig->clear(); BctauPV->clear(); BctauErrPV->clear(); BlxyPV->clear();  
   }

   if(OniaDiTrackReco_){
      Pi1P4->Clear(); Pi2P4->Clear(); KstarP4->Clear(); xP4->Clear(); 

      Pi1D0->clear(); Pi1Dz->clear(); Pi1NHits->clear(); Pi1PixelHits->clear(); Pi1NormChi2->clear();
      Pi2D0->clear(); Pi2Dz->clear(); Pi2NHits->clear(); Pi2PixelHits->clear(); Pi2NormChi2->clear();
      xOindex->clear();
      xVtxCL->clear(); xVtxC2->clear(); 
      xcosAlpha->clear(); xlxySig->clear(); xctauPV->clear(); xctauErrPV->clear(); xlxyPV->clear();
   }

}



// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer::endJob() 
{
   cout << "Total number of events = " << nEvents << endl;
   cout << "============================================================" << endl;
   //tree_->GetDirectory()->cd();
   fOut_->cd();
   tree_->Write();
   fOut_->Close();

}

// ------------ method called when starting to processes a run  ------------
void 
Analyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


#endif
