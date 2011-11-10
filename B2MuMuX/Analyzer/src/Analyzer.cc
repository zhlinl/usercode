// -*- C++ -*-
//
// Package:    B2KstarMuMu
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc B2KstarMuMu/Analyzer/src/Analyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Linlin, Zhang 
//         Created:  Tue Nov. 7 2011
// $Id: Analyzer.h,v 1.4 2011/11/07 20:18:52 zhlinl Exp $
//


#include "../interface/Analyzer.h"

// ------------ method called for each event  ------------
	void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	nEvents++;
	//cout<<">>Debug000"<<endl;

	// Event related infos
	eventNb= iEvent.id().event() ;
	runNb=iEvent.id().run() ;
	lumiBlock= iEvent.luminosityBlock() ;

	if (AnalyzeMC_) MCAna(iEvent);

	// check HLT TriggerReuslts
	this->hltReport(iEvent, iSetup);


	reco::BeamSpot beamSpot;
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
	if ( beamSpotHandle.isValid() ) {beamSpot = *beamSpotHandle; 
		theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}
	else cout << "No beam spot available from EventSetup" << endl;

	Handle<VertexCollection> recVtxs;
	iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
	nPriVtx = recVtxs->size();
	if ( recVtxs->begin() != recVtxs->end() ) {
		thePrimaryV = Vertex(*(recVtxs->begin()));
	}
	else {
		thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
	}

	RefVtx = thePrimaryV.position();
	_thePassedCands.clear();
	this->Onia2MuMu(iEvent,iSetup);

	//cout<<">>thePassedCands.size(): "<<_thePassedCands.size()<<endl;
	for(unsigned int i=0; i<_thePassedCands.size();i++){
		//Onia+1Muon
		if (TriMuReco_) OniaMu(&_thePassedCands[i], iEvent, iSetup);
		//Onia+1track
		if(OniaTrackReco_) OniaTrack(&_thePassedCands[i], iEvent, iSetup );
		//Onia+2track
		if(OniaDiTrackReco_)OniaDitrack( &_thePassedCands[i], iEvent, iSetup );

	}

	tree_->Fill();

	this->resetVariables();

}

void  
Analyzer::MCAna(const edm::Event& iEvent){
	// cout<<"Start analysis the "<<nEvt<<" event"<<endl;
	nGenKstar=0;

	Handle<GenParticleCollection> genparticles;
	iEvent.getByLabel( "genParticles", genparticles );

	for( size_t i = 0; i < genparticles->size(); ++ i ) {
		const GenParticle & p = (*genparticles)[i];

		if (abs(p.pdgId())==313 && ((p.mother())->pdgId() ==511)
				&& ((abs((p.daughter(0))->pdgId()) == 321 && abs((p.daughter(1))->pdgId()) == 211)
					||(abs((p.daughter(0))->pdgId()) == 211 && abs((p.daughter(1))->pdgId()) == 321))
			 ){

			Gen_Kstar_pdgid->push_back(p.pdgId());
			TLorentzVector v1(0.0,0.0,0.0,0.0);
			v1.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
			new((*Gen_Kstar_P4)[nGenKstar])TLorentzVector(v1);

			const Candidate * Kstarmom = p.mother();
			if ( Kstarmom != 0 ) {
				Gen_B_pdgid->push_back(Kstarmom->pdgId());
				TLorentzVector v2(0.0,0.0,0.0,0.0);
				v2.SetPxPyPzE(Kstarmom->px(),Kstarmom->py(),Kstarmom->pz(),Kstarmom->energy());
				new((*Gen_B_P4)[nGenKstar])TLorentzVector(v2);

				const Candidate * dauKstarGen0=p.daughter(0);
				const Candidate * dauKstarGen1=p.daughter(1);
				if(abs(dauKstarGen0->pdgId()) == 321){
					TLorentzVector v21(0.0,0.0,0.0,0.0);
					v21.SetPxPyPzE(dauKstarGen0->px(),dauKstarGen0->py(),dauKstarGen0->pz(),dauKstarGen0->energy());
					Gen_Kaon_pdgid->push_back(dauKstarGen0->pdgId());
					new((*Gen_Kaon_P4)[nGenKstar])TLorentzVector(v21);

					TLorentzVector v22(0.0,0.0,0.0,0.0);
					v22.SetPxPyPzE(dauKstarGen1->px(),dauKstarGen1->py(),dauKstarGen1->pz(),dauKstarGen1->energy());
					Gen_Pion_pdgid->push_back(dauKstarGen1->pdgId());
					new((*Gen_Pion_P4)[nGenKstar])TLorentzVector(v22);
				}
				else{
					TLorentzVector v23(0.0,0.0,0.0,0.0);
					v23.SetPxPyPzE(dauKstarGen0->px(),dauKstarGen0->py(),dauKstarGen0->pz(),dauKstarGen0->energy());
					Gen_Pion_pdgid->push_back(dauKstarGen0->pdgId());
					new((*Gen_Pion_P4)[nGenKstar])TLorentzVector(v23);

					TLorentzVector v24(0.0,0.0,0.0,0.0);
					v24.SetPxPyPzE(dauKstarGen1->px(),dauKstarGen1->py(),dauKstarGen1->pz(),dauKstarGen1->energy());
					Gen_Kaon_pdgid->push_back(dauKstarGen1->pdgId());
					new((*Gen_Kaon_P4)[nGenKstar])TLorentzVector(v24);
				}

				TVector3 b_prodPos,b_endPos;
				b_prodPos.SetXYZ( Kstarmom->vx(), Kstarmom->vy(), Kstarmom->vz());
				b_endPos.SetXYZ( p.vx(), p.vy(), p.vz());
				TVector3 vGenxyz = b_endPos - b_prodPos;   
				new((*Gen_B_V3)[nGenKstar])TVector3(vGenxyz);

				for(size_t k = 0; k < Kstarmom->numberOfDaughters(); ++k) {
					if (abs((Kstarmom->daughter(k))->pdgId())==313) continue;
					const Candidate * dauBGen   = Kstarmom->daughter(k);
					if(abs(dauBGen->pdgId())==14) continue;
					if(dauBGen->pdgId()==13){
						TLorentzVector v3(0.0,0.0,0.0,0.0);
						v3.SetPxPyPzE(dauBGen->px(),dauBGen->py(),dauBGen->pz(),dauBGen->energy());
						new((*Gen_muonPos_P4)[nGenKstar])TLorentzVector(v3);
					}
					else if(dauBGen->pdgId()==-13){
						TLorentzVector v4(0.0,0.0,0.0,0.0);
						v4.SetPxPyPzE(dauBGen->px(),dauBGen->py(),dauBGen->pz(),dauBGen->energy());
						new((*Gen_muonNeg_P4)[nGenKstar])TLorentzVector(v4);
					}
					else if(abs(dauBGen->pdgId()==443)){
						TLorentzVector v5(0.0,0.0,0.0,0.0);
						v5.SetPxPyPzE(dauBGen->px(),dauBGen->py(),dauBGen->pz(),dauBGen->energy());
						new((*Gen_Jpsi_P4)[nGenKstar])TLorentzVector(v5);
					}
				}

			}
			nGenKstar++;
		}
	}

} 

void 
Analyzer::hltReport(const edm::Event & iEvent ,const edm::EventSetup& iSetup){

	Handle< edm::TriggerResults> hltresults;
	try {
		iEvent.getByLabel( tagTriggerResults_, hltresults );
	} catch ( ... ) {
		cout << "Couldn't get handle on HLT Trigger!" << endl;
	}
	if (hltresults.isValid()) {

		int ntrigs=hltresults->size();
		const  edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);

		int nDimuTrigger = HLTBitNames_DoubleMu.size();
		if (nDimuTrigger != int(HLTLastFilterNames_DoubleMu.size())) cout <<"trigger size not equal filter size"<<endl;

		for (int i=0; i<nDimuTrigger; i++){
			DimuTriggerResult[i]= 0;
		}

		for (int it=0; it< ntrigs; it++) {
			string alltrigName = triggerNames_.triggerName(it);
			int hltflag = (*hltresults)[it].accept();
			itriggerflag->push_back(hltflag);
			itriggerNames->push_back(alltrigName);
			for (int itrig=0; itrig<nDimuTrigger; itrig++){
				if (HLTBitNames_DoubleMu[itrig] == alltrigName){

					DimuTriggerResult[itrig]= hltflag;
				}
			}
		}

		for (int itrig=0; itrig<nDimuTrigger; itrig++){

			Dimutriggerflag->push_back(DimuTriggerResult[itrig]);
			DimutriggerNames->push_back(HLTBitNames_DoubleMu[itrig]);
		}
	}else cout << "TriggerResults NOT valid in current event!" << endl;
}


void 
Analyzer::Onia2MuMu(const edm::Event& iEvent ,const edm::EventSetup& iSetup) {
	nonia=0;
	nOniaMu=0;
	nXcand=0;
	nOniaTrack=0;

	try {
		iEvent.getByLabel(muons_, allmuons);
	} catch ( cms::Exception& ex ) {
		edm::LogError("RecoAnalyzer") <<"Error! can't get collection with label : muons";
	}

	vector<pat::Muon>::const_iterator mu1;
	vector<pat::Muon>::const_iterator mu2;   
	vector<pat::Muon>::const_iterator muonP;
	vector<pat::Muon>::const_iterator muonM;

	if( allmuons->size()<2 ) return;
	ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

	for ( mu1 = allmuons->begin(); mu1 != allmuons->end(); mu1++) {
		// both must pass low quality
		if(!lowerPuritySelection_(*mu1)) continue; 
		for( mu2 =mu1+1; mu2 != allmuons->end(); mu2++ ) {
			if (mu1->charge() + mu2->charge() != 0) continue;
			// both must pass low quality
			if(!lowerPuritySelection_(*mu2)) continue; 
			// one must pass tight quality
			if (!(higherPuritySelection_(*mu1) || higherPuritySelection_(*mu2))) continue;

			const pat::Muon* patMuonP(0);const pat::Muon* patMuonM(0);
			if(mu1->charge() > 0){ patMuonP = &(*mu1); patMuonM = &(*mu2);muonP=mu1;muonM=mu2;}
			else if(mu1->charge() < 0){ patMuonP = &(*mu2); patMuonM = &(*mu1);muonP=mu2;muonM=mu1;}

			const reco::Muon *rmuP = dynamic_cast<const reco::Muon *>(muonP->originalObject());
			const reco::Muon *rmuM = dynamic_cast<const reco::Muon *>(muonM->originalObject());


			pat::CompositeCandidate myCand;
			myCand.addDaughter(*muonP,"muonP");
			myCand.addDaughter(*muonM,"muonM");

			myCand.addUserFloat("mupkey",rmuP->track().key());
			myCand.addUserFloat("mumkey",rmuM->track().key());

			// ---- define and set candidate's 4momentum  ----  
			typedef Candidate::LorentzVector LorentzVector;
			LorentzVector jpsi = muonP->p4() + muonM->p4();
			myCand.setP4(jpsi);
			myCand.setCharge(muonP->charge()+muonM->charge());

			// ---- apply the dimuon cut ----
			if(!dimuonSelection_(myCand)) continue;

			bool _selection( 0 );
			int thePassedCats(-1);

			if (patMuonP->isGlobalMuon() && patMuonM->isGlobalMuon() &&
					patMuonP->isTrackerMuon() && patMuonM->isTrackerMuon()   ) {
				if (!_applycuts || (selGlobalMuon(patMuonP) &&
							selGlobalMuon(patMuonM))){
					thePassedCats=0;
					_selection=1;
				}
			}
			if (patMuonP->isGlobalMuon() && patMuonM->isTrackerMuon() &&
					patMuonP->isTrackerMuon()&& !patMuonM->isGlobalMuon() ) {
				if (!_applycuts || (selGlobalMuon(patMuonP) &&
							selTrackerMuon(patMuonM))) {
					thePassedCats=1;
					_selection=1;
				}
			}

			if (patMuonM->isGlobalMuon() && patMuonP->isTrackerMuon() &&
					patMuonM->isTrackerMuon()&& !patMuonP->isGlobalMuon() ) {
				if (!_applycuts || (selGlobalMuon(patMuonM) &&
							selTrackerMuon(patMuonP))) {
					thePassedCats=1;
					_selection=1;
				}
			}

			// tracker + tracker?  
			if (patMuonP->isTrackerMuon() && patMuonM->isTrackerMuon()&&
					!patMuonP->isGlobalMuon() && !patMuonM->isGlobalMuon()) {
				if (!_applycuts || (selTrackerMuon(patMuonP) &&
							selTrackerMuon(patMuonM))) {
					thePassedCats=2;
					_selection=1;
				}
			}
			if(!_selection) continue;

			TrackRef muPtrack=muonP->track();
			TrackRef muMtrack=muonM->track();
			if (muPtrack.isNull() || muMtrack.isNull() ) continue;

			edm::ESHandle<TransientTrackBuilder> theB;
			iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

			TransientTrack muonPTT   = (*theB).build(&(*muPtrack));
			TransientTrack muonMTT   = (*theB).build(&(*muMtrack));

			KinematicParticleFactoryFromTransientTrack pFactory;

			//The mass of a muon and the insignificant mass sigma 
			//to avoid singularities in the covariance matrix.
			ParticleMass muon_mass = 0.10565837; //pdg mass
			float muon_sigma = muon_mass*1.e-6;

			//initial chi2 and ndf before kinematic fits.
			float chi = 0.;
			float ndf = 0.;
			vector<RefCountedKinematicParticle> muonParticles;
			muonParticles.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
			muonParticles.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));

			KinematicParticleVertexFitter fitter;   
			RefCountedKinematicTree psiVertexFitTree;
			psiVertexFitTree = fitter.fit(muonParticles); 


			//-----------------------------
			if (!psiVertexFitTree->isValid()) {
				//	   std::cout << "caught an exception in the psi vertex fit" << std::endl;
				continue; 
			}
			psiVertexFitTree->movePointerToTheTop();
			RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

			if (!psi_vFit_vertex_noMC->vertexIsValid()){
				// cout << "onia fit vertex is not valid" << endl;
				continue;
			}

			double VtxCL=ChiSquaredProbability((double)(psi_vFit_vertex_noMC->chiSquared()),(double)(psi_vFit_vertex_noMC->degreesOfFreedom()));
			if ( VtxCL<0.01 ) continue;

			RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
			KinematicParameters psiP=psi_vFit_noMC->currentState().kinematicParameters();

			psiVertexFitTree->movePointerToTheFirstChild();
			RefCountedKinematicParticle muPCandMC = psiVertexFitTree->currentParticle();
			psiVertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle muMCandMC = psiVertexFitTree->currentParticle();

			KinematicParameters psiMupKP = muPCandMC->currentState().kinematicParameters();
			KinematicParameters psiMumKP = muMCandMC->currentState().kinematicParameters();

			OniaVtxCL->push_back( VtxCL );
			OniaCats->push_back(thePassedCats);

			TLorentzVector jpsiv, muv_p,muv_m;   
			jpsiv.SetPxPyPzE(psiP.momentum().x(), psiP.momentum().y(), psiP.momentum().z(), psiP.energy() );
			muv_p.SetPxPyPzE(psiMupKP.momentum().x(), psiMupKP.momentum().y(), psiMupKP.momentum().z(), psiMupKP.energy() );
			muv_m.SetPxPyPzE(psiMumKP.momentum().x(), psiMumKP.momentum().y(), psiMumKP.momentum().z(), psiMumKP.energy() );

			new((*OniaP4)[nonia])TLorentzVector(jpsiv);
			new((*MuPP4)[nonia])TLorentzVector(muv_p);
			new((*MuMP4)[nonia])TLorentzVector(muv_m);

			_thePassedCands.push_back(myCand );
			myCand.addUserFloat("Onia_index",nonia);

			int ntriggers =  HLTBitNames_DoubleMu.size();

			for (int MatchTrig=0; MatchTrig < ntriggers; MatchTrig++){

				//cout << "result " << DimuTriggerResult[MatchTrig] << " " << HLTBitNames_DoubleMu[MatchTrig] <<endl;
				if (DimuTriggerResult[MatchTrig]!=0){
					const pat::TriggerObjectStandAloneCollection mu1HLTMatches = patMuonP->triggerObjectMatchesByFilter(HLTLastFilterNames_DoubleMu[MatchTrig]  );
					const pat::TriggerObjectStandAloneCollection mu2HLTMatches = patMuonM->triggerObjectMatchesByFilter(HLTLastFilterNames_DoubleMu[MatchTrig]  );
					bool pass1 = mu1HLTMatches.size() > 0;
					bool pass2 = mu2HLTMatches.size() > 0;

					//cout << "Uno " << pass1 << " " << HLTBitNames_DoubleMu[MatchTrig] <<endl;
					//cout << "Due " << pass2 << " " << HLTBitNames_DoubleMu[MatchTrig] <<endl;
					if ((pass1) && (pass2))
						JPsiMuonTrigMatch->push_back(true);

					else
						JPsiMuonTrigMatch->push_back(false);
				}
				else
					JPsiMuonTrigMatch->push_back(false);
			}

			nonia++;

		}
	}
}	  


void  
Analyzer::OniaTrack(const pat::CompositeCandidate* aCand,const edm::Event& iEvent, const edm::EventSetup& iSetup){
	const ParticleMass muon_mass = 0.10565837; //pdg mass
	float muon_sigma = muon_mass*1.e-6;
	const ParticleMass pion_mass = 0.13957018;
	float pion_sigma = pion_mass*1.e-6;



	ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

	const pat::Muon* muonPos = dynamic_cast<const pat::Muon*>(aCand->daughter("muonP"));
	const pat::Muon* muonNeg = dynamic_cast<const pat::Muon*>(aCand->daughter("muonM"));

	TrackRef TrackP = muonPos->track();
	TrackRef TrackN = muonNeg->track();

	TransientTrack track_MuP= (*theB).build(&TrackP);
	TransientTrack track_MuN= (*theB).build(&TrackN);

	//------------------------------------------
	Handle< vector<pat::GenericParticle> > allTracks;
	iEvent.getByLabel(_trackLabelPi, allTracks);
	vector<pat::GenericParticle>::const_iterator iTrack1;

	for (iTrack1=allTracks->begin(); iTrack1 != allTracks->end(); ++iTrack1){
		if (iTrack1->track().key()==aCand->userFloat("mupkey") || iTrack1->track().key()==aCand->userFloat("mumkey")) continue;
		//check track doesn't overlap with the Onia candidate tracks
		if ((iTrack1->track()->chi2()/ iTrack1->track()->ndof()> 10)|| iTrack1->pt() <0.8) continue;


		double OP1DR = deltaR( aCand->eta(), aCand->phi(), iTrack1->eta(),iTrack1->phi() );
		if(OP1DR > 3 ) continue;

		math::XYZTLorentzVector Xp4=aCand->p4()+iTrack1->p4();

		if(Xp4.M() >=20 || Xp4.M() <0) continue;

		TrackRef trk1Ref = iTrack1->track();         
		TransientTrack track_T1= (*theB).build(&trk1Ref);
		if (!selTrack( trk1Ref )) continue;

		//Creating a KinematicParticleFactory

		KinematicParticleFactoryFromTransientTrack pFactory;

		ParticleMass onia_mass=0;
		if ( aCand->mass()>2.8&&aCand->mass()<3.35 ){
			onia_mass = 3.096916;
		}else if ( aCand->mass()>3.5&&aCand->mass()<3.9){
			onia_mass = 3.68609;
		}else if ( aCand->mass()>9.1&&aCand->mass()<9.75 ){
			onia_mass = 9.4603;
		}else if ( aCand->mass()>9.75&&aCand->mass()<10.2 ){
			onia_mass = 10.02326;
		}else if ( aCand->mass()>10.2&&aCand->mass()<10.7 ){
			onia_mass = 10.3552;
		}

		float chi = 0.;
		float ndf = 0.;
		vector<RefCountedKinematicParticle> XParticles;
		XParticles.push_back(pFactory.particle(track_MuP,muon_mass,chi,ndf,muon_sigma));
		XParticles.push_back(pFactory.particle(track_MuN,muon_mass,chi,ndf,muon_sigma));
		XParticles.push_back(pFactory.particle(track_T1,pion_mass,chi,ndf,pion_sigma));
		MultiTrackKinematicConstraint *  onia_c = new  TwoTrackMassKinematicConstraint(onia_mass);
		KinematicConstrainedVertexFitter kcvFitter;
		RefCountedKinematicTree vertexFitTree = kcvFitter.fit(XParticles, onia_c);

		if (!vertexFitTree->isValid()) {
			//std::cout << "caught an exception in the X vertex fit with MC" << std::endl;
			continue;
		}
		vertexFitTree->movePointerToTheTop();
		RefCountedKinematicParticle xCandMC = vertexFitTree->currentParticle();
		RefCountedKinematicVertex xDecayVertexMC = vertexFitTree->currentDecayVertex();
		if (!xDecayVertexMC->vertexIsValid()){
			//cout << "X MC fit vertex is not valid" << endl;
			continue;
		}
		double VtxCL=ChiSquaredProbability((double)(xDecayVertexMC->chiSquared()),(double)(xDecayVertexMC->degreesOfFreedom()));
		if ( VtxCL<0.01 ) continue;


		if ( xDecayVertexMC->chiSquared()<0 || xDecayVertexMC->chiSquared()>10000 ) {
			// cout << " failed chi2 cut in MC fit with chi2 = " << xDecayVertexMC->chiSquared() << endl;
			continue;
		}

		PiD0->push_back(trk1Ref->dxy(RefVtx));      
		PiDz->push_back(trk1Ref->dz(RefVtx));   
		PiNHits->push_back(trk1Ref->numberOfValidHits());
		PiPixelHits->push_back(trk1Ref->hitPattern().numberOfValidPixelHits());
		PiNormChi2->push_back(trk1Ref->chi2()/trk1Ref->ndof());


		vertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle muP =  vertexFitTree->currentParticle();
		vertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle muN = vertexFitTree->currentParticle();
		vertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle pi1 = vertexFitTree->currentParticle();

		TLorentzVector LV_pi1(pi1->currentState().globalMomentum().x(),pi1->currentState().globalMomentum().y(), 
				pi1->currentState().globalMomentum().z(),pi1->currentState().kinematicParameters().energy());
		TLorentzVector LV_P(xCandMC->currentState().globalMomentum().x(),xCandMC->currentState().globalMomentum().y(), 
				xCandMC->currentState().globalMomentum().z(),xCandMC->currentState().kinematicParameters().energy());
		new((*PiP4)[nOniaTrack])TLorentzVector(LV_pi1);
		new((*BP4)[nOniaTrack])TLorentzVector(LV_P);

		BOindex->push_back(aCand->userFloat("Onia_index"));


		BVtxCL->push_back(VtxCL);
		BVtxC2->push_back( xDecayVertexMC->chiSquared() );

		//--distance--
		TVector3 vtx;
		TVector3 pvtx;
		vtx.SetXYZ(xDecayVertexMC->position().x(),xDecayVertexMC->position().y(),0);      
		pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
		TVector3 vdiff = vtx - pvtx;
		//--pt--
		TVector3 pperp(LV_P.Px(),LV_P.Px(),0);
		AlgebraicVector vpperp(3);
		vpperp[0] = pperp.x();
		vpperp[1] = pperp.y();
		vpperp[2] = 0.;

		VertexDistanceXY vdistXY;
		double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
		Measurement1D distXY = vdistXY.distance(Vertex(*xDecayVertexMC), thePrimaryV);
		double lxy=distXY.value();
		double lxyerr=distXY.error();
		double lxySig=lxy/lxyerr;

		double ctauPV =lxy*cosAlpha * xCandMC->currentState().mass()/pperp.Perp();
		GlobalError v1e = (Vertex(*xDecayVertexMC)).error();
		GlobalError v2e = thePrimaryV.error();
		AlgebraicSymMatrix vXYe = v1e.matrix()+ v2e.matrix();
		double ctauErrPV = sqrt(vXYe.similarity(vpperp))* xCandMC->currentState().mass()/(pperp.Perp2());
		double lxyPV = vdiff.Dot(pperp)/pperp.Mag();
		BcosAlpha->push_back(cosAlpha);
		BlxySig->push_back(lxySig);
		BctauPV->push_back( ctauPV );
		BctauErrPV->push_back( ctauErrPV );
		BlxyPV->push_back( lxyPV );
		nOniaTrack++;
	}
}

void  
Analyzer::OniaMu(const pat::CompositeCandidate* aCand,const edm::Event& iEvent, const edm::EventSetup& iSetup){

	//Jpsi mass cut
	// if (aCand->mass()>4.05&&aCand->mass()<2.8) return;
	ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

	const pat::Muon* muonPos = dynamic_cast<const pat::Muon*>(aCand->daughter("muonP"));
	const pat::Muon* muonNeg = dynamic_cast<const pat::Muon*>(aCand->daughter("muonM"));

	TrackRef TrackP = muonPos->track();
	TrackRef TrackN = muonNeg->track();

	TransientTrack track_MuP= (*theB).build(&TrackP);
	TransientTrack track_MuN= (*theB).build(&TrackN);

	//------------------------------------------
	iEvent.getByLabel(muons_, allmuons);
	vector<pat::Muon>::const_iterator muonk;
	if( allmuons->size()<3 ) return;
	for ( muonk = allmuons->begin(); muonk != allmuons->end(); muonk++) {
		if(!lowerPuritySelection_(*muonk)) continue;
		//overlap check
		const reco::Muon *rmu3 = dynamic_cast<const reco::Muon *>(muonk->originalObject());
		if (rmu3->track().key()==aCand->userFloat("mupkey") || rmu3->track().key()==aCand->userFloat("mumkey")) continue;
		const pat::Muon* patMuonk= &(*muonk);
		int thirdmuonCat_;
		bool _selection( 0 );
		if (muonk->isGlobalMuon()&&muonk->isTrackerMuon() &&selGlobalMuon(patMuonk)) { 
			thirdmuonCat_ = 1; 
			_selection=1;
		}
		else if (muonk->isTrackerMuon() && !muonk->isGlobalMuon()&&selTrackerMuon(patMuonk)){  
			thirdmuonCat_ = 2;
			_selection=1;
		}
		if (!_selection) continue; //Muon quallity cut


		// LorentzVector Trimup4 =aCand->p4() + muonk->p4();
		//   if(Trimup4.M() >=3 || Trimup4.M() <20) continue;

		// if (!OniaMuselection_(OniaMuCand)) continue;

		double OMuDR = deltaR(aCand->eta(), aCand->phi(), muonk->eta(), muonk->phi());
		if(OMuDR > 3.14) continue;

		TrackRef mu3track=muonk->track();
		if ( mu3track.isNull()) continue;

		TransientTrack muonMTT= (*theB).build(&(*mu3track));

		//Creating a KinematicParticleFactory

		KinematicParticleFactoryFromTransientTrack pFactory;

		ParticleMass muon_mass = 0.1056583;
		float muon_sigma = muon_mass*1.e-6;
		float chi = 0.;
		float ndf = 0.;
		float jpsiMsigma = 0.00004;

		ParticleMass onia_mass=0;
		if ( aCand->mass()>2.8&&aCand->mass()<3.35 ){
			onia_mass = 3.096916;
		}else if ( aCand->mass()>3.5&&aCand->mass()<3.9){
			onia_mass = 3.68609;
		}else if ( aCand->mass()>9.1&&aCand->mass()<9.75 ){
			onia_mass = 9.4603;
		}else if ( aCand->mass()>9.75&&aCand->mass()<10.2 ){
			onia_mass = 10.02326;
		}else if ( aCand->mass()>10.2&&aCand->mass()<10.7 ){
			onia_mass = 10.3552;
		}

		//making particles
		vector<RefCountedKinematicParticle> allParticles;

		allParticles.push_back(pFactory.particle (track_MuP,muon_mass,chi,ndf,muon_sigma));
		allParticles.push_back(pFactory.particle (track_MuN,muon_mass,chi,ndf,muon_sigma));

		KinematicParticleVertexFitter Fitter;
		RefCountedKinematicTree JpsiTree = Fitter.fit(allParticles);

		if (!JpsiTree->isValid()) {
			//std::cout << "caught an exception in the psi vertex fit" << std::endl;
			continue; 
		}
		KinematicParticleFitter constFitter;
		KinematicConstraint * jpsi_const = new MassKinematicConstraint(onia_mass,jpsiMsigma);
		JpsiTree = constFitter.fit(jpsi_const,JpsiTree);
		if (!JpsiTree->isValid()) {
			// std::cout << "caught an exception in the psi mass constraint fit" << std::endl;
			continue; 
		}

		JpsiTree->movePointerToTheTop();

		RefCountedKinematicParticle Jpsi_branch = JpsiTree->currentParticle();
		vector<RefCountedKinematicParticle> BParticles;

		BParticles.push_back(pFactory.particle (muonMTT,muon_mass,chi,ndf,muon_sigma));
		BParticles.push_back(Jpsi_branch);
		RefCountedKinematicTree vertexFitTree = Fitter.fit(BParticles);

		if (!vertexFitTree->isValid()) {
			//std::cout << "caught an exception in the B fit" << std::endl;
			continue; 
		}

		vertexFitTree->movePointerToTheTop();
		RefCountedKinematicParticle BCandMC = vertexFitTree->currentParticle();
		RefCountedKinematicVertex bVertex = vertexFitTree->currentDecayVertex();

		if (!bVertex->vertexIsValid()) continue;

		//std::cout << "B kinematic fit successful "<<std::endl;

		if ( bVertex->chiSquared()<0 || bVertex->chiSquared()>10000 ) {
			continue;
		}
		double VtxCL=ChiSquaredProbability(BCandMC->chiSquared(), (int)BCandMC->degreesOfFreedom());
		if(VtxCL<0.01)continue;//vertex proberbility preselection

		Mu3D0->push_back(mu3track->dxy(RefVtx));      
		Mu3Dz->push_back(mu3track->dz(RefVtx));   
		Mu3NHits->push_back(mu3track->numberOfValidHits());
		Mu3PixelHits->push_back(mu3track->hitPattern().numberOfValidPixelHits());
		Mu3NormChi2->push_back(mu3track->chi2()/mu3track->ndof());

		vertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle muP =  vertexFitTree->currentParticle();
		vertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle muN = vertexFitTree->currentParticle();
		vertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle mu3 = vertexFitTree->currentParticle();

		TLorentzVector LV_mu3(mu3->currentState().globalMomentum().x(),mu3->currentState().globalMomentum().y(), 
				mu3->currentState().globalMomentum().z(),mu3->currentState().kinematicParameters().energy());
		TLorentzVector LV_P(BCandMC->currentState().globalMomentum().x(),BCandMC->currentState().globalMomentum().y(), 
				BCandMC->currentState().globalMomentum().z(),BCandMC->currentState().kinematicParameters().energy());
		new((*Mu3P4)[nOniaMu])TLorentzVector(LV_mu3);
		new((*TrMP4)[nOniaMu])TLorentzVector(LV_P);

		TrMOindex->push_back(aCand->userFloat("Onia_index"));
		TrMCats->push_back(thirdmuonCat_);

		TrMVtxCL->push_back(VtxCL);
		TrMVtxC2->push_back( bVertex->chiSquared() );

		//--distance--
		TVector3 vtx;
		TVector3 pvtx;
		vtx.SetXYZ(bVertex->position().x(),bVertex->position().y(),0);      
		pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
		TVector3 vdiff = vtx - pvtx;
		//--pt--
		TVector3 pperp(LV_P.Px(),LV_P.Px(),0);
		AlgebraicVector vpperp(3);
		vpperp[0] = pperp.x();
		vpperp[1] = pperp.y();
		vpperp[2] = 0.;

		VertexDistanceXY vdistXY;
		double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
		Measurement1D distXY = vdistXY.distance(Vertex(*bVertex), thePrimaryV);
		double lxy=distXY.value();
		double lxyerr=distXY.error();
		double lxySig=lxy/lxyerr;

		double ctauPV =lxy*cosAlpha * BCandMC->currentState().mass()/pperp.Perp();
		GlobalError v1e = (Vertex(*bVertex)).error();
		GlobalError v2e = thePrimaryV.error();
		AlgebraicSymMatrix vXYe = v1e.matrix()+ v2e.matrix();
		double ctauErrPV = sqrt(vXYe.similarity(vpperp))* BCandMC->currentState().mass()/(pperp.Perp2());
		double lxyPV = vdiff.Dot(pperp)/pperp.Mag();
		TrMcosAlpha->push_back(cosAlpha);
		TrMlxySig->push_back(lxySig);
		TrMctauPV->push_back( ctauPV );
		TrMctauErrPV->push_back( ctauErrPV );
		TrMlxyPV->push_back( lxyPV );
		nOniaMu++;

	}
}

void 
Analyzer::OniaDitrack(const pat::CompositeCandidate* aCand, const edm::Event& iEvent ,const edm::EventSetup& iSetup){
	const ParticleMass muon_mass = 0.10565837; //pdg mass
	float muon_sigma = muon_mass*1.e-6;
	const ParticleMass pion_mass = 0.13957018;
	float pion_sigma = pion_mass*1.e-6;
	const ParticleMass kaon_mass = 0.493677;
	float kaon_sigma = kaon_mass*1.e-6;
	const ParticleMass kstar_mass = 0.89610000;
	float kstar_sigma = kstar_mass*1.e-6;
	//cout<<">>Debug001"<<endl;

	//Jpsi mass cut
	ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);	

	const pat::Muon* muonPos = dynamic_cast<const pat::Muon*>(aCand->daughter("muonP"));
	const pat::Muon* muonNeg = dynamic_cast<const pat::Muon*>(aCand->daughter("muonM"));

	TrackRef TrackP = muonPos->track();
	TrackRef TrackN = muonNeg->track();

	TransientTrack track_MuP= (*theB).build(&TrackP);
	TransientTrack track_MuN= (*theB).build(&TrackN);

	//----------------------------------------------
	Handle< vector<pat::GenericParticle> > allTracks;
	iEvent.getByLabel(_trackLabelPi, allTracks);
	vector<pat::GenericParticle>::const_iterator iTrack1;
	vector<pat::GenericParticle>::const_iterator iTrack2;

	//cout<<">>allTracks->size(): "<<allTracks->size()<<endl;
	for (iTrack1=allTracks->begin(); iTrack1 != allTracks->end(); ++iTrack1){
		if (iTrack1->track().key()==aCand->userFloat("mupkey") || iTrack1->track().key()==aCand->userFloat("mumkey")) continue;
		//check track doesn't overlap with the Onia candidate tracks
		if ((iTrack1->track()->chi2()/ iTrack1->track()->ndof()> 10)|| iTrack1->pt() < PiPt_c) continue;

		for (iTrack2=iTrack1+1; iTrack2 != allTracks->end(); ++iTrack2){

			if(iTrack1->charge()+iTrack2->charge() != 0) continue;

			if ((iTrack2->track()->chi2()/ iTrack2->track()->ndof() > 10)|| iTrack2->pt() < PiPt_c)continue;
			if (iTrack2->track().key()==aCand->userFloat("mupkey") || iTrack2->track().key()==aCand->userFloat("mumkey")) continue;
			// cout <<"mupkey"<<aCand->userFloat("mupkey")<<endl;

			//check track doesn't overlap with the Onia candidate tracks

			double OP1DR = deltaR( aCand->eta(), aCand->phi(), iTrack1->eta(),iTrack1->phi() );
			double OP2DR = deltaR( aCand->eta(), aCand->phi(), iTrack2->eta(),iTrack2->phi() );


			if(OP1DR > OniaPiDR_c || OP2DR > OniaPiDR_c ) continue;

			math::XYZTLorentzVector Xp4=aCand->p4()+iTrack1->p4() + iTrack2->p4();
			math::XYZTLorentzVector Ruop4=iTrack1->p4() + iTrack2->p4();


			//if (((aCand->mass()>4.05||aCand->mass()<2.8)||(Xp4.M() >=5.5||Xp4.M() <0))&&
			//     ((aCand->mass()>11.5||aCand->mass()<8.5)||(Xp4.M() >=13||Xp4.M() <9))
			// )continue; 
			//if(fabs(Ruop4.M()-8.9610000e-01) > 3*0.08) continue;


			TrackRef trk1Ref = iTrack1->track();
			TrackRef trk2Ref = iTrack2->track();

			TransientTrack track_T1= (*theB).build(&trk1Ref);
			TransientTrack track_T2= (*theB).build(&trk2Ref);
			if ( !selTrack( trk1Ref ) || !selTrack( trk2Ref ) ) continue;

			KinematicParticleFactoryFromTransientTrack pFactory;

			float chi = 0.;
			float ndf = 0.;

			ParticleMass onia_mass=0.;
			if ( aCand->mass()>2.8&&aCand->mass()<3.35 ){
				onia_mass = 3.096916;
			}else if ( aCand->mass()>3.5&&aCand->mass()<3.9){
				onia_mass = 3.68609;
			}else if ( aCand->mass()>9.1&&aCand->mass()<9.75 ){
				onia_mass = 9.4603;
			}else if ( aCand->mass()>9.75&&aCand->mass()<10.2 ){
				onia_mass = 10.02326;
			}else if ( aCand->mass()>10.2&&aCand->mass()<10.7 ){
				onia_mass = 10.3552;
			}
			//double Qvalue=Xp4.M()-onia_mass-Ruop4.M();
			//if (Qvalue>0.5) continue;

			vector<RefCountedKinematicParticle> KstarParticles;
			KstarParticles.push_back(pFactory.particle(track_T1,pion_mass,chi,ndf,pion_sigma));
			KstarParticles.push_back(pFactory.particle(track_T2,kaon_mass,chi,ndf,kaon_sigma));
			KinematicParticleVertexFitter KstarCVFitter;
			RefCountedKinematicTree KstarfitTree=KstarCVFitter.fit(KstarParticles);
			if(!KstarfitTree->isValid()) continue;

			KstarfitTree->movePointerToTheTop();
			KinematicParticleFitter csFitterKstar;
			KinematicConstraint * kstar_c = new MassKinematicConstraint(kstar_mass,kstar_sigma);
			KstarfitTree=csFitterKstar.fit(kstar_c,KstarfitTree);

			if(!KstarfitTree->isValid()) continue;
			KstarfitTree->movePointerToTheTop();
			KstarfitTree->movePointerToTheFirstChild();
			RefCountedKinematicParticle track1=KstarfitTree->currentParticle();
			KstarfitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle track2=KstarfitTree->currentParticle();
			TLorentzVector LVtrack1(track1->currentState().globalMomentum().x(),track1->currentState().globalMomentum().y(), 
					track1->currentState().globalMomentum().z(),track1->currentState().kinematicParameters().energy());
			TLorentzVector LVtrack2(track2->currentState().globalMomentum().x(),track2->currentState().globalMomentum().y(), 
					track2->currentState().globalMomentum().z(),track2->currentState().kinematicParameters().energy());
			TLorentzVector LVKstar(LVtrack1+LVtrack2);
		
			//cout<<">>Debug002"<<endl;
			//cout<<">>fabs(LVKstar.M()-8.9610000e-01): "<<fabs(LVKstar.M()-8.9610000e-01)<<endl;
			if( !(fabs(LVKstar.M()-8.9610000e-01)<3.*0.08) ) continue;
			//cout<<">>Debug003"<<endl;

			vector<RefCountedKinematicParticle> XParticles;
			XParticles.push_back(pFactory.particle(track_MuP,muon_mass,chi,ndf,muon_sigma));
			XParticles.push_back(pFactory.particle(track_MuN,muon_mass,chi,ndf,muon_sigma));
			XParticles.push_back(pFactory.particle(track_T1,pion_mass,chi,ndf,pion_sigma));
			XParticles.push_back(pFactory.particle(track_T2,kaon_mass,chi,ndf,kaon_sigma));
			MultiTrackKinematicConstraint *  onia_c = new  TwoTrackMassKinematicConstraint(onia_mass);
			KinematicConstrainedVertexFitter kcvFitter;
			RefCountedKinematicTree vertexFitTree = kcvFitter.fit(XParticles, onia_c);

			if (!vertexFitTree->isValid()) {
				 //std::cout << "caught an exception in the X vertex fit with MC" << std::endl;
				continue;
			}
			//cout<<">>Debug004"<<endl;
			vertexFitTree->movePointerToTheTop();
			RefCountedKinematicParticle xCandMC = vertexFitTree->currentParticle();
			RefCountedKinematicVertex xDecayVertexMC = vertexFitTree->currentDecayVertex();
			if (!xDecayVertexMC->vertexIsValid()){
				  //cout << "X MC fit vertex is not valid" << endl;
				continue;
			}

			if ( xDecayVertexMC->chiSquared()<0 || xDecayVertexMC->chiSquared()>10000 ) {
				 //cout << " failed chi2 cut in MC fit with chi2 = " << xDecayVertexMC->chiSquared() << endl;
				continue;
			}
			//cout<<">>Debug005"<<endl;
			double VtxCL=ChiSquaredProbability(xDecayVertexMC->chiSquared(), (int)xDecayVertexMC->degreesOfFreedom());
			if(VtxCL<0.001)continue;//vertex proberbility preselection
			//cout<<">>Debug006"<<endl;
			Pi1D0->push_back(trk1Ref->dxy(RefVtx));      
			Pi1Dz->push_back(trk1Ref->dz(RefVtx));   
			Pi1NHits->push_back(trk1Ref->numberOfValidHits());
			Pi1PixelHits->push_back(trk1Ref->hitPattern().numberOfValidPixelHits());
			Pi1NormChi2->push_back(trk1Ref->chi2()/trk1Ref->ndof());

			Pi2D0->push_back(trk2Ref->dxy(RefVtx)) ;      
			Pi2Dz->push_back(trk2Ref->dz(RefVtx)) ;   
			Pi2NHits->push_back(trk2Ref->numberOfValidHits()) ;
			Pi2PixelHits->push_back(trk2Ref->hitPattern().numberOfValidPixelHits());
			Pi2NormChi2->push_back(trk2Ref->chi2()/trk2Ref->ndof()) ;



			vertexFitTree->movePointerToTheFirstChild();
			RefCountedKinematicParticle muP =  vertexFitTree->currentParticle();
			vertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle muN = vertexFitTree->currentParticle();
			vertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle pi1 = vertexFitTree->currentParticle();
			vertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle pi2 = vertexFitTree->currentParticle();

			TLorentzVector LV_pi1(pi1->currentState().globalMomentum().x(),pi1->currentState().globalMomentum().y(), 
					pi1->currentState().globalMomentum().z(),pi1->currentState().kinematicParameters().energy());
			TLorentzVector LV_pi2(pi2->currentState().globalMomentum().x(),pi2->currentState().globalMomentum().y(), 
					pi2->currentState().globalMomentum().z(),pi2->currentState().kinematicParameters().energy());

			xOindex->push_back(aCand->userFloat("Onia_index"));
			TLorentzVector LV_P(xCandMC->currentState().globalMomentum().x(),xCandMC->currentState().globalMomentum().y(), 
					xCandMC->currentState().globalMomentum().z(),xCandMC->currentState().kinematicParameters().energy());
			new((*Pi1P4)[nXcand])TLorentzVector(LV_pi1);
			new((*Pi2P4)[nXcand])TLorentzVector(LV_pi2);
			new((*KstarP4)[nXcand])TLorentzVector(LV_pi1+LV_pi2);
			new((*xP4)[nXcand])TLorentzVector(LV_P);

			xVtxCL->push_back( VtxCL);
			xVtxC2->push_back( xDecayVertexMC->chiSquared() );

			//--distance--
			TVector3 vtx;
			TVector3 pvtx;
			vtx.SetXYZ(xDecayVertexMC->position().x(),xDecayVertexMC->position().y(),0);      
			pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
			TVector3 vdiff = vtx - pvtx;
			//--pt--
			TVector3 pperp(LV_P.Px(),LV_P.Px(),0);
			AlgebraicVector vpperp(3);
			vpperp[0] = pperp.x();
			vpperp[1] = pperp.y();
			vpperp[2] = 0.;

			VertexDistanceXY vdistXY;
			double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
			Measurement1D distXY = vdistXY.distance(Vertex(*xDecayVertexMC), thePrimaryV);
			double lxy=distXY.value();
			double lxyerr=distXY.error();
			double lxySig=lxy/lxyerr;

			double ctauPV =lxy*cosAlpha * xCandMC->currentState().mass()/pperp.Perp();
			GlobalError v1e = (Vertex(*xDecayVertexMC)).error();
			GlobalError v2e = thePrimaryV.error();
			AlgebraicSymMatrix vXYe = v1e.matrix()+ v2e.matrix();
			double ctauErrPV = sqrt(vXYe.similarity(vpperp))* xCandMC->currentState().mass()/(pperp.Perp2());
			double lxyPV = vdiff.Dot(pperp)/pperp.Mag();
			xcosAlpha->push_back(cosAlpha);
			xlxySig->push_back(lxySig);
			xctauPV->push_back( ctauPV );
			xctauErrPV->push_back( ctauErrPV );
			xlxyPV->push_back( lxyPV );
			nXcand++;
		}
	}
}

bool
Analyzer::selGlobalMuon(const pat::Muon* aMuon) {

	TrackRef iTrack = aMuon->innerTrack();
	const reco::HitPattern& p = iTrack->hitPattern();
	const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
	const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

	TrackRef gTrack = aMuon->globalTrack();
	const reco::HitPattern& q = gTrack->hitPattern();

	bool trackOK = false;
	// cooler way of cutting on tracks
	if (_applyExpHitcuts) {
		float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
		trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
		// old way of cutting on tracks  
	} else trackOK = (iTrack->found() > 10);

	return (// isMuonInAccept(aMuon) &&
			trackOK &&
			gTrack->chi2()/gTrack->ndof() < 20.0 &&
			q.numberOfValidMuonHits() > 0 &&
			iTrack->chi2()/iTrack->ndof() < 1.8 &&
			aMuon->muonID("TrackerMuonArbitrated") &&
			aMuon->muonID("TMOneStationTight") &&
			p.pixelLayersWithMeasurement() > 1 &&
			fabs(iTrack->dxy(RefVtx)) < 3.0 &&
			fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool 
Analyzer::selTrackerMuon(const pat::Muon* aMuon) {

	TrackRef iTrack = aMuon->innerTrack();
	const reco::HitPattern& p = iTrack->hitPattern();
	const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
	const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

	bool trackOK = false;
	// cooler way of cutting on tracks
	if (_applyExpHitcuts) {
		float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
		trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
		// old way of cutting on tracks  
	} else trackOK = (iTrack->found() > 10);

	return (// isMuonInAccept(aMuon) &&
			trackOK &&
			iTrack->chi2()/iTrack->ndof() < 1.8 &&
			aMuon->muonID("TrackerMuonArbitrated") &&
			aMuon->muonID("TMOneStationTight") &&
			p.pixelLayersWithMeasurement() > 1 &&
			fabs(iTrack->dxy(RefVtx)) < 3.0 &&
			fabs(iTrack->dz(RefVtx)) < 15.0 );
}
bool 
Analyzer::selTrack(TrackRef atrkRef) {
	return (atrkRef->dxy(RefVtx)<3&&atrkRef->numberOfValidHits()>4
			&&atrkRef->hitPattern().numberOfValidPixelHits()>1
			&&atrkRef->chi2()/atrkRef->ndof()<5);   
}
bool
Analyzer::isMuonInAccept(const pat::Muon* aMuon) {
	// *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
	return (fabs(aMuon->eta()) < 2.4 &&
			((fabs(aMuon->eta()) < 1.3 && aMuon->pt() > 3.3) ||
			 (fabs(aMuon->eta()) > 1.3 && fabs(aMuon->eta()) < 2.2 && aMuon->p() > 2.9) ||
			 (fabs(aMuon->eta()) > 2.2 && aMuon->pt() > 0.8)));

	// *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
	// by just returning TRUE
	//  return true;
}



//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
