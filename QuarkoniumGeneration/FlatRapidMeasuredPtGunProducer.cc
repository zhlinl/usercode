/*
 *  $Date: 2012/03/22 13:10:56 $
 *  $Revision: 1.1 $
 *  \author Jean-Roch Vlimant
 */

#include <ostream>

#include "IOMC/ParticleGuns/interface/FlatRapidMeasuredPtGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// #include "FWCore/Utilities/interface/Exception.h"

// #include "CLHEP/Random/RandExpo.h"

using namespace edm;
using namespace std;

FlatRapidMeasuredPtGunProducer::FlatRapidMeasuredPtGunProducer(const ParameterSet& pset) : 
   BaseFlatGunProducer(pset)
{


   ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters") ;
  
   fMinPt = pgun_params.getParameter<double>("MinPt");
   fMaxPt = pgun_params.getParameter<double>("MaxPt");

   fMinRapidity = pgun_params.getUntrackedParameter<double>("MinRapidity",-0.1);
   fMaxRapidity = pgun_params.getUntrackedParameter<double>("MaxRapidity",0.1);

	 fResonance = pgun_params.getUntrackedParameter<double>("Resonance",1);

   produces<HepMCProduct>();
   produces<GenEventInfoProduct>();

   
}


FlatRapidMeasuredPtGunProducer::~FlatRapidMeasuredPtGunProducer()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}

double FlatRapidMeasuredPtGunProducer::getRandom(double MinPt, double MaxPt){

	//realistic pT distribution function
	TF1 *fPT = new TF1("fPT", "[0]*x*pow(1.+(1./([1]-2.))*x*x/[2],-[1])", MinPt, MaxPt);

	if(fResonance == 1){
		fPT->FixParameter(0, 0.1); //Ups(1S) from BPH-11-001
		fPT->FixParameter(1, 3.46);//Ups(1S) from BPH-11-001
		fPT->FixParameter(2, 47.3);//Ups(1S) from BPH-11-001
	}
	else if(fResonance == 2){
		fPT->FixParameter(0, 0.1); //Ups(2S) from BPH-11-001
		fPT->FixParameter(1, 3.27);//Ups(2S) from BPH-11-001
		fPT->FixParameter(2, 65.7);//Ups(2S) from BPH-11-001
	}
	else if(fResonance == 3){
		fPT->FixParameter(0, 0.1); //Ups(3S)
		fPT->FixParameter(1, 3.05);//Ups(3S)
		fPT->FixParameter(2, 80.5);//Ups(3S)
	}
	else if(fResonance == 0){
		fPT->FixParameter(0, 0.1); //Jpsi
		fPT->FixParameter(1, 3.69);//Jpsi
		fPT->FixParameter(2, 12.0);//Jpsi
	}
	else if(fResonance == 4){
		fPT->FixParameter(0, 0.1); //Psi'
		fPT->FixParameter(1, 3.71);//Psi'
		fPT->FixParameter(2, 19.54);//Psi'
	}
	else{
		return fRandomGenerator->fire(MinPt, MaxPt); //flat pT 
	}

	return fPT->GetRandom();
}

void FlatRapidMeasuredPtGunProducer::produce(Event &e, const EventSetup& es) 
{

	if ( fVerbosity > 0 )
	{
		cout << " FlatRapidMeasuredPtGunProducer : Begin New Event Generation" << endl ; 
	}
	// event loop (well, another step in it...)

	// no need to clean up GenEvent memory - done in HepMCProduct
	// 

	// here re-create fEvt (memory)
	//
	fEvt = new HepMC::GenEvent() ;

	// now actualy, cook up the event from PDGTable and gun parameters
	//
	// 1st, primary vertex
	//
	//HepMC::GenVertex* Vtx = new HepMC::GenVertex(CLHEP::HepLorentzVector(0.,0.,0.));
	HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));

	// loop over particles
	//
	int barcode = 1 ;
	for (unsigned int ip=0; ip<fPartIDs.size(); ++ip)
	{

		//the max is to ensure you don't generate at 0
		//the 90% is to get rid of edge effect

		//double pt     =  std::max(0.00001,0.90*fMinPt)+fRandomExpoGenerator->fire(fMeanPt);
		//shoot until in the designated range
		//while (pt<fMinPt || pt>fMaxPt)
		//  {pt = std::max(0.00001,0.90*fMinPt) + fRandomExpoGenerator->fire(fMeanPt);}


		//double pt = fRandomGenerator->fire(fMinPt,fMaxPt);
		double pt = getRandom(fMinPt,fMaxPt);

		double y  = fRandomGenerator->fire(fMinRapidity, fMaxRapidity) ;

		double phi    = fRandomGenerator->fire(fMinPhi, fMaxPhi) ;
		int PartID = fPartIDs[ip] ;
		double mass;
		if(abs(PartID)!= 200553) {
			const HepPDT::ParticleData* 
				PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
			mass   = PData->mass().value() ;
		} else {
			mass = 10.3552;
		}

		double energy = sqrt(pt*pt+mass*mass)*(exp(2*y)+1)/(2*exp(y));
		double pz = energy*(exp(2*y)-1)/(exp(2*y)+1);

		double px     = pt*cos(phi) ;
		double py     = pt*sin(phi) ;

		//CLHEP::Hep3Vector p(px,py,pz) ;
		//HepMC::GenParticle* Part = 
		//    new HepMC::GenParticle(CLHEP::HepLorentzVector(p,energy),PartID,1);
		HepMC::FourVector p(px,py,pz,energy) ;
		HepMC::GenParticle* Part = 
			new HepMC::GenParticle(p,PartID,1);
		Part->suggest_barcode( barcode ) ;
		barcode++ ;
		Vtx->add_particle_out(Part);

		if ( fAddAntiParticle )
		{
			//CLHEP::Hep3Vector ap(-px,-py,-pz) ;
			HepMC::FourVector ap(-px,-py,-pz,energy) ;
			int APartID = -PartID ;
			if ( PartID == 22 || PartID == 23 )
			{
				APartID = PartID ;
			}	  
			//HepMC::GenParticle* APart =
			//   new HepMC::GenParticle(CLHEP::HepLorentzVector(ap,energy),APartID,1);
			HepMC::GenParticle* APart =
				new HepMC::GenParticle(ap,APartID,1);
			APart->suggest_barcode( barcode ) ;
			barcode++ ;
			Vtx->add_particle_out(APart) ;
		}

	}

	fEvt->add_vertex(Vtx) ;
	fEvt->set_event_number(e.id().event()) ;
	fEvt->set_signal_process_id(20) ; 

	if ( fVerbosity > 0 )
	{
		fEvt->print() ;  
	}

	auto_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
	BProduct->addHepMCData( fEvt );
	e.put(BProduct);

	auto_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
	e.put(genEventInfo);

	if ( fVerbosity > 0 )
	{
		// for testing purpose only
		// fEvt->print() ; // prints empty info after it's made into edm::Event
		cout << " FlatRandomPtGunProducer : Event Generation Done " << endl;
	}
}

