#ifndef FlatRapidMeasuredPtGunProducer_H
#define FlatRapidMeasuredPtGunProducer_H

/** \class FlatRapidMeasuredPtGunProducer
 *
 * Generates single particle gun in HepMC format
 * Jean-Roch Vlimant
 ***************************************/

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"
#include "CLHEP/Random/RandExponential.h"
#include "TF1.h"
#include "TRandom3.h"
namespace edm
{
  
  class FlatRapidMeasuredPtGunProducer : public BaseFlatGunProducer
  {
  
  public:
    FlatRapidMeasuredPtGunProducer(const ParameterSet & pset);
    virtual ~FlatRapidMeasuredPtGunProducer();

  private:
   
    virtual void produce(Event & e, const EventSetup& es);
    
  protected :
  
    // data members
    
    double            fMinPt   ;
    double            fMaxPt   ;

    double            fMinRapidity;
    double            fMaxRapidity;
		double              fResonance;
		double getRandom(double MinPt, double MaxPt);


  };
} 

#endif
