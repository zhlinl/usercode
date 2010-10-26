#include "CLHEP/Random/DRand48Engine.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/Randomize.h"

#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/CommonAlignment/interface/AlignableComposite.h"
#include "Alignment/CommonAlignment/interface/AlignableModifier.h"


//__________________________________________________________________________________________________
AlignableModifier::AlignableModifier( void )
{

  theDRand48Engine = new CLHEP::DRand48Engine();

}


//__________________________________________________________________________________________________
AlignableModifier::~AlignableModifier()
{

  delete theDRand48Engine;

}

//__________________________________________________________________________________________________
void AlignableModifier::init_( void )
{

  // Initialize all known parameters (according to ORCA's MisalignmentScenario.cc)
  distribution_ = "";        // Switch for distributions ("fixed","flat","gaussian")
  setError_     = false;     // Apply alignment errors
  setRotations_ = true;      // Apply rotations
  setTranslations_ = true;   // Apply translations
  scale_        = 1.;        // Scale to apply to all movements
  scaleError_   = 1.;        // Scale to apply to alignment errors
  phiX_         = 0.;        // Rotation angle around X [rad]
  phiY_         = 0.;        // Rotation angle around Y [rad]
  phiZ_         = 0.;        // Rotation angle around Z [rad]
  phiXlocal_    = 0.;        // Local rotation angle around X [rad]
  phiYlocal_    = 0.;        // Local rotation angle around Y [rad]
  phiZlocal_    = 0.;        // Local rotation angle around Z [rad]
  dX_           = 0.;        // X displacement [cm]
  dY_           = 0.;        // Y displacement [cm]
  dZ_           = 0.;        // Z displacement [cm]
  dXlocal_      = 0.;        // Local X displacement [cm]
  dYlocal_      = 0.;        // Local Y displacement [cm]
  dZlocal_      = 0.;        // Local Z displacement [cm]
  twist_        = 0.;        // Twist angle [rad]
  shear_        = 0.;        // Shear angle [rad]

  // These are set through 'distribution'
  random_       = true;      // Use random distributions
  gaussian_     = true;      // Use gaussian distribution (otherwise flat)

}

//__________________________________________________________________________________________________
// Return true if given parameter name should be propagated down
const bool AlignableModifier::isPropagated( const std::string& parameterName ) const
{

  if ( parameterName == "distribution"    || 
       parameterName == "setError"        ||
       parameterName == "scaleError"      ||
       parameterName == "setRotations"    ||
       parameterName == "setTranslations" ||
       parameterName == "scale" 
       ) return true;
  
  return false;

}


//__________________________________________________________________________________________________
/// All known parameters and defaults are defined here! Returns true if modification actually applied.
bool AlignableModifier::modify( Alignable* alignable, const edm::ParameterSet& pSet )
{

  // Initialize parameters
  this->init_();
  int rotX_=0, rotY_=0, rotZ_=0; // To check correct backward compatibility


  // Reset counter
  m_modified = 0;
  
  // Retrieve parameters
  std::ostringstream error;
  std::vector<std::string> parameterNames = pSet.getParameterNames();
  for ( std::vector<std::string>::iterator iParam = parameterNames.begin(); 
        iParam != parameterNames.end(); iParam++ ) {
    if  ( (*iParam) == "distribution" ) distribution_ = pSet.getParameter<std::string>( *iParam );
    else if ( (*iParam) == "setError" ) setError_ = pSet.getParameter<bool>( *iParam );
    else if ( (*iParam) == "setRotations") setRotations_ = pSet.getParameter<bool>( *iParam );
    else if ( (*iParam) == "setTranslations") setTranslations_ = pSet.getParameter<bool>( *iParam );
    else if ( (*iParam) == "scale" )    scale_ = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "scaleError" ) scaleError_ = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "phiX" )    phiX_     = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "phiY" )    phiY_     = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "phiZ" )    phiZ_     = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "dX" )      dX_       = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "dY" )      dY_       = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "dZ" )      dZ_       = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "dXlocal" ) dXlocal_  = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "dYlocal" ) dYlocal_  = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "dZlocal" ) dZlocal_  = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "twist" )   twist_    = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "shear" )   shear_    = pSet.getParameter<double>( *iParam );
    else if ( (*iParam) == "localX" ) { phiXlocal_=pSet.getParameter<double>( *iParam ); rotX_++; }
    else if ( (*iParam) == "localY" ) { phiYlocal_=pSet.getParameter<double>( *iParam ); rotY_++; }
    else if ( (*iParam) == "localZ" ) { phiZlocal_=pSet.getParameter<double>( *iParam ); rotZ_++; }
    else if ( (*iParam) == "phiXlocal" ) { phiXlocal_=pSet.getParameter<double>( *iParam ); rotX_++; }
    else if ( (*iParam) == "phiYlocal" ) { phiYlocal_=pSet.getParameter<double>( *iParam ); rotY_++; }
    else if ( (*iParam) == "phiZlocal" ) { phiZlocal_=pSet.getParameter<double>( *iParam ); rotZ_++; }
    else if ( !pSet.existsAs<edm::ParameterSet>(*iParam) ) { // PSets are OK to ignore
      if ( !error.str().length() ) error << "Unknown parameter name(s): ";
      error << " " << *iParam;
    }
  }

  // Check if both 'localN' and 'phiNlocal' have been used
  if ( rotX_==2 ) throw cms::Exception("BadConfig") << "Found both localX and phiXlocal";
  if ( rotY_==2 ) throw cms::Exception("BadConfig") << "Found both localY and phiYlocal";
  if ( rotZ_==2 ) throw cms::Exception("BadConfig") << "Found both localZ and phiZlocal";

  // Check error
  if ( error.str().length() )
    throw cms::Exception("BadConfig") << error.str();

  // Decode distribution
  this->setDistribution( distribution_ );

  //if (scale_) { NO! Different random sequence if only parts scale to zero!

  // Apply displacements
  if ( std::abs(dX_) + std::abs(dY_) + std::abs(dZ_) > 0 && setTranslations_ )
    this->moveAlignable( alignable, random_, gaussian_, scale_*dX_, scale_*dY_, scale_*dZ_ );
  
  // Apply local displacements
  if ( std::abs(dXlocal_) + std::abs(dYlocal_) + std::abs(dZlocal_) > 0 && setTranslations_ )
    this->moveAlignableLocal( alignable, random_, gaussian_, 
			      scale_*dXlocal_, scale_*dYlocal_, scale_*dZlocal_ );
  
  // Apply rotations
  if ( std::abs(phiX_) + std::abs(phiY_) + std::abs(phiZ_) > 0 && setRotations_ )
    this->rotateAlignable( alignable, random_, gaussian_,
			   scale_*phiX_, scale_*phiY_, scale_*phiZ_ );
  
  // Apply local rotations
  if ( std::abs(phiXlocal_) + std::abs(phiYlocal_) + std::abs(phiZlocal_) > 0 && setRotations_ )
    this->rotateAlignableLocal( alignable, random_, gaussian_, 
				scale_*phiXlocal_, scale_*phiYlocal_, scale_*phiZlocal_ );
  
  // Apply twist
  if ( std::abs(twist_) > 0 )
    edm::LogError("NotImplemented") << "Twist is not implemented yet";
  
  // Apply shear
  if ( std::abs(shear_) > 0 )
    edm::LogError("NotImplemented") << "Shear is not implemented yet";
  
  // Apply error - first add scale_ to error
  scaleError_ *= scale_;
  if ( setError_ && scaleError_ ) {
    // Alignment Position Error for flat distribution: 1 sigma
    if ( !gaussian_ ) scaleError_ *= 0.68;
    
    
    // Error on displacement
    if ( std::abs(dX_) + std::abs(dY_) + std::abs(dZ_) > 0 && setTranslations_ )
      this->addAlignmentPositionError( alignable, 
				       scaleError_*dX_, scaleError_*dY_, scaleError_*dZ_ );
    
    // Error on local displacements
    if ( std::abs(dXlocal_) + std::abs(dYlocal_) + std::abs(dZlocal_) > 0 && setTranslations_ )
      this->addAlignmentPositionErrorLocal( alignable,
					    scaleError_*dXlocal_, scaleError_*dYlocal_, 
					    scaleError_*dZlocal_ );
    
    // Error on rotations
    if ( std::abs(phiX_) + std::abs(phiY_) + std::abs(phiZ_) > 0 && setRotations_ )
      this->addAlignmentPositionErrorFromRotation( alignable, 
						   scaleError_*phiX_, scaleError_*phiY_, 
						   scaleError_*phiZ_ );
    
    // Error on local rotations
    if ( std::abs(phiXlocal_) + std::abs(phiYlocal_) + std::abs(phiZlocal_) > 0
	 && setRotations_ )
      this->addAlignmentPositionErrorFromLocalRotation( alignable, 
							scaleError_*phiXlocal_, 
							scaleError_*phiYlocal_, 
							scaleError_*phiZlocal_ );
  }
  // } // end if (scale_)

  return ( m_modified > 0 );
  
}


//__________________________________________________________________________________________________
void AlignableModifier::setDistribution( const std::string& distr )
{

  if ( distr == "fixed" ) random_ = false;
  else if ( distr == "flat" ) {
    random_   = true;
    gaussian_ = false;
  } else if ( distr == "gaussian" ) {
    random_   = true;
    gaussian_ = true;
  }
  
}


//__________________________________________________________________________________________________
/// If 'seed' is zero, asks  RandomNumberGenerator service.
void AlignableModifier::setSeed( const long seed )
{

  long m_seed;

  if ( seed > 0 ) m_seed = seed;
  else {
    edm::Service<edm::RandomNumberGenerator> rng;
    m_seed = rng->mySeed();
  }

  LogDebug("PrintArgs") << "Setting generator seed to " << m_seed;

  theDRand48Engine->setSeed( m_seed );

}

//__________________________________________________________________________________________________
/// If 'random' is false, the given movements are strictly applied. Otherwise, a random
/// number is generated according to a gaussian or a flat distribution depending on 'gaussian'.
void AlignableModifier::moveAlignable( Alignable* alignable, bool random, bool gaussian,
                                       float sigmaX, float sigmaY, float sigmaZ )
{

  
  std::ostringstream message;
 
  // Get movement vector according to arguments
  GlobalVector moveV( sigmaX, sigmaY, sigmaZ ); // Default: fixed
  if ( random ) {
    std::vector<float> randomNumbers;
    message << "random ";
    if (gaussian)	{
      randomNumbers = this->gaussianRandomVector( sigmaX, sigmaY, sigmaZ );
      message << "gaussian ";
    } else 	{
      randomNumbers = this->flatRandomVector( sigmaX, sigmaY, sigmaZ );
      message << "flat ";
    }
    moveV = GlobalVector( randomNumbers[0], randomNumbers[1], randomNumbers[2] );
  }
  
  message << " move with sigma " << sigmaX << " " << sigmaY << " " << sigmaZ;

  LogDebug("PrintArgs") << message.str(); // Arguments

  LogDebug("PrintMovement") << "applied displacement: " << moveV; // Actual movements
  alignable->move(moveV);
  m_modified++;


}

//__________________________________________________________________________________________________
/// If 'random' is false, the given movements are strictly applied. Otherwise, a random
/// number is generated according to a gaussian or a flat distribution depending on 'gaussian'.
void AlignableModifier::moveAlignableLocal( Alignable* alignable, bool random, bool gaussian,
                                            float sigmaX, float sigmaY, float sigmaZ )
{

  
  std::ostringstream message;
 
  // Get movement vector according to arguments
  align::LocalVector moveV( sigmaX, sigmaY, sigmaZ ); // Default: fixed
  if ( random ) {
    std::vector<float> randomNumbers;
    message << "random ";
    if (gaussian) {
      randomNumbers = this->gaussianRandomVector( sigmaX, sigmaY, sigmaZ );
      message << "gaussian ";
    } else {
      randomNumbers = this->flatRandomVector( sigmaX, sigmaY, sigmaZ );
      message << "flat ";
    }
    moveV = align::LocalVector( randomNumbers[0], randomNumbers[1], randomNumbers[2] );
  }
  
  message << " move with sigma " << sigmaX << " " << sigmaY << " " << sigmaZ;

  LogDebug("PrintArgs") << message.str(); // Arguments

  LogDebug("PrintMovement") << "applied local displacement: " << moveV; // Actual movements
  alignable->move( alignable->surface().toGlobal(moveV) );
  m_modified++;


}


//__________________________________________________________________________________________________
/// If 'random' is false, the given rotations are strictly applied. Otherwise, a random
/// number is generated according to a gaussian or a flat distribution depending on 'gaussian'.
void AlignableModifier::rotateAlignable( Alignable* alignable, bool random, bool gaussian,
                                         float sigmaPhiX, float sigmaPhiY, float sigmaPhiZ )
{

  
  std::ostringstream message;

  // Get rotation vector according to arguments
  GlobalVector rotV( sigmaPhiX, sigmaPhiY, sigmaPhiZ ); // Default: fixed
  if ( random ) {
    std::vector<float> randomNumbers;
    message << "random ";
    if (gaussian) {
      randomNumbers = this->gaussianRandomVector( sigmaPhiX, sigmaPhiY, sigmaPhiZ );
      message << "gaussian ";
    } else {
      randomNumbers = flatRandomVector( sigmaPhiX, sigmaPhiY, sigmaPhiZ );
      message << "flat ";
    }
    rotV = GlobalVector( randomNumbers[0], randomNumbers[1], randomNumbers[2] );
  }
  
  message << "global rotation by angles " << sigmaPhiX << " " << sigmaPhiY << " " << sigmaPhiZ;

  LogDebug("PrintArgs") << message.str(); // Arguments

  LogDebug("PrintMovement") << "applied rotation angles: " << rotV; // Actual movements
  if ( std::abs(sigmaPhiX) ) alignable->rotateAroundGlobalX( rotV.x() );
  if ( std::abs(sigmaPhiY) ) alignable->rotateAroundGlobalY( rotV.y() );
  if ( std::abs(sigmaPhiZ) ) alignable->rotateAroundGlobalZ( rotV.z() );
  m_modified++;


}

//__________________________________________________________________________________________________
/// If 'random' is false, the given rotations are strictly applied. Otherwise, a random
/// number is generated according to a gaussian or a flat distribution depending on 'gaussian'.
void 
AlignableModifier::rotateAlignableLocal( Alignable* alignable, bool random, bool gaussian,
                                         float sigmaPhiX, float sigmaPhiY, float sigmaPhiZ )
{

  
  std::ostringstream message;

  // Get rotation vector according to arguments
  align::LocalVector rotV( sigmaPhiX, sigmaPhiY, sigmaPhiZ ); // Default: fixed
  if ( random ) {
    std::vector<float> randomNumbers;
    message << "random ";
    if (gaussian) {
      randomNumbers = this->gaussianRandomVector( sigmaPhiX, sigmaPhiY, sigmaPhiZ );
      message << "gaussian ";
    } else {
      randomNumbers = flatRandomVector( sigmaPhiX, sigmaPhiY, sigmaPhiZ );
      message << "flat ";
    }
    rotV = align::LocalVector( randomNumbers[0], randomNumbers[1], randomNumbers[2] );
  }
  
  message << "local rotation by angles " << sigmaPhiX << " " << sigmaPhiY << " " << sigmaPhiZ;
  
  LogDebug("PrintArgs") << message.str(); // Arguments
  
  LogDebug("PrintMovement") << "applied local rotation angles: " << rotV; // Actual movements
  if ( std::abs(sigmaPhiX) ) alignable->rotateAroundLocalX( rotV.x() );
  if ( std::abs(sigmaPhiY) ) alignable->rotateAroundLocalY( rotV.y() );
  if ( std::abs(sigmaPhiZ) ) alignable->rotateAroundLocalZ( rotV.z() );
  m_modified++;


}


//__________________________________________________________________________________________________
const std::vector<float> 
AlignableModifier::gaussianRandomVector( float sigmaX, float sigmaY, float sigmaZ ) const
{

  // Get absolute value if negative arguments
  if ( sigmaX < 0 ) {
    edm::LogWarning("BadConfig") << " taking absolute value for gaussian sigma_x";
    sigmaX = std::abs(sigmaX);
  }
  if ( sigmaY < 0 ) {
    edm::LogWarning("BadConfig") << " taking absolute value for gaussian sigma_y";
    sigmaY = std::abs(sigmaY);
  }
  if ( sigmaZ < 0 ) {
    edm::LogWarning("BadConfig") << " taking absolute value for gaussian sigma_z";
    sigmaZ = std::abs(sigmaZ);
  }

  // Pass by reference, otherwise pointer is deleted!
  CLHEP::RandGauss aGaussObjX( *theDRand48Engine, 0., sigmaX );
  CLHEP::RandGauss aGaussObjY( *theDRand48Engine, 0., sigmaY );
  CLHEP::RandGauss aGaussObjZ( *theDRand48Engine, 0., sigmaZ );

  std::vector<float> randomVector;
  randomVector.push_back( aGaussObjX.fire() );
  randomVector.push_back( aGaussObjY.fire() );
  randomVector.push_back( aGaussObjZ.fire() );

  return randomVector;

}


//__________________________________________________________________________________________________
const  std::vector<float> 
AlignableModifier::flatRandomVector( float sigmaX,float sigmaY, float sigmaZ ) const
{

  // Get absolute value if negative arguments
  if ( sigmaX < 0 ) {
    edm::LogWarning("BadConfig") << " taking absolute value for flat sigma_x";
    sigmaX = std::abs(sigmaX);
  }
  if ( sigmaY < 0 ) {
    edm::LogWarning("BadConfig") << " taking absolute value for flat sigma_y";
    sigmaY = std::abs(sigmaY);
  }
  if ( sigmaZ < 0 ) {
    edm::LogWarning("BadConfig") << " taking absolute value for flat sigma_z";
    sigmaZ = std::abs(sigmaZ);
  }

  CLHEP::RandFlat aFlatObjX( *theDRand48Engine, -sigmaX, sigmaX );
  CLHEP::RandFlat aFlatObjY( *theDRand48Engine, -sigmaY, sigmaY );
  CLHEP::RandFlat aFlatObjZ( *theDRand48Engine, -sigmaZ, sigmaZ );

  std::vector<float> randomVector;
  randomVector.push_back( aFlatObjX.fire() );
  randomVector.push_back( aFlatObjY.fire() );
  randomVector.push_back( aFlatObjZ.fire() );

  return randomVector;

}



//__________________________________________________________________________________________________
void AlignableModifier::addAlignmentPositionError( Alignable* alignable, 
                                                   float dx, float dy, float dz )
{

  LogDebug("PrintArgs") << "Adding an AlignmentPositionError of size " 
                        << dx << " "  << dy << " "  << dz;

  AlignmentPositionError ape(dx,dy,dz);
  alignable->addAlignmentPositionError( ape, true );

}


//__________________________________________________________________________________________________
void AlignableModifier::addAlignmentPositionErrorLocal( Alignable* alignable, 
                                                        float dx, float dy, float dz )
{

  LogDebug("PrintArgs") << "Adding a local AlignmentPositionError of size " 
                        << dx << " "  << dy << " "  << dz;

  AlgebraicSymMatrix as(3,0); //3x3, zeroed
  as[0][0] = dx*dx; as[1][1] = dy*dy; as[2][2] = dz*dz; //diagonals
  align::RotationType rt = alignable->globalRotation(); //get rotation
  AlgebraicMatrix am(3,3);
  am[0][0]=rt.xx(); am[0][1]=rt.xy(); am[0][2]=rt.xz();
  am[1][0]=rt.yx(); am[1][1]=rt.yy(); am[1][2]=rt.yz();
  am[2][0]=rt.zx(); am[2][1]=rt.zy(); am[2][2]=rt.zz();
  as=as.similarityT(am); //rotate error matrix

  GlobalError ge( as );
  AlignmentPositionError ape( ge );

  alignable->addAlignmentPositionError( ape, true ); // propagate down to components

}



//__________________________________________________________________________________________________
void AlignableModifier::addAlignmentPositionErrorFromRotation( Alignable* alignable, 
                                                               float phiX, float phiY, 
                                                               float phiZ )
{

  align::RotationType rotx( Basic3DVector<float>(1.0, 0.0, 0.0), phiX );
  align::RotationType roty( Basic3DVector<float>(0.0, 1.0, 0.0), phiY );
  align::RotationType rotz( Basic3DVector<float>(0.0, 0.0, 1.0), phiZ );
  align::RotationType rot = rotz * roty * rotx;
  
  this->addAlignmentPositionErrorFromRotation( alignable, rot );

}


//__________________________________________________________________________________________________
void AlignableModifier::addAlignmentPositionErrorFromLocalRotation( Alignable* alignable, 
                                                                    float phiX, float phiY, 
                                                                    float phiZ )
{

  align::RotationType rotx( Basic3DVector<float>(1.0, 0.0, 0.0), phiX );
  align::RotationType roty( Basic3DVector<float>(0.0, 1.0, 0.0), phiY );
  align::RotationType rotz( Basic3DVector<float>(0.0, 0.0, 1.0), phiZ );
  align::RotationType rot = rotz * roty * rotx;
  
  this->addAlignmentPositionErrorFromLocalRotation( alignable, rot );

}


//__________________________________________________________________________________________________
void AlignableModifier::addAlignmentPositionErrorFromRotation( Alignable* alignable, 
                                                               align::RotationType& rotation )
{ 

  LogDebug("PrintArgs") << "Adding an AlignmentPositionError from Rotation" << std::endl 
                        << rotation;

  alignable->addAlignmentPositionErrorFromRotation(rotation, true); // propagate down to components

}


//__________________________________________________________________________________________________
void AlignableModifier::addAlignmentPositionErrorFromLocalRotation( Alignable* alignable, 
                                                                    align::RotationType& rotation )
{ 
  
  LogDebug("PrintArgs") << "Adding an AlignmentPositionError from Local Rotation" << std::endl 
                        << rotation;
  
  // true: propagate down to components
  alignable->addAlignmentPositionErrorFromLocalRotation(rotation, true);
  
}

