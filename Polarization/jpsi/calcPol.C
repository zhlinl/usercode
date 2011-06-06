#include "commonVar.h"
Double_t thisCosTh[jpsi::kNbFrames], thisPhi[jpsi::kNbFrames];
Double_t thisPhi_rad[jpsi::kNbFrames], thisCosPhi[jpsi::kNbFrames];

//=========================================
// calculation of decay angular parameters
//=========================================
void calcPol(TLorentzVector muplus_LAB,
	     TLorentzVector muminus_LAB){
  
  TLorentzVector qqbar_LAB = muplus_LAB + muminus_LAB;
  Double_t rapidity = qqbar_LAB.Rapidity();

  // boost beams and positive muon into the q-qbar rest frame:
  TVector3 LAB_to_QQBAR = -qqbar_LAB.BoostVector();

  TLorentzVector beam1_QQBAR = jpsi::beam1_LAB;
  beam1_QQBAR.Boost( LAB_to_QQBAR );

  TLorentzVector beam2_QQBAR = jpsi::beam2_LAB;
  beam2_QQBAR.Boost( LAB_to_QQBAR );

  TLorentzVector muplus_QQBAR = muplus_LAB;
  muplus_QQBAR.Boost( LAB_to_QQBAR );

  // reference directions in the Jpsi rest frame:

  TVector3 beam1_direction     = beam1_QQBAR.Vect().Unit();
  TVector3 beam2_direction     = beam2_QQBAR.Vect().Unit();
  TVector3 qqbar_direction     = qqbar_LAB.Vect().Unit();
  TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

  // all polarization frames have the same Y axis = the normal to the plane formed by
  // the directions of the colliding hadrons
  TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();
  if ( rapidity < 0. ) Yaxis = -Yaxis; //H: added (5 Dec 2010)

  /////////////////////////////////////////////////////////////////////
  // CS frame

  TVector3 newZaxis = beam1_beam2_bisect;
  TVector3 newYaxis = Yaxis;
  TVector3 newXaxis = newYaxis.Cross( newZaxis );

  TRotation rotation;
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();   // transforms coordinates from the "xyz" system
  // to the "new" (rotated) system having the polarization axis
  // as z axis

  TVector3 muplus_QQBAR_rotated(muplus_QQBAR.Vect());
  
  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[jpsi::CS] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[jpsi::CS] = muplus_QQBAR_rotated.Phi();
  thisPhi[jpsi::CS] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[jpsi::CS] < 0. ) thisPhi[jpsi::CS]= 360. + thisPhi[jpsi::CS];      // phi defined in degrees from 0 to 360
  //  thisPhi[jpsi::CS] += 180.; //H: don't add anything...

  /////////////////////////////////////////////////////////////////////
  // HELICITY frame

  newZaxis = qqbar_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[jpsi::HX] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[jpsi::HX] = muplus_QQBAR_rotated.Phi();
  thisPhi[jpsi::HX] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[jpsi::HX] < 0. ) thisPhi[jpsi::HX] = 360. + thisPhi[jpsi::HX]; // phi defined in degrees from 0 to 360
  //thisPhi[jpsi::HX] += 180.;//H: don't add anything...
 
  /////////////////////////////////////////////////////////////////////
  // GJ1 frame

  newZaxis = beam1_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[jpsi::GJ1] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[jpsi::GJ1] = muplus_QQBAR_rotated.Phi();
  thisPhi[jpsi::GJ1] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[jpsi::GJ1] < 0. ) thisPhi[jpsi::GJ1] = 360. + thisPhi[jpsi::GJ1]; // phi defined in degrees from 0 to 360
  //thisPhi[jpsi::GJ1] += 180.;//H: don't add anything...

  /////////////////////////////////////////////////////////////////////
  // GJ2 frame

  newZaxis = beam2_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[jpsi::GJ2] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[jpsi::GJ2] = muplus_QQBAR_rotated.Phi();
  thisPhi[jpsi::GJ2] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[jpsi::GJ2] < 0. ) thisPhi[jpsi::GJ2] = 360. + thisPhi[jpsi::GJ2]; // phi defined in degrees from 0 to 360
  //thisPhi[jpsi::GJ2] += 180.;//H: don't add anything...
  /////////////////////////////////////////////////////////////////////
  // sGJ frame (symmetrized GJ)

  newZaxis = beam1_direction; if( rapidity < 0. ) newZaxis = beam2_direction;
  newYaxis = Yaxis;

  // try to swith the following line on or off
  //if( rapidity < 0. ) newYaxis = -Yaxis;

  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[jpsi::sGJ] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[jpsi::sGJ] = muplus_QQBAR_rotated.Phi();
  thisPhi[jpsi::sGJ] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[jpsi::sGJ] < 0. ) thisPhi[jpsi::sGJ] = 360. + thisPhi[jpsi::sGJ]; // phi defined in degrees from 0 to 360
  //thisPhi[jpsi::sGJ] += 180.;//H: don't add anything...

  /////////////////////////////////////////////////////////////////////
  // PHX frame ("perpendicular helicity frame" - z axis perpendicular
  // to the CS axis)

  newZaxis = newZaxis = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[jpsi::PHX] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[jpsi::PHX] = muplus_QQBAR_rotated.Phi();
  thisPhi[jpsi::PHX] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[jpsi::PHX] < 0. ) thisPhi[jpsi::PHX] = 360. + thisPhi[jpsi::PHX]; // phi defined in degrees from 0 to 360
  //thisPhi[jpsi::PHX] += 180.;//H: don't add anything...
}
