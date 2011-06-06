#ifndef __jpsi_polarization__
#define __jpsi_polarization__

#include "TLorentzVector.h"
#include "TMath.h"

namespace jpsi{

  //enum {L1DoubleMuOpen, Mu0Track0Jpsi, Mu3Track0Jpsi, DoubleMu0, DoubleMu3};

  // beam energy in GeV
  const double pbeam = 3500.;
  // masses
  const double Mprot = 0.9382720;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );
  const double muMass = 0.105658;
  //rap bins
  Int_t const kNbRapForPTBins = 5;
  /* Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.5, 1.9, 2.3}; */
  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
  //pT bins (optimized to have at least 10.000 entries per bin)
  Int_t const kNbPTMaxBins = 12;
  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 7,8,9,12,12};//all y, y1, y2, y3, y4, y5
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},//all rapidities
    {0., 6., 7., 8., 10., 15., 20., 30.},//mid-rap
    {0., 4., 6., 7., 8., 10., 15., 20., 30.},
    {0., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.}};//most forward
  //the following values are dummy and need to be filled out at end of analysis
  Double_t pTWCentre_rap[kNbRapForPTBins+1][23] = 
    {{0.5, 1.25, 1.75, 1.95, 2.5, 2.8, 3.1, 3.4, 3.9, 4.3, 4.6, 5.0, 5.5, 6.0, 6.5, 7.2, 7.8, 8.5, 9.6, 11.0, 14., 27.},
     {7., 12.},
     {4.5, 6.0, 7.3, 15.},
     {1.0, 2.5, 3.5, 4.5, 5.5, 6.7, 10.},
     {0.5, 1.5, 2.5, 3.5, 5.0, 13.},
     {0.5, 1.5, 2.5, 3.5, 5.0, 13.}};
  //need to be extracted from the real data:
  /* Double_t pTWCentre[kNbPTBins] = {1.527, 3.561, 6.021, 9.540, 15.205, 23.566}; */
  //number of reference frames
  Int_t const kNbFrames = 6;
  Char_t *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  //polarization variables (3rd one was for testing purposes)
  Int_t const kNbPolVar = 3; //cosTheta, phi, cos2Phi
  enum {cosThPol,phiPol,cos2PhiPol};
  //cosTheta 
  Int_t const kNbBinsCosT = 40;
  Double_t cosTMin = -1., cosTMax = 1.;
  //phi for pol. 
  Int_t const kNbBinsPhiPol = 36;
  //  Double_t phiPolMin = 0., phiPolMax = 360.;
  Double_t phiPolMin = -180., phiPolMax = 180.;
  //cos2Phi
  Int_t const kNbBinsCos2Phi = 40;
  Double_t cos2PhiMin = -1., cos2PhiMax = 1.;

  //study the negative and positive rapidity sides separately
  Int_t const kNbRapBins = kNbRapForPTBins;
  /* Double_t rapRange[2*kNbRapBins+1] = {-2.3, -1.9, -1.5, -0.9, 0., 0.9, 1.5, 1.9, 2.3}; */
  Double_t rapRange[2*kNbRapBins+1] = {-2.4, -2.1, -1.6, -1.2, -0.9, 0., 0.9, 1.2, 1.6, 2.1, 2.4};
  
  //phase space limiting cuts:
  Int_t const kNbEtaRegions = 3;
  Double_t etaPS[kNbEtaRegions] = {1.3, 2.2, 2.4}; //pseudo-rap cuts for muons
  Double_t pTMuMin[kNbEtaRegions] = {3.3, 0., 0.8};
  Double_t pMuMin[kNbEtaRegions] = {0., 2.9, 0.};
  Double_t rapYPS = 2.4;
  /* Double_t JpsiCtauMax = 0.100; //100 micron */
  Double_t JpsiCtauMax = 1000.; //effectively no cut on lifetime
  Double_t nSigMass = 2.;
  Double_t nSigBkgLow = 3.5;
  Double_t nSigBkgHigh = 3.5;
  Double_t polMassJpsi[kNbRapForPTBins+1] = {3.092, 3.094, 3.094, 3.092, 3.092, 3.090};//[all rap, rap bin 1-3]
  Double_t sigmaMassJpsi[kNbRapForPTBins+1] = {0.042, 0.024, 0.035, 0.040, 0.045, 0.056};//[all rap, rap bin 1-3]

  //some make up to use the same colour and marker for each pT and rapidity bin
  //in every plotting macro:
  Int_t colour_pT[kNbPTMaxBins+1] = {1, 2, 3, 4, 6, 7, 8, 49, 38, 46, 12, 40};
  Int_t marker_pT[kNbPTMaxBins+1] = {20, 21, 25, 22, 23, 26, 27, 28, 29, 30, 20, 20};
  /* Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 6, 7, 7, 6, 4, 3, 2}; */
  /* Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 29, 20, 24, 30, 23, 25, 24}; */
  /* Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 4, 3, 2}; */
  /* Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 23, 25, 24}; */
  
  Int_t colour_rapForPTBins[kNbRapForPTBins+1] = {1, 30, 4, 2, 3, kMagenta+1};
  Int_t marker_rapForPTBins[kNbRapForPTBins+1] = {20, 21, 25, 20, 22, 29};

  //min number of entries to consider the bin
  //(in acceptance map, etc)
  Int_t minEntriesPerBin = 10;
}

#endif
