#include "TLorentzVector.h"
#include "TMath.h"
import ROOT
from ROOT import TLorentzVector, TMath
from math import sqrt

#enum {L1DoubleMuOpen, Mu0Track0Jpsi, Mu3Track0Jpsi, DoubleMu0, DoubleMu3};

# beam energy in GeV
pbeam = 3500.
# masses
Mprot = 0.9382720
Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot )
beam1_LAB = TLorentzVector ( 0., 0., pbeam, Ebeam )
beam2_LAB = TLorentzVector ( 0., 0., -pbeam, Ebeam )
muMass = 0.105658
#rap bins
kNbRapForPTBins = 5
# Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.5, 1.9, 2.3}; */
rapForPTRange = [[0.,0.9],
                 [0.9,1.2],
                 [1.2,1.6],
                 [1.6,2.1],
                 [2.1,2.4]]
#pT bins (optimized to have at least 10.000 entries per bin)
kNbPTMaxBins = 12;
kNbPTBins = [kNbPTMaxBins, 7,8,9,12,12] #all y, y1, y2, y3, y4, y5
pTRange = [
    [[0.,1.],[1.,2.],[2.,3.],[3.,4.],[4.,5.], [5.,6.], [6.,7.], [7.,8.], [8.,10.], [10.,15.], [15.,20.], [20.,30.]], #all rapidities
    [[0.,6.],[6.,7.],[7.,8.],[8.,10.],[10.,15.],[15.,20.],[20.,30.]],  #mid-rap  
    [[0.,4.],[4.,6.],[6.,7.],[7.,8.],[8.,10.],[10.,15.],[15.,20.],[20.,30.]]] # 
#    [[4.,5.],[5.,6.],[6.,7.],[7.,8.],[8.,10.],[10.,15.],[15.,20.],[20.,30.]], # [0.,4.],
#    [[0.,1.],[1.,2.],[2.,3.],[3.,4.],[4.,5.], [5.,6.], [6.,7.], [7.,8.], [8.,10.], [10.,15.], [15.,20.], [20.,30.]],
#    [[0.,1.],[1.,2.],[2.,3.],[3.,4.],[4.,5.], [5.,6.], [6.,7.], [7.,8.], [8.,10.], [10.,15.], [15.,20.], [20.,30.]]] #//most forward
#the following values are dummy and need to be filled out at end of analysis
pTWCentre_rap = [
    [0.5, 1.25, 1.75, 1.95, 2.5, 2.8, 3.1, 3.4, 3.9, 4.3, 4.6, 5.0, 5.5, 6.0, 6.5, 7.2, 7.8, 8.5, 9.6, 11.0, 14., 27.],
    [7., 12.],
    [4.5, 6.0, 7.3, 15.],
    [1.0, 2.5, 3.5, 4.5, 5.5, 6.7, 10.],
    [0.5, 1.5, 2.5, 3.5, 5.0, 13.],
    [0.5, 1.5, 2.5, 3.5, 5.0, 13.]]
#need to be extracted from the real data:
# Double_t pTWCentre[kNbPTBins] = {1.527, 3.561, 6.021, 9.540, 15.205, 23.566}; */
#number of reference frames
kNbFrames = 6
frameLabel = ["CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"]
CS, HX, PHX, sGJ, GJ1, GJ2 = range(kNbFrames)
#polarization variables (3rd one was for testing purposes)
kNbPolVar = 3 # //cosTheta, phi, cos2Phi
cosThPol,phiPol,cos2PhiPol = range(kNbPolVar) # equivalent to enum {cosThPol, etc...}
#cosTheta 
kNbBinsCosT = 40
cosTMin = -1.0
cosTMax = 1.0
#phi for pol. 
kNbBinsPhiPol = 36
phiPolMin = 0.
phiPolMax = 360.
#cos2Phi
kNbBinsCos2Phi = 40
cos2PhiMin = -1.0
cos2PhiMax = 1.0

#study the negative and positive rapidity sides separately
kNbRapBins = kNbRapForPTBins;
# Double_t rapRange[2*kNbRapBins+1] = {-2.3, -1.9, -1.5, -0.9, 0., 0.9, 1.5, 1.9, 2.3}; */
rapRange = [-2.4, -2.1, -1.6, -1.2, -0.9, 0., 0.9, 1.2, 1.6, 2.1, 2.4]
  
#phase space limiting cuts:
etaPS = [1.3,2.2,2.4] # //pseudo-rap cuts for muons
pTMuMin = [3.3,0.0,0.8]
pMuMin = [0.0,2.9,0.0]
rapYPS = 2.4
# JpsiCtauMax = 0.100 #//100 micron */
JpsiCtauMax = 1000. #effectively no cut on lifetime
nSigMass = 3.5
nSigBkgLow = 4.0
nSigBkgHigh = 4.0
polMassJpsi = [3.092, 3.094, 3.094, 3.092, 3.092, 3.090]  #//[all rap, rap bin 1-3]
sigmaMassJpsi = [0.042, 0.024, 0.035, 0.040, 0.045, 0.056] #//[all rap, rap bin 1-3]

#some make up to use the same colour and marker for each pT and rapidity bin
#in every plotting macro:
# Int_t colour_pT[kNbPTBins+1] = {1, 2, 3, 4, 6, 7, 8}; */
# Int_t marker_pT[kNbPTBins+1] = {20, 21, 22, 23, 20, 21, 29}; */
# Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 6, 7, 7, 6, 4, 3, 2}; */
# Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 29, 20, 24, 30, 23, 25, 24}; */
# Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 4, 3, 2}; */
# Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 23, 25, 24}; */
  
colour_rapForPTBins = [1, 30, 4, 2, 3, ROOT.kMagenta+1];
marker_rapForPTBins = [20, 21, 25, 20, 22, 29];
  
#min number of entries to consider the bin
#(in acceptance map, etc)
minEntriesPerBin = 10

