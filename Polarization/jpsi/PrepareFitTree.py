#!/usr/bin/env python

import ROOT
import commonVar as jpsi
from optparse import OptionParser

from ROOT import TFile, TTree, TBranch, gROOT, TChain, AddressOf

def main(options,args):
	gROOT.Reset()

	#define trigger unprescaled periods
	jpsi_Trig_unP = {#'HLT_Mu0_Track0_Jpsi'            :[133446,141882],
			 #'HLT_Mu0_TkMu0_Jpsi'             :[140116,144114],
			 #'HLT_Mu0_TkMu0_OST_Jpsi'         :[146428,148058],
			 #'HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2':[148819,149182],
			 #'HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3':[149291,149442]}
			 'HLT_DoubleMu0'                  :[133446,147116],
			 'HLT_DoubleMu0_Quarkonium_v1'    :[147196,149442]}
	
   	#load calcPol routine
	gROOT.ProcessLine('.L calcPol.C+')

	inputs = TChain(options.treeName)
	for arg in args:
		print 'Adding: ',arg
		inputs.Add(arg)
        
	out = TFile.Open(options.output,'RECREATE')
	outTree = inputs.CloneTree(0)    
	# polarization stuff
	gROOT.ProcessLine('struct theBranches { TLorentzVector* JpsiP; TLorentzVector* muPosP; TLorentzVector* muNegP; Double_t JpsiMass; Double_t JpsiPt; Double_t JpsiRap; Double_t muPosPt; Double_t muPosEta; Double_t muPosPhi; Double_t muNegPt; Double_t muNegEta; Double_t muNegPhi; Double_t costh_CS; Double_t phi_CS; Double_t costh_HX; Double_t phi_HX; Double_t costh_PHX; Double_t phi_PHX; Double_t costh_sGJ; Double_t phi_sGJ; Double_t costh_GJ1; Double_t phi_GJ1; Double_t costh_GJ2; Double_t phi_GJ2; }')
	# trigger stuff
	gROOT.ProcessLine('struct triggers { Int_t runNb; Int_t HLT_Mu0_Track0_Jpsi; Int_t HLT_Mu0_TkMu0_Jpsi; Int_t HLT_Mu0_TkMu0_OST_Jpsi; Int_t HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2; Int_t HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3; Int_t HLT_DoubleMu0; Int_t HLT_DoubleMu0_Quarkonium_v1; }')

	branches = ROOT.theBranches()

	thetriggers = ROOT.triggers()

	inputs.SetBranchAddress('runNb',AddressOf(thetriggers,'runNb'))

	inputs.SetBranchAddress('HLT_Mu0_Track0_Jpsi',AddressOf(thetriggers,'HLT_Mu0_Track0_Jpsi'))
	inputs.SetBranchAddress('HLT_Mu0_TkMu0_Jpsi',AddressOf(thetriggers,'HLT_Mu0_TkMu0_Jpsi'))
	inputs.SetBranchAddress('HLT_Mu0_TkMu0_OST_Jpsi',AddressOf(thetriggers,'HLT_Mu0_TkMu0_OST_Jpsi'))
	inputs.SetBranchAddress('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2',AddressOf(thetriggers,'HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2'))
	inputs.SetBranchAddress('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3',AddressOf(thetriggers,'HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3'))
	inputs.SetBranchAddress('HLT_DoubleMu0',AddressOf(thetriggers,'HLT_DoubleMu0'))
	inputs.SetBranchAddress('HLT_DoubleMu0_Quarkonium_v1',AddressOf(thetriggers,'HLT_DoubleMu0_Quarkonium_v1'))

	outTree.SetBranchAddress('runNb',AddressOf(thetriggers,'runNb'))

	outTree.SetBranchAddress('HLT_Mu0_Track0_Jpsi',AddressOf(thetriggers,'HLT_Mu0_Track0_Jpsi'))
	outTree.SetBranchAddress('HLT_Mu0_TkMu0_Jpsi',AddressOf(thetriggers,'HLT_Mu0_TkMu0_Jpsi'))
	outTree.SetBranchAddress('HLT_Mu0_TkMu0_OST_Jpsi',AddressOf(thetriggers,'HLT_Mu0_TkMu0_OST_Jpsi'))
	outTree.SetBranchAddress('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2',AddressOf(thetriggers,'HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2'))
	outTree.SetBranchAddress('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3',AddressOf(thetriggers,'HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3'))
	outTree.SetBranchAddress('HLT_DoubleMu0',AddressOf(thetriggers,'HLT_DoubleMu0'))
	outTree.SetBranchAddress('HLT_DoubleMu0_Quarkonium_v1',AddressOf(thetriggers,'HLT_DoubleMu0_Quarkonium_v1'))

	inputs.SetBranchAddress('JpsiP',AddressOf(branches,'JpsiP'))
	inputs.SetBranchAddress('muPosP',AddressOf(branches,'muPosP'))
	inputs.SetBranchAddress('muNegP',AddressOf(branches,'muNegP'))

	phiFolded = [0,0,0,0,0,0]
	thetaAdjusted = [0,0,0,0,0,0]

	outTree.Branch("JpsiMass",AddressOf(branches,'JpsiMass'),"JpsiMass/D");
	outTree.Branch("JpsiPt",AddressOf(branches,'JpsiPt'),"JpsiPt/D");
	outTree.Branch("JpsiRap",AddressOf(branches,'JpsiRap'),"JpsiRap/D");
    
	outTree.Branch("muPosPt",AddressOf(branches,'muPosPt'),"muPosPt/D");
	outTree.Branch("muPosEta",AddressOf(branches,'muPosEta'),"muPosEta/D");
	outTree.Branch("muPosPhi",AddressOf(branches,'muPosPhi'),"muPosPhi/D");
	outTree.Branch("muNegPt",AddressOf(branches,'muNegPt'),"muNegPt/D");
	outTree.Branch("muNegEta",AddressOf(branches,'muNegEta'),"muNegEta/D");
	outTree.Branch("muNegPhi",AddressOf(branches,'muNegPhi'),"muNegPhi/D");
	
	outTree.Branch("costh_CS",AddressOf(branches,'costh_CS'),"costh_CS/D");
	outTree.Branch("phi_CS",AddressOf(branches,'phi_CS'),"phi_CS/D");
	
	outTree.Branch("costh_HX",AddressOf(branches,'costh_HX'),"costh_HX/D");
	outTree.Branch("phi_HX",AddressOf(branches,'phi_HX'),"phi_HX/D");
	
	outTree.Branch("costh_PHX",AddressOf(branches,'costh_PHX'),"costh_PHX/D");
	outTree.Branch("phi_PHX",AddressOf(branches,'phi_PHX'),"phi_PHX/D");
	
	outTree.Branch("costh_sGJ",AddressOf(branches,'costh_sGJ'),"costh_sGJ/D");
	outTree.Branch("phi_sGJ",AddressOf(branches,'phi_sGJ'),"phi_sGJ/D");
    
	outTree.Branch("costh_GJ1",AddressOf(branches,'costh_GJ1'),"costh_GJ1/D");
	outTree.Branch("phi_GJ1",AddressOf(branches,'phi_GJ1'),"phi_GJ1/D");
	
	outTree.Branch("costh_GJ2",AddressOf(branches,'costh_GJ2'),"costh_GJ2/D");
	outTree.Branch("phi_GJ2",AddressOf(branches,'phi_GJ2'),"phi_GJ2/D");
    
	nEvents=inputs.GetEntries()
	print 'There are ',nEvents,' events to process!'

	FiducialCutCountPos = 0
	FiducialCutCountNeg = 0

	for j in range(0,nEvents):
		nb = inputs.GetEntry(j)        
		
		if branches.muPosP.Pt() == branches.muNegP.Pt():
			continue
		
		passesTrigger = True
		
		if not options.isMC:
			passesTrigger=False
			for trig in jpsi_Trig_unP.keys():
				if thetriggers.runNb >= jpsi_Trig_unP[trig][0] and thetriggers.runNb <= jpsi_Trig_unP[trig][-1] and getattr(thetriggers,trig) == 1:
					passesTrigger = True		
		
		if not passesTrigger:
			continue
		
		branches.muPosPt = branches.muPosP.Pt()
		branches.muPosEta = branches.muPosP.Eta()
		branches.muPosPhi = branches.muPosP.Phi()
		
		branches.muNegPt = branches.muNegP.Pt()
		branches.muNegEta = branches.muNegP.Eta()
		branches.muNegPhi = branches.muNegP.Phi()
		
		#(a) on the positive muon
		if((abs(branches.muPosEta) < jpsi.etaPS[0] and branches.muPosPt < jpsi.pTMuMin[0]) or 
		   (abs(branches.muPosEta) > jpsi.etaPS[0] and abs(branches.muPosEta) < jpsi.etaPS[1] and branches.muPosP.P() < jpsi.pMuMin[1]) or
		   (abs(branches.muPosEta) > jpsi.etaPS[1] and abs(branches.muPosEta) < jpsi.etaPS[2] and branches.muPosPt < jpsi.pTMuMin[2])):
			FiducialCutCountPos += 1
			continue
        	#(b) on the negative muon
        	if((abs(branches.muNegEta) < jpsi.etaPS[0] and branches.muNegPt < jpsi.pTMuMin[0]) or
		   (abs(branches.muNegEta) > jpsi.etaPS[0] and abs(branches.muNegEta) < jpsi.etaPS[1] and branches.muNegP.P() < jpsi.pMuMin[1]) or
		   (abs(branches.muNegEta) > jpsi.etaPS[1] and abs(branches.muNegEta) < jpsi.etaPS[2] and branches.muNegPt < jpsi.pTMuMin[2])):
			FiducialCutCountNeg += 1
			continue
		
		ROOT.calcPol(branches.muPosP,branches.muNegP)

		#set up KLUDGE to make DoubleMu0 and DoubleMu0_Quarkonium to be the same as they have the same acceptance map
		if thetriggers.HLT_DoubleMu0_Quarkonium_v1 == 1:
			thetriggers.HLT_DoubleMu0 = 1

		if options.isFolded:
			for frame in range(len(jpsi.frameLabel)):
				phiFolded[frame] = ROOT.thisPhi[frame]
				thetaAdjusted[frame] = ROOT.thisCosTh[frame]
				if ROOT.thisPhi[frame] >= -90 and ROOT.thisPhi[frame] < 0.:
					phiFolded[frame] *= -1
				elif ROOT.thisPhi[frame] > 90 and ROOT.thisPhi[frame] < 180:
					phiFolded[frame] = 180 - ROOT.thisPhi[frame]
					thetaAdjusted[frame] *= -1
				elif ROOT.thisPhi[frame] > -180 and ROOT.thisPhi[frame] < -90:
					phiFolded[frame] = 180 + ROOT.thisPhi[frame]
					thetaAdjusted[frame] *= -1
			branches.costh_CS = thetaAdjusted[0]
			branches.phi_CS = phiFolded[0]
						
			branches.costh_HX = thetaAdjusted[1]
			branches.phi_HX = phiFolded[1]
						
			branches.costh_PHX = thetaAdjusted[2]
			branches.phi_PHX = phiFolded[2]
						
			branches.costh_sGJ = thetaAdjusted[3]
			branches.phi_sGJ = phiFolded[3]
						
			branches.costh_GJ1 = thetaAdjusted[4]
			branches.phi_GJ1 = phiFolded[4]
				
			branches.costh_GJ2 = thetaAdjusted[5]
			branches.phi_GJ2 = phiFolded[5]
		else:
			branches.costh_CS = ROOT.thisCosTh[0]
			branches.phi_CS = ROOT.thisPhi[0]
						
			branches.costh_HX = ROOT.thisCosTh[1]
			branches.phi_HX = ROOT.thisPhi[1]
						
			branches.costh_PHX = ROOT.thisCosTh[2]
			branches.phi_PHX = ROOT.thisPhi[2]
						
			branches.costh_sGJ = ROOT.thisCosTh[3]
			branches.phi_sGJ = ROOT.thisPhi[3]
						
			branches.costh_GJ1 = ROOT.thisCosTh[4]
			branches.phi_GJ1 = ROOT.thisPhi[4]
				
			branches.costh_GJ2 = ROOT.thisCosTh[5]
			branches.phi_GJ2 = ROOT.thisPhi[5]
			
		branches.JpsiMass = branches.JpsiP.M()
		branches.JpsiPt = branches.JpsiP.Pt()
		branches.JpsiRap = branches.JpsiP.Rapidity()       
		
		outTree.Fill()

	out.Write()
	out.Close()
        
if __name__ == '__main__':
	parser = OptionParser(description='%prog : J/Psi Polarization Fitter.',
	                      usage='polarizationFit.py --output=foo inputData.root ...')
	parser.add_option('--output',dest='output',help='name of output file')
	parser.add_option('--treeName',dest='treeName',default='data',help='Name of the input TTree.')
	parser.add_option('--mc',dest='isMC',default=False,action='store_true',help='Is data MC?')
	parser.add_option('--notFolded',dest='isFolded',default=True,action='store_false',help='Create non-folded polarization variables.')
	(options,args) = parser.parse_args()

	miss_options = False

	if options.output is None:
		print 'Need to specify --output'
		miss_options=True
        
	if miss_options:
		exit(1)

	main(options,args)
