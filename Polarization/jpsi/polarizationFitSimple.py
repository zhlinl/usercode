#!/usr/bin/env python

# implementation of the full polarization fit using a RooSimultaneous
# should be much more stable than what we had previously
# nuisance parameters will be fit without polarization first and
# then constrained during polarization fit.

import os
from optparse import OptionParser
import ROOT
import commonVar as jpsi

from math import sqrt

from ROOT import RooFit, RooWorkspace, RooArgSet, RooArgList, RooRealVar, gROOT, RooSimWSTool
from ROOT import RooCategory, RooFormulaVar, RooDataSet, RooAddPdf, RooDataHist, RooHistPdf
from ROOT import TFile, TTree, TChain, RooGExpModel, RooAddition, RooMinuit, RooHistFunc
from ROOT import RooSuperCategory, RooDataWeightedAverage, RooMinimizer
from ROOT import RooStats,RooProduct


def main(options,args):
    gROOT.Reset()

    #load our super special Polarization PDF and beta distribution wrapper
    gROOT.ProcessLine('.L RooPolarizationPdf.cxx+')
    gROOT.ProcessLine('.L RooPolarizationConstraint.cxx+')
    gROOT.ProcessLine('.L betaWrapper.cxx+')
    
    #setup integration
    intConf = ROOT.RooAbsReal.defaultIntegratorConfig()
    #intConf.Print('v')
    intConf.setPrintEvalCounter(True)
    #intConf.method1D().setLabel('RooAdaptiveGaussKronrodIntegrator1D')
    intConf.setEpsAbs(1e-13)
    intConf.setEpsRel(1e-13)
    print intConf.epsAbs()
    print intConf.epsRel()
    #intConf.method2D().setLabel('RooIntegrator2D')
    #intConf.method2D().setLabel()
    
    #intConf.methodND().setLabel('RooMCIntegrator')

    #intConf.getConfigSection('RooIntegrator1D').find('extrapolation').setLabel('None')
    #intConf.getConfigSection('RooIntegrator1D').find('maxSteps').setVal(500)
    #intConf.getConfigSection('RooIntegrator1D').find('minSteps').setVal(10)
    #intConf.getConfigSection('RooIntegrator1D').find('fixSteps').setVal(100)

    intConf.getConfigSection('RooIntegrator1D').Print('v')

    #exit(1)

    infile = TFile.Open(options.workspaceName+'.root')
    

    theWS = infile.Get(options.workspaceName)

    #save the polarization PDF code in the RooWorkspace
    
    theWS.importClassCode(ROOT.RooPolarizationPdf.Class(),True)
    theWS.importClassCode(ROOT.RooPolarizationConstraint.Class(),True)
    
    #buildDataAndCategories(theWS,options,args)    

    #buildMassAndLifetimePDF(theWS)

    theWS.Print("v")

    if options.fitFrame is not None:
        buildPolarizationPDF(theWS,options)
        
        theWS.Print('v')

        ROOT.RooMsgService.instance().Print()

        doPolFit(theWS,options)

    output = TFile.Open(options.workspaceName+'-'+options.fitFrame+'.root','RECREATE')
    output.cd()
    theWS.Write()
    output.Close()
    infile.Close()

def doPolFit(ws,options):
    rap_bins = range(1,len(jpsi.pTRange))
    pt_bins = None
    
    if options.testBin is not None:
        rap_bins = [int(options.testBin.split(',')[0])]
        pt_bins  = [int(options.testBin.split(',')[1])-1]
        
    for rap_bin in rap_bins:
        if options.testBin is None:
            pt_bins = range(len(jpsi.pTRange[rap_bin]))
        for pt_bin in pt_bins:

            sigMaxMass = jpsi.polMassJpsi[rap_bin] + jpsi.nSigMass*jpsi.sigmaMassJpsi[rap_bin]
            sigMinMass = jpsi.polMassJpsi[rap_bin] - jpsi.nSigMass*jpsi.sigmaMassJpsi[rap_bin]

            sbHighMass = jpsi.polMassJpsi[rap_bin] + jpsi.nSigBkgHigh*jpsi.sigmaMassJpsi[rap_bin]
            sbLowMass  = jpsi.polMassJpsi[rap_bin] - jpsi.nSigBkgLow*jpsi.sigmaMassJpsi[rap_bin]

            jPsiMass = ws.var('JpsiMass')
            jPsicTau = ws.var('Jpsict')

            jPsiMass.setRange('mlfit_promptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('mlfit_nonPromptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('mlfit_leftMassSideBand',2.7,sbLowMass)
            jPsiMass.setRange('mlfit_rightMassSideBandt',sbHighMass,3.5)

            jPsiMass.setRange('NormalizationRangeFormlfit_promptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('NormalizationRangeFormlfit_nonPromptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('NormalizationRangeFormlfit_leftMassSideBand',2.7,sbLowMass)
            jPsiMass.setRange('NormalizationRangeFormlfit_rightMassSideBandt',sbHighMass,3.5)   

            jPsicTau.setRange('NormalizationRangeFormlfit_promptSignal',-1,.1)
            jPsicTau.setRange('NormalizationRangeFormlfit_nonPromptSignal',.1,2.5)
            jPsicTau.setRange('NormalizationRangeFormlfit_leftMassSideBand',-1,2.5)
            jPsicTau.setRange('NormalizationRangeFormlfit_rightMassSideBand',-1,2.5)
            
            jPsicTau.setRange('mlfit_promptSignal',-1,.1)
            jPsicTau.setRange('mlfit_nonPromptSignal',.1,2.5)
            jPsicTau.setRange('mlfit_leftMassSideBand',-1,2.5)
            jPsicTau.setRange('mlfit_rightMassSideBand',-1,2.5)

            #load data
            data = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

            ws.var('costh_'+options.fitFrame).setBins(50)
            ws.var('phi_'+options.fitFrame).setBins(90)

            datahist = RooDataHist('datahist'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   'datahist',
                                   RooArgSet(ws.var('costh_'+options.fitFrame),
                                             ws.var('phi_'+options.fitFrame),
                                             ws.cat('mlRegion')),
                                   data
                                   )
            
            
            #reset parameters            
            ws.loadSnapshot('snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

            
            constrs = RooArgSet()

            constrs.add(ws.pdf('BetaFunPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            constrs.add(ws.pdf('BetaFunNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            constrs.add(ws.pdf('BetaFunBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))

            constrs.add(ws.pdf('BetaFunPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            constrs.add(ws.pdf('BetaFunNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            constrs.add(ws.pdf('BetaFunBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
                                
            if options.polOnly:
                NLLs.removeAll()
                
            PPdf = ws.pdf('PPdf'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))

            PPdf.Print("v")
            ws.pdf('promptPolExtended'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)).Print("v")
            ws.pdf('nonpromptPolExtended'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)).Print("v")
            ws.pdf('promptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)).Print("v")
            ws.pdf('nonpromptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)).Print("v")
                                
            #ws.var('lambda_phi_'+options.fitFrame+'_p').setVal(0)
            ws.var('lambda_theta_'+options.fitFrame+'_p').setVal(0)
            ws.var('lambda_thetaphi_'+options.fitFrame+'_p').setVal(0)
            
            #ws.var('lambda_phi_'+options.fitFrame+'_np').setVal(0)
            ws.var('lambda_theta_'+options.fitFrame+'_np').setVal(0)
            ws.var('lambda_thetaphi_'+options.fitFrame+'_np').setVal(0)                

            #ws.var('nNonPromptPol').setVal(1000)
            #ws.var('nPromptPol').setVal(1000)
                
            PolNLL = PPdf.createNLL(data,
                                    ROOT.RooFit.Range('mlfit'),
                                    ROOT.RooFit.SplitRange(False),                                    
                                    ROOT.RooFit.ConditionalObservables(RooArgSet(ws.var('JpsictErr'))),
                                    #ROOT.RooFit.ExternalConstraints(constrs),
                                    ROOT.RooFit.NumCPU(3),
                                    ROOT.RooFit.Verbose(False))
            
            simNLL = RooAddition('add','add',RooArgSet(PolNLL))

            minuit = RooMinuit(PolNLL)
            minuit.setStrategy(2)
            minuit.setPrintEvalErrors(-1)
            minuit.setErrorLevel(.5)
            #minuit.setEvalErrorWall(True)
            #minuit.setVerbose(True)

            minuit2 = RooMinimizer(PolNLL)
            minuit2.setStrategy(2)
            minuit2.setPrintEvalErrors(-1)
            minuit2.setErrorLevel(.5)
            #minuit2.setEvalErrorWall(True)
            #minuit2.setVerbose(True)

            bfgs = RooMinimizer(PolNLL)
            bfgs.setStrategy(2)
            bfgs.setPrintEvalErrors(-1)
            bfgs.setErrorLevel(.5)
            #bfgs.setEvalErrorWall(True)
            #bfgs.setVerbose(True)

            minuit.hesse()
            minuit.migrad()
            minuit.migrad()
            #minuit.improve()
            #minuit2.minimize('Minuit2','Scan')
            #minuit2.minimize('Minuit2','Combined')
            #bfgs.minimize('GSLMultiMin','ConjugatePR')
            #bfgs.minimize('GSLMultiMin','BFGS2')            
            minuit.hesse()
            #minuit.hesse()
        
            #minuit.migrad()
            #minuit.migrad()
            #minuit.hesse()

            #minuit.simplex()
            #minuit.minimize('Minuit2','migrad')
            #minuit.minimize('Minuit2','migrad')
            #minuit.minimize('GSLMultiMin','BFGS2')
            #minuit.minimize('Minuit2','hesse')

            minosVars = RooArgSet(ws.var('lambda_theta_'+options.fitFrame+'_p'),
                                  ws.var('lambda_thetaphi_'+options.fitFrame+'_p'))

            if options.lambdaPhiSub == 'lambda_phi':
                minosVars.add(ws.var('lambda_phi_'+options.fitFrame+'_p'))
            if options.lambdaPhiSub == 'F_invar':
                minosVars.add(ws.var('F_invar_p'))
            if options.lambdaPhiSub == 'lambda_tilde':
                minosVars.add(ws.var('lambda_tilde_p'))
            
            if options.doNonPrompt:                
                minosVars.add(ws.var('lambda_theta_'+options.fitFrame+'_np'))
                minosVars.add(ws.var('lambda_thetaphi_'+options.fitFrame+'_np'))

                if options.lambdaPhiSub == 'lambda_phi':
                    minosVars.add(ws.var('lambda_phi_'+options.fitFrame+'_np'))
                if options.lambdaPhiSub == 'F_invar':
                    minosVars.add(ws.var('F_invar_np'))
                if options.lambdaPhiSub == 'lambda_tilde':
                    minosVars.add(ws.var('lambda_tilde_np'))

                

            minuit.setProfile(True)
            minuit.minos( minosVars )
            
            fitresult = minuit.save('pol_fitresult_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
            getattr(ws,'import')(fitresult)           

            ws.saveSnapshot('pol_snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),ws.allVars())

def doMLFit(ws,options):
    rap_bins = range(1,len(jpsi.pTRange))
    pt_bins = None
    
    if options.testBin is not None:
        rap_bins = [int(options.testBin.split(',')[0])]
        pt_bins  = [int(options.testBin.split(',')[1])-1]
        
    for rap_bin in rap_bins:
        if options.testBin is None:
            pt_bins = range(len(jpsi.pTRange[rap_bin]))
        for pt_bin in pt_bins:

            sigMaxMass = jpsi.polMassJpsi[rap_bin] + jpsi.nSigMass*jpsi.sigmaMassJpsi[rap_bin]
            sigMinMass = jpsi.polMassJpsi[rap_bin] - jpsi.nSigMass*jpsi.sigmaMassJpsi[rap_bin]

            sbHighMass = jpsi.polMassJpsi[rap_bin] + jpsi.nSigBkgHigh*jpsi.sigmaMassJpsi[rap_bin]
            sbLowMass  = jpsi.polMassJpsi[rap_bin] - jpsi.nSigBkgLow*jpsi.sigmaMassJpsi[rap_bin]

            jPsiMass = ws.var('JpsiMass')
            jPsicTau = ws.var('Jpsict')
            
            jPsiMass.setRange('mlfit_promptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('mlfit_nonPromptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('mlfit_leftMassSideBand',2.7,sbLowMass)
            jPsiMass.setRange('mlfit_rightMassSideBandt',sbHighMass,3.5)
            
            jPsiMass.setRange('NormalizationRangeFormlfit_promptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('NormalizationRangeFormlfit_nonPromptSignal',sigMinMass,sigMaxMass)
            jPsiMass.setRange('NormalizationRangeFormlfit_leftMassSideBand',2.7,sbLowMass)
            jPsiMass.setRange('NormalizationRangeFormlfit_rightMassSideBandt',sbHighMass,3.5)   

            jPsicTau.setRange('NormalizationRangeFormlfit_promptSignal',-1,.1)
            jPsicTau.setRange('NormalizationRangeFormlfit_nonPromptSignal',.1,2.5)
            jPsicTau.setRange('NormalizationRangeFormlfit_leftMassSideBand',-1,2.5)
            jPsicTau.setRange('NormalizationRangeFormlfit_rightMassSideBand',-1,2.5)
            
            jPsicTau.setRange('mlfit_promptSignal',-1,.1)
            jPsicTau.setRange('mlfit_nonPromptSignal',.1,2.5)
            jPsicTau.setRange('mlfit_leftMassSideBand',-1,2.5)
            jPsicTau.setRange('mlfit_rightMassSideBand',-1,2.5)

            #reset parameters
            ws.var('CBn').setVal(10)
            ws.var('CBalpha').setVal(2.1)
            ws.var('CBmass').setVal(3.096)
            ws.var('CBsigma').setVal(.02)
            ws.var('bkgLambda').setVal(-.5)
            
            ws.var('bkgTauSSDL').setVal(.5)
            #ws.var('bkgTauFDL').setVal(.5)
            ws.var('bkgTauDSDL').setVal(.5)
            ws.var('fBkgSSDL').setVal(.5)
            ws.var('fBkgLR').setVal(.5)

            ws.var('bkgTauSSDR').setVal(.5)
            #ws.var('bkgTauFDR').setVal(.5)
            ws.var('bkgTauDSDR').setVal(.5)
            ws.var('fBkgSSDR').setVal(.5)
            #ws.var('fBkgFDR').setVal(.25)
            
            ws.var('nPromptSignal').setVal(5000)
            ws.var('nPromptL').setVal(100)
            ws.var('nPromptR').setVal(100)
            ws.var('nNonPromptSignal').setVal(500)
            ws.var('nNonPromptL').setVal(50)
            ws.var('nNonPromptR').setVal(50)
            ws.var('nBackgroundSignal').setVal(100)
            ws.var('nBackgroundL').setVal(50)
            ws.var('nBackgroundR').setVal(50)
            
            ws.var('nonPromptTau').setVal(.32)
            ws.var('promptMean').setVal(0)
            ws.var('ctResolution').setVal(1)

                        
            LPdf = ws.pdf('LPdf')
            MPdf = ws.pdf('MPdf')            
            
            MLPdf = ws.pdf('MLPdf')
                
            data = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

            NLLs = RooArgSet()

            MassNLL = MPdf.createNLL(data,
                                     ROOT.RooFit.Range('mlfit'),
                                     ROOT.RooFit.SplitRange(False),                                    
                                     ROOT.RooFit.ConditionalObservables(RooArgSet(ws.var('JpsictErr'))),
                                     ROOT.RooFit.NumCPU(2))
                                     

            CTauNLL = LPdf.createNLL(data,
                                     ROOT.RooFit.Range('mlfit'),
                                     ROOT.RooFit.SplitRange(False),                                    
                                     ROOT.RooFit.ConditionalObservables(RooArgSet(ws.var('JpsictErr'))),
                                     ROOT.RooFit.NumCPU(2))
                
            NLLs.add(MassNLL)
            NLLs.add(CTauNLL)

            massMinuit = RooMinuit(MassNLL)
            massMinuit.setStrategy(2)
            massMinuit.setPrintEvalErrors(-1)
            
            massMinuit.migrad() #minimize the mass first to get a first good guess on many parameters
            
            simNLL = RooAddition('add','add',NLLs)

            minuit = RooMinuit(simNLL)
            minuit.setStrategy(2)
            minuit.setPrintEvalErrors(-1)
            minuit.setEvalErrorWall(False)
            minuit.setVerbose(False)
            
            minuit.simplex()
            minuit.migrad()
            minuit.migrad()
            minuit.hesse()
            
            fitresult = minuit.save('ml_fitresult_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
            getattr(ws,'import')(fitresult)           

            ws.saveSnapshot('ml_snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),ws.allVars())
            

def buildDataAndCategories(ws,options,args):
    #Get the input data
    inputData = TChain(options.treeName,'The input data')
    for arg in args:
        print 'Adding data from: ',arg
        inputData.Add(arg)

    foldname = ''
    phirange = [0,90]
    
    if not options.folded:
        foldname=''
        phirange = [-180,180]
    
    #variables necessary for j/psi mass,lifetime,polarization fit
    jPsiMass      = RooRealVar('JpsiMass','M [GeV]',2.7,3.5)
    jPsiRap       = RooRealVar('JpsiRap','#nu',-2.3,2.3)
    jPsiPt        = RooRealVar("JpsiPt","pT [GeV]",0,40);
    jPsicTau      = RooRealVar('Jpsict','l_{J/#psi} [mm]',-1,2.5)
    jPsicTauError = RooRealVar('JpsictErr','Error on l_{J/#psi} [mm]',0,2)
    jPsiHXcosth   = None
    jPsiHXphi     = None

    jPsicTau.setBins(10000,"cache")
        
    if options.fitFrame is not None:
        jPsiHXcosth   = RooRealVar('costh_'+options.fitFrame+foldname,'cos(#theta)_{'+options.fitFrame+'}',-1,1)
        jPsiHXphi     = RooRealVar('phi_'+options.fitFrame+foldname,'#phi_{'+options.fitFrame+'}',phirange[0],phirange[1])
    else:
        jPsiHXcosth   = RooRealVar('costh_CS'+foldname,'cos(#theta)_{CS}',-1,1)
        jPsiHXphi     = RooRealVar('phi_CS'+foldname,'#phi_{CS}',phirange[0],phirange[1])
    
    #vars needed for on the fly calc of polarization variables
    jPsimuPosPx = RooRealVar('muPosPx','+ Muon P_{x} [GeV]',0)
    jPsimuPosPy = RooRealVar('muPosPy','+ Muon P_{y} [GeV]',0)
    jPsimuPosPz = RooRealVar('muPosPz','+ Muon P_{z} [GeV]',0)
    jPsimuNegPx = RooRealVar('muNegPx','- Muon P_{x} [GeV]',0)
    jPsimuNegPy = RooRealVar('muNegPy','- Muon P_{y} [GeV]',0)
    jPsimuNegPz = RooRealVar('muNegPz','- Muon P_{z} [GeV]',0)

    #create RooArgSet for eventual dataset creation
    dataVars = RooArgSet(jPsiMass,jPsiRap,jPsiPt,
                         jPsicTau,jPsicTauError,
                         jPsimuPosPx,jPsimuPosPy,jPsimuPosPz)
    
    #add trigger requirement if specified
    if options.triggerName:
        trigger = RooRealVar(options.triggerName,'Passes Trigger',0.5,1.5)
        dataVars.add(trigger)
    
    dataVars.add(jPsimuNegPx)
    dataVars.add(jPsimuNegPy)
    dataVars.add(jPsimuNegPz)
    dataVars.add(jPsiHXcosth)
    dataVars.add(jPsiHXphi)
    
    
    redVars = RooArgSet(jPsiMass,jPsiRap,jPsiPt,
                        jPsicTau,jPsicTauError)
    redVars.add(jPsiHXcosth)
    redVars.add(jPsiHXphi)
    fitVars = redVars.Clone()    

    ### HERE IS WHERE THE BIT FOR CALCULATING POLARIZATION VARS GOES

    ctauStates = RooCategory('ctauRegion','Cut Region in lifetime')
    ctauStates.defineType('prompt',0)
    ctauStates.defineType('nonPrompt',1)

    massStates = RooCategory('massRegion','Cut Region in mass')
    massStates.defineType('signal',1)
    massStates.defineType('separation',0)
    massStates.defineType('leftMassSideBand',-2)
    massStates.defineType('rightMassSideBand',-1)

    states = RooCategory('mlRegion','Cut Region in mass')
    states.defineType('nonPromptSignal',2)
    states.defineType('promptSignal',1)
    states.defineType('separation',0)
    states.defineType('leftMassSideBand',-2)
    states.defineType('rightMassSideBand',-1)

    #define corresponding ranges in roorealvars
    #mass is a little tricky since the sidebands change definitions in each rap bin
    #define the names here and change as we do the fits
    #jPsiMass.setRange('NormalizationRangeFormlfit_promptSignal',2.7,3.5)
    #jPsiMass.setRange('NormalizationRangeFormlfit_nonPromptSignal',2.7,3.5)
    #jPsiMass.setRange('NormalizationRangeFormlfit_leftMassSideBand',2.7,3.1)
    #jPsiMass.setRange('NormalizationRangeFormlfit_rightMassSideBand',3.1,3.5)

    #want the prompt fit only done in prompt region
    #non-prompt only in non-prompt region
    #background over entire cTau range
    #jPsicTau.setRange('NormalizationRangeFormlfit_promptSignal',-1,.1)
    #jPsicTau.setRange('NormalizationRangeFormlfit_nonPromptSignal',.1,2.5)
    #jPsicTau.setRange('NormalizationRangeFormlfit_leftMassSideBand',-1,2.5)
    #jPsicTau.setRange('NormalizationRangeFormlfit_rightMassSideBand',-1,2.5)

    #redVars.add(ctauStates)
    #redVars.add(massStates)
    #redVars.add(states)
    fitVars.add(ctauStates)
    fitVars.add(massStates)
    fitVars.add(states)
    
    fullData = RooDataSet('fullData','The Full Data From the Input ROOT Trees',
                          dataVars,
                          ROOT.RooFit.Import(inputData))    

    for rap_bin in range(1,len(jpsi.pTRange)):
        yMin  = jpsi.rapForPTRange[rap_bin-1][0]
        yMax  = jpsi.rapForPTRange[rap_bin-1][-1]
        for pt_bin in range(len(jpsi.pTRange[rap_bin])):

            ptMin = jpsi.pTRange[rap_bin][pt_bin][0]
            ptMax = jpsi.pTRange[rap_bin][pt_bin][-1]               

            sigMaxMass = jpsi.polMassJpsi[rap_bin] + jpsi.nSigMass*jpsi.sigmaMassJpsi[rap_bin]
            sigMinMass = jpsi.polMassJpsi[rap_bin] - jpsi.nSigMass*jpsi.sigmaMassJpsi[rap_bin]

            sbHighMass = jpsi.polMassJpsi[rap_bin] + jpsi.nSigBkgHigh*jpsi.sigmaMassJpsi[rap_bin]
            sbLowMass  = jpsi.polMassJpsi[rap_bin] - jpsi.nSigBkgLow*jpsi.sigmaMassJpsi[rap_bin]

            ctauNonPrompt = .1
            
            massFun = RooFormulaVar('massRegion','Function that returns the mass state.',
                                     '('+jPsiMass.GetName()+' < '+str(sigMaxMass)+' && '+jPsiMass.GetName()+' > '+str(sigMinMass)+
                                     ') - ('+jPsiMass.GetName()+' > '+str(sbHighMass)+')'+
                                     '-2*('+jPsiMass.GetName()+' < '+str(sbLowMass)+')',
                                     RooArgList(jPsiMass,jPsicTau))
            
            ctauFun = RooFormulaVar('ctauRegion','Function that returns the ctau state.',
                                     '('+jPsicTau.GetName()+' > '+str(ctauNonPrompt)+')',
                                     RooArgList(jPsiMass,jPsicTau))

            mlFun = RooFormulaVar('mlRegion','Function that returns the mass and lifetime state.',
                                  '('+jPsiMass.GetName()+' < '+str(sigMaxMass)+' && '+jPsiMass.GetName()+' > '+str(sigMinMass)+
                                  ') + ('+jPsiMass.GetName()+' < '+str(sigMaxMass)+' && '+jPsiMass.GetName()+' > '+
                                  str(sigMinMass)+' && '+jPsicTau.GetName()+' > '+str(ctauNonPrompt)+
                                  ') - ('+jPsiMass.GetName()+' > '+str(sbHighMass)+')'+
                                  '-2*('+jPsiMass.GetName()+' < '+str(sbLowMass)+')',
                                  RooArgList(jPsiMass,jPsicTau))
            

            cutStringPt = '('+jPsiPt.GetName()+' > '+str(ptMin)+' && '+jPsiPt.GetName()+' < '+str(ptMax)+')'
            cutStringY  = '( abs('+jPsiRap.GetName()+') > '+str(yMin)+' && abs('+jPsiRap.GetName()+') < '+str(yMax)+')'
            #cutStringM1 = '('+jPsiMass.GetName()+' < '+str(sigMinMass)+' && '+jPsiMass.GetName()+' > '+str(sbLowMass)+')'
            #cutStringM2 = '('+jPsiMass.GetName()+' < '+str(sbHighMass)+' && '+jPsiMass.GetName()+' > '+str(sigMaxMass)+')'
            #cutStringMT = '!('+cutStringM1+' || '+cutStringM2+')'
            cutString   = cutStringPt+' && '+cutStringY #+' && '+cutStringMT

            print cutString

            #get the reduced dataset we'll do the fit on
            binData = fullData.reduce(ROOT.RooFit.SelectVars(redVars),
                                      ROOT.RooFit.Cut(cutString),
                                      ROOT.RooFit.Name('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                      ROOT.RooFit.Title('Data For Fitting'))

            binDataWithCategory = RooDataSet('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                             'Data For Fitting',
                                             fitVars)
            #categorize
            binData.addColumn(ctauStates)
            binData.addColumn(massStates)
            binData.addColumn(states)
            for ev in range(binData.numEntries()):
                args = binData.get(ev)
                                
                jPsiMass.setVal(args.find(jPsiMass.GetName()).getVal())
                jPsiRap.setVal(args.find(jPsiRap.GetName()).getVal())
                jPsiPt.setVal(args.find(jPsiPt.GetName()).getVal())
                jPsicTau.setVal(args.find(jPsicTau.GetName()).getVal())                    
                jPsicTauError.setVal(args.find(jPsicTauError.GetName()).getVal())
            
                jPsiHXcosth.setVal(args.find(jPsiHXcosth.GetName()).getVal())
                jPsiHXphi.setVal(args.find(jPsiHXphi.GetName()).getVal())

                massStates.setIndex(int(massFun.getVal()))
                ctauStates.setIndex(int(ctauFun.getVal()))
                states.setIndex(int(mlFun.getVal()))
                
                binDataWithCategory.add(fitVars)
            

            getattr(ws,'import')(binDataWithCategory)

#    states.Print("v")
#    getattr(ws,'import')(states)
    
def buildMassAndLifetimePDF(ws):
    
    #define mass shape
    ws.factory('RooCBShape::massCBShape(JpsiMass,CBmass[3.1,3.05,3.15],CBsigma[0.02,0.0001,1],CBalpha[1,.0001,5],CBn[10,.0001,50])')
    #ws.factory('RooGaussian::massCBShape(JpsiMass,CBmass[3.1,3.05,3.15],CBsigma[0.02,0.0001,1])')
    ws.factory('RooExponential::bkgMassShape(JpsiMass,bkgLambda[0,-5,5])')
    #ws.factory('RooChebychev::bkgMassShape(JpsiMass,{p0[1,-10,10],p1[1,-10,10],p2[1,-10,10]})')    
    #lifetime
    ws.factory('RooGaussModel::promptLifetimeRaw(Jpsict,promptMean[0,-1,1],ctResolution[1,.001,5],1,JpsictErr)')    
    ws.pdf('promptLifetimeRaw').advertiseFlatScaleFactorIntegral(True)
    
    ws.factory('RooDecay::nonPromptSSDRaw(Jpsict,nonPromptTau[.3,.01,3],promptLifetimeRaw,RooDecay::SingleSided)')
    
    ws.factory('RooDecay::backgroundSSDRawL(Jpsict,bkgTauSSDL[.3,0,3],promptLifetimeRaw,RooDecay::SingleSided)')
    #ws.factory('RooDecay::backgroundFDRaw(Jpsict,bkgTauFD[.3,.0001,3],promptLifetimeRaw,RooDecay::Flipped)')
    ws.factory('RooDecay::backgroundDSDRawL(Jpsict,bkgTauDSDL[.3,0,3],promptLifetimeRaw,RooDecay::DoubleSided)')    
    ws.factory('SUM::backgroundRawL(fBkgSSDL[.5,0,1]*backgroundSSDRawL,backgroundDSDRawL)')
    
    ws.factory('RooDecay::backgroundSSDRawR(Jpsict,bkgTauSSDR[.3,0,3],promptLifetimeRaw,RooDecay::SingleSided)')
    #ws.factory('RooDecay::backgroundFDRaw(Jpsict,bkgTauFD[.3,.0001,3],promptLifetimeRaw,RooDecay::Flipped)')
    ws.factory('RooDecay::backgroundDSDRawR(Jpsict,bkgTauDSDR[.3,0,3],promptLifetimeRaw,RooDecay::DoubleSided)')    
    ws.factory('SUM::backgroundRawR(fBkgSSDR[.5,0,1]*backgroundSSDRawR,backgroundDSDRawR)')

    ws.factory('SUM::backgroundRaw(fBkgLR[.5]*backgroundRawL,backgroundDSDRawR)')

    ws.factory('PROD::promptRawMass(massCBShape)') #,promptLifetimeRaw|JpsictErr
    ws.factory('PROD::nonpromptRawMass(massCBShape)') #,nonPromptSSDRaw|JpsictErr
    ws.factory('PROD::backgroundRawMass(bkgMassShape)') #,backgroundRaw|JpsictErr
    
    ws.factory('PROD::promptRawMassLife(massCBShape,promptLifetimeRaw|JpsictErr)') #
    ws.factory('PROD::nonpromptRawMassLife(massCBShape,nonPromptSSDRaw|JpsictErr)') #
    ws.factory('PROD::backgroundRawMassLife(bkgMassShape,backgroundRaw|JpsictErr)') #
    ws.factory('PROD::backgroundRawMassLifeL(bkgMassShape,backgroundRawL|JpsictErr)') #
    ws.factory('PROD::backgroundRawMassLifeR(bkgMassShape,backgroundRawR|JpsictErr)') #

    ws.factory('PROD::promptRawCTau(promptLifetimeRaw|JpsictErr)') 
    ws.factory('PROD::nonpromptRawCTau(nonPromptSSDRaw|JpsictErr)') 
    ws.factory('PROD::backgroundRawCTauL(backgroundRawL|JpsictErr)') 
    ws.factory('PROD::backgroundRawCTauR(backgroundRawR|JpsictErr)') 
    ws.factory('PROD::backgroundRawCTau(backgroundRaw|JpsictErr)')
    
    #extended pdfs    
    ws.factory('RooExtendPdf::promptExtMass(promptRawMass,sum::nPrompt(nPromptSignal[1000,0,1000000],nPromptL[1000,0,1000000],nPromptR[1000,0,1000000]))')
    ws.factory('RooExtendPdf::nonpromptExtMass(nonpromptRawMass,sum::nNonPrompt(nNonPromptSignal[500,0,1000000],nNonPromptL[500,0,1000000],nNonPromptR[500,0,1000000]))')
    ws.factory('RooExtendPdf::backgroundExtMass(backgroundRawMass,sum::nBackground(nBackgroundSignal[100,0,1000000],nBackgroundL[100,0,1000000],nBackgroundR[100,0,1000000]))')

    ws.factory('RooExtendPdf::promptExtCTauL(promptRawCTau,nPromptL)')
    ws.factory('RooExtendPdf::promptExtCTauR(promptRawCTau,nPromptR)')
    ws.factory('RooExtendPdf::nonpromptExtCTauL(nonpromptRawCTau,nNonPromptL)')
    ws.factory('RooExtendPdf::nonpromptExtCTauR(nonpromptRawCTau,nNonPromptR)')
    ws.factory('RooExtendPdf::backgroundExtCTauL(backgroundRawCTauL,nBackgroundL)')
    ws.factory('RooExtendPdf::backgroundExtCTauR(backgroundRawCTauR,nBackgroundR)')

    ws.factory('RooExtendPdf::promptExtCTau(promptRawCTau,nPromptSignal)')
    ws.factory('RooExtendPdf::nonpromptExtCTau(nonpromptRawCTau,nNonPromptSignal)')
    ws.factory('RooExtendPdf::backgroundExtCTau(backgroundRawCTau,nBackgroundSignal)')

    #prod pdfs of mass and lifetime... should work....
    #ws.factory('RooExtendPdf::promptExtMassCTauL(promptRawMassLife,nPromptL)')    
    #ws.factory('RooExtendPdf::nonpromptExtMassCTauL(nonpromptRawMassLife,nNonPromptL)')    
    #ws.factory('RooExtendPdf::backgroundExtMassCTauL(backgroundRawCTauL,nBackgroundL)')

    #ws.factory('RooExtendPdf::promptExtMassCTauR(promptRawMassLife,nPromptR)')
    #ws.factory('RooExtendPdf::nonpromptExtMassCTauR(nonpromptRawMassLife,nNonPromptR)')
    #ws.factory('RooExtendPdf::backgroundExtMassCTauR(backgroundRawCTauR,nBackgroundR)')

    #ws.factory('RooExtendPdf::promptExtMassCTauinP(promptRawMassLife,nPromptSignalinP[1000,0,1000000])')
    #ws.factory('RooExtendPdf::nonpromptExtMassCTauinP(nonpromptRawMassLife,nNonPromptSignalinP[500,0,1000000])')
    #ws.factory('RooExtendPdf::backgroundExtMassCTauinP(backgroundRawMassLife,nBackgroundSignalinP[50,0,1000000])')

    #ws.factory('RooExtendPdf::promptExtMassCTauinNP(promptRawMassLife,nPromptSignalinNP[500,0,1000000])')
    #ws.factory('RooExtendPdf::nonpromptExtMassCTauinNP(nonpromptRawMassLife,nNonPromptSignalinNP[1000,0,1000000])')
    #ws.factory('RooExtendPdf::backgroundExtMassCTauinNP(backgroundRawMassLife,nBackgroundSignalinNP[50,0,1000000])')
    
    
    #final AddPdfs
    promptMArgList = RooArgList(ws.pdf('promptExtMass'))
    nonPromptMArgList = RooArgList(ws.pdf('promptExtMass'),
                                   ws.pdf('nonpromptExtMass'))
    
    if options.doNonPrompt:
        promptMArgList.add(ws.pdf('nonpromptExtMass'))
    if options.doBackground:        
        promptMArgList.add(ws.pdf('backgroundExtMass'))
        nonPromptMArgList.add(ws.pdf('backgroundExtMass'))
    
    promptM = RooAddPdf('promptMass','prompt',promptMArgList)
    nonpromptM = RooAddPdf('nonpromptMass','nonprompt',nonPromptMArgList)
    backgroundM = RooAddPdf('backgroundMass','background',RooArgList(ws.pdf('promptExtMass'),
                                                                     ws.pdf('nonpromptExtMass'),
                                                                     ws.pdf('backgroundExtMass')))
    
    promptCTArgList = RooArgList(ws.pdf('promptExtCTau'))
    nonPromptCTArgList = RooArgList(ws.pdf('promptExtCTau'),
                                    ws.pdf('nonpromptExtCTau'))
    
    if options.doNonPrompt:
        promptCTArgList.add(ws.pdf('nonpromptExtCTau'))        
    if options.doBackground:        
        promptCTArgList.add(ws.pdf('backgroundExtCTau'))
        nonPromptCTArgList.add(ws.pdf('backgroundExtCTau'))

    promptCT = RooAddPdf('promptCTau','prompt',promptCTArgList)
    nonpromptCT = RooAddPdf('nonpromptCTau','nonprompt',nonPromptCTArgList)
    backgroundCT = RooAddPdf('backgroundCTau','background',RooArgList(ws.pdf('promptExtCTau'),
                                                                      ws.pdf('nonpromptExtCTau'),
                                                                      ws.pdf('backgroundExtCTau')))
    #playing around with signal in left sideband
    backgroundCTL = RooAddPdf('backgroundCTauL','background',RooArgList(ws.pdf('promptExtCTauL'),
                                                                        ws.pdf('nonpromptExtCTauL'),
                                                                        ws.pdf('backgroundExtCTauL')))
    backgroundCTR = RooAddPdf('backgroundCTauR','background',RooArgList(ws.pdf('promptExtCTauR'),
                                                                        ws.pdf('nonpromptExtCTaR'),
                                                                        ws.pdf('backgroundExtCTauR')))   

    promptMLArgList = RooArgList(ws.pdf('promptExtMassCTauinP'),
                                 ws.pdf('nonpromptExtMassCTauinP'),
                                 ws.pdf('backgroundExtMassCTauinP'))
    nonPromptMLArgList = RooArgList(ws.pdf('promptExtMassCTauinNP'),
                                    ws.pdf('nonpromptExtMassCTauinNP'),
                                    ws.pdf('backgroundExtMassCTauinNP'))
    bkgMLArgListL = RooArgList(ws.pdf('promptExtMassCTauL'),
                               ws.pdf('nonpromptExtMassCTauL'),
                               ws.pdf('backgroundExtMassCTauL'))
    bkgMLArgListR = RooArgList(ws.pdf('promptExtMassCTauR'),
                               ws.pdf('nonpromptExtMassCTauR'),
                               ws.pdf('backgroundExtMassCTauR'))

    promptMCT = RooAddPdf('promptMCT','prompt',promptMLArgList)
    nonpromptMCT = RooAddPdf('nonpromptMCT','nonprompt',nonPromptMLArgList)
    backgroundMCTL = RooAddPdf('backgroundMCTL','backgroundL',bkgMLArgListL)
    backgroundMCTR = RooAddPdf('backgroundMCTR','backgroundR',bkgMLArgListR)


    getattr(ws,'import')(promptMCT)
    getattr(ws,'import')(nonpromptMCT)
    getattr(ws,'import')(backgroundMCTL)
    getattr(ws,'import')(backgroundMCTR)

    getattr(ws,'import')(promptM)
    getattr(ws,'import')(nonpromptM)
    getattr(ws,'import')(backgroundM)

    getattr(ws,'import')(promptCT)
    getattr(ws,'import')(nonpromptCT)
    getattr(ws,'import')(backgroundCT)
    getattr(ws,'import')(backgroundCTL)
    getattr(ws,'import')(backgroundCTR)

    massConfig = 'SIMUL::MPdf(ctauRegion,'
    ctConfig = 'SIMUL::LPdf(massRegion,'

    mlConfig = 'SIMUL::MLPdf(mlRegion,'
    
    if options.doPrompt:
        massConfig += 'prompt=promptMass,'
        ctConfig += 'signal=promptCTau,'
        mlConfig += 'promptSignal=promptMCT,'
    if options.doNonPrompt:
        massConfig += 'nonPrompt=nonpromptMass'
        mlConfig += 'nonPromptSignal=nonpromptMCT,'
    if options.doBackground:               
        ctConfig += 'leftMassSideBand=backgroundCTauL,rightMassSideBand=backgroundCTauR'
        mlConfig += 'leftMassSideBand=backgroundMCTL,rightMassSideBand=backgroundMCTR'

    massConfig += ')'
    ctConfig += ')'
    mlConfig += ')'
    
    #simultaneous
    ws.factory(massConfig)
    ws.factory(ctConfig)
    ws.factory(mlConfig)

def buildPolarizationPDF(ws,options):

    foldname = ''
    if not options.folded:
        foldname=''
        
    accMapsP = TFile.Open(options.acceptanceMap.split(',')[0])
    recoEffMapsP = TFile.Open(options.recoEfficiencyMap.split(',')[0])
    trigEffMapsP = TFile.Open(options.trigEfficiencyMap.split(',')[0])

    accMapsNP = TFile.Open(options.acceptanceMap.split(',')[-1])
    recoEffMapsNP = TFile.Open(options.recoEfficiencyMap.split(',')[-1])
    trigEffMapsNP = TFile.Open(options.trigEfficiencyMap.split(',')[-1])

    ws.factory('lambda_theta_'+options.fitFrame+'_p[0,-1,1]')
    
    ws.factory('lambda_thetaphi_'+options.fitFrame+'_p[0,-1,1]')
    ws.var('lambda_thetaphi_'+options.fitFrame+'_p').setMin(-sqrt(2)/2)
    ws.var('lambda_thetaphi_'+options.fitFrame+'_p').setMax(sqrt(2)/2)    

    ws.factory('lambda_theta_'+options.fitFrame+'_np[0,-1,1]')    
    ws.factory('lambda_thetaphi_'+options.fitFrame+'_np[0,-1,1]')
    ws.var('lambda_thetaphi_'+options.fitFrame+'_np').setMin(-sqrt(2)/2)
    ws.var('lambda_thetaphi_'+options.fitFrame+'_np').setMax(sqrt(2)/2)

    if options.lambdaPhiSub == 'lambda_phi':
        ws.factory('lambda_phi_'+options.fitFrame+'_p[0,-1,1]')
        ws.factory('lambda_phi_'+options.fitFrame+'_np[0,-1,1]')
    if options.lambdaPhiSub == 'F_invar':
        ws.factory('F_invar_p[.3333,0,1]')
        ws.factory('F_invar_np[.3333,0,1]')
        ws.factory('expr::lambda_phi_'+options.fitFrame+'_p(".5*(@0*(3+@1)-1-@1)",F_invar_p,lambda_theta_'+options.fitFrame+'_p)')
        ws.factory('expr::lambda_phi_'+options.fitFrame+'_np(".5*(@0*(3+@1)-1-@1)",F_invar_np,lambda_theta_'+options.fitFrame+'_np)')
    if options.lambdaPhiSub == 'lambda_tilde':
        ws.factory('lambda_tilde_p[0,-1,100]')
        ws.factory('lambda_tilde_np[0,-1,100]')
        ws.factory('expr::lambda_phi_'+options.fitFrame+'_p("(@0-@1)/(3+@0)",lambda_tilde_p,lambda_theta_'+options.fitFrame+'_p)')
        ws.factory('expr::lambda_phi_'+options.fitFrame+'_np("(@0-@1)/(3+@0)",lambda_tilde_np,lambda_theta_'+options.fitFrame+'_np)')
    
    ws.defineSet('POI',
                 RooArgSet(ws.var('lambda_theta_'+options.fitFrame+'_p'),                           
                           ws.var('lambda_thetaphi_'+options.fitFrame+'_p'),
                           ws.var('lambda_theta_'+options.fitFrame+'_np'),                           
                           ws.var('lambda_thetaphi_'+options.fitFrame+'_np')),
                 True)

    ws.factory('RooUniform::prior_lth_'+options.fitFrame+'_p(lambda_theta_'+options.fitFrame+'_p)')    
    ws.factory('RooUniform::prior_lthphi_'+options.fitFrame+'_p(lambda_thetaphi_'+options.fitFrame+'_p)')

    ws.factory('RooUniform::prior_lth_'+options.fitFrame+'_np(lambda_theta_'+options.fitFrame+'_np)')    
    ws.factory('RooUniform::prior_lthphi_'+options.fitFrame+'_np(lambda_thetaphi_'+options.fitFrame+'_np)')

    
    PriorPOIList= RooArgList(ws.pdf('prior_lth_'+options.fitFrame+'_p'),
                             ws.pdf('prior_lthphi_'+options.fitFrame+'_p'),
                             ws.pdf('prior_lth_'+options.fitFrame+'_np'),                             
                             ws.pdf('prior_lthphi_'+options.fitFrame+'_np'))

    if options.lambdaPhiSub == 'lambda_phi':
        ws.set('POI').add(ws.var('lambda_phi_'+options.fitFrame+'_p'))
        ws.set('POI').add(ws.var('lambda_phi_'+options.fitFrame+'_np'))

        ws.factory('RooUniform::prior_lphi_'+options.fitFrame+'_p(lambda_phi_'+options.fitFrame+'_p)')
        ws.factory('RooUniform::prior_lphi_'+options.fitFrame+'_np(lambda_phi_'+options.fitFrame+'_np)')
        
        PriorPOIList.add(ws.pdf('prior_lphi_'+options.fitFrame+'_p'))
        PriorPOIList.add(ws.pdf('prior_lphi_'+options.fitFrame+'_np'))        
        
    if options.lambdaPhiSub == 'F_invar':
        ws.set('POI').add(ws.var('F_invar_p'))
        ws.set('POI').add(ws.var('F_invar_np'))

        ws.factory('RooUniform::prior_F_'+options.fitFrame+'_p(F_invar_p)')
        ws.factory('RooUniform::prior_F_'+options.fitFrame+'_np(F_invar_np)')
        
        PriorPOIList.add(ws.pdf('prior_F_'+options.fitFrame+'_p'))
        PriorPOIList.add(ws.pdf('prior_F_'+options.fitFrame+'_np'))
        
    if options.lambdaPhiSub == 'lambda_tilde':
        ws.set('POI').add(ws.var('lambda_tilde_p'))
        ws.set('POI').add(ws.var('lambda_tilde_np'))
        
        ws.factory('RooUniform::prior_ltilde_'+options.fitFrame+'_p(lambda_tilde_p)')
        ws.factory('RooUniform::prior_ltilde_'+options.fitFrame+'_np(lambda_tilde_np)')
        
        PriorPOIList.add(ws.pdf('prior_ltilde_'+options.fitFrame+'_p'))
        PriorPOIList.add(ws.pdf('prior_ltilde_'+options.fitFrame+'_np'))        
                 

    priorPOI = ROOT.RooProdPdf('PriorPOI','PriorPOI',PriorPOIList)

    getattr(ws,'import')(priorPOI)

    ws.factory('nPromptinP[500,0,1000000]')
    ws.factory('sum::nPromptinNP(1*nPromptSignal,-1*nPromptinP)')

    ws.factory('nNonPromptinNP[100,0,1000000]')
    ws.factory('sum::nNonPromptinP(1*nNonPromptSignal,-1*nNonPromptinNP)')
    
    ws.factory('nBackgroundinP[50,0,1000000]')
    ws.factory('sum::nBackgroundinNP(1*nBackgroundSignal,-1*nBackgroundinP)')

    ws.factory('RooPolarizationConstraint::promptConstr('+
               'lambda_theta_'+options.fitFrame+'_p,'+
               'lambda_phi_'+options.fitFrame+'_p,'+
               'lambda_thetaphi_'+options.fitFrame+'_p)')
    ws.factory('RooPolarizationConstraint::nonPromptConstr('+
               'lambda_theta_'+options.fitFrame+'_np,'+
               'lambda_phi_'+options.fitFrame+'_np,'+
               'lambda_thetaphi_'+options.fitFrame+'_np)')

    rap_bins = range(1,len(jpsi.pTRange))
    pt_bins = None
    
    if options.testBin is not None:
        rap_bins = [int(options.testBin.split(',')[0])]
        pt_bins  = [int(options.testBin.split(',')[1])-1]
      
    #for each rap/pt cell make a unique simultaneous fit of prompt,non-prompt, background
    for rap_bin in rap_bins:
        if options.testBin is None:
            pt_bins = range(len(jpsi.pTRange[rap_bin]))
        for pt_bin in pt_bins:
            #load in proper set of mass and lifetime fitted parameters
            ws.loadSnapshot('snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

            #prepare data-weighted average of lifetime description for integral-ly goodness
            dsigproj = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).reduce(RooArgSet(ws.var('JpsictErr')),
                                                                                   'massRegion==massRegion::signal')
            dsigproj.SetName('projdata_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
            getattr(ws,'import')(dsigproj)
                                                                                   
            #dsigproj.Print("v")
            ws.var('Jpsict').setRange('fNPIntegral',-1,.1) #change to rapidity dependent thingy later

            prInt = ws.pdf('promptLifetimeRaw').createIntegral(RooArgSet(ws.var('Jpsict')),
                                                               RooArgSet(ws.var('Jpsict')),
                                                               'fNPIntegral')
            
            fPInP = RooDataWeightedAverage('fPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                           'Fraction of Prompt in Prompt Region',
                                           prInt,
                                           ws.data('projdata_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                           RooArgSet(ws.var('JpsictErr')),1,False,False)
            
            prInP  = RooFormulaVar('nPromptInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   'nPromptInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   '@0*@1',
                                   RooArgList(fPInP,ws.var('nPromptSignal')))
            prInNP = RooFormulaVar('nPromptInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   'nPromptInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   '(1-@0)*@1',
                                   RooArgList(fPInP,ws.var('nPromptSignal')))

            

#            getattr(ws,'import')(prInP,ROOT.RooFit.RecycleConflictNodes())
#            getattr(ws,'import')(prInNP,ROOT.RooFit.RecycleConflictNodes())
            

            npInt = ws.pdf('nonPromptSSDRaw').createIntegral(RooArgSet(ws.var('Jpsict')),
                                                             RooArgSet(ws.var('Jpsict')),
                                                             'fNPIntegral')

            #npInt.Print("v")
            #npInt.getVal()

            fNPInP = RooDataWeightedAverage('fNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                            'Fraction of Prompt in Prompt Region',
                                            npInt,
                                            ws.data('projdata_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                            RooArgSet(ws.var('JpsictErr')),1,False,False)

            #fNPInP.Print("v")
            
            npInP  = RooFormulaVar('nNonPromptInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   'nNonPromptInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   '@0*@1',
                                   RooArgList(fNPInP,ws.var('nNonPromptSignal')))
            npInNP = RooFormulaVar('nNonPromptInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   'nNonPromptInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                   '(1-@0)*@1',
                                   RooArgList(fNPInP,ws.var('nNonPromptSignal')))

#            getattr(ws,'import')(npInP,ROOT.RooFit.RecycleConflictNodes())
#            getattr(ws,'import')(npInNP,ROOT.RooFit.RecycleConflictNodes())

            bkgInt = ws.pdf('backgroundRaw').createIntegral(RooArgSet(ws.var('Jpsict')),
                                                            RooArgSet(ws.var('Jpsict')),
                                                            'fNPIntegral')

            fBkgInP = RooDataWeightedAverage('fBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                             'Fraction of Prompt in Prompt Region',
                                             bkgInt,
                                             ws.data('projdata_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                             RooArgSet(ws.var('JpsictErr')),1,False,False)
            
            bkgInP  = RooFormulaVar('nBackgroundInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                    'nBackgroundInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                    '@0*@1',
                                    RooArgList(fBkgInP,ws.var('nBackgroundSignal')))
            bkgInNP = RooFormulaVar('nBackgroundInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                    'nBackgroundInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                    '(1-@0)*@1',
                                    RooArgList(fBkgInP,ws.var('nBackgroundSignal')))

#            getattr(ws,'import')(bkgInP,ROOT.RooFit.RecycleConflictNodes())
#            getattr(ws,'import')(bkgInNP,ROOT.RooFit.RecycleConflictNodes())
            
            fPinP_actual  = RooFormulaVar('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                          'True fPinP',
                                          '@0/(@0+@1+@2)', #
                                          RooArgList(prInP,
                                                     npInP,
                                                     bkgInP))
            fPinNP_actual = RooFormulaVar('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                          'True fPinNP',
                                          '0',#@0/(@0+@1+@2)
                                          RooArgList(prInNP,
                                                     npInNP,
                                                     bkgInNP))
            
            fNPinP_actual  = RooFormulaVar('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                           'True fNPinP',
                                           '@1/(@0+@1+@2)',#
                                           RooArgList(prInP,
                                                      npInP,
                                                      bkgInP))            
            fNPinNP_actual = RooFormulaVar('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                           'True fNPinNP',
                                           '@1/(@0+@1+@2)',#
                                           RooArgList(prInNP,
                                                      npInNP,
                                                      bkgInNP))
            
            fBkginP_actual  = RooFormulaVar('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                            'True fBkginP',
                                            '@2/(@0+@1+@2)',
                                            RooArgList(prInP,
                                                       npInP,
                                                       bkgInP))
            fBkginNP_actual = RooFormulaVar('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                            'True fBkginNP',
                                            '@2/(@0+@1+@2)',
                                            RooArgList(prInNP,
                                                       npInNP,
                                                       bkgInNP))
            #get some doubles
            nPinP = prInP.getVal()
            nNPinP = npInP.getVal()
            nBGinP = bkgInP.getVal()
            
            nPinPerr = ws.var('nPromptSignal').getError()*fPinP_actual.getVal()
            nNPinPerr = ws.var('nNonPromptSignal').getError()*fNPinP_actual.getVal()
            nBGinPerr = ws.var('nBackgroundSignal').getError()*fBkginP_actual.getVal()
            
            nPinNP = prInNP.getVal()
            nNPinNP = npInNP.getVal()
            nBGinNP = bkgInNP.getVal()
            
            nPinNPerr = ws.var('nPromptSignal').getError()*fPinNP_actual.getVal()
            nNPinNPerr = ws.var('nNonPromptSignal').getError()*fNPinNP_actual.getVal()
            nBGinNPerr = ws.var('nBackgroundSignal').getError()*fBkginNP_actual.getVal()

            PinP = fPinP_actual.getVal()
            NPinP = fNPinP_actual.getVal()
            BGinP = fBkginP_actual.getVal()

            PinNP = fPinNP_actual.getVal()
            NPinNP = fNPinNP_actual.getVal()
            BGinNP = fBkginNP_actual.getVal()

            #create errors, corrected alphas and betas for each fraction 
            errPinP  = ( (1/(nPinP+nNPinP+nBGinP) + nPinP/pow(nPinP+nNPinP+nBGinP,2))*nPinPerr +
                         nPinP/pow(nPinP+nNPinP+nBGinP,2)*nNPinPerr + nPinP/pow(nPinP+nNPinP+nBGinP,2)*nBGinPerr )
            
            alphaPinP = PinP*((PinP*(1-PinP))/errPinP/errPinP - 1) 
            betaPinP = (1-PinP)*((PinP*(1-PinP))/errPinP/errPinP - 1) + 1
            
            errNPinP = ( (1/(nPinP+nNPinP+nBGinP) + nNPinP/pow(nPinP+nNPinP+nBGinP,2))*nNPinPerr +
                         nNPinP/pow(nPinP+nNPinP+nBGinP,2)*nPinPerr + nNPinP/pow(nPinP+nNPinP+nBGinP,2)*nBGinPerr )

            alphaNPinP = NPinP*((NPinP*(1-NPinP))/errNPinP/errNPinP - 1) 
            betaNPinP = (1-NPinP)*((NPinP*(1-NPinP))/errNPinP/errNPinP - 1) + 1
            
            errBGinP = ( (1/(nPinP+nNPinP+nBGinP) + nBGinP/pow(nPinP+nNPinP+nBGinP,2))*nBGinPerr +
                         nBGinP/pow(nPinP+nNPinP+nBGinP,2)*nPinPerr + nBGinP/pow(nPinP+nNPinP+nBGinP,2)*nNPinPerr )
            
            alphaBGinP = BGinP*((BGinP*(1-BGinP))/errBGinP/errBGinP - 1) 
            betaBGinP = (1-BGinP)*((BGinP*(1-BGinP))/errBGinP/errBGinP - 1) + 1
            
            errPinNP  = ( (1/(nPinNP+nNPinNP+nBGinNP) + nPinNP/pow(nPinNP+nNPinNP+nBGinNP,2))*nPinNPerr +
                          nPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nNPinNPerr + nPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nBGinNPerr )
            
            alphaPinNP = PinNP*((PinNP*(1-PinNP))/errPinNP/errPinNP - 1) 
            betaPinNP = (1-PinNP)*((PinNP*(1-PinNP))/errPinNP/errPinNP - 1) + 1
            
            errNPinNP = ( (1/(nPinNP+nNPinNP+nBGinNP) + nNPinNP/pow(nPinNP+nNPinNP+nBGinNP,2))*nNPinNPerr +
                          nNPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nPinNPerr + nNPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nBGinNPerr )
            
            alphaNPinNP = NPinNP*((NPinNP*(1-NPinNP))/errNPinNP/errNPinNP - 1) 
            betaNPinNP = (1-NPinNP)*((NPinNP*(1-NPinNP))/errNPinNP/errNPinNP - 1) + 1
            
            errBGinNP = ( (1/(nPinNP+nNPinNP+nBGinNP) + nBGinNP/pow(nPinNP+nNPinNP+nBGinNP,2))*nBGinNPerr +
                          nBGinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nNPinNPerr + nBGinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nPinNPerr )
        
            alphaBGinNP = BGinNP*((BGinNP*(1-BGinNP))/errBGinNP/errBGinNP - 1)
            betaBGinNP = (1-BGinNP)*((BGinNP*(1-BGinNP))/errBGinNP/errBGinNP - 1) + 1

            #P enriched region
            ws.factory('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[.5]')
            ws.var('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(fPinP_actual.getVal())
            ws.var('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setError(errPinP)
            ws.var('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            ws.factory('alphaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            ws.factory('betaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            if(alphaPinP > 0 and betaPinP > 0):
                ws.var('alphaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(alphaPinP)
                ws.var('betaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(betaPinP)
            ws.var('alphaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            ws.var('betaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
                        
            ws.factory('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[.5]')
            ws.var('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(fNPinP_actual.getVal())
            ws.var('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setError(errNPinP)
            ws.var('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            ws.factory('alphaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            ws.factory('betaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            if(alphaNPinP > 0 and betaNPinP > 0):
                ws.var('alphaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(alphaNPinP)
                ws.var('betaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(betaNPinP)
            ws.var('alphaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            ws.var('betaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            
            ws.factory('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[.5]')
            ws.var('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(fBkginP_actual.getVal())
            ws.var('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setError(errBGinP)
            ws.var('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            ws.factory('alphaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            ws.factory('betaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            if(alphaNPinP > 0 and betaNPinP > 0):
                ws.var('alphaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(alphaBGinP)
                ws.var('betaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(betaBGinP)
            ws.var('alphaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            ws.var('betaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            
            #NP enriched region
            ws.factory('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[.5]')
            ws.var('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(fPinNP_actual.getVal())
            ws.var('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setError(errPinNP)
            ws.var('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            ws.factory('alphaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            ws.factory('betaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            if(alphaPinNP > 0 and betaPinNP > 0):
                ws.var('alphaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(alphaPinNP)
                ws.var('betaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(betaPinNP)
            ws.var('alphaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            ws.var('betaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            
            ws.factory('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[.5]')
            ws.var('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(fNPinNP_actual.getVal())
            ws.var('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setError(errNPinNP)
            ws.var('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            ws.factory('alphaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            ws.factory('betaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            if(alphaNPinNP > 0 and betaNPinNP > 0):
                ws.var('alphaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(alphaNPinNP)
                ws.var('betaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(betaNPinNP)
            ws.var('alphaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            ws.var('betaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            
            ws.factory('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[.5]')
            ws.var('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(fBkginNP_actual.getVal())
            ws.var('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setError(errBGinNP)
            ws.var('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            ws.factory('alphaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            ws.factory('betaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+'[1]')
            if(alphaBGinNP > 0 and betaBGinNP > 0):
                ws.var('alphaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(alphaBGinNP)
                ws.var('betaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).setVal(betaBGinNP)
            ws.var('alphaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()
            ws.var('betaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).Print()

            #create constraint function
            #prompt constraints
            x_PinP = RooRealVar('x_PinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                'x_PinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                fPinP_actual.getVal(),
                                max(fPinP_actual.getVal()-3*errPinP,0),
                                min(fPinP_actual.getVal()+3*errPinP,1))            
            BetaPinP = ROOT.makeBetaPdf('BetaFunPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                        x_PinP,
                                        ws.var('alphaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                        ws.var('betaPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))

            getattr(ws,'import')(BetaPinP)

            x_NPinP = RooRealVar('x_NPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                 'x_NPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                 fNPinP_actual.getVal(),
                                 max(fNPinP_actual.getVal()-3*errNPinP,0),
                                 min(fNPinP_actual.getVal()+3*errNPinP,1))            
            BetaNPinP = ROOT.makeBetaPdf('BetaFunNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                         x_NPinP,
                                         ws.var('alphaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                         ws.var('betaNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            
            getattr(ws,'import')(BetaNPinP)

            x_BGinP = RooRealVar('x_BkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                 'x_BkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                 fBkginP_actual.getVal(),
                                 max(fBkginP_actual.getVal()-3*errBGinP,0),
                                 min(fBkginP_actual.getVal()+3*errBGinP,1))            
            BetaBGinP = ROOT.makeBetaPdf('BetaFunBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                         x_BGinP,
                                         ws.var('alphaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                         ws.var('betaBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            
            getattr(ws,'import')(BetaBGinP)
            
            #non-prompt constraints
            x_PinNP = RooRealVar('x_PinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                 'x_PinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                 fPinNP_actual.getVal(),
                                 max(fPinNP_actual.getVal()-3*errPinNP,0),
                                 min(fPinNP_actual.getVal()+3*errPinNP,1))            
            BetaPinNP = ROOT.makeBetaPdf('BetaFunPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                         x_PinNP,
                                         ws.var('alphaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                         ws.var('betaPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            BetaPinNP.Print()
            getattr(ws,'import')(BetaPinNP)
            
            x_NPinNP = RooRealVar('x_NPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                  'x_NPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                  fNPinNP_actual.getVal(),
                                  max(fNPinNP_actual.getVal()-3*errNPinNP,0),
                                  min(fNPinNP_actual.getVal()+3*errNPinNP,1))            
            BetaNPinNP = ROOT.makeBetaPdf('BetaFunNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                         x_NPinNP,
                                         ws.var('alphaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                         ws.var('betaNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            BetaNPinNP.Print()
            getattr(ws,'import')(BetaNPinNP)
            
            x_BGinNP = RooRealVar('x_BkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                  'x_BkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                  fBkginNP_actual.getVal(),
                                  max(fBkginNP_actual.getVal()-3*errBGinNP,0),
                                  min(fBkginNP_actual.getVal()+3*errBGinNP,1))            
            BetaBGinNP = ROOT.makeBetaPdf('BetaFunBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                          x_BGinNP,
                                          ws.var('alphaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                          ws.var('betaBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))

            getattr(ws,'import')(BetaBGinNP)
            
            ListPriorNuis = RooArgList(ws.pdf('BetaFunPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                      ws.pdf('BetaFunNPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                      ws.pdf('BetaFunBkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                      ws.pdf('BetaFunPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                      ws.pdf('BetaFunNPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                      ws.pdf('BetaFunBkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            
            priorNuis = ROOT.RooProdPdf('PriorNuis_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                        'PriorNuis_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                                        ListPriorNuis)

            getattr(ws,'import')(priorNuis)
            
            #getattr(ws,'import')(fPinP_actual,ROOT.RooFit.RecycleConflictNodes())
            #getattr(ws,'import')(fPinNP_actual,ROOT.RooFit.RecycleConflictNodes())
            #getattr(ws,'import')(fNPinP_actual,ROOT.RooFit.RecycleConflictNodes())
            #getattr(ws,'import')(fNPinNP_actual,ROOT.RooFit.RecycleConflictNodes())
            #getattr(ws,'import')(fBkginP_actual,ROOT.RooFit.RecycleConflictNodes())
            #getattr(ws,'import')(fBkginNP_actual,ROOT.RooFit.RecycleConflictNodes())
            
            #make acceptance maps for prompt
            accMapPHist = accMapsP.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            accMapPHist.SetName('hAcc2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            recoEffMapPHist = recoEffMapsP.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            recoEffMapPHist.SetName('hRecoEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            trigEffMapPHist = trigEffMapsP.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            trigEffMapPHist.SetName('hTrigEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            
            getattr(ws,'import')(accMapPHist)
            getattr(ws,'import')(recoEffMapPHist)
            getattr(ws,'import')(trigEffMapPHist)
            
            accMap = RooDataHist('accMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                 'acceptance map',
                                 RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                 ROOT.RooFit.Import(ws.obj('hAcc2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(accMap)

            recoEffMap = RooDataHist('recoEffMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'reconstruction efficiency map',
                                     RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ROOT.RooFit.Import(ws.obj('hRecoEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(recoEffMap)

            trigEffMap = RooDataHist('trigEffMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'acceptance map',
                                     RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ROOT.RooFit.Import(ws.obj('hTrigEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(trigEffMap)
                        
            accMapFunc = RooHistFunc('accMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'acceptance map',
                                     RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ws.data('accMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(accMapFunc)

            recoEffMapFunc = RooHistFunc('recoEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                         'reconstruction efficiency map',
                                         RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                         ws.data('recoEffMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(recoEffMapFunc)

            trigEffMapFunc = RooHistFunc('trigEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                         'trigger efficiency map',
                                         RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                         ws.data('trigEffMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(trigEffMapFunc)

            accEffMapFuncTemp = RooProduct('accEffMapFuncPTemp'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                           'acceptance*efficiency map',
                                           RooArgSet(ws.function('accMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)),
                                                     ws.function('recoEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)),
                                                     ws.function('trigEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))))

            accEffMapHist = accEffMapFuncTemp.createHistogram('hAccEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin),
                                                              ws.var('costh_'+options.fitFrame+foldname),
                                                              ROOT.RooFit.Binning(accMapPHist.GetNbinsX()),
                                                              ROOT.RooFit.YVar(ws.var('phi_'+options.fitFrame+foldname),
                                                                               ROOT.RooFit.Binning(accMapPHist.GetNbinsY())),
                                                              ROOT.RooFit.Scaling(False))
            accEffMapHist.SetName('hAccEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            getattr(ws,'import')(accEffMapHist)

            accEffMapData = RooDataHist('accEffMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                        'acceptance*efficiency map',
                                        RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                        ROOT.RooFit.Import(ws.obj('hAccEff2D_P_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(accEffMapData)
            
            accEffMapFunc = RooHistFunc('accEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                        'acceptance*efficiency map',
                                        RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                        ws.data('accEffMapP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(accEffMapFunc)
            
            ws.factory('expr::accEffMapFuncPTheta'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@2(@0,@1)*@0*@0",'+
                       'costh_'+options.fitFrame+foldname+','
                       'phi_'+options.fitFrame+foldname+','
                       'accEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            ws.factory('expr::accEffMapFuncPPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@2(@0,@1)*(1.-@0*@0)*cos(@1*2.*pi/180.)",'+
                       'costh_'+options.fitFrame+foldname+','
                       'phi_'+options.fitFrame+foldname+','
                       'accEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            ws.factory('expr::accEffMapFuncPThetaPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@2(@0,@1)*2.*sqrt(1-@0*@0)*@0*cos(@1*pi/180.)",'+
                       'costh_'+options.fitFrame+foldname+','
                       'phi_'+options.fitFrame+foldname+','
                       'accEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')

            #make acceptance maps for non-prompt
            accMapNPHist = accMapsNP.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            accMapNPHist.SetName('hAcc2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            recoEffMapNPHist = recoEffMapsNP.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            recoEffMapNPHist.SetName('hRecoEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            trigEffMapNPHist = trigEffMapsNP.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            trigEffMapNPHist.SetName('hTrigEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            
            getattr(ws,'import')(accMapNPHist)
            getattr(ws,'import')(recoEffMapNPHist)
            getattr(ws,'import')(trigEffMapNPHist)
            
            accMap = RooDataHist('accMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                 'acceptance map',
                                 RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                 ROOT.RooFit.Import(ws.obj('hAcc2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(accMap)

            recoEffMap = RooDataHist('recoEffMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'reconstruction efficiency map',
                                     RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ROOT.RooFit.Import(ws.obj('hRecoEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(recoEffMap)

            trigEffMap = RooDataHist('trigEffMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'acceptance map',
                                     RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ROOT.RooFit.Import(ws.obj('hTrigEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(trigEffMap)
                        
            accMapFunc = RooHistFunc('accMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'acceptance map',
                                     RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ws.data('accMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(accMapFunc)

            recoEffMapFunc = RooHistFunc('recoEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                         'reconstruction efficiency map',
                                         RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                         ws.data('recoEffMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(recoEffMapFunc)

            trigEffMapFunc = RooHistFunc('trigEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                         'trigger efficiency map',
                                         RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                         ws.data('trigEffMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(trigEffMapFunc)

            accEffMapFuncTemp = RooProduct('accEffMapFuncNPTemp'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                           'acceptance*efficiency map',
                                           RooArgSet(ws.function('accMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)),
                                                     ws.function('recoEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)),
                                                     ws.function('trigEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))))

            accEffMapHist = accEffMapFuncTemp.createHistogram('hAccEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin),
                                                              ws.var('costh_'+options.fitFrame+foldname),
                                                              ROOT.RooFit.Binning(accMapNPHist.GetNbinsX()),
                                                              ROOT.RooFit.YVar(ws.var('phi_'+options.fitFrame+foldname),
                                                                               ROOT.RooFit.Binning(accMapNPHist.GetNbinsY())),
                                                              ROOT.RooFit.Scaling(False))
            accEffMapHist.SetName('hAccEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            getattr(ws,'import')(accEffMapHist)

            accEffMapData = RooDataHist('accEffMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                        'acceptance*efficiency map',
                                        RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                        ROOT.RooFit.Import(ws.obj('hAccEff2D_NP_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin)),False))
            getattr(ws,'import')(accEffMapData)
            
            accEffMapFunc = RooHistFunc('accEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                        'acceptance*efficiency map',
                                        RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                        ws.data('accEffMapNP'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),0)
            getattr(ws,'import')(accEffMapFunc)
                        
            ws.factory('expr::accEffMapFuncNPTheta'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@2(@0,@1)*@0*@0",'+
                       'costh_'+options.fitFrame+foldname+','
                       'phi_'+options.fitFrame+foldname+','
                       'accEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            ws.factory('expr::accEffMapFuncNPPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@2(@0,@1)*(1.-@0*@0)*cos(@1*2.*pi/180.)",'+
                       'costh_'+options.fitFrame+foldname+','
                       'phi_'+options.fitFrame+foldname+','
                       'accEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            ws.factory('expr::accEffMapFuncNPThetaPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@2(@0,@1)*2.*sqrt(1-@0*@0)*@0*cos(@1*pi/180.)",'+
                       'costh_'+options.fitFrame+foldname+','
                       'phi_'+options.fitFrame+foldname+','
                       'accEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            # END ACCEPTANCE MAPS ::::: THIS NEEDS TO BE TURNED INTO A FUNCTION

            #make a *real* hist pdf :-)
            ws.factory('RooGenericPdf::accMap'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@0",{accMapFunc'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+'})')
            
            #create datahist for L and R sidebands
            ws.var('costh_'+options.fitFrame).setBins(45,'bkg') #set bins for background shape!
            ws.var('phi_'+options.fitFrame).setBins(40,'bkg') #set bins for background shape!
            bkgL = RooDataHist('bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1),
                               'bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1),
                               RooArgSet(ws.var('costh_'+options.fitFrame),ws.var('phi_'+options.fitFrame)),
                               'bkg')
            bkgR = RooDataHist('bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1),
                               'bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1),
                               RooArgSet(ws.var('costh_'+options.fitFrame),ws.var('phi_'+options.fitFrame)),
                               'bkg')

            getattr(ws,'import')(bkgL)
            getattr(ws,'import')(bkgR)
            
            #ws.data('bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)).reset()
            #ws.data('bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)).reset()

            #fill them
            ws.data('bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)).add(ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                                                                             'massRegion == massRegion::leftMassSideBand')
            ws.data('bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)).add(ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                                                                             'massRegion == massRegion::rightMassSideBand')

            

            #make histpdfs and combination
            ws.factory('RooHistPdf::bkgShape'+options.fitFrame+'L_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '({costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+'},'+
                       'bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)+',1)')
            ws.factory('RooHistPdf::bkgShape'+options.fitFrame+'R_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '({costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+'},'+
                       'bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)+',1)')
            ws.factory('RooGenericPdf::bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@0*@1+(1-@0)*@2",{fBkgLR,bkgShape'+options.fitFrame+'L_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       ',bkgShape'+options.fitFrame+'R_'+str(rap_bin)+'_'+str(pt_bin+1)+'})')

            #test new polarization pdf
            ws.factory('RooPolarizationPdf::basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+
                       ',lambda_theta_'+options.fitFrame+'_p,lambda_phi_'+options.fitFrame+'_p,lambda_thetaphi_'+options.fitFrame+'_p,'+
                       'accEffMapFuncP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                       'accEffMapFuncPTheta'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                       'accEffMapFuncPPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                       'accEffMapFuncPThetaPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            
            ws.factory('RooPolarizationPdf::basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+
                       ',lambda_theta_'+options.fitFrame+'_np,lambda_phi_'+options.fitFrame+'_np,lambda_thetaphi_'+options.fitFrame+'_np,'+
                       'accEffMapFuncNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                       'accEffMapFuncNPTheta'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                       'accEffMapFuncNPPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                       'accEffMapFuncNPThetaPhi'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')


            prF = ws.pdf('basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1))
            npF = ws.pdf('basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1))

            print prF.getVal(RooArgSet(ws.var('costh_'+options.fitFrame),ws.var('phi_'+options.fitFrame)))
            print npF.getVal(RooArgSet(ws.var('costh_'+options.fitFrame),ws.var('phi_'+options.fitFrame)))

            print prF.createIntegral(RooArgSet(ws.var('costh_'+options.fitFrame),ws.var('phi_'+options.fitFrame))).getVal()
            print npF.createIntegral(RooArgSet(ws.var('costh_'+options.fitFrame),ws.var('phi_'+options.fitFrame))).getVal()

            #should add back *(3+@2) ?
            #sin(2theta) = 2sin(theta)cos(theta) = 2*sqrt(1 - cos(theta)*cos(theta))*cos(theta)
            #ws.factory('RooGenericPdf::basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '("(1+@2*pow(@0,2.)+@3*(1-pow(@0,2.))*cos(2*@1*pi/180)+@4*sin(2*acos(@0))*cos(@1*pi/180))"'+
            #           ',{costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+',lambda_theta_'+options.fitFrame+
            #           '_p,lambda_phi_'+options.fitFrame+'_p,lambda_thetaphi_'+options.fitFrame+'_p})')
            
            #ws.factory('RooGenericPdf::basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '("(1+@2*@0*@0+@3*(1-@0*@0)*cos(2*@1*pi/180)+@4*sin(2*acos(@0))*cos(@1*pi/180))"'+
            #           ',{costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+',lambda_theta_'+options.fitFrame+
            #           '_np,lambda_phi_'+options.fitFrame+'_np,lambda_thetaphi_'+options.fitFrame+'_np})')    

            ws.factory('PROD::polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +'(basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +',promptConstr)') #add in constraint functions as multiplicative penalties to PDF
            
            ws.factory('PROD::polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +'(basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +',nonPromptConstr)') #add in constraint function as multiplicative penalties to PDF
            
            #ws.factory('RooExtendPdf::promptPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '(polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nPromptInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+')')

            #ws.factory('RooExtendPdf::promptPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '(polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nPromptInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+')')

            #ws.factory('RooExtendPdf::nonPromptPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '(polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nNonPromptInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+')')
            
            #ws.factory('RooExtendPdf::nonPromptPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '(polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nNonPromptInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+')')
            
            #ws.factory('RooExtendPdf::bkgPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '(bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',nBackgroundInPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+')')

            #ws.factory('RooExtendPdf::bkgPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '(bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',nBackgroundInNonPrompt_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)+')')

            promptPOLArgList = RooArgList(ws.pdf('polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)))
            promptPOLFracList = RooArgList(ws.var('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))

            nonPromptPOLArgList = RooArgList(ws.pdf('polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)),
                                             ws.pdf('polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)))
            nonPromptPOLFracList = RooArgList(ws.var('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                              ws.var('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
                        
            if options.doNonPrompt:
                promptPOLArgList.add(ws.pdf('polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)))
                promptPOLFracList.add(ws.var('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            if options.doBackground:        
                promptPOLArgList.add(ws.pdf('bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)))
                nonPromptPOLArgList.add(ws.pdf('bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)))
                promptPOLFracList.add(ws.var('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
                nonPromptPOLFracList.add(ws.var('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            
            promptP = RooAddPdf('promptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                'prompt',
                                promptPOLArgList,
                                promptPOLFracList
                                )
            
            nonpromptP = RooAddPdf('nonpromptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                   'nonprompt',
                                   nonPromptPOLArgList,
                                   nonPromptPOLFracList
                                   )

            getattr(ws,'import')(promptP)
            getattr(ws,'import')(nonpromptP)   

            dataSize = str(ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).numEntries())
            halfSize = str(ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).numEntries()/2.0)
            
            ws.factory('RooExtendPdf::promptPolExtended'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(promptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',nPromptPol['+halfSize+',0,'+dataSize+'])')

            ws.factory('RooExtendPdf::nonpromptPolExtended'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(nonpromptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',nNonPromptPol['+halfSize+',0,'+dataSize+'])')

            params = RooArgSet()

            params.add(ws.var('nPromptPol'))
            params.add(ws.var('nNonPromptPol'))
            params.add(ws.var('x_PinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            params.add(ws.var('x_NPinP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            params.add(ws.var('x_BkginP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            params.add(ws.var('x_PinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            params.add(ws.var('x_NPinNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            params.add(ws.var('x_BkginNP_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)))
            ws.defineSet('param_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),
                         params,
                         True)   

            polConfig = 'SIMUL::PPdf'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+'(mlRegion,'
             
            if options.doPrompt:
                polConfig += 'promptSignal=promptPolExtended'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',' #
            if options.doNonPrompt:
                polConfig += 'nonPromptSignal=nonpromptPolExtended'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1) #
            #if options.doBackground:       
            #    polConfig += ('prompt;leftMassSideBand=bkgPolExtL'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','
            #                  'prompt;rightMassSideBand=bkgPolExtR'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))
            #    polConfig += ('nonPrompt;leftMassSideBand=bkgPolExtL'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','
            #                  'nonPrompt;rightMassSideBand=bkgPolExtR'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))
                 

            polConfig += ')'
            
            #simultaneous
            ws.factory(polConfig)

            print "made pdf!"

            #ws.pdf('polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)).setNormValueCaching(3)
            #ws.pdf('polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)).setNormValueCaching(3)
            #normSet = RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname))
            #print 'HEY',ws.pdf('polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)).getVal(normSet)
            #print 'HEY',ws.pdf('polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)).getVal(normSet)

    accMapsP.Close()
    recoEffMapsP.Close()
    trigEffMapsP.Close()

    accMapsNP.Close()
    recoEffMapsNP.Close()
    trigEffMapsNP.Close()

if __name__ == '__main__':
    parser = OptionParser(description='%prog : J/Psi Polarization Fitter.',
                          usage='polarizationFit.py --workspaceName=foo inputData.root ...')
    parser.add_option('--workspaceName',dest='workspaceName',help='The name of your RooWorkspace')
    parser.add_option('--triggerName',dest='triggerName',help='Name of the Trigger you want to use.')
    parser.add_option('--treeName',dest='treeName',default='data',help='Name of the input TTree.')
    parser.add_option('--acceptanceMap',dest='acceptanceMap',help='Name of the file containing the acceptance maps.')
    parser.add_option('--recoEfficiencyMap',dest='recoEfficiencyMap',help='Name of the file containing the reco efficiency maps.')
    parser.add_option('--trigEfficiencyMap',dest='trigEfficiencyMap',help='Name of the file containing the trigger efficiency maps.')
    parser.add_option('--fitFrame',dest='fitFrame',help='CS,HX,PHX,sGJ,GJ1,GJ2')
    parser.add_option('--testBin',dest='testBin',help='order pair of rapidity and pt bin number for testing.')
    parser.add_option('--notfolded',dest='folded',default=True,action='store_false',help='use folded phi')
    parser.add_option('--lambdaPhiSub',dest='lambdaPhiSub',default='lambda_phi',help='substitute something for lambda phi')
    

    #define options to turn off pieces of final pdf
    parser.add_option('--noPrompt',dest='doPrompt',default=True,action='store_false',help='turn off prompt part of fit')
    parser.add_option('--noNonPrompt',dest='doNonPrompt',default=True,action='store_false',help='turn off nonprompt fit')
    parser.add_option('--noBackground',dest='doBackground',default=True,action='store_false',help='turn off background fit')
    parser.add_option('--polarizationOnly',dest='polOnly',default=False,action='store_true',help='turn off mass and lifetime fits')
    
    (options,args) = parser.parse_args()

    miss_options = False

    if options.workspaceName is None:
        print 'Need to specify --workspaceName'
        miss_options=True
#    if options.triggerName is None:
#        print 'Need to specify --triggerName'
#        miss_options=True
    if options.acceptanceMap is None and options.fitFrame is not None:
        print 'Need to specify --acceptanceMap'
        miss_options=True
    if options.recoEfficiencyMap is None and options.fitFrame is not None:
        print 'Need to specify --recoEfficiencyMap'
        miss_options=True
    if options.trigEfficiencyMap is None and options.fitFrame is not None:
        print 'Need to specify --trigEfficiencyMap'
        miss_options=True
        
    if miss_options:
        exit(1)

    main(options,args)
