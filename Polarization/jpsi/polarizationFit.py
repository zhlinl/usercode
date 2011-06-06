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
from ROOT import RooSuperCategory, RooMinimizer

def main(options,args):
    gROOT.Reset()

    #load our super special Polarization PDF
    gROOT.ProcessLine('.L RooPolarizationPdf.cxx+')
    gROOT.ProcessLine('.L RooPolarizationConstraint.cxx+')
    
    #setup integration
    intConf = ROOT.RooAbsReal.defaultIntegratorConfig()
    #intConf.Print('v')
#    intConf.method1D().setLabel('RooAdaptiveGaussKronrodIntegrator1D')
    intConf.setEpsAbs(1e-13)
    intConf.setEpsRel(1e-13)
    print intConf.epsAbs()
    print intConf.epsRel()
#    intConf.method2D().setLabel('RooMCIntegrator')
#    intConf.methodND().setLabel('RooMCIntegrator')

    output = TFile.Open(options.workspaceName+'.root','RECREATE')

    theWS = RooWorkspace(options.workspaceName,1)

    #save the polarization PDF code in the RooWorkspace
    theWS.importClassCode('RooPolarization*',True)
    
    buildDataAndCategories(theWS,options,args)    

    buildMassAndLifetimePDF(theWS)

#    if options.fitFrame is not None:
#        buildPolarizationPDF(theWS,options)

    #root is stupid
    output.cd()
    
    theWS.Print('v')

    ROOT.RooMsgService.instance().Print()

    doFit(theWS,options)
    
    theWS.Write()
    output.Close()

def doFit(ws,options):
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

            jPsiMass.setRange('mlfit_prompt',2.7,3.5)
            jPsiMass.setRange('mlfit_nonPrompt',2.7,3.5)            

            jPsiMass.setRange('NormalizationRangeFormlfit_prompt',2.7,3.5)
            jPsiMass.setRange('NormalizationRangeFormlfit_nonPrompt',2.7,3.5)

            jPsicTau.setRange('mlfit_signal',-1,2.5)
            jPsicTau.setRange('mlfit_leftMassSideBand',-1,2.5)
            jPsicTau.setRange('mlfit_rightMassSideBand',-1,2.5)

            jPsicTau.setRange('NormalizationRangeFormlfit_signal',-1,2.5)
            jPsicTau.setRange('NormalizationRangeFormlfit_leftMassSideBand',-1,2.5)
            jPsicTau.setRange('NormalizationRangeFormlfit_rightMassSideBand',-1,2.5)

            #jPsicTau.setRange('NormalizationRangeFormlfit_promptSignal',-1,.1)
            #jPsicTau.setRange('NormalizationRangeFormlfit_nonPromptSignal',.1,2.5)
            #jPsicTau.setRange('NormalizationRangeFormlfit_leftMassSideBand',-1,2.5)
            #jPsicTau.setRange('NormalizationRangeFormlfit_rightMassSideBand',-1,2.5)
            
            #jPsicTau.setRange('mlfit_promptSignal',-1,.1)
            #jPsicTau.setRange('mlfit_nonPromptSignal',.1,2.5)
            #jPsicTau.setRange('mlfit_leftMassSideBand',-1,2.5)
            #jPsicTau.setRange('mlfit_rightMassSideBand',-1,2.5)

            #reset parameters
            ws.var('CBn').setVal(.5)
            ws.var('CBalpha').setVal(.5)
            ws.var('CBmass').setVal(3.1)
            ws.var('CBsigma').setVal(.02)
            ws.var('bkgLambda').setVal(0)
            
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
            
            #ws.var('nPrompt').setVal(5000)
            #ws.var('nNonPrompt').setVal(500)
            #ws.var('nBackground').setVal(100)
            #ws.var('nBackgroundL').setVal(50)
            #ws.var('nBackgroundR').setVal(50)
            
            ws.var('nonPromptTau').setVal(.5)
            ws.var('promptMean').setVal(0)
            ws.var('ctResolution').setVal(1)

                        
            LPdf = ws.pdf('LPdf')
            MPdf = ws.pdf('MPdf')            
            
                
            data = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

            NLLs = RooArgSet()

            MassNLL = MPdf.createNLL(data,
                                     ROOT.RooFit.Range('mlfit'),
                                     ROOT.RooFit.SplitRange(True),                                    
                                     ROOT.RooFit.ConditionalObservables(RooArgSet(ws.var('JpsictErr'))),
                                     ROOT.RooFit.NumCPU(2))
                                     

            CTauNLL = LPdf.createNLL(data,
                                     ROOT.RooFit.Range('mlfit'),
                                     ROOT.RooFit.SplitRange(True),                                    
                                     ROOT.RooFit.ConditionalObservables(RooArgSet(ws.var('JpsictErr'))),
                                     ROOT.RooFit.NumCPU(2))
                
            NLLs.add(MassNLL)
            NLLs.add(CTauNLL)

            simNLL = RooAddition('add','add',NLLs)

            minuit = RooMinuit(simNLL)
            minuit.setStrategy(2)
            minuit.setPrintEvalErrors(-1)

            minuit.simplex()
            minuit.migrad()
            minuit.migrad()
            minuit.hesse()
            
            fitresult = minuit.save('polfitresult_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
            getattr(ws,'import')(fitresult)           

            ws.saveSnapshot('snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1),ws.allVars())
            

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
    jPsiVprob     = RooRealVar('JpsiVprob','',.01,1)
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

    dataVars.add(jPsiVprob)
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
    ws.factory('PROD::backgroundRawMass(bkgMassShape)') #,backgroundDSDRaw|JpsictErr

    ws.factory('PROD::promptRawCTau(promptLifetimeRaw|JpsictErr)') 
    ws.factory('PROD::nonpromptRawCTau(nonPromptSSDRaw|JpsictErr)') 
    ws.factory('PROD::backgroundRawCTauL(backgroundRawL|JpsictErr)') 
    ws.factory('PROD::backgroundRawCTauR(backgroundRawR|JpsictErr)') 
    ws.factory('PROD::backgroundRawCTau(backgroundRaw|JpsictErr)')
        
    #extended pdfs    
    ws.factory('RooExtendPdf::promptExtMass(promptRawMass,sum::nPrompt(nPromptSignal[1000,0,1000000]))') #,nPromptL[1000,0,1000000],nPromptR[1000,0,1000000]
    ws.factory('RooExtendPdf::nonpromptExtMass(nonpromptRawMass,sum::nNonPrompt(nNonPromptSignal[500,0,1000000]))')
    #[500,0,1000000],,nNonPromptL[500,0,1000000],nNonPromptR[500,0,1000000]
    ws.factory('RooExtendPdf::backgroundExtMass(backgroundRawMass,sum::nBackground(nBackgroundSignal[100,0,1000000],nBackgroundL[100,0,1000000],nBackgroundR[100,0,1000000]))') #

    ws.factory('RooExtendPdf::promptExtCTauL(promptRawCTau,0)') #nPromptL[1,0,1000000]
    ws.factory('RooExtendPdf::promptExtCTauR(promptRawCTau,0)')#nPromptR[1,0,1000000]
    ws.factory('RooExtendPdf::nonpromptExtCTauL(nonpromptRawCTau,0)') #nNonPromptL[1,0,1000000]
    ws.factory('RooExtendPdf::nonpromptExtCTauR(nonpromptRawCTau,0)') #nNonPromptR[1,0,1000000]
    ws.factory('RooExtendPdf::backgroundExtCTauL(backgroundRawCTauL,nBackgroundL)')
    ws.factory('RooExtendPdf::backgroundExtCTauR(backgroundRawCTauR,nBackgroundR)')

    ws.factory('RooExtendPdf::promptExtCTau(promptRawCTau,nPromptSignal)')
    ws.factory('RooExtendPdf::nonpromptExtCTau(nonpromptRawCTau,nNonPromptSignal)')
    ws.factory('RooExtendPdf::backgroundExtCTau(backgroundRawCTau,nBackgroundSignal)')
    
    
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
                                                                        ws.pdf('nonpromptExtCTauR'),
                                                                        ws.pdf('backgroundExtCTauR')))   

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
    
    if options.doPrompt:
        massConfig += 'prompt=promptMass,'
        ctConfig += 'signal=promptCTau,'
    if options.doNonPrompt:
        massConfig += 'nonPrompt=promptMass'
    if options.doBackground:               
        ctConfig += 'leftMassSideBand=backgroundCTauL,rightMassSideBand=backgroundCTauR'

    massConfig += ')'
    ctConfig += ')'
    
    #simultaneous
    ws.factory(massConfig)
    ws.factory(ctConfig)

def buildPolarizationPDF(ws,options):

    foldname = ''
    if not options.folded:
        foldname=''
        
    accMaps = TFile.Open(options.acceptanceMap)
    #effMaps = TFile.Open(options.efficiencyMap)

    ws.factory('lambda_theta_'+options.fitFrame+'_p[0,-1,1]')
    ws.factory('lambda_phi_'+options.fitFrame+'_p[0,-1,1]')
    ws.factory('lambda_thetaphi_'+options.fitFrame+'_p[0,-1,1]')

    ws.factory('lambda_theta_'+options.fitFrame+'_np[0,-1,1]')
    ws.factory('lambda_phi_'+options.fitFrame+'_np[0,-1,1]')
    ws.factory('lambda_thetaphi_'+options.fitFrame+'_np[0,-1,1]')

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
      
    #for each rap/pt cell make a unique simultaneous fit of prompt,non-prompt, background
    for rap_bin in range(1,len(jpsi.pTRange)):
        for pt_bin in range(len(jpsi.pTRange[rap_bin])):

            accMapHist = accMaps.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
            for effMap in options.efficiencyMap.split(','):
                effMaps = TFile.Open(effMap)
                effMapHist = effMaps.Get('hAcc2D_'+options.fitFrame+'_pT'+str(pt_bin+1)+'_rap'+str(rap_bin))
                accMapHist.Multiply(effMapHist)
                effMaps.Close()
            
            getattr(ws,'import')(accMapHist)

            accMap = RooDataHist('accMap'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1),
                                 'acceptance map',
                                 RooArgList(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                 ROOT.RooFit.Import(accMapHist,False))
            getattr(ws,'import')(accMap)
                        
            accMapFunc = RooHistFunc('accMapFunc'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                     'acceptance map',
                                     RooArgSet(ws.var('costh_'+options.fitFrame+foldname),ws.var('phi_'+options.fitFrame+foldname)),
                                     ws.data('accMap'+options.fitFrame+'Data_'+str(rap_bin)+'_'+str(pt_bin+1)),1)
            getattr(ws,'import')(accMapFunc)

            #make a *real* hist pdf :-)
            ws.factory('RooGenericPdf::accMap'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '("@0",{accMapFunc'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+'})')
            
            #create datahist for L and R sidebands
            getattr(ws,'import')(accMap.Clone('bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)))
            getattr(ws,'import')(accMap.Clone('bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)))
            ws.data('bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)).reset()
            ws.data('bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)).reset()

            #fill them
            ws.data('bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)).add(ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                                                                             'massRegion == massRegion::leftMassSideBand')
            ws.data('bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)).add(ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)),
                                                                                             'massRegion == massRegion::rightMassSideBand')

            #make histpdfs and combination
            ws.factory('RooHistPdf::bkgShape'+options.fitFrame+'L_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '({costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+'},'+
                       'bkgShape'+options.fitFrame+'DataL_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            ws.factory('RooHistPdf::bkgShape'+options.fitFrame+'R_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '({costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+'},'+
                       'bkgShape'+options.fitFrame+'DataR_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            ws.factory('SUM::bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(fBkgLR*bkgShape'+options.fitFrame+'L_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       ',bkgShape'+options.fitFrame+'R_'+str(rap_bin)+'_'+str(pt_bin+1)+')')

            #test new polarization pdf
            ws.factory('RooPolarizationPdf::basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+
                       ',lambda_theta_'+options.fitFrame+'_p,lambda_phi_'+options.fitFrame+'_p,lambda_thetaphi_'+options.fitFrame+'_p)')
            ws.factory('RooPolarizationPdf::basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+
                       ',lambda_theta_'+options.fitFrame+'_np,lambda_phi_'+options.fitFrame+'_np,lambda_thetaphi_'+options.fitFrame+'_np)')

            #should add back *(3+@2) ?
            #sin(2theta) = 2sin(theta)cos(theta) = 2*sqrt(1 - cos(theta)*cos(theta))*cos(theta)
            #ws.factory('RooGenericPdf::basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '("1+@2*pow(@0,2.)+@3*(1-pow(@0,2.))*cos(2*@1*pi/180)+@4*sin(2*acos(@0))*cos(@1*pi/180)"'+
            #           ',{costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+',lambda_theta_'+options.fitFrame+
            #           '_p,lambda_phi_'+options.fitFrame+'_p,lambda_thetaphi_'+options.fitFrame+'_p})')
            
            #ws.factory('RooGenericPdf::basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+
            #           '("1+@2*@0*@0+@3*(1-@0*@0)*cos(2*@1*pi/180)+@4*sin(2*acos(@0))*cos(@1*pi/180)"'+
            #           ',{costh_'+options.fitFrame+foldname+',phi_'+options.fitFrame+foldname+',lambda_theta_'+options.fitFrame+
            #           '_np,lambda_phi_'+options.fitFrame+'_np,lambda_thetaphi_'+options.fitFrame+'_np})')    

            ws.factory('PROD::polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +'(basePolPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +',accMap'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            
            ws.factory('PROD::polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +'(basePolPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)
                       +',accMap'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
            
            ws.factory('RooExtendPdf::promptPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nPromptinP)')

            ws.factory('RooExtendPdf::promptPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(polPdf'+options.fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nPromptinNP)')

            ws.factory('RooExtendPdf::nonPromptPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nNonPromptinNP)')
            
            ws.factory('RooExtendPdf::nonPromptPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(polPdf'+options.fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1)+',nNonPromptinP)')
            
            ws.factory('RooExtendPdf::bkgPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',nBackgroundinP)')

            ws.factory('RooExtendPdf::bkgPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                       '(bkgShape'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+',nBackgroundinNP)')

            promptPOLArgList = RooArgList(ws.pdf('promptPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)))

            nonPromptPOLArgList = RooArgList(ws.pdf('promptPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)),
                                             ws.pdf('nonPromptPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))
                                             )
            
            if options.doNonPrompt:
                promptPOLArgList.add(ws.pdf('nonPromptPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)))
            if options.doBackground:        
                promptPOLArgList.add(ws.pdf('bkgPolExtinP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)))
                nonPromptPOLArgList.add(ws.pdf('bkgPolExtinNP'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)))
            
            promptP = RooAddPdf('promptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                'prompt',
                                promptPOLArgList
                                )
            
            nonpromptP = RooAddPdf('nonpromptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1),
                                   'nonprompt',
                                   nonPromptPOLArgList
                                   )
    
            getattr(ws,'import')(promptP)
            getattr(ws,'import')(nonpromptP)

            polConfig = 'SIMUL::PPdf'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+'(mlRegion,'
             
            if options.doPrompt:
                polConfig += 'promptSignal=promptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','
            if options.doNonPrompt:
                polConfig += 'nonPromptSignal=nonpromptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','
            #if options.doBackground:       
            #    polConfig += ('prompt;leftMassSideBand=bkgPolExtL'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','
            #                  'prompt;rightMassSideBand=bkgPolExtR'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))
            #    polConfig += ('nonPrompt;leftMassSideBand=bkgPolExtL'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','
            #                  'nonPrompt;rightMassSideBand=bkgPolExtR'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1))
                 

            polConfig += ')'
            
            #simultaneous
            ws.factory(polConfig)

    accMaps.Close()

if __name__ == '__main__':
    parser = OptionParser(description='%prog : J/Psi Polarization Fitter.',
                          usage='polarizationFit.py --workspaceName=foo inputData.root ...')
    parser.add_option('--workspaceName',dest='workspaceName',help='The name of your RooWorkspace')
    parser.add_option('--triggerName',dest='triggerName',help='Name of the Trigger you want to use.')
    parser.add_option('--treeName',dest='treeName',default='data',help='Name of the input TTree.')
    parser.add_option('--acceptanceMap',dest='acceptanceMap',help='Name of the file containing the acceptance maps.')
    parser.add_option('--efficiencyMap',dest='efficiencyMap',help='Name of the file containing the efficiency maps.')
    parser.add_option('--fitFrame',dest='fitFrame',help='CS,HX,PHX,sGJ,GJ1,GJ2')
    parser.add_option('--testBin',dest='testBin',help='order pair of rapidity and pt bin number for testing.')
    parser.add_option('--notfolded',dest='folded',default=True,action='store_false',help='use folded phi')

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
    if options.efficiencyMap is None and options.fitFrame is not None:
        print 'Need to specify --efficiencyMap'
        miss_options=True
        
    if miss_options:
        exit(1)

    main(options,args)
