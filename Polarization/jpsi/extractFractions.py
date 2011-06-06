#!/usr/bin/env python

# implementation of the full polarization fit using a RooSimultaneous
# should be much more stable than what we had previously
# nuisance parameters will be fit without polarization first and
# then constrained during polarization fit.

import os
from optparse import OptionParser
import ROOT
import commonVar as jpsi

from math import sqrt,fabs,exp
from random import random

from ROOT import RooFit, RooWorkspace, RooArgSet, RooArgList, RooRealVar,gROOT
from ROOT import RooCategory, RooFormulaVar, RooDataSet, RooAddPdf, TGraphAsymmErrors
from ROOT import TFile, TTree, TChain, RooGExpModel, RooAddition, RooMinuit, TVectorD
from ROOT import TH1F, TH2F
from ROOT import RooPlot



def main(options,args):
    gROOT.Reset()
    gROOT.Reset()

    #setup integration
    intConf = ROOT.RooAbsReal.defaultIntegratorConfig()
    #intConf.Print('v')
    #intConf.method1D().setLabel('RooAdaptiveGaussKronrodIntegrator1D')

    #for i in range(1000):
    #    print unitarityTriplet()
    
    for f in args:

        output = TFile.Open(f,'UPDATE')
        output.cd()

        print f[:f.find('.root')].split('-')[0]
        
        ws = output.Get(f[:f.find('.root')].split('-')[0])

        for rap_bin in range(1,len(jpsi.pTRange)):                    
            for pt_bin in range(len(jpsi.pTRange[rap_bin])):

                if options.testBin is not None:
                    if rap_bin != int(options.testBin[0]) or pt_bin+1 != int(options.testBin[1]):
                        continue

                ws.loadSnapshot('pol_snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                
                PinP  = ws.var('fPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).getVal()
                NPinP = ws.var('fNPinP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).getVal()
                BGinP = ws.var('fBkginP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).getVal()

                PinNP  = ws.var('fPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).getVal()
                NPinNP = ws.var('fNPinNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).getVal()
                BGinNP = ws.var('fBkginNP_actual_rap'+str(rap_bin)+'_pt'+str(pt_bin+1)).getVal()

                nPinP = ws.var('nPromptSignal').getVal()*PinP
                nNPinP = ws.var('nNonPromptSignal').getVal()*NPinP
                nBGinP = ws.var('nBackgroundSignal').getVal()*BGinP

                nPinPerr = ws.var('nPromptSignal').getError()*PinP
                nNPinPerr = ws.var('nNonPromptSignal').getError()*NPinP
                nBGinPerr = ws.var('nBackgroundSignal').getError()*BGinP

                nPinNP = ws.var('nPromptSignal').getVal()*PinNP
                nNPinNP = ws.var('nNonPromptSignal').getVal()*NPinNP
                nBGinNP = ws.var('nBackgroundSignal').getVal()*BGinNP

                nPinNPerr = ws.var('nPromptSignal').getError()*PinNP
                nNPinNPerr = ws.var('nNonPromptSignal').getError()*NPinNP
                nBGinNPerr = ws.var('nBackgroundSignal').getError()*BGinNP


                errPinP  = ( (1/(nPinP+nNPinP+nBGinP) + nPinP/pow(nPinP+nNPinP+nBGinP,2))*nPinPerr +
                             nPinP/pow(nPinP+nNPinP+nBGinP,2)*nNPinPerr + nPinP/pow(nPinP+nNPinP+nBGinP,2)*nBGinPerr )

                alphaPinP = 1 #PinP*((PinP*(1-PinP))/errPinP/errPinP - 1)
                betaPinP = 1 #(alphaPinP*(1-PinP) + 2*(PinP-.5))/PinP
                
                errNPinP = ( (1/(nPinP+nNPinP+nBGinP) + nNPinP/pow(nPinP+nNPinP+nBGinP,2))*nNPinPerr +
                             nNPinP/pow(nPinP+nNPinP+nBGinP,2)*nPinPerr + nNPinP/pow(nPinP+nNPinP+nBGinP,2)*nBGinPerr )

                alphaNPinP = 1 #NPinP*((NPinP*(1-NPinP))/errNPinP/errNPinP - 1)
                betaNPinP = 1 #(alphaNPinP*(1-NPinP) + 2*(NPinP-.5))/NPinP
                
                errBGinP = ( (1/(nPinP+nNPinP+nBGinP) + nBGinP/pow(nPinP+nNPinP+nBGinP,2))*nBGinPerr +
                             nBGinP/pow(nPinP+nNPinP+nBGinP,2)*nPinPerr + nBGinP/pow(nPinP+nNPinP+nBGinP,2)*nNPinPerr )

                alphaBGinP = 1 #BGinP*((BGinP*(1-BGinP))/errBGinP/errBGinP - 1)
                betaBGinP = 1 #(alphaBGinP*(1-BGinP) + 2*(BGinP-.5))/BGinP
                
                errPinNP  = ( (1/(nPinNP+nNPinNP+nBGinNP) + nPinNP/pow(nPinNP+nNPinNP+nBGinNP,2))*nPinNPerr +
                             nPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nNPinNPerr + nPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nBGinNPerr )

                alphaPinNP = 1 #PinNP*((PinNP*(1-PinNP))/errPinNP/errPinNP - 1)
                betaPinNP = 1 #(alphaPinNP*(1-PinNP) + 2*(PinNP-.5))/PinNP
                
                errNPinNP = ( (1/(nPinNP+nNPinNP+nBGinNP) + nNPinNP/pow(nPinNP+nNPinNP+nBGinNP,2))*nNPinNPerr +
                             nNPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nPinNPerr + nNPinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nBGinNPerr )

                alphaNPinNP = 1 #NPinNP*((NPinNP*(1-NPinNP))/errNPinNP/errNPinNP - 1)
                betaNPinNP = 1 #(1-NPinNP)*((NPinNP*(1-NPinNP))/errNPinNP/errNPinNP - 1)
                #if NPinNP*(1-NPinNP)/errNPinNP/errNPinNP <= 1 and NPinNP < .1:
                #    
                #    alphaNPinNP = 1
                #    betaNPinNP = (1-NPinNP)/NPinNP
                #elif  NPinNP*(1-NPinNP)/errNPinNP/errNPinNP <= 1 and NPinNP > .9:
                #    
                #    betaNPinNP = 1
                #    alphaNPinNP = NPinNP/(1-NPinNP)
                
                errBGinNP = ( (1/(nPinNP+nNPinNP+nBGinNP) + nBGinNP/pow(nPinNP+nNPinNP+nBGinNP,2))*nBGinNPerr +
                             nBGinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nNPinNPerr + nBGinNP/pow(nPinNP+nNPinNP+nBGinNP,2)*nPinNPerr )

                alphaBGinNP = 1 #BGinNP*((BGinNP*(1-BGinNP))/errBGinNP/errBGinNP - 1)
                betaBGinNP = 1 #(alphaBGinNP*(1-BGinNP) + 2*(BGinNP-.5))/BGinNP

                print 'rap = ',rap_bin,' pt = ', pt_bin+1
                print 'fPinP = ', PinP,' +- ', errPinP
                print 'fNPinP = ', NPinP,' +- ', errNPinP
                print 'fBGinP = ', BGinP,' +- ', errBGinP
                
                print 'fPinNP = ', PinNP,' +- ', errPinNP
                print 'fNPinNP = ', NPinNP,' +- ', errNPinNP
                print 'fBGinNP = ', BGinNP,' +- ', errBGinNP

                print
                
                print ' P in P alpha = ',alphaPinP,' beta = ',betaPinP
                print 'NP in P alpha = ',alphaNPinP,' beta = ',betaNPinP
                print 'BG in P alpha = ',alphaBGinP,' beta = ',betaBGinP

                print ' P in NP alpha = ',alphaPinNP,' beta = ',betaPinNP
                print 'NP in NP alpha = ',alphaNPinNP,' beta = ',betaNPinNP
                print 'BG in NP alpha = ',alphaBGinNP,' beta = ',betaBGinNP

                print

                d = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                dpr = d.reduce('mlRegion == mlRegion::promptSignal')
                dnp = d.reduce('mlRegion == mlRegion::nonPromptSignal')

                print 'N_Events (l < .1 mm): ',dpr.sumEntries()
                print 'N_Events (l > .1 mm): ',dnp.sumEntries()


        output.Close()


if __name__ == '__main__':
    parser = OptionParser(description='%prog : J/Psi Polarization Fitter.',
                          usage='polarizationFit.py --workspaceName=foo inputData.root ...')
    
    parser.add_option('--treeName',dest='treeName',default='data',help='Name of the input TTree.')
    parser.add_option('--plotPol',dest='plotPol',default=False,action='store_true',help='Make plots of polarization fit')
    parser.add_option('--testBin',dest='testBin',help='Only plot one bin.')
    parser.add_option('--points',dest='points',help='number of likelihood points to test')
    parser.add_option('--fitFrame',dest='fitFrame',help='HX,CS,....')
    (options,args) = parser.parse_args()

    miss_options = False

    if options.testBin is not None:
        options.testBin = options.testBin.split(',')
    if options.points is not None:
        options.points = int(options.points)
        
    if miss_options:
        exit(1)

    main(options,args)
