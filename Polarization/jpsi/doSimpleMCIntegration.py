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

    #load our super special Polarization PDF
    gROOT.ProcessLine('.L RooPolarizationPdf.cxx+')
    gROOT.ProcessLine('.L RooPolarizationConstraint.cxx+')

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

    #intConf.getConfigSection('RooIntegrator1D').Print('v')

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

                th_p_PostProb = TH1F(options.fitFrame+'postprob_lth_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                      'postprob_lth_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                      100,-1,1)
                phi_p_PostProb = TH1F(options.fitFrame+'postprob_lphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                      'postprob_lphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                      100,-1,1)
                thphi_p_PostProb = TH1F(options.fitFrame+'postprob_lthphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                        'postprob_lthphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                        100,-sqrt(2)/2,sqrt(2)/2)
                ltilde_p_PostProb = TH1F(options.fitFrame+'postprob_ltilde_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                         'postprob_ltilde_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                         100,-1,1)
                
                th_np_PostProb = TH1F(options.fitFrame+'postprob_lth_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                      'postprob_lth_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                      100,-1,1)
                phi_np_PostProb = TH1F(options.fitFrame+'postprob_lphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                       'postprob_lphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                       100,-1,1)
                thphi_np_PostProb = TH1F(options.fitFrame+'postprob_lthphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                         'postprob_lthphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                         100,-sqrt(2)/2,sqrt(2)/2)
                ltilde_np_PostProb = TH1F(options.fitFrame+'postprob_ltilde_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                          'postprob_ltilde_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                          100,-1,1)
                
                th_p_NLL = TH1F(options.fitFrame+'NLL_lth_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                'NLL_lth_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                100,-1,1)
                phi_p_NLL = TH1F(options.fitFrame+'NLL_lphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                 'NLL_lphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                 100,-1,1)
                thphi_p_NLL = TH1F(options.fitFrame+'NLL_lthphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                   'NLL_lthphi_p_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                   100,-sqrt(2)/2,sqrt(2)/2)
                
                th_np_NLL = TH1F(options.fitFrame+'NLL_lth_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                 'NLL_lth_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                 100,-1,1)
                phi_np_NLL = TH1F(options.fitFrame+'NLL_lphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                  'NLL_lphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                  100,-1,1)
                thphi_np_NLL = TH1F(options.fitFrame+'NLL_lthphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                    'NLL_lthphi_np_rap'+str(rap_bin)+'_pT'+str(pt_bin+1),
                                    100,-sqrt(2)/2,sqrt(2)/2)

                ws.loadSnapshot('pol_snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

                data = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                PPdf = ws.factory('SIMUL::NoExtPPdf'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+
                                  '(mlRegion,promptSignal=promptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+','+
                                  'nonPromptSignal=nonpromptPol'+options.fitFrame+'_'+str(rap_bin)+'_'+str(pt_bin+1)+')')
                                  

                theNLL = PPdf.createNLL(data,
                                        ROOT.RooFit.Range('mlfit'),
                                        ROOT.RooFit.SplitRange(False),                                    
                                        ROOT.RooFit.ConditionalObservables(RooArgSet(ws.var('JpsictErr'))),
                                        #ROOT.RooFit.ExternalConstraints(constrs),
                                        ROOT.RooFit.NumCPU(3),
                                        ROOT.RooFit.Verbose(False))

                weight = theNLL.getVal()

                offset = weight

                print 'best fit : ',weight-offset, exp(-(weight-offset))
                
                th_p_PostProb.Fill(ws.var('lambda_theta_'+options.fitFrame+'_p').getVal(),exp(-(weight-offset)))
                phi_p_PostProb.Fill(ws.var('lambda_phi_'+options.fitFrame+'_p').getVal(),exp(-(weight-offset)))
                thphi_p_PostProb.Fill(ws.var('lambda_thetaphi_'+options.fitFrame+'_p').getVal(),exp(-(weight-offset)))
                ltilde = (ws.var('lambda_theta_'+options.fitFrame+'_p').getVal()+3*ws.var('lambda_phi_'+options.fitFrame+'_p').getVal())/(1-ws.var('lambda_phi_'+options.fitFrame+'_p').getVal())
                ltilde_p_PostProb.Fill(ltilde,exp(-(weight-offset)))

                
                th_np_PostProb.Fill(ws.var('lambda_theta_'+options.fitFrame+'_np').getVal(),exp(-(weight-offset)))
                phi_np_PostProb.Fill(ws.var('lambda_phi_'+options.fitFrame+'_np').getVal(),exp(-(weight-offset)))
                thphi_np_PostProb.Fill(ws.var('lambda_thetaphi_'+options.fitFrame+'_np').getVal(),exp(-(weight-offset)))
                ltilde = (ws.var('lambda_theta_'+options.fitFrame+'_np').getVal()+3*ws.var('lambda_phi_'+options.fitFrame+'_np').getVal())/(1-ws.var('lambda_phi_'+options.fitFrame+'_np').getVal())
                ltilde_np_PostProb.Fill(ltilde,exp(-(weight-offset)))

                for i in range(options.points):

                    th_p,phi_p,thphi_p = unitarityTriplet()
                    th_np,phi_np,thphi_np = unitarityTriplet()
                    
                    
                    ws.var('lambda_theta_'+options.fitFrame+'_p').setVal(th_p)
                    ws.var('lambda_phi_'+options.fitFrame+'_p').setVal(phi_p)
                    ws.var('lambda_thetaphi_'+options.fitFrame+'_p').setVal(thphi_p)
                    
                    ws.var('lambda_theta_'+options.fitFrame+'_np').setVal(th_np)
                    ws.var('lambda_phi_'+options.fitFrame+'_np').setVal(phi_np)
                    ws.var('lambda_thetaphi_'+options.fitFrame+'_np').setVal(thphi_np)    

                    weight = theNLL.getVal()

                    print i,' : ',weight-offset, exp(-(weight-offset))
                    print 'Prompt   : ',th_p,phi_p,thphi_p
                    print 'NonPrompt: ',th_np,phi_np,thphi_np

                    
                    th_p_PostProb.Fill(th_p,exp(-(weight-offset)))
                    phi_p_PostProb.Fill(phi_p,exp(-(weight-offset)))
                    thphi_p_PostProb.Fill(thphi_p,exp(-(weight-offset)))
                    ltilde = (th_p+3*phi_p)/(1-phi_p)
                    ltilde_p_PostProb.Fill(ltilde,exp(-(weight-offset)))
                    
                    
                    th_np_PostProb.Fill(th_np,exp(-(weight-offset)))
                    phi_np_PostProb.Fill(phi_np,exp(-(weight-offset)))
                    thphi_np_PostProb.Fill(thphi_np,exp(-(weight-offset)))
                    ltilde = (th_np+3*phi_np)/(1-phi_np)
                    ltilde_np_PostProb.Fill(ltilde,exp(-(weight-offset)))

                    #th_p_bin = th_p_NLL.FindBin(th_p)
                    #phi_p_bin = th_p_NLL.FindBin(phi_p)
                    #thphi_p_bin = th_p_NLL.FindBin(th_p)
                    
                    #th_np_bin = th_p_NLL.FindBin(th_np)
                    #phi_np_bin = th_p_NLL.FindBin(phi_np)
                    #thphi_np_bin = th_p_NLL.FindBin(thphi_np)

                    #th_p_NLL.Fill(th_p,weight)
                    #phi_p_NLL.Fill(phi_p,weight)
                    #thphi_p_NLL.Fill(thphi_p,weight)
                    
                    #th_np_NLL.Fill(th_np,weight)
                    #phi_np_NLL.Fill(phi_np,weight)
                    #thphi_np_NLL.Fill(thphi_np,weight)

                th_p_PostProb.Write()
                phi_p_PostProb.Write()
                thphi_p_PostProb.Write()
                ltilde_p_PostProb.Write()
                
                th_np_PostProb.Write()
                phi_np_PostProb.Write()
                thphi_np_PostProb.Write()
                ltilde_np_PostProb.Write()
                
                
                #th_p_NLL.Write()
                #phi_p_NLL.Write()
                #thphi_p_NLL.Write()
                
                #th_np_NLL.Write()
                #phi_np_NLL.Write()
                #thphi_np_NLL.Write()
                

        output.Close()

#bound returns 1 if unitarity is preserved, 0 otherwise
def unitarityTriplet(bound_th=1,bound_ph=1,bound_thph=sqrt(2)/2):

    th = bound_th*(2*random() - 1)
    phi = bound_ph*(2*random() - 1)
    thphi = bound_thph*(2*random() - 1)

    while not goodUnitarity(th,phi,thphi):
        #print 'triplet not good!'
        th = bound_th*(2*random() - 1)
        phi = bound_ph*(2*random() - 1)
        thphi = bound_thph*(2*random() - 1)
    
    return th,phi,thphi

def goodUnitarity(th,phi,thphi):
    if fabs(phi) > .5 *(1+th):
        return False
    if fabs(thphi) > .5*(1-phi):
        return False
    if th*th + 2*thphi*thphi > 1:
        return False
    if phi <= -1.0/3.0 and (pow(1.+2*phi,2) + 2*pow(thphi,2)) > 1:
        return  False
    return True

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
