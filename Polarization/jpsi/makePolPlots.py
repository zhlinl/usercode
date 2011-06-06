#!/usr/bin/env python

# implementation of the full polarization fit using a RooSimultaneous
# should be much more stable than what we had previously
# nuisance parameters will be fit without polarization first and
# then constrained during polarization fit.

import os
from optparse import OptionParser
import ROOT
import commonVar as jpsi

from math import sqrt,fabs

from ROOT import RooFit, RooWorkspace, RooArgSet, RooArgList, RooRealVar,gROOT
from ROOT import RooCategory, RooFormulaVar, RooDataSet, RooAddPdf, TGraphAsymmErrors
from ROOT import TFile, TTree, TChain, RooGExpModel, RooAddition, RooMinuit, TVectorD
from ROOT import RooPlot


def main(options,args):
    gROOT.Reset()

    #load our super special Polarization PDF
    #gROOT.ProcessLine('.L RooPolarizationPdf.cxx+')
    #gROOT.ProcessLine('.L RooPolarizationConstraint.cxx+')

    #setup integration
    intConf = ROOT.RooAbsReal.defaultIntegratorConfig()
    #intConf.setEpsAbs(1e-13)
    #intConf.setEpsRel(1e-13)
    #intConf.Print('v')
    intConf.method1D().setLabel('RooGaussKronrodIntegrator1D')
    #intConf.method2D().setLabel('RooIntegrator2D')
    #intConf.method2D().setLabel()
    
    #intConf.methodND().setLabel('RooMCIntegrator')

    #intConf.getConfigSection('RooIntegrator1D').find('extrapolation').setLabel('None')
    #intConf.getConfigSection('RooIntegrator1D').find('maxSteps').setVal(500)
    #intConf.getConfigSection('RooIntegrator1D').find('minSteps').setVal(10)
    #intConf.getConfigSection('RooIntegrator1D').find('fixSteps').setVal(100)
    
    for f in args:

        output = TFile.Open(f,'UPDATE')
        output.cd()

        print f[:f.find('.root')].split('-')[0]
        
        ws = output.Get(f[:f.find('.root')].split('-')[0])

        for rap_bin in range(1,len(jpsi.pTRange)):        
            ptMean = TVectorD(len(jpsi.pTRange[rap_bin]))
            errptMeanl = TVectorD(len(jpsi.pTRange[rap_bin]))
            errptMeanh = TVectorD(len(jpsi.pTRange[rap_bin]))
            bFraction = TVectorD(len(jpsi.pTRange[rap_bin]))
            bFractionErrl = TVectorD(len(jpsi.pTRange[rap_bin]))
            bFractionErrh = TVectorD(len(jpsi.pTRange[rap_bin]))

            ptMean_p = TVectorD(len(jpsi.pTRange[rap_bin]))
            errptMeanl_p = TVectorD(len(jpsi.pTRange[rap_bin]))
            errptMeanh_p = TVectorD(len(jpsi.pTRange[rap_bin]))

            ptMean_np = TVectorD(len(jpsi.pTRange[rap_bin]))
            errptMeanl_np = TVectorD(len(jpsi.pTRange[rap_bin]))
            errptMeanh_np = TVectorD(len(jpsi.pTRange[rap_bin]))

            th_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_th_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_th_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            
            phi_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_phi_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_phi_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            
            thphi_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_thphi_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_thphi_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))

            thtilde_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_thtilde_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_thtilde_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))

            f_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_f_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_f_p_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            
            th_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_th_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_th_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            
            phi_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_phi_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_phi_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            
            thphi_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_thphi_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_thphi_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))

            thtilde_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_thtilde_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_thtilde_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))

            f_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errl_f_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            errh_f_np_Pt = TVectorD(len(jpsi.pTRange[rap_bin]))
            
            for pt_bin in range(len(jpsi.pTRange[rap_bin])):

                if options.testBin is not None:
                    if rap_bin != int(options.testBin[0]) or pt_bin+1 != int(options.testBin[1]):
                        continue
                
                #datasets
                d = ws.data('data_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                dsig = d.reduce('massRegion == massRegion::signal')
                dpr = d.reduce('mlRegion == mlRegion::promptSignal')
                dnp = d.reduce('mlRegion == mlRegion::nonPromptSignal')
                dbkg = d.reduce('massRegion == massRegion::rightMassSideBand | massRegion == massRegion::leftMassSideBand')
                dbkgl = d.reduce('massRegion == massRegion::leftMassSideBand')
                dbkgr = d.reduce('massRegion == massRegion::rightMassSideBand')

                #make plots of mass and lifetime for signal and background
                ctauFrameSig = ws.var('Jpsict').frame()
                ctauFrameBkg = ws.var('Jpsict').frame()
                massFrame = ws.var('JpsiMass').frame()

                ctauFrameSig.SetName('ctausig_plot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                ctauFrameBkg.SetName('ctaubkg_plot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                massFrame.SetName('mass_plot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

                ws.loadSnapshot('pol_snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                
                dsig.plotOn(ctauFrameSig)
                dbkg.plotOn(ctauFrameBkg)
                d.plotOn(massFrame)

                ws.pdf('MPdf').plotOn(massFrame,
                                      RooFit.ProjWData(d))
                
                ws.pdf('MPdf').plotOn(massFrame,
                                      RooFit.Components('backgroundExtMass'),
                                      RooFit.LineStyle(7),
                                      RooFit.LineColor(ROOT.kPink+3),
                                      RooFit.LineWidth(2),
                                      RooFit.ProjWData(d))
                
                ws.pdf('MPdf').plotOn(massFrame,
                                      RooFit.Components('nonpromptExtMass'),
                                      RooFit.LineStyle(2),
                                      RooFit.LineColor(ROOT.kRed),
                                      RooFit.LineWidth(2),
                                      RooFit.ProjWData(d))

                ws.pdf('MPdf').plotOn(massFrame,
                                      RooFit.Components('promptExtMass'),
                                      RooFit.LineStyle(5),
                                      RooFit.LineColor(ROOT.kBlue),
                                      RooFit.LineWidth(2),
                                      RooFit.ProjWData(d))
                

                ws.pdf('LPdf').plotOn(ctauFrameSig,
                                      RooFit.ProjWData(dsig))
                
                ws.pdf('LPdf').plotOn(ctauFrameSig,
                                      RooFit.Components('promptExtCTau'),
                                       RooFit.LineStyle(5),
                                      RooFit.LineColor(ROOT.kBlue),
                                      RooFit.LineWidth(2),
                                      RooFit.ProjWData(dsig))

                ws.pdf('LPdf').plotOn(ctauFrameSig,
                                      RooFit.Components('nonpromptExtCTau'),
                                      RooFit.LineStyle(2),
                                      RooFit.LineColor(ROOT.kRed),
                                      RooFit.LineWidth(2),
                                      RooFit.ProjWData(dsig))

                ws.pdf('LPdf').plotOn(ctauFrameSig,
                                      RooFit.Components('backgroundExtCTau'),
                                      RooFit.LineStyle(7),
                                      RooFit.LineColor(ROOT.kPink+3),
                                      RooFit.LineWidth(2),
                                      RooFit.ProjWData(dsig))                
                
                ws.pdf('LPdf').plotOn(ctauFrameBkg,
                                      RooFit.ProjWData(dbkg))

                #make polarization plots
                if options.plotPol:
                    dataVars = d.get()
                    
                    phiFramePr = None
                    costhFramePr = None
                    phiFrameNP = None
                    costhFrameNP = None

                    
                    phiFrameNPPeda = None
                    costhFrameNPPeda = None
                    
                    phiFrameLPeda = None
                    costhFrameLPeda = None
                    phiFrameRPeda = None
                    costhFrameRPeda = None
                    
                    fitFrame = None
                    di =  dataVars.createIterator()
                    
                    while (di.Next()):
                        name = di.GetName()                    
                        if 'phi' in name:
                            di.setBins(25)
                            phiFramePr = di.frame()
                            phiFrameNP = di.frame()
                            phiFramePrPeda = di.frame()
                            phiFrameNPPeda = di.frame()
                            di.setBins(12)
                            phiFrameLPeda = di.frame()
                            phiFrameRPeda = di.frame()
                            fitFrame = name.split('_')[-1]
                        if 'costh' in name:
                            di.setBins(25)
                            costhFramePr = di.frame()
                            costhFrameNP = di.frame()
                            costhFramePrPeda = di.frame()
                            costhFrameNPPeda = di.frame()
                            di.setBins(12)
                            costhFrameLPeda = di.frame()
                            costhFrameRPeda = di.frame()

                    if phiFramePr is None:
                        print 'ARRGG PHI'
                    if costhFramePr is None:
                        print 'ARRGG COSTH'
                    
                    dpr.plotOn(phiFramePr)
                    dpr.plotOn(costhFramePr)
                    
                    dnp.plotOn(phiFrameNP)
                    dnp.plotOn(costhFrameNP)

                    #prompt region phi
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFramePr,
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFramePr,
                                                     RooFit.Components('polPdf'+fitFrame+'Prompt_'+str(rap_bin)+
                                                                       '_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(5),
                                                     RooFit.LineColor(ROOT.kBlue),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFramePr,
                                                     RooFit.Components('polPdf'+fitFrame+'NonPrompt_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(2),
                                                     RooFit.LineColor(ROOT.kRed),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFramePr,
                                                     RooFit.Components('bkgShape'+fitFrame+'_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(7),
                                                     RooFit.LineColor(ROOT.kPink+3),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    #prompt region costh
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFramePr,
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFramePr,
                                                     RooFit.Components('polPdf'+fitFrame+'Prompt_'+str(rap_bin)+
                                                                       '_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(5),
                                                     RooFit.LineColor(ROOT.kBlue),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFramePr,
                                                     RooFit.Components('polPdf'+fitFrame+'NonPrompt_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(2),
                                                     RooFit.LineColor(ROOT.kRed),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFramePr,
                                                     RooFit.Components('bkgShape'+fitFrame+'_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(7),
                                                     RooFit.LineColor(ROOT.kPink+3),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                        

                    #non-prompt region phi
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFrameNP,
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFrameNP,
                                                     RooFit.Components('polPdf'+fitFrame+'Prompt_'+str(rap_bin)+
                                                                       '_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(5),
                                                     RooFit.LineColor(ROOT.kBlue),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFrameNP,
                                                     RooFit.Components('polPdf'+fitFrame+'NonPrompt_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(2),
                                                     RooFit.LineColor(ROOT.kRed),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(phiFrameNP,
                                                     RooFit.Components('bkgShape'+fitFrame+'_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(7),
                                                     RooFit.LineColor(ROOT.kPink+3),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    #non-prompt region costh
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFrameNP,
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFrameNP,
                                                     RooFit.Components('polPdf'+fitFrame+'Prompt_'+str(rap_bin)+
                                                                       '_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(5),
                                                     RooFit.LineColor(ROOT.kBlue),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFrameNP,
                                                     RooFit.Components('polPdf'+fitFrame+'NonPrompt_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(2),
                                                     RooFit.LineColor(ROOT.kRed),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                    ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                           '_'+str(pt_bin+1)).plotOn(costhFrameNP,
                                                     RooFit.Components('bkgShape'+fitFrame+'_'+
                                                                       str(rap_bin)+'_'+str(pt_bin+1)),
                                                     RooFit.LineStyle(7),
                                                     RooFit.LineColor(ROOT.kPink+3),
                                                     RooFit.LineWidth(2),
                                                     RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))                    

                    if options.pedagogical:

                        dpr.plotOn(phiFramePrPeda)
                        dpr.plotOn(costhFramePrPeda)
                        
                        dnp.plotOn(phiFrameNPPeda)
                        dnp.plotOn(costhFrameNPPeda)

                        dbkgl.plotOn(phiFrameLPeda)
                        dbkgl.plotOn(costhFrameLPeda)
                        
                        dbkgr.plotOn(phiFrameRPeda)                        
                        dbkgr.plotOn(costhFrameRPeda)

                        #prompt region phi
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))
                        #prompt region costh
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr))

                        #non-prompt region phi
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))

                        #non-prompt region costh
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp))
                        
                        ws.var('lambda_theta_'+fitFrame+'_p').setVal(0.0)
                        ws.var('lambda_phi_'+fitFrame+'_p').setVal(0.0)
                        ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)
                        
                        ws.var('lambda_theta_'+fitFrame+'_np').setVal(0.0)
                        ws.var('lambda_phi_'+fitFrame+'_np').setVal(0.0)
                        ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)
                    
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameLPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameLPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameRPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameRPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))                        
                        
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))
                        
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(1),
                                                         RooFit.LineWidth(1))

                        if fitFrame == "CS":
                            ws.var('lambda_theta_'+fitFrame+'_p').setVal(1.0)
                            ws.var('lambda_phi_'+fitFrame+'_p').setVal(0.0)
                            ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)
                        elif fitFrame == "HX":
                            ws.var('lambda_theta_'+fitFrame+'_p').setVal(0.0)
                            ws.var('lambda_phi_'+fitFrame+'_p').setVal(0.5)
                            ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)
                        
                        ws.var('lambda_theta_'+fitFrame+'_np').setVal(0.0)
                        ws.var('lambda_phi_'+fitFrame+'_np').setVal(0.0)
                        ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)
                    

                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameLPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameLPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameRPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameRPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))   
                        
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))

                        ws.var('lambda_theta_'+fitFrame+'_p').setVal(0.0)
                        ws.var('lambda_phi_'+fitFrame+'_p').setVal(0.0)
                        ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)

                        if fitFrame == "CS":
                            ws.var('lambda_theta_'+fitFrame+'_np').setVal(1.0)
                            ws.var('lambda_phi_'+fitFrame+'_np').setVal(0.0)
                            ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)
                        elif fitFrame == "HX":
                            ws.var('lambda_theta_'+fitFrame+'_np').setVal(0.0)
                            ws.var('lambda_phi_'+fitFrame+'_np').setVal(0.5)
                            ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)
                        
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(7),
                                                         RooFit.LineWidth(1))


                        if fitFrame == "CS":
                            ws.var('lambda_theta_'+fitFrame+'_p').setVal(-1.0)
                            ws.var('lambda_phi_'+fitFrame+'_p').setVal(0.0)
                            ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)
                        elif fitFrame == "HX":
                            ws.var('lambda_theta_'+fitFrame+'_p').setVal(0.0)
                            ws.var('lambda_phi_'+fitFrame+'_p').setVal(-0.5)
                            ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)
                        
                        ws.var('lambda_theta_'+fitFrame+'_np').setVal(0.0)
                        ws.var('lambda_phi_'+fitFrame+'_np').setVal(0.0)
                        ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)

                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameLPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameLPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameRPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))
                        ws.pdf('polPdf'+fitFrame+'Prompt_'+
                               str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameRPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))   
                        
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFramePrPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dpr),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))

                        ws.var('lambda_theta_'+fitFrame+'_p').setVal(0.0)
                        ws.var('lambda_phi_'+fitFrame+'_p').setVal(0.0)
                        ws.var('lambda_thetaphi_'+fitFrame+'_p').setVal(0.0)

                        if fitFrame == "CS":
                            ws.var('lambda_theta_'+fitFrame+'_np').setVal(-1.0)
                            ws.var('lambda_phi_'+fitFrame+'_np').setVal(0.0)
                            ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)
                        elif fitFrame == "HX":
                            ws.var('lambda_theta_'+fitFrame+'_np').setVal(0.0)
                            ws.var('lambda_phi_'+fitFrame+'_np').setVal(-0.5)
                            ws.var('lambda_thetaphi_'+fitFrame+'_np').setVal(0.0)
                        
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(phiFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))
                        ws.pdf('PPdf'+fitFrame+'_'+str(rap_bin)+
                               '_'+str(pt_bin+1)).plotOn(costhFrameNPPeda,
                                                         RooFit.ProjWData(RooArgSet(ws.cat('mlRegion')),dnp),
                                                         RooFit.LineColor(ROOT.kRed),
                                                         RooFit.LineStyle(2),
                                                         RooFit.LineWidth(1))
                        phiFramePrPeda.SetName('polPhi'+fitFrame+'Prompt_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        costhFramePrPeda.SetName('polCosTh'+fitFrame+'Prompt_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        phiFrameNPPeda.SetName('polPhi'+fitFrame+'NonPrompt_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        costhFrameNPPeda.SetName('polCosTh'+fitFrame+'NonPrompt_peda_'+str(rap_bin)+'_'+str(pt_bin+1))

                        phiFrameLPeda.SetName('polPhi'+fitFrame+'L_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        costhFrameLPeda.SetName('polCosTh'+fitFrame+'L_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        phiFrameRPeda.SetName('polPhi'+fitFrame+'R_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        costhFrameRPeda.SetName('polCosTh'+fitFrame+'R_peda_'+str(rap_bin)+'_'+str(pt_bin+1))
                        
                        phiFramePrPeda.Write()
                        costhFramePrPeda.Write()
                        phiFrameNPPeda.Write()
                        costhFrameNPPeda.Write()

                        phiFrameLPeda.Write()
                        costhFrameLPeda.Write()
                        phiFrameRPeda.Write()
                        costhFrameRPeda.Write()
                        #end pedagogical

                    phiFramePr.SetName('polPhi'+fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1))
                    costhFramePr.SetName('polCosTh'+fitFrame+'Prompt_'+str(rap_bin)+'_'+str(pt_bin+1))
                    phiFrameNP.SetName('polPhi'+fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1))
                    costhFrameNP.SetName('polCosTh'+fitFrame+'NonPrompt_'+str(rap_bin)+'_'+str(pt_bin+1))                    
                                        
                    phiFramePr.Write()
                    costhFramePr.Write()
                    phiFrameNP.Write()
                    costhFrameNP.Write()                    

                    fr = ws.obj('pol_fitresult_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))

                    #th_ph_ErrPlot_p = ROOT.RooPlot(ws.var('lambda_theta_'+fitFrame+'_p'),
                    #                               ws.var('lambda_phi_'+fitFrame+'_p'),-1.,1,-1.,1.)
                   
                    #th_thph_ErrPlot_p = ROOT.RooPlot(ws.var('lambda_theta_'+fitFrame+'_p'),
                    #                                 ws.var('lambda_thetaphi_'+fitFrame+'_p'),-1.,1,-1.,1.)
                    
                    #thph_ph_ErrPlot_p = ROOT.RooPlot(ws.var('lambda_phi_'+fitFrame+'_p'),
                    #                                 ws.var('lambda_thetaphi_'+fitFrame+'_p'),-1.,1,-1.,1.)

                    #th_ph_ErrPlot_np = ROOT.RooPlot(ws.var('lambda_theta_'+fitFrame+'_np'),
                    #                                ws.var('lambda_phi_'+fitFrame+'_np'),-1.,1,-1.,1.)
                   
                    #th_thph_ErrPlot_np = ROOT.RooPlot(ws.var('lambda_theta_'+fitFrame+'_np'),
                    #                                  ws.var('lambda_thetaphi_'+fitFrame+'_np'),-1.,1,-1.,1.)
                    
                    #thph_ph_ErrPlot_np = ROOT.RooPlot(ws.var('lambda_phi_'+fitFrame+'_np'),
                    #                                  ws.var('lambda_thetaphi_'+fitFrame+'_np'),-1.,1,-1.,1.)

                    #fr.plotOn(th_ph_ErrPlot_p,
                    #          ws.var('lambda_theta_'+fitFrame+'_p'),
                    #          ws.var('lambda_phi_'+fitFrame+'_p'),'ME12')
                    #fr.plotOn(th_thph_ErrPlot_p,
                    #          ws.var('lambda_theta_'+fitFrame+'_p'),
                    #          ws.var('lambda_thetaphi_'+fitFrame+'_p'),'ME12')
                    #fr.plotOn(thph_ph_ErrPlot_p,
                    #          ws.var('lambda_phi_'+fitFrame+'_p'),
                    #          ws.var('lambda_thetaphi_'+fitFrame+'_p'),'ME12')

                    #fr.plotOn(th_ph_ErrPlot_np,
                    #          ws.var('lambda_theta_'+fitFrame+'_np'),
                    #          ws.var('lambda_phi_'+fitFrame+'_np'),'ME12')
                    #fr.plotOn(th_thph_ErrPlot_np,
                    #          ws.var('lambda_theta_'+fitFrame+'_np'),
                    #          ws.var('lambda_thetaphi_'+fitFrame+'_np'),'ME12')
                    #fr.plotOn(thph_ph_ErrPlot_np,
                    #          ws.var('lambda_phi_'+fitFrame+'_np'),
                    #          ws.var('lambda_thetaphi_'+fitFrame+'_np'),'ME12')

                    #th_ph_ErrPlot_p.SetName('th_ph_p_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    #th_thph_ErrPlot_p.SetName('th_thph_p_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    #thph_ph_ErrPlot_p.SetName('thph_ph_p_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    
                    #th_ph_ErrPlot_np.SetName('th_ph_np_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    #th_thph_ErrPlot_np.SetName('th_thph_np_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    #thph_ph_ErrPlot_np.SetName('thph_ph_np_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    
                    #th_ph_ErrPlot_p.Write()
                    #th_thph_ErrPlot_p.Write()
                    #thph_ph_ErrPlot_p.Write()
                    
                    #th_ph_ErrPlot_np.Write()
                    #th_thph_ErrPlot_np.Write()
                    #thph_ph_ErrPlot_np.Write()

                    

                    #reload snapshot
                    ws.loadSnapshot('pol_snapshot_rap'+str(rap_bin)+'_pt'+str(pt_bin+1))
                    
                    #do prompt
                    ptMeanThisBin = 0
                    for i in range(dpr.numEntries()):
                        set = dpr.get(i)
                        ptMeanThisBin += set.find('JpsiPt').getVal()
                    ptMeanThisBin /= dpr.numEntries()
                    ptMean_p[pt_bin] = ptMeanThisBin
                    errptMeanl_p[pt_bin] = ptMeanThisBin - jpsi.pTRange[rap_bin][pt_bin][0]
                    errptMeanh_p[pt_bin] = jpsi.pTRange[rap_bin][pt_bin][-1] - ptMeanThisBin

                    th_p_Pt[pt_bin] = ws.var('lambda_theta_'+fitFrame+'_p').getVal()                    
                    thphi_p_Pt[pt_bin] = ws.var('lambda_thetaphi_'+fitFrame+'_p').getVal()

                    errl_th_p_Pt[pt_bin] = -ws.var('lambda_theta_'+fitFrame+'_p').getErrorLo()
                    errl_thphi_p_Pt[pt_bin] = -ws.var('lambda_thetaphi_'+fitFrame+'_p').getErrorLo()

                    errh_th_p_Pt[pt_bin] = ws.var('lambda_theta_'+fitFrame+'_p').getErrorHi()
                    errh_thphi_p_Pt[pt_bin] = ws.var('lambda_thetaphi_'+fitFrame+'_p').getErrorHi()

                    #other thing to plot, depending on mode
                    if ws.set('POI').find('lambda_phi_'+fitFrame+'_p'):
                        phi_p_Pt[pt_bin] = ws.var('lambda_phi_'+fitFrame+'_p').getVal()
                        errl_phi_p_Pt[pt_bin] = -ws.var('lambda_phi_'+fitFrame+'_p').getErrorLo()
                        errh_phi_p_Pt[pt_bin] = ws.var('lambda_phi_'+fitFrame+'_p').getErrorHi()                        
                    if ws.set('POI').find('lambda_tilde_p'):
                        thtilde_p_Pt[pt_bin] = ws.var('lambda_tilde_p').getVal()                    
                        errl_thtilde_p_Pt[pt_bin] = -ws.var('lambda_tilde_p').getErrorLo()
                        errh_thtilde_p_Pt[pt_bin] = ws.var('lambda_tilde_p').getErrorHi()
                    if ws.set('POI').find('F_invar_p'):
                        f_p_Pt[pt_bin] = ws.var('F_invar_p').getVal()                    
                        errl_f_p_Pt[pt_bin] = -ws.var('F_invar_p').getErrorLo()
                        errh_f_p_Pt[pt_bin] = ws.var('F_invar_p').getErrorHi()
                    

                    #do non-prompt
                    ptMeanThisBin = 0
                    for i in range(dnp.numEntries()):
                        set = dnp.get(i)
                        ptMeanThisBin += set.find('JpsiPt').getVal()
                    ptMeanThisBin /= dnp.numEntries()
                    ptMean_np[pt_bin] = ptMeanThisBin
                    errptMeanl_np[pt_bin] = ptMeanThisBin - jpsi.pTRange[rap_bin][pt_bin][0]
                    errptMeanh_np[pt_bin] = jpsi.pTRange[rap_bin][pt_bin][-1] - ptMeanThisBin

                    th_np_Pt[pt_bin] = ws.var('lambda_theta_'+fitFrame+'_np').getVal()                                        
                    thphi_np_Pt[pt_bin] = ws.var('lambda_thetaphi_'+fitFrame+'_np').getVal()

                    errl_th_np_Pt[pt_bin] = -ws.var('lambda_theta_'+fitFrame+'_np').getErrorLo()
                    errl_thphi_np_Pt[pt_bin] = -ws.var('lambda_thetaphi_'+fitFrame+'_np').getErrorLo()

                    errh_th_np_Pt[pt_bin] = ws.var('lambda_theta_'+fitFrame+'_np').getErrorHi()
                    errh_thphi_np_Pt[pt_bin] = ws.var('lambda_thetaphi_'+fitFrame+'_np').getErrorHi()

                    #other thing to plot, depending on mode
                    if ws.set('POI').find('lambda_phi_'+fitFrame+'_np'):
                        phi_np_Pt[pt_bin] = ws.var('lambda_phi_'+fitFrame+'_np').getVal()
                        errl_phi_np_Pt[pt_bin] = -ws.var('lambda_phi_'+fitFrame+'_np').getErrorLo()
                        errh_phi_np_Pt[pt_bin] = ws.var('lambda_phi_'+fitFrame+'_np').getErrorHi()                        
                    if ws.set('POI').find('lambda_tilde_np'):
                        thtilde_np_Pt[pt_bin] = ws.var('lambda_tilde_np').getVal()                    
                        errl_thtilde_np_Pt[pt_bin] = -ws.var('lambda_tilde_np').getErrorLo()
                        errh_thtilde_np_Pt[pt_bin] = ws.var('lambda_tilde_np').getErrorHi()
                    if ws.set('POI').find('F_invar_np'):
                        f_np_Pt[pt_bin] = ws.var('F_invar_np').getVal()                    
                        errl_f_np_Pt[pt_bin] = -ws.var('F_invar_np').getErrorLo()
                        errh_f_np_Pt[pt_bin] = ws.var('F_invar_np').getErrorHi()

                #make b-fraction diagnostic plots
                ptMeanThisBin = 0
                rapMeanThisBin = 0
                for i in range(dsig.numEntries()):
                    set = dsig.get(i)
                    ptMeanThisBin += set.find('JpsiPt').getVal()
                    rapMeanThisBin += fabs(set.find('JpsiRap').getVal())
                rapMeanThisBin /=dsig.numEntries()
                print rapMeanThisBin
                ptMeanThisBin /= dsig.numEntries()
                ptMean[pt_bin] = ptMeanThisBin
                errptMeanl[pt_bin] = ptMeanThisBin - jpsi.pTRange[rap_bin][pt_bin][0]
                errptMeanh[pt_bin] = jpsi.pTRange[rap_bin][pt_bin][-1] - ptMeanThisBin
                                

                nNP = ws.var('nNonPromptSignal').getVal()
                nP = ws.var('nPromptSignal').getVal()

                errNP = ws.var('nNonPromptSignal').getError()
                errP = ws.var('nPromptSignal').getError()

                bFraction[pt_bin] = (nNP/(nP+nNP))
                bFractionErrl[pt_bin] = ((1.0/(nNP+nP)+nNP/((nNP+nP)*(nNP+nP)))*errNP + nNP/((nNP+nP)*(nNP+nP))*errP)
                bFractionErrh[pt_bin] = bFractionErrl[pt_bin]

                ctauFrameSig.Write()
                ctauFrameBkg.Write()
                massFrame.Write()

            if options.plotPol:
                lth_p_pt = TGraphAsymmErrors(ptMean_p,th_p_Pt,errptMeanl_p,errptMeanh_p,errl_th_p_Pt,errh_th_p_Pt)
                lth_p_pt.SetName('lth_p_rap'+str(rap_bin))
                lth_p_pt.Write()
                
                lphi_p_pt = TGraphAsymmErrors(ptMean_p,phi_p_Pt,errptMeanl_p,errptMeanh_p,errl_phi_p_Pt,errh_phi_p_Pt)
                lphi_p_pt.SetName('lphi_p_rap'+str(rap_bin))
                lphi_p_pt.Write()
                
                lthphi_p_pt = TGraphAsymmErrors(ptMean_p,thphi_p_Pt,errptMeanl_p,errptMeanh_p,errl_thphi_p_Pt,errh_thphi_p_Pt)
                lthphi_p_pt.SetName('lthphi_p_rap'+str(rap_bin))
                lthphi_p_pt.Write()
                
                lthtilde_p_pt = TGraphAsymmErrors(ptMean_p,thtilde_p_Pt,errptMeanl_p,errptMeanh_p,errl_thtilde_p_Pt,errh_thtilde_p_Pt)
                lthtilde_p_pt.SetName('lthtilde_p_rap'+str(rap_bin))
                lthtilde_p_pt.Write()

                f_p_pt = TGraphAsymmErrors(ptMean_p,f_p_Pt,errptMeanl_p,errptMeanh_p,errl_f_p_Pt,errh_f_p_Pt)
                f_p_pt.SetName('f_p_rap'+str(rap_bin))
                f_p_pt.Write()
                
                lth_np_pt = TGraphAsymmErrors(ptMean_np,th_np_Pt,errptMeanl_np,errptMeanh_np,errl_th_np_Pt,errh_th_np_Pt)
                lth_np_pt.SetName('lth_np_rap'+str(rap_bin))
                lth_np_pt.Write()
                
                lphi_np_pt = TGraphAsymmErrors(ptMean_np,phi_np_Pt,errptMeanl_np,errptMeanh_np,errl_phi_np_Pt,errh_phi_np_Pt)
                lphi_np_pt.SetName('lphi_np_rap'+str(rap_bin))
                lphi_np_pt.Write()

                lthphi_np_pt = TGraphAsymmErrors(ptMean_np,thphi_np_Pt,errptMeanl_np,errptMeanh_np,errl_thphi_np_Pt,errh_thphi_np_Pt)
                lthphi_np_pt.SetName('lthphi_np_rap'+str(rap_bin))
                lthphi_np_pt.Write()            
                
                lthtilde_p_pt = TGraphAsymmErrors(ptMean_np,thtilde_np_Pt,errptMeanl_np,errptMeanh_np,errl_thtilde_np_Pt,errh_thtilde_np_Pt)
                lthtilde_p_pt.SetName('lthtilde_np_rap'+str(rap_bin))
                lthtilde_p_pt.Write()

                f_np_pt = TGraphAsymmErrors(ptMean_np,f_np_Pt,errptMeanl_np,errptMeanh_np,errl_f_np_Pt,errh_f_np_Pt)
                f_np_pt.SetName('f_np_rap'+str(rap_bin))
                f_np_pt.Write()
                
            bFracGraph = TGraphAsymmErrors(ptMean,bFraction,errptMeanl,errptMeanh,bFractionErrl,bFractionErrh)
            bFracGraph.SetName('bFrac_rap'+str(rap_bin))
            bFracGraph.Write()


    output.Close()


if __name__ == '__main__':
    parser = OptionParser(description='%prog : J/Psi Polarization Fitter.',
                          usage='polarizationFit.py --workspaceName=foo inputData.root ...')
    
    parser.add_option('--treeName',dest='treeName',default='data',help='Name of the input TTree.')
    parser.add_option('--plotPol',dest='plotPol',default=False,action='store_true',help='Make plots of polarization fit')
    parser.add_option('--testBin',dest='testBin',help='Only plot one bin.')
    parser.add_option('--pedagogical',dest='pedagogical',default=False,action='store_true',help='make +1,0,-1 lambda theta and lambda phi plots')
    (options,args) = parser.parse_args()

    miss_options = False

    if options.testBin is not None:
        options.testBin = options.testBin.split(',')
        
    if miss_options:
        exit(1)

    main(options,args)
