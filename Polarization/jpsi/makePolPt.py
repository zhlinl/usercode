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

from ROOT import RooFit, RooWorkspace, RooArgSet, RooArgList, RooRealVar,gROOT
from ROOT import RooCategory, RooFormulaVar, RooDataSet, RooAddPdf, TGraphAsymmErrors
from ROOT import TFile, TTree, TChain, RooGExpModel, RooAddition, RooMinuit, TVectorD
from ROOT import RooPlot, TList


def main(options,args):
    gROOT.Reset()

    #load our super special Polarization PDF
    gROOT.ProcessLine('.L RooPolarizationPdf.cxx+')
    gROOT.ProcessLine('.L RooPolarizationConstraint.cxx+')

    #setup integration
    intConf = ROOT.RooAbsReal.defaultIntegratorConfig()
    #intConf.Print('v')
    #intConf.method1D().setLabel('RooAdaptiveGaussKronrodIntegrator1D')

    output = TFile.Open('PolVsPt'+options.fitFrame+'.root','recreate')
    output.cd()
    
    lth_pt_p = TGraphAsymmErrors()
    lphi_pt_p = TGraphAsymmErrors()
    lthphi_pt_p = TGraphAsymmErrors()
    lthtilde_pt_p = TGraphAsymmErrors()
    f_pt_p = TGraphAsymmErrors()

    lth_pt_p_list = TList()
    lphi_pt_p_list = TList()
    lthphi_pt_p_list = TList()
    
    lth_pt_p.SetName('lth_pt_p')
    lphi_pt_p.SetName('lphi_pt_p')
    lthphi_pt_p.SetName('lthphi_pt_p')
    lthtilde_pt_p.SetName('lthtilde_pt_p')
    f_pt_p.SetName('f_pt_p')

    lth_pt_np = TGraphAsymmErrors()
    lphi_pt_np = TGraphAsymmErrors()
    lthphi_pt_np = TGraphAsymmErrors()
    lthtilde_pt_np = TGraphAsymmErrors()
    f_pt_np = TGraphAsymmErrors()
    
    lth_pt_np_list = TList()
    lphi_pt_np_list = TList()
    lthphi_pt_np_list = TList()
    
    lth_pt_np.SetName('lth_pt_np')
    lphi_pt_np.SetName('lphi_pt_np')
    lthphi_pt_np.SetName('lthphi_pt_np')
    lthtilde_pt_np.SetName('lthtilde_pt_np')
    f_pt_np.SetName('f_pt_np')

    total_nonzero_points = 0

    for f in args:
        
        infile = TFile.Open(f)
        
        #ws = output.Get(f[:f.find('.root')].split('-')[0])
        frame = f[:f.find('.root')].split('-')[1]
        rapBin = f[:f.find('.root')].split('-')[2].split('_')[0]
        ptBin = f[:f.find('.root')].split('-')[2].split('_')[1]

        lth_p = infile.Get('lth_p_rap'+rapBin)

        for point in range(lth_p.GetN()):
            x = ROOT.Double(0)
            y = ROOT.Double(0)
            exl = float(0)
            exh = float(0)
            eyl = float(0)
            eyh = float(0)
            lth_p.GetPoint(point,x,y)
            if x != 0.0 and y != 0.0:

                infile.Get('lth_p_rap'+rapBin).GetPoint(point,x,y)
                lth_pt_p.SetPoint(total_nonzero_points,x,y)
                
                infile.Get('lphi_p_rap'+rapBin).GetPoint(point,x,y)
                lphi_pt_p.SetPoint(total_nonzero_points,x,y)
                
                infile.Get('lthphi_p_rap'+rapBin).GetPoint(point,x,y)
                lthphi_pt_p.SetPoint(total_nonzero_points,x,y)

                infile.Get('lthtilde_p_rap'+rapBin).GetPoint(point,x,y)
                lthtilde_pt_p.SetPoint(total_nonzero_points,x,y)

                infile.Get('f_p_rap'+rapBin).GetPoint(point,x,y)
                f_pt_p.SetPoint(total_nonzero_points,x,y)

                exh = infile.Get('lth_p_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lth_p_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lth_p_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lth_p_rap'+rapBin).GetErrorYlow(point)
                lth_pt_p.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)
                
                exh = infile.Get('lphi_p_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lphi_p_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lphi_p_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lphi_p_rap'+rapBin).GetErrorYlow(point)
                lphi_pt_p.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)
                
                exh = infile.Get('lthphi_p_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lthphi_p_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lthphi_p_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lthphi_p_rap'+rapBin).GetErrorYlow(point)
                lthphi_pt_p.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                exh = infile.Get('lthtilde_p_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lthtilde_p_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lthtilde_p_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lthtilde_p_rap'+rapBin).GetErrorYlow(point)
                lthtilde_pt_p.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                exh = infile.Get('f_p_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('f_p_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('f_p_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('f_p_rap'+rapBin).GetErrorYlow(point)
                f_pt_p.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                infile.Get('lth_np_rap'+rapBin).GetPoint(point,x,y)
                lth_pt_np.SetPoint(total_nonzero_points,x,y)
                
                infile.Get('lphi_np_rap'+rapBin).GetPoint(point,x,y)
                lphi_pt_np.SetPoint(total_nonzero_points,x,y)
                
                infile.Get('lthphi_np_rap'+rapBin).GetPoint(point,x,y)
                lthphi_pt_np.SetPoint(total_nonzero_points,x,y)

                infile.Get('lthtilde_np_rap'+rapBin).GetPoint(point,x,y)
                lthtilde_pt_np.SetPoint(total_nonzero_points,x,y)

                infile.Get('f_np_rap'+rapBin).GetPoint(point,x,y)
                f_pt_np.SetPoint(total_nonzero_points,x,y)

                exh = infile.Get('lth_np_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lth_np_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lth_np_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lth_np_rap'+rapBin).GetErrorYlow(point)
                lth_pt_np.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)
                
                exh = infile.Get('lphi_np_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lphi_np_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lphi_np_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lphi_np_rap'+rapBin).GetErrorYlow(point)
                lphi_pt_np.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                exh = infile.Get('lthphi_np_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lthphi_np_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lthphi_np_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lthphi_np_rap'+rapBin).GetErrorYlow(point)
                lthphi_pt_np.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                exh = infile.Get('lthtilde_np_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('lthtilde_np_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('lthtilde_np_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('lthtilde_np_rap'+rapBin).GetErrorYlow(point)
                lthtilde_pt_np.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                exh = infile.Get('f_np_rap'+rapBin).GetErrorXhigh(point)
                exl = infile.Get('f_np_rap'+rapBin).GetErrorXlow(point)
                eyh = infile.Get('f_np_rap'+rapBin).GetErrorYhigh(point)
                eyl = infile.Get('f_np_rap'+rapBin).GetErrorYlow(point)
                f_pt_np.SetPointError(total_nonzero_points,exl,exh,eyl,eyh)

                total_nonzero_points += 1

    output.cd()

       
    lth_pt_p.Write()
    lphi_pt_p.Write()
    lthphi_pt_p.Write()
    lthtilde_pt_p.Write()
    f_pt_p.Write()
    
    lth_pt_np.Write()
    lphi_pt_np.Write()
    lthphi_pt_np.Write()
    lthtilde_pt_np.Write()
    f_pt_np.Write()
    
    output.Close()


if __name__ == '__main__':
    parser = OptionParser(description='%prog : J/Psi Polarization Fitter.',
                          usage='polarizationFit.py --workspaceName=foo inputData.root ...')
    
    parser.add_option('--treeName',dest='treeName',default='data',help='Name of the input TTree.')
    parser.add_option('--plotPol',dest='plotPol',default=False,action='store_true',help='Make plots of polarization fit')
    parser.add_option('--testBin',dest='testBin',help='Only plot one bin.')
    parser.add_option('--ptBins',dest='ptBins',help='which pt bins to do')
    parser.add_option('--rapBins',dest='rapBins',help='which rapidity bins to do')
    parser.add_option('--fitFrame',dest='fitFrame',help='fit frame label')
    
    (options,args) = parser.parse_args()

    miss_options = False

    

    if options.testBin is not None:
        options.testBin = options.testBin.split(',')
        
    if miss_options:
        exit(1)

    main(options,args)
