from ROOT import TFile, gROOT, gStyle, TGraphAsymmErrors, gPad, TMultiGraph, \
TGraph, TGraphErrors, kFullCircle, kBlue, TLatex, TLegend, TCanvas

def normalizeGraph(graph, binwidth):
    xsecSum = 0.
    for pt in range(1, graph.GetN()):
        xsecSum += graph.GetY()[pt]*binwidth[pt]
    normGraph = TGraphAsymmErrors(graph.GetN()-1);
    for pt in range(0, normGraph.GetN()):
        normGraph.SetPoint(pt, graph.GetX()[pt+1], graph.GetY()[pt+1]/xsecSum)
        normGraph.SetPointEYhigh(pt, graph.GetErrorYhigh(pt+1)/xsecSum)
        normGraph.SetPointEYlow(pt, graph.GetErrorYlow(pt+1)/xsecSum)
    return normGraph

gROOT.ProcessLine('.L ../xsection/setTDRStyle_modified.C')
from ROOT import setTDRStyle
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetLineScalePS(6.)
gStyle.SetPadLeftMargin(0.18)
gStyle.SetTitleYOffset(1.35)
gStyle.SetTextFont(42)
bins1s = [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,17.,20.,30.]
bins2s = [0., 2., 4., 6., 9.,12.,16.,20.,30.]
bins3s = [0., 3., 6., 9., 14., 20., 30.]

Tevatron = TFile('TevatronResults.root')

l = TLatex()
l.SetNDC()
l.SetTextFont(42)
l.SetTextAlign(31)
l.SetTextSize(0.06)

for resonance in range(1, 4):
    fname = 'xsection_%is_rap_0_2.root' % resonance
    print fname
    cmsfile = TFile(fname)
    mg1 = cmsfile.Get('mg1')
    c = TCanvas('c_Y%iS' % resonance, 'c_Y%iS' % resonance, 800, 800)

    unnormCMS = mg1.GetListOfGraphs().At(1)
    if (resonance == 1):
        bins = bins1s
        D0 = Tevatron.Get('D0_Y%iS' % resonance)
        CDF = Tevatron.Get('CDF_Y%iS' % resonance)
    elif (resonance == 2):
        bins = bins2s
        D0 = None
        CDF = Tevatron.Get('CDF_Y%iS' % resonance)
    elif (resonance == 3):
        bins = bins3s
        D0 = None
        CDF = Tevatron.Get('CDF_Y%iS' % resonance)
    binwidths = []
    for i in range(0, unnormCMS.GetN()):
        if (i > 0):
            binwidths.append(bins[i]-bins[i-1])
        else:
            binwidths.append(30.0)
    #print binwidths
    normCMS = normalizeGraph(unnormCMS, binwidths)
    normCMS.SetMarkerStyle(kFullCircle)
    normCMS.SetMarkerSize(1.4)
    normCMS.SetMarkerColor(kBlue+2)
    normCMS.SetLineColor(kBlue+2)
    normCMS.Draw('apz')
    normCMS.SetMaximum(0.3)
    leg = TLegend(0.4, 1. - gPad.GetTopMargin() - 0.03 - 0.18,
                  1. - gPad.GetRightMargin() - 0.02,
                  1. - gPad.GetTopMargin() - 0.03,
                  '', 'NDC')
    leg.SetMargin(0.15)
    #leg.Dump()
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(normCMS, 'CMS, |y|<2, #sqrt{s}=7 TeV', 'pe')
    if (D0):
        D0.Draw('pz')
        leg.AddEntry(D0, 'D#oslash, |y|<1.8, #sqrt{s}=1.96 TeV', 'pe')
    else:
        leg.SetY1(1. - gPad.GetTopMargin() - 0.03 - 0.12)
    if (CDF):
        CDF.Draw('pz')
        leg.AddEntry(CDF, 'CDF, |y|<0.4, #sqrt{s}=1.8 TeV', 'pe')
    #leg.Dump()
    leg.Draw('same')
    l.DrawLatex(1.0 - gPad.GetRightMargin() - 0.04,
            leg.GetY1() - 0.06,
            '#varUpsilon(%iS)' % resonance)
    gPad.SetLogy()
    gPad.Update()
    normCMS.GetXaxis().SetTitle('p_{T} (GeV/c)')
    normCMS.GetYaxis().SetTitle('(d#sigma/dp_{T})/#sigma_{TOT} (GeV/c)^{-1}')
    normCMS.GetXaxis().SetLimits(0., 30.)
    gPad.Modified()
    gPad.Print('TevatronCompare%iS.eps' % resonance)
    gPad.Print('TevatronCompare%iS.png' % resonance)
    gPad.Print('TevatronCompare%iS.pdf' % resonance)
    
