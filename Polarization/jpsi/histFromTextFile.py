import ROOT,os

from optparse import OptionParser
from ROOT import TH2F

def main(options,args):

    out = ROOT.TFile.Open(options.output,'RECREATE')
    
    for fs in args:
        f = open(fs,'r')

        ptBin = fs.split('.txt')[0].split('pT')[-1][0]
        yBin = fs.split('.txt')[0].split('rap')[-1][0]
        frame = fs.split('_')[1]

        print ptBin,yBin,frame
        
        outHist = TH2F('hAcc2D_Onia_'+frame+'_pT'+ptBin+'_rap'+yBin,'converted text file',
                       int(options.costhBins),-1,1,
                       int(options.phiBins),0,360)

        for line in f:
            (z,costh,phi) = line.split()

            outHist.SetBinContent(outHist.FindBin(float(costh),
                                                  float(phi)),
                                  float(z))
        out.Add(outHist.Clone(),True)
            
        out.Write()

        #for line in f:
        #    z,costh,phi = line.split()

if __name__ == "__main__":

    parser = OptionParser(description="takes a list of formatted files and turns them into TH2Fs")

    parser.add_option("--output",dest="output",help="Name of output file")
    parser.add_option("--phiBins",dest="phiBins",help="Number of bins in phi")
    parser.add_option("--costhBins",dest="costhBins",help="Number of bins in cos(theta)")

    (options,args) = parser.parse_args()
    
    main(options,args)
