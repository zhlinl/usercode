import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from YZheng.UpsilonAna.tagAndProbe import *

process.TnP = Template.clone(
    InputFileNames = cms.vstring("../EventSelection/MuOnia.root"),
    OutputFileName = cms.string("effTQ.root"),
    Efficiencies = cms.PSet(
        TQ = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("TQ","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                #pt = cms.vdouble(2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 20.0),
                pt = cms.vdouble(2.5, 50.0),
                abseta = cms.vdouble(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4),
                #abseta = cms.vdouble(0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4),
                TM = cms.vstring("pass"),
                L1DoubleMuOpen_Tight03 = cms.vstring("pass"),
                tag_L1DoubleMuOpen_Tight03 = cms.vstring("pass"),
                pair_drM2 = cms.vdouble(0.6, 5),
                # run = cms.vdouble(0, 141900),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ),
        TQint = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("TQ","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(2.5, 50.0),
                abseta = cms.vdouble(0.0, 2.4),
                TM = cms.vstring("pass"),
                L1DoubleMuOpen_Tight03 = cms.vstring("pass"),
                tag_L1DoubleMuOpen_Tight03 = cms.vstring("pass"),
                pair_drM2 = cms.vdouble(0.6, 5),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ),
    ),
)

process.fit = cms.Path(
    process.TnP
)

