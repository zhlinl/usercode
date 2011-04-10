import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from zhlinl.UpsilonAna.tagAndProbe import *

process.TnP = Template.clone(
    InputFileNames = cms.vstring("/uscms_data/d2/zgecse/Onia2MuMu-v5/Upsilon-v8/CMSSW_3_8_1/src/YZheng/UpsilonAna/test/EventSelection/MuOnia.root"),
    OutputFileName = cms.string("effTrig.root"),
    Efficiencies = cms.PSet(
        L1DoubleMuOpen_Tight = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("L1DoubleMuOpen_Tight03","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 50.0),
                #pt = cms.vdouble(5.0, 20.0),
                # abseta = cms.vdouble(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4),
                abseta = cms.vdouble(0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4),
                TQ = cms.vstring("pass"),
                tag_MuX_TrackY_Jpsi_MU = cms.vstring("pass"),
                pair_drM2 = cms.vdouble(0.6, 5),
                TM = cms.vstring("pass"),
                # run = cms.vdouble(0, 141900),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ),
    ),
)

process.fit = cms.Path(
    process.TnP
)

