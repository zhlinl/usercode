import FWCore.ParameterSet.Config as cms

process = cms.Process("SELECTION")

from YZheng.UpsilonAna.selection_cff import *
selection(process, GlobalTag="START38_V8::All", MC=True, SelectionTrigger="hltL1MuOpenL1Filtered0")

process.source.fileNames = cms.untracked.vstring(
    '/store/user/zgecse/BpToJPsiMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/Summer10-START36_V9_S09-v1-Onia2MuMu-v5/bef246ddbd6b4e7664a45838bf80a640/onia2MuMuPAT_9_1_xqy.root',
)   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

